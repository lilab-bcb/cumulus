version 1.0

workflow starsolo_count {
    input {
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (local paths or URIs)
        String input_fastqs_directories
        String assay = ""
        String genome
        String read1_fastq_pattern = "_S*_L*_R1_001.fastq.gz"
        String read2_fastq_pattern = "_S*_L*_R2_001.fastq.gz"
        String barcode_read
        String output_directory
        File acronym_file
        String? outSAMtype
        String? soloType
        File? soloCBwhitelist
        Int? soloCBstart
        Int? soloCBlen
        Int? soloUMIstart
        Int? soloUMIlen
        Int? soloBarcodeReadLength
        Int? soloBarcodeMate
        String? soloCBposition
        String? soloUMIposition
        String? soloAdapterSequence
        Int? soloAdapterMismatchesNmax
        String? soloCBmatchWLtype
        String? soloInputSAMattrBarcodeSeq
        String? soloInputSAMattrBarcodeQual
        String? soloStrand
        String? soloFeatures
        String? soloMultiMappers
        String? soloUMIdedup
        String? soloUMIfiltering
        String? soloCellFilter
        String? soloOutFormatFeaturesGeneField3
        String docker_registry
        String star_version
        String zones
        String memory
        Int num_cpu
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String awsQueueArn
        String backend
    }

    Map[String, String] acronym2uri = read_map(acronym_file)
    File genome_file = if sub(genome, "^.+\\.(tgz|gz)$", "PATH") == "PATH" then genome else acronym2uri[genome]

    String whitelist_uri = if assay != '' then acronym2uri[assay] else acronym2uri["Nil"]

    call run_starsolo {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            read1_fastq_pattern = read1_fastq_pattern,
            read2_fastq_pattern = read2_fastq_pattern,
            barcode_read = barcode_read,
            assay = assay,
            genome = genome_file,
            output_directory = output_directory,
            outSAMtype = outSAMtype,
            soloType = soloType,
            soloCBwhitelist = if defined(soloCBwhitelist) then select_first([soloCBwhitelist]) else whitelist_uri,
            soloCBstart = soloCBstart,
            soloCBlen = soloCBlen,
            soloUMIstart = soloUMIstart,
            soloUMIlen = soloUMIlen,
            soloBarcodeReadLength = soloBarcodeReadLength,
            soloBarcodeMate = soloBarcodeMate,
            soloCBposition = soloCBposition,
            soloUMIposition = soloUMIposition,
            soloAdapterSequence =soloAdapterSequence,
            soloAdapterMismatchesNmax = soloAdapterMismatchesNmax,
            soloCBmatchWLtype = soloCBmatchWLtype,
            soloInputSAMattrBarcodeSeq = soloInputSAMattrBarcodeSeq,
            soloInputSAMattrBarcodeQual = soloInputSAMattrBarcodeQual,
            soloStrand = soloStrand,
            soloFeatures = soloFeatures,
            soloMultiMappers = soloMultiMappers,
            soloUMIdedup = soloUMIdedup,
            soloUMIfiltering = soloUMIfiltering,
            soloCellFilter = soloCellFilter,
            soloOutFormatFeaturesGeneField3 = soloOutFormatFeaturesGeneField3,
            docker_registry = docker_registry,
            version = star_version,
            zones = zones,
            memory = memory,
            num_cpu = num_cpu,
            disk_space = disk_space,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    output {
        File monitoringLog = run_starsolo.monitoringLog
        File starsoloLog = run_starsolo.starsoloLog
        String output_count_directory = run_starsolo.output_count_directory
    }
}

task run_starsolo {
    input {
        String sample_id
        String input_fastqs_directories
        String read1_fastq_pattern
        String read2_fastq_pattern
        String barcode_read
        String assay
        File genome
        String output_directory
        String? outSAMtype
        String? soloType
        File soloCBwhitelist
        Int? soloCBstart
        Int? soloCBlen
        Int? soloUMIstart
        Int? soloUMIlen
        Int? soloBarcodeReadLength
        Int? soloBarcodeMate
        String? soloCBposition
        String? soloUMIposition
        String? soloAdapterSequence
        Int? soloAdapterMismatchesNmax
        String? soloCBmatchWLtype
        String? soloInputSAMattrBarcodeSeq
        String? soloInputSAMattrBarcodeQual
        String? soloStrand
        String? soloFeatures
        String? soloMultiMappers
        String? soloUMIdedup
        String? soloUMIfiltering
        String? soloCellFilter
        String? soloOutFormatFeaturesGeneField3
        String docker_registry
        String version
        String zones
        String memory
        Int num_cpu
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        mkdir genome_ref
        tar -zxf "~{genome}" -C genome_ref --strip-components 1
        rm "~{genome}"

        mkdir result

        python <<CODE
        import os, re
        from fnmatch import fnmatch
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT

        def set_up_input_fastq_files(l, folder, pattern):
            file_list = [f for f in os.listdir(folder) if fnmatch(f, "*" + pattern)]
            file_list.sort()
            for f in file_list:
                l.append(folder + '/' + f)

        def generate_args_list(args_dict):
            res_list = list()
            for k, v in args_dict.items():
                if isinstance(v, list):
                    res_list.extend([k] + v)
                else:
                    res_list.extend([k, v])
            return res_list

        r1_list = list()
        r2_list = list()
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            target = "~{sample_id}_" + str(i)
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)

                set_up_input_fastq_files(r1_list, target, '~{read1_fastq_pattern}')
                set_up_input_fastq_files(r2_list, target, '~{read2_fastq_pattern}')

            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{read1_fastq_pattern}' , target + '/']
                print(' '.join(call_args))
                check_call(call_args)
                set_up_input_fastq_files(r1_list, target, '~{read1_fastq_pattern}')

                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{read2_fastq_pattern}' , target + '/']
                print(' '.join(call_args))
                check_call(call_args)
                set_up_input_fastq_files(r2_list, target, '~{read2_fastq_pattern}')

        assert len(r1_list) == len(r2_list)
        file_ext = '.fastq.gz' if os.path.splitext("~{read1_fastq_pattern}")[-1] == '.gz' else '.fastq'

        def remove_extra_space(s):
            return re.sub(' +', ' ', s.strip())

        call_args = ['STAR', '--genomeDir', 'genome_ref', '--runThreadN', '~{num_cpu}', '--outFileNamePrefix', 'result/']

        barcode_read = '~{barcode_read}'
        args_dict = dict()

        if '~{assay}' in ['tenX_v2', 'tenX_v3', 'ShareSeq', 'tenX_5p', 'tenX_5p_pe', 'tenX_multiome']:
            args_dict['--soloType'] = 'CB_UMI_Simple'
            args_dict['--soloCBmatchWLtype'] = '1MM_multi_Nbase_pseudocounts'
            args_dict['--soloUMIfiltering'] = 'MultiGeneUMI_CR'
            args_dict['--soloUMIdedup'] = '1MM_CR'
            args_dict['--outFilterScoreMin'] = '30'
            args_dict['--outSAMtype'] = ['BAM', 'SortedByCoordinate']
            args_dict['--outSAMattributes'] = ['CR', 'UR', 'CY', 'UY', 'CB', 'UB']

            if '~{assay}' in ['tenX_v3', 'tenX_multiome']:
                args_dict['--soloCBstart'] = '1'
                args_dict['--soloCBlen'] = '16'
                args_dict['--soloUMIstart'] = '17'
                args_dict['--soloUMIlen'] = '12'
                args_dict['--clipAdapterType'] = 'CellRanger4'
                barcode_read = 'read1'
            elif '~{assay}' in ['tenX_v2', 'tenX_5p', 'tenX_5p_pe']:
                args_dict['--soloCBstart'] = '1'
                args_dict['--soloCBlen'] = '16'
                args_dict['--soloUMIstart'] = '17'
                args_dict['--soloUMIlen'] = '10'

                if '~{assay}' == 'tenX_v2':
                    args_dict['--clipAdapterType'] = 'CellRanger4'
                    barcode_read = 'read1'
                elif '~{assay}' == 'tenX_5p':
                    barcode_read = 'read1'
                    args_dict['--soloStrand'] = 'Reverse'
                else:
                    args_dict['--soloBarcodeMate'] = '1'
                    args_dict['--clip5pNbases'] = ['39', '0']
                    barcode_read = 'read2'
            elif '~{assay}' == 'ShareSeq':
                args_dict['--soloCBstart'] = '1'
                args_dict['--soloCBlen'] = '24'
                args_dict['--soloUMIstart'] = '25'
                args_dict['--soloUMIlen'] = '10'
                barcode_read = 'read2'
        elif '~{assay}' in ['SeqWell', 'DropSeq']:
            args_dict['--soloType'] = 'CB_UMI_Simple'
            args_dict['--soloCBstart'] = '1'
            args_dict['--soloCBlen'] = '12'
            args_dict['--soloUMIstart'] = '13'
            args_dict['--soloUMIlen'] = '8'
            args_dict['--outSAMtype'] = ['BAM', 'SortedByCoordinate']
            args_dict['--outSAMattributes'] = ['CR', 'UR', 'CY', 'UY', 'CB', 'UB']
            barcode_read = 'read1'
        else:
            args_dict['--outSAMattributes'] = ['BAM', 'Unsorted']

        if file_ext == '.fastq.gz':
            args_dict['--readFilesCommand'] = 'zcat'

        if barcode_read == 'read1':
            args_dict['--readFilesIn'] = [','.join(r2_list), ','.join(r1_list)]
        else:
            args_dict['--readFilesIn'] = [','.join(r1_list), ','.join(r2_list)]

        if '~{outSAMtype}' != '':
            args_dict['--outSAMtype'] = remove_extra_space('~{outSAMtype}').split(' ')
        if '~{soloType}' != '':
            args_dict['--soloType'] = '~{soloType}'

        if '~{soloCBwhitelist}' != '' and os.path.basename('~{soloCBwhitelist}') != 'null':
            fn_tup = os.path.splitext("~{soloCBwhitelist}")
            if fn_tup[1] == '.gz':
                with open(fn_tup[0], 'w') as fp:
                    check_call(['zcat', "~{soloCBwhitelist}"], stdout=fp)
                whitelist_file = fn_tup[0]
            else:
                whitelist_file = '~{soloCBwhitelist}'
            args_dict['--soloCBwhitelist'] = whitelist_file
        else:
            args_dict['--soloCBwhitelist'] = 'None'

        if '~{soloCBstart}' != '':
            args_dict['--soloCBstart'] = '~{soloCBstart}'
        if '~{soloCBlen}' != '':
            args_dict['--soloCBlen'] = '~{soloCBlen}'
        if '~{soloUMIstart}' != '':
            args_dict['--soloUMIstart'] = '~{soloUMIstart}'
        if '~{soloUMIlen}' != '':
            args_dict['--soloUMIlen'] = '~{soloUMIlen}'
        if '~{soloBarcodeReadLength}' != '':
            args_dict['--soloBarcodeReadLength'] = '~{soloBarcodeReadLength}'
        if '~{soloBarcodeMate}' != '':
            args_dict['--soloBarcodeMate'] = '~{soloBarcodeMate}'
        if '~{soloType}' == 'CB_UMI_Complex':
            if '~{soloCBposition}' != '':
                args_dict['--soloCBposition'] = remove_extra_space('~{soloCBposition}').split(' ')
            if '~{soloUMIposition}' != '':
                args_dict['--soloUMIposition'] = remove_extra_space('~{soloUMIposition}').split(' ')
        if '~{soloAdapterSequence}' != '':
            args_dict['--soloAdapterSequence'] = '~{soloAdapterSequence}'
        if '~{soloAdapterMismatchesNmax}' != '':
            args_dict['--soloAdapterMismatchesNmax'] = '~{soloAdapterMismatchesNmax}'
        if '~{soloCBmatchWLtype}' != '':
            args_dict['--soloCBmatchWLtype'] = '~{soloCBmatchWLtype}'
        if '~{soloInputSAMattrBarcodeSeq}' != '':
            args_dict['--soloInputSAMattrBarcodeSeq'] = remove_extra_space('~{soloInputSAMattrBarcodeSeq}').split(' ')
        if '~{soloInputSAMattrBarcodeQual}' != '':
            args_dict['--soloInputSAMattrBarcodeQual'] = remove_extra_space('~{soloInputSAMattrBarcodeQual}').split(' ')
        if '~{soloStrand}' != '':
            args_dict['--soloStrand'] = '~{soloStrand}'
        if '~{soloFeatures}' != '':
            feature_list = remove_extra_space('~{soloFeatures}').split(' ')
            if ('Velocyto' in feature_list) and ('Gene' not in feature_list):
                feature_list.append('Gene')
            args_dict['--soloFeatures'] = feature_list
        if '~{soloMultiMappers}' != '':
            args_dict['--soloMultiMappers'] = remove_extra_space('~{soloMultiMappers}').split(' ')
        if '~{soloUMIdedup}' != '':
            args_dict['--soloUMIdedup'] = remove_extra_space('~{soloUMIdedup}').split(' ')
        if '~{soloUMIfiltering}' != '':
            args_dict['--soloUMIfiltering'] = remove_extra_space('~{soloUMIfiltering}').split(' ')
        if '~{soloCellFilter}' != '':
            args_dict['--soloCellFilter'] = remove_extra_space('~{soloCellFilter}').split(' ')
        if '~{soloOutFormatFeaturesGeneField3}' != '':
            args_dict['--soloOutFormatFeaturesGeneField3'] = remove_extra_space('~{soloOutFormatFeaturesGeneField3}').split(' ')

        call_args += generate_args_list(args_dict)

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m result "~{output_directory}/~{sample_id}"
    }

    output {
        File monitoringLog = 'monitoring.log'
        File starsoloLog = 'result/Log.out'
        String output_count_directory = '~{output_directory}/~{sample_id}'
    }

    runtime {
        docker: "~{docker_registry}/starsolo:~{version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
        queueArn: awsQueueArn
    }
}
