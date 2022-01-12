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
        String backend
    }

    Map[String, String] acronym2uri = read_map(acronym_file)
    File genome_file = if sub(genome, "^.+\\.(tgz|gz)$", "PATH") == "PATH" then genome else acronym2uri[genome]

    String whitelist_uri = if assay != '' then acronym2uri[assay] else 'null'

    call run_starsolo {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            read1_fastq_pattern = read1_fastq_pattern,
            read2_fastq_pattern = read2_fastq_pattern,
            assay = assay,
            genome = genome_file,
            output_directory = output_directory,
            outSAMtype = outSAMtype,
            soloType = soloType,
            soloCBwhitelist = if whitelist_uri != 'null' then whitelist_uri else soloCBwhitelist,
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
            backend = backend
    }

    output {
        File monitoringLog = run_starsolo.monitoringLog
        String output_folder = run_starsolo.output_folder
    }
}

task run_starsolo {
    input {
        String sample_id
        String input_fastqs_directories
        String read1_fastq_pattern
        String read2_fastq_pattern
        String assay
        File genome
        String output_directory
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
        String version
        String zones
        String memory
        Int num_cpu
        Int disk_space
        Int preemptible
        Int awsMaxRetries
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

        call_args = ['STAR', '--genomeDir', 'genome_ref', '--runThreadN', '~{num_cpu}']

        if '~{assay}' in ['tenX_v2', 'tenX_v3', 'ShareSeq']:
            call_args.extend(['--soloType', 'CB_UMI_Simple', '--soloCBmatchWLtype', '1MM_multi_Nbase_pseudocounts', '--soloUMIfiltering', 'MultiGeneUMI_CR', \
                              '--soloUMIdedup', '1MM_CR', '--clipAdapterType', 'CellRanger4', '--outFilterScoreMin', '30', \
                              '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMattributes', 'CR', 'UR', 'CY', 'UY', 'CB', 'UB'])
            if '~{assay}' == 'tenX_v3':
                call_args.extend(['--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '12'])
            elif '~{assay}' == 'tenX_v2':
                call_args.extend(['--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '10'])
            elif '~{assay}' == 'ShareSeq':
                call_args.extend(['--soloCBstart', '1', '--soloCBlen', '24', '--soloUMIstart', '25', '--soloUMIlen', '10'])
        elif '~{assay}' in ['SeqWell', 'DropSeq']:
            call_args.extend(['--soloType', 'CB_UMI_Simple', '--soloCBstart', '1', '--soloCBlen', '12', '--soloUMIstart', '13', '--soloUMIlen', '8', \
                              '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMattributes', 'CR', 'UR', 'CY', 'UY', 'CB', 'UB'])
        else:
            call_args.extend(['--outSAMtype', 'BAM', 'Unsorted'])

        if file_ext == '.fastq.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        call_args.extend(['--readFilesIn', ','.join(r2_list), ','.join(r1_list)])
        call_args.extend(['--outFileNamePrefix', 'result/~{sample_id}_'])

        if '~{outSAMtype}' != '':
            call_args.extend(['--outSAMtype'] + remove_extra_space('~{outSAMtype}').split(' '))
        if '~{soloType}' != '':
            call_args.extend(['--soloType', '~{soloType}'])

        if '~{soloCBwhitelist}' != '':
            fn_tup = os.path.splitext("~{soloCBwhitelist}")
            if fn_tup[1] == '.gz':
                with open(fn_tup[0], 'w') as fp:
                    check_call(['zcat', "~{soloCBwhitelist}"], stdout=fp)
                whitelist_file = fn_tup[0]
            else:
                whitelist_file = '~{soloCBwhitelist}'
            call_args.extend(['--soloCBwhitelist', whitelist_file])
        else:
            call_args.extend(['--soloCBwhitelist', 'None'])

        if '~{soloCBstart}' != '':
            call_args.extend(['--soloCBstart', '~{soloCBstart}'])
        if '~{soloCBlen}' != '':
            call_args.extend(['--soloCBlen', '~{soloCBlen}'])
        if '~{soloUMIstart}' != '':
            call_args.extend(['--soloUMIstart', '~{soloUMIstart}'])
        if '~{soloUMIlen}' != '':
            call_args.extend(['--soloUMIlen', '~{soloUMIlen}'])
        if '~{soloBarcodeReadLength}' != '':
            call_args.extend(['--soloBarcodeReadLength', '~{soloBarcodeReadLength}'])
        if '~{soloBarcodeMate}' != '':
            call_args.extend(['--soloBarcodeMate', '~{soloBarcodeMate}'])
        if '~{soloType}' == 'CB_UMI_Complex':
            if '~{soloCBposition}' != '':
                call_args.extend(['--soloCBposition'] + remove_extra_space('~{soloCBposition}').split(' '))
            if '~{soloUMIposition}' != '':
                call_args.extend(['--soloUMIposition'] + remove_extra_space('~{soloUMIposition}').split(' '))
        if '~{soloAdapterSequence}' != '':
            call_args.extend(['--soloAdapterSequence', '~{soloAdapterSequence}'])
        if '~{soloAdapterMismatchesNmax}' != '':
            call_args.extend(['--soloAdapterMismatchesNmax', '~{soloAdapterMismatchesNmax}'])
        if '~{soloCBmatchWLtype}' != '':
            call_args.extend(['--soloCBmatchWLtype', '~{soloCBmatchWLtype}'])
        if '~{soloInputSAMattrBarcodeSeq}' != '':
            call_args.extend(['--soloInputSAMattrBarcodeSeq'] + remove_extra_space('~{soloInputSAMattrBarcodeSeq}').split(' '))
        if '~{soloInputSAMattrBarcodeQual}' != '':
            call_args.extend(['--soloInputSAMattrBarcodeQual'] + remove_extra_space('~{soloInputSAMattrBarcodeQual}').split(' '))
        if '~{soloStrand}' != '':
            call_args.extend(['--soloStrand', '~{soloStrand}'])
        if '~{soloFeatures}' != '':
            feature_list = remove_extra_space('~{soloFeatures}').split(' ')
            if ('Velocyto' in feature_list) and ('Gene' not in feature_list):
                feature_list.append('Gene')
            call_args.extend(['--soloFeatures'] + feature_list)
        if '~{soloMultiMappers}' != '':
            call_args.extend(['--soloMultiMappers'] + remove_extra_space('~{soloMultiMappers}').split(' '))
        if '~{soloUMIdedup}' != '':
            call_args.extend(['--soloUMIdedup'] + remove_extra_space('~{soloUMIdedup}').split(' '))
        if '~{soloUMIfiltering}' != '':
            call_args.extend(['--soloUMIfiltering'] + remove_extra_space('~{soloUMIfiltering}').split(' '))
        if '~{soloCellFilter}' != '':
            call_args.extend(['--soloCellFilter'] + remove_extra_space('~{soloCellFilter}').split(' '))
        if '~{soloOutFormatFeaturesGeneField3}' != '':
            call_args.extend(['--soloOutFormatFeaturesGeneField3'] + remove_extra_space('~{soloOutFormatFeaturesGeneField3}').split(' '))

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m result "~{output_directory}/~{sample_id}"
    }

    output {
        File monitoringLog = 'monitoring.log'
        String output_folder = '~{output_directory}/~{sample_id}'
    }

    runtime {
        docker: "~{docker_registry}/starsolo:~{version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
