version 1.0

workflow starsolo_per_sample {
    input {
        String sample_id
        String r1_fastqs
        String r2_fastqs
        String? preset
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
        String star_version
        String zones
        String memory
        Int num_cpu
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    call run_starsolo {
        input:
            sample_id = sample_id,
            r1_fastqs = r1_fastqs,
            r2_fastqs = r2_fastqs,
            preset = preset,
            genome = genome,
            output_directory = output_directory,
            outSAMtype = outSAMtype,
            soloType = soloType,
            soloCBwhitelist = soloCBwhitelist,
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
        String r1_fastqs
        String r2_fastqs
        String? preset
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
        from subprocess import check_call

        def remove_extra_space(s):
            return re.sub(' +', ' ', s.strip())

        r1_list = '~{r1_fastqs}'.split(',')
        r2_list = '~{r2_fastqs}'.split(',')
        assert len(r1_list) == len(r2_list)

        def fetch_and_rename_fastq(in_file, out_file):
            call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', in_file, out_file]
            print(' '.join(call_args))
            check_call(call_args)

        for i in range(len(r1_list)):
            file_ext = '.fastq.gz' if os.path.splitext(r1_list[i])[-1] == '.gz' else '.fastq'
            fetch_and_rename_fastq(r1_list[i], 'R1_'+str(i)+file_ext)
            fetch_and_rename_fastq(r2_list[i], 'R2_'+str(i)+file_ext)

        call_args = ['STAR', '--genomeDir', 'genome_ref', '--runThreadN', '~{num_cpu}']

        if '~{preset}' in ['tenX_v2', 'tenX_v3']:
            call_args.extend(['--soloType', 'CB_UMI_Simple', '--soloCBmatchWLtype', '1MM_multi_Nbase_pseudocounts', '--soloUMIfiltering', 'MultiGeneUMI_CR', \
                              '--soloUMIdedup', '1MM_CR', '--clipAdapterType', 'CellRanger4', '--outFilterScoreMin', '30', \
                              '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMattributes', 'CR', 'UR', 'CY', 'UY', 'CB', 'UB'])
            if '~{preset}' == 'tenX_v3':
                call_args.extend(['--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '12'])
            elif '~{preset}' == 'tenX_v2':
                call_args.extend(['--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '10'])
        elif '~{preset}' in ['SeqWell', 'DropSeq']:
            call_args.extend(['--soloType', 'CB_UMI_Simple', '--soloCBwhitelist', 'None', '--soloCBstart', '1', '--soloCBlen', '12', '--soloUMIstart', '13', '--soloUMIlen', '8', \
                              '--outSAMtype', 'BAM', 'SortedByCoordinate', '--outSAMattributes', 'CR', 'UR', 'CY', 'UY', 'CB', 'UB'])
        else:
            call_args.extend(['--outSAMtype', 'BAM', 'Unsorted'])

        if file_ext == '.fastq.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        def set_up_readfiles(prefix, n_files, f_ext):
            return ','.join([prefix+str(n)+f_ext for n in range(n_files)])

        call_args.extend(['--readFilesIn', set_up_readfiles('R2_', len(r2_list), file_ext), set_up_readfiles('R1_', len(r1_list), file_ext)])
        call_args.extend(['--outFileNamePrefix', 'result/~{sample_id}_'])

        if '~{outSAMtype}' != '':
            call_args.extend(['--outSAMtype'] + remove_extra_space('~{outSAMtype}').split(' '))
        if '~{soloType}' != '':
            call_args.extend(['--soloType', '~{soloType}'])
        if '~{soloCBwhitelist}' != '':
            call_args.extend(['--soloCBwhitelist'] + remove_extra_space('~{soloCBwhitelist}').split(' '))
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
