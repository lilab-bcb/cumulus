version 1.0

task run_starsolo_per_sample {
    input {
        String sample_id
        String r1_fastqs
        String r2_fastqs
        String preset
        File genome
        String output_directory
        String soloType
        File? preset_whitelist
        File? soloCBwhitelist
        Int? soloCBstart
        Int? soloCBlen
        Int? soloUMIstart
        Int? soloUMIlen
        Int? soloBarcodeReadLength
        Int? soloBarcodeMate
        Int? soloCBposition
        Int? soloUMIposition
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
        String? soloOutFileNames
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
        import os
        from subprocess import check_call

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

        call_args = ['STAR', '--soloType', '~{soloType}', '--genomeDir', 'genome_ref', '--runThreadN', '~{num_cpu}', '--outSAMtype', 'BAM', 'SortedByCoordinate', '--clipAdapterType', 'CellRanger4', '--outFilterScoreMin', '30']

        white_list = 'None'
        if '~{preset_whitelist}' != '':
            white_list = '~{preset_whitelist}'
        else:
            if '~{preset}' == 'custom' and '~{soloCBwhitelist}' != '':
                white_list = '~{soloCBwhitelist}'

        if '~{preset}' in ['tenX_v2', 'tenX_v3']:
            call_args.extend(['--soloType', 'CB_UMI_Simple', '--soloCBmatchWLtype', '1MM_multi_Nbase_pseudocounts', '--soloUMIfiltering', 'MultiGeneUMI_CR', '--soloUMIdedup', '1MM_CR'])
            if '~{preset}' == 'tenX_v3':
                call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '12'])
            elif '~{preset}' == 'tenX_v2':
                call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '10'])
        elif '~{preset}' in ['SeqWell', 'DropSeq']:
            call_args.extend(['--soloType', 'CB_UMI_Simple', '--soloCBwhitelist', 'None', '--soloCBstart', '1', '--soloCBlen', '12', '--soloUMIstart', '13', '--soloUMIlen', '8'])
        elif '~{preset}' == 'custom':
            call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '~{soloCBstart}', '--soloCBlen', '~{soloCBlen}', '--soloUMIstart', '~{soloUMIstart}', '--soloUMIlen', '~{soloUMIlen}'])

        if file_ext == '.fastq.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        def set_up_readfiles(prefix, n_files, f_ext):
            return ','.join([prefix+str(n)+f_ext for n in range(n_files)])

        call_args.extend(['--readFilesIn', set_up_readfiles('R2_', len(r2_list), file_ext), set_up_readfiles('R1_', len(r1_list), file_ext)])
        call_args.extend(['--outFileNamePrefix', 'result/~{sample_id}_'])

        if '~{soloType}' != '':
            call_args.extend(['--soloType', '~{soloType}'])
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
                call_args.extend(['--soloCBposition'] + '~{soloCBposition}'.split(','))
            if '~{soloUMIposition}' != '':
                call_args.extend(['--soloUMIposition', '~{soloUMIposition}'])
        if '~{soloAdapterSequence}' != '':
            call_args.extend(['--soloAdapterSequence', '~{soloAdapterSequence}'])
        if '~{soloAdapterMismatchesNmax}' != '':
            call_args.extend(['--soloAdapterMismatchesNmax', '~{soloAdapterMismatchesNmax}'])
        if '~{soloCBmatchWLtype}' != '':
            call_args.extend(['--soloCBmatchWLtype', '~{soloCBmatchWLtype}'])
        if '~{soloInputSAMattrBarcodeSeq}' != '':
            call_args.extend(['--soloInputSAMattrBarcodeSeq'] + '~{soloInputSAMattrBarcodeSeq}'.split(','))
        if '~{soloInputSAMattrBarcodeQual}' != '':
            call_args.extend(['--soloInputSAMattrBarcodeQual'] + '~{soloInputSAMattrBarcodeQual}'.split(','))
        if '~{soloStrand}' != '':
            call_args.extend(['--soloStrand', '~{soloStrand}'])
        if '~{soloFeatures}' != '':
            feature_list = '~{soloFeatures}'.split(',')
            if ('Velocyto' in feature_list) and ('Gene' not in feature_list):
                feature_list.append('Gene')
            call_args.extend(['--soloFeatures'] + feature_list)
        if '~{soloMultiMappers}' != '':
            call_args.extend(['--soloMultiMappers'] + '~{soloMultiMappers}'.split(','))
        if '~{soloUMIdedup}' != '':
            call_args.extend(['--soloUMIdedup'] + '~{soloUMIdedup}'.split(','))
        if '~{soloUMIfiltering}' != '':
            call_args.extend(['--soloUMIfiltering'] + '~{soloUMIfiltering}'.split(','))
        if '~{soloOutFileNames}' != '':
            call_args.extend(['--soloOutFileNames'] + '~{soloOutFileNames}'.split(','))
        if '~{soloCellFilter}' != '':
            call_args.extend(['--soloCellFilter'] + '~{soloCellFilter}'.split(','))
        if '~{soloOutFormatFeaturesGeneField3}' != '':
            call_args.extend(['--soloOutFormatFeaturesGeneField3'] + '~{soloOutFormatFeaturesGeneField3}'.split(','))

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
