version 1.0

task run_starsolo_per_sample {
    input {
        String sample_id
        String r1_fastqs
        String r2_fastqs
        String chemistry
        File genome
        String output_directory
        String solo_type
        File? preset_whitelist
        File? CBwhitelist
        Int? CBstart
        Int? CBlen
        Int? UMIstart
        Int? UMIlen
        Int? BarcodeReadLength
        Int? BarcodeMate
        Int? CBposition
        Int? UMIposition
        String? adapterSequence
        Int? adapterMismatchesNmax
        String? CBmatchWLtype
        String? inputSAMattrBarcodeSeq
        String? inputSAMattrBarcodeQual
        String? strand
        String? features
        String? multiMappers
        String? UMIdedup
        String? UMIfiltering
        String? outFileNames
        String? cellFilter
        String? outFormatFeaturesGeneField3
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

        call_args = ['STAR', '--soloType', '~{solo_type}', '--genomeDir', 'genome_ref', '--runThreadN', '~{num_cpu}', '--outSAMtype', 'BAM', 'SortedByCoordinate', '--clipAdapterType', 'CellRanger4', '--outFilterScoreMin', '30']

        white_list = 'None'
        if '~{preset_whitelist}' != '':
            white_list = '~{preset_whitelist}'
        else:
            if '~{chemistry}' == 'custom' and '~{CBwhitelist}' != '':
                white_list = '~{CBwhitelist}'

        if '~{chemistry}' in ['tenX_v2', 'tenX_v3']:
            call_args.extend(['--soloCBmatchWLtype', '1MM_multi_Nbase_pseudocounts', '--soloUMIfiltering', 'MultiGeneUMI_CR', '--soloUMIdedup', '1MM_CR'])
            if '~{chemistry}' == 'tenX_v3':
                call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '12'])
            elif '~{chemistry}' == 'tenX_v2':
                call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '1', '--soloCBlen', '16', '--soloUMIstart', '17', '--soloUMIlen', '10'])
        elif '~{chemistry}' in ['SeqWell', 'DropSeq']:
            call_args.extend(['--soloCBwhitelist', 'None', '--soloCBstart', '1', '--soloCBlen', '12', '--soloUMIstart', '13', '--soloUMIlen', '8'])
        elif '~{chemistry}' == 'custom':
            call_args.extend(['--soloCBwhitelist', white_list, '--soloCBstart', '~{CBstart}', '--soloCBlen', '~{CBlen}', '--soloUMIstart', '~{UMIstart}', '--soloUMIlen', '~{UMIlen}'])

        if file_ext == '.fastq.gz':
            call_args.extend(['--readFilesCommand', 'zcat'])

        def set_up_readfiles(prefix, n_files, f_ext):
            return ','.join([prefix+str(n)+f_ext for n in range(n_files)])

        call_args.extend(['--readFilesIn', set_up_readfiles('R2_', len(r2_list), file_ext), set_up_readfiles('R1_', len(r1_list), file_ext)])
        call_args.extend(['--outFileNamePrefix', 'result/~{sample_id}_'])

        if '~{BarcodeReadLength}' != '':
            call_args.extend(['--soloBarcodeReadLength', '~{BarcodeReadLength}'])
        if '~{BarcodeMate}' != '':
            call_args.extend(['--soloBarcodeMate', '~{BarcodeMate}'])
        if '~{solo_type}' == 'CB_UMI_Complex':
            if '~{CBposition}' != '':
                call_args.extend(['--soloCBposition'] + '~{CBposition}'.split(','))
            if '~{UMIposition}' != '':
                call_args.extend(['--soloUMIposition', '~{UMIposition}'])
        if '~{adapterSequence}' != '':
            call_args.extend(['--soloAdapterSequence', '~{adapterSequence}'])
        if '~{adapterMismatchesNmax}' != '':
            call_args.extend(['--soloAdapterMismatchesNmax', '~{adapterMismatchesNmax}'])
        if '~{CBmatchWLtype}' != '':
            call_args.extend(['--soloCBmatchWLtype', '~{CBmatchWLtype}'])
        if '~{inputSAMattrBarcodeSeq}' != '':
            call_args.extend(['--soloInputSAMattrBarcodeSeq'] + '~{inputSAMattrBarcodeSeq}'.split(','))
        if '~{inputSAMattrBarcodeQual}' != '':
            call_args.extend(['--soloInputSAMattrBarcodeQual'] + '~{inputSAMattrBarcodeQual}'.split(','))
        if '~{strand}' != '':
            call_args.extend(['--soloStrand', '~{strand}'])
        if '~{features}' != '':
            feature_list = '~{features}'.split(',')
            if ('Velocyto' in feature_list) and ('Gene' not in feature_list):
                feature_list.append('Gene')
            call_args.extend(['--soloFeatures'] + '~{features}'.split(','))
        if '~{multiMappers}' != '':
            call_args.extend(['--soloMultiMappers'] + '~{multiMappers}'.split(','))
        if '~{UMIdedup}' != '':
            call_args.extend(['--soloUMIdedup'] + '~{UMIdedup}'.split(','))
        if '~{UMIfiltering}' != '':
            call_args.extend(['--soloUMIfiltering'] + '~{UMIfiltering}'.split(','))
        if '~{outFileNames}' != '':
            call_args.extend(['--soloOutFileNames'] + '~{outFileNames}'.split(','))
        if '~{cellFilter}' != '':
            call_args.extend(['--soloCellFilter'] + '~{cellFilter}'.split(','))
        if '~{outFormatFeaturesGeneField3}' != '':
            call_args.extend(['--soloOutFormatFeaturesGeneField3'] + '~{outFormatFeaturesGeneField3}'.split(','))

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} result "~{output_directory}/~{sample_id}"
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
