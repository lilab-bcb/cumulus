version 1.0

workflow cellranger_count {
    input {
        # Sample ID
        String sample_id
        # A comma-separated list of input sample names
        String? input_samples
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # A comma-separated list of input data types
        String? input_data_types
        # A comma-separated list of input feature barcode files
        String? input_fbf
        # CellRanger output directory, gs url
        String output_directory

        # GRCh38, hg19, mm10, GRCh38_and_mm10, GRCh38_premrna, mm10_premrna, GRCh38_premrna_and_mm10_premrna or a URL to a tar.gz file
        String genome
        # Index TSV file
        File acronym_file

        # Target panel CSV for targeted gene expression analysis
        File? target_panel

        # chemistry of the channel
        String chemistry = "auto"
        # Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells.
        Int? force_cells
        # Expected number of recovered cells. Mutually exclusive with force_cells
        Int? expect_cells
        # If count reads mapping to intronic regions
        Boolean include_introns = true
        # If generate bam outputs
        Boolean no_bam = false
        # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
        Boolean secondary = false

        # cellranger version
        String cellranger_version
        # Which docker registry to use: cumulusprod (default) or quay.io/cumulus
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per cellranger job
        Int num_cpu = 32
        # Memory string, e.g. 120G
        String memory = "120G"
        # Disk space in GB
        Int disk_space = 500
        # Number of preemptible tries
        Int preemptible = 2
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend = "gcp"
    }

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    # If reference is a url
    Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

    File genome_file = (if is_url then genome else acronym2gsurl[genome])

    call run_cellranger_count {
        input:
            sample_id = sample_id,
            input_samples = input_samples,
            input_fastqs_directories = input_fastqs_directories,
            input_data_types = input_data_types,
            input_fbf = input_fbf,
            output_directory = output_directory,
            genome_file = genome_file,
            target_panel = target_panel,
            chemistry = chemistry,
            force_cells = force_cells,
            expect_cells = expect_cells,
            include_introns = include_introns,
            no_bam = no_bam,
            secondary = secondary,
            cellranger_version = cellranger_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    output {
        String output_count_directory = run_cellranger_count.output_count_directory
        String output_metrics_summary = run_cellranger_count.output_metrics_summary
        String output_web_summary = run_cellranger_count.output_web_summary
        File monitoringLog = run_cellranger_count.monitoringLog
    }
}

task run_cellranger_count {
    input {
        String sample_id
        String? input_samples
        String input_fastqs_directories
        String? input_data_types
        String? input_fbf
        String output_directory
        File genome_file
        File? target_panel
        String chemistry
        Int? force_cells
        Int? expect_cells
        Boolean include_introns
        Boolean no_bam
        Boolean secondary
        String cellranger_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &
        mkdir -p genome_dir
        tar xf ~{genome_file} -C genome_dir --strip-components 1

        python <<CODE
        import re
        import os
        import sys
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT
        from packaging import version

        samples = data_types = fbfs = None
        fastqs_dirs = []

        if '~{input_samples}' != '':
            assert version.parse('~{cellranger_version}') >= version.parse('3.0.2')
            samples = '~{input_samples}'.split(',')
            data_types = '~{input_data_types}'.split(',')
            fbfs = '~{input_fbf}'.split(',')

            target_panel = set()
            feature_file = set()
            for dtype, fbf in zip(data_types, fbfs):
                if dtype == 'rna':
                    target_panel.add(fbf)
                else:
                    feature_file.add(fbf)

            def _locate_file(file_set, keyword):
                if len(file_set) > 1:
                    print("Detected multiple " + keyword + " files!", file = sys.stderr)
                    sys.exit(1)
                if len(file_set) == 0 or list(file_set)[0] == 'null':
                    return ''
                file_loc = list(file_set)[0]
                call_args = ['strato', 'cp', '-m', file_loc, '.']
                print(' '.join(call_args))
                check_call(call_args)
                return os.path.abspath(os.path.basename(file_loc))

            target_panel = _locate_file(target_panel, 'target panel')
            feature_file = _locate_file(feature_file, 'feature reference')
            assert feature_file != ''

            with open('libraries.csv', 'w') as fout:
                fout.write('fastqs,sample,library_type\n')
                for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
                    directory = re.sub('/+$', '', directory) # remove trailing slashes
                    target = samples[i] + "_" + str(i)
                    try:
                        call_args = ['strato', 'exists', directory + '/' + samples[i] + '/']
                        print(' '.join(call_args))
                        check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                        call_args = ['strato', 'sync', '-m', directory + '/' + samples[i], target]
                        print(' '.join(call_args))
                        check_call(call_args)
                    except CalledProcessError:
                        if not os.path.exists(target):
                            os.mkdir(target)
                        call_args = ['strato', 'cp', '-m', directory + '/' + samples[i] + '_S*_L*_*_001.fastq.gz' , target]
                        print(' '.join(call_args))
                        check_call(call_args)
                    feature_type = ''
                    if data_types[i] == 'rna':
                        feature_type = 'Gene Expression'
                    elif data_types[i] == 'crispr':
                        feature_type = 'CRISPR Guide Capture'
                    elif data_types[i] == 'citeseq':
                        feature_type = 'Antibody Capture'
                    elif data_types[i] == 'hashing':
                        feature_type = 'Custom'
                    if feature_type == '':
                        print("Do not expect " + data_types[i] + " in a cellranger count (Feature Barcode) run!", file = sys.stderr)
                        sys.exit(1)
                    fout.write(os.path.abspath(target) + ',' + samples[i] + ',' + feature_type + '\n')
        else:
            for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
                directory = re.sub('/+$', '', directory) # remove trailing slashes
                target = '~{sample_id}_' + str(i)
                try:
                    call_args = ['strato', 'exists', directory + '/~{sample_id}/']
                    print(' '.join(call_args))
                    check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                    call_args = ['strato', 'sync', '-m', directory + '/~{sample_id}', target]
                    print(' '.join(call_args))
                    check_call(call_args)
                except CalledProcessError:
                    if not os.path.exists(target):
                        os.mkdir(target)
                    call_args = ['strato', 'cp', '-m', directory + '/~{sample_id}' + '_S*_L*_*_001.fastq.gz' , target]
                    print(' '.join(call_args))
                    check_call(call_args)
                fastqs_dirs.append(target)

        mem_size = re.findall(r"\d+", "~{memory}")[0]
        call_args = ['cellranger', 'count', '--id=results', '--transcriptome=genome_dir', '--chemistry=~{chemistry}', '--jobmode=local', '--localcores=~{num_cpu}', '--localmem='+mem_size]

        if samples is None: # not Feature Barcode
            call_args.extend(['--sample=~{sample_id}', '--fastqs=' + ','.join(fastqs_dirs)])
            if ('~{target_panel}' != '') and (os.path.basename('~{target_panel}') != 'null'):
                assert version.parse('~{cellranger_version}') >= version.parse('4.0.0')
                call_args.append('--target-panel=~{target_panel}')
        else:
            call_args.extend(['--libraries=libraries.csv', '--feature-ref=' + feature_file])
            if target_panel != '':
                assert version.parse('~{cellranger_version}') >= version.parse('4.0.0')
                call_args.append('--target-panel=' + target_panel)

        if '~{force_cells}' != '':
            call_args.append('--force-cells=~{force_cells}')
        if '~{expect_cells}' != '':
            call_args.append('--expect-cells=~{expect_cells}')
        if '~{include_introns}' == 'true':
            assert version.parse('~{cellranger_version}') >= version.parse('5.0.0')
            if version.parse('~{cellranger_version}') >= version.parse('7.0.0'):
                call_args.extend(['--include-introns', '~{include_introns}'])
            else:
                call_args.append('--include-introns')
        elif version.parse('~{cellranger_version}') >= version.parse('7.0.0'):
            call_args.extend(['--include-introns', '~{include_introns}'])

        # For generating BAM output
        if version.parse('~{cellranger_version}') >= version.parse('8.0.0'):
            if '~{no_bam}' == 'false':
                call_args.append('--create-bam=true')
            else:
                call_args.append('--create-bam=false')
        else:
            if '~{no_bam}' == 'true':
                assert version.parse('~{cellranger_version}') >= version.parse('5.0.0')
                call_args.append('--no-bam')

        if '~{secondary}' != 'true':
            call_args.append('--nosecondary')
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync -m results/outs "~{output_directory}/~{sample_id}"
    }

    output {
        String output_count_directory = "~{output_directory}/~{sample_id}"
        String output_metrics_summary = "~{output_directory}/~{sample_id}/metrics_summary.csv"
        String output_web_summary = "~{output_directory}/~{sample_id}/web_summary.html"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger:~{cellranger_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
