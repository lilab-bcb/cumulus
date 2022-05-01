version 1.0

workflow cellranger_multi {
    input {
        # Link ID
        String link_id
        # A comma-separated list of input sample names
        String input_samples
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # A comma-separated list of input data types
        String input_data_types
        # A comma-separated list of input feature barcode files
        String input_fbf
        # CellRanger multi output directory
        String output_directory

        # Keywords or a URL to a tar.gz file
        String genome

        # CMO set CSV file, delaring CMO constructs and associated barcodes
        File? cmo_set

        # Index TSV file
        File acronym_file

        # Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells.
        Int? force_cells
        # Expected number of recovered cells. Mutually exclusive with force_cells
        Int? expect_cells
        # If count reads mapping to intronic regions
        Boolean include_introns = false
        # If generate bam outputs
        Boolean no_bam = false
        # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
        Boolean secondary = false

        # CellRanger version, must be >= 6.0.0
        String cellranger_version
        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with Cromwell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per cellranger job
        Int num_cpu = 32
        # Memory string, e.g. 45G
        String memory = "120G"
        # Disk space in GB
        Int disk_space = 500
        # Number of preemptible tries
        Int preemptible = 2
        # Max number of retries for AWS instance
        Int awsMaxRetries = 5
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend = "gcp"
    }

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    # If reference is a url
    Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

    File genome_file = (if is_url then genome else acronym2gsurl[genome])

    call run_cellranger_multi {
        input:
            link_id = link_id,
            input_samples = input_samples,
            input_fastqs_directories = input_fastqs_directories,
            input_data_types = input_data_types,
            input_fbf = input_fbf,
            output_directory = output_directory,
            genome_file = genome_file,
            cmo_set = cmo_set,
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
            awsMaxRetries = awsMaxRetries,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    output {
        String output_multi_directory = run_cellranger_multi.output_multi_directory
    }
}

task run_cellranger_multi {
    input {
        String link_id
        String input_samples
        String input_fastqs_directories
        String input_data_types
        String input_fbf
        String output_directory
        File genome_file
        File? cmo_set
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
        Int awsMaxRetries
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        mkdir -p genome_dir
        tar xf ~{genome_file} -C genome_dir --strip-components 1

        python <<CODE
        import re
        import os
        import sys
        from subprocess import check_call, CalledProcessError, STDOUT, DEVNULL

        samples = '~{input_samples}'.split(',')
        data_types = '~{input_data_types}'.split(',')
        fbfs = '~{input_fbf}'.split(',')

        target_panel = set()
        cmo_file = set()
        feature_file = set()
        for dtype, fbf in zip(data_types, fbfs):
            if dtype == 'rna':
                target_panel.add(fbf)
            elif dtype == 'cmo':
                cmo_file.add(fbf)
            else:
                feature_file.add(fbf)

        def _locate_file(file_set, keyword):
            if len(file_set) > 1:
                print("Detected multiple " + keyword + " files!", file = sys.stderr)
                sys.exit(1)
            if len(file_set) == 0 or list(file_set)[0] == 'null':
                return ''
            file_loc = list(file_set)[0]
            call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', file_loc, '.']
            print(' '.join(call_args))
            check_call(call_args)
            return os.path.abspath(os.path.basename(file_loc))

        target_panel = _locate_file(target_panel, 'target panel')
        cmo_file = _locate_file(cmo_file, 'CMO sample')
        feature_file = _locate_file(feature_file, 'feature reference')

        with open('multi.csv', 'w') as fout:
            fout.write('[gene-expression]\n')
            fout.write('reference,' + os.path.abspath('genome_dir') + '\n')
            if '~{cmo_set}' != '':
                fout.write('cmo-set,~{cmo_set}\n')
            if target_panel != '':
                fout.write('target-panel,' + target_panel + '\n')
            if '~{force_cells}' != '':
                fout.write('force-cells,~{force_cells}\n')
            if '~{expect_cells}' != '':
                fout.write('expect-cells,~{expect_cells}\n')
            if '~{include_introns}' == 'true':
                fout.write('include-introns,true\n')
            if '~{secondary}' == 'false':
                fout.write('no-secondary,true\n')
            if '~{no_bam}' == 'true':
                fout.write('no-bam,true\n')

            if feature_file != '':
                fout.write('\n[feature]\nreference,' + feature_file + '\n')

            fout.write('\n[libraries]\nfastq_id,fastqs,feature_types\n')
            for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
                directory = re.sub('/+$', '', directory) # remove trailing slashes
                target = samples[i] + "_" + str(i)
                try:
                    call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/' + samples[i] + '/']
                    print(' '.join(call_args))
                    check_call(call_args, stderr=STDOUT, stdout=DEVNULL)
                    call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', '-r', directory + '/' + samples[i], target]
                    print(' '.join(call_args))
                    check_call(call_args)
                except CalledProcessError:
                    if not os.path.exists(target):
                        os.mkdir(target)
                    call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/' + samples[i] + '_S*_L*_*_001.fastq.gz' , target]
                    print(' '.join(call_args))
                    check_call(call_args)
                feature_type = ''
                if data_types[i] == 'rna':
                    feature_type = 'Gene Expression'
                elif data_types[i] == 'crispr':
                    feature_type = 'CRISPR Guide Capture'
                elif data_types[i] == 'citeseq':
                    feature_type = 'Antibody Capture'
                elif data_types[i] == 'cmo':
                    feature_type = 'Multiplexing Capture'
                if feature_type == '':
                    print("Do not expect " + data_types[i] + " in a cellranger multi run!", file = sys.stderr)
                    sys.exit(1)
                fout.write(samples[i] + ',' + os.path.abspath(target) + ',' +  feature_type + '\n')

            if cmo_file == '':
                print("Cannot locate a CMO sample file!", file = sys.stderr)
                sys.exit(1)
            fout.write('\n[samples]\nsample_id,cmo_ids\n')
            with open(cmo_file) as fin:
                for line in fin:
                    fout.write(line)

        call_args = ['cellranger', 'multi', '--id=results', '--csv=multi.csv', '--jobmode=local']
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m --ionice results/outs "~{output_directory}"/~{link_id}
    }

    output {
        String output_multi_directory = "~{output_directory}/~{link_id}"
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
        maxRetries: if backend == "aws" then awsMaxRetries else 0
        queueArn: awsQueueArn
    }
}
