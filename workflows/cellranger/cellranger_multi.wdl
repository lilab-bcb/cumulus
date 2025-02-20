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
        # Keyword to a preset probe set CSV file
        String probe_set = "null"

        # CMO set CSV file, delaring CMO constructs and associated barcodes
        File? cmo_set

        # Index TSV file
        File acronym_file

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

        # CellRanger version
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
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend = "gcp"
    }

    Map[String, String] acronym2uri = read_map(acronym_file)
    # If reference is a URI
    Boolean is_genome_uri = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"
    File genome_file = (if is_genome_uri then genome else acronym2uri[genome])

    # If probe set is specified
    Boolean is_probeset_uri = sub(probe_set, "^.+\\.csv$", "URL") == "URL"
    File probe_set_file = (if probe_set != "null" then (if is_probeset_uri then probe_set else acronym2uri[probe_set]) else acronym2uri["null_file"])

    call run_cellranger_multi {
        input:
            link_id = link_id,
            input_samples = input_samples,
            input_fastqs_directories = input_fastqs_directories,
            input_data_types = input_data_types,
            input_fbf = input_fbf,
            output_directory = output_directory,
            genome_file = genome_file,
            probe_set_file = probe_set_file,
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
        File probe_set_file
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
        from subprocess import check_call, CalledProcessError, STDOUT, DEVNULL
        from packaging import version

        samples = '~{input_samples}'.split(',')
        data_types = '~{input_data_types}'.split(',')
        fbfs = '~{input_fbf}'.split(',')

        target_panel = set()
        cmo_file = set()
        frp_file = set()
        feature_file = set()
        for dtype, fbf in zip(data_types, fbfs):
            if dtype == 'rna':
                target_panel.add(fbf)
            elif dtype == 'cmo':
                cmo_file.add(fbf)
            elif dtype == 'frp':
                frp_file.add(fbf)
            else:
                feature_file.add(fbf)

        def _locate_file(file_set, keyword):
            if len(file_set) > 1:
                print("Detected multiple " + keyword + " files!", file = sys.stderr)
                sys.exit(1)
            if len(file_set) == 0 or list(file_set)[0] == 'null':
                return ''
            file_loc = list(file_set)[0]
            call_args = ['strato', 'cp', file_loc, '.']
            print(' '.join(call_args))
            check_call(call_args)
            return os.path.abspath(os.path.basename(file_loc))

        def is_null_file(filename):
            return filename == "" or os.path.basename(filename) == "null"

        target_panel = _locate_file(target_panel, 'target panel')
        cmo_file = _locate_file(cmo_file, 'CMO sample')
        feature_file = _locate_file(feature_file, 'feature reference')
        frp_file = _locate_file(frp_file, 'FRP sample')

        with open('multi.csv', 'w') as fout:
            fout.write('[gene-expression]\n')
            fout.write('reference,' + os.path.abspath('genome_dir') + '\n')

            if is_null_file('~{probe_set_file}'):
                if '~{include_introns}' == 'true':
                    fout.write('include-introns,true\n')
                elif version.parse('~{cellranger_version}') >= version.parse('7.0.0'):
                    fout.write('include-introns,false\n')
            else:
                if version.parse('~{cellranger_version}') >= version.parse('7.0.0'):
                    fout.write('probe-set,~{probe_set_file}\n')
                else:
                    print("Fixed RNA Profiling only works in Cell Ranger v7.0.0+!", file=sys.stderr)
                    sys.exit(1)
            if '~{cmo_set}' != '':
                fout.write('cmo-set,~{cmo_set}\n')
            if target_panel != '':
                fout.write('target-panel,' + target_panel + '\n')
            if '~{force_cells}' != '':
                fout.write('force-cells,~{force_cells}\n')
            if '~{expect_cells}' != '':
                fout.write('expect-cells,~{expect_cells}\n')
            if '~{secondary}' == 'false':
                fout.write('no-secondary,true\n')
            if version.parse('~{cellranger_version}') >= version.parse('8.0.0'):
                if '~{no_bam}' == 'false':
                    fout.write('create-bam,true\n')
                else:
                    fout.write('create-bam,false\n')
            else:
                if '~{no_bam}' == 'true':
                    fout.write('no-bam,true\n')

            if feature_file != '':
                fout.write('\n[feature]\nreference,' + feature_file + '\n')

            fout.write('\n[libraries]\nfastq_id,fastqs,feature_types\n')
            for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
                directory = re.sub('/+$', '', directory) # remove trailing slashes
                target = samples[i] + "_" + str(i)
                try:
                    call_args = ['strato', 'exists', directory + '/' + samples[i] + '/']
                    print(' '.join(call_args))
                    check_call(call_args, stderr=STDOUT, stdout=DEVNULL)
                    call_args = ['strato', 'cp', '-r', directory + '/' + samples[i], target]
                    print(' '.join(call_args))
                    check_call(call_args)
                except CalledProcessError:
                    if not os.path.exists(target):
                        os.mkdir(target)
                    call_args = ['strato', 'cp', directory + '/' + samples[i] + '_S*_L*_*_001.fastq.gz' , target]
                    print(' '.join(call_args))
                    check_call(call_args)
                feature_type = ''
                if data_types[i] in ['rna', 'frp']:
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

            def write_csv_wise(full_columns, fin, fout):
                lines = fin.readlines()
                columns = full_columns[0:len(lines[0].split(','))]
                fout.write("\n[samples]\n" + ",".join(columns) + "\n")
                for l in lines:
                    fout.write(l)

            has_cmo = False
            if cmo_file != '':
                with open(cmo_file, 'r') as fin:
                    write_csv_wise(['sample_id', 'cmo_ids', 'description'], fin, fout)
                has_cmo = True

            has_frp = False
            if 'frp' in data_types:
                if frp_file != '':
                    if has_cmo:
                        raise Exception("Cannot have both CMO and FRP sample files!")
                    with open(frp_file, 'r') as fin:
                        write_csv_wise(['sample_id', 'probe_barcode_ids', 'description', 'expect_cells', 'force_cells'], fin, fout)
                has_frp = True

            if (not has_cmo) and (not has_frp):
                raise Exception("Cannot locate CMO or FRP sample file!")

        mem_size = re.findall(r"\d+", "~{memory}")[0]
        call_args = ['cellranger', 'multi', '--id=results', '--csv=multi.csv', '--jobmode=local', '--localcores=~{num_cpu}', '--localmem='+mem_size]
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --ionice results/outs "~{output_directory}/~{link_id}"
    }

    output {
        String output_multi_directory = "~{output_directory}/~{link_id}"
        File multi_csv = "multi.csv"
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
