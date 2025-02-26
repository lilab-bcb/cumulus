version 1.0

workflow cellranger_multi {
    input {
        # Link ID
        String link_id
        # A comma-separated list of input sample names
        String input_samples
        # A comma-separated list of input FASTQs directories
        String input_fastqs_directories
        # A comma-separated list of input data types
        String input_data_types
        # A comma-separated list of input auxiliary files
        String input_aux
        # CellRanger multi output directory
        String output_directory

        # Keywords or a URL to a tar.gz file
        String genome
        String vdj_ref = "null"
        # Flex probe set CSV file
        File probe_set_file

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
    Boolean is_genome_uri = sub(genome, "^.+\\.(tgz|gz)$", "URI") == "URI"
    File genome_file = (if is_genome_uri then genome else acronym2uri[genome])

    # If vdj reference is a URI
    Boolean is_vdj_ref_uri = sub(vdj_ref, "^.+\\.(tgz|gz)$", "URI") == "URI"
    File vdj_ref_file = (if vdj_ref != "null" then (if is_vdj_ref_uri then vdj_ref else acronym2uri[vdj_ref]) else acronym2uri["null_file"])

    call run_cellranger_multi {
        input:
            link_id = link_id,
            input_samples = input_samples,
            input_fastqs_directories = input_fastqs_directories,
            input_data_types = input_data_types,
            input_aux = input_aux,
            output_directory = output_directory,
            genome_file = genome_file,
            vdj_ref_file = vdj_ref_file,
            probe_set_file = probe_set_file,
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
        String input_aux
        String output_directory
        File genome_file
        File vdj_ref_file
        File probe_set_file
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
        auxs = '~{input_aux}'.split(',')

        rna_file = set()
        vdj_file = set()
        cmo_file = set()
        flex_file = set()
        feature_file = set()
        has_hto = False
        has_cmo = False
        has_flex = False
        for dtype, aux in zip(data_types, auxs):
            if dtype == 'rna':  # OCM, HTO, and CMO cases
                rna_file.add(aux)
            elif dtype == 'cmo':
                cmo_file.add(aux)
                has_cmo = True
            elif dtype == 'frp':
                flex_file.add(aux)
                has_flex = True
            else:  # hashing, citeseq, adt
                feature_file.add(aux)
                if dtype in ['hashing', 'adt']:
                    has_hto = True

        def is_null_file(filename):
            return filename == "" or os.path.basename(filename) == "null"

        def _locate_file(file_set, keyword):
            if len(file_set) > 1:
                print("Detected multiple " + keyword + " files!", file = sys.stderr)
                sys.exit(1)
            if len(file_set) == 0 or is_null_file(list(file_set)[0]) or list(file_set)[0] == 'null':
                return ''
            file_loc = list(file_set)[0]
            call_args = ['strato', 'cp', file_loc, '.']
            print(' '.join(call_args))
            check_call(call_args)
            return os.path.abspath(os.path.basename(file_loc))

        if has_cmo:
            rna_file = _locate_file(rna_file, 'CMO sample')
        elif has_hto:
            rna_file = _locate_file(rna_file, 'HTO sample')
        else:
            rna_file = _locate_file(rna_file, 'OCM sample')
        vdj_file = _locate_file(vdj_file, 'VDJ inner-enrichment-primers')
        cmo_file = _locate_file(cmo_file, 'CMO cmo-set')
        feature_file = _locate_file(feature_file, 'feature reference')
        flex_file = _locate_file(flex_file, 'Flex sample')

        with open('multi.csv', 'w') as fout:
            #############################
            # [gene-expression] section #
            #############################
            fout.write('[gene-expression]\n')
            fout.write('reference,' + os.path.abspath('genome_dir') + '\n')

            if is_null_file('~{probe_set_file}'):  # GEX case
                if '~{include_introns}' == 'false':
                    fout.write('include-introns,false\n')
            else:   # Flex case
                fout.write('probe-set,~{probe_set_file}\n')

            # For CellPlex
            if not is_null_file(cmo_file):
                fout.write('cmo-set,' + cmo_file + '\n')

            # Other options
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

            #################
            # [vdj] section #
            #################
            if not is_null_file('~{vdj_ref_file}'):
                fout.write('\n[vdj]\nreference,~{vdj_ref_file}\n')
                if not is_null_file(vdj_file):
                    fout.write('inner-enrichment-primers,' + vdj_file + '\n')

            #####################
            # [feature] section #
            #####################
            if feature_file != '':
                fout.write('\n[feature]\nreference,' + feature_file + '\n')

            #######################
            # [libraries] section #
            #######################
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
                elif data_types[i] == 'vdj':
                    feature_type = 'VDJ'
                elif data_types[i] == 'vdj_t':
                    feature_type = 'VDJ-T'
                elif data_types[i] == 'vdj_b':
                    feature_type = 'VDJ-B'
                elif data_types[i] == 'vdj_t_gd':
                    feature_type = 'VDJ-T-GD'
                elif data_types[i] == 'crispr':
                    feature_type = 'CRISPR Guide Capture'
                elif data_types[i] in ['citeseq', 'hashing', 'adt']:
                    feature_type = 'Antibody Capture'
                elif data_types[i] == 'cmo':
                    feature_type = 'Multiplexing Capture'
                if feature_type == '':
                    print("Do not expect " + data_types[i] + " in a cellranger multi run!", file = sys.stderr)
                    sys.exit(1)
                fout.write(samples[i] + ',' + os.path.abspath(target) + ',' +  feature_type + '\n')

            #####################
            # [samples] section #
            #####################
            def write_csv_wise(full_columns, fin, fout):
                lines = fin.readlines()
                fout.write("\n[samples]\n")
                if 'sample_id' not in lines[0].strip().split(','):
                    # Add header
                    columns = full_columns[0:len(lines[0].split(','))]
                    fout.write(",".join(columns) + "\n")
                for l in lines:
                    fout.write(l)

            has_ocm = False
            if rna_file != '':
                with open(rna_file, 'r') as fin:
                    if has_hto:  # HTO
                        write_csv_wise(['sample_id', 'hashtag_ids', 'description'], fin, fout)
                    elif has_cmo:  # CMO
                        write_csv_wise(['sample_id', 'cmo_ids', 'description'], fin, fout)
                    else:  # OCM
                        write_csv_wise(['sample_id', 'ocm_barcode_ids', 'description'], fin, fout)
                        has_ocm = True

            if has_flex:
                if flex_file != '':  # Multiplex Flex case
                    with open(flex_file, 'r') as fin:
                        write_csv_wise(['sample_id', 'probe_barcode_ids', 'description'], fin, fout)

            if (not has_ocm) and (not has_hto) and (not has_cmo) and (not has_flex):
                raise Exception("Cannot locate OCM, HTO, CMO, or Flex sample file!")

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
