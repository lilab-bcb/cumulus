version 1.0

workflow cellranger_arc_count {
    input {
        # Link ID
        String link_id
        # A comma-separated list of input sample names
        String input_samples
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # A comma-separated list of input data types
        String input_data_types
        # CellRanger ARC output directory
        String output_directory

        # Keywords or a URL to a tar.gz file
        String genome

        # Index TSV file
        File acronym_file

        # Disable counting of intronic reads. In this mode, only reads that are exonic and compatible with annotated splice junctions in the reference are counted.
        # Note: using this mode will reduce the UMI counts in the feature-barcode matrix
        Boolean gex_exclude_introns = false
        # Do not generate bam files
        Boolean no_bam = false
        # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
        Boolean secondary = false
        # Cell caller override: define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode.
        # Note: this option must be specified in conjunction with `min-gex-count`.
        # With `--min-atac-count=X` and `--min-gex-count=Y` a barcode is defined as a cell if it contains at least X ATAC counts AND at least Y GEX UMI counts
        Int? min_atac_count
        # Cell caller override: define the minimum number of GEX UMI counts for a cell barcode.
        # Note: this option must be specified in conjunction with `min-atac-count`.
        # With `--min-atac-count=X` and `--min-gex-count=Y` a barcode is defined as a cell if it contains at least X ATAC counts AND at least Y GEX UMI counts
        Int? min_gex_count
        # Override peak caller: specify peaks to use in downstream analyses from supplied 3-column BED file.
        # The supplied peaks file must be sorted by position and not contain overlapping peaks; comment lines beginning with `#` are allowed
        File? peaks

        # CellRanger ARC version
        String cellranger_arc_version
        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with Cromwell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per cellranger job
        Int num_cpu = 64
        # Memory string, e.g. 160G
        String memory = "160G"
        # Disk space in GB
        Int disk_space = 700
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

    call run_cellranger_arc_count {
        input:
            link_id = link_id,
            input_samples = input_samples,
            input_fastqs_directories = input_fastqs_directories,
            input_data_types = input_data_types,
            output_directory = output_directory,
            genome_file = genome_file,
            gex_exclude_introns = gex_exclude_introns,
            no_bam = no_bam,
            secondary = secondary,
            min_atac_count = min_atac_count,
            min_gex_count = min_gex_count,
            peaks = peaks,
            cellranger_arc_version = cellranger_arc_version,
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
        String output_count_directory = run_cellranger_arc_count.output_count_directory
        String output_metrics_summary = run_cellranger_arc_count.output_metrics_summary
        String output_web_summary = run_cellranger_arc_count.output_web_summary
    }
}

task run_cellranger_arc_count {
    input {
        String link_id
        String input_samples
        String input_fastqs_directories
        String input_data_types
        String output_directory
        File genome_file
        Boolean gex_exclude_introns
        Boolean no_bam
        Boolean secondary
        Int? min_atac_count
        Int? min_gex_count
        File? peaks
        String cellranger_arc_version
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
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT
        from packaging import version

        samples = '~{input_samples}'.split(',')
        data_types = '~{input_data_types}'.split(',')
        with open('libraries.csv', 'w') as fout:
            fout.write('fastqs,sample,library_type\n')
            for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
                directory = re.sub('/+$', '', directory) # remove trailing slashes
                target = samples[i] + '_' + str(i)
                try:
                    call_args = ['strato', 'exists', directory + '/' + samples[i] + '/']
                    print(' '.join(call_args))
                    check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                    call_args = ['strato', 'sync', directory + '/' + samples[i], target]
                    print(' '.join(call_args))
                    check_call(call_args)
                except CalledProcessError:
                    if not os.path.exists(target):
                        os.mkdir(target)
                    call_args = ['strato', 'cp', directory + '/' + samples[i] + '_S*_L*_*_001.fastq.gz' , target]
                    print(' '.join(call_args))
                    check_call(call_args)
                fout.write(os.path.abspath(target) + ',' + samples[i] + ',' + ('Gene Expression' if data_types[i] == 'rna' else 'Chromatin Accessibility') + '\n')

        mem_size = re.findall(r"\d+", "~{memory}")[0]
        call_args = ['cellranger-arc', 'count', '--id=results', '--libraries=libraries.csv', '--reference=genome_dir', '--jobmode=local', '--localcores=~{num_cpu}', '--localmem='+mem_size]
        if '~{gex_exclude_introns}' == 'true':
            call_args.append('--gex-exclude-introns')

        # For generating BAM output
        if version.parse("~{cellranger_arc_version}") >= version.parse("2.1.0"):
            if '~{no_bam}' == 'false':
                call_args.append('--create-bam=true')
            else:
                call_args.append('--create-bam=false')
        else:
            if '~{no_bam}' == 'true':
                call-args.append('--no-bam')

        if '~{secondary}' != 'true':
            call_args.append('--nosecondary')

        if '~{min_atac_count}' != '':
            call_args.extend(['--min-atac-count', '~{min_atac_count}'])
        if '~{min_gex_count}' != '':
            call_args.extend(['--min-gex-count', '~{min_gex_count}'])
        if '~{peaks}' != '':
            call_args.extend(['--peaks', '~{peaks}'])
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync results/outs "~{output_directory}/~{link_id}"
    }

    output {
        String output_count_directory = "~{output_directory}/~{link_id}"
        String output_metrics_summary = "~{output_directory}/~{link_id}/summary.csv"
        String output_web_summary = "~{output_directory}/~{link_id}/web_summary.html"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger-arc:~{cellranger_arc_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
