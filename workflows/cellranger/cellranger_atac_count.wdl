version 1.0

workflow cellranger_atac_count {
    input {
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # cellRanger-atac output directory, gs url
        String output_directory

        # Keywords or a URL to a tar.gz file
        String genome

        # Index TSV file
        File acronym_file

        # Chemistry to use. If auto full analysis will be made if ARC-v1, just the individual ATAC library from a
        # multiome assay will be used
        String chemistry = "auto"

        # Force pipeline to use this number of cells, bypassing the cell detection algorithm
        Int? force_cells
        # Choose the algorithm for dimensionality reduction prior to clustering and tsne: 'lsa' (default), 'plsa', or 'pca'.
        String? dim_reduce
        # A BED file to override peak caller
        File? peaks

        # cellranger-atac version
        String cellranger_atac_version
        # Which docker registry to use
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per cellranger job
        Int num_cpu = 64
        # Memory string, e.g. 57.6G
        String memory = "57.6G"
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

    call run_cellranger_atac_count {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            output_directory = sub(output_directory, "/+$", ""),
            genome_file = genome_file,
            force_cells = force_cells,
            dim_reduce = dim_reduce,
            peaks = peaks,
            chemistry = chemistry,
            cellranger_atac_version = cellranger_atac_version,
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
        String output_count_directory = run_cellranger_atac_count.output_count_directory
        String output_metrics_summary = run_cellranger_atac_count.output_metrics_summary
        String output_web_summary = run_cellranger_atac_count.output_web_summary
        File monitoringLog = run_cellranger_atac_count.monitoringLog
    }
}


task run_cellranger_atac_count {
    input {
        String sample_id
        String input_fastqs_directories
        String output_directory
        File genome_file
        Int? force_cells
        String? dim_reduce
        File? peaks
        String chemistry
        String cellranger_atac_version
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

        fastqs = []
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            target = '~{sample_id}_' + str(i)
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)
            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}' + '_S*_L*_*_001.fastq.gz' , target]
                print(' '.join(call_args))
                check_call(call_args)
            fastqs.append(target)

        call_args = ['cellranger-atac', 'count', '--id=results', '--reference=genome_dir', '--fastqs=' + ','.join(fastqs), '--sample=~{sample_id}', '--jobmode=local']
        if '~{force_cells}' != '':
            call_args.append('--force-cells=~{force_cells}')
        if '~{dim_reduce}' != '':
            call_args.append('--dim-reduce=~{dim_reduce}')
        if '~{peaks}' != '':
            assert version.parse('~{cellranger_atac_version}') >= version.parse('2.0.0')
            call_args.append('--peaks=~{peaks}')
        if '~{chemistry}'.casefold() == 'ARC-v1'.casefold():
            call_args.append('--chemistry=ARC-v1')
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m results/outs "~{output_directory}/~{sample_id}"
    }

    output {
        String output_count_directory = "~{output_directory}/~{sample_id}"
        String output_metrics_summary = "~{output_directory}/~{sample_id}/summary.csv"
        String output_web_summary = "~{output_directory}/~{sample_id}/web_summary.html"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger-atac:~{cellranger_atac_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
