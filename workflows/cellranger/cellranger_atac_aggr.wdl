version 1.0

workflow cellranger_atac_aggr {
    input {
        # Aggregate ID
        String aggr_id
        # A comma-separated list of input atac count result directories (gs urls), note that each directory should contain fragments.tsv.gz and singlecell.csv
        String input_counts_directories
        # CellRanger-atac output directory, gs url
        String output_directory
        # Index TSV file
        File acronym_file
        # Backend
        String backend = "gcp"

        # Keywords or a URL to a tar.gz file
        String genome

        # Sample normalization MODE: none (default), depth, signal
        String normalize = "none"
        # Perform secondary analysis (dimensionality reduction, clustering and visualization). Default: false
        Boolean secondary = false
        # Chose the algorithm for dimensionality reduction prior to clustering and tsne: 'lsa' (default), 'plsa', or 'pca'.
        String dim_reduce = "lsa"
        # A BED file to override peak caller
        File? peaks

        # 2.0.0, 1.2.0, 1.1.0
        String cellranger_atac_version = "2.0.0"
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

        # Which docker registry to use: cumulusprod (default) or quay.io/cumulus
        String docker_registry = "cumulusprod"
    }

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    # If reference is a url
    Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

    File genome_file = (if is_url then genome else acronym2gsurl[genome])

    call run_cellranger_atac_aggr {
        input:
            aggr_id = aggr_id,
            input_counts_directories = input_counts_directories,
            output_directory = sub(output_directory, "/+$", ""),
            genome_file = genome_file,
            normalize = normalize,
            secondary = secondary,
            dim_reduce = dim_reduce,
            peaks = peaks,
            cellranger_atac_version = cellranger_atac_version,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            docker_registry = docker_registry,
            backend = backend
    }

    output {
        String output_aggr_directory = run_cellranger_atac_aggr.output_aggr_directory
        String output_metrics_summary = run_cellranger_atac_aggr.output_metrics_summary
        String output_web_summary = run_cellranger_atac_aggr.output_web_summary
        File monitoringLog = run_cellranger_atac_aggr.monitoringLog
    }
}

task run_cellranger_atac_aggr {
    input {
        String aggr_id
        String input_counts_directories
        String output_directory
        File genome_file
        String normalize
        Boolean secondary
        String dim_reduce
        File? peaks
        String cellranger_atac_version
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        String docker_registry
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
        from subprocess import check_call
        from packaging import version

        counts = []
        with open('aggr.csv', 'w') as fout:
            fout.write('library_id,fragments,cells\n')
            libs_seen = set()
            current_dir = os.getcwd()
            for i, directory in enumerate('~{input_counts_directories}'.split(',')):
                directory = re.sub('/+$', '', directory) # remove trailing slashes

                library_id = os.path.basename(directory)
                if library_id in libs_seen:
                    raise Exception("Found duplicated library id " + library_id + "!")
                libs_seen.add(library_id)

                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', '-r', directory, current_dir]
                print(' '.join(call_args))
                check_call(call_args)
                counts.append(library_id)
                fout.write(library_id + "," + current_dir + '/' + library_id + "/fragments.tsv.gz," + current_dir + '/' + library_id + "/singlecell.csv\n")

        call_args = ['cellranger-atac', 'aggr', '--id=results', '--reference=genome_dir', '--csv=aggr.csv', '--normalize=~{normalize}', '--jobmode=local']
        if '~{secondary}' != 'true':
            call_args.append('--nosecondary')
        else:
            call_args.append('--dim-reduce=~{dim_reduce}')
        if '~{peaks}' != '':
            assert version.parse('~{cellranger_atac_version}') >= version.parse('2.0.0')
            call_args.append('--peaks=~{peaks}')
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m results/outs "~{output_directory}/~{aggr_id}"
    }

    output {
        String output_aggr_directory = "~{output_directory}/~{aggr_id}"
        String output_metrics_summary = "~{output_directory}/~{aggr_id}/summary.csv"
        String output_web_summary = "~{output_directory}/~{aggr_id}/web_summary.html"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger-atac:~{cellranger_atac_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: "~{num_cpu}"
        preemptible: "~{preemptible}"
    }
}
