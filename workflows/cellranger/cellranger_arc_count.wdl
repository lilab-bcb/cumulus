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

        # If generate bam outputs
        Boolean no_bam = false

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
            no_bam = no_bam,
            cellranger_arc_version = cellranger_arc_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
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
        Boolean no_bam
        String cellranger_arc_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
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

        samples = '~{input_samples}'.split(',')
        data_types = '~{input_data_types}'.split(',')
        with open('libraries.csv', 'w') as fout:
            fout.write('fastqs,sample,library_type\n')
            for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
                directory = re.sub('/+$', '', directory) # remove trailing slashes
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', '-r', directory + '/' + samples[i], '.']
                print(' '.join(call_args))
                check_call(call_args)
                fastqs = samples[i] + '_' + str(i)
                call_args = ['mv', samples[i], fastqs]
                print(' '.join(call_args))
                check_call(call_args)
                fout.write(os.path.abspath(fastqs) + ',' + samples[i] + ',' + ('Gene Expression' if data_types[i] == 'rna' else 'Chromatin Accessibility') + '\n')

        call_args = ['cellranger-arc', 'count', '--id=results', '--libraries=libraries.csv', '--reference=genome_dir', '--jobmode=local']
        if '~{no_bam}' is 'true':
            call_args.append('--no-bam')
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m results/outs "~{output_directory}"/~{link_id}
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
        cpu: "~{num_cpu}"
        preemptible: "~{preemptible}"
    }
}
