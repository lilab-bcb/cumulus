version 1.0

workflow cellranger_vdj {
    input {
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # CellRanger output directory, gs url
        String output_directory

        # GRCh38_vdj, GRCm38_vdj or a URL to a tar.gz file
        String genome

        # Do not align reads to reference V(D)J sequences before de novo assembly. Default: false
        Boolean denovo = false
        
        # Force the analysis to be carried out for a particular chain type. The accepted values are:
        #   "auto" for autodetection based on TR vs IG representation (default),
        #   "TR" for T cell receptors,
        #   "IG" for B cell receptors,
        # Use this in rare cases when automatic chain detection fails.
        String chain = "auto"

        # cellranger version
        String cellranger_version
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

        # Which docker registry to use: cumulusprod (default) or quay.io/cumulus
        String docker_registry
    }

    File acronym_file = "gs://regev-lab/resources/cellranger/index.tsv"
    # File acronym_file = "index.tsv"
    Map[String, String] acronym2gsurl = read_map(acronym_file)
    # If reference is a url
    Boolean is_url = sub(genome, "^.+\\.(tgz|gz)$", "URL") == "URL"

    File genome_file = (if is_url then genome else acronym2gsurl[genome])

    call run_cellranger_vdj {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            output_directory = sub(output_directory, "/+$", ""),
            genome_file = genome_file,
            denovo = denovo,
            chain = chain,
            cellranger_version = cellranger_version,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            docker_registry = docker_registry
    }

    output {
        String output_vdj_directory = run_cellranger_vdj.output_vdj_directory
        String output_metrics_summary = run_cellranger_vdj.output_metrics_summary
        String output_web_summary = run_cellranger_vdj.output_web_summary
        File monitoringLog = run_cellranger_vdj.monitoringLog
    }
}

task run_cellranger_vdj {
    input {
        String sample_id
        String input_fastqs_directories
        String output_directory
        File genome_file
        Boolean denovo
        String chain
        String cellranger_version
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        String docker_registry
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        mkdir -p ref_dir
        tar xf ~{genome_file} -C ref_dir --strip-components 1

        python <<CODE
        import re
        from subprocess import check_call

        fastqs = []
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes 
            call_args = ['gsutil', '-q', '-m', 'cp', '-r', directory + '/~{sample_id}', '.']
            # call_args = ['cp', '-r', directory + '/~{sample_id}', '.']
            print(' '.join(call_args))
            check_call(call_args)
            call_args = ['mv', '~{sample_id}', '~{sample_id}_' + str(i)]
            print(' '.join(call_args))
            check_call(call_args)
            fastqs.append('~{sample_id}_' + str(i))

        call_args = ['cellranger', 'vdj', '--chain=~{chain}', '--id=results', '--reference=ref_dir', '--fastqs=' + ','.join(fastqs), '--sample=~{sample_id}', '--jobmode=local']
        if '~{denovo}' is not 'false':
            call_args.append('--denovo')
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -q -m rsync -d -r results/outs ~{output_directory}/~{sample_id}
        # cp -r results/outs ~{output_directory}/~{sample_id}
    }

    output {
        String output_vdj_directory = "~{output_directory}/~{sample_id}"
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
        cpu: "~{num_cpu}"
        preemptible: "~{preemptible}"
    }
}
