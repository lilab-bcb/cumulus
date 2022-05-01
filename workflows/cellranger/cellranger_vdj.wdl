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

        # Index TSV file
        File acronym_file

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

    call run_cellranger_vdj {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            output_directory = sub(output_directory, "/+$", ""),
            genome_file = genome_file,
            denovo = denovo,
            chain = chain,
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
        mkdir -p ref_dir
        tar xf ~{genome_file} -C ref_dir --strip-components 1

        python <<CODE
        import re
        import os
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT

        fastqs = []
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            target = '~{sample_id}_' + str(i)
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', '-r', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)
            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}' + '_S*_L*_*_001.fastq.gz' , target]
                print(' '.join(call_args))
                check_call(call_args)
            fastqs.append(target)

        call_args = ['cellranger', 'vdj', '--chain=~{chain}', '--id=results', '--reference=ref_dir', '--fastqs=' + ','.join(fastqs), '--sample=~{sample_id}', '--jobmode=local']
        if '~{denovo}' is not 'false':
            call_args.append('--denovo')
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend "~{backend}" -m results/outs "~{output_directory}/~{sample_id}"
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
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
        queueArn: awsQueueArn
    }
}
