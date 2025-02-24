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
        String reference

        # vdj, vdj_t, vdj_b, or vdj_tgd
        String data_type

        # Index TSV file
        File acronym_file

        # Do not align reads to reference V(D)J sequences before de novo assembly. Default: false
        Boolean denovo = false

        # If inner enrichment primers other than those provided in the 10x kits are used, they need to be specified here as a textfile with one primer per line. Disable secondary analysis, e.g. clustering
        # A cloud URI to the text file
        File inner_enrichment_primers

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
    Boolean is_url = sub(reference, "^.+\\.(tgz|gz)$", "URL") == "URL"

    File ref_file = (if is_url then reference else acronym2gsurl[reference])

    call run_cellranger_vdj {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            output_directory = sub(output_directory, "/+$", ""),
            ref_file = ref_file,
            data_type = data_type,
            denovo = denovo,
            inner_enrichment_primers = inner_enrichment_primers,
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
        File ref_file
        Boolean denovo
        String data_type
        File inner_enrichment_primers
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
        mkdir -p ref_dir
        tar xf ~{ref_file} -C ref_dir --strip-components 1

        python <<CODE
        import re
        import os
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT

        def is_null_file(filename):
            return filename == "" or os.path.basename(filename) == "null"

        fastqs = []
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            target = '~{sample_id}_' + str(i)
            try:
                call_args = ['strato', 'exists', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                call_args = ['strato', 'cp', '-r', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)
            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', directory + '/~{sample_id}' + '_S*_L*_*_001.fastq.gz' , target]
                print(' '.join(call_args))
                check_call(call_args)
            fastqs.append(target)

        mem_size = re.findall(r"\d+", "~{memory}")[0]
        call_args = ['cellranger', 'vdj', '--id=results', '--reference=ref_dir', '--fastqs=' + ','.join(fastqs), '--sample=~{sample_id}', '--jobmode=local', '--localcores=~{num_cpu}', '--localmem='+mem_size]
        chain = 'auto'
        if '~{data_type}' == 'vdj_t':
            chain = 'TR'
        elif '~{data_type}' == 'vdj_t_gd':
            chain = 'TR'
            assert not is_null_file('~{inner_enrichment_primers}'), "Sample '~{sample_id}' of vdj_t_gd DataType doesn't have associated inner_enrichment_primers!"
        else:
            chain = 'IG'
        call_args.append('--chain=' + chain)
        if '~{denovo}' != 'false':
            call_args.append('--denovo')
        if not is_null_file('~{inner_enrichment_primers}'):
            call_args.extend(['--inner-enrichment-primers', '~{inner_enrichment_primers}'])
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync results/outs "~{output_directory}/~{sample_id}"
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
        queueArn: awsQueueArn
    }
}
