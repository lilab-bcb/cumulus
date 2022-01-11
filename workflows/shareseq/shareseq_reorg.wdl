version 1.0

workflow shareseq_reorg {
    input {
        # Sample_id
        String sample_id
        # Type of data, choose from ['atac','gex']
        String type
        # R1 Fastq pattern
        String r1_fastq_pattern
        # R2 Fastq pattern
        String r2_fastq_pattern
        # Index Fastq pattern
        String index_fastq_pattern
        # Input FASTQ directory, gs url
        String input_fastqs_directories
        # Shareseq reorg output directory, gs url
        String output_directory

        # shareseq_reorg version
        String shareseq_reorg_version = "0.1.0"
        # Which docker registry to use
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per reorg job
        Int num_cpu = 4
        # Memory string, e.g. 120G
        String memory = "120G"
        # Disk space in GB
        Int disk_space = 1500
        # Number of preemptible tries
        Int preemptible = 2
        # Max number of retries for AWS instance
        Int awsMaxRetries = 5
        # Backend
        String backend = "gcp"
    }

    call run_shareseq_reorg {
        input:
            sample_id = sample_id,
            type = type,
            input_fastqs_directories = sub(input_fastqs_directories, "/+$", ""),
            r1_fastq_pattern = r1_fastq_pattern,
            r2_fastq_pattern = r2_fastq_pattern,
            index_fastq_pattern = index_fastq_pattern,
            output_directory = sub(output_directory, "/+$", ""),
            shareseq_reorg_version = shareseq_reorg_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend
    }

    output {
        String output_reorg_directory = run_shareseq_reorg.output_reorg_directory
        File monitoringLog = run_shareseq_reorg.monitoringLog
    }
}

task run_shareseq_reorg {
    input {
        String sample_id
        String type
        String r1_fastq_pattern
        String r2_fastq_pattern
        String index_fastq_pattern
        String input_fastqs_directories
        String output_directory
        String shareseq_reorg_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        
        mkdir -p _out_reorg

        python <<CODE
        import os, re
        from fnmatch import fnmatch
        from subprocess import check_call, CalledProcessError, DEVNULL, STDOUT

        target_dirs = []
        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            target = "~{sample_id}_" + str(i)
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/~{sample_id}']
                print(' '.join(call_args))
                check_call(call_args, stdout=DEVNULL, stderr=STDOUT)
                call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)
                target_dirs.append(target)

            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{r1_fastq_pattern}' , target + '/']
                print(' '.join(call_args))
                check_call(call_args)

                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{r2_fastq_pattern}' , target + '/']
                print(' '.join(call_args))
                check_call(call_args)

                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}~{index_fastq_pattern}' , target + '/']
                print(' '.join(call_args))
                check_call(call_args)
                target_dirs.append(target)               


        call_args = ['shareseq_reorg_barcodes', '/indices/shareseq_barcode_index.csv', '/indices/shareseq_flanking_sequence.csv',
                     '~{sample_id}', '~{type}', target_dirs, '_out_reorg',
                     '--r1-pattern', '~{r1_fastq_pattern}', '--r2-pattern', '~{r2_fastq_pattern}',
                     '--r3-pattern', '~{index_fastq_pattern}']
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        strato sync --backend ~{backend} -m _out_reorg ~{output_directory}/fastqs_reorg/~{sample_id}
    }

    output {
        String output_reorg_directory = "~{output_directory}/fastqs_reorg/~{sample_id}"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/shareseq_reorg:~{shareseq_reorg_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
