version 1.0

workflow cellranger_vdj_create_reference {
    input {
        # Output directory, gs URL
        String output_directory

        # Input genome reference in either FASTA or FASTA.gz format
        File input_fasta
        # Input gene annotation file in either GTF or GTF.gz format
        File input_gtf
        # Genome reference name. New reference will be stored in a folder named genome
        String genome
        # reference version string
        String ref_version = ""

        # Which docker registry to use
        String docker_registry = "quay.io/cumulus"
        # 6.1.1, 6.0.2, 6.0.1
        String cellranger_version = "6.1.1"

        # Disk space in GB
        Int disk_space = 100
        # Google cloud zones
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Memory string
        String memory = "32G"

        # Number of preemptible tries
        Int preemptible = 2
        # Max number of retries for AWS instance
        Int awsMaxRetries = 5
        # Backend
        String backend = "gcp"
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    call run_cellranger_vdj_create_reference {
        input:
            docker_registry = docker_registry,
            cellranger_version = cellranger_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend,
            input_fasta = input_fasta,
            input_gtf = input_gtf,
            output_dir = output_directory_stripped,
            genome = genome,
            ref_version = ref_version
    }

    output {
        File output_reference = run_cellranger_vdj_create_reference.output_reference
    }
}

task run_cellranger_vdj_create_reference {
    input {
        String docker_registry
        String cellranger_version
        Int disk_space
        String zones
        String memory
        Int preemptible
        Int awsMaxRetries
        String backend
        File input_fasta
        File input_gtf
        String output_dir
        String genome
        String ref_version
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        import os
        from subprocess import check_call

        # Unzip fa if needed.
        fa_file = '~{input_fasta}'
        root, ext = os.path.splitext(fa_file)
        if ext == '.gz':
            call_args = ['gunzip', '-f', fa_file]
            print(' '.join(call_args))
            check_call(call_args)
            fa_file = root

        # Unzip gtf if needed.
        gtf_file = '~{input_gtf}'
        root, ext = os.path.splitext(gtf_file)
        if ext == '.gz':
            call_args = ['gunzip', '-f', gtf_file]
            print(' '.join(call_args))
            check_call(call_args)
            gtf_file = root

        call_args = ['cellranger', 'mkvdjref', '--genome=~{genome}', '--fasta=' + fa_file, '--genes=' + gtf_file]

        if '~{ref_version}' is not '':
            call_args.append('--ref-version=~{ref_version}')

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        tar -czf ~{genome}.tar.gz ~{genome}
        strato cp --backend ~{backend} -m ~{genome}.tar.gz "~{output_dir}"
    }

    output {
        File output_reference = "~{genome}.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger:~{cellranger_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
