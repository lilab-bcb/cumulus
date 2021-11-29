version 1.0

workflow chromap_create_reference {
    input {
        # Which docker registry to use
        String docker_registry = "quay.io/cumulus"
        # cellranger-atac version: 2.0.0, 1.2.0, 1.1.0
        String chromap_version = "0.1.3"

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

        # Genome name
        String genome
        # GSURL for input fasta file
        File input_fasta
        # Kmer
        Int? kmer
        # Minimum window size
        Int? mini_win_size
        
        # Output directory, gs URL
        String output_directory
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    call run_chromap_create_reference {
        input:
            docker_registry = docker_registry,
            chromap_version = chromap_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend,
            genome = genome,
            kmer = kmer,
            mini_win_size = mini_win_size,
            input_fasta = input_fasta,
            output_dir = output_directory_stripped
    }

    output {
        File output_reference = run_chromap_create_reference.output_reference
    }

}

task run_chromap_create_reference {
    input {
        String docker_registry
        String chromap_version
        Int disk_space
        String zones
        String memory
        Int preemptible
        Int awsMaxRetries
        String backend
        String genome
        File input_fasta
        Int? kmer
        Int? mini_win_size
        String output_dir
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call, CalledProcessError
        import sys

        call_args = ['chromap', '-i']

        if '~{kmer}' !=  '':
            call_args.extend(['-k','~{kmer}'])
        if '~{mini_win_size}' !=  '':
            call_args.extend(['-w','~{mini_win_size}'])

        call_args.extend(['-r', '~{input_fasta}', '-o', 'ref.index'])
        print(' '.join(call_args))
        check_call(call_args)    
        CODE

        strato cp --backend ~{backend} -m input_fasta ~{output_dir}
        strato cp --backend ~{backend} -m ref.index ~{output_dir}
        tar -czf ~{genome}.tar.gz ~{output_dir}
    }

    output {
        File output_reference = "~{genome}.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/chromap:~{chromap_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}



