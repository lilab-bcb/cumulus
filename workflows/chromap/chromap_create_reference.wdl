version 1.0

workflow chromap_create_reference {
    input {
        # Which docker registry to use
        String docker_registry = "quay.io/cumulus"
        # chromap version
        String chromap_version = "0.2.6"

        # Disk space in GB
        Int disk_space = 100
        # Google cloud zones
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Memory string
        String memory = "80G"

        # Number of preemptible tries
        Int preemptible = 2
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend = "gcp"

        # Genome name, the reference package will be named as genome.tar.gz
        String genome
        # URL for input fasta file
        File input_fasta
        # Kmer
        Int? kmer
        # Minimum window size
        Int? mini_win_size

        # Output directory, URL
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
            awsQueueArn = awsQueueArn,
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
        String awsQueueArn
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
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        mkdir -p ~{genome}

        chromap -i ~{"-k" + kmer} ~{"-w" + mini_win_size} -r ~{input_fasta} -o ~{genome}/ref.index

        mv ~{input_fasta} ~{genome}/ref.fa
        tar -czf ~{genome}.tar.gz ~{genome}
        strato cp -m ~{genome}.tar.gz "~{output_dir}"/
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
        queueArn: awsQueueArn
    }
}
