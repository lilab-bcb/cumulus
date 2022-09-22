version 1.0

workflow starsolo_create_reference {
    input {
        File input_fasta
        File input_gtf
        String genome
        String output_directory

        String docker_registry = "quay.io/cumulus"
        String star_version = "2.7.10a"
        Int num_cpu = 32
        Int disk_space = 100
        String memory = "80G"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String backend = "gcp"
        Int preemptible = 2
        Int awsMaxRetries = 5
        String awsQueueArn = ""
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    call run_starsolo_create_reference {
        input:
            input_fasta = input_fasta,
            input_gtf = input_gtf,
            genome = genome,
            output_directory = output_directory_stripped,
            docker_registry = docker_registry,
            version = star_version,
            num_cpu = num_cpu,
            disk_space = disk_space,
            memory = memory,
            zones = zones,
            backend = backend,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            awsQueueArn = awsQueueArn
    }

    output {
        File output_reference = run_starsolo_create_reference.output_reference
    }
}

task run_starsolo_create_reference {
    input {
        File input_fasta
        File input_gtf
        String genome
        String output_directory
        String docker_registry
        String version
        Int num_cpu
        Int disk_space
        String memory
        String zones
        String backend
        Int preemptible
        Int awsMaxRetries
        String awsQueueArn
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call

        call_args = ['STAR', '--runMode', 'genomeGenerate', '--runThreadN', '~{num_cpu}', '--genomeDir', 'starsolo-ref', '--genomeFastaFiles', '~{input_fasta}', '--sjdbGTFfile', '~{input_gtf}']

        mem_digit = ""
        for c in "~{memory}":
            if c.isdigit():
                mem_digit += c
            else:
                break
        mem_bytes = str(int(mem_digit) * 10**9)
        call_args.extend(['--limitGenomeGenerateRAM', mem_bytes])

        check_call(call_args)
        CODE

        tar -czf ~{genome}-starsolo.tar.gz starsolo-ref
        strato cp --backend ~{backend} -m ~{genome}-starsolo.tar.gz "~{output_directory}"/
    }

    output {
        File output_reference = "~{genome}-starsolo.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/starsolo:~{version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
        queueArn: awsQueueArn
    }
}
