version 1.0

workflow bam2fastq {
    input {
        File input_bam
        String sample_id
        String output_directory

        String docker_registry = "quay.io/cumulus"
        String samtools_version = "1.11"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int num_cpu = 4
        Int memory = 200
        Int disk_space = 200
        Int preemptible = 2
    }

    call convert_to_fastq {
        input:
            input_bam = input_bam,
            sample_id = sample_id,
            output_directory = output_directory,
            docker_registry = docker_registry,
            samtools_version = samtools_version,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible
    }

    output {
        Array[String] output_fastqs = convert_to_fastq.output_fastqs
        File monitorLog = convert_to_fastq.monitorLog
    }
}

task convert_to_fastq {
    input {
        File input_bam
        String sample_id
        String output_directory

        String docker_registry
        String samtools_version
        String zones
        Int num_cpu
        Int memory
        Int disk_space
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        samtools fastq -n -0 ~{sample_id}_sam.fastq -T on,QX,OQ,CR,CY ~{input_bam}
        python /software/process_fastq.py ~{sample_id}_sam.fastq ~{sample_id}
        gzip ~{sample_id}_R1.fastq
        gzip ~{sample_id}_R2.fastq

        gsutil -q -m cp ~{sample_id}_R1.fastq.gz ~{sample_id}_R2.fastq.gz ~{output_directory}/
    }

    output {
        Array[String] output_fastqs = ["~{sample_id}_R1.fastq.gz", "~{sample_id}_R2.fastq.gz"]
        File monitorLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/samtools:~{samtools_version}"
        zones: zones
        memory: "~{memory}G"
        disks: "local_disk ~{disk_space} HDD"
        cpu: "~{num_cpu}"
        preemptible: "~{preemptible}"
    }
}
