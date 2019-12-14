workflow smartseq2_create_reference {
    File fasta
    File gtf
    # Output directory, gs URL
    String output_directory
    String genome
    String? smartseq2_version = "1.0.0"
    String? zones = "us-central1-b"
    Int? cpu = 8
    String? memory = "7.2G"
    # disk space in GB
    Int? extra_disk_space = 15
    # Number of preemptible tries 
    Int? preemptible = 2
    String? docker_registry = "cumulusprod/"

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "/+$", "")

    call rsem_prepare_reference {
        input:
            fasta=fasta,
            gtf=gtf,
            output_dir = output_directory_stripped,
            genome = genome,
            smartseq2_version=smartseq2_version,
            zones=zones,
            preemptible=preemptible,
            cpu=cpu,
            memory=memory,
            disk_space=extra_disk_space,
            docker_registry=docker_registry
    }
}

task rsem_prepare_reference {
    File fasta
    File gtf
    String output_dir
    String genome
    String smartseq2_version
    String zones
    Int preemptible
    Int cpu
    String memory
    Int disk_space
    String docker_registry

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        mkdir ${genome}
        rsem-prepare-reference --gtf ${gtf} --bowtie2 -p ${cpu} ${fasta} ${genome}/${genome}
        tar -czf ${genome}.tar.gz ${genome}

        gsutil -m cp ${genome}.tar.gz ${output_dir}
        # mkdir -p ${output_dir}
        # cp ${genome}.tar.gz ${output_dir}
    }

    output {
        File reference = "${genome}.tar.gz"
    }

    runtime {
        disks: "local-disk " + ceil(disk_space + 8*size(fasta,"GB") + size(gtf,"GB")) + " HDD"
        docker: "${docker_registry}smartseq2:${smartseq2_version}"
        zones: zones
        preemptible: "${preemptible}"
        cpu:"${cpu}"
        memory:"${memory}"
    }
}

