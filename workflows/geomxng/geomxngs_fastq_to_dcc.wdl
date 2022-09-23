version 1.0

workflow geomxngs_fastq_to_dcc {
    input {
        File ini
        String fastq_directory
        String output_directory
        File? fastq_rename
        String docker_registry = "gcr.io/broad-cumulus"
        String geomxngs_version = "2.3.3.10"
        Int cpu = 4
        Int disk_space = 500
        String memory = "64GB"
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f"
        String backend = "gcp"
        String aws_queue_arn = ""
        Boolean delete_fastq_directory = false
        Int max_retries = 0
    }


    output {
        String geomxngs_output = geomxngs_task.geomxngs_output
        String dcc_zip = geomxngs_task.dcc_zip
    }

    call geomxngs_task {
        input:
            fastq_directory = fastq_directory,
            delete_fastq_directory=delete_fastq_directory,
            ini = ini,
            output_directory = output_directory,
            fastq_rename = fastq_rename,
            docker_registry = docker_registry,
            geomxngs_version = geomxngs_version,
            cpu = cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            zones = zones,
            backend = backend,
            aws_queue_arn = aws_queue_arn,
            max_retries=max_retries
    }

}

task geomxngs_task {
    input {
        String fastq_directory
        String output_directory
        File ini
        File? fastq_rename
        String memory
        Int cpu
        Int preemptible
        Int disk_space
        String docker_registry
        String geomxngs_version
        String zones
        String backend
        String aws_queue_arn
        Boolean delete_fastq_directory
        Int max_retries
    }
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")
    String local_ini = basename(ini)
    # update cpus in ini fie

    command {
        set -e
        monitor_script.sh > monitoring.log &

        python /software/scripts/update-cpu.py --ini ~{ini} --cpu ~{cpu} --out ~{local_ini}

        mkdir fastqs
        strato sync --backend ~{backend} -m ~{fastq_directory} fastqs/

        if [[ '~{fastq_rename}' != '' ]]; then
            python /software/scripts/rename-fastqs.py --fastqs fastqs --rename ~{fastq_rename}
        fi

        geomxngspipeline --in=fastqs --out=results --ini=~{local_ini} ~{true="--check-illumina-naming=false" false="" geomxngs_version!="2.2.0.2"}
        strato sync --backend ~{backend} -m results ~{output_directory_stripped}
        if [[ '~{delete_fastq_directory}' = 'true' ]]; then
            python /software/scripts/delete-url.py --backend ~{backend} --url ~{fastq_directory}
        fi

    }

    output {
        File monitoring_log = "monitoring.log"
        String geomxngs_output = "${output_directory_stripped}"
        String local_dcc_zip = glob("results/*.zip")[0]
        String dcc_zip = output_directory_stripped  + "/" + basename(local_dcc_zip)
    }

    runtime {
        docker: "~{docker_registry}/geomxngspipeline:~{geomxngs_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: cpu
        preemptible: preemptible
        queueArn: aws_queue_arn
        maxRetries:max_retries
    }

}
