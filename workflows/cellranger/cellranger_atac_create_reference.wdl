version 1.0

workflow cellranger_atac_create_reference {
    input {
        String docker_registry = "cumulusprod"
        String cellranger_atac_version = '1.1.0'
        Int disk_space = 100
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String memory = "32G"

        File config_json
        # Output directory, gs URL
        String output_directory
        String genome
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    call run_cellranger_atac_create_reference {
        input:
            docker_registry = docker_registry,
            cellranger_atac_version = cellranger_atac_version,
            disk_space = disk_space,
            preemptible = preemptible,
            zones = zones,
            memory = memory,
            config_json = config_json,
            output_dir = output_directory_stripped,
            genome = genome
    }

}

task run_cellranger_atac_create_reference {
    input {
        String docker_registry
        String cellranger_atac_version
        Int disk_space
        String zones
        String memory
        Int preemptible

        File config_json
        String output_dir
        String genome
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        cellranger-atac mkref ~{genome} --config ~{config_json}
        tar -czf ~{genome}.tar.gz ~{genome}
        gsutil -m cp ~{genome}.tar.gz "~{output_dir}"
        # mkdir -p "~{output_dir}"
        # cp ~{genome}.tar.gz "~{output_dir}"
    }

    output {
        File output_reference = "~{genome}.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger-atac:~{cellranger_atac_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: "~{preemptible}"
    }
}
