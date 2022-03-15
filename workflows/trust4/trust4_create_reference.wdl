version 1.0

workflow trust4_create_reference {
    input {
        # Which docker registry to use
        String docker_registry = "quay.io/cumulus"
        # trust4 version
        String trust4_version = "master"

        # Disk space in GB
        Int disk_space = 100
        # Google cloud zones
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Memory string
        String memory = "80G"

        # Number of preemptible tries
        Int preemptible = 2
        # Max number of retries for AWS instance
        Int awsMaxRetries = 5
        # Backend
        String backend = "gcp"

        # Reference FASTA file
        File reference_fasta
        # Annotation GTF file
        File annotation_gtf
        # gene name list of interest
        File gene_name_list
        # Species name
        String species
        
        # Output directory, URL
        String output_directory
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    call run_trust4_create_reference {
        input:
            docker_registry = docker_registry,
            trust4_version = trust4_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend,
            reference_fasta = reference_fasta,
            annotation_gtf = annotation_gtf,
            gene_name_list = gene_name_list,
            species = species,           
            output_dir = output_directory_stripped
    }

    output {
        File bcrtcr = run_trust4_create_reference.bcrtcr
        File imgt = run_trust4_create_reference.imgt
    }

}

task run_trust4_create_reference {
    input {
        String docker_registry
        String trust4_version
        Int disk_space
        String zones
        String memory
        Int preemptible
        Int awsMaxRetries
        String backend
        File reference_fasta
        File annotation_gtf
        File gene_name_list
        String species
        String output_dir
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        
        mkdir -p ~{species}

        perl BuildDatabaseFa.pl ~{reference_fasta} ~{annotation_gtf} ~{gene_name_list} > ~{species}/bcrtcr.fa
        perl BuildImgtAnnot.pl ~{species} > ~{species}/IMGT+C.fa
        
        strato sync --backend ~{backend} -m ~{species} ~{output_dir}/
    }

    output {
        File bcrtcr = "~{species}/bcrtcr.fa"
        File imgt = "~{species}/IMGT+C.fa"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/trust4:~{trust4_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}



