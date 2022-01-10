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
        # Number of cpus per bcl2fastq job
        Int num_cpu = 32
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
        String output_directory = run_shareseq_reorg.output_directory
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

        shareseq_reorg_barcodes /indices/shareseq_barcode_index.csv /indices/shareseq_flanking_sequence.csv \
                                ~{sample_id} ~{type} ~{input_fastqs_directories} _out_reorg \
                                ~{'--r1-pattern '+ r1_fastq_pattern} ~{'--r2-pattern '+ r2_fastq_pattern} \
                                ~{'--r3-pattern '+ index_fastq_pattern}

        strato sync --backend ~{backend} -m _out_reorg ~{output_directory}/~{sample_id}_fastqs_reorg
    }

    output {
        String output_reorg_directory = "~{output_directory}/~{sample_id}_fastqs_reorg"
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
