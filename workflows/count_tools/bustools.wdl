version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/alexandria:kallisto-bustools_count/versions/2/plain-WDL/descriptor" as kbc
# import "../../../kallisto-bustools_workflow/WDL/kallisto-bustools_count.wdl" as kbc

workflow bustools {
    input {
        String sample_id
        File r1_fastq
        File r2_fastq
        String genome_url
        String chemistry
        String output_directory
        Int num_cpu
        Boolean output_loom
        Boolean output_h5ad  

        String docker_registry = "cumulusprod"
        String bustools_version = '0.24.4'
        Int disk_space = 500
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int memory = 120
    }
    
    String chemistry_str = if chemistry == 'tenX_v3' then '10XV3' else if chemistry == 'tenX_v2' then '10XV2' else ''

    call set_up_resources as src {
        input:
            output_directory = output_directory,
            genome = genome_url,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }

    call kbc.kallisto_bustools_count as kallisto_bustools_count {
        input:
            docker = "shaleklab/kallisto-bustools:" + bustools_version,
            number_cpu_threads = num_cpu,
            task_memory_GB = memory,
            preemptible = preemptible,
            zones = zones,
            boot_disk_size_GB = 12,
            disks = "local-disk " + disk_space + " HDD",
            bucket = src.bucket,
            output_path = src.output_path,
            index = src.output_index,
            T2G_mapping = src.output_t2g,
            technology = chemistry_str,
            R1_fastq = r1_fastq,
            sample_name = sample_id,
            R2_fastq = r2_fastq,
            use_lamanno = false,
            bustools_filter = false,
            loom = output_loom,
            h5ad = output_h5ad,
            delete_bus_files = false
    }

    output {
        String output_folder = kallisto_bustools_count.count_output_path
    }
}

task set_up_resources {
    input {
        String output_directory
        File genome
        String docker_registry
        String zones
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp

        mkdir genome_ref
        tar -zxf "~{genome}" -C genome_ref --strip-components 1
        rm "~{genome}"

        python <<CODE
        import os

        dir_list = '~{output_directory}'[5:].split('/')
        bucket = 'gs://' + dir_list[0]
        output_path = '/'.join(dir_list[1:])

        with open('bucket.txt', 'w') as fo1, open('output_path.txt', 'w') as fo2:
            fo1.write(bucket)
            fo2.write(output_path)

        CODE
    }

    output {
        String bucket = read_string('bucket.txt')
        String output_path = read_string('output_path.txt')
        File output_index = 'genome_ref/transcriptome.idx'
        File output_t2g = 'genome_ref/transcripts_to_genes.txt'
        File output_fasta_file = 'genome_ref/cdna.fa'
    }

    runtime {
        docker: "~{docker_registry}/count"
        zones: zones
        preemptible: "~{preemptible}"
    }
}