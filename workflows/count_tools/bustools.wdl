import "https://api.firecloud.org/ga4gh/v1/tools/alexandria:kallisto-bustools_count/versions/1/plain-WDL/descriptor" as kbc
# import "../../../kallisto-bustools_workflow/WDL/kallisto-bustools_count.wdl" as kbc

workflow bustools {
    String sample_id
    File r1_fastq
    File r2_fastq
    String genome_url
    String chemistry
    String output_directory
    Int num_cpu
    Boolean output_loom
    Boolean output_h5ad


    String? docker = "shaleklab/kallisto-bustools"
    String? bustools_version = '0.24.4'
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    String chemistry_str = if chemistry == 'tenX_v3' then '10XV3' else '10XV2'

    call set_up_resources as src {
        input:
            output_directory = output_directory,
            genome_url = genome_url + '/bustools.tar.gz',
            docker_registry = docker_registry,
            kb_python_version = kb_python_version,
            zones = zones,
            preemptible = preemptible
    }

    call kbc.kallisto_bustools_count as kallisto_bustools_count {
        input:
            docker = docker + ':' + kb_python_version,
            number_cpu_threads = num_cpu,
            task_memory_GB = memory / 5 * 6,
            program_memory_GB = memory,
            preemptible = preemptible,
            zones = zones,
            boot_disk_size_GB = 12,
            bucket = src.output_loc[0],
            output_path = src.output_loc[1],
            index = src.output_index,
            T2G_mapping = src.output_t2g,
            technology = chemistry_str,
            R1_fastq = r1_fastq,
            sample_name = sample_id,
            R2_fastq = r2_fastq,
            use_lamanno = false,
            loom = output_loom,
            h5ad = output_h5ad,
            delete_bus_files = false
    }

    output {
        String output_folder = kallisto_bustools_count.count_output_path
    }
}

task set_up_resources {
    String output_directory
    String genome_url
    String docker_registry
    String kb_python_version
    String zones
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp

        gsutil -q -m cp ${genome_url} genome.tar.gz
        tar -zxvf genome.tar.gz
        rm genome.tar.gz

        python <<CODE
        import os

        dir_list = '${output_directory}'[5:].split('/')
        bucket = 'gs://' + dir_list[0]
        output_path = '/'.join(dir_list[1:])

        with open('output_location.tsv', 'w') as fo:
            fo.write(bucket + '\t' + output_path + '\n')

        CODE
    }

    output {
        Array[String] output_loc = read_tsv('output_location.tsv')
        File output_index = 'bustools-ref/transcriptome.idx'
        File output_t2g = 'bustools-ref/transcripts_to_genes.txt'
        File output_fasta_file = 'cdna.fa'

    }

    runtime {
        docker: "${docker_registry}/kallisto-bustools:${bustools_version}"
        zones: zones
        preemptible: "${preemptible}"
    }
}