version 1.0

import "https://raw.githubusercontent.com/lilab-bcb/skylab/cumulus/pipelines/optimus/Optimus.wdl" as opm
# import "../../../skylab/pipelines/optimus/Optimus.wdl" as opm

workflow optimus_count {
    input {
        String sample_id
        File r1_fastq
        File r2_fastq
        File? i1_fastq
        String genome_url
        String chemistry
        String output_directory
        Boolean output_loom

        String docker_registry = "cumulusprod"
        String version = 'optimus_v1.4.0'
        Int disk_space = 100
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int memory = 32
    }
    
    File wl_index_file = "gs://regev-lab/resources/count_tools/whitelist_index.tsv"
    # File wl_index_file = "whitelist_index.tsv"
    Map[String, String] wl_index2gsurl = read_map(wl_index_file)
    String whitelist_url = wl_index2gsurl[chemistry]
    
    call get_reference as ref {
        input:
            genome_url = genome_url,
            docker_registry = docker_registry,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    call opm.Optimus as Optimus {
        input:
            version = version,
            sample_id = sample_id,
            r1_fastq = [r1_fastq],
            r2_fastq = [r2_fastq],
            i1_fastq = if defined(i1_fastq) then [i1_fastq] else '',
            tar_star_reference = ref.star_gz,
            annotations_gtf = ref.gtf,
            ref_genome_fasta = ref.fasta,
            whitelist = whitelist_url,
            chemistry = chemistry,
            fastq_suffix = '.gz',
            output_loom = output_loom
    }

    call organize_result {
        input:
            output_directory = output_directory,
            sample_id = sample_id,
            results = [Optimus.bam, Optimus.matrix, Optimus.matrix_row_index, Optimus.matrix_col_index, Optimus.cell_metrics, Optimus.gene_metrics, Optimus.cell_calls],
            zarr_files = Optimus.zarr_output_files,
            loom_file = Optimus.loom_output_file,
            docker_registry = docker_registry,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    output {
        String output_folder = organize_result.output_folder
    }
}

task get_reference {
    input {
        String genome_url
        String docker_registry
        Int disk_space
        String zones
        Int memory
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp

        gsutil -q -m cp ${genome_url} optimus.tar.gz
        # cp ${genome_url} optimus.tar.gz
        tar -zxvf optimus.tar.gz
        rm optimus.tar.gz
    }

    output {
        File star_gz = 'optimus-ref/star.tar/gz'
        File fasta = 'optimus-ref/genome.fa'
        File gtf = 'optimus-ref/genes.gtf'
    }

    runtime {
        docker: "${docker_registry}/count"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ${disk_space} HDD"
        cpu: 1
        preemptible: "${preemptible}"
    }

}

task organize_result {
    input {
        String output_directory
        String sample_id
        Array[File] results
        Array[File] zarr_files
        File? loom_file
        String docker_registry
        Int disk_space
        String zones
        Int memory
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp

        mkdir output
        cp ${sep=" " results} output

        mkdir output/zarr_output
        cp ${sep=" " zarr_files} output/zarr_output

        cp ${loom_file} output

        gsutil -q -m rsync -r output ${output_directory}/${sample_id}
        # mkdir -p ${output_directory}/${sample_id}
        # cp -r output/* ${output_directory}/${sample_id}
    }

    output {
        String output_folder = "${output_directory}/${sample_id}"
    }

    runtime {
        docker: "${docker_registry}/count"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ${disk_space} HDD"
        cpu: 2
        preemptible: "${preemptible}"
    }
}