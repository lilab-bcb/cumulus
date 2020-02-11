import "https://raw.githubusercontent.com/HumanCellAtlas/skylab/optimus_v1.4.0_terra/pipelines/optimus/Optimus.wdl" as opm


workflow optimus_count {
    String sample_id
    File r1_fastq
    File r2_fastq
    File i1_fastq
    String genome_url
    String chemistry

    
    Int? star_align_cpu = 32

    String? docker_registry = "cumulusprod"
    String? version = 'optimus_v1.4.0'
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    File wl_index_file = "gs://regev-lab/resources/count_tools/whitelist_index.tsv"
    # File wl_index_file = "whitelist_index.tsv"
    Map[String, String] wl_index2gsurl = read_map(wl_index_file)
    String whitelist_url = wl_index2gsurl[chemistry]
    
    run get_reference as ref {
        input:
            genome_url = genome_url,
            docker_registry = docker_registry,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    run opm.Optimus as optimus {
        input:
            version = version,
            sample_id = sample_id,
            r1_fastq = r1_fastq,
            r2_fastq = r2_fastq,
            i1_fastq = i1_fastq,
            tar_star_reference = ref.star_gz,
            annotations_gtf = ref.gtf,
            ref_genome_fasta = ref.fasta,
            whitelist = whitelist_url,
            chemistry = chemistry,
            fastq_suffix = '.gz',
            output_loom = output_loom,
            optimus.StarAlign = star_align_cpu,
            optimus.AttachBarcodes = preemptible,
            optimus.AttachBarcodesNoIndex = preemptible,
            optimus.CalculateCellMetrics = preemptible,
            optimus.CalculateGeneMetrics = preemptible,
            optimus.CellSortBam = preemptible,
            optimus.CorrectUMItools = preemptible,
            optimus.CreateSparseCountMatrix = preemptible,
            optimus.FastqToUBam = preemptible,
            optimus.GeneSortBam = preemptible,
            optimus.MergeCellMetrics = preemptible,
            optimus.MergeCountFiles = preemptible,
            optimus.MergeGeneMetrics = preemptible,
            optimus.MergeSorted = preemptible,
            optimus.ModifyGtf = preemptible,
            optimus.OptimusZarrConversion = preemptible,
            optimus.OptimusZarrToLoom = preemptible,
            optimus.PreCountSort = preemptible,
            optimus.PreMergeSort = preemptible,
            optimus.PreUMISort = preemptible,
            optimus.RunEmptyDrops = preemptible,
            optimus.SplitBamByCellBarcode = preemptible,
            optimus.StarAlign = preemptible,
            optimus.TagGenes = preemptible
    }

    output {
        File bam = optimus.bam
        File matrix = optimus.matrix
        File matrix_row_index = optimus.matrix_row_index
        File matrix_col_index = optimus.matrix_col_index
        File cell_metrics = optimus.cell_metrics
        File gene_metrics = optimus.gene_metrics
        File cell_calls = optimus.cell_calls
        File? loom_output_file = optimus.loom_output_file
    }
}

task get_reference {
    String genome_url
    String docker_registry
    Int disk_space
    String zones
    Int memory
    Int preemptible

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