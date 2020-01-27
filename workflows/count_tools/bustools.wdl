workflow bustools {
    File ref_fasta_file
    File ref_gtf_file
    String generate_index_directory
    Int? num_cpu_ref = 1

    File input_sample_sheet
    String chemistry
    String output_directory
    Boolean? output_loom = false
    Boolean? output_h5ad = false
    Int? num_cpu_count = 32


    String? docker_registry = "cumulusprod"
    String? kb_python_version = '0.24'
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    call run_kb_ref {
        input:
            input_fasta_file = ref_fasta_file,
            input_gtf_file = ref_gtf_file,
            output_directory = generate_index_directory,
            num_cpu = num_cpu_ref,
            docker_registry = docker_registry,
            kb_python_version = kb_python_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    call run_kb_count {
        input:
            input_sample_sheet = input_sample_sheet,
            index = run_kb_ref.output_index,
            t2g = run_kb_ref.output_t2g,
            chemistry = chemistry,
            output_loom = output_loom,
            output_h5ad = output_h5ad,
            output_directory = output_directory,
            num_cpu = num_cpu_count,
            docker_registry = docker_registry,
            kb_python_version = kb_python_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File output_mtx_file = run_kb_count.output_mtx_file
        File output_genes_txt = run_kb_count.output_genes_txt
        File output_barcodes_txt = run_kb_count.output_barcodes_txt
        File? output_loom_file = run_kb_count.output_loom_file
        File? output_h5ad_file = run_kb_count.output_h5ad_file
    }
}

task run_kb_ref {
    String docker_registry
    String kb_python_version
    Int disk_space
    String zones
    Int memory
    Int preemptible
    Int num_cpu

    File input_fasta_file
    File input_gtf_file
    String output_directory

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        kb ref -i transcriptome.idx -g transcripts_to_genes.txt -f1 cdna.fa \
        ${input_fasta_file} ${input_gtf_file}
        gsutil cp transcriptome.idx transcripts_to_genes.txt cdna.fa ${output_directory}
        # mkdir -p ${output_directory}
        # cp transcriptome.idx transcripts_to_genes.txt cdna.fa ${output_directory}
    }

    output {
        File output_index = 'transcriptome.idx'
        File output_t2g = 'transcripts_to_genes.txt'
        File output_fasta_file = 'cdna.fa'
    }

    runtime {
        docker: "${docker_registry}/bustools:${kb_python_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${preemptible}"
    }
}

task run_kb_count {
    String docker_registry
    String kb_python_version
    Int disk_space
    String zones
    Int memory
    Int preemptible

    File input_sample_sheet
    File index
    File t2g
    String chemistry
    Boolean output_loom
    Boolean output_h5ad
    String output_directory
    Int num_cpu

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv("${input_sample_sheet}")
        r1_fastqs = df['r1_fastq'].values
        r2_fastqs = df['r2_fastq'].values
        input_fastqs = [item for tup in list(zip(r1_fastqs, r2_fastqs)) for item in tup]

        call_args = ['kb', 'count', '-i', '${index}', '-g', '${t2g}', '-x', '${chemistry}', '-o', 'output', '-t', '${num_cpu}', '-m', '${memory}']
        if '${output_loom}' is 'true':
            call_args.append('--loom')
        if '${output_h5ad}' is 'true':
            call_args.append('--h5ad')
        call_args.extend(input_fastqs)

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil cp -r output ${output_directory} 
        # mkdir -p ${output_directory}
        # cp -r output ${output_directory}
    }

    output {
        File output_mtx_file = 'output/counts_unfiltered/cells_x_genes.mtx'
        File output_genes_txt = 'output/counts_unfiltered/cells_x_genes.genes.txt'
        File output_barcodes_txt = 'output/counts_unfiltered/cells_x_genes.barcodes.txt'
        File? output_loom_file = 'output/counts_unfiltered/adata.loom'
        File? output_h5ad_file = 'output/counts_unfiltered/adata.h5ad'
    }


    runtime {
        docker: "${docker_registry}/bustools:${kb_python_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${preemptible}"
    }
}