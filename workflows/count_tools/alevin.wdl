workflow alevin {
    File input_sample_sheet
    File ref_genome_fa
    File ref_transcriptome_fa
    File ref_genome_anno_fa
    String resource_directory
    String index_directory
    String output_directory
    Int? num_cpu = 32
    String? library_type = "ISR"
    String protocol

    String? docker_registry = "cumulusprod"
    String? alevin_version = '1.1'
    Int? disk_space = 100
    Int? preemptible = 2
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    Int? memory = 32

    call set_up_resources as src {
        input:
            input_genome_fa = ref_genome_fa,
            input_transcriptome_fa = ref_transcriptome_fa,
            input_genome_anno_fa = ref_genome_anno_fa,
            output_directory = resource_directory,
            docker_registry = docker_registry,
            alevin_version = alevin_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    call run_salmon_index {
        input:
            input_gentrome_fa = src.output_gentrome_fa,
            input_decoys_txt = src.output_decoys_txt,
            num_cpu = num_cpu,
            output_directory = index_directory,
            docker_registry = docker_registry,
            alevin_version = alevin_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    call run_alevin {
        input:
            input_sample_sheet = input_sample_sheet,
            library_type = library_type,
            protocol = protocol,
            index_directory = run_salmon_index.output_salmon_index,
            num_cpu = num_cpu,
            output_directory = output_directory,
            ref_tgMap_tsv = src.output_t2g_tsv,
            docker_registry = docker_registry,
            alevin_version = alevin_version,
            disk_space = disk_space,
            zones = zones,
            memory = memory,
            preemptible = preemptible
    }

    output {
        File output_mat_gz = run_alevin.output_mat_gz
        File output_tier_mat_gz = run_alevin.output_tier_mat_gz
        File output_genes_txt = run_alevin.output_genes_txt
        File output_cells_txt = run_alevin.output_cells_txt
    }

}

task set_up_resources {
    File input_genome_fa
    File input_transcriptome_fa
    File input_genome_anno_fa
    String output_directory

    String docker_registry
    String alevin_version
    Int disk_space
    String zones
    Int memory
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        make_trascripts_fa.sh ${input_genome_fa} ${input_transcriptome_fa}
        make_tgmap.sh ${input_genome_anno_fa}

        gsutil -m cp gentrome.fa.gz decoys.txt txp2gene.tsv ${output_directory}
        # mkdir -p ${output_directory}
        # cp gentrome.fa.gz decoys.txt txp2gene.tsv ${output_directory}
    }

    output {
        File output_gentrome_fa = 'gentrome.fa.gz'
        File output_decoys_txt = 'decoys.txt'
        File output_t2g_tsv = 'txp2gene.tsv'
    }

    runtime {
        docker: "${docker_registry}/alevin:${alevin_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: 1
        preemptible: "${preemptible}"
    }
}

task run_salmon_index {
    File input_gentrome_fa
    File input_decoys_txt
    Int num_cpu
    String output_directory

    String docker_registry
    String alevin_version
    Int disk_space
    String zones
    Int memory
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        salmon index -t ${input_gentrome_fa} -d ${input_decoys_txt} -p ${num_cpu} -i salmon_index --gencode 

        gsutil -m cp -r salmon_index/* ${output_directory}
        # mkdir -p ${output_directory}
        # cp -r salmon_index/* ${output_directory}
    }

    output {
        File output_salmon_index = 'salmon_index'
    }

    runtime {
        docker: "${docker_registry}/alevin:${alevin_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${preemptible}"
    }
}

task run_alevin {
    File input_sample_sheet
    String library_type
    String protocol
    String index_directory
    Int num_cpu
    String output_directory
    File ref_tgMap_tsv

    String docker_registry
    String alevin_version
    Int disk_space
    String zones
    Int memory
    Int preemptible

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        import pandas as pd

        df = pd.read_csv("${input_sample_sheet}")
        r1_fastqs = df['R1_fastq'].values
        r2_fastqs = df['R2_fastq'].values

        call_args = ['salmon', 'alevin', '-l', '${library_type}', '-1']
        call_args.extend(r1_fastqs)
        call_args.append('-2')
        call_args.extend(r2_fastqs)
        call_args.extend(['-i', '${index_directory}'])

        if '${protocol}' is 'dropseq':
            call_args.append('--dropseq')
        elif '${protocol}' is 'V2':
            call_args.append('--chromium')
        else:
            call_args.append('--chromiumV3')

        call_args.extend(['-p', '${num_cpu}', '-o', 'alevin_output', '--tgMap', '${ref_tgMap_tsv}'])

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        gsutil -m cp -r alevin_output/* ${output_directory}
        # mkdir -p ${output_directory}
        # cp -r alevin_output/* ${output_directory}
    }

    output {
        File output_mat_gz = 'alevin_output/alevin/quants_mat.gz'
        File output_tier_mat_gz = 'alevin_output/alevin/quants_tier_mat.gz'
        File output_genes_txt = 'alevin_output/alevin/quants_mat_cols.txt'
        File output_cells_txt = 'alevin_output/alevin/quants_mat_rows.txt'
    }

    runtime {
        docker: "${docker_registry}/alevin:${alevin_version}"
        zones: zones
        memory: "${memory}G"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_cpu}"
        preemptible: "${preemptible}"
    }
}