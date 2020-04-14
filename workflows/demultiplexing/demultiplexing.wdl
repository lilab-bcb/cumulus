version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:demuxEM/versions/2/plain-WDL/descriptor" as dem
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:souporcell/versions/5/plain-WDL/descriptor" as soc
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:demuxlet/versions/3/plain-WDL/descriptor" as dmx

#import "demuxEM.wdl" as dem
#import "souporcell.wdl" as soc
#import "demuxlet.wdl" as dmx

workflow demultiplexing {
    input {
        File input_sample_sheet
        String output_directory
        String genome
        String demultiplexing_algorithm = "souporcell"
        Int min_num_genes = 500

        String docker_registry = "cumulusprod"
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"

        # For demuxEM
        # The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse.
        Float? demuxEM_alpha_on_samples
        # Only demultiplex cells/nuclei with at least <demuxEM_min_num_umis> of UMIs. [default: 100]
        Int? demuxEM_min_num_umis
        # Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown. [default: 10.0]
        Float? demuxEM_min_signal_hashtag
        # The random seed used in the KMeans algorithm to separate empty ADT droplets from others. [default: 0]
        Int? demuxEM_random_state
        # Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc. [default: true]
        Boolean demuxEM_generate_diagnostic_plots = true
        # Generate violin plots using gender-specific genes (e.g. Xist). <demuxEM_generate_gender_plot> is a comma-separated list of gene names.
        String? demuxEM_generate_gender_plot
        String demuxEM_version = "0.1.1"
        Int demuxEM_num_cpu = 8
        Int demuxEM_disk_space = 20
        Int demuxEM_memory = 10

        # For souporcell
        Int souporcell_num_clusters = 1
        Boolean souporcell_de_novo_mode = true
        String souporcell_rename_donors = ""
        String souporcell_version = "2020.03"
        Int souporcell_num_cpu = 32
        Int souporcell_disk_space = 500
        Int souporcell_memory = 120

        # For demuxlet
        String demuxlet_version = "0.1b"
        Int demuxlet_memory = 10
        Int demuxlet_disk_space = 2

    }

    File ref_index_file = "gs://regev-lab/resources/cellranger/index.tsv"
    # File ref_index_file = "index.tsv"
    Map[String, String] ref_index2gsurl = read_map(ref_index_file)
    String genome_url = ref_index2gsurl[genome]

    call generate_demux_config as Config {
        input:
            input_sample_sheet = input_sample_sheet,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }

    if (Config.hashing_ids[0] != '') {
        scatter (hashing_id in Config.hashing_ids) {
            call dem.demuxEM as demuxEM {
                input:
                    sample_id = hashing_id,
                    output_directory = output_directory,
                    input_rna = Config.id2rna[hashing_id],
                    input_adt_csv = Config.id2tag[hashing_id],
                    genome = genome,
                    alpha_on_samples = demuxEM_alpha_on_samples,
                    min_num_genes = min_num_genes,
                    min_num_umis = demuxEM_min_num_umis,
                    min_signal_hashtag = demuxEM_min_signal_hashtag,
                    random_state = demuxEM_random_state,
                    generate_diagnostic_plots = demuxEM_generate_diagnostic_plots,
                    generate_gender_plot = demuxEM_generate_gender_plot,
                    docker_registry = docker_registry,
                    demuxEM_version = demuxEM_version,
                    zones = zones,
                    num_cpu = demuxEM_num_cpu,
                    memory = demuxEM_memory,
                    disk_space = demuxEM_disk_space,
                    preemptible = preemptible
            }
        }
    }

    if (Config.pooling_ids[0] != '') {
        scatter (pooling_id in Config.pooling_ids) {
            if (demultiplexing_algorithm == "souporcell") {
                call soc.souporcell as souporcell {
                    input:
                        sample_id = pooling_id,
                        output_directory = output_directory,
                        input_rna = Config.id2rna[pooling_id],
                        input_bam = Config.id2tag[pooling_id],
                        genome_url = genome_url,
                        ref_genotypes_url = Config.id2genotype[pooling_id],
                        de_novo_mode = souporcell_de_novo_mode,
                        min_num_genes = min_num_genes,
                        num_clusters = souporcell_num_clusters,
                        donor_rename = souporcell_rename_donors,
                        souporcell_version = souporcell_version,
                        docker_registry = docker_registry,
                        num_cpu = souporcell_num_cpu,
                        disk_space = souporcell_disk_space,
                        memory = souporcell_memory,
                        zones = zones,
                        preemptible = preemptible
                }
            }

            if (demultiplexing_algorithm == "demuxlet") {
                call dmx.demuxlet as demuxlet {
                    input:
                        sample_id = pooling_id,
                        output_directory = output_directory,
                        input_rna = Config.id2rna[pooling_id],
                        input_bam = Config.id2tag[pooling_id],
                        ref_genotypes = Config.id2genotype[pooling_id],
                        min_num_genes = min_num_genes,
                        docker_registry = docker_registry,
                        demuxlet_version = demuxlet_version,
                        extra_disk_space = demuxlet_disk_space,
                        memory = demuxlet_memory,
                        zones = zones,
                        preemptible = preemptible
                }
            }
        }
    }

    output {
        Array[String] output_folders = select_all(flatten(select_all([demuxEM.output_folder, souporcell.output_folder, demuxlet.output_folder])))
        Array[File] output_zarr_files = select_all(flatten(select_all([demuxEM.output_zarr, souporcell.output_zarr, demuxlet.output_zarr])))
    }
}

task generate_demux_config {
    input {
        File input_sample_sheet
        String docker_registry
        String zones
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import re, sys
        import pandas as pd

        df = pd.read_csv('~{input_sample_sheet}', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()
        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['OUTNAME'].str.contains(regex_pat)):
            print('OUTNAME must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        tag_key = 'TagFile' if 'TagFile' in df.columns else 'ADT'
        assert tag_key in ['TagFile', 'ADT']
            
        with open('hashing_ids.txt', 'w') as fo_hashing, open('pooling_ids.txt', 'w') as fo_pooling, open('id2rna.txt', 'w') as fo_rnas, open('id2tag.txt', 'w') as fo_tags, open('id2genotype.txt', 'w') as fo_genotypes:
            
            for idx, row in df.iterrows():
                if row['TYPE'] == 'genetic-pooling':
                    fo_pooling.write(row['OUTNAME'] + '\n')
                else:
                    assert row['TYPE'] in ['cell-hashing', 'nucleus-hashing']
                    fo_hashing.write(row['OUTNAME'] + '\n')

                fo_rnas.write(row['OUTNAME'] + '\t' + row['RNA'] + '\n')
                fo_tags.write(row['OUTNAME'] + '\t' + row[tag_key] + '\n')

                if 'Genotype' in df.columns:
                    fo_genotypes.write(row['OUTNAME'] + '\t' + row['Genotype'] + '\n')
                else:
                    fo_genotypes.write(row['OUTNAME'] + '\tnull\n')
        CODE
    }

    output {
        Array[String] hashing_ids = read_lines('hashing_ids.txt')
        Array[String] pooling_ids = read_lines('pooling_ids.txt')
        Map[String, String] id2rna = read_map('id2rna.txt')
        Map[String, String] id2tag = read_map('id2tag.txt')
        Map[String, String] id2genotype = read_map('id2genotype.txt')
    }

    runtime {
        docker: "~{docker_registry}/config"
        zones: zones
        preemptible: "~{preemptible}"
    }
}