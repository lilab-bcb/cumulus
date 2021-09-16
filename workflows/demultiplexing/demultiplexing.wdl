version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:demuxEM/versions/6/plain-WDL/descriptor" as dem
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:souporcell/versions/6/plain-WDL/descriptor" as soc
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:popscle/versions/6/plain-WDL/descriptor" as psc

workflow demultiplexing {
    input {
        # Input CSV file describing metadata of RNA and hashtag/genetic data pairing.
        File input_sample_sheet
        # This is the output directory (gs url + path) for all results. There will be one folder per RNA-hashtag/genetic data pair under this directory.
        String output_directory
        # Reference genome name
        String genome = ""
        # demultiplexing algorithm to use for genetic-pooling data, choosing from 'souporcell' or 'popscle' (demuxlet/freemuxlet)
        String demultiplexing_algorithm = "souporcell"
        # Only demultiplex cells/nuclei with at least <min_num_genes> expressed genes
        Int min_num_genes = 100

        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry = "quay.io/cumulus"
        # Number of preemptible tries
        Int preemptible = 2
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Backend
        String backend = "gcp"

        # For demuxEM
        # The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse. [default: 0.0]
        Float? demuxEM_alpha_on_samples
        # Only demultiplex cells/nuclei with at least <demuxEM_min_num_umis> of UMIs. [default: 100]
        Int? demuxEM_min_num_umis
        # Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown. [default: 10.0]
        Float? demuxEM_min_signal_hashtag
        # The random seed used in the KMeans algorithm to separate empty ADT droplets from others. [default: 0]
        Int? demuxEM_random_state
        # Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc. [default: true]
        Boolean demuxEM_generate_diagnostic_plots = true
        # Generate violin plots using gender-specific genes (e.g. Xist). <demuxEM_generate_gender_plot> is a comma-separated list of gene names
        String? demuxEM_generate_gender_plot
        # Reference Index TSV
        File ref_index_file = "gs://regev-lab/resources/cellranger/index.tsv"
        # DemuxEM version
        String demuxEM_version = "0.1.6"
        # Number of CPUs used
        Int demuxEM_num_cpu = 8
        # Disk space in GB
        Int demuxEM_disk_space = 20
        # Memory in GB
        Int demuxEM_memory = 10

        # For souporcell
        # Number of expected clusters when doing clustering
        Int souporcell_num_clusters = 1
        # If true, run souporcell in de novo mode without reference genotypes; and if a reference genotype vcf file is provided in the sample sheet, use it only for matching the cluster labels computed by souporcell. If false, run souporcell with --known_genotypes option using the reference genotype vcf file specified in sample sheet, and souporcell_rename_donors is required in this case.
        Boolean souporcell_de_novo_mode = true
        # Users can provide a common variants list in VCF format for Souporcell to use, instead of calling SNPs de novo
        File? souporcell_common_variants
        # Skip remap step. Only recommended in non denovo mode or common variants are provided
        Boolean souporcell_skip_remap = false
        # A comma-separated list of donor names for renaming clusters achieved by souporcell
        String souporcell_rename_donors = ""
        # Souporcell version to use. Available versions: "2020.07", "2021.03", "2020.03"
        String souporcell_version = "2020.07"
        # Number of CPUs to request for souporcell per pair
        Int souporcell_num_cpu = 32
        # Disk space (integer) in GB needed for souporcell per pair
        Int souporcell_disk_space = 500
        # Memory size (integer) in GB needed for souporcell per pair
        Int souporcell_memory = 120

        # For popscle (demuxlet/freemuxlet)
        # Minimum mapping quality to consider (lower MQ will be ignored) [default: 20]
        Int? popscle_min_MQ
        # Minimum distance to the tail (lower will be ignored) [default: 0]
        Int? popscle_min_TD
        # Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB  [default: "CB"]
        String? popscle_tag_group
        # Tag representing UMIs. For 10x genomics, use UB  [default: "UB"]
        String? popscle_tag_UMI
        # Default is 0, means to use demuxlet, if this number > 0, use freemuxlet
        Int popscle_num_samples = 0
        # FORMAT field to extract the genotype, likelihood, or posterior from
        String popscle_field = "GT"
        # Grid of alpha to search for [default: "0.1,0.2,0.3,0.4,0.5"]
        String? popscle_alpha
        # Popscle version. Available versions: "2021.05", "0.1b"
        String popscle_version = "0.1b"
        # A comma-separated list of donor names for renaming clusters achieved by freemuxlet
        String? popscle_rename_donors
        # Number of CPUs used for popscle per pair
        Int popscle_num_cpu = 1
        # Memory size in GB needed for popscle per pair
        Int popscle_memory = 120
        # Extra disk space (integer) in GB needed for popscle per pair
        Int popscle_extra_disk_space = 100

        # Version of config docker image to use. This docker is used for parsing the input sample sheet for downstream execution. Available options: "0.2", "0.1"
        String config_version = "0.2"
    }

    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    call generate_demux_config as Config {
        input:
            input_sample_sheet = input_sample_sheet,
            config_version = config_version,
            docker_registry = docker_registry,
            zones = zones,
            preemptible = preemptible
    }

    if (length(Config.hashing_ids) > 0) {
        scatter (hashing_id in Config.hashing_ids) {
            call dem.demuxEM as demuxEM {
                input:
                    sample_id = hashing_id,
                    output_directory = output_directory_stripped,
                    input_rna_h5 = Config.id2rna[hashing_id],
                    input_hto_csv = Config.id2tag[hashing_id],
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
                    preemptible = preemptible,
                    backend = backend
            }
        }
    }

    if (length(Config.pooling_ids) > 0) {
        scatter (pooling_id in Config.pooling_ids) {
            if (demultiplexing_algorithm == "souporcell") {
                Map[String, String] ref_index2gsurl = read_map(ref_index_file)
                String genome_url = ref_index2gsurl[genome]

                call soc.souporcell as souporcell {
                    input:
                        sample_id = pooling_id,
                        output_directory = output_directory_stripped,
                        input_rna = Config.id2rna[pooling_id],
                        input_bam = Config.id2tag[pooling_id],
                        genome_url = genome_url,
                        ref_genotypes_url = Config.id2genotype[pooling_id],
                        common_variants = souporcell_common_variants,
                        skip_remap = souporcell_skip_remap,
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
                        preemptible = preemptible,
                        backend = backend
                }
            }

            if (demultiplexing_algorithm == "popscle") {
                call psc.popscle as popscle {
                    input:
                        sample_id = pooling_id,
                        output_directory = output_directory_stripped,
                        input_rna = Config.id2rna[pooling_id],
                        input_bam = Config.id2tag[pooling_id],
                        ref_genotypes = Config.id2genotype[pooling_id],
                        min_num_genes = min_num_genes,
                        min_MQ = popscle_min_MQ,
                        min_TD = popscle_min_TD,
                        tag_group = popscle_tag_group,
                        tag_UMI = popscle_tag_UMI,
                        field = popscle_field,
                        alpha = popscle_alpha,
                        nsample = popscle_num_samples,
                        donor_rename = popscle_rename_donors,
                        docker_registry = docker_registry,
                        popscle_version = popscle_version,
                        extra_disk_space = popscle_extra_disk_space,
                        num_cpu = popscle_num_cpu,
                        memory = popscle_memory,
                        zones = zones,
                        preemptible = preemptible,
                        backend = backend
                }
            }
        }
    }

    output {
        Array[String] output_folders = select_all(flatten(select_all([demuxEM.output_folder, souporcell.output_folder, popscle.output_folder])))
        Array[File] output_zarr_files = select_all(flatten(select_all([demuxEM.output_zarr, souporcell.output_zarr, popscle.output_zarr])))
    }
}

task generate_demux_config {
    input {
        File input_sample_sheet
        String config_version
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

        with open('hashing_ids.txt', 'w') as fo_hashing, open('pooling_ids.txt', 'w') as fo_pooling, open('id2rna.txt', 'w') as fo_rnas, open('id2tag.txt', 'w') as fo_tags, open('id2genotype.txt', 'w') as fo_genotypes:

            for idx, row in df.iterrows():
                if row['TYPE'] in ['genetic-pooling', 'genetic_pooling']:
                    fo_pooling.write(row['OUTNAME'] + '\n')
                else:
                    if row['TYPE'] not in ['cell-hashing', 'nucleus-hashing', 'cell_hashing', 'nucleus_hashing']:
                        print("Warning: demultiplexing type " + row['TYPE'] + " is neither cell-hashing nor nucleus-hashing. But we still assume hashing data and use DemuxEM for demultiplexing!")
                    fo_hashing.write(row['OUTNAME'] + '\n')

                fo_rnas.write(row['OUTNAME'] + '\t' + row['RNA'] + '\n')
                fo_tags.write(row['OUTNAME'] + '\t' + row[tag_key] + '\n')

                if 'Genotype' in df.columns and (not pd.isnull(row['Genotype'])):
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
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: "~{preemptible}"
    }
}
