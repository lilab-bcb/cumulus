version 1.0

import "cellranger_count.wdl" as crc
import "cellranger_multi.wdl" as crmulti
import "cellranger_vdj.wdl" as crv
import "../cumulus/cumulus_adt.wdl" as ca
import "cellranger_atac_count.wdl" as crac
import "cellranger_arc_count.wdl" as crarc

workflow cellranger_workflow {
    input {
        # Columns: Sample, Reference, Flowcell, [Chemistry, DataType, FeatureBarcodeFile, Link]).
        File input_csv_file
        # Output directory, AWS or GCP URI
        String output_directory

        # For cellranger count

        # Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells.
        Int? force_cells
        # Expected number of recovered cells. Mutually exclusive with force_cells
        Int? expect_cells
        # If count reads mapping to intronic regions
        Boolean include_introns = true
        # If generate bam outputs. This is also a spaceranger argument.
        Boolean no_bam = false
        # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false.
        Boolean secondary = false
        # Cell annotation model to use. Valid model names: auto, human_pca_v1_beta, mouse_pca_v1_beta. Default: Do not run.
        String cell_annotation_model = ""

        # For vdj

        # Do not align reads to reference V(D)J sequences before de novo assembly. Default: false
        Boolean vdj_denovo = false

        # Force the analysis to be carried out for a particular chain type. The accepted values are:
        #   "auto" for auto detection based on TR vs IG representation (default),
        #   "TR" for T cell receptors,
        #   "IG" for B cell receptors,
        # Use this in rare cases when automatic chain detection fails.
        String vdj_chain = "auto"
        # If inner enrichment primers other than those provided in the 10x kits are used, they need to be specified here as a textfile with one primer per line. Disable secondary analysis, e.g. clustering
        # A cloud URI to the text file
        String vdj_inner_enrichment_primers = ""

        # For extracting ADT count

        # Barcode start position at Read 2 (0-based coordinate) for CRISPR
        Int? crispr_barcode_pos
        # scaffold sequence for CRISPR, default is ""
        String? scaffold_sequence
        # maximum hamming distance in feature barcodes (change default to 2)
        Int max_mismatch = 2
        # minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for data type crispr
        Float min_read_ratio = 0.1

        # For atac

        # For atac, choose the algorithm for dimensionality reduction prior to clustering and tsne: 'lsa' (default), 'plsa', or 'pca'.
        String? atac_dim_reduce = "lsa"
        # A BED file to override peak caller
        File? peaks

        # For arc

        # Disable counting of intronic reads.
        Boolean arc_gex_exclude_introns = false
        # Cell caller override: define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode.
        Int? arc_min_atac_count
        # Cell caller override: define the minimum number of GEX UMI counts for a cell barcode.
        Int? arc_min_gex_count

        # For multi

        # CMO set CSV file, delaring CMO constructs and associated barcodes
        File? cmo_set

        # Index TSV file
        File acronym_file = "gs://cumulus-ref/resources/cellranger/index.tsv"

        # 9.0.1, 9.0.0, 8.0.1, 8.0.0, 7.2.0, 7.1.0, 7.0.1, 7.0.0
        String cellranger_version = "9.0.1"
        String cumulus_feature_barcoding_version = "0.11.4"
        # 2.1.0, 2.0.0
        String cellranger_atac_version = "2.1.0"
        # 2.0.2.strato, 2.0.2.custom-max-cell, 2.0.2, 2.0.1, 2.0.0
        String cellranger_arc_version = "2.0.2.strato"
        # config version
        String config_version = "0.3"

        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry = "quay.io/cumulus"
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Backend
        String backend = "gcp"
        # Number of cpus per cellranger and spaceranger job
        Int num_cpu = 32
        # Memory string
        String memory = "120G"

        # Number of cpus for cellranger-atac count
        Int atac_num_cpu = 64
        # Memory string for cellranger-atac count
        String atac_memory = "57.6G"

        # Number of cpus for cumulus-adt
        Int feature_num_cpu = 4
        # Optional memory string for cumulus_adt
        String feature_memory = "32G"

        # Number of cpus for cellranger-arc count
        Int arc_num_cpu = 64
        # Memory string for cellranger-arc count
        String arc_memory = "160G"

        # Optional disk space needed for cell ranger count.
        Int count_disk_space = 500
        # Optional disk space needed for cell ranger multi.
        Int multi_disk_space = 1500
        # Optional disk space needed for cell ranger vdj.
        Int vdj_disk_space = 500
        # Optional disk space needed for cumulus_adt
        Int feature_disk_space = 100
        # Optional disk space needed for cellranger-atac count
        Int atac_disk_space = 500
        # Optional disk space needed for cellranger-arc count
        Int arc_disk_space = 700

        # Number of preemptible tries
        Int preemptible = 2
        # Arn string of AWS queue to use
        String awsQueueArn = ""
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    String docker_registry_stripped = sub(docker_registry, "/+$", "")

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    String null_file = acronym2gsurl["null_file"]


    call generate_count_config {
        input:
            input_csv_file = input_csv_file,
            output_dir = output_directory_stripped,
            config_version = config_version,
            docker_registry = docker_registry_stripped,
            zones = zones,
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend,
            null_file = null_file
    }

    if (length(generate_count_config.sample_ids) > 0) {
        scatter (sample_id in generate_count_config.sample_ids) {
            call crc.cellranger_count as cellranger_count {
                input:
                    sample_id = sample_id,
                    input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                    output_directory = output_directory_stripped,
                    genome = generate_count_config.sample2genome[sample_id],
                    target_panel = generate_count_config.sample2fbf[sample_id],
                    chemistry = generate_count_config.sample2chemistry[sample_id],
                    include_introns = include_introns,
                    acronym_file = acronym_file,
                    no_bam = no_bam,
                    secondary = secondary,
                    cell_annotation_model = cell_annotation_model,
                    force_cells = force_cells,
                    expect_cells = expect_cells,
                    cellranger_version = cellranger_version,
                    docker_registry = docker_registry_stripped,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = count_disk_space,
                    preemptible = preemptible,
                    backend = backend,
                    awsQueueArn = awsQueueArn
            }
        }

        call collect_summaries {
            input:
                summaries = cellranger_count.output_metrics_summary,
                sample_ids = cellranger_count.output_count_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }
    }

    if (length(generate_count_config.sample_vdj_ids) > 0) {
        scatter (sample_id in generate_count_config.sample_vdj_ids) {
            call crv.cellranger_vdj as cellranger_vdj {
                input:
                    sample_id = sample_id,
                    input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                    output_directory = output_directory_stripped,
                    genome = generate_count_config.sample2genome[sample_id],
                    acronym_file = acronym_file,
                    denovo = vdj_denovo,
                    chain = vdj_chain,
                    inner_enrichment_primers = vdj_inner_enrichment_primers,
                    cellranger_version = cellranger_version,
                    docker_registry = docker_registry_stripped,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = vdj_disk_space,
                    preemptible = preemptible,
                    backend = backend,
                    awsQueueArn = awsQueueArn
            }
        }

        call collect_summaries as collect_summaries_vdj {
            input:
                summaries = cellranger_vdj.output_metrics_summary,
                sample_ids = cellranger_vdj.output_vdj_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }
    }

    if (length(generate_count_config.sample_feature_ids) > 0) {
        scatter (sample_id in generate_count_config.sample_feature_ids) {
            call ca.cumulus_adt as cumulus_adt {
                input:
                    sample_id = sample_id,
                    input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                    output_directory = output_directory_stripped,
                    chemistry = generate_count_config.sample2chemistry[sample_id],
                    data_type = generate_count_config.sample2datatype[sample_id],
                    feature_barcode_file = generate_count_config.sample2fbf[sample_id],
                    crispr_barcode_pos = crispr_barcode_pos,
                    scaffold_sequence = scaffold_sequence,
                    max_mismatch = max_mismatch,
                    min_read_ratio = min_read_ratio,
                    cumulus_feature_barcoding_version = cumulus_feature_barcoding_version,
                    docker_registry = docker_registry_stripped,
                    acronym_file = acronym_file,
                    zones = zones,
                    num_cpu = feature_num_cpu,
                    memory = feature_memory,
                    disk_space = feature_disk_space,
                    preemptible = preemptible,
                    backend = backend,
                    awsQueueArn = awsQueueArn
            }
        }
    }

    if (length(generate_count_config.sample_atac_ids) > 0) {
        scatter (sample_id in generate_count_config.sample_atac_ids) {
            call crac.cellranger_atac_count as cellranger_atac_count {
                input:
                    sample_id = sample_id,
                    input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                    output_directory = output_directory_stripped,
                    genome = generate_count_config.sample2genome[sample_id],
                    acronym_file = acronym_file,
                    force_cells = force_cells,
                    dim_reduce = atac_dim_reduce,
                    peaks = peaks,
                    chemistry = generate_count_config.sample2chemistry[sample_id],
                    cellranger_atac_version = cellranger_atac_version,
                    docker_registry = docker_registry_stripped,
                    zones = zones,
                    num_cpu = atac_num_cpu,
                    memory = atac_memory,
                    disk_space = atac_disk_space,
                    preemptible = preemptible,
                    backend = backend,
                    awsQueueArn = awsQueueArn
            }
        }

        call collect_summaries as collect_summaries_atac {
            input:
                summaries = cellranger_atac_count.output_metrics_summary,
                sample_ids = cellranger_atac_count.output_count_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }
    }

    if (length(generate_count_config.link_arc_ids) > 0) {
        scatter (link_id in generate_count_config.link_arc_ids) {
            call crarc.cellranger_arc_count as cellranger_arc_count {
                input:
                    link_id = link_id,
                    input_samples = generate_count_config.link2sample[link_id],
                    input_fastqs_directories = generate_count_config.sample2dir[link_id],
                    input_data_types = generate_count_config.sample2datatype[link_id],
                    output_directory = output_directory_stripped,
                    acronym_file = acronym_file,
                    genome = generate_count_config.sample2genome[link_id],
                    gex_exclude_introns = arc_gex_exclude_introns,
                    no_bam = no_bam,
                    min_atac_count = arc_min_atac_count,
                    min_gex_count = arc_min_gex_count,
                    peaks = peaks,
                    cellranger_arc_version = cellranger_arc_version,
                    docker_registry = docker_registry_stripped,
                    zones = zones,
                    num_cpu = arc_num_cpu,
                    memory = arc_memory,
                    disk_space = arc_disk_space,
                    preemptible = preemptible,
                    backend = backend,
                    awsQueueArn = awsQueueArn
            }
        }

        call collect_summaries as collect_summaries_arc {
            input:
                summaries = cellranger_arc_count.output_metrics_summary,
                sample_ids = cellranger_arc_count.output_count_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }
    }

    if (length(generate_count_config.link_multi_ids) > 0) {
        scatter (link_id in generate_count_config.link_multi_ids) {
            call crmulti.cellranger_multi as cellranger_multi {
                input:
                    link_id = link_id,
                    input_samples = generate_count_config.link2sample[link_id],
                    input_fastqs_directories = generate_count_config.sample2dir[link_id],
                    input_data_types = generate_count_config.sample2datatype[link_id],
                    input_fbf = generate_count_config.sample2fbf[link_id],
                    output_directory = output_directory_stripped,
                    acronym_file = acronym_file,
                    genome = generate_count_config.sample2genome[link_id],
                    probe_set = generate_count_config.sample2probeset[link_id],
                    cmo_set = cmo_set,
                    include_introns = include_introns,
                    no_bam = no_bam,
                    secondary = secondary,
                    force_cells = force_cells,
                    expect_cells = expect_cells,
                    cellranger_version = cellranger_version,
                    docker_registry = docker_registry_stripped,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = multi_disk_space,
                    preemptible = preemptible,
                    backend = backend,
                    awsQueueArn = awsQueueArn
            }
        }
    }

    if (length(generate_count_config.link_fbc_ids) > 0) {
        scatter (link_id in generate_count_config.link_fbc_ids) {
            call crc.cellranger_count as cellranger_count_fbc {
                input:
                    sample_id = link_id,
                    input_samples = generate_count_config.link2sample[link_id],
                    input_fastqs_directories = generate_count_config.sample2dir[link_id],
                    input_data_types = generate_count_config.sample2datatype[link_id],
                    input_fbf = generate_count_config.sample2fbf[link_id],
                    output_directory = output_directory_stripped,
                    acronym_file = acronym_file,
                    genome = generate_count_config.sample2genome[link_id],
                    include_introns = include_introns,
                    no_bam = no_bam,
                    secondary = secondary,
                    force_cells = force_cells,
                    expect_cells = expect_cells,
                    cellranger_version = cellranger_version,
                    docker_registry = docker_registry_stripped,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = count_disk_space,
                    preemptible = preemptible,
                    backend = backend,
                    awsQueueArn = awsQueueArn
            }
        }

        call collect_summaries as collect_summaries_fbc {
            input:
                summaries = cellranger_count_fbc.output_metrics_summary,
                sample_ids = cellranger_count_fbc.output_count_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }
    }

    output {
        Map[String, Array[String]?] count_outputs = {
            "gex": cellranger_count.output_count_directory,
            "vdj": cellranger_vdj.output_vdj_directory,
            "adt": cumulus_adt.output_count_directory,
            "atac": cellranger_atac_count.output_count_directory,
            "arc": cellranger_arc_count.output_count_directory,
            "multi": cellranger_multi.output_multi_directory,
            "fbc": cellranger_count_fbc.output_count_directory
        }
        File count_matrix = generate_count_config.count_matrix
    }
}


task generate_count_config {
    input {
        File input_csv_file
        String output_dir
        Array[String]? fastq_dirs
        Array[String]? fastq_dirs_atac
        Array[String]? fastq_dirs_arc
        String config_version
        String docker_registry
        String zones
        Int preemptible
        String awsQueueArn
        String backend
        String null_file
    }

    command {
        set -e
        export TMPDIR=/tmp

        python /software/check_uri.py "~{backend}" "~{output_dir}"

        python <<CODE
        import os
        import re
        import sys
        import pandas as pd
        from collections import defaultdict

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False)
        df.columns = df.columns.str.strip()

        if 'DataType' not in df.columns:
            df['DataType'] = 'rna'
        else:
            df.loc[df['DataType'].isna(), 'DataType'] = 'rna'

        if 'Chemistry' not in df.columns:
            df['Chemistry'] = 'auto'
        else:
            df.loc[df['Chemistry'].isna(), 'Chemistry'] = 'auto'

        for c in df.columns:
            df[c] = df[c].str.strip()

        for idx, row in df.iterrows():
            row['Flowcell'] = re.sub('/+$', '', row['Flowcell'])
            if row['DataType'] not in ['rna', 'vdj', 'adt', 'citeseq', 'cmo', 'crispr', 'atac', 'hashing', 'frp']:
                print("Unknown DataType " + row['DataType'] + " is detected!", file = sys.stderr)
                sys.exit(1)
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)

        with open('sample_ids.txt', 'w') as fo1, open('sample_vdj_ids.txt', 'w') as fo2, open('sample_feature_ids.txt', 'w') as fo3, open('sample_atac_ids.txt', 'w') as fo4, \
             open('sample2dir.txt', 'w') as foo1, open('sample2datatype.txt', 'w') as foo2, open('sample2genome.txt', 'w') as foo3, \
             open('sample2chemistry.txt', 'w') as foo4, open('sample2fbf.txt', 'w') as foo5, open('count_matrix.csv', 'w') as foo6, \
             open('link_arc_ids.txt', 'w') as fo5, open('link_multi_ids.txt', 'w') as fo6, open('link_fbc_ids.txt', 'w') as fo7, \
             open('link2sample.txt', 'w') as foo7, open('sample2probeset.txt', 'w') as fo_s2probeset:

            n_ref = n_chem = n_fbf = n_link = n_probeset = 0 # this mappings can be empty
            foo6.write('Sample,Location,Bam,BamIndex,Barcodes,Reference,Chemistry\n') # count_matrix.csv
            datatype2fo = dict([('rna', fo1), ('vdj', fo2), ('adt', fo3), ('citeseq', fo3), ('hashing', fo3), ('cmo', fo3), ('crispr', fo3), ('atac', fo4)])

            multiomics = defaultdict(set)
            link2sample = defaultdict(list)
            link2dir = defaultdict(list)
            link2dt = defaultdict(list)
            link2fbf = defaultdict(list)
            link2ref = defaultdict(set)
            link2probset = defaultdict(set)

            for sample_id in df['Sample'].unique():
                df_local = df.loc[df['Sample'] == sample_id]
                if df_local['DataType'].unique().size > 1:
                    print('Detected multiple DataType values for sample ' + sample_id + '!', file = sys.stderr)
                    sys.exit(1)

                datatype = df_local['DataType'].iat[0]

                dirs = df_local['Flowcell'].values

                reference = 'null'
                if datatype in ['rna', 'vdj', 'atac', 'frp']:
                    if df_local['Reference'].unique().size > 1:
                        print("Detected multiple references for sample " + sample_id + "!", file = sys.stderr)
                        sys.exit(1)
                    reference = df_local['Reference'].iat[0]

                probeset = 'null'
                if datatype == 'frp':
                    if 'ProbeSet' not in df_local.columns:
                        probeset = 'FRP_human_probe_v1.0.1'
                    else:
                        if df_local['ProbeSet'].unique().size > 1:
                            print("Detected multiple probe sets for sample " + sample_id + "!", file = sys.stderr)
                            sys.exit(1)
                        probeset = df_local['ProbeSet'].iat[0]

                feature_barcode_file = 'null'
                if datatype in ['rna', 'adt', 'citeseq', 'hashing', 'cmo', 'crispr', 'frp']:
                    has_fbf = ('FeatureBarcodeFile' in df_local.columns) and isinstance(df_local['FeatureBarcodeFile'].iat[0], str) and (df_local['FeatureBarcodeFile'].iat[0] != '')
                    if has_fbf:
                        if df_local['FeatureBarcodeFile'].unique().size > 1:
                            print("Detected multiple feature barcode or target panel files for sample " + sample_id + "!", file = sys.stderr)
                            sys.exit(1)
                        feature_barcode_file = df_local['FeatureBarcodeFile'].iat[0]
                    else:
                        if datatype in ['adt', 'citeseq', 'hashing', 'cmo', 'crispr']:
                            print("Please specify one feature barcode file for sample " + sample_id + "!", file = sys.stderr)
                            sys.exit(1)
                        feature_barcode_file = '~{null_file}'

                if ('Link' in df_local) or (datatype == 'frp'): # if multiomics or FRP (Fixed RNA Profiling)
                    if 'Link' in df_local:
                        if df_local['Link'].unique().size > 1:
                            print("Detected multiple Link values for sample '" + sample_id + "'!", file = sys.stderr)
                            sys.exit(1)
                        link = df_local['Link'].iat[0]
                    else:
                        link = sample_id

                    if pd.notnull(link) and (link != ''):
                        multiomics[link].add(datatype)
                        size = dirs.size
                        link2sample[link].extend([sample_id] * size)
                        link2dir[link].extend(list(dirs))
                        link2dt[link].extend([datatype] * size)
                        if feature_barcode_file == '~{null_file}':
                            feature_barcode_file = 'null'
                        link2fbf[link].extend([feature_barcode_file] * size)
                        if reference != 'null':
                            link2ref[link].add(reference)
                        if probeset != 'null':
                            link2probset[link].add(probeset)
                        continue

                datatype2fo[datatype].write(sample_id + '\n')

                foo1.write(sample_id + '\t' + ','.join(dirs) + '\n')
                foo2.write(sample_id + '\t' + datatype + '\n')

                if reference != 'null':
                    foo3.write(sample_id + '\t' + reference + '\n')
                    n_ref += 1

                if feature_barcode_file != 'null':
                    foo5.write(sample_id + '\t' + feature_barcode_file + '\n')
                    n_fbf += 1

                if datatype in ['rna', 'adt', 'citeseq', 'hashing', 'cmo', 'crispr', 'atac']:
                    if df_local['Chemistry'].unique().size > 1:
                        print("Detected multiple chemistry strings for sample " + sample_id + "!", file = sys.stderr)
                        sys.exit(1)
                    chemistry = df_local['Chemistry'].iat[0]
                    if (chemistry in ['auto', 'threeprime']) and (datatype in ['adt', 'citeseq', 'hashing', 'cmo', 'crispr']):
                        chemistry = 'SC3Pv3' # default is different
                    foo4.write(sample_id + '\t' + chemistry + '\n')
                    n_chem += 1

                if datatype == 'rna':
                    prefix = '~{output_dir}/' + sample_id
                    bam = prefix + '/possorted_genome_bam.bam'
                    bai = prefix + '/possorted_genome_bam.bam.bai'
                    count_matrix = prefix + '/filtered_feature_bc_matrix.h5' # assume cellranger version >= 3.0.0
                    barcodes = prefix + '/filtered_feature_bc_matrix/barcodes.tsv.gz'
                    foo6.write(sample_id + ',' + count_matrix + ',' + bam + ',' + bai + ',' + barcodes + ',' + reference + ',' + chemistry + '\n')

            for link_id in multiomics.keys():
                n_link += 1

                ref_set = link2ref.get(link_id, set())
                if len(ref_set) > 1:
                    print("Link '" + link_id + "' contains multiple references!", file = sys.stderr)
                    sys.exit(1)
                if len(ref_set) == 1:
                    n_ref += 1
                    foo3.write(link_id + '\t' + list(ref_set)[0] + '\n')

                probeset_set = link2probset.get(link_id, set())
                if len(probeset_set) > 1:
                    print("Link '" + link_id + "' contains multiple probe sets!", file = sys.stderr)
                    sys.exit(1)
                if len(probeset_set) == 1:
                    n_probeset += 1
                    fo_s2probeset.write(link_id + '\t' + list(probeset_set)[0] + '\n')

                foo1.write(link_id + '\t' + ','.join(link2dir[link_id]) + '\n')
                foo2.write(link_id + '\t' + ','.join(link2dt[link_id]) + '\n')
                foo7.write(link_id + '\t' + ','.join(link2sample[link_id]) + '\n')

                if 'atac' in multiomics[link_id]:
                    if multiomics[link_id] != set(['atac', 'rna']):
                        print("CellRanger ARC only works with ATAC+RNA data! Link '" + link_id + "' contains " + ', '.join(list(multiomics[link_id])) + '.', file = sys.stderr)
                        sys.exit(1)
                    fo5.write(link_id + '\n')
                elif 'cmo' in multiomics[link_id]:
                    if 'rna' not in multiomics[link_id]:
                        print("CellRanger multi expect RNA modality!", file = sys.stderr)
                        sys.exit(1)
                    if not multiomics[link_id].issubset(set(['rna', 'cmo', 'crispr', 'citeseq'])):
                        print("CellRanger multi only works with RNA/CMO/CRISPR/CITESEQ data! Link '" + link_id + "' contains " + ', '.join(list(multiomics[link_id])) + '.', file = sys.stderr)
                        sys.exit(1)
                    fo6.write(link_id + '\n')
                    foo5.write(link_id + '\t' + ','.join(link2fbf[link_id]) + '\n')
                    n_fbf += 1
                elif 'frp' in multiomics[link_id]:
                    fo6.write(link_id + '\n')
                    foo5.write(link_id + '\t' + ','.join(link2fbf[link_id]) + '\n')
                    n_fbf += 1
                else:
                    if not multiomics[link_id].issubset(set(['rna', 'crispr', 'citeseq', 'hashing'])):
                        print("CellRanger count only works with RNA/CRISPR/CITESEQ/HASHING data! Link '" + link_id + "' contains " + ', '.join(list(multiomics[link_id])) + '.', file = sys.stderr)
                        sys.exit(1)
                    fo7.write(link_id + '\n')
                    foo5.write(link_id + '\t' + ','.join(link2fbf[link_id]) + '\n')
                    n_fbf += 1

            if n_ref == 0:
                foo3.write('null\tnull\n')
            if n_probeset == 0:
                fo_s2probeset.write('null\tnull\n')
            if n_chem == 0:
                foo4.write('null\tnull\n')
            if n_fbf == 0:
                foo5.write('null\tnull\n')
            if n_link == 0:
                foo7.write('null\tnull\n')
        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Array[String] sample_vdj_ids = read_lines('sample_vdj_ids.txt')
        Array[String] sample_feature_ids = read_lines('sample_feature_ids.txt')
        Array[String] sample_atac_ids = read_lines('sample_atac_ids.txt')
        Array[String] link_arc_ids = read_lines('link_arc_ids.txt')
        Array[String] link_multi_ids = read_lines('link_multi_ids.txt')
        Array[String] link_fbc_ids = read_lines('link_fbc_ids.txt')
        Map[String, String] sample2dir = read_map('sample2dir.txt')
        Map[String, String] sample2datatype = read_map('sample2datatype.txt')
        Map[String, String] sample2genome = read_map('sample2genome.txt')
        Map[String, String] sample2probeset = read_map('sample2probeset.txt')
        Map[String, String] sample2chemistry = read_map('sample2chemistry.txt')
        Map[String, String] sample2fbf = read_map('sample2fbf.txt')
        Map[String, String] link2sample = read_map('link2sample.txt')
        File count_matrix = "count_matrix.csv"
    }

    runtime {
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}

task collect_summaries {
    input {
        Array[File] summaries
        Array[String] sample_ids
        String config_version
        String docker_registry
        String zones
        Int preemptible
        String awsQueueArn
        String backend
    }

    command {
        python <<CODE
        import pandas as pd
        import os
        import xlsxwriter
        summaries = pd.read_csv('~{write_lines(summaries)}', header = None)
        sample_ids = pd.read_csv('~{write_lines(sample_ids)}', header = None).applymap(lambda x: os.path.basename(x))
        df_info = pd.concat([summaries, sample_ids], axis = 1)
        df_info.columns = ['summary', 'sample_id']
        dfs = []
        for idx, row in df_info.iterrows():
            df = pd.read_csv(row['summary'], header = 0)
            df.index = [row['sample_id']]
            dfs.append(df)
        df_res = pd.concat(dfs)
        df_res.index.name = "Sample"
        writer = pd.ExcelWriter('summaries.xlsx', engine = 'xlsxwriter')
        df_res.to_excel(writer, sheet_name = "summaries")
        writer.save()
        CODE
    }

    output {
        File metrics_summaries = "summaries.xlsx"
    }

    runtime {
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
