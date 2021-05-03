version 1.0

import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_mkfastq/versions/10/plain-WDL/descriptor" as crm
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_count/versions/7/plain-WDL/descriptor" as crc
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_vdj/versions/8/plain-WDL/descriptor" as crv
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cumulus_adt/versions/8/plain-WDL/descriptor" as ca
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_atac_mkfastq/versions/5/plain-WDL/descriptor" as cram
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_atac_count/versions/6/plain-WDL/descriptor" as crac

workflow cellranger_workflow {
    input {
        # 5 - 8 columns (Sample, Reference, Flowcell, Lane, Index, [Chemistry, DataType, FeatureBarcodeFile]). gs URL
        File input_csv_file
        # Output directory, gs URL
        String output_directory

        # If run mkfastq
        Boolean run_mkfastq = true
        # If run count
        Boolean run_count = true

        # for mkfastq

        # Whether to delete input_bcl_directory, default: false
        Boolean delete_input_bcl_directory = false
        # Number of allowed mismatches per index
        Int? mkfastq_barcode_mismatches
        # Only demultiplex samples identified by an i7-only sample index, ignoring dual-indexed samples.  Dual-indexed samples will not be demultiplexed.
        Boolean mkfastq_filter_single_index = false
        # Override the read lengths as specified in RunInfo.xml
        String? mkfastq_use_bases_mask

        # For cellranger count

        # Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. This is also a cellranger-atac argument.
        Int? force_cells
        # Expected number of recovered cells. Mutually exclusive with force_cells
        Int? expect_cells
        # If count reads mapping to intronic regions
        Boolean include_introns = false
        # If generate bam outputs. This is also a spaceranger argument.
        Boolean no_bam = false
        # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false. This is also a spaceranger argument.
        Boolean secondary = false

        # For vdj

        # Do not align reads to reference V(D)J sequences before de novo assembly. Default: false
        Boolean vdj_denovo = false

        # Force the analysis to be carried out for a particular chain type. The accepted values are:
        #   "auto" for auto detection based on TR vs IG representation (default),
        #   "TR" for T cell receptors,
        #   "IG" for B cell receptors,
        # Use this in rare cases when automatic chain detection fails.
        String vdj_chain = "auto"

        # For extracting ADT count

        # scaffold sequence for Perturb-seq, default is "", which for Perturb-seq means barcode starts at position 0 of read 2
        String scaffold_sequence = ""
        # maximum hamming distance in feature barcodes
        Int max_mismatch = 3
        # minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for data type crispr
        Float min_read_ratio = 0.1

        # For atac

        # For atac, choose the algorithm for dimensionality reduction prior to clustering and tsne: 'lsa' (default), 'plsa', or 'pca'.
        String? atac_dim_reduce

        # 6.0.1, 6.0.0, 5.0.1, 5.0.0, 4.0.0, 3.1.0, 3.0.2, 2.2.0
        String cellranger_version = "6.0.1"
        # 0.5.0, 0.4.0, 0.3.0, 0.2.0
        String cumulus_feature_barcoding_version = "0.5.0"
        # 1.2.0, 1.1.0
        String cellranger_atac_version = "1.2.0"
        # 0.2
        String config_version = "0.2"

        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry = "quay.io/cumulus"
        # cellranger/cellranger-atac mkfastq registry, default to gcr.io/broad-cumulus
        String mkfastq_docker_registry = "gcr.io/broad-cumulus"
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Number of cpus per cellranger and spaceranger job
        Int num_cpu = 32
        # Memory string
        String memory = "120G"

        # Number of cpus for cellranger-atac count
        Int atac_num_cpu = 64
        # Memory string for cellranger-atac count
        String atac_memory = "57.6G"

        # Optional memory string for cumulus_adt
        String feature_memory = "32G"

        # Optional disk space for mkfastq.
        Int mkfastq_disk_space = 1500
        # Optional disk space needed for cell ranger count.
        Int count_disk_space = 500
        # Optional disk space needed for cell ranger vdj.
        Int vdj_disk_space = 500
        # Optional disk space needed for cumulus_adt
        Int feature_disk_space = 100
        # Optional disk space needed for cellranger-atac count
        Int atac_disk_space = 500

        # Number of preemptible tries
        Int preemptible = 2
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    String docker_registry_stripped = sub(docker_registry, "/+$", "")
    String mkfastq_docker_registry_stripped = sub(mkfastq_docker_registry, "/+$", "")

    if (run_mkfastq) {
        call generate_bcl_csv {
            input:
                input_csv_file = input_csv_file,
                config_version = config_version,
                zones = zones,
                preemptible = preemptible,
                docker_registry = docker_registry_stripped
        }

        if (length(generate_bcl_csv.bcl_csv_rna) > 0) {
            scatter (bcl_csv in generate_bcl_csv.bcl_csv_rna) {
                String rna_key = basename(bcl_csv)
                call crm.cellranger_mkfastq as cellranger_mkfastq {
                    input:
                        input_bcl_directory = generate_bcl_csv.inpdirs[rna_key],
                        input_csv_file = bcl_csv,
                        output_directory = output_directory_stripped,
                        delete_input_bcl_directory = delete_input_bcl_directory,
                        barcode_mismatches = mkfastq_barcode_mismatches,
                        filter_single_index = mkfastq_filter_single_index,
                        use_bases_mask = mkfastq_use_bases_mask,
                        cellranger_version = cellranger_version,
                        docker_registry = mkfastq_docker_registry_stripped,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        disk_space = mkfastq_disk_space,
                        preemptible = preemptible
                }
            }
        }

        if (length(generate_bcl_csv.bcl_csv_atac) > 0) {
            scatter (bcl_csv in generate_bcl_csv.bcl_csv_atac) {
                String atac_key = basename(bcl_csv)
                call cram.cellranger_atac_mkfastq as cellranger_atac_mkfastq {
                    input:
                        input_bcl_directory = generate_bcl_csv.inpdirs[atac_key],
                        input_csv_file = bcl_csv,
                        output_directory = output_directory_stripped,
                        delete_input_bcl_directory = delete_input_bcl_directory,
                        barcode_mismatches = mkfastq_barcode_mismatches,
                        cellranger_atac_version = cellranger_atac_version,
                        docker_registry = mkfastq_docker_registry_stripped,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        disk_space = mkfastq_disk_space,
                        preemptible = preemptible
                }
            }
        }
    }

    if (run_count) {
        call generate_count_config {
            input:
                input_csv_file = input_csv_file,
                output_dir = output_directory_stripped,
                fastq_dirs = cellranger_mkfastq.output_fastqs_flowcell_directory,
                fastq_dirs_atac = cellranger_atac_mkfastq.output_fastqs_flowcell_directory,
                config_version = config_version,
                zones = zones,
                preemptible = preemptible,
                docker_registry = docker_registry_stripped
        }

        if (generate_count_config.sample_ids[0] != '') {
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
                        no_bam = no_bam,
                        secondary = secondary,
                        force_cells = force_cells,
                        expect_cells = expect_cells,
                        cellranger_version = cellranger_version,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        disk_space = count_disk_space,
                        preemptible = preemptible,
                        docker_registry = docker_registry_stripped
                }
            }

            call collect_summaries {
                input:
                    summaries = cellranger_count.output_metrics_summary,
                    sample_ids = cellranger_count.output_count_directory,
                    config_version = config_version,
                    zones = zones,
                    preemptible = preemptible,
                    docker_registry = docker_registry_stripped
            }
        }

        if (generate_count_config.sample_vdj_ids[0] != '') {
            scatter (sample_id in generate_count_config.sample_vdj_ids) {
                call crv.cellranger_vdj as cellranger_vdj {
                    input:
                        sample_id = sample_id,
                        input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                        output_directory = output_directory_stripped,
                        genome = generate_count_config.sample2genome[sample_id],
                        denovo = vdj_denovo,
                        chain = vdj_chain,
                        cellranger_version = cellranger_version,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        disk_space = vdj_disk_space,
                        preemptible = preemptible,
                        docker_registry = docker_registry_stripped
                }
            }

            call collect_summaries as collect_summaries_vdj {
                input:
                    summaries = cellranger_vdj.output_metrics_summary,
                    sample_ids = cellranger_vdj.output_vdj_directory,
                    config_version = config_version,
                    zones = zones,
                    preemptible = preemptible,
                    docker_registry = docker_registry_stripped
            }
        }

        if (generate_count_config.sample_feature_ids[0] != '') {
            scatter (sample_id in generate_count_config.sample_feature_ids) {
                call ca.cumulus_adt as cumulus_adt {
                    input:
                        sample_id = sample_id,
                        input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                        output_directory = output_directory_stripped,
                        chemistry = generate_count_config.sample2chemistry[sample_id],
                        data_type = generate_count_config.sample2datatype[sample_id],
                        feature_barcode_file = generate_count_config.sample2fbf[sample_id],
                        scaffold_sequence = scaffold_sequence,
                        max_mismatch = max_mismatch,
                        min_read_ratio = min_read_ratio,
                        cumulus_feature_barcoding_version = cumulus_feature_barcoding_version,
                        zones = zones,
                        memory = feature_memory,
                        disk_space = feature_disk_space,
                        preemptible = preemptible,
                        docker_registry = docker_registry_stripped
                }
            }
        }

        if (generate_count_config.sample_atac_ids[0] != '') {
            scatter (sample_id in generate_count_config.sample_atac_ids) {
                call crac.cellranger_atac_count as cellranger_atac_count {
                    input:
                        sample_id = sample_id,
                        input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                        output_directory = output_directory_stripped,
                        genome = generate_count_config.sample2genome[sample_id],
                        force_cells = force_cells,
                        dim_reduce = atac_dim_reduce,
                        cellranger_atac_version = cellranger_atac_version,
                        zones = zones,
                        num_cpu = atac_num_cpu,
                        memory = atac_memory,
                        disk_space = atac_disk_space,
                        preemptible = preemptible,
                        docker_registry = docker_registry_stripped
                }
            }

            call collect_summaries as collect_summaries_atac {
                input:
                    summaries = cellranger_atac_count.output_metrics_summary,
                    sample_ids = cellranger_atac_count.output_count_directory,
                    config_version = config_version,
                    zones = zones,
                    preemptible = preemptible,
                    docker_registry = docker_registry_stripped
            }
        }
    }

    output {
        String? count_matrix = generate_count_config.count_matrix
    }
}

task generate_bcl_csv {
    input {
        File input_csv_file
        String config_version
        String zones
        Int preemptible
        String docker_registry
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import os
        import re
        import sys
        import pandas as pd

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False)

        if 'DataType' not in df.columns:
            df['DataType'] = 'rna'
        else:
            df.loc[df['DataType'].isna(), 'DataType'] = 'rna'

        for c in df.columns:
            df[c] = df[c].str.strip()

        for idx, row in df.iterrows():
            row['Flowcell'] = re.sub('/+$', '', row['Flowcell'])
            if row['DataType'] not in ['rna', 'vdj', 'adt', 'cmo', 'crispr', 'atac']:
                print("Unknown DataType " + row['DataType'] + " is detected!", file = sys.stderr)
                sys.exit(1)
            if row['DataType'] in ['vdj', 'adt', 'cmo', 'crispr']:
                row['DataType'] = 'rna'
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)

        with open('inpdirs.txt', 'w') as fo:
            for input_dir in df['Flowcell'].unique():
                run_id = os.path.basename(input_dir)
                flowcell_df = df.loc[df['Flowcell'] == input_dir]
                for datatype in flowcell_df['DataType'].unique():
                    bcl_df = flowcell_df.loc[flowcell_df['DataType'] == datatype, ['Lane', 'Sample', 'Index']]
                    bcl_file = run_id + '_' + datatype + '_bcl.csv'
                    bcl_df.to_csv(bcl_file, index = False)
                    fo.write(bcl_file + '\t' + input_dir + '\n')
        CODE
    }

    output {
        Map[String, String] inpdirs = read_map('inpdirs.txt')
        Array[File] bcl_csv_rna = glob('*_rna_bcl.csv')
        Array[File] bcl_csv_atac = glob('*_atac_bcl.csv')
    }

    runtime {
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: "~{preemptible}"
    }
}

task generate_count_config {
    input {
        File input_csv_file
        String output_dir
        Array[String]? fastq_dirs
        Array[String]? fastq_dirs_atac
        String config_version
        String zones
        Int preemptible
        String docker_registry
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import os
        import re
        import sys
        import pandas as pd

        null_file = 'gs://regev-lab/resources/cellranger/null' # null file

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False)

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
            if row['DataType'] not in ['rna', 'vdj', 'adt', 'cmo', 'crispr', 'atac']:
                print("Unknown DataType " + row['DataType'] + " is detected!", file = sys.stderr)
                sys.exit(1)
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)

        def parse_fastq_dirs(dirs_str):
            r2f = dict()
            if dirs_str == '':
                return r2f
            dirs = dirs_str.split(',')
            for dir in dirs:
                run_id = dir.split('/')[-3].rpartition('_')[0]
                r2f[run_id] = dir
            return r2f

        r2f = parse_fastq_dirs('~{sep="," fastq_dirs}')
        ar2f = parse_fastq_dirs('~{sep="," fastq_dirs_atac}')

        with open('sample_ids.txt', 'w') as fo1, open('sample_vdj_ids.txt', 'w') as fo2, open('sample_feature_ids.txt', 'w') as fo3, open('sample_atac_ids.txt', 'w') as fo4, \
             open('sample2dir.txt', 'w') as foo1, open('sample2datatype.txt', 'w') as foo2, open('sample2genome.txt', 'w') as foo3, \
             open('sample2chemistry.txt', 'w') as foo4, open('sample2fbf.txt', 'w') as foo5, open('count_matrix.csv', 'w') as foo6:

            n_ref = n_chem = n_fbf = 0 # this mappings can be empty
            foo6.write('Sample,Location,Bam,BamIndex,Barcodes,Reference,Chemistry\n') # count_matrix.csv
            datatype2fo = dict([('rna', fo1), ('vdj', fo2), ('adt', fo3), ('cmo', fo3) ('crispr', fo3), ('atac', fo4)])
            datatype2r2f = dict([('rna', r2f), ('vdj', r2f), ('adt', r2f), ('cmo', r2f), ('crispr', r2f), ('atac', ar2f)])
            for sample_id in df['Sample'].unique():
                df_local = df.loc[df['Sample'] == sample_id]
                if df_local['DataType'].unique().size > 1:
                    print('Detected multiple DataType values for sample ' + sample_id + '!', file = sys.stderr)
                    sys.exit(1)

                datatype = df_local['DataType'].iat[0]
                datatype2fo[datatype].write(sample_id + '\n')

                rid2fdir = datatype2r2f[datatype]
                if len(rid2fdir) > 0:
                    dirs = df_local['Flowcell'].map(lambda x: rid2fdir[os.path.basename(x)]).values # if also run mkfastq
                else:
                    dirs = df_local['Flowcell'].values # if start from count step
                foo1.write(sample_id + '\t' + ','.join(dirs) + '\n')

                foo2.write(sample_id + '\t' + datatype + '\n')

                if datatype in ['rna', 'vdj', 'atac']:
                    if df_local['Reference'].unique().size > 1:
                        print("Detected multiple references for sample " + sample_id + "!", file = sys.stderr)
                        sys.exit(1)
                    reference = df_local['Reference'].iat[0]
                    foo3.write(sample_id + '\t' + reference + '\n')
                    n_ref += 1

                if datatype in ['rna', 'adt', 'cmo', 'crispr']:
                    if df_local['Chemistry'].unique().size > 1:
                        print("Detected multiple chemistry strings for sample " + sample_id + "!", file = sys.stderr)
                        sys.exit(1)
                    chemistry = df_local['Chemistry'].iat[0]
                    if chemistry == 'auto' and datatype in ['adt', 'cmo', 'crispr']:
                        chemistry = 'SC3Pv3' # default is different
                    foo4.write(sample_id + '\t' + chemistry + '\n')
                    n_chem += 1

                if datatype in ['rna', 'adt', 'cmo', 'crispr']:
                    has_fbf = ('FeatureBarcodeFile' in df_local.columns) and isinstance(df_local['FeatureBarcodeFile'].iat[0], str) and (df_local['FeatureBarcodeFile'].iat[0] != '')
                    if has_fbf:
                        if df_local['FeatureBarcodeFile'].unique().size > 1:
                            print("Detected multiple feature barcode or target panel files for sample " + sample_id + "!", file = sys.stderr)
                            sys.exit(1)
                        feature_barcode_file = df_local['FeatureBarcodeFile'].iat[0]
                    else:
                        if datatype in ['adt', 'cmo', 'crispr']:
                            print("Please specify one feature barcode file for sample " + sample_id + "!", file = sys.stderr)
                            sys.exit(1)
                        feature_barcode_file = null_file
                    foo5.write(sample_id + '\t' + feature_barcode_file + '\n')
                    n_fbf += 1

                if datatype == 'rna':
                    prefix = '~{output_dir}/' + sample_id
                    bam = prefix + '/possorted_genome_bam.bam'
                    bai = prefix + '/possorted_genome_bam.bam.bai'
                    count_matrix = prefix + '/filtered_feature_bc_matrix.h5' # assume cellranger version >= 3.0.0
                    barcodes = prefix + '/filtered_feature_bc_matrix/barcodes.tsv.gz'
                    foo6.write(sample_id + ',' + count_matrix + ',' + bam + ',' + bai + ',' + barcodes + ',' + reference + ',' + chemistry + '\n')

            if n_ref == 0:
                foo3.write('null\tnull\n')
            if n_chem == 0:
                foo4.write('null\tnull\n')
            if n_fbf == 0:
                foo5.write('null\tnull\n')
        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Array[String] sample_vdj_ids = read_lines('sample_vdj_ids.txt')
        Array[String] sample_feature_ids = read_lines('sample_feature_ids.txt')
        Array[String] sample_atac_ids = read_lines('sample_atac_ids.txt')
        Map[String, String] sample2dir = read_map('sample2dir.txt')
        Map[String, String] sample2datatype = read_map('sample2datatype.txt')
        Map[String, String] sample2genome = read_map('sample2genome.txt')
        Map[String, String] sample2chemistry = read_map('sample2chemistry.txt')
        Map[String, String] sample2fbf = read_map('sample2fbf.txt')
        File count_matrix = "count_matrix.csv"
    }

    runtime {
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: "~{preemptible}"
    }
}

task collect_summaries {
    input {
        Array[File] summaries
        Array[String] sample_ids
        String config_version
        String zones
        Int preemptible
        String docker_registry
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
        preemptible: "~{preemptible}"
    }
}
