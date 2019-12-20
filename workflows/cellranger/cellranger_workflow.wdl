import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_atac_count/versions/3/plain-WDL/descriptor" as crac
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_atac_mkfastq/versions/2/plain-WDL/descriptor" as cram
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_count/versions/3/plain-WDL/descriptor" as crc
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_mkfastq/versions/2/plain-WDL/descriptor" as crm
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cellranger_vdj/versions/3/plain-WDL/descriptor" as crv
import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cumulus_adt/versions/4/plain-WDL/descriptor" as ca

workflow cellranger_workflow {
    # 5 - 8 columns (Sample, Reference, Flowcell, Lane, Index, [Chemistry, DataType, FeatureBarcodeFile]). gs URL
    File input_csv_file
    # Output directory, gs URL
    String output_directory
    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "/+$", "")

    # If run cellranger mkfastq
    Boolean run_mkfastq = true
    # If run cellranger count
    Boolean run_count = true

    # for mkfastq

    # Whether to delete input_bcl_directory, default: false
    Boolean? delete_input_bcl_directory = false
    # Number of allowed mismatches per index
    Int? mkfastq_barcode_mismatches

    # common to cellranger count/vdj and cellranger-atac count

    # Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells.
    Int? force_cells

    # For count

    # Expected number of recovered cells. Mutually exclusive with force_cells
    Int? expect_cells
    # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false
    Boolean? secondary = false

    # For vdj

    # Do not align reads to reference V(D)J sequences before de novo assembly. Default: false
    Boolean? vdj_denovo = false
    # Force the web summary HTML and metrics summary CSV to only report on a particular chain type. The accepted values are: auto for autodetection based on TR vs IG representation, TR for T cell receptors, IG for B cell receptors, all for all chain types.
    String? vdj_chain 

    # For extracting ADT count

    # scaffold sequence for Perturb-seq, set to "" for adt
    String? scaffold_sequence = ""
    # maximum hamming distance in feature barcodes
    Int? max_mismatch = 3
    # minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for data type crispr
    Float? min_read_ratio = 0.1

    # For atac

    # For atac, choose the algorithm for dimensionality reduction prior to clustering and tsne: 'lsa' (default), 'plsa', or 'pca'.
    String? atac_dim_reduce


    # 3.1.0, 3.0.2, 2.2.0 
    String? cellranger_version = "3.1.0"
    # 1.1.0
    String? cellranger_atac_version = "1.1.0"
    # cumulus version, default to "0.12.0"
    String? cumulus_version = "0.12.0"
    String? docker_registry = "cumulusprod"
    String? docker_registry_stripped = sub(docker_registry, "/+$", "")
    String? cellranger_mkfastq_docker_registry = "gcr.io/broad-cumulus"
    String? cellranger_mkfastq_docker_registry_stripped = sub(cellranger_mkfastq_docker_registry, "/+$", "")
    # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    # Number of cpus per cellranger job
    Int? num_cpu = 32
    # Memory string
    String? memory = "120G"

    # Number of cpus for cellranger-atac count
    Int? atac_num_cpu = 64
    # Memory string for cellranger-atac count
    String? atac_memory = "57.6G"

    # Optional memory string for cumulus_adt
    String? feature_memory = "32G"
    # Optional disk space for mkfastq.
    Int? mkfastq_disk_space = 1500
    # Optional disk space needed for cell ranger count.
    Int? count_disk_space = 500    
    # Optional disk space needed for cell ranger vdj.
    Int? vdj_disk_space = 500    
    # Optional disk space needed for cumulus_adt
    Int? feature_disk_space = 100
    # Optional disk space needed for cellranger-atac count
    Int? atac_disk_space = 500
    # Number of preemptible tries 
    Int? preemptible = 2

    if (run_mkfastq) {
        call generate_bcl_csv {
            input:
                input_csv_file = input_csv_file,
                output_dir = output_directory_stripped,
                cellranger_version = cellranger_version,
                zones = zones,
                preemptible = preemptible,
                docker_registry = docker_registry_stripped
        }

        if (generate_bcl_csv.run_ids[0] != '') {
            scatter (run_id in generate_bcl_csv.run_ids) {
                call crm.cellranger_mkfastq as cellranger_mkfastq {
                    input:
                        input_bcl_directory = generate_bcl_csv.inpdirs[run_id],
                        input_csv_file = generate_bcl_csv.bcls[run_id],
                        output_directory = output_directory_stripped,
                        delete_input_bcl_directory = delete_input_bcl_directory,
                        barcode_mismatches = mkfastq_barcode_mismatches,
                        cellranger_version = cellranger_version,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        docker_registry = cellranger_mkfastq_docker_registry_stripped,
                        disk_space = mkfastq_disk_space,
                        preemptible = preemptible
                }
            }
        }

        if (generate_bcl_csv.run_ids_atac[0] != '') {
            scatter (run_id in generate_bcl_csv.run_ids_atac) {
                call cram.cellranger_atac_mkfastq as cellranger_atac_mkfastq {
                    input:
                        input_bcl_directory = generate_bcl_csv.inpdirs[run_id],
                        input_csv_file = generate_bcl_csv.bcls[run_id],
                        output_directory = output_directory_stripped,
                        delete_input_bcl_directory = delete_input_bcl_directory,
                        barcode_mismatches = mkfastq_barcode_mismatches,
                        cellranger_atac_version = cellranger_atac_version,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        docker_registry = cellranger_mkfastq_docker_registry_stripped,
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
                run_ids = generate_bcl_csv.run_ids,
                fastq_dirs = cellranger_mkfastq.output_fastqs_flowcell_directory,
                run_ids_atac = generate_bcl_csv.run_ids_atac,
                fastq_dirs_atac = cellranger_atac_mkfastq.output_fastqs_flowcell_directory,
                cellranger_version = cellranger_version,
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
                        chemistry = generate_count_config.sample2chemistry[sample_id],
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
                    cellranger_version = cellranger_version,
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
                        force_cells = force_cells,
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
                    cellranger_version = cellranger_version,
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
                        cumulus_version = cumulus_version,
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
                    cellranger_version = cellranger_version,
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
    File input_csv_file
    String output_dir
    String cellranger_version
    String zones
    Int preemptible
    String docker_registry

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import os
        import re
        import sys
        import pandas as pd 
        from subprocess import check_call

        df = pd.read_csv('${input_csv_file}', header = 0, dtype=str, index_col=False)
        for c in df.columns:
            df[c] = df[c].str.strip()
        df['Flowcell'] = df['Flowcell'].map(lambda x: re.sub('/+$', '', x)) # remove trailing slashes
        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)
        with open('run_ids.txt', 'w') as fo1, open('inpdirs.txt', 'w') as fo2, open('bcls.txt', 'w') as fo3, open('run_ids_atac.txt', 'w') as fo4:
            for input_dir in df['Flowcell'].unique():
                run_id = os.path.basename(input_dir)
                bcl_df = df.loc[df['Flowcell'] == input_dir, ['Lane', 'Sample', 'Index']]
                bcl_df.to_csv(run_id + '_bcl.csv', index = False)
                call_args = ['gsutil', 'cp', run_id + '_bcl.csv', '${output_dir}/']
                # call_args = ['cp', run_id + '_bcl.csv', '${output_dir}/']
                print(' '.join(call_args))
                check_call(call_args)
                if 'DataType' in df.columns and df.loc[df['Flowcell'] == input_dir, 'DataType'].iat[0] == 'atac':
                    fo4.write(run_id + '\n')
                else:
                    fo1.write(run_id + '\n')
                fo2.write(run_id + '\t' + input_dir + '\n')
                fo3.write(run_id + '\t${output_dir}/' + run_id + '_bcl.csv\n')
        CODE
    }

    output {
        Array[String] run_ids = read_lines('run_ids.txt')
        Array[String] run_ids_atac = read_lines('run_ids_atac.txt')
        Map[String, String] inpdirs = read_map('inpdirs.txt')
        Map[String, String] bcls = read_map('bcls.txt')
    }

    runtime {
        docker: "${docker_registry}/cellranger:${cellranger_version}"
        zones: zones
        preemptible: "${preemptible}"
    }
}

task generate_count_config {
    File input_csv_file
    String output_dir
    Array[String]? run_ids
    Array[String]? fastq_dirs
    Array[String]? run_ids_atac
    Array[String]? fastq_dirs_atac
    String cellranger_version
    String zones
    Int preemptible
    String docker_registry

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import os
        import re
        import sys
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv('${input_csv_file}', header = 0, dtype=str)
        for c in df.columns:
            df[c] = df[c].str.strip()

        df['Flowcell'] = df['Flowcell'].map(lambda x: re.sub('/+$', '', x)) # remove trailing slashes
        run_ids = '${sep="," run_ids}'.split(',')
        fastq_dirs = '${sep="," fastq_dirs}'.split(',')

        run_ids.extend('${sep="," run_ids_atac}'.split(','))
        fastq_dirs.extend('${sep="," fastq_dirs_atac}'.split(','))

        rid2fdir = dict()
        for run_id, fastq_dir in zip(run_ids, fastq_dirs):
            if run_id is not '':
                rid2fdir[run_id] = fastq_dir

        with open('sample_ids.txt', 'w') as fo1, open('sample2dir.txt', 'w') as fo2, open('sample2genome.txt', 'w') as fo3, open('sample2chemistry.txt', 'w') as fo4, \
             open('count_matrix.csv', 'w') as fo5, open('sample_vdj_ids.txt', 'w') as fo6, open('sample_feature_ids.txt', 'w') as fo7, \
             open('sample2datatype.txt', 'w') as fo8, open('sample2fbf.txt', 'w') as fo9, open('sample_atac_ids.txt', 'w') as fo10:

            fo5.write('Sample,Location,Bam,BamIndex,Barcodes,Reference,Chemistry\n')

            n_ref = n_chem = n_fbf = 0

            for sample_id in df['Sample'].unique():
                df_local = df.loc[df['Sample'] == sample_id]
                
                data_type = 'rna'
                if 'DataType' in df_local.columns:
                    assert df_local['DataType'].unique().size == 1
                    data_type = df_local['DataType'].iat[0]

                if data_type == 'rna':
                    fo1.write(sample_id + '\n')
                elif data_type == 'vdj':
                    fo6.write(sample_id + '\n')
                elif data_type == 'adt' or data_type == 'crispr':
                    fo7.write(sample_id + '\n')
                elif data_type == 'atac':
                    fo10.write(sample_id + '\n')
                else:
                    print('Invalid data type: ' + data_type + '!')
                    assert False

                dirs = df_local['Flowcell'].map(lambda x: x if len(rid2fdir) == 0 else rid2fdir[os.path.basename(x)]).values
                fo2.write(sample_id + '\t' + ','.join(dirs) + '\n')
                
                if data_type == 'rna' or data_type == 'vdj' or data_type == 'atac':
                    assert df_local['Reference'].unique().size == 1
                    reference = df_local['Reference'].iat[0]
                    fo3.write(sample_id + '\t' + reference + '\n')
                    n_ref += 1

                if data_type == 'rna' or data_type == 'adt' or data_type == 'crispr':
                    chemistry = 'auto'
                    if 'Chemistry' in df_local.columns:
                        assert df_local['Chemistry'].unique().size == 1
                        chemistry = df_local['Chemistry'].iat[0]
                    if chemistry == 'auto' and data_type != 'rna':
                        chemistry = 'SC3Pv3'
                    fo4.write(sample_id + '\t' + chemistry + '\n')
                    n_chem += 1

                if data_type == 'adt' or data_type == 'crispr':
                    assert 'FeatureBarcodeFile' in df_local.columns
                    assert df_local['FeatureBarcodeFile'].unique().size == 1
                    feature_barcode_file = df_local['FeatureBarcodeFile'].iat[0]
                    assert feature_barcode_file != ''
                    fo8.write(sample_id + '\t' + data_type + '\n')
                    fo9.write(sample_id + '\t' + feature_barcode_file + '\n')
                    n_fbf += 1

                if data_type == 'rna':
                    prefix = '${output_dir}/' + sample_id
                    bam = prefix + '/possorted_genome_bam.bam'
                    bai = prefix + '/possorted_genome_bam.bam.bai'
                    if int('${cellranger_version}'.split('.')[0]) >= 3:
                        count_matrix = prefix + '/filtered_feature_bc_matrix.h5'
                        barcodes = prefix + '/filtered_feature_bc_matrix/barcodes.tsv.gz'
                    else:
                        count_matrix = prefix + '/filtered_gene_bc_matrices_h5.h5'
                        barcodes = prefix + '/filtered_gene_bc_matrices/' + os.path.splitext(os.path.basename(reference))[0] + '/barcodes.tsv'
                    fo5.write(sample_id + ',' + count_matrix + ',' + bam + ',' + bai + ',' + barcodes + ',' + reference + ',' + chemistry + '\n')

            if n_ref == 0:
                fo3.write('null\tnull\n')
            if n_chem == 0:
                fo4.write('null\tnull\n')
            if n_fbf == 0:
                fo8.write('null\tnull\n')
                fo9.write('null\tnull\n')
        CODE

        gsutil -m cp count_matrix.csv ${output_dir}/
        #mkdir -p ${output_dir}
        #cp count_matrix.csv ${output_dir}/
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] sample2dir = read_map('sample2dir.txt')
        Map[String, String] sample2genome = read_map('sample2genome.txt')
        Map[String, String] sample2chemistry = read_map('sample2chemistry.txt')
        String count_matrix = "${output_dir}/count_matrix.csv"
        Array[String] sample_vdj_ids = read_lines('sample_vdj_ids.txt')
        Array[String] sample_feature_ids = read_lines('sample_feature_ids.txt')
        Map[String, String] sample2datatype = read_map('sample2datatype.txt')
        Map[String, String] sample2fbf = read_map('sample2fbf.txt')
        Array[String] sample_atac_ids = read_lines('sample_atac_ids.txt')
    }

    runtime {
        docker: "${docker_registry}/cellranger:${cellranger_version}"
        zones: zones
        preemptible: "${preemptible}"
    }
}

task collect_summaries {
    Array[File] summaries
    Array[String] sample_ids
    String cellranger_version
    String zones
    Int preemptible
    String docker_registry

    command {
        python <<CODE
        import pandas as pd
        import os
        import xlsxwriter
        summaries = pd.read_csv('${write_lines(summaries)}', header = None)
        sample_ids = pd.read_csv('${write_lines(sample_ids)}', header = None).applymap(lambda x: os.path.basename(x))
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
        docker: "${docker_registry}/cellranger:${cellranger_version}"
        zones: zones
        preemptible: "${preemptible}"
    }
}
