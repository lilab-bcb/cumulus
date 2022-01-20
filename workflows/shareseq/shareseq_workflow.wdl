version 1.0

import "shareseq_mkfastq.wdl" as shm
import "shareseq_reorg.wdl" as shr
import "../starsolo/starsolo_count.wdl" as sc
import "../chromap/chromap.wdl" as cm

workflow shareseq_workflow {
    input {
        # 4 - 6 columns (Sample, Lane, Index, Flowcell, [optional for mkfastq: Reference, Type]). gs URL
        File input_csv_file
        # Output directory, gs URL
        String output_directory
        # If run mkfastq
        Boolean run_mkfastq = true
        # If run count
        Boolean run_count = true

        # for mkfastq

        # Dual-index Paired-end flowcell workflow type, choosing from 'auto', 'forward' and 'reverse'. 'auto' means automatically determine the workflow type.
        String dual_index_type = "auto"
        # Whether to delete input_bcl_directory, default: false
        Boolean delete_input_bcl_directory = false
        # Number of allowed mismatches per index
        Int? barcode_mismatches
        # Override the read lengths as specified in RunInfo.xml
        String? use_bases_mask

        # for fastq patterns - shareseq reorg

        # read1 fastq pattern
        String read1_fastq_pattern = "_S*_L*_R1_001.fastq.gz"
        # read2 fastq pattern
        String read2_fastq_pattern = "_S*_L*_R3_001.fastq.gz"
        # index fastq pattern
        String index_fastq_pattern = "_S*_L*_R2_001.fastq.gz"        

        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry = "quay.io/cumulus"
        # shareseq mkfastq registry, default to gcr.io/broad-cumulus
        String mkfastq_docker_registry = "gcr.io/broad-cumulus"
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Backend
        String backend = "gcp"


        # SHARE-Seq index TSV
        File acronym_file = "gs://gred-cumulus-output/ref-data/shareseq/index.tsv"

        # 0.2
        String config_version = "0.2"
        # 0.1.0
        String shareseq_mkfastq_version = "0.1.0"
        # 0.1.0
        String shareseq_reorg_version = "0.1.0"
        # 2.7.9a
        String star_version = "2.7.9a"
        # 0.1.4
        String chromap_version = "0.1.4"

        # Number of cpus shareseq_mkfastq 
        Int shareseq_mkfastq_num_cpu = 32
        # Memory string
        String shareseq_mkfastq_memory = "120G"

        # num cpu sharseq_reorg
        Int sharseq_reorg_num_cpu = 4
        # memory sharseq_reorg
        String shareseq_reorg_memory = "8G"

        # num cpu starsolo
        Int starsolo_num_cpu = 32
        # memory starsolo
        String starsolo_memory = "120G"

        # num cpu chromap
        Int chromap_num_cpu = 8
        # memory chromap
        String chromap_memory = "64G"

        # Optional disk space for mkfastq.
        Int mkfastq_disk_space = 1500
        # Optional disk space for shareseq_reorg.
        Int shareseq_reorg_disk_space = 500
        # Optional disk space for starsolo.
        Int starsolo_disk_space = 500
        # Optional disk space for chromap.
        Int chromap_disk_space = 200

        # Number of preemptible tries
        Int preemptible = 2
        # Max number of retries for AWS instance
        Int awsMaxRetries = 5
                     
    }
    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    String docker_registry_stripped = sub(docker_registry, "/+$", "")
    String mkfastq_docker_registry_stripped = sub(mkfastq_docker_registry, "/+$", "")

    Map[String, String] acronym2gsurl = read_map(acronym_file)

    if (run_mkfastq) {
        call generate_mkfastq_input {
            input:
                input_csv_file = input_csv_file,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsMaxRetries = awsMaxRetries,
                backend = backend   
        }

        if (length(generate_mkfastq_input.bcl_csv) > 0) {
            scatter (bcl_csv in generate_mkfastq_input.bcl_csv) {
                String bcl_key = basename(bcl_csv)
                call shm.shareseq_mkfastq as shareseq_mkfastq {
                    input:
                        input_bcl_directory = generate_mkfastq_input.inpdirs[bcl_key],
                        input_csv_file = bcl_csv,
                        output_directory = output_directory_stripped,
                        delete_input_bcl_directory = delete_input_bcl_directory,
                        dual_index_type = dual_index_type,
                        barcode_mismatches = barcode_mismatches,
                        use_bases_mask = use_bases_mask,
                        shareseq_mkfastq_version = shareseq_mkfastq_version,
                        docker_registry = mkfastq_docker_registry_stripped,
                        zones = zones,
                        num_cpu = shareseq_mkfastq_num_cpu,
                        memory = shareseq_mkfastq_memory,
                        disk_space = mkfastq_disk_space,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend
                }
            }
        }
    }

    if (run_count) {
        call generate_count_config {
            input:
                input_csv_file = input_csv_file,
                fastq_dirs = shareseq_mkfastq.output_fastqs_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsMaxRetries = awsMaxRetries,
                backend = backend
        }

        if (length(generate_count_config.sample_gex_ids) > 0) {
            scatter (sample_id in generate_count_config.sample_gex_ids) {
                call shr.shareseq_reorg as shareseq_reorg_gex {
                    input:
                        sample_id = sample_id,
                        type = 'gex',
                        input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                        r1_fastq_pattern = read1_fastq_pattern,
                        r2_fastq_pattern = read2_fastq_pattern,
                        index_fastq_pattern = index_fastq_pattern,
                        output_directory = output_directory_stripped,
                        shareseq_reorg_version = shareseq_reorg_version,
                        docker_registry = docker_registry,
                        zones = zones,
                        num_cpu = sharseq_reorg_num_cpu,
                        memory = shareseq_reorg_memory,
                        disk_space = shareseq_reorg_disk_space,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend
                }    
                call sc.starsolo_count as starsolo_count {
                    input:
                        sample_id = sample_id,
                        input_fastqs_directories = shareseq_reorg_gex.output_reorg_directory,
                        read1_fastq_pattern = '_f*_R1.fastq.gz',
                        read2_fastq_pattern = '_f*_R2.fastq.gz',
                        assay = 'ShareSeq',
                        barcode_read = 'read2',
                        genome = generate_count_config.sample2genome[sample_id],
                        acronym_file = acronym2gsurl['starsolo'],
                        output_directory = output_directory_stripped,
                        docker_registry = docker_registry,
                        star_version = star_version,
                        zones = zones,
                        memory = starsolo_memory,
                        num_cpu = starsolo_num_cpu,
                        disk_space = starsolo_disk_space,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend
                }
            }
        }

        if (length(generate_count_config.sample_atac_ids) > 0) {
            scatter (sample_id in generate_count_config.sample_atac_ids) {
                call shr.shareseq_reorg as shareseq_reorg_atac {
                    input:
                        sample_id = sample_id,
                        type = 'atac',
                        input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                        r1_fastq_pattern = read1_fastq_pattern,
                        r2_fastq_pattern = read2_fastq_pattern,
                        index_fastq_pattern = index_fastq_pattern,
                        output_directory = output_directory_stripped,
                        shareseq_reorg_version = shareseq_reorg_version,
                        docker_registry = docker_registry,
                        zones = zones,
                        num_cpu = sharseq_reorg_num_cpu,
                        memory = shareseq_reorg_memory,
                        disk_space = shareseq_reorg_disk_space,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend   
                }    
                call cm.chromap as chromap {
                    input:
                        chromap_version = chromap_version,
                        read1_fastq_pattern = '_f*_R1.fastq.gz',
                        read2_fastq_pattern = '_f*_R2.fastq.gz',
                        barcode_fastq_pattern = '_f*_I1.fastq.gz',
                        preset = 'atac',
                        acronym_file = acronym2gsurl['chromap'],
                        sample_id = sample_id,
                        output_directory = output_directory_stripped,
                        input_fastqs_directories = shareseq_reorg_atac.output_reorg_directory,
                        genome = generate_count_config.sample2genome[sample_id],
                        disk_space = chromap_disk_space,
                        docker_registry = docker_registry_stripped,
                        zones = zones,
                        num_cpu = chromap_num_cpu,
                        memory = chromap_memory,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend
                }
            }
        }
    }    
}

task generate_mkfastq_input {
    input {
        File input_csv_file
        String config_version
        String docker_registry
        String zones
        Int preemptible
        Int awsMaxRetries
        String backend
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

        for c in df.columns:
            df[c] = df[c].str.strip()

        mkfastq_cols = ["Lane","Sample","Index","Flowcell"]
        for c in mkfastq_cols: 
            if c not in df.columns:
                print("To run SHARE-seq mkfastq, following column is required: " + c + ". Please correct sample sheet.", file = sys.stderr)
                sys.exit(1)

        for idx, row in df.iterrows():
            row['Flowcell'] = re.sub('/+$', '', row['Flowcell'])
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)

        with open('inpdirs.txt', 'w') as fo:
            for input_dir in df['Flowcell'].unique():
                run_id = os.path.basename(input_dir)
                flowcell_df = df.loc[df['Flowcell'] == input_dir]
                bcl_file = run_id + '_bcl.csv'
                header = ["Lane","Sample","Index"]
                flowcell_df.to_csv(bcl_file, index = False,columns = header)
                fo.write(bcl_file + '\t' + input_dir + '\n')
        CODE
    }

    output {
        Map[String, String] inpdirs = read_map('inpdirs.txt')
        Array[File] bcl_csv = glob('*_bcl.csv')
    }

    runtime {
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}

task generate_count_config {
    input {
        File input_csv_file
        Array[String]? fastq_dirs
        String config_version
        String docker_registry
        String zones
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import os
        import re
        import sys
        import pandas as pd
        from collections import defaultdict

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False)

        for c in df.columns:
            df[c] = df[c].str.strip()

        count_cols = ["Sample","Flowcell","Type","Reference"]
        for c in count_cols: 
            if c not in df.columns:
                print("To run SHARE-seq count (ATAC/GEX), following column is required: " + c + ". Please correct sample sheet.", file = sys.stderr)
                sys.exit(1)

        for idx, row in df.iterrows():
            row['Flowcell'] = re.sub('/+$', '', row['Flowcell'])
            if row['Type'] not in ['gex', 'atac']:
                print("Unknown Type " + row['Type'] + " is detected!", file = sys.stderr)
                sys.exit(1)
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)
            if pd.isnull(row['Reference']):
                print('Reference cannot be empty. No reference provided for sample ' + row['Sample'] + ".")
                sys.exit(1)

        def parse_fastq_dirs(dirs_str):
            r2f = dict()
            if dirs_str == '':
                return r2f
            dirs = dirs_str.split(',')
            for dir in dirs:
                run_id = dir.split('/')[-2].rpartition('_')[0]
                print(run_id)
                r2f[run_id] = dir
            return(r2f)

        r2f = parse_fastq_dirs('~{sep="," fastq_dirs}')

        with open('sample_gex_ids.txt', 'w') as fo1, open('sample_atac_ids.txt', 'w') as fo2, \
             open('sample2dir.txt', 'w') as foo1, open('sample2genome.txt', 'w') as foo2:

            for sample_id in df['Sample'].unique():
                df_local = df.loc[df['Sample'] == sample_id]
                if df_local['Type'].unique().size > 1:
                    print('Detected multiple Type values for sample ' + sample_id + '!', file = sys.stderr)
                    sys.exit(1)

                datatype = df_local['Type'].iat[0]

                if len(r2f) > 0:
                    dirs = df_local['Flowcell'].map(lambda x: r2f[os.path.basename(x)]).values # if also run mkfastq
                else:
                    dirs = df_local['Flowcell'].values # if start from count step

                if datatype in ['gex', 'atac']:
                    if df_local['Reference'].unique().size > 1:
                        print("Detected multiple references for sample " + sample_id + "!", file = sys.stderr)
                        sys.exit(1)
                    reference = df_local['Reference'].iat[0]

                foo1.write(sample_id + '\t' + ','.join(dirs) + '\n')
                foo2.write(sample_id + '\t' + reference + '\n')
                
        CODE
    }

    output {
        Array[String] sample_gex_ids = read_lines('sample_gex_ids.txt')
        Array[String] sample_atac_ids = read_lines('sample_atac_ids.txt')
        Map[String, String] sample2dir = read_map('sample2dir.txt')
        Map[String, String] sample2genome = read_map('sample2genome.txt')
    }

    runtime {
        docker: "~{docker_registry}/config:~{config_version}"
        zones: zones
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}