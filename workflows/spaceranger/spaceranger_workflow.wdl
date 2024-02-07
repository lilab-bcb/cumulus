version 1.0

import "spaceranger_mkfastq.wdl" as srm
import "spaceranger_count.wdl" as src

workflow spaceranger_workflow {
    input {
        # 5 - 15 columns (Sample, Reference, Flowcell, Lane, Index, [ProbeSet, Image, DarkImage, ColorizedImage, CytaImage, Slide, Area, SlideFile, LoupeAlignment, TargetPanel])
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

        # Referece index TSV
        File acronym_file = "gs://regev-lab/resources/cellranger/index.tsv"

        # For spaceranger count

        # Use with automatic image alignment to specify that images may not be in canonical orientation with the hourglass in the top left corner of the image. The automatic fiducial alignment will attempt to align any rotation or mirroring of the image.
        Boolean reorient_images = true
        # Whether to filter the probe set using the "included" column of the probe set CSV. Default: true
        Boolean filter_probes = true

        # Index of DAPI channel (1-indexed) of fluorescence image, only used in the CytaAssist case, with dark background image.
        Int? dapi_index
        # Use this option if the slide serial number and area identifier have been lost. Choose from visium-1, visium-2 and visium-2-large.
        String? unknown_slide

        # If generate bam outputs. This is also a spaceranger argument.
        Boolean no_bam = false
        # Perform secondary analysis of the gene-barcode matrix (dimensionality reduction, clustering and visualization). Default: false. This is also a spaceranger argument.
        Boolean secondary = false

        # Hard trim the input Read 1 to this length before analysis
        Int? r1_length
        # Hard trim the input Read 2 to this length before analysis
        Int? r2_length

        # Space Ranger version: 2.1.1, 2.0.1, 2.0.0, 1.3.1, 1.3.0
        String spaceranger_version = "2.1.1"
        # Config version
        String config_version = "0.3"

        # Which docker registry to use: quay.io/cumulus (default) or cumulusprod
        String docker_registry = "quay.io/cumulus"
        # spaceranger mkfastq registry, default to gcr.io/broad-cumulus
        String mkfastq_docker_registry = "gcr.io/broad-cumulus"
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Number of cpus per spaceranger and spaceranger job
        Int num_cpu = 32
        # Memory string
        String memory = "120G"

        # Optional disk space for mkfastq.
        Int mkfastq_disk_space = 1500
        # Optional disk space needed for count.
        Int count_disk_space = 500

        # Number of preemptible tries
        Int preemptible = 2
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend = "gcp"
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    String docker_registry_stripped = sub(docker_registry, "/+$", "")
    String mkfastq_docker_registry_stripped = sub(mkfastq_docker_registry, "/+$", "")

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    String null_file = acronym2gsurl["null_file"]

    if (run_mkfastq) {
        call generate_bcl_csv {
            input:
                input_csv_file = input_csv_file,
                output_dir = output_directory_stripped,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }

        if (length(generate_bcl_csv.bcl_csv) > 0) {
            scatter (bcl_csv in generate_bcl_csv.bcl_csv) {
                String key = basename(bcl_csv)
                call srm.spaceranger_mkfastq as spaceranger_mkfastq {
                    input:
                        input_bcl_directory = generate_bcl_csv.inpdirs[key],
                        input_csv_file = bcl_csv,
                        output_directory = output_directory_stripped,
                        delete_input_bcl_directory = delete_input_bcl_directory,
                        barcode_mismatches = mkfastq_barcode_mismatches,
                        spaceranger_version = spaceranger_version,
                        docker_registry = mkfastq_docker_registry_stripped,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        disk_space = mkfastq_disk_space,
                        preemptible = preemptible,
                        awsQueueArn = awsQueueArn,
                        backend = backend
                }
            }
        }
    }

    if (run_count) {
        call generate_count_config {
            input:
                input_csv_file = input_csv_file,
                output_dir = output_directory_stripped,
                fastq_dirs = spaceranger_mkfastq.output_fastqs_flowcell_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                null_file = null_file,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }

        scatter (sample_row in generate_count_config.sample_table) {
            call src.spaceranger_count as spaceranger_count {
                input:
                    sample_id = sample_row[0],
                    input_fastqs_directories = sample_row[1],
                    output_directory = output_directory_stripped,
                    acronym_file = acronym_file,
                    genome = sample_row[2],
                    probe_set = sample_row[3],
                    filter_probes = filter_probes,
                    image = sample_row[4],
                    darkimagestr = sample_row[5],
                    colorizedimage = sample_row[6],
                    cytaimage = sample_row[7],
                    dapi_index = dapi_index,
                    slide = sample_row[8],
                    area = sample_row[9],
                    slidefile = sample_row[10],
                    reorient_images = reorient_images,
                    loupe_alignment = sample_row[11],
                    target_panel = sample_row[12],
                    no_bam = no_bam,
                    secondary = secondary,
                    r1_length = r1_length,
                    r2_length = r2_length,
                    spaceranger_version = spaceranger_version,
                    docker_registry = docker_registry_stripped,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = count_disk_space,
                    preemptible = preemptible,
                    awsQueueArn = awsQueueArn,
                    backend = backend
            }
        }

        call collect_summaries {
            input:
                summaries = spaceranger_count.output_metrics_summary,
                sample_ids = spaceranger_count.output_count_directory,
                config_version = config_version,
                docker_registry = docker_registry_stripped,
                zones = zones,
                preemptible = preemptible,
                awsQueueArn = awsQueueArn,
                backend = backend
        }
    }

    output {
        Array[String]? fastq_outputs = spaceranger_mkfastq.output_fastqs_flowcell_directory
        Array[String]? count_outputs = spaceranger_count.output_count_directory
        File? metrics_summaries = collect_summaries.metrics_summaries
    }
}

task generate_bcl_csv {
    input {
        File input_csv_file
        String output_dir
        String config_version
        String docker_registry
        String zones
        Int preemptible
        String awsQueueArn
        String backend
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

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False)

        df.columns = df.columns.str.strip()
        for c in df.columns:
            df[c] = df[c].str.strip()

        for idx, row in df.iterrows():
            row['Flowcell'] = re.sub('/+$', '', row['Flowcell'])
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)

        with open('inpdirs.txt', 'w') as fo:
            for input_dir in df['Flowcell'].unique():
                run_id = os.path.basename(input_dir)
                bcl_df = df.loc[df['Flowcell'] == input_dir, ['Lane', 'Sample', 'Index']]
                bcl_file = run_id + '_bcl.csv'
                bcl_df.to_csv(bcl_file, index = False)
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
        queueArn: awsQueueArn
    }
}

task generate_count_config {
    input {
        File input_csv_file
        String output_dir
        Array[String]? fastq_dirs
        String config_version
        String docker_registry
        String zones
        String null_file
        Int preemptible
        String awsQueueArn
        String backend
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

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False)

        df.columns = df.columns.str.strip()
        for c in df.columns:
            df[c] = df[c].str.strip()

        for idx, row in df.iterrows():
            row['Flowcell'] = re.sub('/+$', '', row['Flowcell'])
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)

        for c in ['Image', 'ColorizedImage', 'CytaImage', 'SlideFile', 'LoupeAlignment', 'TargetPanel']:
            if c not in df.columns:
                df[c] = '~{null_file}'
            else:
                df.loc[df[c].isna(), c] = '~{null_file}'

        for c in ['DarkImage', 'Slide', 'Area', 'ProbeSet']:
            if c not in df.columns:
                df[c] = ''
            else:
                df.loc[df[c].isna(), c] = ''

        def parse_fastq_dirs(dirs_str):
            r2f = dict()
            if dirs_str == '':
                return r2f
            dirs = dirs_str.split(',')
            for dir in dirs:
                run_id = dir.split('/')[-3].rpartition('_')[0]
                r2f[run_id] = dir
            return r2f

        rid2fdir = parse_fastq_dirs('~{sep="," fastq_dirs}')

        for sample_id in df['Sample'].unique():
            df_local = df.loc[df['Sample'] == sample_id]
            if len(rid2fdir) > 0:
                dirs = df_local['Flowcell'].map(lambda x: rid2fdir[os.path.basename(x)]).values # if also run mkfastq
            else:
                dirs = df_local['Flowcell'].values # if start from count step
            out_str = df_local['Sample'].iat[0] + '\t' + ','.join(dirs) + '\t' + df_local['Reference'].iat[0]
            for c in ['ProbeSet', 'Image', 'DarkImage', 'ColorizedImage', 'CytaImage', 'Slide', 'Area', 'SlideFile', 'LoupeAlignment', 'TargetPanel']:
                out_str += '\t' + df_local[c].iat[0]
            print(out_str)
        CODE
    }

    output {
        Array[Array[String]] sample_table = read_tsv(stdout())
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
