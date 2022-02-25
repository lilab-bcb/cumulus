version 1.0

import "starsolo_count.wdl" as ssc

workflow starsolo_workflow {
    input {
        # Input CSV sample sheet describing metadata of each sample
        File input_csv_file
        # URL of output directory
        String output_directory
        String read1_fastq_pattern = "_S*_L*_R1_001.fastq.gz"
        String read2_fastq_pattern = "_S*_L*_R2_001.fastq.gz"
        # Which read contains cell barcodes and UMIs. Either "read1" or "read2"
        String barcode_read = "read1"
        # Type of SAM/BAM output
        String? outSAMtype
        # Type of single-cell RNA-seq
        String? soloType
        # Cell barcode start position (1-based coordinate)
        Int? soloCBstart
        # Cell barcode length
        Int? soloCBlen
        # UMI start position (1-based coordinate)
        Int? soloUMIstart
        # UMI length
        Int? soloUMIlen
        # Cell barcode white list
        String soloCBwhitelist = ""
        # Length of the barcode read
        Int? soloBarcodeReadLength
        # Identifies which read mate contains the barcode (CB+UMI) sequence
        Int? soloBarcodeMate
        # position of Cell Barcode(s) on the barcode read. Presently only works with --soloType CB_UMI_Complex, and barcodes are assumed to be on Read2.
        String? soloCBposition
        # position of the UMI on the barcode read, same as CBposition
        String? soloUMIposition
        # adapter sequence to anchor barcodes.
        String? soloAdapterSequence
        # maximum number of mismatches allowed in adapter sequence.
        Int? soloAdapterMismatchesNmax
        # matching the Cell Barcodes to the WhiteList
        String? soloCBmatchWLtype
        # when inputting reads from a SAM file (--readsFileType SAM SE/PE), these SAM attributes mark the barcode sequence (in proper order).
        # For instance, for 10X CellRanger or STARsolo BAMs, use --soloInputSAMattrBarcodeSeq CR UR .
        # This parameter is required when running STARsolo with input from SAM.
        String? soloInputSAMattrBarcodeSeq
        # when inputting reads from a SAM file (--readsFileType SAM SE/PE), these SAM attributes mark the barcode qualities (in proper order).
        # For instance, for 10X CellRanger or STARsolo BAMs, use --soloInputSAMattrBarcodeQual CY UY .
        # If this parameter is '-' (default), the quality 'H' will be assigned to all bases.
        String? soloInputSAMattrBarcodeQual
        # strandedness of the solo libraries
        String? soloStrand
        # genomic features for which the UMI counts per Cell Barcode are collected
        String? soloFeatures
        # counting method for reads mapping to multiple genes
        String? soloMultiMappers
        # type of UMI deduplication (collapsing) algorithm
        String? soloUMIdedup
        # type of UMI filtering (for reads uniquely mapping to genes)
        String? soloUMIfiltering
        # cell filtering type and parameters
        String? soloCellFilter
        # field 3 in the Gene features.tsv file. If "-", then no 3rd field is output.
        String? soloOutFormatFeaturesGeneField3
        # Number of CPUs to request per sample
        Int num_cpu = 32
        # STAR version to use. Currently support: 2.7.9a
        String star_version = "2.7.9a"
        # Docker registry, default to quay.io/cumulus
        String docker_registry = "quay.io/cumulus"
        # Reference Index TSV
        File acronym_file = "gs://regev-lab/resources/starsolo/index.tsv"
        # Disk space in GB needed per sample
        Int disk_space = 500
        # Google cloud zones to consider for execution
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Memory size string represent memory needed for count per sample
        String memory = "120G"
        # Number of maximum preemptible tries allowed
        Int preemptible = 2
        # Number of maximum retries when running on AWS
        Int awsMaxRetries = 5
        # backend choose from "gcp", "aws", "local"
        String backend = "gcp"
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")

    call generate_count_config {
        input:
            input_csv_file = input_csv_file,
            docker_registry = docker_registry,
            zones = zones,
            star_version = star_version,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend
    }

    if (length(generate_count_config.sample_ids) > 0) {
        scatter (sample_id in generate_count_config.sample_ids) {
            call ssc.starsolo_count as starsolo_count {
                input:
                    sample_id = sample_id,
                    input_fastqs_directories = generate_count_config.sample2dir[sample_id],
                    genome = generate_count_config.sample2genome[sample_id],
                    assay = generate_count_config.sample2assay[sample_id],
                    acronym_file = acronym_file,
                    read1_fastq_pattern = read1_fastq_pattern,
                    read2_fastq_pattern = read2_fastq_pattern,
                    barcode_read = barcode_read,
                    outSAMtype = outSAMtype,
                    soloType = soloType,
                    soloCBstart = soloCBstart,
                    soloCBlen = soloCBlen,
                    soloUMIstart = soloUMIstart,
                    soloUMIlen = soloUMIlen,
                    soloCBwhitelist = soloCBwhitelist,
                    soloBarcodeReadLength = soloBarcodeReadLength,
                    soloBarcodeMate = soloBarcodeMate,
                    soloCBposition = soloCBposition,
                    soloUMIposition = soloUMIposition,
                    soloAdapterSequence = soloAdapterSequence,
                    soloAdapterMismatchesNmax = soloAdapterMismatchesNmax,
                    soloCBmatchWLtype = soloCBmatchWLtype,
                    soloInputSAMattrBarcodeSeq = soloInputSAMattrBarcodeSeq,
                    soloInputSAMattrBarcodeQual = soloInputSAMattrBarcodeQual,
                    soloStrand = soloStrand,
                    soloFeatures = soloFeatures,
                    soloMultiMappers = soloMultiMappers,
                    soloUMIdedup = soloUMIdedup,
                    soloUMIfiltering = soloUMIfiltering,
                    soloCellFilter = soloCellFilter,
                    soloOutFormatFeaturesGeneField3 = soloOutFormatFeaturesGeneField3,
                    output_directory = output_directory_stripped,
                    docker_registry = docker_registry,
                    star_version = star_version,
                    zones = zones,
                    num_cpu = num_cpu,
                    memory = memory,
                    disk_space = disk_space,
                    preemptible = preemptible,
                    awsMaxRetries = awsMaxRetries,
                    backend = backend
            }
        }
    }

    output {
        Array[File]? starsoloLogs = starsolo_count.starsoloLog
        String output_folder = output_directory
    }

}

task generate_count_config {
    input {
        File input_csv_file
        String zones
        String star_version
        String docker_registry
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        import re, sys
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv('~{input_csv_file}', header = 0, dtype = str, index_col = False)
        df.columns = df.columns.str.strip()
        for c in df.columns:
            df[c] = df[c].str.strip()

        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        if 'Assay' not in df.columns:
            df['Assay'] = 'tenX_v3'
        else:
            df.loc[df['Assay'].isna(), 'Assay'] = 'tenX_v3'

        for idx, row in df.iterrows():
            row['Location'] = re.sub('/+$', '', row['Location'])
            if re.search('[^a-zA-Z0-9_-]', row['Sample']) is not None:
                print('Sample must contain only alphanumeric characters, hyphens, and underscores.', file = sys.stderr)
                print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &', file = sys.stderr)
                sys.exit(1)

        with open('sample_ids.txt', 'w') as fo1, open('sample2dir.txt', 'w') as fo2, open('sample2genome.txt', 'w') as fo3, open('sample2assay.txt', 'w') as fo4:
            for sample_id in df['Sample'].unique():
                df_local = df.loc[df['Sample'] == sample_id]

                dirs = df_local['Location'].values

                if df_local['Reference'].unique().size > 1:
                    print("Detected multiple references for sample " + sample_id + "!", file = sys.stderr)
                    sys.exit(1)
                reference = df_local['Reference'].iat[0]

                if df_local['Assay'].unique().size > 1:
                    print("Detected multiple assay names for sample " + sample_id + "!", file = sys.stderr)
                assay = df_local['Assay'].iat[0]

                fo1.write(sample_id + '\n')
                fo2.write(sample_id + '\t' + ','.join(dirs) + '\n')
                fo3.write(sample_id + '\t' + reference + '\n')
                fo4.write(sample_id + '\t' + assay + '\n')
        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] sample2dir = read_map('sample2dir.txt')
        Map[String, String] sample2genome = read_map('sample2genome.txt')
        Map[String, String] sample2assay = read_map('sample2assay.txt')
    }

    runtime {
        docker: "~{docker_registry}/starsolo:~{star_version}"
        zones: zones
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
