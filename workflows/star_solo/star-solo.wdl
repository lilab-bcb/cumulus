version 1.0

import "starsolo_per_sample.wdl" as stps

workflow starsolo {
    input {
        # Input TSV sample sheet describing metadata of each sample
        File input_tsv_file
        # Genome reference
        String genome
        # Preset options, choosing from tenX_v3 (for 10X V3 chemistry), tenX_v2 (for 10X V2 chemistry), DropSeq, SeqWell and custom
        String preset = ""
        # Cell barcode start position (1-based coordinate)
        String? soloType
        Int? soloCBstart
        # Cell barcode length
        Int? soloCBlen
        # UMI start position (1-based coordinate)
        Int? soloUMIstart
        # UMI length
        Int? soloUMIlen
        # Cell barcode white list
        File? soloCBwhitelist
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
        # file names for STARsolo output
        String? soloOutFileNames
        # cell filtering type and parameters
        String? soloCellFilter
        # field 3 in the Gene features.tsv file. If "-", then no 3rd field is output.
        String? soloOutFormatFeaturesGeneField3
        # URL of output directory
        String output_directory
        # Number of CPUs to request per sample
        Int num_cpu = 32
        # STAR version to use. Currently support: 2.7.9a, 2.7.6a.
        String star_version = "2.7.9a"
        # Docker registry, default to quay.io/cumulus
        String docker_registry = "quay.io/cumulus"
        # Reference Index TSV
        File ref_index_file = "gs://regev-lab/resources/count_tools/ref_index.tsv"
        # Whitelist Index TSV
        File whitelist_index_file = "gs://regev-lab/resources/count_tools/whitelist_index.tsv"
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

    Map[String, String] ref_index2url = read_map(ref_index_file)
    String genome_url = if sub(genome, "^(gs|s3)://.*", "URL")=="URL" then genome else ref_index2url[genome] + '/starsolo.tar.gz'

    Map[String, String] wl_index2url = read_map(whitelist_index_file)

    call generate_count_config {
        input:
            input_tsv_file = input_tsv_file,
            docker_registry = docker_registry,
            zones = zones,
            star_version = star_version,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend
    }

    if (length(generate_count_config.sample_ids) > 0) {
        scatter (sample_id in generate_count_config.sample_ids) {
            call stps.run_starsolo_per_sample as run_starsolo_per_sample {
                input:
                    sample_id = sample_id,
                    r1_fastqs = generate_count_config.id2r1[sample_id],
                    r2_fastqs = generate_count_config.id2r2[sample_id],
                    preset = preset,
                    genome = genome_url,
                    soloType = soloType,
                    soloCBstart = soloCBstart,
                    soloCBlen = soloCBlen,
                    soloUMIstart = soloUMIstart,
                    soloUMIlen = soloUMIlen,
                    soloCBwhitelist = if soloCBwhitelist != '' then soloCBwhitelist else wl_index2url[preset],
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
                    soloOutFileNames = soloOutFileNames,
                    soloCellFilter = soloCellFilter,
                    soloOutFormatFeaturesGeneField3 = soloOutFormatFeaturesGeneField3,
                    output_directory = output_directory_stripped,
                    docker_registry = docker_registry,
                    version = star_version,
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
        String output_folder = output_directory
    }

}

task generate_count_config {
    input {
        File input_tsv_file
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
        import re
        import pandas as pd
        from subprocess import check_call

        df = pd.read_csv('~{input_tsv_file}', sep = '\t', header = 0, dtype = str, index_col = False)
        for c in df.columns:
            df[c] = df[c].str.strip()

        regex_pat = re.compile('[^a-zA-Z0-9_-]')
        if any(df['Sample'].str.contains(regex_pat)):
            print('Sample must contain only alphanumeric characters, hyphens, and underscores.')
            print('Examples of common characters that are not allowed are the space character and the following: ?()[]/\=+<>:;"\',*^| &')
            sys.exit(1)

        with open('sample_ids.txt', 'w') as fo1, open('sample_r1.tsv', 'w') as fo2, open('sample_r2.tsv', 'w') as fo3:
            for idx, row in df.iterrows():
                fo1.write(row['Sample'] + '\n')

                if 'Flowcells' in df.columns: # Fetch R1 and R2 fastqs automatically.
                    input_dir_list = list(map(lambda x: x.strip(), row['Flowcells'].split(',')))
                    r1_list = []
                    r2_list = []
                    for directory in input_dir_list:
                        directory = re.sub('/+$', '', directory)

                        call_args = ['gsutil', 'ls', directory]
                        # call_args = ['ls', directory]
                        with open('list_dir.txt', 'w') as tmp_fo:
                            check_call(call_args, stdout=tmp_fo)

                        with open('list_dir.txt', 'r') as tmp_fin:
                            f_list = tmp_fin.readlines()
                            f_list = list(map(lambda s: s.strip(), f_list))

                        r1_files = [f for f in f_list if re.match('.*_R1_.*.fastq.gz', f)]
                        r2_files = [f for f in f_list if re.match('.*_R2_.*.fastq.gz', f)]
                        r1_files.sort()
                        r2_files.sort()
                        # r1_files = list(map(lambda s: directory+'/'+s, r1_files))
                        # r2_files = list(map(lambda s: directory+'/'+s, r2_files))

                        r1_list.extend(r1_files)
                        r2_list.extend(r2_files)

                else:  # R1 and R2 fastqs specified in sample sheet.
                    r1_list = list(map(lambda s: s.strip(), row['R1'].split(',')))
                    r2_list = list(map(lambda s: s.strip(), row['R2'].split(',')))

                fo2.write(row['Sample'] + '\t' + ','.join(r1_list) + '\n')
                fo3.write(row['Sample'] + '\t' + ','.join(r2_list) + '\n')
        CODE
    }

    output {
        Array[String] sample_ids = read_lines('sample_ids.txt')
        Map[String, String] id2r1 = read_map('sample_r1.tsv')
        Map[String, String] id2r2 = read_map('sample_r2.tsv')
    }

    runtime {
        docker: "~{docker_registry}/starsolo:~{star_version}"
        zones: zones
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
