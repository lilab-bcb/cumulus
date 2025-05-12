version 1.0

workflow cumulus_adt {
    input {
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # Output directory, gs url
        String output_directory

        # genome reference
        String genome
        # 10x genomics chemistry
        String chemistry

        # data type, either hashing, citeseq, cmo, crispr or adt
        String data_type

        # feature barcodes in csv format
        File feature_barcode_file

        # Barcode start position at Read 2 (0-based coordinate) for CRISPR
        Int? crispr_barcode_pos
        # scaffold sequence for CRISPR, default is ""
        String? scaffold_sequence

        # maximum hamming distance in feature barcodes, change default to 2.
        Int max_mismatch = 2
        # PCR chimeric filitering: minimum read count ratio cutoff (non-inclusive) to justify a feature given a cell barcode and UMI combination, only used for crispr
        Float read_ratio_cutoff = 0.5

        # cumulus_feature_barcoding version
        String cumulus_feature_barcoding_version
        # Which docker registry to use: quay.io/cumulus (default), or cumulusprod
        String docker_registry = "quay.io/cumulus"

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per adt job
        Int num_cpu = 3
        # Memory string, e.g. 32G
        String memory = "32G"
        # Disk space in GB
        Int disk_space = 100

        # Number of preemptible tries
        Int preemptible = 2
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend = "gcp"
    }

    call run_generate_count_matrix_ADTs {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            output_directory = output_directory,
            genome = genome,
            chemistry = chemistry,
            data_type = data_type,
            feature_barcodes = feature_barcode_file,
            crispr_barcode_pos = crispr_barcode_pos,
            scaffold_sequence = scaffold_sequence,
            max_mismatch = max_mismatch,
            read_ratio_cutoff = read_ratio_cutoff,
            cumulus_feature_barcoding_version = cumulus_feature_barcoding_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    output {
        String output_count_directory = run_generate_count_matrix_ADTs.output_count_directory
        File monitoringLog = run_generate_count_matrix_ADTs.monitoringLog
    }
}

task run_generate_count_matrix_ADTs {
    input {
        String sample_id
        String input_fastqs_directories
        String output_directory
        String genome
        String chemistry
        String data_type
        File feature_barcodes
        Int? crispr_barcode_pos
        String? scaffold_sequence
        Int max_mismatch
        Float read_ratio_cutoff
        String cumulus_feature_barcoding_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        python <<CODE
        import re
        from subprocess import check_call, CalledProcessError, STDOUT, DEVNULL
        import os
        import glob

        fastqs = []

        for i, directory in enumerate('~{input_fastqs_directories}'.split(',')):
            directory = re.sub('/+$', '', directory) # remove trailing slashes
            target = '~{sample_id}_' + str(i)
            try:
                call_args = ['strato', 'exists', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stderr=STDOUT, stdout=DEVNULL)
                call_args = ['strato', 'sync', '-m', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)
            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '-m', directory + '/~{sample_id}' + '_S*_L*_*_001.fastq.gz' , target]
                print(' '.join(call_args))
                check_call(call_args)
            fastqs.append(target)

        call_args = ['generate_feature_barcode_matrices', '/software/barcode_list', '~{data_type}', '~{feature_barcodes}', ','.join(fastqs), '~{sample_id}', '-p', '~{num_cpu}', '--max-mismatch-feature', '~{max_mismatch}']
        if '~{genome}' != '':
            call_args.extend(['--genome', '~{genome}'])
        if '~{chemistry}' != '':
            call_args.extend(['--chemistry', '~{chemistry}'])
        if '~{data_type}' == 'crispr':
            call_args.extend(['--read-ratio-cutoff', '~{read_ratio_cutoff}'])
            if '~{scaffold_sequence}' != '':
                call_args.extend(['--scaffold-sequence', '~{scaffold_sequence}'])
            if '~{crispr_barcode_pos}' != '':
                call_args.extend(['--barcode-pos', '~{crispr_barcode_pos}'])
        print(' '.join(call_args))
        check_call(call_args)

        CODE

        strato cp -m "~{sample_id}".*h5 "~{sample_id}".report.txt "~{output_directory}/~{sample_id}/"
    }

    output {
        String output_count_directory = "~{output_directory}/~{sample_id}"
        String output_text_summary = "~{output_directory}/~{sample_id}/~{sample_id}.report.txt"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cumulus_feature_barcoding:~{cumulus_feature_barcoding_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        queueArn: awsQueueArn
    }
}
