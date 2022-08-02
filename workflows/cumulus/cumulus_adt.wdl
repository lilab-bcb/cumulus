version 1.0

workflow cumulus_adt {
    input {
        # Sample ID
        String sample_id
        # A comma-separated list of input FASTQs directories (gs urls)
        String input_fastqs_directories
        # Output directory, gs url
        String output_directory

        # 10x genomics chemistry
        String chemistry

        # data type, either adt or crispr
        String data_type

        # feature barcodes in csv format
        File feature_barcode_file

        # Index TSV file
        File acronym_file

        # Barcode start position at Read 2 (0-based coordinate) for CRISPR
        Int? crispr_barcode_pos
        # scaffold sequence for CRISPR, default is ""
        String? scaffold_sequence

        # maximum hamming distance in feature barcodes, change default to 2.
        Int max_mismatch = 2
        # minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for crispr
        Float min_read_ratio = 0.1

        # cumulus_feature_barcoding version
        String cumulus_feature_barcoding_version
        # Which docker registry to use: quay.io/cumulus (default), or cumulusprod
        String docker_registry = "quay.io/cumulus"

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per adt job
        Int num_cpu = 4
        # Memory string, e.g. 32G
        String memory = "32G"
        # Disk space in GB
        Int disk_space = 100

        # Number of preemptible tries
        Int preemptible = 2
        # Number of maximum retries when running on AWS
        Int awsMaxRetries = 5
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend = "gcp"
    }

    Map[String, String] acronym2gsurl = read_map(acronym_file)
    File cell_barcode_file = acronym2gsurl[chemistry]


    call run_generate_count_matrix_ADTs {
        input:
            sample_id = sample_id,
            input_fastqs_directories = input_fastqs_directories,
            output_directory = output_directory,
            chemistry = chemistry,
            data_type = data_type,
            cell_barcodes = cell_barcode_file,
            feature_barcodes = feature_barcode_file,
            crispr_barcode_pos = crispr_barcode_pos,
            scaffold_sequence = scaffold_sequence,
            max_mismatch = max_mismatch,
            min_read_ratio = min_read_ratio,
            cumulus_feature_barcoding_version = cumulus_feature_barcoding_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    call run_generate_ridge_plot {
        input:
            sample_id = sample_id,
            output_count_directory = run_generate_count_matrix_ADTs.output_count_directory,
            cumulus_feature_barcoding_version = cumulus_feature_barcoding_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            awsQueueArn = awsQueueArn,
            backend = backend
    }

    output {
        String output_count_directory = run_generate_count_matrix_ADTs.output_count_directory
        String ridge_plot_pdf = run_generate_ridge_plot.ridge_plot_pdf
        File monitoringLog = run_generate_count_matrix_ADTs.monitoringLog
    }
}
task run_generate_count_matrix_ADTs {
    input {
        String sample_id
        String input_fastqs_directories
        String output_directory
        String chemistry
        String data_type
        File cell_barcodes
        File feature_barcodes
        Int? crispr_barcode_pos
        String? scaffold_sequence
        Int max_mismatch
        Float min_read_ratio
        String cumulus_feature_barcoding_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
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
                call_args = ['strato', 'exists', '--backend', '~{backend}', directory + '/~{sample_id}/']
                print(' '.join(call_args))
                check_call(call_args, stderr=STDOUT, stdout=DEVNULL)
                call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', directory + '/~{sample_id}', target]
                print(' '.join(call_args))
                check_call(call_args)
            except CalledProcessError:
                if not os.path.exists(target):
                    os.mkdir(target)
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', directory + '/~{sample_id}' + '_S*_L*_*_001.fastq.gz' , target]
                print(' '.join(call_args))
                check_call(call_args)
            fastqs.append(target)

        call_args = ['generate_count_matrix_ADTs', '~{cell_barcodes}', '~{feature_barcodes}', ','.join(fastqs), '~{sample_id}', '-p', '~{num_cpu}', '--max-mismatch-feature', '~{max_mismatch}']
        if '~{data_type}' == 'crispr':
            call_args.extend(['--feature', 'crispr'])
            if '~{scaffold_sequence}' != '':
                call_args.extend(['--scaffold-sequence', '~{scaffold_sequence}'])
            if '~{crispr_barcode_pos}' != '':
                call_args.extend(['--barcode-pos', '~{crispr_barcode_pos}'])
            if '~{chemistry}' == 'SC3Pv3' and '~{crispr_barcode_pos}' == '' and '~{scaffold_sequence}' == '':
                call_args.append('--convert-cell-barcode')
        else:
            call_args.extend(['--feature', 'antibody'])
            if '~{data_type}' == 'cmo':
                call_args.append('--convert-cell-barcode')
        if '~{chemistry}' == 'SC3Pv3':
            call_args.extend(['--max-mismatch-cell', '0', '--umi-length', '12'])
        elif '~{chemistry}' == 'multiome':
            call_args.extend(['--max-mismatch-cell', '1', '--umi-length', '12'])
        else:
            call_args.extend(['--max-mismatch-cell', '1', '--umi-length', '10'])
        print(' '.join(call_args))
        check_call(call_args)
     
        report_file = '~{sample_id}'+'.report.txt'
        call_args = ['convert_to_pdf.py', report_file, '~{sample_id}'+'.report.pdf']
        print(' '.join(call_args))
        check_call(call_args)

        CODE

        if [ -f "~{sample_id}".stat.csv.gz ]
        then
            filter_chimeric_reads ~{data_type} ~{feature_barcodes} "~{sample_id}.stat.csv.gz" ~{min_read_ratio} ~{sample_id}
        fi

        if [ -f "~{sample_id}".crispr.stat.csv.gz ]
        then
            filter_chimeric_reads ~{data_type} ~{feature_barcodes} "~{sample_id}.crispr.stat.csv.gz" ~{min_read_ratio} ~{sample_id}
        fi

        strato cp --backend ~{backend} -m "~{sample_id}".*csv* "~{sample_id}".report.txt "~{sample_id}".report.pdf "~{output_directory}/~{sample_id}/"

        if [ -f "~{sample_id}".umi_count.pdf ]
        then
            strato cp --backend ~{backend} "~{sample_id}".umi_count.pdf "~{output_directory}/~{sample_id}/"
        fi
    }

    output {
        String output_count_directory = "~{output_directory}/~{sample_id}"
        String output_text_summary = "~{output_directory}/~{sample_id}/~{sample_id}.report.pdf"
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
        maxRetries: if backend == "aws" then awsMaxRetries else 0
        queueArn: awsQueueArn
    }
}

task run_generate_ridge_plot {
    input {
        String sample_id
        String output_count_directory
        String cumulus_feature_barcoding_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String awsQueueArn
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        python <<CODE

        from glob import glob

        call_args = ["strato", "cp", "--backend", "~{backend}", "~{output_count_directory}/~{sample_id}/~{sample_id}.*.csv" , "."]
        print(' '.join(call_args))
        check_call(call_args)

        call_args = ["strato", "cp", "--backend", "~{backend}", "~{output_count_directory}/~{sample_id}/~{sample_id}.report.pdf" , "."]
        print(' '.join(call_args))
        check_call(call_args)

        for count_matrix_file in glob("~{sample_id}.*.csv"):
            output_plot_pdf = count_matrix+"_ridge_plot.pdf"
            count_matrix = os.path.splitext(os.path.basename(count_matrix_file))[0]
            call_args = ["generate_ridgeplots.py", count_matrix_file, count_matrix, output_plot_pdf]
            print(' '.join(call_args))
            check_call(call_args)

        report_pdf = "~{sample_id}.report.pdf"
        ridgeplot_pdf_list = glob("*_ridge_plot.pdf")
        files_to_merge = ",".join([report_pdf] + ridgeplot_pdf_list)

        call_args = ["merge_pdf.py", files_to_merge, "~{sample_id}"+".report_ridgeplot.pdf"]
        print(' '.join(call_args))
        check_call(call_args)

        CODE
        
        strato cp --backend "~{backend}" -m "~{sample_id}"_ridge_plot.pdf "~{output_count_directory}/~{sample_id}/"

    }

    output {
        String ridge_plot_pdf = "~{output_count_directory}/~{sample_id}/~{sample_id}_ridge_plot.pdf"
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
        maxRetries: if backend == "aws" then awsMaxRetries else 0
        queueArn: awsQueueArn
    }
}
