version 1.0

workflow shareseq_mkfastq {
    input {
        # Input BCL directory, gs url
        String input_bcl_directory
        # 3 column CSV file (Lane, Sample, Index)
        File input_csv_file
        # Shareseq output directory, gs url
        String output_directory

        # Whether to delete input bcl directory. If false, you should delete this folder yourself so as to not incur storage charges.
        Boolean delete_input_bcl_directory = false
        # Number of allowed mismatches per index
        Int? barcode_mismatches
        # Override the read lengths as specified in RunInfo.xml
        String? use_bases_mask

        # shareseqdemux version
        String shareseqdemux_version = "0.1.0"
        # Which docker registry to use
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per bcl2fastq job
        Int num_cpu = 32
        # Memory string, e.g. 120G
        String memory = "120G"
        # Disk space in GB
        Int disk_space = 1500
        # Number of preemptible tries
        Int preemptible = 2
        # Max number of retries for AWS instance
        Int awsMaxRetries = 5
        # Backend
        String backend = "gcp"
    }

    call run_shareseq_mkfastq {
        input:
            input_bcl_directory = sub(input_bcl_directory, "/+$", ""),
            input_csv_file = input_csv_file,
            output_directory = sub(output_directory, "/+$", ""),
            delete_input_bcl_directory = delete_input_bcl_directory,
            barcode_mismatches = barcode_mismatches,
            shareseqdemux_version = shareseqdemux_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend
    }

    output {
        String output_fastqs_directory = run_shareseq_mkfastq.output_fastqs_directory
        File monitoringLog = run_shareseq_mkfastq.monitoringLog
    }
}

task run_shareseq_mkfastq {
    input {
        String input_bcl_directory
        File input_csv_file
        String output_directory
        Boolean delete_input_bcl_directory
        Int? barcode_mismatches
        String? use_bases_mask
        String shareseqdemux_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    String run_id = basename(input_bcl_directory)

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        strato sync --backend "~{backend}" -m ~{input_bcl_directory} ~{run_id}
        shareseq2bcl ~{input_csv_file} ~{run_id} _bcl_sample_sheet.csv
        
        echo "this is the file being written"
        cat _bcl_sample_sheet.csv

        python <<CODE
        from subprocess import check_call, check_output, CalledProcessError

        call_args = ['bcl2fastq', '-o', '_out', '--sample-sheet' ,'_bcl_sample_sheet.csv', '--runfolder-dir', '~{run_id}', '--create-fastq-for-index-reads']
        if '~{barcode_mismatches}' != '':
            call_args.apppend('--barcode-mismatches ~{barcode_mismatches}')
        if '~{use_bases_mask}' != '':
            call_args.apppend('--use-bases-mask ~{use_bases_mask}')
        print(' '.join(call_args))
        check_call(call_args)

        call_args = ['strato', 'sync', '--backend', '~{backend}', '-m', '_out', '~{output_directory}/~{run_id}_fastqs']
        print(' '.join(call_args))
        check_call(call_args)

        if '~{delete_input_bcl_directory}' is 'true':
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', '~{output_directory}/~{run_id}_fastqs/input_samplesheet.csv']
                print(' '.join(call_args))
                check_output(call_args)
                call_args = ['strato', 'rm', '--backend', '~{backend}', '-m', '-r', '~{input_bcl_directory}']
                print(' '.join(call_args))
                check_call(call_args)
                print('~{input_bcl_directory} is deleted!')
            except CalledProcessError:
                print("Failed to delete BCL directory.")
        CODE
    }

    output {
        String output_fastqs_directory = "~{output_directory}/~{run_id}_fastqs"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/shareseqdemux:~{shareseqdemux_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
