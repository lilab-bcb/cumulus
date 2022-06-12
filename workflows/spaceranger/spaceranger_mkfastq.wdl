version 1.0

workflow spaceranger_mkfastq {
    input {
        # Input BCL directory, gs url
        String input_bcl_directory
        # 3 column CSV file (Lane, Sample, Index)
        File input_csv_file
        # spaceranger output directory, gs url
        String output_directory

        # Whether to delete input bcl directory. If false, you should delete this folder yourself so as to not incur storage charges
        Boolean delete_input_bcl_directory = false
        # Number of allowed mismatches per index
        Int? barcode_mismatches

        # spaceranger version
        String spaceranger_version
        # Which docker registry to use
        String docker_registry

        # Google cloud zones, default to "us-central1-b", which is consistent with CromWell's genomics.default-zones attribute
        String zones = "us-central1-b"
        # Number of cpus per spaceranger job
        Int num_cpu = 32
        # Memory string
        String memory = "120G"
        # Disk space in GB
        Int disk_space = 1500
        # Number of preemptible tries
        Int preemptible = 2
        # Number of maximum retries when running on AWS
        Int awsMaxRetries = 5
        # Arn string of AWS queue
        String awsQueueArn = ""
        # Backend
        String backend
    }

    call run_spaceranger_mkfastq {
        input:
            input_bcl_directory = input_bcl_directory,
            input_csv_file = input_csv_file,
            output_directory = output_directory,
            delete_input_bcl_directory = delete_input_bcl_directory,
            barcode_mismatches = barcode_mismatches,
            spaceranger_version = spaceranger_version,
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
        String output_fastqs_directory = run_spaceranger_mkfastq.output_fastqs_directory
        String output_fastqs_flowcell_directory = run_spaceranger_mkfastq.output_fastqs_flowcell_directory
        File monitoringLog = run_spaceranger_mkfastq.monitoringLog
    }
}

task run_spaceranger_mkfastq {
    input {
        String input_bcl_directory
        File input_csv_file
        String output_directory
        Boolean delete_input_bcl_directory
        Int? barcode_mismatches
        String spaceranger_version
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

    String run_id = basename(input_bcl_directory)

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &
        strato sync --backend ~{backend} -m ~{input_bcl_directory} ~{run_id}

        python <<CODE
        import os
        import glob
        import sys
        import pandas as pd
        import subprocess

        from packaging import version

        mkfastq_args = ['spaceranger', 'mkfastq', '--id=results', '--run=~{run_id}', '--csv=~{input_csv_file}', '--jobmode=local']
        if version.parse("~{spaceranger_version}") < version.parse("1.3.0"):
            mkfastq_args.append('--qc')
        barcode_mismatches = '~{barcode_mismatches}'
        if barcode_mismatches != '':
            mkfastq_args += ['--barcode-mismatches', barcode_mismatches]
        p = subprocess.run(mkfastq_args)

        if p.returncode != 0: # ./MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/chnk0-u8d92d5526b/_stderr
            if os.path.exists('results/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/'):
                output_dirs = os.listdir('results/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/')
                for output in output_dirs:
                    if output.startswith('chnk0-'):
                        break
                with open(os.path.join('results/MAKE_FASTQS_CS/MAKE_FASTQS/BCL2FASTQ_WITH_SAMPLESHEET/fork0/', output, '_stderr'), "r") as error_in:
                    for line in error_in:
                        print(line, file=sys.stderr)
            sys.exit(1)

        with open("output_fastqs_flowcell_directory.txt", "w") as fout:
            flowcell = [name for name in os.listdir('results/outs/fastq_path') if name != 'Reports' and name != 'Stats' and os.path.isdir('results/outs/fastq_path/' + name)][0]
            fout.write('~{output_directory}/~{run_id}_spatialfastqs/fastq_path/' + flowcell + '\n')

        CODE

        strato sync --backend ~{backend} -m results/outs "~{output_directory}/~{run_id}_spatialfastqs"

        python <<CODE
        from subprocess import check_call, check_output, CalledProcessError
        if '~{delete_input_bcl_directory}' == 'true':
            try:
                call_args = ['strato', 'exists', '--backend', '~{backend}', '~{output_directory}/~{run_id}_spatialfastqs/input_samplesheet.csv']
                print(' '.join(call_args))
                check_output(call_args)
                call_args = ['strato', 'rm', '-m', '-r', '~{input_bcl_directory}']
                print(' '.join(call_args))
                check_call(call_args)
                print('~{input_bcl_directory} is deleted!')
            except CalledProcessError:
                print("Failed to delete BCL directory.")
        CODE
    }

    output {
        String output_fastqs_directory = "~{output_directory}/~{run_id}_spatialfastqs"
        String output_fastqs_flowcell_directory = read_lines("output_fastqs_flowcell_directory.txt")[0]
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/spaceranger:~{spaceranger_version}"
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
