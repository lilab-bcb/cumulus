version 1.0

workflow cellranger_vdj_create_reference {
    input {
        String docker_registry = "cumulusprod"
        String cellranger_version = '4.0.0'
        Int disk_space = 100
        Int preemptible = 2
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String memory = "32G"

        File input_fasta
        File input_gtf
        # Output directory, gs URL
        String output_directory
        String genome
        String ref_version = ""
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "/+$", "")

    call run_cellranger_vdj_create_reference {
        input:
            docker_registry = docker_registry,
            cellranger_version = cellranger_version,
            disk_space = disk_space,
            preemptible = preemptible,
            zones = zones,
            memory = memory,
            input_fasta = input_fasta,
            input_gtf = input_gtf,
            output_dir = output_directory_stripped,
            genome = genome,
            ref_version = ref_version
    }
}

task run_cellranger_vdj_create_reference {
    input {
        String docker_registry
        String cellranger_version
        Int disk_space
        Int preemptible
        String zones
        String memory
        File input_fasta
        File input_gtf
        String output_dir
        String genome
        String ref_version
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        import os
        from subprocess import check_call

        # Unzip fa if needed.
        fa_file = '~{input_fasta}'
        root, ext = os.path.splitext(fa_file)
        if ext == '.gz':
            call_args = ['gunzip', '-f', fa_file]
            print(' '.join(call_args))
            check_call(call_args)
            fa_file = root

        # Unzip gtf if needed.
        gtf_file = '~{input_gtf}'
        root, ext = os.path.splitext(gtf_file)
        if ext == '.gz':
            call_args = ['gunzip', '-f', gtf_file]
            print(' '.join(call_args))
            check_call(call_args)
            gtf_file = root

        call_args = ['cellranger', 'mkvdjref', '--genome=~{genome}', '--fasta=' + fa_file, '--genes=' + gtf_file]

        if '~{ref_version}' is not '':
            call_args.append('--ref-version=~{ref_version}')

        print(' '.join(call_args))
        check_call(call_args)
        CODE

        tar -czf ~{genome}.tar.gz ~{genome}
        gsutil -m cp ~{genome}.tar.gz ~{output_dir}
        # mkdir -p ~{output_dir}
        # cp ~{genome}.tar.gz ~{output_dir}
    }

    output {
        File output_reference = "~{genome}.tar.gz"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cellranger:~{cellranger_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: "~{preemptible}"
    }
}
