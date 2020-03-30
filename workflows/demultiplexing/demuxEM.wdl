version 1.0

workflow demuxEM {
    input {
        String sample_id
        String output_directory
        String input_rna
        String input_adt_csv
        String genome
        Int min_num_genes
        Float? alpha_on_samples
        Int? min_num_umis
        Float? min_signal_hashtag
        Int? random_state
        Boolean generate_diagnostic_plots
        String? generate_gender_plot

        String docker_registry = "cumulusprod"
        String demuxEM_version = "0.1"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        Int num_cpu = 8
        Int memory = 10
        Int disk_space = 20
        Int preemptible = 2
    }

    call run_demuxEM {
        input:
            sample_id = sample_id,
            output_directory = output_directory,
            input_rna = input_rna,
            input_adt_csv = input_adt_csv,
            genome = genome,
            alpha_on_samples = alpha_on_samples,
            min_num_genes = min_num_genes,
            min_num_umis = min_num_umis,
            min_signal_hashtag = min_signal_hashtag,
            random_state = random_state,
            generate_diagnostic_plots = generate_diagnostic_plots,
            generate_gender_plot = generate_gender_plot,
            docker_registry = docker_registry,
            demuxEM_version = demuxEM_version,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible
    }

    output {
        String output_folder = run_demuxEM.output_folder
        File output_h5sc = run_demuxEM.output_h5sc
        File monitoringLog = run_demuxEM.monitoringLog
    }
}

task run_demuxEM {
    input {
        String sample_id
        String output_directory
        File input_rna
        File input_adt_csv
        String genome
        Int min_num_genes
        Float? alpha_on_samples
        Int? min_num_umis
        Float? min_signal_hashtag
        Int? random_state
        Boolean generate_diagnostic_plots
        String? generate_gender_plot

        String docker_registry
        String demuxEM_version
        String zones
        Int num_cpu
        Int memory
        Int disk_space
        Int preemptible
    }

    command {
        set -e
        export TMPDIR=/tmp
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call
        call_args = ['demuxEM', '~{input_adt_csv}', '~{input_rna}', '~{sample_id}', '-p', '~{num_cpu}']
        if '~{genome}' is not '':
            call_args.extend(['--genome', '~{genome}'])
        if '~{alpha_on_samples}' is not '':
            call_args.extend(['--alpha_on_samples', '~{alpha_on_samples}'])
        if '~{min_num_genes}' is not '':
            call_args.extend(['--min-num-genes', '~{min_num_genes}'])
        if '~{min_num_umis}' is not '':
            call_args.extend(['--min-num-umis', '~{min_num_umis}'])
        if '~{min_signal_hashtag}' is not '':
            call_args.extend(['--min-signal-hashtag', '~{min_signal_hashtag}'])
        if '~{random_state}' is not '':
            call_args.extend(['--random-state', '~{random_state}'])
        if '~{generate_diagnostic_plots}' is 'true':
            call_args.append('--generate-diagnostic-plots')
        if '~{generate_gender_plot}' is not '':
            call_args.extend(['--generate-gender-plot', '~{generate_gender_plot}'])
        print(' '.join(call_args))
        check_call(call_args)
        CODE

        mkdir result
        cp ~{sample_id}_demux.h5sc ~{sample_id}_ADTs.h5ad ~{sample_id}_demux.h5ad ~{sample_id}.*.pdf result
        gsutil -q -m rsync -r result ~{output_directory}/~{sample_id}
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp result/* ~{output_directory}/~{sample_id}
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_h5sc = "result/~{sample_id}_demux.h5sc"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/demuxem:~{demuxEM_version}"
        zones: zones
        memory: "~{memory}G"
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
    }
}