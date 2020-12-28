version 1.0

workflow demuxEM {
    input {
        String sample_id
        String output_directory
        # Input raw RNA expression matrix in 10x hdf5 format.
        String input_rna_h5
        # Input HTO (antibody tag) count matrix in CSV format.
        String input_hto_csv
        # Reference genome name. If not provided, we will infer it from the expression matrix file.
        String genome
        # Only demultiplex cells/nuclei with at least <number> of expressed genes. [default: 100]
        Int min_num_genes
        # The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse. [default: 0.0]
        Float? alpha_on_samples
        # Only demultiplex cells/nuclei with at least <number> of UMIs. [default: 100]
        Int? min_num_umis
        # Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown. [default: 10.0]
        Float? min_signal_hashtag
        # The random seed used in the KMeans algorithm to separate empty ADT droplets from others. [default: 0]
        Int? random_state
        # Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc.
        Boolean generate_diagnostic_plots
        # Generate violin plots using gender-specific genes (e.g. Xist). <gene> is a comma-separated list of gene names.
        String? generate_gender_plot

        String docker_registry = "quay.io/cumulus"
        String demuxEM_version = "0.1.5"
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
            input_rna_h5 = input_rna_h5,
            input_hto_csv = input_hto_csv,
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
        File output_zarr = run_demuxEM.output_zarr
        File monitoringLog = run_demuxEM.monitoringLog
    }
}

task run_demuxEM {
    input {
        String sample_id
        String output_directory
        File input_rna_h5
        File input_hto_csv
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
        call_args = ['demuxEM', '~{input_rna_h5}', '~{input_hto_csv}', '~{sample_id}', '-p', '~{num_cpu}']
        if '~{genome}' is not '':
            call_args.extend(['--genome', '~{genome}'])
        if '~{alpha_on_samples}' is not '':
            call_args.extend(['--alpha-on-samples', '~{alpha_on_samples}'])
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
        cp ~{sample_id}_demux.zarr.zip ~{sample_id}.out.demuxEM.zarr.zip ~{sample_id}.*.pdf result
        gsutil -q -m rsync -r result ~{output_directory}/~{sample_id}
        # mkdir -p ~{output_directory}/~{sample_id}
        # cp result/* ~{output_directory}/~{sample_id}
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr = "result/~{sample_id}_demux.zarr.zip"
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
