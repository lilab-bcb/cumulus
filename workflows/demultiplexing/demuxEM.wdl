version 1.0

workflow demuxEM {
    input {
        # Sample ID
        String sample_id
        # Output directory (gs url + path)
        String output_directory
        # Input raw RNA expression matrix in 10x hdf5 format
        String input_rna_h5
        # Input HTO (antibody tag) count matrix in CSV format
        String input_hto_csv
        # Reference genome name. If not provided, we will infer it from the expression matrix file
        String genome
        # Only demultiplex cells/nuclei with at least <number> of expressed genes
        Int min_num_genes
        # The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse
        Float? alpha_on_samples
        # Only demultiplex cells/nuclei with at least <number> of UMIs
        Int? min_num_umis
        # Any cell/nucleus with less than <count> hashtags from the signal will be marked as unknown
        Float? min_signal_hashtag
        # The random seed used in the KMeans algorithm to separate empty ADT droplets from others
        Int? random_state
        # Generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells etc.
        Boolean generate_diagnostic_plots
        # Generate violin plots using gender-specific genes (e.g. Xist). <gene> is a comma-separated list of gene names
        String? generate_gender_plot
        # Which docker registry to use
        String docker_registry
        # DemuxEM version
        String demuxEM_version
        # Google cloud zones
        String zones
        # Number of CPUs used
        Int num_cpu
        # Memory size string for demuxEM
        String memory
        # Disk space in GB
        Int disk_space
        # Number of preemptible tries
        Int preemptible
        # Number of maximum retries when running on AWS
        Int awsMaxRetries
        # Backend
        String backend
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
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend
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
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
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
        cp "~{sample_id}_demux".zarr.zip "~{sample_id}".out.demuxEM.zarr.zip "~{sample_id}".*.pdf result
        strato sync --backend ~{backend} -m result "~{output_directory}/~{sample_id}"
    }

    output {
        String output_folder = "~{output_directory}/~{sample_id}"
        File output_zarr = "result/~{sample_id}_demux.zarr.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/demuxem:~{demuxEM_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
