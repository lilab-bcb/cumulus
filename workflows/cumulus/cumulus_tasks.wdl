version 1.0

workflow cumulus_tasks {
}

task run_cumulus_aggregate_matrices {
    input {
        File input_count_matrix_csv
        String output_directory
        String output_name
        String pegasus_version
        String docker_registry
        String zones
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
        String? restrictions
        String? attributes
        String? default_reference
        Boolean? select_only_singlets
        String? remap_singlets
        String? subset_singlets
        Int? minimum_number_of_genes
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        python <<CODE
        from subprocess import check_call

        call_args = ['pegasus', 'aggregate_matrix', '~{input_count_matrix_csv}', '~{output_name}.aggr']
        if '~{restrictions}' != '':
            ress = '~{restrictions}'.split(';')
            for res in ress:
                call_args.extend(['--restriction', res])
        if '~{attributes}' != '':
            call_args.extend(['--attributes', '~{attributes}'])
        if '~{default_reference}' != '':
            call_args.extend(['--default-reference', '~{default_reference}'])
        if '~{select_only_singlets}' == 'true':
            call_args.append('--select-only-singlets')
        if '~{remap_singlets}' != '':
            call_args.extend(['--remap-singlets', '~{remap_singlets}'])
        if '~{subset_singlets}' != '':
            call_args.extend(['--subset-singlets', '~{subset_singlets}'])
        if '~{minimum_number_of_genes}' != '':
            call_args.extend(['--min-genes', '~{minimum_number_of_genes}'])

        print(' '.join(call_args))
        check_call(call_args)

        output_directory = '~{output_directory}'
        if output_directory != '':
            dest = output_directory + '/' + '~{output_name}' + '/'
            call_args = ['strato', 'cp', '--backend', '~{backend}', '~{output_name}.aggr.zarr.zip', dest]
            print(' '.join(call_args))
            check_call(call_args)
        CODE
    }

    output {
        File output_zarr = "~{output_name}.aggr.zarr.zip"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cumulus:~{pegasus_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}

task run_cumulus_cluster {
    input {
        File input_file
        String output_directory
        String output_name
        String pegasus_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
        String? channel
        String? black_list
        Int? min_genes_before_filtration
        Boolean? select_singlets
        String? remap_singlets
        String? subset_singlets
        String? genome
        String? focus
        String? append
        Boolean? output_filtration_results
        Boolean? plot_filtration_results
        String? plot_filtration_figsize
        Boolean? output_h5ad
        Boolean? output_loom
        Int? min_genes
        Int? max_genes
        Int? min_umis
        Int? max_umis
        String? mito_prefix
        Float? percent_mito
        Float? gene_percent_cells
        Float? counts_per_cell_after
        String? select_hvf_flavor
        Int? select_hvf_ngenes
        Boolean? no_select_hvf
        Boolean? plot_hvf
        Boolean? correct_batch_effect
        String? correction_method
        String? batch_group_by
        Int? inmf_lambda
        Int? random_state
        String? gene_signature_set
        Int? nPC
        Boolean? run_nmf
        Int? nmf_n
        Int? knn_K
        Boolean? knn_full_speed
        Boolean? run_diffmap
        Int? diffmap_ndc
        Int? diffmap_maxt
        Boolean? run_louvain
        Float? louvain_resolution
        String? louvain_class_label
        Boolean? run_leiden
        Float? leiden_resolution
        Int? leiden_niter
        String? leiden_class_label
        Boolean? run_spectral_louvain
        String? spectral_louvain_basis
        Float? spectral_louvain_resolution
        String? spectral_louvain_class_label
        Boolean? run_spectral_leiden
        String? spectral_leiden_basis
        Float? spectral_leiden_resolution
        String? spectral_leiden_class_label
        Boolean? run_tsne
        Float? tsne_perplexity
        String? tsne_initialization
        Boolean? run_umap
        Int? umap_K
        Float? umap_min_dist
        Float? umap_spread
        Boolean? run_fle
        Int? fle_K
        Float? fle_target_change_per_node
        Int? fle_target_steps
        Float? net_down_sample_fraction
        Boolean? run_net_umap
        String? net_umap_out_basis
        Boolean? run_net_fle
        String? net_fle_out_basis
        Boolean? infer_doublets
        Float? expected_doublet_rate
        String? doublet_cluster_attribute
        Boolean? citeseq
        Boolean? citeseq_umap
        String? citeseq_umap_exclude
    }

    Boolean is_url = defined(gene_signature_set) && sub(select_first([gene_signature_set]), "^gs://.+", "URL") == "URL"

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        if [ ~{is_url} == true ]; then
            strato cp --backend ~{backend} ~{gene_signature_set} gene_signature.gmt
        fi

        python <<CODE
        from subprocess import check_call

        call_args = ['pegasus', 'cluster', '~{input_file}', '~{output_name}', '-p', '~{num_cpu}']
        if '~{channel}' != '':
            call_args.extend(['--channel', '~{channel}'])
        if '~{black_list}' != '':
            call_args.extend(['--black-list', '~{black_list}'])
        if '~{min_genes_before_filtration}' != '':
            call_args.extend(['--min-genes-before-filtration', '~{min_genes_before_filtration}'])
        if '~{select_singlets}' == 'true':
            call_args.append('--select-singlets')
        if '~{remap_singlets}' != '':
            call_args.extend(['--remap-singlets', '~{remap_singlets}'])
        if '~{subset_singlets}' != '':
            call_args.extend(['--subset-singlets', '~{subset_singlets}'])
        if '~{genome}' != '':
            call_args.extend(['--genome', '~{genome}'])
        if '~{focus}' != '':
            call_args.extend(['--focus', '~{focus}'])
        if '~{append}' != '':
            call_args.extend(['--append', '~{append}'])
        if '~{output_filtration_results}' == 'true':
            call_args.append('--output-filtration-results')
        if '~{plot_filtration_results}' == 'true':
            call_args.append('--plot-filtration-results')
        if '~{plot_filtration_figsize}' != '':
            call_args.extend(['--plot-filtration-figsize', '~{plot_filtration_figsize}'])
        if '~{output_h5ad}' == 'true':
            call_args.append('--output-h5ad')
        if '~{output_loom}' == 'true':
            call_args.append('--output-loom')
        if '~{correct_batch_effect}' == 'true':
            call_args.append('--correct-batch-effect')
            if '~{correction_method}' != '':
                call_args.extend(['--correction-method', '~{correction_method}'])
            if '~{batch_group_by}' != '':
                call_args.extend(['--batch-group-by', '~{batch_group_by}'])
            if '~{correction_method}' == 'inmf' and '~{inmf_lambda}' != '':
                call_args.extend(['--inmf-lambda', '~{inmf_lambda}'])
        if '~{min_genes}' != '':
            call_args.extend(['--min-genes', '~{min_genes}'])
        if '~{max_genes}' != '':
            call_args.extend(['--max-genes', '~{max_genes}'])
        if '~{min_umis}' != '':
            call_args.extend(['--min-umis', '~{min_umis}'])
        if '~{max_umis}' != '':
            call_args.extend(['--max-umis', '~{max_umis}'])
        if '~{mito_prefix}' != '':
            call_args.extend(['--mito-prefix', '~{mito_prefix}'])
        if '~{percent_mito}' != '' :
            call_args.extend(['--percent-mito', '~{percent_mito}'])
        if '~{gene_percent_cells}' != '':
            call_args.extend(['--gene-percent-cells', '~{gene_percent_cells}'])
        if '~{counts_per_cell_after}' != '':
            call_args.extend(['--counts-per-cell-after', '~{counts_per_cell_after}'])
        if '~{random_state}' != '':
            call_args.extend(['--random-state', '~{random_state}'])
        if '~{gene_signature_set}' != '':
            if '~{is_url}' == 'true':
                call_args.extend(['--calc-signature-scores', 'gene_signature.gmt'])
            else:
                call_args.extend(['--calc-signature-scores', '~{gene_signature_set}'])
        if '~{no_select_hvf}' == 'true':
            call_args.append('--no-select-hvf')
        if '~{select_hvf_flavor}' != '':
            call_args.extend(['--select-hvf-flavor', '~{select_hvf_flavor}'])
        if '~{select_hvf_ngenes}' != '':
            call_args.extend(['--select-hvf-ngenes', '~{select_hvf_ngenes}'])
        if '~{plot_hvf}' == 'true':
            call_args.append('--plot-hvf')
        if '~{nPC}' != '':
            call_args.extend(['--pca-n', '~{nPC}'])
        if '~{run_nmf}' != '':
            call_args.append('--nmf')
        if '~{nmf_n}' != '':
            call_args.extend(['--nmf-n', '~{nmf_n}'])
        if '~{knn_K}' != '':
            call_args.extend(['--knn-K', '~{knn_K}'])
        if '~{knn_full_speed}' == 'true':
            call_args.append('--knn-full-speed')
        if '~{run_diffmap}' == 'true':
            call_args.append('--diffmap')
        if '~{diffmap_ndc}' != '':
            call_args.extend(['--diffmap-ndc', '~{diffmap_ndc}'])
        if '~{diffmap_maxt}' != '':
            call_args.extend(['--diffmap-maxt', '~{diffmap_maxt}'])
        if '~{run_louvain}' == 'true':
            call_args.append('--louvain')
        if '~{louvain_resolution}' != '':
            call_args.extend(['--louvain-resolution', '~{louvain_resolution}'])
        if '~{louvain_class_label}' != '':
            call_args.extend(['--louvain-class-label', '~{louvain_class_label}'])
        if '~{run_leiden}' == 'true':
            call_args.append('--leiden')
        if '~{leiden_resolution}' != '':
            call_args.extend(['--leiden-resolution', '~{leiden_resolution}'])
        if '~{leiden_niter}' != '':
            call_args.extend(['--leiden-niter', '~{leiden_niter}'])
        if '~{leiden_class_label}' != '':
            call_args.extend(['--leiden-class-label', '~{leiden_class_label}'])
        if '~{run_spectral_louvain}' == 'true':
            call_args.append('--spectral-louvain')
        if '~{spectral_louvain_basis}' != '':
            call_args.extend(['--spectral-louvain-basis', '~{spectral_louvain_basis}'])
        if '~{spectral_louvain_resolution}' != '':
            call_args.extend(['--spectral-louvain-resolution', '~{spectral_louvain_resolution}'])
        if '~{spectral_louvain_class_label}' != '':
            call_args.extend(['--spectral-louvain-class-label', '~{spectral_louvain_class_label}'])
        if '~{run_spectral_leiden}' == 'true':
            call_args.append('--spectral-leiden')
        if '~{spectral_leiden_basis}' != '':
            call_args.extend(['--spectral-leiden-basis', '~{spectral_leiden_basis}'])
        if '~{spectral_leiden_resolution}' != '':
            call_args.extend(['--spectral-leiden-resolution', '~{spectral_leiden_resolution}'])
        if '~{spectral_leiden_class_label}' != '':
            call_args.extend(['--spectral-leiden-class-label', '~{spectral_leiden_class_label}'])
        if '~{run_tsne}' == 'true':
            call_args.append('--tsne')
        if '~{tsne_perplexity}' != '':
            call_args.extend(['--tsne-perplexity', '~{tsne_perplexity}'])
        if '~{tsne_initialization}' != '':
            call_args.extend(['--tsne-initialization', '~{tsne_initialization}'])
        if '~{run_umap}' == 'true':
            call_args.append('--umap')
        if '~{umap_K}' != '':
            call_args.extend(['--umap-K', '~{umap_K}'])
        if '~{umap_min_dist}' != '':
            call_args.extend(['--umap-min-dist', '~{umap_min_dist}'])
        if '~{umap_spread}' != '':
            call_args.extend(['--umap-spread', '~{umap_spread}'])
        if '~{run_fle}' == 'true':
            call_args.append('--fle')
        if '~{fle_K}' != '':
            call_args.extend(['--fle-K', '~{fle_K}'])
        if '~{fle_target_change_per_node}' != '':
            call_args.extend(['--fle-target-change-per-node', '~{fle_target_change_per_node}'])
        if '~{fle_target_steps}' != '':
            call_args.extend(['--fle-target-steps', '~{fle_target_steps}'])
        if '~{net_down_sample_fraction}' != '':
            call_args.extend(['--net-down-sample-fraction', '~{net_down_sample_fraction}'])
        if '~{run_net_umap}' == 'true':
            call_args.append('--net-umap')
        if '~{net_umap_out_basis}' != '':
            call_args.extend(['--net-umap-out-basis', '~{net_umap_out_basis}'])
        if '~{run_net_fle}' == 'true':
            call_args.append('--net-fle')
        if '~{net_fle_out_basis}' != '':
            call_args.extend(['--net-fle-out-basis', '~{net_fle_out_basis}'])
        if '~{infer_doublets}' == 'true':
            call_args.append('--infer-doublets')
        if '~{expected_doublet_rate}' != '':
            call_args.extend(['--expected-doublet-rate', '~{expected_doublet_rate}'])
        if '~{doublet_cluster_attribute}' != '':
            call_args.extend(['--dbl-cluster-attr', '~{doublet_cluster_attribute}'])
        if '~{citeseq}' == 'true':
            call_args.append('--citeseq')
        if '~{citeseq_umap}' == 'true':
            call_args.append('--citeseq-umap')
        if '~{citeseq_umap_exclude}' != '':
            call_args.extend(['--citeseq-umap-exclude', '~{citeseq_umap_exclude}'])
        print(' '.join(call_args))
        check_call(call_args)

        import glob
        output_directory = '~{output_directory}'
        if output_directory != '':
            dest = output_directory + '/' + '~{output_name}' + '/'
            files = ['~{output_name}.zarr.zip', '~{output_name}.log']
            if ('~{output_h5ad}' == 'true') or ('~{citeseq}' == 'true'):
                files.extend(glob.glob('~{output_name}.*.h5ad'))
            if '~{output_filtration_results}' == 'true':
                files.extend(glob.glob('~{output_name}.*.filt.xlsx'))
            if ('~{no_select_hvf}' != 'true') and ('~{plot_hvf}' == 'true'):
                files.extend(glob.glob('~{output_name}.*.hvf.pdf'))
            if '~{plot_filtration_results}' == 'true':
                files.extend(glob.glob('~{output_name}.*.filt.*.pdf'))
            if '~{infer_doublets}' == 'true':
                files.extend(glob.glob('~{output_name}.*.dbl.png'))
            if '~{output_loom}' == 'true':
                files.extend(glob.glob('~{output_name}.*.loom'))
            for file in files:
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', file, dest]
                print(' '.join(call_args))
                check_call(call_args)
        CODE
    }

    output {
        File output_zarr = "~{output_name}.zarr.zip"
        Array[File] output_h5ad_file = glob("~{output_name}.*.h5ad")
        Array[File] output_filt_xlsx = glob("~{output_name}.*.filt.xlsx")
        Array[File] output_filt_plot = glob("~{output_name}.*.filt.*.pdf")
        Array[File] output_hvf_plot = glob("~{output_name}.*.hvf.pdf")
        Array[File] output_dbl_plot = glob("~{output_name}.*.dbl.png")
        Array[File] output_loom_file = glob("~{output_name}.*.loom")
        File output_log = "~{output_name}.log"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cumulus:~{pegasus_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}

task run_cumulus_cirro_output {
    input {
        File input_h5ad
        String output_directory
        String output_name
        String pegasus_version
        String docker_registry
        String zones
        String memory
        Int disk_space
        Int num_cpu
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        cirro prepare_data --out "~{output_name}".cirro ~{input_h5ad}

        python <<CODE

        from subprocess import check_call

        call_args=['strato', 'cp', '--backend', '~{backend}', '-m', '-r', "~{output_name}".cirro]

        if output_directory != '':
            output_directory=~{output_directory}+'/'
            call_args.append(output_directory)
            print(' '.join(call_args))
            check_call(call_args)
        CODE

    }

    output {
        String output_cirro_folder = "~{output_directory}/~{output_name}.cirro"
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cumulus:~{pegasus_version}"
        zones: zones
        memory: memory
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}

task run_cumulus_de_analysis {
    input {
        File input_h5ad
        String output_directory
        String output_name
        String pegasus_version
        String docker_registry
        String zones
        Int num_cpu
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
        String? labels
        Boolean? t_test
        Boolean? fisher
        Float? alpha

        Boolean? annotate_cluster
        String? annotate_de_test
        String? organism
        Float? minimum_report_score

        Boolean? find_markers_lightgbm
        Boolean? remove_ribo
        Float? min_gain
        Int? random_state
    }

    Boolean is_url = defined(organism) && sub(select_first([organism]), "^.+\\.json", "URL") == "URL"

    command {
        set -e
        export TMPDIR=/tmp
        export BACKEND=~{backend}
        monitor_script.sh > monitoring.log &

        if [ ~{is_url} == true ]; then
            strato cp --backend ~{backend} ~{organism} markers.json
        fi

        python <<CODE
        from subprocess import check_call

        call_args = ['mv', '-f', '~{input_h5ad}', '~{output_name}.h5ad']
        print(' '.join(call_args))
        check_call(call_args)
        call_args = ['pegasus', 'de_analysis', '~{output_name}.h5ad', '~{output_name}.de.xlsx', '-p', '~{num_cpu}']
        if '~{labels}' == '':
            call_args.extend(['--labels', 'louvain_labels'])
        else:
            call_args.extend(['--labels', '~{labels}'])
        if '~{alpha}' != '':
            call_args.extend(['--alpha', '~{alpha}'])
        if '~{t_test}' == 'true':
            call_args.append('--t')
        if '~{fisher}' == 'true':
            call_args.append('--fisher')
        print(' '.join(call_args))
        check_call(call_args)
        if '~{find_markers_lightgbm}' == 'true':
            call_args = ['pegasus', 'find_markers', '~{output_name}.h5ad', '~{output_name}.markers.xlsx', '-p', '~{num_cpu}']
            if '~{labels}' != '':
                call_args.extend(['--labels', '~{labels}'])
            if '~{remove_ribo}' == 'true':
                call_args.append('--remove-ribo')
            if '~{min_gain}' != '':
                call_args.extend(['--min-gain', '~{min_gain}'])
            if '~{random_state}' != '':
                call_args.extend(['--random-state', '~{random_state}'])
            print(' '.join(call_args))
            check_call(call_args)
        if '~{annotate_cluster}' == 'true':
            call_args = ['pegasus', 'annotate_cluster', '~{output_name}.h5ad', '~{output_name}.anno.txt']
            if '~{organism}' != '':
                if '~{is_url}' == 'true':
                    call_args.extend(['--markers', 'markers.json'])
                else:
                    call_args.extend(['--markers', '~{organism}'])
            if '~{annotate_de_test}' != '':
                call_args.extend(['--de-test', '~{annotate_de_test}'])
            if '~{alpha}' != '':
                call_args.extend(['--de-alpha', '~{alpha}'])
            if '~{minimum_report_score}' != '':
                call_args.extend(['--minimum-report-score', '~{minimum_report_score}'])
            print(' '.join(call_args))
            check_call(call_args)

        output_directory = '~{output_directory}'
        if output_directory != '':
            dest = output_directory + '/'
            files = ['~{output_name}.h5ad', '~{output_name}.de.xlsx']
            if '~{find_markers_lightgbm}' == 'true':
                files.append('~{output_name}.markers.xlsx')
            if '~{annotate_cluster}' == 'true':
                files.append('~{output_name}.anno.txt')
            for file in files:
                call_args = ['strato', 'cp', '--backend', '~{backend}', file, dest]
                print(' '.join(call_args))
                check_call(call_args)
        CODE
    }

    output {
        File output_de_h5ad = "~{output_name}.h5ad"
        File output_de_xlsx = "~{output_name}.de.xlsx"
        Array[File] output_markers_xlsx = glob("~{output_name}.markers.xlsx")
        Array[File] output_anno_file = glob("~{output_name}.anno.txt")
        File monitoringLog = "monitoring.log"
    }

    runtime {
        docker: "~{docker_registry}/cumulus:~{pegasus_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: num_cpu
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}

task run_cumulus_plot {
    input {
        File input_h5ad
        String output_directory
        String output_name
        String pegasus_version
        String docker_registry
        String zones
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
        String? plot_composition
        String? plot_tsne
        String? plot_umap
        String? plot_fle
        String? plot_net_umap
        String? plot_net_fle
        String? plot_citeseq_umap
        String? plot_nmf
        Int nmf_n
    }

    command {
        set -e
        export TMPDIR=/tmp

        python <<CODE
        from subprocess import check_call
        if '~{plot_composition}' != '':
            pairs = '~{plot_composition}'.split(',')
            for pair in pairs:
                lab, attr = pair.split(':')
                call_args = ['pegasus', 'plot', 'compo', '--groupby', lab, '--condition', attr, '--style', 'normalized', '~{input_h5ad}', '~{output_name}.' + lab + '.' + attr + '.composition.pdf']
                print(' '.join(call_args))
                check_call(call_args)
        if '~{plot_tsne}' != '':
            call_args = ['pegasus', 'plot', 'scatter',  '--basis', 'tsne', '--attributes', '~{plot_tsne}', '~{input_h5ad}', '~{output_name}.tsne.pdf']
            print(' '.join(call_args))
            check_call(call_args)
        if '~{plot_umap}' != '':
            call_args = ['pegasus', 'plot', 'scatter', '--basis', 'umap', '--attributes', '~{plot_umap}', '~{input_h5ad}', '~{output_name}.umap.pdf']
            print(' '.join(call_args))
            check_call(call_args)
        if '~{plot_fle}' != '':
            call_args = ['pegasus', 'plot', 'scatter', '--basis', 'fle', '--attributes', '~{plot_fle}', '~{input_h5ad}', '~{output_name}.fle.pdf']
            print(' '.join(call_args))
            check_call(call_args)
        if '~{plot_net_umap}' != '':
            call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_umap', '--attributes', '~{plot_net_umap}', '~{input_h5ad}', '~{output_name}.net.umap.pdf']
            print(' '.join(call_args))
            check_call(call_args)
        if '~{plot_net_fle}' != '':
            call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_fle', '--attributes', '~{plot_net_fle}', '~{input_h5ad}', '~{output_name}.net.fle.pdf']
            print(' '.join(call_args))
            check_call(call_args)
        if '~{plot_citeseq_umap}' != '':
            call_args = ['pegasus', 'plot', 'scatter', '--basis', 'citeseq_umap', '--attributes', '~{plot_citeseq_umap}', '~{input_h5ad}', '~{output_name}.citeseq.umap.pdf']
            print(' '.join(call_args))
            check_call(call_args)
        if '~{plot_nmf}' != '':
            for i in range(int('~{nmf_n}')):
                cur_factor = 'X_' + '~{plot_nmf}' + '@' + str(i)
                call_args = ['pegasus', 'plot', 'scatter', '--basis', 'umap', '--attributes', cur_factor, '~{input_h5ad}', '~{output_name}.~{plot_nmf}_'+str(i)+'.umap.pdf']
                print(' '.join(call_args))
                check_call(call_args)

                call_args = ['pegasus', 'plot', 'wordcloud', '--factor', str(i), '~{input_h5ad}', '~{output_name}.word_cloud_'+str(i)+'.pdf']
                print(' '.join(call_args))
                check_call(call_args)

        import glob
        output_directory = '~{output_directory}'
        if output_directory != '':
            dest = output_directory + '/'
            files = glob.glob('*.pdf')
            files.extend(glob.glob('*.html'))

            for file in files:
                call_args = ['strato', 'cp', '--backend', '~{backend}', '-m', file, dest]
                print(' '.join(call_args))
                check_call(call_args)
        CODE
    }

    output {
        Array[File] output_pdfs = glob("*.pdf")
        Array[File] output_htmls = glob("*.html")
    }

    runtime {
        docker: "~{docker_registry}/cumulus:~{pegasus_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}

task run_cumulus_scp_output {
    input {
        File input_h5ad
        String output_directory
        String output_name
        Boolean output_dense
        String pegasus_version
        String docker_registry
        String zones
        String memory
        Int disk_space
        Int preemptible
        Int awsMaxRetries
        String backend
    }

    command {
        set -e
        export TMPDIR=/tmp
        pegasus scp_output ~{true='--dense' false='' output_dense} ~{input_h5ad} "~{output_name}"

        python <<CODE

        from subprocess import check_call

        call_args=['strato', 'cp', '--backend', '~{backend}']

        if output_directory != '':
            output_directory='~{output_directory}' + '/'
            call_args.extend(['~{output_name}.scp.*', output_directory])
            print(' '.join(call_args))
            check_call(call_args)

        CODE
    }

    output {
        Array[File] output_scp_files = glob("~{output_name}.scp.*")
    }

    runtime {
        docker: "~{docker_registry}/cumulus:~{pegasus_version}"
        zones: zones
        memory: memory
        bootDiskSizeGb: 12
        disks: "local-disk ~{disk_space} HDD"
        cpu: 1
        preemptible: preemptible
        maxRetries: if backend == "aws" then awsMaxRetries else 0
    }
}
