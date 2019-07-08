workflow scCloud_tasks {
}

task run_scCloud_aggregate_matrices {
	File input_count_matrix_csv
	String output_name
	String sccloud_version
	String zones
	String memory
	Int disk_space
	Int preemptible
	String? restrictions
	String? attributes
	Boolean? select_only_singlets
	Int? minimum_number_of_genes
	String? dropseq_genome

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call
		call_args = ['scCloud', 'aggregate_matrix', '${input_count_matrix_csv}', '${output_name}', '--google-cloud']
		# call_args = ['scCloud', 'aggregate_matrix', '${input_count_matrix_csv}', '${output_name}']
		if '${restrictions}' is not '':
			ress = '${restrictions}'.split(';')
			for res in ress:
				call_args.extend(['--restriction', res])
		if '${attributes}' is not '':
			call_args.extend(['--attributes', '${attributes}'])
		if '${select_only_singlets}' is 'true':
			call_args.append('--select-only-singlets')
		if '${minimum_number_of_genes}' is not '':
			call_args.extend(['--minimum-number-of-genes', '${minimum_number_of_genes}'])
		if '${dropseq_genome}' is not '':
			call_args.extend(['--dropseq-genome', '${dropseq_genome}'])

		print(' '.join(call_args))
		check_call(call_args)
		CODE
	}

	output {
		File output_10x_h5 = '${output_name}_10x.h5'
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_scCloud_cluster {
	File input_10x_file
	String output_name
	String sccloud_version
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	String? genome
	Boolean? cite_seq
	Float? cite_seq_capping
	Boolean? output_filtration_results
	Boolean? plot_filtration_results
	String? plot_filtration_figsize
	Boolean? output_seurat_compatible
	Boolean? output_loom
	Boolean? output_parquet
	Boolean? correct_batch_effect
	String? batch_group_by
	Int? min_genes
	Int? max_genes
	Int? min_umis
	Int? max_umis
	String? mito_prefix
	Float? percent_mito
	Float? gene_percent_cells
	Int? min_genes_on_raw
	Float? counts_per_cell_after
	Int? random_state
	Boolean? run_uncentered_pca
	Boolean? no_variable_gene_selection
	Boolean? no_submat_to_dense
	Int? nPC
	Int? nDC
	Float? diffmap_alpha
	Int? diffmap_K
	Boolean? diffmap_full_speed
	Boolean? run_louvain
	Float? louvain_resolution
	String? louvain_class_label
	Boolean? run_leiden
	Float? leiden_resolution
	String? leiden_class_label
	Boolean? run_approximated_louvain
	String? approx_louvain_basis
	Float? approx_louvain_resolution
	String? approx_louvain_class_label
	Boolean? run_approximated_leiden
	String? approx_leiden_basis
	Float? approx_leiden_resolution
	String? approx_leiden_class_label
	Boolean? run_tsne
	Boolean? run_fitsne
	Float? tsne_perplexity
	Boolean? run_umap
	Int? umap_K
	Float? umap_min_dist
	Float? umap_spread
	Boolean? run_fle
	Int? fle_K
	Float? fle_target_change_per_node
	Int? fle_target_steps
	Boolean? fle_3D
	Float? net_down_sample_fraction
	Boolean? net_ds_full_speed
	Boolean? run_net_tsne
	String? net_tsne_out_basis
	Boolean? run_net_fitsne
	String? net_fitsne_out_basis
	Boolean? run_net_umap
	String? net_umap_out_basis
	Boolean? run_net_fle
	Boolean? net_fle_ds_full_speed
	String? net_fle_out_basis

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['scCloud', 'cluster', '${input_10x_file}', '${output_name}', '-p', '${num_cpu}']
		if '${genome}' is not '':
			call_args.extend(['--genome', '${genome}'])
		if '${cite_seq}' is 'true':
			call_args.append('--cite-seq')
		if '${cite_seq_capping}' is not '':
			call_args.extend(['--cite-seq-capping', '${cite_seq_capping}'])
		if '${output_filtration_results}' is 'true':
			call_args.append('--output-filtration-results')
		if '${plot_filtration_results}' is 'true':
			call_args.append('--plot-filtration-results')
		if '${plot_filtration_figsize}' is not '':
			call_args.extend(['--plot-filtration-figsize', '${plot_filtration_figsize}'])
		if '${output_seurat_compatible}' is 'true':
			call_args.append('--output-seurat-compatible')
		if '${output_loom}' is 'true':
			call_args.append('--output-loom')
		if '${correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
			if '${batch_group_by}' is not '':
				call_args.extend(['--batch-group-by', '${batch_group_by}'])
		if '${min_genes}' is not '':
			call_args.extend(['--min-genes', '${min_genes}'])
		if '${max_genes}' is not '':
			call_args.extend(['--max-genes', '${max_genes}'])
		if '${min_umis}' is not '':
			call_args.extend(['--min-umis', '${min_umis}'])
		if '${max_umis}' is not '':
			call_args.extend(['--max-umis', '${max_umis}'])
		if '${mito_prefix}' is not '':
			call_args.extend(['--mito-prefix', '${mito_prefix}'])
		if '${percent_mito}' is not '' :
			call_args.extend(['--percent-mito', '${percent_mito}'])
		if '${gene_percent_cells}' is not '':
			call_args.extend(['--gene-percent-cells', '${gene_percent_cells}'])
		if '${min_genes_on_raw}' is not '':
			call_args.extend(['--min-genes-on-raw', '${min_genes_on_raw}'])
		if '${counts_per_cell_after}' is not '':
			call_args.extend(['--counts-per-cell-after', '${counts_per_cell_after}'])
		if '${random_state}' is not '':
			call_args.extend(['--random-state', '${random_state}'])
		if '${run_uncentered_pca}' is 'true':
			call_args.append('--run-uncentered-pca')
		if '${no_variable_gene_selection}' is 'true':
			call_args.append('--no-variable-gene-selection')
		if '${no_submat_to_dense}' is 'true':
			call_args.append('--no-submat-to-dense')
		if '${nPC}' is not '':
			call_args.extend(['--nPC', '${nPC}'])
		if '${nDC}' is not '':
			call_args.extend(['--nDC', '${nDC}'])
		if '${diffmap_alpha}' is not '':
			call_args.extend(['--diffmap-alpha', '${diffmap_alpha}'])
		if '${diffmap_K}' is not '':
			call_args.extend(['--diffmap-K', '${diffmap_K}'])
		if '${diffmap_full_speed}' is 'true':
			call_args.append('--diffmap-full-speed')
		if '${run_louvain}' is 'true':
			call_args.append('--run-louvain')
		if '${louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '${louvain_resolution}'])
		if '${louvain_class_label}' is not '':
			call_args.extend(['--louvain-class-label', '${louvain_class_label}'])
		if '${run_leiden}' is 'true':
			call_args.append('--run-leiden')
		if '${leiden_resolution}' is not '':
			call_args.extend(['--leiden-resolution', '${leiden_resolution}'])
		if '${leiden_class_label}' is not '':
			call_args.extend(['--leiden-class-label', '${leiden_class_label}'])
		if '${run_approximated_louvain}' is 'true':
			call_args.append('--run-approximated-louvain')
		if '${approx_louvain_basis}' is not '':
			call_args.extend(['--approx-louvain-basis', '${approx_louvain_basis}'])
		if '${approx_louvain_resolution}' is not '':
			call_args.extend(['--approx-louvain-resolution', '${approx_louvain_resolution}'])
		if '${approx_louvain_class_label}' is not '':
			call_args.extend(['--approx-louvain-class-label', '${approx_louvain_class_label}'])
		if '${run_approximated_leiden}' is 'true':
			call_args.append('--run-approximated-leiden')
		if '${approx_leiden_basis}' is not '':
			call_args.extend(['--approx-leiden-basis', '${approx_leiden_basis}'])
		if '${approx_leiden_resolution}' is not '':
			call_args.extend(['--approx-leiden-resolution', '${approx_leiden_resolution}'])
		if '${approx_leiden_class_label}' is not '':
			call_args.extend(['--approx-leiden-class-label', '${approx_leiden_class_label}'])
		if '${run_tsne}' is 'true':
			call_args.append('--run-tsne')
		if '${run_fitsne}' is 'true':
			call_args.append('--run-fitsne')
		if '${tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '${tsne_perplexity}'])
		if '${run_umap}' is 'true':
			call_args.append('--run-umap')
		if '${umap_K}' is not '':
			call_args.extend(['--umap-K', '${umap_K}'])
		if '${umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '${umap_min_dist}'])
		if '${umap_spread}' is not '':
			call_args.extend(['--umap-spread', '${umap_spread}'])
		if '${run_fle}' is 'true':
			call_args.append('--run-fle')
		if '${fle_K}' is not '':
			call_args.extend(['--fle-K', '${fle_K}'])
		if '${fle_target_change_per_node}' is not '':
			call_args.extend(['--fle-target-change-per-node', '${fle_target_change_per_node}'])
		if '${fle_target_steps}' is not '':
			call_args.extend(['--fle-target-steps', '${fle_target_steps}'])
		if '${fle_3D}' is 'true':
			call_args.append('--fle-3D')
		if '${net_down_sample_fraction}' is not '':
			call_args.extend(['--net-down-sample-fraction', '${net_down_sample_fraction}'])
		if '${net_ds_full_speed}' is 'true':
			call_args.append('--net-ds-full-speed')
		if '${run_net_tsne}' is 'true':
			call_args.append('--run-net-tsne')
		if '${net_tsne_out_basis}' is not '':
			call_args.extend(['--net-tsne-out-basis', '${net_tsne_out_basis}'])
		if '${run_net_fitsne}' is 'true':
			call_args.append('--run-net-fitsne')
		if '${net_fitsne_out_basis}' is not '':
			call_args.extend(['--net-fitsne-out-basis', '${net_fitsne_out_basis}'])
		if '${run_net_umap}' is 'true':
			call_args.append('--run-net-umap')
		if '${net_umap_out_basis}' is not '':
			call_args.extend(['--net-umap-out-basis', '${net_umap_out_basis}'])
		if '${run_net_fle}' is 'true':
			call_args.append('--run-net-fle')
		if '${net_fle_ds_full_speed}' is 'true':
			call_args.append('--net-fle-ds-full-speed')
		if '${net_fle_out_basis}' is not '':
			call_args.extend(['--net-fle-out-basis', '${net_fle_out_basis}'])
		print(' '.join(call_args))
		check_call(call_args)
		if '${output_parquet}' is 'true':
			call_args = ['scCloud', 'parquet', '${output_name}.h5ad', '${output_name}', '-p', '${num_cpu}']
			print(' '.join(call_args))
			check_call(call_args)			
		CODE
	}

	output {
		File output_h5ad = "${output_name}.h5ad"
		Array[File] output_seurat_h5ad = glob("${output_name}.seurat.h5ad")
		Array[File] output_filt_xlsx = glob("${output_name}.filt.xlsx")
		Array[File] output_filt_plot = glob("${output_name}.filt.*.pdf")
		Array[File] output_loom_file = glob("${output_name}.loom")
		Array[File] output_parquet_file = glob("${output_name}.parquet")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_scCloud_de_analysis {
	File input_h5ad
	String output_name
	String sccloud_version
	String zones	
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	String? labels
	Float? alpha
	Boolean? fisher
	Boolean? mwu
	Boolean? roc

	Boolean? annotate_cluster
	String? organism
	Float? minimum_report_score

	Boolean? find_markers_lightgbm
	Boolean? remove_ribo
	Float? min_gain
	Int? random_state

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['mv', '-f', '${input_h5ad}', '${output_name}.h5ad']
		print(' '.join(call_args))
		check_call(call_args)			
		call_args = ['scCloud', 'de_analysis', '${output_name}.h5ad', '${output_name}.de.xlsx', '-p', '${num_cpu}']
		if '${labels}' is not '':
			call_args.extend(['--labels', '${labels}'])
		if '${alpha}' is not '':
			call_args.extend(['--alpha', '${alpha}'])
		if '${fisher}' is 'true':
			call_args.append('--fisher')
		if '${mwu}' is 'true':
			call_args.append('--mwu')
		if '${roc}' is 'true':
			call_args.append('--roc')
		print(' '.join(call_args))
		check_call(call_args)
		if '${find_markers_lightgbm}' is 'true':
			call_args = ['scCloud', 'find_markers', '${output_name}.h5ad', '${output_name}.markers.xlsx', '-p', '${num_cpu}']
			if '${labels}' is not '':
				call_args.extend(['--labels', '${labels}'])
			if '${remove_ribo}' is 'true':
				call_args.append('--remove-ribo')
			if '${min_gain}' is not '':
				call_args.extend(['--min-gain', '${min_gain}'])
			if '${random_state}' is not '':
				call_args.extend(['--random-state', '${random_state}'])
			print(' '.join(call_args))
			check_call(call_args)
		if '${annotate_cluster}' is 'true':
			call_args = ['scCloud', 'annotate_cluster', '${output_name}.h5ad', '${output_name}' + '.anno.txt']
			if '${organism}' is not '':
				call_args.extend(['--json-file', '${organism}'])
			if '${minimum_report_score}' is not '':
				call_args.extend(['--minimum-report-score', '${minimum_report_score}'])
			print(' '.join(call_args))
			check_call(call_args)			
		CODE
	}

	output {
		File output_de_h5ad = "${output_name}.h5ad"
		File output_de_xlsx = "${output_name}.de.xlsx"
		Array[File] output_markers_xlsx = glob("${output_name}.markers.xlsx")
		Array[File] output_anno_file = glob("${output_name}.anno.txt")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_scCloud_plot {
	File input_h5ad
	String output_name
	String sccloud_version
	String zones
	String memory
	Int disk_space
	Int preemptible
	String? plot_composition
	String? plot_tsne
	String? plot_fitsne
	String? plot_umap
	String? plot_fle
	String? plot_diffmap
	String? plot_citeseq_fitsne
	String? plot_net_tsne
	String? plot_net_fitsne
	String? plot_net_umap
	String? plot_net_fle

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call
		if '${plot_composition}' is not '':
			pairs = '${plot_composition}'.split(',')
			for pair in pairs:
				lab, attr = pair.split(':')
				call_args = ['scCloud', 'plot', 'composition', '--cluster-labels', lab, '--attribute', attr, '--style', 'normalized', '--not-stacked', '${input_h5ad}', '${output_name}.' + lab + '.' + attr + '.composition.pdf']
				print(' '.join(call_args))
				check_call(call_args)
		if '${plot_tsne}' is not '':
			call_args = ['scCloud', 'plot', 'scatter',  '--basis', 'tsne', '--attributes', '${plot_tsne}', '${input_h5ad}', '${output_name}.tsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_fitsne}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'fitsne', '--attributes', '${plot_fitsne}', '${input_h5ad}', '${output_name}.fitsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_umap}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'umap', '--attributes', '${plot_umap}', '${input_h5ad}', '${output_name}.umap.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_fle}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'fle', '--attributes', '${plot_fle}', '${input_h5ad}', '${output_name}.fle.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_diffmap}' is not '':
			attrs = '${plot_diffmap}'.split(',')
			for attr in attrs:
				call_args = ['scCloud', 'iplot', '--attribute', attr, 'diffmap_pca', '${input_h5ad}', '${output_name}.' + attr + '.diffmap_pca.html']
				print(' '.join(call_args))
				check_call(call_args)
		if '${plot_citeseq_fitsne}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'citeseq_fitsne', '--attributes', '${plot_citeseq_fitsne}', '${input_h5ad}', '${output_name}.citeseq.fitsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_net_tsne}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'net_tsne', '--attributes', '${plot_net_tsne}', '${input_h5ad}', '${output_name}.net.tsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_net_fitsne}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'net_fitsne', '--attributes', '${plot_net_fitsne}', '${input_h5ad}', '${output_name}.net.fitsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_net_umap}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'net_umap', '--attributes', '${plot_net_umap}', '${input_h5ad}', '${output_name}.net.umap.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_net_fle}' is not '':
			call_args = ['scCloud', 'plot', 'scatter', '--basis', 'net_fle', '--attributes', '${plot_net_fle}', '${input_h5ad}', '${output_name}.net.fle.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		CODE
	}

	output {
		Array[File] output_pdfs = glob("*.pdf")
		Array[File] output_htmls = glob("*.html")
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_scCloud_scp_output {
	File input_h5ad
	String output_name
	Boolean output_dense
	String sccloud_version
	String zones
	String memory
	Int disk_space
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp
		scCloud scp_output ${true='--dense' false='' output_dense} ${input_h5ad} ${output_name}
	}

	output {
		Array[File] output_scp_files = glob("${output_name}.scp.*")
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_scCloud_subcluster {
	File input_h5ad
	String output_name
	String sccloud_version
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	String subset_selections 
	Boolean? correct_batch_effect
	Boolean? output_loom
	Boolean? output_parquet
	Int? random_state
	Boolean? run_uncentered_pca
	Boolean? no_variable_gene_selection
	Boolean? no_submat_to_dense
	Int? nPC
	Int? nDC
	Float? diffmap_alpha
	Int? diffmap_K
	Boolean? diffmap_full_speed
	String? calculate_pseudotime
	Boolean? run_louvain
	Float? louvain_resolution
	String? louvain_class_label
	Boolean? run_leiden
	Float? leiden_resolution
	String? leiden_class_label
	Boolean? run_approximated_louvain
	String? approx_louvain_basis
	Float? approx_louvain_resolution
	String? approx_louvain_class_label
	Boolean? run_approximated_leiden
	String? approx_leiden_basis
	Float? approx_leiden_resolution
	String approx_leiden_class_label
	Boolean? run_tsne
	Boolean? run_fitsne
	Float? tsne_perplexity
	Boolean? run_umap
	Int? umap_K
	Float? umap_min_dist
	Float? umap_spread
	Boolean? run_fle
	Int? fle_K
	Float? fle_target_change_per_node
	Int? fle_target_steps
	Boolean? fle_3D
	Float? net_down_sample_fraction
	Boolean? net_ds_full_speed
	Boolean? run_net_tsne
	String? net_tsne_out_basis
	Boolean? run_net_fitsne
	String? net_fitsne_out_basis
	Boolean? run_net_umap
	String? net_umap_out_basis
	Boolean? run_net_fle
	Boolean? net_fle_ds_full_speed
	String? net_fle_out_basis

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['scCloud', 'subcluster', '${input_h5ad}', '${output_name}', '-p', '${num_cpu}']
		if '${subset_selections}' is not '':
			sels = '${subset_selections}'.split(';')
			for sel in sels:
				call_args.extend(['--subset-selection', sel])
		if '${correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
		if '${output_loom}' is 'true':
			call_args.append('--output-loom')
		if '${random_state}' is not '':
			call_args.extend(['--random-state', '${random_state}'])
		if '${run_uncentered_pca}' is 'true':
			call_args.append('--run-uncentered-pca')
		if '${no_variable_gene_selection}' is 'true':
			call_args.append('--no-variable-gene-selection')
		if '${no_submat_to_dense}' is 'true':
			call_args.append('--no-submat-to-dense')
		if '${nPC}' is not '':
			call_args.extend(['--nPC', '${nPC}'])
		if '${nDC}' is not '':
			call_args.extend(['--nDC', '${nDC}'])
		if '${diffmap_alpha}' is not '':
			call_args.extend(['--diffmap-alpha', '${diffmap_alpha}'])
		if '${diffmap_K}' is not '':
			call_args.extend(['--diffmap-K', '${diffmap_K}'])
		if '${diffmap_full_speed}' is 'true':
			call_args.append('--diffmap-full-speed')
		if '${calculate_pseudotime}' is not '':
			call_args.extend(['--calculate-pseudotime', '${calculate_pseudotime}'])
		if '${run_louvain}' is 'true':
			call_args.append('--run-louvain')
		if '${louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '${louvain_resolution}'])
		if '${louvain_class_label}' is not '':
			call_args.extend(['--louvain-class-label', '${louvain_class_label}'])
		if '${run_leiden}' is 'true':
			call_args.append('--run-leiden')
		if '${leiden_resolution}' is not '':
			call_args.extend(['--leiden-resolution', '${leiden_resolution}'])
		if '${leiden_class_label}' is not '':
			call_args.extend(['--leiden-class-label', '${leiden_class_label}'])
		if '${run_approximated_louvain}' is 'true':
			call_args.append('--run-approximated-louvain')
		if '${approx_louvain_basis}' is not '':
			call_args.extend(['--approx-louvain-basis', '${approx_louvain_basis}'])
		if '${approx_louvain_resolution}' is not '':
			call_args.extend(['--approx-louvain-resolution', '${approx_louvain_resolution}'])
		if '${approx_louvain_class_label}' is not '':
			call_args.extend(['--approx-louvain-class-label', '${approx_louvain_class_label}'])
		if '${run_approximated_leiden}' is 'true':
			call_args.append('--run-approximated-leiden')
		if '${approx_leiden_basis}' is not '':
			call_args.extend(['--approx-leiden-basis', '${approx_leiden_basis}'])
		if '${approx_leiden_resolution}' is not '':
			call_args.extend(['--approx-leiden-resolution', '${approx_leiden_resolution}'])
		if '${approx_leiden_class_label}' is not '':
			call_args.extend(['--approx-leiden-class-label', '${approx_leiden_class_label}'])
		if '${run_tsne}' is 'true':
			call_args.append('--run-tsne')
		if '${run_fitsne}' is 'true':
			call_args.append('--run-fitsne')
		if '${tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '${tsne_perplexity}'])
		if '${run_umap}' is 'true':
			call_args.append('--run-umap')
		if '${umap_K}' is not '':
			call_args.extend(['--umap-K', '${umap_K}'])
		if '${umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '${umap_min_dist}'])
		if '${umap_spread}' is not '':
			call_args.extend(['--umap-spread', '${umap_spread}'])
		if '${run_fle}' is 'true':
			call_args.append('--run-fle')
		if '${fle_K}' is not '':
			call_args.extend(['--fle-K', '${fle_K}'])
		if '${fle_target_change_per_node}' is not '':
			call_args.extend(['--fle-target-change-per-node', '${fle_target_change_per_node}'])
		if '${fle_target_steps}' is not '':
			call_args.extend(['--fle-target-steps', '${fle_target_steps}'])
		if '${fle_3D}' is 'true':
			call_args.append('--fle-3D')
		if '${net_down_sample_fraction}' is not '':
			call_args.extend(['--net-down-sample-fraction', '${net_down_sample_fraction}'])
		if '${net_ds_full_speed}' is 'true':
			call_args.append('--net-ds-full-speed')
		if '${run_net_tsne}' is 'true':
			call_args.append('--run-net-tsne')
		if '${net_tsne_out_basis}' is not '':
			call_args.extend(['--net-tsne-out-basis', '${net_tsne_out_basis}'])
		if '${run_net_fitsne}' is 'true':
			call_args.append('--run-net-fitsne')
		if '${net_fitsne_out_basis}' is not '':
			call_args.extend(['--net-fitsne-out-basis', '${net_fitsne_out_basis}'])
		if '${run_net_umap}' is 'true':
			call_args.append('--run-net-umap')
		if '${net_umap_out_basis}' is not '':
			call_args.extend(['--net-umap-out-basis', '${net_umap_out_basis}'])
		if '${run_net_fle}' is 'true':
			call_args.append('--run-net-fle')
		if '${net_fle_ds_full_speed}' is 'true':
			call_args.append('--net-fle-ds-full-speed')
		if '${net_fle_out_basis}' is not '':
			call_args.extend(['--net-fle-out-basis', '${net_fle_out_basis}'])
		print(' '.join(call_args))
		check_call(call_args)
		if '${output_parquet}' is 'true':
			call_args = ['scCloud', 'parquet', '${output_name}.h5ad', '${output_name}', '-p', '${num_cpu}']
			print(' '.join(call_args))
			check_call(call_args)
		CODE
	}

	output {
		File output_h5ad = "${output_name}.h5ad"
		Array[File] output_loom_file = glob("${output_name}.loom")
		Array[File] output_parquet_file = glob("${output_name}.parquet")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task organize_results {
	String output_name
	String sccloud_version
	String zones
	Int disk_space
	Int preemptible
	File? output_10x_h5
	File? output_h5ad
	Array[File]? output_seurat_h5ad
	Array[File]? output_filt_xlsx
	Array[File]? output_filt_plot
	Array[File]? output_loom_file
	Array[File]? output_parquet_file
	File? output_de_h5ad
	File? output_de_xlsx
	Array[File]? output_markers_xlsx
	Array[File]? output_anno_file
	Array[File]? output_pdfs
	Array[File]? output_htmls
	Array[File]? output_scp_files

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import os
		from subprocess import check_call
		dest = os.path.dirname('${output_name}') + '/'

		# check_call(['mkdir', '-p', dest])
		
		files = ['${output_10x_h5}', '${sep=" " output_seurat_h5ad}', '${sep=" " output_filt_xlsx}', '${sep=" " output_loom_file}', '${sep=" " output_parquet_file}', '${output_de_xlsx}', '${sep=" " output_markers_xlsx}', '${sep=" " output_anno_file}']
		files.append('${output_h5ad}' if '${output_de_h5ad}' is '' else '${output_de_h5ad}')
		files.extend('${sep="," output_filt_plot}'.split(','))
		files.extend('${sep="," output_pdfs}'.split(','))
		files.extend('${sep="," output_htmls}'.split(','))
		files.extend('${sep="," output_scp_files}'.split(','))
		for file in files:
			if file is not '':
				# call_args = ['cp', file, dest]
				call_args = ['gsutil', '-q', 'cp', file, dest]
				print(' '.join(call_args))
				check_call(call_args)
		CODE
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: "30 GB"
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task generate_hashing_cite_seq_tasks {
	File input_sample_sheet
	String sccloud_version
	String zones
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		import pandas as pd 
		from subprocess import check_call

		df = pd.read_csv('${input_sample_sheet}', header = 0, index_col = 0)
		with open('hashing.txt', 'w') as fo1, open('cite_seq.txt', 'w') as fo2, open('id2rna.txt', 'w') as fo3, open('id2adt.txt', 'w') as fo4:
			for outname, row in df.iterrows():
				if row['TYPE'] == 'cite-seq':
					fo2.write(outname + '\n')
				else:
					assert row['TYPE'] in ['cell-hashing', 'nuclei-hashing']
					fo1.write(outname + '\n')
				fo3.write(outname + '\t' + row['RNA'] + '\n')
				fo4.write(outname + '\t' + row['ADT'] + '\n')
		CODE
	}

	output {
		Array[String] hashing_ids = read_lines('hashing.txt')
		Array[String] cite_seq_ids = read_lines('cite_seq.txt')
		Map[String, String] id2rna = read_map('id2rna.txt')
		Map[String, String] id2adt = read_map('id2adt.txt')
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		preemptible: preemptible
	}
}

task run_scCloud_demuxEM {
	File input_adt_csv
	File input_raw_gene_bc_matrices_h5
	String output_dir
	String output_name
	String sccloud_version
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	String? genome
	Float? alpha_on_samples
	Int? min_num_genes
	Int? min_num_umis
	Float? min_signal_hashtag
	Int? random_state
	Boolean? generate_diagnostic_plots
	String? generate_gender_plot

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['scCloud', 'demuxEM', '${input_adt_csv}', '${input_raw_gene_bc_matrices_h5}', '${output_name}', '-p', '${num_cpu}']
		if '${genome}' is not '':
			call_args.extend(['--genome', '${genome}'])
		if '${alpha_on_samples}' is not '':
			call_args.extend(['--alpha_on_samples', '${alpha_on_samples}'])
		if '${min_num_genes}' is not '':
			call_args.extend(['--min-num-genes', '${min_num_genes}'])
		if '${min_num_umis}' is not '':
			call_args.extend(['--min-num-umis', '${min_num_umis}'])
		if '${min_signal_hashtag}' is not '':
			call_args.extend(['--min-signal-hashtag', '${min_signal_hashtag}'])
		if '${random_state}' is not '':
			call_args.extend(['--random-state', '${random_state}'])
		if '${generate_diagnostic_plots}' is 'true':
			call_args.append('--generate-diagnostic-plots')
		if '${generate_gender_plot}' is not '':
			call_args.extend(['--generate-gender-plot', '${generate_gender_plot}'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE

		gsutil -q cp ${output_name}_demux_10x.h5 ${output_dir}/${output_name}/
		gsutil -q cp ${output_name}_ADTs.h5ad ${output_dir}/${output_name}/
		gsutil -q cp ${output_name}_demux.h5ad ${output_dir}/${output_name}/
		gsutil -q -m cp ${output_name}.*.pdf ${output_dir}/${output_name}/
		# mkdir -p ${output_dir}/${output_name}
		# cp ${output_name}_demux_10x.h5 ${output_dir}/${output_name}/
		# cp ${output_name}_ADTs.h5ad ${output_dir}/${output_name}/
		# cp ${output_name}_demux.h5ad ${output_dir}/${output_name}/
		# cp ${output_name}.*.pdf ${output_dir}/${output_name}/
	}

	output {
		String output_folder = "${output_dir}/${output_name}"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_scCloud_merge_rna_adt {
	File input_raw_gene_bc_matrices_h5
	File input_adt_csv
	File? antibody_control_csv
	String output_dir
	String output_name
	String sccloud_version
	String zones
	String memory
	Int disk_space
	Int preemptible

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['scCloud', 'merge_rna_adt', '${input_raw_gene_bc_matrices_h5}', '${input_adt_csv}', '${output_name}_merged_10x.h5']
		if '${antibody_control_csv}' is not '':
			call_args.extend(['--antibody-control-csv', '${antibody_control_csv}'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE

		gsutil -q cp ${output_name}_merged_10x.h5 ${output_dir}/${output_name}/
		# mkdir -p ${output_dir}/${output_name}
		# cp ${output_name}_merged_10x.h5 ${output_dir}/${output_name}/
	}

	output {
		String output_folder = "${output_dir}/${output_name}"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "regevlab/sccloud-${sccloud_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}
