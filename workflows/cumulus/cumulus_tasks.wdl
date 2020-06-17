version 1.0

workflow cumulus_tasks {
}

task run_cumulus_aggregate_matrices {
	input {
		File input_count_matrix_csv
		String output_directory
		String output_name
		String cumulus_version
		String zones
		String memory
		Int disk_space
		Int preemptible
		String? restrictions
		String? attributes
		String? default_reference
		Boolean? select_only_singlets
		Int? minimum_number_of_genes
		String docker_registry
	}

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call

		call_args = ['pegasusio', 'aggregate_matrix', '~{input_count_matrix_csv}', '~{output_name}.aggr']
		if '~{restrictions}' is not '':
			ress = '~{restrictions}'.split(';')
			for res in ress:
				call_args.extend(['--restriction', res])
		if '~{attributes}' is not '':
			call_args.extend(['--attributes', '~{attributes}'])
		if '~{default_reference}' is not '':
			call_args.extend(['--default-reference', '~{default_reference}'])
		if '~{select_only_singlets}' is 'true':
			call_args.append('--select-only-singlets')
		if '~{minimum_number_of_genes}' is not '':
			call_args.extend(['--min-genes', '~{minimum_number_of_genes}'])

		print(' '.join(call_args))
		check_call(call_args)

		dest = '~{output_directory}' + '/' + '~{output_name}' + '/'
		# check_call(['mkdir', '-p', dest])
		# call_args = ['cp', '~{output_name}.aggr.zarr.zip', dest]
		call_args = ['gsutil', 'cp', '~{output_name}.aggr.zarr.zip', dest]
		print(' '.join(call_args))
		check_call(call_args)
		CODE
	}

	output {
		File output_zarr = '~{output_name}.aggr.zarr.zip'
	}

	runtime {
		docker: "~{docker_registry}/cumulus:~{cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_cumulus_cluster {
	input {
		File input_file
		String output_directory
		String output_name
		String cumulus_version
		String zones
		Int num_cpu
		String memory
		Int disk_space
		Int preemptible
		String? channel
		String? black_list
		Int? min_genes_before_filtration
		Boolean? select_singlets
		String? remap_singlets
		String? subset_singlets
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
		Boolean? correct_batch_effect
		String? correction_method
		String? batch_group_by
		Int? random_state
		File? gene_signature_file
		Int? nPC
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
		Float? net_down_sample_fraction
		Boolean? run_net_tsne
		String? net_tsne_out_basis
		Boolean? run_net_umap
		String? net_umap_out_basis
		Boolean? run_net_fle
		String? net_fle_out_basis
		String docker_registry
	}

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call

		call_args = ['pegasus', 'cluster', '~{input_file}', '~{output_name}', '-p', '~{num_cpu}']
		if '~{channel}' is not '':
			call_args.extend(['--channel', '~{channel}'])
		if '~{black_list}' is not '':
			call_args.extend(['--black-list', '~{black_list}'])
		if '~{min_genes_before_filtration}' is not '':
			call_args.extend(['--min-genes-before-filtration', '~{min_genes_before_filtration}'])
		if '~{select_singlets}' is 'true':
			call_args.append('--select-singlets')
		if '~{remap_singlets}' is not '':
			call_args.extend(['--remap-singlets', '~{remap_singlets}'])
		if '~{subset_singlets}' is not '':
			call_args.extend(['--subset-singlets', '~{subset_singlets}'])
		if '~{focus}' is not '':
			call_args.extend(['--focus', '~{focus}'])
		if '~{append}' is not '':
			call_args.extend(['--append', '~{append}'])
		if '~{output_filtration_results}' is 'true':
			call_args.append('--output-filtration-results')
		if '~{plot_filtration_results}' is 'true':
			call_args.append('--plot-filtration-results')
		if '~{plot_filtration_figsize}' is not '':
			call_args.extend(['--plot-filtration-figsize', '~{plot_filtration_figsize}'])
		if '~{output_h5ad}' is 'true':
			call_args.append('--output-h5ad')
		if '~{output_loom}' is 'true':
			call_args.append('--output-loom')
		if '~{correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
			if '~{correction_method}' is not '':
				call_args.extend(['--correction-method', '~{correction_method}'])
			if '~{batch_group_by}' is not '':
				call_args.extend(['--batch-group-by', '~{batch_group_by}'])
		if '~{min_genes}' is not '':
			call_args.extend(['--min-genes', '~{min_genes}'])
		if '~{max_genes}' is not '':
			call_args.extend(['--max-genes', '~{max_genes}'])
		if '~{min_umis}' is not '':
			call_args.extend(['--min-umis', '~{min_umis}'])
		if '~{max_umis}' is not '':
			call_args.extend(['--max-umis', '~{max_umis}'])
		if '~{mito_prefix}' is not '':
			call_args.extend(['--mito-prefix', '~{mito_prefix}'])
		if '~{percent_mito}' is not '' :
			call_args.extend(['--percent-mito', '~{percent_mito}'])
		if '~{gene_percent_cells}' is not '':
			call_args.extend(['--gene-percent-cells', '~{gene_percent_cells}'])
		if '~{counts_per_cell_after}' is not '':
			call_args.extend(['--counts-per-cell-after', '~{counts_per_cell_after}'])
		if '~{random_state}' is not '':
			call_args.extend(['--random-state', '~{random_state}'])
		if '~{gene_signature_file}' is not '':
			call_args.extend(['--calc-signature-scores', '~{gene_signature_file}'])
		if '~{no_select_hvf}' is 'true':
			call_args.append('--no-select-hvf')
		if '~{select_hvf_flavor}' is not '':
			call_args.extend(['--select-hvf-flavor', '~{select_hvf_flavor}'])
		if '~{select_hvf_ngenes}' is not '':
			call_args.extend(['--select-hvf-ngenes', '~{select_hvf_ngenes}'])
		if '~{nPC}' is not '':
			call_args.extend(['--pca-n', '~{nPC}'])
		if '~{knn_K}' is not '':
			call_args.extend(['--knn-K', '~{knn_K}'])
		if '~{knn_full_speed}' is 'true':
			call_args.append('--knn-full-speed')
		if '~{run_diffmap}' is 'true':
			call_args.extend(['--diffmap', '--diffmap-to-3d'])
		if '~{diffmap_ndc}' is not '':
			call_args.extend(['--diffmap-ndc', '~{diffmap_ndc}'])
		if '~{diffmap_maxt}' is not '':
			call_args.extend(['--diffmap-maxt', '~{diffmap_maxt}'])
		if '~{run_louvain}' is 'true':
			call_args.append('--louvain')
		if '~{louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '~{louvain_resolution}'])
		if '~{louvain_class_label}' is not '':
			call_args.extend(['--louvain-class-label', '~{louvain_class_label}'])
		if '~{run_leiden}' is 'true':
			call_args.append('--leiden')
		if '~{leiden_resolution}' is not '':
			call_args.extend(['--leiden-resolution', '~{leiden_resolution}'])
		if '~{leiden_niter}' is not '':
			call_args.extend(['--leiden-niter', '~{leiden_niter}'])
		if '~{leiden_class_label}' is not '':
			call_args.extend(['--leiden-class-label', '~{leiden_class_label}'])
		if '~{run_spectral_louvain}' is 'true':
			call_args.append('--spectral-louvain')
		if '~{spectral_louvain_basis}' is not '':
			call_args.extend(['--spectral-louvain-basis', '~{spectral_louvain_basis}'])
		if '~{spectral_louvain_resolution}' is not '':
			call_args.extend(['--spectral-louvain-resolution', '~{spectral_louvain_resolution}'])
		if '~{spectral_louvain_class_label}' is not '':
			call_args.extend(['--spectral-louvain-class-label', '~{spectral_louvain_class_label}'])
		if '~{run_spectral_leiden}' is 'true':
			call_args.append('--spectral-leiden')
		if '~{spectral_leiden_basis}' is not '':
			call_args.extend(['--spectral-leiden-basis', '~{spectral_leiden_basis}'])
		if '~{spectral_leiden_resolution}' is not '':
			call_args.extend(['--spectral-leiden-resolution', '~{spectral_leiden_resolution}'])
		if '~{spectral_leiden_class_label}' is not '':
			call_args.extend(['--spectral-leiden-class-label', '~{spectral_leiden_class_label}'])
		if '~{run_tsne}' is 'true':
			call_args.append('--tsne')
		if '~{run_fitsne}' is 'true':
			call_args.append('--fitsne')
		if '~{tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '~{tsne_perplexity}'])
		if '~{run_umap}' is 'true':
			call_args.append('--umap')
		if '~{umap_K}' is not '':
			call_args.extend(['--umap-K', '~{umap_K}'])
		if '~{umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '~{umap_min_dist}'])
		if '~{umap_spread}' is not '':
			call_args.extend(['--umap-spread', '~{umap_spread}'])
		if '~{run_fle}' is 'true':
			call_args.append('--fle')
		if '~{fle_K}' is not '':
			call_args.extend(['--fle-K', '~{fle_K}'])
		if '~{fle_target_change_per_node}' is not '':
			call_args.extend(['--fle-target-change-per-node', '~{fle_target_change_per_node}'])
		if '~{fle_target_steps}' is not '':
			call_args.extend(['--fle-target-steps', '~{fle_target_steps}'])
		if '~{net_down_sample_fraction}' is not '':
			call_args.extend(['--net-down-sample-fraction', '~{net_down_sample_fraction}'])
		if '~{run_net_tsne}' is 'true':
			call_args.append('--net-tsne')
		if '~{net_tsne_out_basis}' is not '':
			call_args.extend(['--net-tsne-out-basis', '~{net_tsne_out_basis}'])
		if '~{run_net_umap}' is 'true':
			call_args.append('--net-umap')
		if '~{net_umap_out_basis}' is not '':
			call_args.extend(['--net-umap-out-basis', '~{net_umap_out_basis}'])
		if '~{run_net_fle}' is 'true':
			call_args.append('--net-fle')
		if '~{net_fle_out_basis}' is not '':
			call_args.extend(['--net-fle-out-basis', '~{net_fle_out_basis}'])
		print(' '.join(call_args))
		check_call(call_args)

		import glob
		dest = '~{output_directory}' + '/' + '~{output_name}' + '/'
		# check_call(['mkdir', '-p', dest])
		files = ['~{output_name}.zarr.zip', '~{output_name}.log']
		if '~{output_h5ad}' is 'true':
			files.extend(glob.glob('~{output_name}.*.h5ad'))
		if '~{output_filtration_results}' is 'true':
			files.extend(glob.glob('~{output_name}.*.filt.xlsx'))
		if '~{plot_filtration_results}' is 'true':
			files.extend(glob.glob('~{output_name}.*.filt.*.pdf'))
		if '~{output_loom}' is 'true':
			files.extend(glob.glob('~{output_name}.*.loom'))
		for file in files:
			# call_args = ['cp', file, dest]
			call_args = ['gsutil', '-m', 'cp', file, dest]
			print(' '.join(call_args))
			check_call(call_args)
		CODE
	}

	output {
		File output_zarr = "~{output_name}.zarr.zip"
		Array[File] output_h5ad = glob("~{output_name}.*.h5ad")
		Array[File] output_filt_xlsx = glob("~{output_name}.*.filt.xlsx")
		Array[File] output_filt_plot = glob("~{output_name}.*.filt.*.pdf")
		Array[File] output_loom_file = glob("~{output_name}.*.loom")
		File output_log = "~{output_name}.log"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "~{docker_registry}/cumulus:~{cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_cumulus_cirro_output {
	input {
		File input_h5ad
		String output_directory
		String output_name
		String docker_registry
		String cumulus_version
		String zones
		String memory
		Int disk_space
		Int num_cpu
		Int preemptible
	}

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python /software/prepare_data.py --out ~{output_name}.cirro ~{input_h5ad}
		gsutil -q -m cp -r ~{output_name}.cirro ~{output_directory}/
		# mkdir -p ~{output_directory}/
		# cp -r ~{output_name}.cirro ~{output_directory}/
	}

	output {
		String output_cirro_folder = "~{output_directory}/~{output_name}.cirro"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "~{docker_registry}/cumulus:~{cumulus_version}"
		zones: zones
		memory: memory
		disks: "local-disk ~{disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_cumulus_de_analysis {
	input {
		File input_h5ad
		String output_directory
		String output_name
		String cumulus_version
		String zones
		Int num_cpu
		String memory
		Int disk_space
		Int preemptible
		String? labels
		Boolean? auc
		Boolean? t_test
		Boolean? fisher
		Boolean? mwu
		Float? alpha

		Boolean? annotate_cluster
		String? annotate_de_test
		String? organism
		Float? minimum_report_score

		Boolean? find_markers_lightgbm
		Boolean? remove_ribo
		Float? min_gain
		Int? random_state
		String docker_registry
	}

	Boolean is_url = defined(organism) && sub(select_first([organism]), "^.+\\.json", "URL") == "URL"

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		if [ ~{is_url} == true ]; then
			gsutil -q cp ~{organism} markers.json
		fi

		python <<CODE
		from subprocess import check_call

		call_args = ['mv', '-f', '~{input_h5ad}', '~{output_name}.h5ad']
		print(' '.join(call_args))
		check_call(call_args)
		call_args = ['pegasus', 'de_analysis', '~{output_name}.h5ad', '~{output_name}.de.xlsx', '-p', '~{num_cpu}']
		if '~{labels}' is not '':
			call_args.extend(['--labels', '~{labels}'])
		if '~{alpha}' is not '':
			call_args.extend(['--alpha', '~{alpha}'])
		if '~{auc}' is 'true':
			call_args.append('--auc')
		if '~{t_test}' is 'true':
			call_args.append('--t')
		if '~{fisher}' is 'true':
			call_args.append('--fisher')
		if '~{mwu}' is 'true':
			call_args.append('--mwu')
		print(' '.join(call_args))
		check_call(call_args)
		if '~{find_markers_lightgbm}' is 'true':
			call_args = ['pegasus', 'find_markers', '~{output_name}.h5ad', '~{output_name}.markers.xlsx', '-p', '~{num_cpu}']
			if '~{labels}' is not '':
				call_args.extend(['--labels', '~{labels}'])
			if '~{remove_ribo}' is 'true':
				call_args.append('--remove-ribo')
			if '~{min_gain}' is not '':
				call_args.extend(['--min-gain', '~{min_gain}'])
			if '~{random_state}' is not '':
				call_args.extend(['--random-state', '~{random_state}'])
			print(' '.join(call_args))
			check_call(call_args)
		if '~{annotate_cluster}' is 'true':
			call_args = ['pegasus', 'annotate_cluster', '~{output_name}.h5ad', '~{output_name}.anno.txt']
			if '~{organism}' is not '':
				if '~{is_url}' is 'true':
					call_args.extend(['--marker-file', 'markers.json'])
				else:
					call_args.extend(['--marker-file', '~{organism}'])
			if '~{annotate_de_test}' is not '':
				call_args.extend(['--de-test', '~{annotate_de_test}'])
			if '~{alpha}' is not '':
				call_args.extend(['--de-alpha', '~{alpha}'])
			if '~{minimum_report_score}' is not '':
				call_args.extend(['--minimum-report-score', '~{minimum_report_score}'])
			print(' '.join(call_args))
			check_call(call_args)

		dest = '~{output_directory}' + '/'
		# check_call(['mkdir', '-p', dest])
		files = ['~{output_name}.h5ad', '~{output_name}.de.xlsx']
		if '~{find_markers_lightgbm}' is 'true':
			files.append('~{output_name}.markers.xlsx')
		if '~{annotate_cluster}' is 'true':
			files.append('~{output_name}.anno.txt')
		for file in files:
			# call_args = ['cp', file, dest]
			call_args = ['gsutil', 'cp', file, dest]
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
		docker: "~{docker_registry}/cumulus:~{cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_cumulus_plot {
	input {
		File input_h5ad
		String output_directory
		String output_name
		String cumulus_version
		String zones
		String memory
		Int disk_space
		Int preemptible
		String? plot_composition
		String? plot_tsne
		String? plot_fitsne
		String? plot_umap
		String? plot_fle
		String? plot_net_tsne
		String? plot_net_umap
		String? plot_net_fle
		String docker_registry
	}

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call
		if '~{plot_composition}' is not '':
			pairs = '~{plot_composition}'.split(',')
			for pair in pairs:
				lab, attr = pair.split(':')
				call_args = ['pegasus', 'plot', 'composition', '--cluster-labels', lab, '--attribute', attr, '--style', 'normalized', '--not-stacked', '~{input_h5ad}', '~{output_name}.' + lab + '.' + attr + '.composition.pdf']
				print(' '.join(call_args))
				check_call(call_args)
		if '~{plot_tsne}' is not '':
			call_args = ['pegasus', 'plot', 'scatter',  '--basis', 'tsne', '--attributes', '~{plot_tsne}', '~{input_h5ad}', '~{output_name}.tsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '~{plot_fitsne}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'fitsne', '--attributes', '~{plot_fitsne}', '~{input_h5ad}', '~{output_name}.fitsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '~{plot_umap}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'umap', '--attributes', '~{plot_umap}', '~{input_h5ad}', '~{output_name}.umap.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '~{plot_fle}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'fle', '--attributes', '~{plot_fle}', '~{input_h5ad}', '~{output_name}.fle.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '~{plot_net_tsne}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_tsne', '--attributes', '~{plot_net_tsne}', '~{input_h5ad}', '~{output_name}.net.tsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '~{plot_net_umap}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_umap', '--attributes', '~{plot_net_umap}', '~{input_h5ad}', '~{output_name}.net.umap.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '~{plot_net_fle}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_fle', '--attributes', '~{plot_net_fle}', '~{input_h5ad}', '~{output_name}.net.fle.pdf']
			print(' '.join(call_args))
			check_call(call_args)

		import glob
		dest = '~{output_directory}' + '/'
		# check_call(['mkdir', '-p', dest])
		files = glob.glob('*.pdf')
		files.extend(glob.glob('*.html'))

		for file in files:
			# call_args = ['cp', file, dest]
			call_args = ['gsutil', '-m', 'cp', file, dest]
			print(' '.join(call_args))
			check_call(call_args)
		CODE
	}

	output {
		Array[File] output_pdfs = glob("*.pdf")
		Array[File] output_htmls = glob("*.html")
	}

	runtime {
		docker: "~{docker_registry}/cumulus:~{cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_cumulus_scp_output {
	input {
		File input_h5ad
		String output_directory
		String output_name
		Boolean output_dense
		String cumulus_version
		String zones
		String memory
		Int disk_space
		Int preemptible
		String docker_registry
	}

	command {
		set -e
		export TMPDIR=/tmp
		pegasus scp_output ~{true='--dense' false='' output_dense} ~{input_h5ad} ~{output_name}
		# mkdir -p ~{output_directory} ; cp ~{output_name}.scp.* ~{output_directory}/
		gsutil -m cp ~{output_name}.scp.* ~{output_directory}/
	}

	output {
		Array[File] output_scp_files = glob("~{output_name}.scp.*")
	}

	runtime {
		docker: "~{docker_registry}/cumulus:~{cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_cumulus_subcluster {
	input {
		File input_h5ad
		String output_directory
		String output_name
		String cumulus_version
		String zones
		Int num_cpu
		String memory
		Int disk_space
		Int preemptible
		String subset_selections
		Boolean? correct_batch_effect
		String? correction_method
		String? batch_group_by
		Boolean? output_loom
		String? select_hvf_flavor
		Int? select_hvf_ngenes
		Boolean? no_select_hvf
		Int? random_state
		Int? nPC
		Int? knn_K
		Boolean? knn_full_speed
		Boolean? run_diffmap
		Int? diffmap_ndc
		Int? diffmap_maxt
		String? calculate_pseudotime
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
		Float? net_down_sample_fraction
		Boolean? run_net_tsne
		String? net_tsne_out_basis
		Boolean? run_net_umap
		String? net_umap_out_basis
		Boolean? run_net_fle
		String? net_fle_out_basis
		String docker_registry
	}

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['pegasus', 'subcluster', '~{input_h5ad}', '~{output_name}', '-p', '~{num_cpu}']
		if '~{subset_selections}' is not '':
			sels = '~{subset_selections}'.split(';')
			for sel in sels:
				call_args.extend(['--subset-selection', sel])
		if '~{correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
			if '~{correction_method}' is not '':
				call_args.extend(['--correction-method', '~{correction_method}'])
			if '~{batch_group_by}' is not '':
				call_args.extend(['--batch-group-by', '~{batch_group_by}'])
		if '~{output_loom}' is 'true':
			call_args.append('--output-loom')
		if '~{select_hvf_flavor}' is not '':
			call_args.extend(['--select-hvf-flavor', '~{select_hvf_flavor}'])
		if '~{select_hvf_ngenes}' is not '':
			call_args.extend(['--select-hvf-ngenes', '~{select_hvf_ngenes}'])
		if '~{no_select_hvf}' is 'true':
			call_args.append('--no-select-hvf')
		if '~{random_state}' is not '':
			call_args.extend(['--random-state', '~{random_state}'])
		if '~{nPC}' is not '':
			call_args.extend(['--pca-n', '~{nPC}'])
		if '~{knn_K}' is not '':
			call_args.extend(['--knn-K', '~{knn_K}'])
		if '~{knn_full_speed}' is 'true':
			call_args.append('--knn-full-speed')
		if '~{run_diffmap}' is 'true':
			call_args.append('--diffmap')
		if '~{diffmap_ndc}' is not '':
			call_args.extend(['--diffmap-ndc', '~{diffmap_ndc}'])
		if '~{diffmap_maxt}' is not '':
			call_args.extend(['--diffmap-maxt', '~{diffmap_maxt}'])
		if '~{calculate_pseudotime}' is not '':
			call_args.extend(['--calculate-pseudotime', '~{calculate_pseudotime}'])
		if '~{run_louvain}' is 'true':
			call_args.append('--louvain')
		if '~{louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '~{louvain_resolution}'])
		if '~{louvain_class_label}' is not '':
			call_args.extend(['--louvain-class-label', '~{louvain_class_label}'])
		if '~{run_leiden}' is 'true':
			call_args.append('--leiden')
		if '~{leiden_resolution}' is not '':
			call_args.extend(['--leiden-resolution', '~{leiden_resolution}'])
		if '~{leiden_niter}' is not '':
			call_args.extend(['--leiden-niter', '~{leiden_niter}'])
		if '~{leiden_class_label}' is not '':
			call_args.extend(['--leiden-class-label', '~{leiden_class_label}'])
		if '~{run_spectral_louvain}' is 'true':
			call_args.append('--spectral-louvain')
		if '~{spectral_louvain_basis}' is not '':
			call_args.extend(['--spectral-louvain-basis', '~{spectral_louvain_basis}'])
		if '~{spectral_louvain_resolution}' is not '':
			call_args.extend(['--spectral-louvain-resolution', '~{spectral_louvain_resolution}'])
		if '~{spectral_louvain_class_label}' is not '':
			call_args.extend(['--spectral-louvain-class-label', '~{spectral_louvain_class_label}'])
		if '~{run_spectral_leiden}' is 'true':
			call_args.append('--spectral-leiden')
		if '~{spectral_leiden_basis}' is not '':
			call_args.extend(['--spectral-leiden-basis', '~{spectral_leiden_basis}'])
		if '~{spectral_leiden_resolution}' is not '':
			call_args.extend(['--spectral-leiden-resolution', '~{spectral_leiden_resolution}'])
		if '~{spectral_leiden_class_label}' is not '':
			call_args.extend(['--spectral-leiden-class-label', '~{spectral_leiden_class_label}'])
		if '~{run_tsne}' is 'true':
			call_args.append('--tsne')
		if '~{run_fitsne}' is 'true':
			call_args.append('--fitsne')
		if '~{tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '~{tsne_perplexity}'])
		if '~{run_umap}' is 'true':
			call_args.append('--umap')
		if '~{umap_K}' is not '':
			call_args.extend(['--umap-K', '~{umap_K}'])
		if '~{umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '~{umap_min_dist}'])
		if '~{umap_spread}' is not '':
			call_args.extend(['--umap-spread', '~{umap_spread}'])
		if '~{run_fle}' is 'true':
			call_args.append('--fle')
		if '~{fle_K}' is not '':
			call_args.extend(['--fle-K', '~{fle_K}'])
		if '~{fle_target_change_per_node}' is not '':
			call_args.extend(['--fle-target-change-per-node', '~{fle_target_change_per_node}'])
		if '~{fle_target_steps}' is not '':
			call_args.extend(['--fle-target-steps', '~{fle_target_steps}'])
		if '~{net_down_sample_fraction}' is not '':
			call_args.extend(['--net-down-sample-fraction', '~{net_down_sample_fraction}'])
		if '~{run_net_tsne}' is 'true':
			call_args.append('--net-tsne')
		if '~{net_tsne_out_basis}' is not '':
			call_args.extend(['--net-tsne-out-basis', '~{net_tsne_out_basis}'])
		if '~{run_net_umap}' is 'true':
			call_args.append('--net-umap')
		if '~{net_umap_out_basis}' is not '':
			call_args.extend(['--net-umap-out-basis', '~{net_umap_out_basis}'])
		if '~{run_net_fle}' is 'true':
			call_args.append('--net-fle')
		if '~{net_fle_out_basis}' is not '':
			call_args.extend(['--net-fle-out-basis', '~{net_fle_out_basis}'])
		print(' '.join(call_args))
		check_call(call_args)

		dest = '~{output_directory}' + '/' + '~{output_name}' + '/'
		# check_call(['mkdir', '-p', dest])
		files = ['~{output_name}.h5ad', '~{output_name}.log']
		if '~{output_loom}' is 'true':
			files.append('~{output_name}.loom')
		for file in files:
			# call_args = ['cp', file, dest]
			call_args = ['gsutil', 'cp', file, dest]
			print(' '.join(call_args))
			check_call(call_args)
		CODE
	}

	output {
		File output_h5ad = "~{output_name}.h5ad"
		File output_log = "~{output_name}.log"
		Array[File] output_loom_file = glob("~{output_name}.loom")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "~{docker_registry}/cumulus:~{cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ~{disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}
