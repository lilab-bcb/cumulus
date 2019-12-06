workflow cumulus_tasks {
}

task run_cumulus_aggregate_matrices {
	File input_count_matrix_csv
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

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call
		call_args = ['pegasus', 'aggregate_matrix', '${input_count_matrix_csv}', '${output_name}']
		if '${restrictions}' is not '':
			ress = '${restrictions}'.split(';')
			for res in ress:
				call_args.extend(['--restriction', res])
		if '${attributes}' is not '':
			call_args.extend(['--attributes', '${attributes}'])
		if '${default_reference}' is not '':
			call_args.extend(['--default-reference', '${default_reference}'])
		if '${select_only_singlets}' is 'true':
			call_args.append('--select-only-singlets')
		if '${minimum_number_of_genes}' is not '':
			call_args.extend(['--minimum-number-of-genes', '${minimum_number_of_genes}'])

		print(' '.join(call_args))
		check_call(call_args)

		import os
		dest = os.path.dirname('${output_name}') + '/'
		# check_call(['mkdir', '-p', dest])
		# call_args = ['cp', '${output_name}.h5sc', dest]
		call_args = ['gsutil', 'cp', '${output_name}.h5sc', dest]
		print(' '.join(call_args))
		check_call(call_args)
		CODE
	}

	output {
		File output_h5sc = '${output_name}.h5sc'
	}

	runtime {
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_cumulus_cluster {
	File input_file
	String output_name
	String cumulus_version
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	String? considered_refs
	String? channel
	String? black_list
	Int? min_genes_on_raw
	Boolean? select_singlets
	Boolean? cite_seq
	Float? cite_seq_capping
	Boolean? output_filtration_results
	Boolean? plot_filtration_results
	String? plot_filtration_figsize
	Boolean? output_seurat_compatible
	Boolean? output_loom
	Boolean? output_parquet
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
	String? batch_group_by
	Int? random_state
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

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['pegasus', 'cluster', '${input_file}', '${output_name}', '-p', '${num_cpu}']
		if '${considered_refs}' is not '':
			call_args.extend(['--considered-refs', '${considered_refs}'])
		if '${channel}' is not '':
			call_args.extend(['--channel', '${channel}'])
		if '${black_list}' is not '':
			call_args.extend(['--black-list', '${black_list}'])
		if '${min_genes_on_raw}' is not '':
			call_args.extend(['--min-genes-on-raw', '${min_genes_on_raw}'])
		if '${select_singlets}' is 'true':
			call_args.append('--select-singlets')
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
		if '${counts_per_cell_after}' is not '':
			call_args.extend(['--counts-per-cell-after', '${counts_per_cell_after}'])
		if '${random_state}' is not '':
			call_args.extend(['--random-state', '${random_state}'])
		if '${no_select_hvf}' is 'true':
			call_args.append('--no-select-hvf')
		if '${select_hvf_flavor}' is not '':
			call_args.extend(['--select-hvf-flavor', '${select_hvf_flavor}'])
		if '${select_hvf_ngenes}' is not '':
			call_args.extend(['--select-hvf-ngenes', '${select_hvf_ngenes}'])
		if '${nPC}' is not '':
			call_args.extend(['--nPC', '${nPC}'])
		if '${knn_K}' is not '':
			call_args.extend(['--knn-K', '${knn_K}'])
		if '${knn_full_speed}' is 'true':
			call_args.append('--knn-full-speed')
		if '${run_diffmap}' is 'true':
			call_args.append('--diffmap')
		if '${diffmap_ndc}' is not '':
			call_args.extend(['--diffmap-ndc', '${diffmap_ndc}'])
		if '${diffmap_maxt}' is not '':
			call_args.extend(['--diffmap-maxt', '${diffmap_maxt}'])
		if '${run_louvain}' is 'true':
			call_args.append('--louvain')
		if '${louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '${louvain_resolution}'])
		if '${louvain_class_label}' is not '':
			call_args.extend(['--louvain-class-label', '${louvain_class_label}'])
		if '${run_leiden}' is 'true':
			call_args.append('--leiden')
		if '${leiden_resolution}' is not '':
			call_args.extend(['--leiden-resolution', '${leiden_resolution}'])
		if '${leiden_niter}' is not '':
			call_args.extend(['--leiden-niter', '${leiden_niter}'])
		if '${leiden_class_label}' is not '':
			call_args.extend(['--leiden-class-label', '${leiden_class_label}'])
		if '${run_spectral_louvain}' is 'true':
			call_args.append('--spectral-louvain')
		if '${spectral_louvain_basis}' is not '':
			call_args.extend(['--spectral-louvain-basis', '${spectral_louvain_basis}'])
		if '${spectral_louvain_resolution}' is not '':
			call_args.extend(['--spectral-louvain-resolution', '${spectral_louvain_resolution}'])
		if '${spectral_louvain_class_label}' is not '':
			call_args.extend(['--spectral-louvain-class-label', '${spectral_louvain_class_label}'])
		if '${run_spectral_leiden}' is 'true':
			call_args.append('--spectral-leiden')
		if '${spectral_leiden_basis}' is not '':
			call_args.extend(['--spectral-leiden-basis', '${spectral_leiden_basis}'])
		if '${spectral_leiden_resolution}' is not '':
			call_args.extend(['--spectral-leiden-resolution', '${spectral_leiden_resolution}'])
		if '${spectral_leiden_class_label}' is not '':
			call_args.extend(['--spectral-leiden-class-label', '${spectral_leiden_class_label}'])
		if '${run_tsne}' is 'true':
			call_args.append('--tsne')
		if '${run_fitsne}' is 'true':
			call_args.append('--fitsne')
		if '${tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '${tsne_perplexity}'])
		if '${run_umap}' is 'true':
			call_args.append('--umap')
		if '${umap_K}' is not '':
			call_args.extend(['--umap-K', '${umap_K}'])
		if '${umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '${umap_min_dist}'])
		if '${umap_spread}' is not '':
			call_args.extend(['--umap-spread', '${umap_spread}'])
		if '${run_fle}' is 'true':
			call_args.append('--fle')
		if '${fle_K}' is not '':
			call_args.extend(['--fle-K', '${fle_K}'])
		if '${fle_target_change_per_node}' is not '':
			call_args.extend(['--fle-target-change-per-node', '${fle_target_change_per_node}'])
		if '${fle_target_steps}' is not '':
			call_args.extend(['--fle-target-steps', '${fle_target_steps}'])
		if '${net_down_sample_fraction}' is not '':
			call_args.extend(['--net-down-sample-fraction', '${net_down_sample_fraction}'])
		if '${run_net_tsne}' is 'true':
			call_args.append('--net-tsne')
		if '${net_tsne_out_basis}' is not '':
			call_args.extend(['--net-tsne-out-basis', '${net_tsne_out_basis}'])
		if '${run_net_umap}' is 'true':
			call_args.append('--net-umap')
		if '${net_umap_out_basis}' is not '':
			call_args.extend(['--net-umap-out-basis', '${net_umap_out_basis}'])
		if '${run_net_fle}' is 'true':
			call_args.append('--net-fle')
		if '${net_fle_out_basis}' is not '':
			call_args.extend(['--net-fle-out-basis', '${net_fle_out_basis}'])
		print(' '.join(call_args))
		check_call(call_args)
		if '${output_parquet}' is 'true':
			call_args = ['pegasus', 'parquet', '${output_name}.h5ad', '${output_name}', '-p', '${num_cpu}']
			print(' '.join(call_args))
			check_call(call_args)

		import os
		dest = os.path.dirname('${output_name}') + '/'
		# check_call(['mkdir', '-p', dest])
		files = ['${output_name}.h5ad', '${output_name}.log']
		if '${output_seurat_compatible}' is 'true':
			files.append('${output_name}.seurat.h5ad')
		if '${output_filtration_results}' is 'true':
			files.append('${output_name}.filt.xlsx')
		if '${plot_filtration_results}' is 'true':
			files.append('${output_name}.filt.*.pdf')
		if '${output_loom}' is 'true':
			files.append('${output_name}.loom')
		if '${output_parquet}' is 'true':
			files.append('${output_name}.parquet')
		for file in files:
			# call_args = ['cp', file, dest]
			call_args = ['gsutil', '-m', 'cp', file, dest]
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
		File output_log = "${output_name}.log"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_cumulus_de_analysis {
	File input_h5ad
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

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['mv', '-f', '${input_h5ad}', '${output_name}.h5ad']
		print(' '.join(call_args))
		check_call(call_args)			
		call_args = ['pegasus', 'de_analysis', '${output_name}.h5ad', '${output_name}.de.xlsx', '-p', '${num_cpu}']
		if '${labels}' is not '':
			call_args.extend(['--labels', '${labels}'])
		if '${alpha}' is not '':
			call_args.extend(['--alpha', '${alpha}'])
		if '${auc}' is 'true':
			call_args.append('--auc')
		if '${t_test}' is 'true':
			call_args.append('--t')
		if '${fisher}' is 'true':
			call_args.append('--fisher')
		if '${mwu}' is 'true':
			call_args.append('--mwu')
		print(' '.join(call_args))
		check_call(call_args)
		if '${find_markers_lightgbm}' is 'true':
			call_args = ['pegasus', 'find_markers', '${output_name}.h5ad', '${output_name}.markers.xlsx', '-p', '${num_cpu}']
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
			call_args = ['pegasus', 'annotate_cluster', '${output_name}.h5ad', '${output_name}' + '.anno.txt']
			if '${organism}' is not '':
				call_args.extend(['--marker-file', '${organism}'])
			if '${annotate_de_test}' is not '':
				call_args.extend(['--de-test', '${annotate_de_test}'])
			if '${alpha}' is not '':
				call_args.extend(['--de-alpha', '${alpha}'])
			if '${minimum_report_score}' is not '':
				call_args.extend(['--minimum-report-score', '${minimum_report_score}'])
			print(' '.join(call_args))
			check_call(call_args)

		import os
		dest = os.path.dirname('${output_name}') + '/'
		# check_call(['mkdir', '-p', dest])
		files = ['${output_name}.h5ad', '${output_name}.de.xlsx']
		if '${find_markers_lightgbm}' is 'true':
			files.append('${output_name}.markers.xlsx')
		if '${annotate_cluster}' is 'true':
			files.append('${output_name}.anno.txt')
		for file in files:
			# call_args = ['cp', file, dest]
			call_args = ['gsutil', 'cp', file, dest]
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
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_cumulus_plot {
	File input_h5ad
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
	String? plot_diffmap
	String? plot_citeseq_fitsne
	String? plot_net_tsne
	String? plot_net_umap
	String? plot_net_fle
    String docker_registry

	command {
		set -e
		export TMPDIR=/tmp

		python <<CODE
		from subprocess import check_call
		if '${plot_composition}' is not '':
			pairs = '${plot_composition}'.split(',')
			for pair in pairs:
				lab, attr = pair.split(':')
				call_args = ['pegasus', 'plot', 'composition', '--cluster-labels', lab, '--attribute', attr, '--style', 'normalized', '--not-stacked', '${input_h5ad}', '${output_name}.' + lab + '.' + attr + '.composition.pdf']
				print(' '.join(call_args))
				check_call(call_args)
		if '${plot_tsne}' is not '':
			call_args = ['pegasus', 'plot', 'scatter',  '--basis', 'tsne', '--attributes', '${plot_tsne}', '${input_h5ad}', '${output_name}.tsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_fitsne}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'fitsne', '--attributes', '${plot_fitsne}', '${input_h5ad}', '${output_name}.fitsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_umap}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'umap', '--attributes', '${plot_umap}', '${input_h5ad}', '${output_name}.umap.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_fle}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'fle', '--attributes', '${plot_fle}', '${input_h5ad}', '${output_name}.fle.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_diffmap}' is not '':
			attrs = '${plot_diffmap}'.split(',')
			for attr in attrs:
				call_args = ['pegasus', 'iplot', '--attribute', attr, 'diffmap_pca', '${input_h5ad}', '${output_name}.' + attr + '.diffmap_pca.html']
				print(' '.join(call_args))
				check_call(call_args)
		if '${plot_citeseq_fitsne}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'citeseq_fitsne', '--attributes', '${plot_citeseq_fitsne}', '${input_h5ad}', '${output_name}.citeseq.fitsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_net_tsne}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_tsne', '--attributes', '${plot_net_tsne}', '${input_h5ad}', '${output_name}.net.tsne.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_net_umap}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_umap', '--attributes', '${plot_net_umap}', '${input_h5ad}', '${output_name}.net.umap.pdf']
			print(' '.join(call_args))
			check_call(call_args)
		if '${plot_net_fle}' is not '':
			call_args = ['pegasus', 'plot', 'scatter', '--basis', 'net_fle', '--attributes', '${plot_net_fle}', '${input_h5ad}', '${output_name}.net.fle.pdf']
			print(' '.join(call_args))
			check_call(call_args)

		import os
		import glob
		dest = os.path.dirname('${output_name}') + '/'
		# check_call(['mkdir', '-p', dest])
		files = []
		if len(glob.glob('*.pdf')) > 0:
			files.append('*.pdf')
		if len(glob.glob('*.html')) > 0:
			files.append('*.html')
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
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_cumulus_scp_output {
	File input_h5ad
	String output_name
	Boolean output_dense
	String cumulus_version
	String zones
	String memory
	Int disk_space
	Int preemptible
	String docker_registry

	command {
		set -e
		export TMPDIR=/tmp
		export DIRNAME=`dirname ${output_name}`
		pegasus scp_output ${true='--dense' false='' output_dense} ${input_h5ad} ${output_name}
		# mkdir -p ${DIRNAME} ; cp ${output_name}.scp.* ${DIRNAME}
		gsutil -m cp ${output_name}.scp.* ${DIRNAME}
	}

	output {
		Array[File] output_scp_files = glob("${output_name}.scp.*")
	}

	runtime {
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}

task run_cumulus_subcluster {
	File input_h5ad
	String output_name
	String cumulus_version
	String zones
	Int num_cpu
	String memory
	Int disk_space
	Int preemptible
	String subset_selections 
	Boolean? correct_batch_effect
	String? batch_group_by
	Boolean? output_loom
	Boolean? output_parquet
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

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['pegasus', 'subcluster', '${input_h5ad}', '${output_name}', '-p', '${num_cpu}']
		if '${subset_selections}' is not '':
			sels = '${subset_selections}'.split(';')
			for sel in sels:
				call_args.extend(['--subset-selection', sel])
		if '${correct_batch_effect}' is 'true':
			call_args.append('--correct-batch-effect')
			if '${batch_group_by}' is not '':
				call_args.extend(['--batch-group-by', '${batch_group_by}'])
		if '${output_loom}' is 'true':
			call_args.append('--output-loom')
		if '${select_hvf_flavor}' is not '':
			call_args.extend(['--select-hvf-flavor', '${select_hvf_flavor}'])
		if '${select_hvf_ngenes}' is not '':
			call_args.extend(['--select-hvf-ngenes', '${select_hvf_ngenes}'])
		if '${no_select_hvf}' is 'true':
			call_args.append('--no-select-hvf')
		if '${random_state}' is not '':
			call_args.extend(['--random-state', '${random_state}'])
		if '${nPC}' is not '':
			call_args.extend(['--nPC', '${nPC}'])
		if '${knn_K}' is not '':
			call_args.extend(['--knn-K', '${knn_K}'])
		if '${knn_full_speed}' is 'true':
			call_args.append('--knn-full-speed')
		if '${run_diffmap}' is 'true':
			call_args.append('--diffmap')
		if '${diffmap_ndc}' is not '':
			call_args.extend(['--diffmap-ndc', '${diffmap_ndc}'])
		if '${diffmap_maxt}' is not '':
			call_args.extend(['--diffmap-maxt', '${diffmap_maxt}'])
		if '${calculate_pseudotime}' is not '':
			call_args.extend(['--calculate-pseudotime', '${calculate_pseudotime}'])
		if '${run_louvain}' is 'true':
			call_args.append('--louvain')
		if '${louvain_resolution}' is not '':
			call_args.extend(['--louvain-resolution', '${louvain_resolution}'])
		if '${louvain_class_label}' is not '':
			call_args.extend(['--louvain-class-label', '${louvain_class_label}'])
		if '${run_leiden}' is 'true':
			call_args.append('--leiden')
		if '${leiden_resolution}' is not '':
			call_args.extend(['--leiden-resolution', '${leiden_resolution}'])
		if '${leiden_niter}' is not '':
			call_args.extend(['--leiden-niter', '${leiden_niter}'])
		if '${leiden_class_label}' is not '':
			call_args.extend(['--leiden-class-label', '${leiden_class_label}'])
		if '${run_spectral_louvain}' is 'true':
			call_args.append('--spectral-louvain')
		if '${spectral_louvain_basis}' is not '':
			call_args.extend(['--spectral-louvain-basis', '${spectral_louvain_basis}'])
		if '${spectral_louvain_resolution}' is not '':
			call_args.extend(['--spectral-louvain-resolution', '${spectral_louvain_resolution}'])
		if '${spectral_louvain_class_label}' is not '':
			call_args.extend(['--spectral-louvain-class-label', '${spectral_louvain_class_label}'])
		if '${run_spectral_leiden}' is 'true':
			call_args.append('--spectral-leiden')
		if '${spectral_leiden_basis}' is not '':
			call_args.extend(['--spectral-leiden-basis', '${spectral_leiden_basis}'])
		if '${spectral_leiden_resolution}' is not '':
			call_args.extend(['--spectral-leiden-resolution', '${spectral_leiden_resolution}'])
		if '${spectral_leiden_class_label}' is not '':
			call_args.extend(['--spectral-leiden-class-label', '${spectral_leiden_class_label}'])
		if '${run_tsne}' is 'true':
			call_args.append('--tsne')
		if '${run_fitsne}' is 'true':
			call_args.append('--fitsne')
		if '${tsne_perplexity}' is not '':
			call_args.extend(['--tsne-perplexity', '${tsne_perplexity}'])
		if '${run_umap}' is 'true':
			call_args.append('--umap')
		if '${umap_K}' is not '':
			call_args.extend(['--umap-K', '${umap_K}'])
		if '${umap_min_dist}' is not '':
			call_args.extend(['--umap-min-dist', '${umap_min_dist}'])
		if '${umap_spread}' is not '':
			call_args.extend(['--umap-spread', '${umap_spread}'])
		if '${run_fle}' is 'true':
			call_args.append('--fle')
		if '${fle_K}' is not '':
			call_args.extend(['--fle-K', '${fle_K}'])
		if '${fle_target_change_per_node}' is not '':
			call_args.extend(['--fle-target-change-per-node', '${fle_target_change_per_node}'])
		if '${fle_target_steps}' is not '':
			call_args.extend(['--fle-target-steps', '${fle_target_steps}'])
		if '${net_down_sample_fraction}' is not '':
			call_args.extend(['--net-down-sample-fraction', '${net_down_sample_fraction}'])
		if '${run_net_tsne}' is 'true':
			call_args.append('--net-tsne')
		if '${net_tsne_out_basis}' is not '':
			call_args.extend(['--net-tsne-out-basis', '${net_tsne_out_basis}'])
		if '${run_net_umap}' is 'true':
			call_args.append('--net-umap')
		if '${net_umap_out_basis}' is not '':
			call_args.extend(['--net-umap-out-basis', '${net_umap_out_basis}'])
		if '${run_net_fle}' is 'true':
			call_args.append('--net-fle')
		if '${net_fle_out_basis}' is not '':
			call_args.extend(['--net-fle-out-basis', '${net_fle_out_basis}'])
		print(' '.join(call_args))
		check_call(call_args)
		if '${output_parquet}' is 'true':
			call_args = ['pegasus', 'parquet', '${output_name}.h5ad', '${output_name}', '-p', '${num_cpu}']
			print(' '.join(call_args))
			check_call(call_args)

		import os
		dest = os.path.dirname('${output_name}') + '/'
		# check_call(['mkdir', '-p', dest])
		files = ['${output_name}.h5ad', '${output_name}.log']
		if '${output_loom}' is 'true':
			files.append('${output_name}.loom')
		if '${output_parquet}' is 'true':
			files.append('${output_name}.parquet')
		for file in files:
			# call_args = ['cp', file, dest]
			call_args = ['gsutil', 'cp', file, dest]
			print(' '.join(call_args))
			check_call(call_args)
		CODE
	}

	output {
		File output_h5ad = "${output_name}.h5ad"
		File output_log = "${output_name}.log"
		Array[File] output_loom_file = glob("${output_name}.loom")
		Array[File] output_parquet_file = glob("${output_name}.parquet")
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task generate_hashing_cite_seq_tasks {
	File input_sample_sheet
	String cumulus_version
	String zones
	Int preemptible
	String docker_registry

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
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		preemptible: preemptible
	}
}

task run_cumulus_demuxEM {
	File input_adt_csv
	File input_raw_gene_bc_matrices_h5
	String output_dir
	String output_name
	String cumulus_version
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
    String docker_registry

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['pegasus', 'demuxEM', '${input_adt_csv}', '${input_raw_gene_bc_matrices_h5}', '${output_name}', '-p', '${num_cpu}']
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

		gsutil -q cp ${output_name}_demux.h5sc ${output_dir}/${output_name}/
		gsutil -q cp ${output_name}_ADTs.h5ad ${output_dir}/${output_name}/
		gsutil -q cp ${output_name}_demux.h5ad ${output_dir}/${output_name}/
		gsutil -q -m cp ${output_name}.*.pdf ${output_dir}/${output_name}/
		# mkdir -p ${output_dir}/${output_name}
		# cp ${output_name}_demux.h5sc ${output_dir}/${output_name}/
		# cp ${output_name}_ADTs.h5ad ${output_dir}/${output_name}/
		# cp ${output_name}_demux.h5ad ${output_dir}/${output_name}/
		# cp ${output_name}.*.pdf ${output_dir}/${output_name}/
	}

	output {
		String output_folder = "${output_dir}/${output_name}"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: num_cpu
		preemptible: preemptible
	}
}

task run_cumulus_merge_rna_adt {
	File input_raw_gene_bc_matrices_h5
	File input_adt_csv
	File? antibody_control_csv
	String output_dir
	String output_name
	String cumulus_version
	String zones
	String memory
	Int disk_space
	Int preemptible
	String docker_registry

	command {
		set -e
		export TMPDIR=/tmp
		monitor_script.sh > monitoring.log &

		python <<CODE
		from subprocess import check_call
		call_args = ['pegasus', 'merge_rna_adt', '${input_raw_gene_bc_matrices_h5}', '${input_adt_csv}', '${output_name}']
		if '${antibody_control_csv}' is not '':
			call_args.extend(['--antibody-control-csv', '${antibody_control_csv}'])
		print(' '.join(call_args))
		check_call(call_args)
		CODE

		gsutil -q cp ${output_name}.h5sc ${output_dir}/${output_name}/
		# mkdir -p ${output_dir}/${output_name}
		# cp ${output_name}.h5sc ${output_dir}/${output_name}/
	}

	output {
		String output_folder = "${output_dir}/${output_name}"
		File monitoringLog = "monitoring.log"
	}

	runtime {
		docker: "${docker_registry}cumulus:${cumulus_version}"
		zones: zones
		memory: memory
		bootDiskSizeGb: 12
		disks: "local-disk ${disk_space} HDD"
		cpu: 1
		preemptible: preemptible
	}
}
