import "https://api.firecloud.org/ga4gh/v1/tools/cumulus:cumulus_tasks/versions/3/plain-WDL/descriptor" as tasks
# import "cumulus_tasks.wdl" as tasks

workflow cumulus_subcluster {
	File input_h5ad
	# google bucket, subdirectory name and results name prefix
	String output_name
	# Specify which cells will be included in the subcluster analysis. This field contains one or more <subset_selection> strings separated by ';'. Each <subset_selection> string takes the format of ‘attr:value,…,value’, which means select cells with attr in the values. If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings.
	String subset_selections

	# cumulus version, default to "0.11.0"
	String? cumulus_version = "0.11.0"
	# Docker registry to use
	String? docker_registry = "cumulusprod/"
	# Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	String? zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	# Number of cpus per cumulus job
	Int? num_cpu = 64
	# Memory size string
	String? memory = "200G"
	# Total disk space
	Int? disk_space = 100
	# Number of preemptible tries 
	Int? preemptible = 2


	String out_name = basename(output_name)


	# for subcluster

	# If correct batch effects [default: false]
	Boolean? correct_batch_effect
	# Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either 'attr', or 'attr1+attr2+...+attrn', or 'attr=value11,...,value1n_1;value21,...,value2n_2;...;valuem1,...,valuemn_m'. In the first form, 'attr' should be an existing sample attribute, and groups are defined by 'attr'. In the second form, 'attr1',...,'attrn' are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute 'attr' has a value among valuei1,...,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.
	String? batch_group_by
	# If output loom-formatted file [default: false]
	Boolean? output_loom
	# If output parquet-formatted file [default: false]
	Boolean? output_parquet
	# Highly variable feature selection method. <flavor> can be 'pegasus' or 'Seurat'. [default: pegasus]
	String? select_hvf_flavor
	# Select top <nfeatures> highly variable features. If <flavor> is 'Seurat' and <nfeatures> is 'None', select HVGs with z-score cutoff at 0.5. [default: 2000]
	Int? select_hvf_ngenes
	# Do not select highly variable features.
	Boolean? no_select_hvf
	# Random number generator seed. [default: 0]
	Int? random_state
	# Number of principal components. [default: 50]
	Int? nPC
	# Number of nearest neighbors for building kNN graph. [default: 100]
	Int? knn_K
	# For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads. [default: false]
	Boolean? knn_full_speed
	# Calculate diffusion maps.
	Boolean? run_diffmap
	# Number of diffusion components. [default: 50]
	Int? diffmap_ndc
	# Maximum time stamp to search for the knee point. [default: 5000]
	Float? diffmap_maxt
	# Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes.
	String? calculate_pseudotime
	# Run louvain clustering algorithm.
	Boolean? run_louvain = true
	# Resolution parameter for the louvain clustering algorithm. [default: 1.3]
	Float? louvain_resolution
	# Louvain cluster label name in AnnData. [default: louvain_labels]
	String? louvain_class_label
	# Run leiden clustering algorithm.
	Boolean? run_leiden
	# Resolution parameter for the leiden clustering algorithm. [default: 1.3]
	Float? leiden_resolution
	# Number of iterations of running the Leiden algorithm. If <niter> is negative, run Leiden iteratively until no improvement. [default: -1]
	Int? leiden_niter
	# Leiden cluster label name in AnnData. [default: leiden_labels]
	String? leiden_class_label
	# Run spectral louvain clustering algorithm.
	Boolean? run_spectral_louvain
	# Basis used for KMeans clustering. Can be 'pca', or 'diffmap'. [default: diffmap]
	String? spectral_louvain_basis
	# Resolution parameter for louvain. [default: 1.3]
	Float? spectral_louvain_resolution
	# Spectral louvain label name in AnnData. [default: spectral_louvain_labels]
	String? spectral_louvain_class_label
	# Run spectral leiden clustering algorithm.
	Boolean? run_spectral_leiden
	# Basis used for KMeans clustering. Can be 'pca', or 'diffmap'. [default: diffmap]
	String? spectral_leiden_basis
	# Resolution parameter for leiden. [default: 1.3]
	Float? spectral_leiden_resolution
	# Spectral leiden label name in AnnData. [default: spectral_leiden_labels]
	String? spectral_leiden_class_label
	# Run multi-core tSNE for visualization.
	Boolean? run_tsne
	# tSNE’s perplexity parameter. [default: 30]
	Float? tsne_perplexity
	# Run FItSNE for visualization.
	Boolean? run_fitsne = true
	# Run umap for visualization.
	Boolean? run_umap
	# K neighbors for umap. [default: 15]
	Int? umap_K
	# Umap parameter. [default: 0.1]
	Float? umap_min_dist
	# Umap parameter. [default: 1.0]
	Float? umap_spread
	# Run force-directed layout embedding.
	Boolean? run_fle
	# K neighbors for building graph for FLE. [default: 50]
	Int? fle_K
	# Target change per node to stop forceAtlas2. [default: 2.0]
	Float? fle_target_change_per_node
	# Maximum number of iterations before stopping the forceAtlas2 algoritm. [default: 5000]
	Int? fle_target_steps
	# Down sampling fraction for net-related visualization. [default: 0.1]
	Float? net_down_sample_fraction
	# Run net tSNE for visualization.
	Boolean? run_net_tsne
	# Output basis for net-tSNE. [default: net_tsne]
	String? net_tsne_out_basis
	# Run net UMAP.
	Boolean? run_net_umap
	# Output basis for net-UMAP. [default: net_umap]
	String? net_umap_out_basis
	# Run net FLE.
	Boolean? run_net_fle
	# Output basis for net-FLE. [default: net_fle]
	String? net_fle_out_basis


	# for de_analysis and annotate_cluster

	# If perform de analysis
	Boolean perform_de_analysis = true
	# Specify the cluster labels used for differential expression analysis. [default: louvain_labels]
	String? cluster_labels
	# Control false discovery rate at <alpha>. [default: 0.05]
	Float? alpha
	# Calculate area under ROC (AUROC) and area under Precision-Recall (AUPR).
	Boolean? auc = true
	# Calculate Welch's t-test.
	Boolean? t_test = true
	# Calculate Fisher’s exact test.
	Boolean? fisher = true
	# Calculate Mann-Whitney U test.
	Boolean? mwu

	# If also detect markers using LightGBM
	Boolean? find_markers_lightgbm
	# Remove ribosomal genes with either RPL or RPS as prefixes
	Boolean? remove_ribo
	# Only report genes with a feature importance score (in gain) of at least <gain>. [default: 1.0]
	Float? min_gain

	# If also annotate cell types for clusters based on DE results.
	Boolean? annotate_cluster
	# DE test to use to infer cell types, could be either "t", "fisher", or "mwu". [default: t]
	String? annotate_de_test
	# Organism, could either be "human_immune" or "mouse_immune" or "mouse_brain" [default: human_immune]
	String? organism
	# Minimum cell type score to report a potential cell type. [default: 0.5]
	Float? minimum_report_score


	# for plot

	# Takes the format of "label:attr,label:attr,...,label:attr". If non-empty, generate composition plot for each "label:attr" pair. "label" refers to cluster labels and "attr" refers to sample conditions.
	String? plot_composition
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored tSNEs side by side. 
	String? plot_tsne
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored FItSNEs side by side.
	String? plot_fitsne
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored UMAPs side by side.
	String? plot_umap
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored FLEs side by side.
	String? plot_fle
	# Takes the format of "attr,attr,...,attr". If non-empty, generate attr colored 3D interactive plot. The 3 coordinates are the first 3 PCs of all diffusion components.
	String? plot_diffmap
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored tSNEs side by side based on net tSNE result.
	String? plot_net_tsne
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored UMAPs side by side based on net UMAP result.
	String? plot_net_umap
	# Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored FLEs side by side based on net FLE result.
	String? plot_net_fle


	# for scp_output

	# If generate outputs required by single cell portal
	Boolean generate_scp_outputs = false
	# Organism, could either be "human_immune", "mouse_immune", "human_brain", "mouse_brain" or a JSON file describing the markers. [default: human_immune]
	Boolean output_dense = false


	call tasks.run_cumulus_subcluster as subcluster {
		input:
			input_h5ad = input_h5ad,
			output_name = out_name,
			subset_selections = subset_selections,
			correct_batch_effect = correct_batch_effect,
			batch_group_by = batch_group_by,
			output_loom = output_loom,
			output_parquet = output_parquet,
			select_hvf_flavor = select_hvf_flavor,
			select_hvf_ngenes = select_hvf_ngenes,
			no_select_hvf = no_select_hvf,
			random_state = random_state,
			nPC = nPC,
			knn_K = knn_K,
			knn_full_speed = knn_full_speed,
			run_diffmap = run_diffmap,
			diffmap_ndc = diffmap_ndc,
			diffmap_maxt = diffmap_maxt,
			calculate_pseudotime = calculate_pseudotime,
			run_louvain = run_louvain,
			louvain_resolution = louvain_resolution,
			louvain_class_label = louvain_class_label,
			run_leiden = run_leiden,
			leiden_resolution = leiden_resolution,
			leiden_niter = leiden_niter,
			leiden_class_label = leiden_class_label,
			run_spectral_louvain = run_spectral_louvain,
			spectral_louvain_basis = spectral_louvain_basis,
			spectral_louvain_resolution = spectral_louvain_resolution,
			spectral_louvain_class_label = spectral_louvain_class_label,
			run_spectral_leiden = run_spectral_leiden,
			spectral_leiden_basis = spectral_leiden_basis,
			spectral_leiden_resolution = spectral_leiden_resolution,
			spectral_leiden_class_label = spectral_leiden_class_label,
			run_tsne = run_tsne,
			tsne_perplexity = tsne_perplexity,
			run_fitsne = run_fitsne,
			run_umap = run_umap,
			umap_K = umap_K,
			umap_min_dist = umap_min_dist,
			umap_spread = umap_spread,
			run_fle = run_fle,
			fle_K = fle_K,
			fle_target_change_per_node = fle_target_change_per_node,
			fle_target_steps = fle_target_steps,
			net_down_sample_fraction = net_down_sample_fraction,
			run_net_tsne = run_net_tsne,
			net_tsne_out_basis = net_tsne_out_basis,
			run_net_umap = run_net_umap,
			net_umap_out_basis = net_umap_out_basis,
			run_net_fle = run_net_fle,
			net_fle_out_basis = net_fle_out_basis,
			cumulus_version = cumulus_version,
			zones = zones,
			num_cpu = num_cpu,
			memory = memory,
			disk_space = disk_space,
			preemptible = preemptible,
			docker_registry = docker_registry
	}

	if (perform_de_analysis) {
		call tasks.run_cumulus_de_analysis as de_analysis {
			input:
				input_h5ad = subcluster.output_h5ad,
				output_name = out_name,
				labels = cluster_labels,
				alpha = alpha,
				auc = auc,
				t_test = t_test,
				fisher = fisher,
				mwu = mwu,
				find_markers_lightgbm = find_markers_lightgbm,
				remove_ribo = remove_ribo,
				min_gain = min_gain,
				random_state = random_state,
				annotate_cluster = annotate_cluster,
				annotate_de_test = annotate_de_test,
				organism = organism,
				minimum_report_score = minimum_report_score,
				cumulus_version = cumulus_version,
				zones = zones,
				num_cpu = num_cpu,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible,
				docker_registry = docker_registry
		}
	}

	if (defined(plot_composition) || defined(plot_tsne) || defined(plot_fitsne) || defined(plot_umap) || defined(plot_fle) || defined(plot_diffmap) || defined(plot_net_tsne) || defined(plot_net_umap) || defined(plot_net_fle)) {
		call tasks.run_cumulus_plot as plot {
			input:
				input_h5ad = subcluster.output_h5ad,
				output_name = out_name,
				plot_composition = plot_composition,
				plot_tsne = plot_tsne,
				plot_fitsne = plot_fitsne,
				plot_umap = plot_umap,
				plot_fle = plot_fle,
				plot_diffmap = plot_diffmap,
				plot_net_tsne = plot_net_tsne,
				plot_net_umap = plot_net_umap,
				plot_net_fle = plot_net_fle,
				cumulus_version = cumulus_version,
				zones = zones,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible,
				docker_registry = docker_registry
		}
	}

	if (generate_scp_outputs) {
		call tasks.run_cumulus_scp_output as scp_output {
			input:
				input_h5ad = subcluster.output_h5ad,
				output_name = out_name,
				output_dense = output_dense,
				cumulus_version = cumulus_version,
				zones = zones,
				memory = memory,
				disk_space = disk_space,
				preemptible = preemptible,
				docker_registry = docker_registry
		}
	}

	call tasks.organize_results {
		input:
			output_name = output_name,
			output_h5ad = subcluster.output_h5ad,
			output_loom_file = subcluster.output_loom_file,
			output_parquet_file = subcluster.output_parquet_file,
			output_de_h5ad = de_analysis.output_de_h5ad,
			output_de_xlsx = de_analysis.output_de_xlsx,
			output_markers_xlsx = de_analysis.output_markers_xlsx,
			output_anno_file = de_analysis.output_anno_file,
			output_pdfs = plot.output_pdfs,
			output_htmls = plot.output_htmls,
			output_scp_files = scp_output.output_scp_files,
			cumulus_version = cumulus_version,
			zones = zones,
			disk_space = disk_space,
			preemptible = preemptible,
			docker_registry = docker_registry
	}
}
