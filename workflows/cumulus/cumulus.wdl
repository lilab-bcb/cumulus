version 1.0

import "cumulus_tasks.wdl" as tasks

workflow cumulus {
    input {
        # Input file: can be either a csv-formatted file containing information of each scRNA-Seq run or a single input file
        File input_file
        # Google bucket and directory name.
        String output_directory
        # Results name prefix and subdirectory name.
        String output_name

        # Pegasus version, default to "1.4.3"
        String pegasus_version = "1.4.3"
        # Docker registry to use
        String docker_registry = "quay.io/cumulus"
        # Google cloud zones, default to "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        String zones = "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
        # Number of cpus per cumulus job
        Int num_cpu = 32
        # Memory size string
        String memory = "200G"
        # Total disk space
        Int disk_space = 100
        # Number of preemptible tries

        # Backend
        String backend = "gcp"
        Int preemptible = 2
        # Number of maximum retries when running on AWS
        Int awsMaxRetries = 5

        # If sample count matrix is in either DGE, mtx, csv, tsv or loom format and there is no Reference column in the csv_file, use default_reference as the reference string.
        String? default_reference


        # for aggregate_matrices

        # Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value. Multiple restrictions are separated by ';'
        String? restrictions
        # Specify a comma-separated list of outputted attributes. These attributes should be column names in the csv file
        String? attributes
        # If we have demultiplexed data, turning on this option will make cumulus only include barcodes that are predicted as singlets
        Boolean select_only_singlets = false
        # Remap singlet names using <remap_string>, where <remap_string> takes the format "new_name_i:old_name_1,old_name_2;new_name_ii:old_name_3;...". For example, if we hashed 5 libraries from 3 samples sample1_lib1, sample1_lib2, sample2_lib1, sample2_lib2 and sample3, we can remap them to 3 samples using this string: "sample1:sample1_lib1,sample1_lib2;sample2:sample2_lib1,sample2_lib2". After that, original singlet names will be kept in metadate field with key name 'assignment.orig'.
        String? remap_singlets
        # If select singlets, only select singlets in the <subset_string>, which takes the format "name1,name2,...". Note that if --remap-singlets is specified, subsetting happens after remapping. For example, we can only select singlets from sampe 1 and 3 using "sample1,sample3".
        String? subset_singlets
        # Only keep barcodes with at least this number of expressed genes
        Int minimum_number_of_genes = 100
        # If inputs are dropseq data
        Boolean is_dropseq = false

        # for cluster

        # Specify the cell barcode attribute to represent different samples.
        String? channel
        # Specify cell barcode attributes to be popped out.
        String? black_list
        # If input are raw 10x matrix, which include all barcodes, perform a pre-filtration step to keep the data size small. In the pre-filtration step, only keep cells with at least <number> of genes. [default: 100]
        Int? min_genes_before_filtration
        # Focus analysis on Unimodal data with <keys>. <keys> is a comma-separated list of keys. If None, the self._selected will be the focused one.
        String? focus
        # Append Unimodal data <key> to any <keys> in --focus.
        String? append
        # If write cell and gene filtration results as a spreadsheet. [default: true]
        Boolean output_filtration_results = true
        # If plot filtration results as PDF files. [default: true]
        Boolean plot_filtration_results = true
        # Figure size for filtration plots. <figsize> is a comma-separated list of two numbers, the width and height of the figure (e.g. 6,4).
        String? plot_filtration_figsize
        # Output seurat-compatible h5ad file. Caution: File size might be large, do not turn this option on for large data sets. [default: false]
        Boolean output_h5ad = true
        # If output loom-formatted file [default: false]
        Boolean? output_loom
        # Only keep cells with at least <number> of genes. [default: 500]
        Int? min_genes
        # Only keep cells with less than <number> of genes. [default: 6000]
        Int? max_genes
        # Only keep cells with at least <number> of UMIs. [default: None]
        Int? min_umis
        # Only keep cells with less than <number> of UMIs. [default: None]
        Int? max_umis
        # Prefix for mitochondrial genes. [default: GRCh38:MT-; mm10:mt-]
        String? mito_prefix
        # Only keep cells with mitochondrial percent less than <percent>%. [default: 20]
        Float? percent_mito
        # Only use genes that are expressed in at <percent>% of cells to select variable genes. [default: 0.05]
        Float? gene_percent_cells
        # Total counts per cell after normalization. [default: 1e5]
        Float? counts_per_cell_after
        # Highly variable feature selection method. <flavor> can be 'pegasus' or 'Seurat'. [default: pegasus]
        String? select_hvf_flavor
        # Select top <nfeatures> highly variable features. If <flavor> is 'Seurat' and <nfeatures> is 'None', select HVGs with z-score cutoff at 0.5. [default: 2000]
        Int? select_hvf_ngenes
        # Plot highly variable feature selection.
        Boolean? plot_hvf
        # Do not select highly variable features. [default: false]
        Boolean? no_select_hvf
        # If correct batch effects [default: false]
        Boolean? correct_batch_effect
        # Batch correction method, can be either 'L/S' for location/scale adjustment algorithm (Li and Wong. The analysis of Gene Expression Data 2003),
        # 'harmony' for Harmony (Korsunsky et al. Nature Methods 2019),
        # or 'inmf' for integrative NMF (Yang and Michailidis Bioinformatics 2016, Welch et al. Cell 2019, Gao et al. Natuer Biotechnology 2021). [default: harmony]
        String? correction_method
        # Batch correction assumes the differences in gene expression between channels are due to batch effects. However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. In this case, we will only perform batch correction for channels within each group. This option defines the groups. If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>. <expression> takes the form of either ‘attr’, or ‘attr1+attr2+…+attrn’, or ‘attr=value11,…,value1n_1;value21,…,value2n_2;…;valuem1,…,valuemn_m’. In the first form, ‘attr’ should be an existing sample attribute, and groups are defined by ‘attr’. In the second form, ‘attr1’,…,’attrn’ are n existing sample attributes and groups are defined by the Cartesian product of these n attributes. In the last form, there will be m + 1 groups. A cell belongs to group i (i > 0) if and only if its sample attribute ‘attr’ has a value among valuei1,…,valuein_i. A cell belongs to group 0 if it does not belong to any other groups.
        String? batch_group_by
        # Coefficient of regularization for iNMF. [default: 5.0]
        Int? inmf_lambda
        # Random number generator seed. [default: 0]
        Int? random_state
        # Calculate signature scores for gene sets in <sig_list>. <sig_list> is a comma-separated list of strings. Each string should either be a <GMT_file> or one of 'cell_cycle_human', 'cell_cycle_mouse', 'gender_human', 'gender_mouse', 'mitochondrial_genes_human', 'mitochondrial_genes_mouse', 'ribosomal_genes_human' and 'ribosomal_genes_mouse'
        String? calc_signature_scores
        # Number of PCs. [default: 50]
        Int? nPC
        # Compute nonnegative matrix factorization (NMF) on highly variable features.
        Boolean? run_nmf
        # Number of NMF components. IF iNMF is used for batch correction, this parameter also sets iNMF number of components. [default: 20]
        Int nmf_n = 20
        # Number of neighbors used for constructing affinity matrix. [default: 100]
        Int? knn_K
        # For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads. [default: false]
        Boolean? knn_full_speed
        # Calculate diffusion map.
        Boolean? run_diffmap
        # Number of diffusion components. [default: 100]
        Int? diffmap_ndc
        # Maximum time stamp to search for the knee point. [default: 5000]
        Int? diffmap_maxt
        # Run louvain clustering algorithm.
        Boolean run_louvain = true
        # Resolution parameter for the louvain clustering algorithm. [default: 1.3]
        Float? louvain_resolution
        # Louvain cluster label name in AnnData. [default: louvain_labels]
        String? louvain_class_label
        # Run leiden clustering algorithm.
        Boolean? run_leiden
        # Resolution parameter for the leiden clustering algorithm. [default: 1.3]
        Float? leiden_resolution
        # Number of iterations of running the Leiden algorithm. If negative, run Leiden iteratively until no improvement. [default: -1]
        Int? leiden_niter
        # Leiden cluster label name in AnnData. [default: leiden_labels]
        String? leiden_class_label
        # Run spectral louvain clustering algorithm.
        Boolean? run_spectral_louvain
        # Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. [default: diffmap]
        String? spectral_louvain_basis
        # Resolution parameter for louvain. [default: 1.3]
        Float? spectral_louvain_resolution
        # Spectral louvain label name in AnnData. [default: spectral_louvain_labels]
        String? spectral_louvain_class_label
        # Run spectral leiden clustering algorithm.
        Boolean? run_spectral_leiden
        # Basis used for KMeans clustering. Can be 'pca' or 'diffmap'. [default: diffmap]
        String? spectral_leiden_basis
        # Resolution parameter for leiden. [default: 1.3]
        Float? spectral_leiden_resolution
        # Approximated leiden label name in AnnData. [default: spectral_louvain_labels]
        String? spectral_leiden_class_label
        # Run FIt-SNE package to compute t-SNE embeddings for visualization.
        Boolean? run_tsne
        # tSNE’s perplexity parameter. [default: 30]
        Float? tsne_perplexity
        # <choice> can be either 'random' or 'pca'. 'random' refers to random initialization. 'pca' refers to PCA initialization as described in (CITE Kobak et al. 2019) [default: pca]
        String? tsne_initialization
        # Run umap for visualization.
        Boolean run_umap = true
        # Run umap on diffusion components. [default: 15]
        Int? umap_K
        # Umap parameter. [default: 0.5]
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
        Boolean? run_net_umap
        # Output basis for net-UMAP. [default: net_umap]
        String? net_umap_out_basis
        # Run net FLE.
        Boolean? run_net_fle
        # Output basis for net-FLE. [default: net_fle]
        String? net_fle_out_basis
        # Infer doublets using the method described in https://github.com/klarman-cell-observatory/pegasus/raw/master/doublet_detection.pdf. Obs attribute 'doublet_score' stores Scrublet-like doublet scores and attribute 'demux_type' stores 'doublet/singlet' assignments.
        Boolean? infer_doublets
        # The expected doublet rate per sample. By default, calculate the expected rate based on number of cells from the 10x multiplet rate table.
        Float? expected_doublet_rate
        # <attr> refers to a cluster attribute containing cluster labels (e.g. 'louvain_labels'). Doublet clusters will be marked based on <attr> with the following criteria: passing the Fisher's exact test and having >= 50% of cells identified as doublets. By default, the first computed cluster attribute in the list of leiden, louvain, spectral_ledein and spectral_louvain is used.
        String? doublet_cluster_attribute
        # Input data contain both RNA and CITE-Seq modalities. This will set --focus to be the RNA modality and --append to be the CITE-Seq modality. In addition, 'ADT-' will be added in front of each antibody name to avoid name conflict with genes in the RNA modality.
        Boolean? citeseq
        # For high quality cells kept in the RNA modality, generate a UMAP based on their antibody expression.
        Boolean? citeseq_umap
        # <list> is a comma-separated list of antibodies to be excluded from the UMAP calculation (e.g. Mouse-IgG1,Mouse-IgG2a).
        String? citeseq_umap_exclude

        # for de_analysis and annotate_cluster

        # If perform de analysis
        Boolean perform_de_analysis = true
        # Specify the cluster labels used for differential expression analysis. [default: louvain_labels]
        String? cluster_labels
        # Control false discovery rate at <alpha>. [default: 0.05]
        Float? alpha
        # Calculate Welch's t-test.
        Boolean t_test = false
        # Calculate Fisher’s exact test.
        Boolean fisher = false

        # If also detect markers using LightGBM
        Boolean? find_markers_lightgbm
        # Remove ribosomal genes with either RPL or RPS as prefixes
        Boolean? remove_ribo
        # Only report genes with a feature importance score (in gain) of at least <gain>. [default: 1.0]
        Float? min_gain

        # If also annotate cell types for clusters based on DE results.
        Boolean? annotate_cluster
        # Organism, could either be "human_immune", "mouse_immune", "human_brain", "mouse_brain", "human_lung" or a JSON file describing the markers. [default: human_immune]
        String? organism
        # DE test to use to infer cell types, could be either "mwu", "t", or "fisher". [default: mwu]
        String? annotate_de_test
        # Minimum cell type score to report a potential cell type. [default: 0.5]
        Float? minimum_report_score


        # for plot

        # Takes the format of "label:attr,label:attr,...,label:attr". If non-empty, generate composition plot for each "label:attr" pair. "label" refers to cluster labels and "attr" refers to sample conditions.
        String? plot_composition
        # Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored tSNEs side by side.
        String? plot_tsne
        # Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored UMAPs side by side.
        String? plot_umap
        # Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored FLEs side by side.
        String? plot_fle
        # Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored tSNEs side by side based on net tSNE result.
        String? plot_net_umap
        # Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored FLEs side by side based on net FLE result.
        String? plot_net_fle
        # Takes the format of "attr,attr,...,attr". If non-empty, plot attr colored UMAPs side by side based on CITE-Seq UMAP result.
        String? plot_citeseq_umap

        # for cirro_output
        # If generate Cirrocumulus inputs
        Boolean generate_cirro_inputs = false

        # for scp_output

        # If generate outputs required by single cell portal.
        Boolean generate_scp_outputs = false
        # Output dense expression matrix instead.
        Boolean output_dense = false
    }

    # Output directory, with trailing slashes stripped
    String output_directory_stripped = sub(output_directory, "[/\\s]+$", "")
    # If input file is a sample sheet in csv format.
    Boolean is_sample_sheet = sub(input_file, "^.+\\.csv$", "CSV") == "CSV"

    if (is_sample_sheet) {
        call tasks.run_cumulus_aggregate_matrices as aggregate_matrices {
            input:
                input_count_matrix_csv = input_file,
                output_directory = output_directory_stripped,
                output_name = output_name,
                restrictions = restrictions,
                attributes = attributes,
                default_reference = default_reference,
                select_only_singlets = select_only_singlets,
                remap_singlets = remap_singlets,
                subset_singlets = subset_singlets,
                minimum_number_of_genes = minimum_number_of_genes,
                pegasus_version = pegasus_version,
                docker_registry = docker_registry,
                zones = zones,
                memory = memory,
                disk_space = disk_space,
                preemptible = preemptible,
                awsMaxRetries = awsMaxRetries,
                backend = backend
        }
    }



    call tasks.run_cumulus_cluster as cluster {
        input:
            input_file = if is_sample_sheet then select_first([aggregate_matrices.output_zarr]) else input_file,
            output_directory = output_directory_stripped,
            output_name = output_name,
            channel = channel,
            black_list = black_list,
            min_genes_before_filtration = min_genes_before_filtration,
            select_singlets = if is_sample_sheet then false else select_only_singlets,
            remap_singlets = if is_sample_sheet then "" else remap_singlets,
            subset_singlets = if is_sample_sheet then "" else subset_singlets,
            genome = default_reference,
            focus = focus,
            append = append,
            output_filtration_results = output_filtration_results,
            plot_filtration_results = plot_filtration_results,
            plot_filtration_figsize = plot_filtration_figsize,
            output_h5ad = output_h5ad,
            output_loom = output_loom,
            min_genes = min_genes,
            max_genes = max_genes,
            min_umis = min_umis,
            max_umis = max_umis,
            mito_prefix = mito_prefix,
            percent_mito = percent_mito,
            gene_percent_cells = gene_percent_cells,
            counts_per_cell_after = counts_per_cell_after,
            select_hvf_flavor = select_hvf_flavor,
            select_hvf_ngenes = select_hvf_ngenes,
            no_select_hvf = no_select_hvf,
            plot_hvf = plot_hvf,
            correct_batch_effect = correct_batch_effect,
            correction_method = correction_method,
            batch_group_by = batch_group_by,
            inmf_lambda = inmf_lambda,
            random_state = random_state,
            gene_signature_set = calc_signature_scores,
            nPC = nPC,
            run_nmf = run_nmf,
            nmf_n = nmf_n,
            knn_K = knn_K,
            knn_full_speed = knn_full_speed,
            run_diffmap = run_diffmap,
            diffmap_ndc = diffmap_ndc,
            diffmap_maxt = diffmap_maxt,
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
            tsne_initialization = tsne_initialization,
            run_umap = run_umap,
            umap_K = umap_K,
            umap_min_dist = umap_min_dist,
            umap_spread = umap_spread,
            run_fle = run_fle,
            fle_K = fle_K,
            fle_target_change_per_node = fle_target_change_per_node,
            fle_target_steps = fle_target_steps,
            net_down_sample_fraction = net_down_sample_fraction,
            run_net_umap = run_net_umap,
            net_umap_out_basis = net_umap_out_basis,
            run_net_fle = run_net_fle,
            net_fle_out_basis = net_fle_out_basis,
            infer_doublets = infer_doublets,
            expected_doublet_rate = expected_doublet_rate,
            doublet_cluster_attribute = doublet_cluster_attribute,
            citeseq = citeseq,
            citeseq_umap = citeseq_umap,
            citeseq_umap_exclude = citeseq_umap_exclude,
            pegasus_version = pegasus_version,
            docker_registry = docker_registry,
            zones = zones,
            num_cpu = num_cpu,
            memory = memory,
            disk_space = disk_space,
            preemptible = preemptible,
            awsMaxRetries = awsMaxRetries,
            backend = backend
    }

    if (length(cluster.output_h5ad_file) > 0) {
        scatter (focus_h5ad in cluster.output_h5ad_file) {
            String focus_prefix = basename(focus_h5ad, ".h5ad")

            if (perform_de_analysis) {
                call tasks.run_cumulus_de_analysis as de_analysis {
                    input:
                        input_h5ad = focus_h5ad,
                        output_directory = output_directory_stripped + '/' + output_name,
                        output_name = focus_prefix,
                        labels = cluster_labels,
                        alpha = alpha,
                        t_test = t_test,
                        fisher = fisher,
                        find_markers_lightgbm = find_markers_lightgbm,
                        remove_ribo = remove_ribo,
                        min_gain = min_gain,
                        random_state = random_state,
                        annotate_cluster = annotate_cluster,
                        annotate_de_test = annotate_de_test,
                        organism = organism,
                        minimum_report_score = minimum_report_score,
                        pegasus_version = pegasus_version,
                        zones = zones,
                        num_cpu = num_cpu,
                        memory = memory,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        docker_registry = docker_registry,
                        backend = backend
                }
            }

            Boolean do_plot = defined(plot_composition) || defined(plot_tsne) || defined(plot_umap) || defined(plot_fle) || defined(plot_net_umap) || defined(plot_net_fle)

            String plot_nmf = if (defined(correct_batch_effect) && defined(correction_method) && correction_method == 'inmf') then 'inmf' else (if (defined(run_nmf) && run_nmf) then 'nmf' else '')

            if (do_plot) {
                call tasks.run_cumulus_plot as plot {
                    input:
                        input_h5ad = focus_h5ad,
                        output_directory = output_directory_stripped + '/' + output_name,
                        output_name = focus_prefix,
                        plot_composition = plot_composition,
                        plot_tsne = plot_tsne,
                        plot_umap = plot_umap,
                        plot_fle = plot_fle,
                        plot_net_umap = plot_net_umap,
                        plot_net_fle = plot_net_fle,
                        plot_citeseq_umap = plot_citeseq_umap,
                        plot_nmf = plot_nmf,
                        nmf_n = nmf_n,
                        pegasus_version = pegasus_version,
                        docker_registry = docker_registry,
                        zones = zones,
                        memory = memory,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend
                }
            }

            if (generate_cirro_inputs) {
                call tasks.run_cumulus_cirro_output as cirro_output {
                    input:
                        input_h5ad = focus_h5ad,
                        output_directory = output_directory_stripped + '/' + output_name,
                        output_name = focus_prefix,
                        pegasus_version = pegasus_version,
                        docker_registry = docker_registry,
                        zones = zones,
                        memory = memory,
                        disk_space = disk_space,
                        num_cpu = num_cpu,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend
                }
            }

            if (generate_scp_outputs) {
                call tasks.run_cumulus_scp_output as scp_output {
                    input:
                        input_h5ad = focus_h5ad,
                        output_directory = output_directory_stripped + '/' + output_name,
                        output_name = focus_prefix,
                        output_dense = output_dense,
                        pegasus_version = pegasus_version,
                        docker_registry = docker_registry,
                        zones = zones,
                        memory = memory,
                        disk_space = disk_space,
                        preemptible = preemptible,
                        awsMaxRetries = awsMaxRetries,
                        backend = backend
                }
            }
        }
    }


    output {
        File? output_aggr_zarr = aggregate_matrices.output_zarr
        File output_zarr = cluster.output_zarr
        File output_cluster_log = cluster.output_log
        Array[File] output_h5ad_file = cluster.output_h5ad_file
        Array[File] output_filt_xlsx = cluster.output_filt_xlsx
        Array[File] output_filt_plot = cluster.output_filt_plot
        Array[File] output_hvf_plot = cluster.output_hvf_plot
        Array[File] output_dbl_plot = cluster.output_dbl_plot
        Array[File] output_loom_file = cluster.output_loom_file
        Array[File?]? output_de_h5ad = de_analysis.output_de_h5ad
        Array[File?]? output_de_xlsx =  de_analysis.output_de_xlsx
        Array[Array[File]?]? output_markers_xlsx =  de_analysis.output_markers_xlsx
        Array[Array[File]?]? output_anno_file =  de_analysis.output_anno_file
        Array[Array[File]?]? output_pdfs = plot.output_pdfs
        Array[Array[File]?]? output_htmls = plot.output_htmls
        Array[String?]? output_cirro_folder = cirro_output.output_cirro_folder
        Array[Array[File]?]? output_scp_files= scp_output.output_scp_files
    }
}
