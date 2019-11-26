Run Cumulus for sc/snRNA-Seq data analysis
---------------------------------------------------------------------------------


Run Cumulus analysis
~~~~~~~~~~~~~~~~~~~~~~~

Prepare Input Data
^^^^^^^^^^^^^^^^^^^^^^^

Case One: Sample Sheet
+++++++++++++++++++++++++

Follow the steps below to run **cumulus** on Terra_.

#. Create a sample sheet, **count_matrix.csv**, which describes the metadata for each 10x channel. The sample sheet should at least contain 2 columns --- *Sample* and *Location*. *Sample* refers to sample names and *Location* refers to the location of the channel-specific count matrix in either of

  - 10x format with v2 chemistry. For example, ``gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5``.
  - 10x format with v3 chemistry. For example, ``gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_feature_bc_matrices.h5``.
  - Drop-seq format. For example, ``gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_2/sample_2.umi.dge.txt.gz``. 

You are free to add any other columns and these columns will be used in selecting channels for futher analysis. In the example below, we have *Source*, which refers to the tissue of origin, *Platform*, which refers to the sequencing platform, *Donor*, which refers to the donor ID, and *Reference*, which refers to the reference genome.

	Example::

		Sample,Source,Platform,Donor,Reference,Location
		sample_1,bone_marrow,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5
		sample_2,bone_marrow,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_2/filtered_gene_bc_matrices_h5.h5
		sample_3,pbmc,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_3/filtered_feature_bc_matrices.h5
		sample_4,pbmc,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_4/filtered_feature_bc_matrices.h5

	If you ran **cellranger_workflow** ahead, you should already obtain a template **count_matrix.csv** file that you can modify from **generate_count_config**'s outputs.

#. Upload your sample sheet to the workspace.  

	Example::
	
		gsutil cp /foo/bar/projects/my_count_matrix.csv gs://fc-e0000000-0000-0000-0000-000000000000/

	where ``/foo/bar/projects/my_count_matrix.csv`` is the path to your sample sheet in local machine, and ``gs://fc-e0000000-0000-0000-0000-000000000000/`` is the location on Google bucket to hold it.

#. Import *cumulus* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *cumulus* workflow is under ``Broad Methods Repository`` with name "**cumulus/cumulus**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cumulus* workflow in the drop-down menu.

#. In your workspace, open ``cumulus`` in ``WORKFLOWS`` tab. Select ``Process single workflow from files`` as below

	.. image:: images/single_workflow.png

   and click the ``SAVE`` button.


Case Two: Single File
+++++++++++++++++++++++

Alternatively, if you only have one single count matrix for analysis, you can go without sample sheets. **Cumulus** currently supports the following formats:

* 10x genomics v2/v3 formats (hdf5 or mtx);
* HCA DCP mtx and loom formats;
* Drop-seq dge formats.

Simply upload your data to the Google Bucket of your workspace, and specify its URL in ``input_file`` field of Cumulus' global inputs (see below). Notice that for dge and loom files, the ``genome`` field in global inputs is required.

In this case, the **aggregate_matrices** step will be skipped.


---------------------------------

Cumulus steps:
^^^^^^^^^^^^^^

**Cumulus** processes single cell data in the following steps:

#. **aggregate_matrices** (optional). When given a CSV format sample sheet, this step aggregates channel-specific count matrices into one big count matrix. Users can specify which channels they want to analyze and which sample attributes they want to import to the count matrix in this step. Otherwise, if a single count matrix file is given, skip this step.

#. **cluster**. This is the main analysis step. In this step, **Cumulus** performs low quality cell filtration, highly variable gene selection, batch correction, dimension reduction, diffusion map calculation, graph-based clustering and 2D visualization calculation (e.g. t-SNE/UMAP/FLE).

#. **de_analysis**. This step is optional. In this step, **Cumulus** can calculate potential markers for each cluster by performing a variety of differential expression (DE) analysis. The available DE tests include Welch's t test, Fisher's exact test, and Mann-Whitney U test. **Cumulus** can also calculate the area under ROC (AUROC) curve values for putative markers. If ``find_markers_lightgbm`` is on, **Cumulus** will try to identify cluster-specific markers by training a LightGBM classifier. If the samples are human or mouse immune cells, **Cumulus** can also optionally annotate putative cell types for each cluster based on known markers.

#. **plot**. This step is optional. In this step, **Cumulus** can generate 6 types of figures based on the **cluster** step results:

	- **composition** plots which are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects. 
	- **tsne**, **fitsne**, and **net_tsne**: t-SNE like plots based on different algorithms, respectively. Users can specify cell attributes (e.g. cluster labels, conditions) for coloring side-by-side.
	- **umap** and **net_umap**: UMAP like plots based on different algorithms, respectively. Users can specify cell attributes (e.g. cluster labels, conditions) for coloring side-by-side.
	- **fle** and **net_fle**: FLE (Force-directed Layout Embedding) like plots based on different algorithms, respectively. Users can specify cell attributes (e.g. cluster labels, conditions) for coloring side-by-side.
	- **diffmap** plots which are 3D interactive plots showing the diffusion maps. The 3 coordinates are the first 3 PCs of all diffusion components. 
	- If input is CITE-Seq data, there will be **citeseq_fitsne** plots which are FIt-SNE plots based on epitope expression.

#. **organize_results**. Copy analysis results from execution environment to destination location on Google bucket.

In the following sections, we will first introduce global inputs and then introduce the WDL inputs and outputs for each step separately. But please note that you need to set inputs from all steps simultaneously in the Terra WDL.

Note that we will make the required inputs/outputs bold and all other inputs/outputs are optional.

---------------------------------

global inputs
^^^^^^^^^^^^^

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_file**
	  - Input CSV sample sheet describing metadata of each 10x channel, or a single input count matrix file
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_count_matrix.csv"
	  - 
	* - **output_name**
	  - This is the prefix for all output files. It should contain the google bucket url, subdirectory name and output name prefix
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results"
	  - 
	* - genome
	  - A string contains comma-separated genome names. Cumulus will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
	  - "GRCh38"
	  - 
	* - cumulus_version
	  - cumulus version to use. Versions available: 0.10.0.
	  - "0.10.0"
	  - "0.10.0"
	* - docker_registry
	  - Docker registry to use. Options:

	  	- "cumulusprod/" for Docker Hub images; 

	  	- "quay.io/cumulus/" for backup images on Red Hat registry.
	  - "cumulusprod/"
	  - "cumulusprod/"
	* - zones
	  - Google cloud zones to consider for execution.
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-east1-d us-west1-a us-west1-b"
	* - num_cpu
	  - Number of CPUs per Cumulus job
	  - 32
	  - 64
	* - memory
	  - Memory size string
	  - "200G"
	  - "200G"
	* - disk_space
	  - Total disk space in GB
	  - 100
	  - 100
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

---------------------------------

aggregate_matrices
^^^^^^^^^^^^^^^^^^

aggregate_matrices inputs
+++++++++++++++++++++++++

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - restrictions
	  - Select channels that satisfy all restrictions. Each restriction takes the format of name:value,...,value. Multiple restrictions are separated by ';'
	  - "Source:bone_marrow;Platform:NextSeq"
	  - 
	* - attributes
	  - Specify a comma-separated list of outputted attributes. These attributes should be column names in the count_matrix.csv file
	  - "Source,Platform,Donor"
	  - 
	* - minimum_number_of_genes
	  - Only keep barcodes with at least this number of expressed genes
	  - 100
	  - 100
	* - is_dropseq
	  - If inputs are dropseq data
	  - true
	  - false

aggregate_matrices output
+++++++++++++++++++++++++

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_h5sc
	  - File
	  - Aggregated count matrix in Cumulus hdf5 (h5sc) format

---------------------------------

cluster
^^^^^^^

cluster inputs
++++++++++++++

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - channel
	  - Specify the cell barcode attribute to represent different samples.
	  - "Donor" 
	  - 
	* - black_list
	  - Cell barcode attributes in black list will be poped out. Format is "attr1,attr2,...,attrn".
	  - "attr1,attr2,attr3""
	  - 
	* - min_genes_on_raw
	  - If input are raw 10x matrix, which include all barcodes, perform a pre-filtration step to keep the data size small. In the pre-filtration step, only keep cells with at least <min_genes_on_raw> of genes
	  - 100
	  - 100
	* - cite_seq
	  - | Data are CITE-Seq data. cumulus will perform analyses on RNA count matrix first. 
	    | Then it will attach the ADT matrix to the RNA matrix with all antibody names changing to 'AD-' + antibody_name. 
	    | Lastly, it will embed the antibody expression using FIt-SNE (the basis used for plotting is 'citeseq_fitsne')
	  - false
	  - false
	* - cite_seq_capping
	  - For CITE-Seq surface protein expression, make all cells with expression > <percentile> to the value at <percentile> to smooth outlier. Set <percentile> to 100.0 to turn this option off.
	  - 10.0
	  - 99.99
	* - select_only_singlets
	  - If we have demultiplexed data, turning on this option will make cumulus only include barcodes that are predicted as singlets
	  - false
	  - false
	* - output_filtration_results
	  - If write cell and gene filtration results to a spreadsheet
	  - true
	  - true
	* - plot_filtration_results
	  - If plot filtration results as PDF files
	  - true
	  - true
	* - plot_filtration_figsize
	  - Figure size for filtration plots. <figsize> is a comma-separated list of two numbers, the width and height of the figure (e.g. 6,4)
	  - 6,4
	  -
	* - output_seurat_compatible
	  - Generate Seurat-compatible h5ad file. Caution: File size might be large, do not turn this option on for large data sets.
	  - false
	  - false
	* - output_loom
	  - If generate loom-formatted file
	  - false
	  - false
	* - output_parquet
	  - If generate parquet-formatted file
	  - false
	  - false
	* - min_genes
	  - Only keep cells with at least <min_genes> of genes
	  - 500
	  - 500
	* - max_genes
	  - Only keep cells with less than <max_genes> of genes
	  - 6000
	  - 6000
	* - min_umis
	  - Only keep cells with at least <min_umis> of UMIs
	  - 100
	  - 100
	* - max_umis
	  - Only keep cells with less than <max_umis> of UMIs
	  - 600000
	  - 600000
	* - mito_prefix
	  - Prefix of mitochondrial gene names. This is to identify mitochondrial genes.
	  - "mt-"
	  - "MT-"
	* - percent_mito
	  - Only keep cells with mitochondrial ratio less than <percent_mito>% of total counts
	  - 50
	  - 10.0
	* - gene_percent_cells
	  - Only use genes that are expressed in at <gene_percent_cells>% of cells to select variable genes
	  - 50
	  - 0.05
	* - counts_per_cell_after
	  - Total counts per cell after normalization, before transforming the count matrix into Log space.
	  - 1e5
	  - 1e5
	* - select_hvf_flavor
	  - Highly variable feature selection method. Options: 

	  	- "pegasus": New selection method proposed in Pegasus, the analysis module of Cumulus workflow.
	  	- "Seurat": Conventional selection method used by Seurat and SCANPY.
	  - "pegasus"
	  - "pegasus"
	* - select_hvf_ngenes
	  - Select top <select_hvf_ngenes> highly variable features. If <select_hvf_flavor> is "Seurat" and <select_hvf_ngenes> is "None", select HVGs with z-score cutoff at 0.5.
	  - 2000
	  - 2000
	* - no_select_hvf
	  - Do not select highly variable features.
	  - false
	  - false
	* - correct_batch_effect
	  - If correct batch effects
	  - false
	  - false
	* - batch_group_by
	  - | Batch correction assumes the differences in gene expression between channels are due to batch effects. 
	    | However, in many cases, we know that channels can be partitioned into several groups and each group is biologically different from others. 
	    | In this case, we will only perform batch correction for channels within each group. This option defines the groups. 
	    | If <expression> is None, we assume all channels are from one group. Otherwise, groups are defined according to <expression>.
	    | <expression> takes the form of either ‘attr’, or ‘attr1+attr2+…+attrn’, or ‘attr=value11,…,value1n_1;value21,…,value2n_2;…;valuem1,…,valuemn_m’.
	    | In the first form, ‘attr’ should be an existing sample attribute, and groups are defined by ‘attr’.
	    | In the second form, ‘attr1’,…,’attrn’ are n existing sample attributes and groups are defined by the Cartesian product of these n attributes.
	    | In the last form, there will be m + 1 groups. 
	    | A cell belongs to group i (i > 0) if and only if its sample attribute ‘attr’ has a value among valuei1,…,valuein_i. 
	    | A cell belongs to group 0 if it does not belong to any other groups
	  - "Donor"
	  - None
	* - random_state
	  - Random number generator seed
	  - 0
	  - 0
	* - nPC
	  - Number of principal components
	  - 50
	  - 50
	* - knn_K
	  - Number of nearest neighbors used for constructing affinity matrix.
	  - 50
	  - 100
	* - knn_full_speed
	  - For the sake of reproducibility, we only run one thread for building kNN indices. Turn on this option will allow multiple threads to be used for index building. However, it will also reduce reproducibility due to the racing between multiple threads.
	  - false
	  - false
	* - run_diffmap
	  - Whether to calculate diffusion map or not. It will be automatically set to ``true`` when input **run_fle** or **run_net_fle** is set.
	  - false
	  - false
	* - diffmap_ndc
	  - Number of diffusion components
	  - 100
	  - 100
	* - diffmap_maxt
	  - Maximum time stamp in diffusion map computation to search for the knee point.
	  - 5000
	  - 5000
	* - run_louvain
	  - Run Louvain clustering algorithm
	  - true
	  - true
	* - louvain_resolution
	  - Resolution parameter for the Louvain clustering algorithm
	  - 1.3
	  - 1.3
	* - louvain_class_label
	  - Louvain cluster label name in analysis result.
	  - "louvain_labels"
	  - "louvain_labels"
	* - run_leiden
	  - Run Leiden clustering algorithm.
	  - false
	  - false
	* - leiden_resolution
	  - Resolution parameter for the Leiden clustering algorithm.
	  - 1.3
	  - 1.3
	* - leiden_niter
	  - Number of iterations of running the Leiden algorithm. If negative, run Leiden iteratively until no improvement.
	  - 2
	  - -1
	* - leiden_class_label
	  - Leiden cluster label name in analysis result.
	  - "leiden_labels"
	  - "leiden_labels"
	* - run_spectral_louvain
	  - Run Spectral Louvain clustering algorithm
	  - false
	  - false
	* - spectral_louvain_basis
	  - Basis used for KMeans clustering. Use diffusion map by default. If diffusion map is not calculated, use PCA coordinates. Users can also specify "pca" to directly use PCA coordinates.
	  - "diffmap"
	  - "diffmap"
	* - spectral_louvain_resolution
	  - Resolution parameter for louvain.
	  - 1.3
	  - 1.3
	* - spectral_louvain_class_label
	  - Spectral louvain label name in analysis result.
	  - "spectral_louvain_labels"
	  - "spectral_louvain_labels"
	* - run_spectral_leiden
	  - Run Spectral Leiden clustering algorithm.
	  - false
	  - false
	* - spectral_leiden_basis
	  - Basis used for KMeans clustering. Use diffusion map by default. If diffusion map is not calculated, use PCA coordinates. Users can also specify "pca" to directly use PCA coordinates.
	  - "diffmap"
	  - "diffmap"
	* - spectral_leiden_resolution
	  - Resolution parameter for leiden.
	  - 1.3
	  - 1.3
	* - spectral_leiden_class_label
	  - Spectral leiden label name in analysis result.
	  - "spectral_leiden_labels"
	  - "spectral_leiden_labels"
	* - run_tsne
	  - Run multi-core t-SNE for visualization
	  - false
	  - false
	* - tsne_perplexity
	  - t-SNE’s perplexity parameter, also used by FIt-SNE.
	  - 30
	  - 30
	* - run_fitsne
	  - Run FIt-SNE for visualization
	  - true
	  - true
	* - run_umap
	  - Run UMAP for visualization
	  - false
	  - false
	* - umap_K
	  - K neighbors for UMAP.
	  - 15
	  - 15
	* - umap_min_dist
	  - UMAP parameter.
	  - 0.5
	  - 0.5
	* - umap_spread
	  - UMAP parameter.
	  - 1.0
	  - 1.0
	* - run_fle
	  - Run force-directed layout embedding (FLE) for visualization
	  - false
	  - false
	* - fle_K
	  - Number of neighbors for building graph for FLE
	  - 50
	  - 50
	* - fle_target_change_per_node
	  - Target change per node to stop FLE.
	  - 2.0
	  - 2.0
	* - fle_target_steps
	  - Maximum number of iterations before stopping the algoritm
	  - 5000
	  - 5000
	* - net_down_sample_fraction
	  - Down sampling fraction for net-related visualization
	  - 0.1
	  - 0.1
	* - run_net_tsne
	  - Run Net tSNE for visualization
	  - false
	  - false
	* - net_tsne_out_basis
	  - Basis name for Net t-SNE coordinates in analysis result
	  - "net_tsne"
	  - "net_tsne"
	* - run_net_umap
	  - Run Net UMAP for visualization
	  - false
	  - false
	* - net_umap_out_basis
	  - Basis name for Net UMAP coordinates in analysis result
	  - "net_umap"
	  - "net_umap"
	* - run_net_fle
	  - Run Net FLE for visualization
	  - false
	  - false
	* - net_fle_out_basis
	  - Basis name for Net FLE coordinates in analysis result.
	  - "net_fle"
	  - "net_fle"

cluster outputs
+++++++++++++++

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - **output_h5ad**
	  - File
	  - | Output file in h5ad format (output_name.h5ad).
	    | To load this file in Python, you need to first install `Pegasus <https://pegasus.readthedocs.io/en/latest/installation.html>`_ on your local machine. Then use ``import pegasus as pg; data = pg.read_input('output_name.h5ad')`` in Python interpreter.
	    | The log-normalized expression matrix is stored in ``data.X`` as a CSR-format sparse matrix.
	    | The ``obs`` field contains cell related attributes, including clustering results.
	    | For example, ``data.obs_names`` records cell barcodes; ``data.obs['Channel']`` records the channel each cell comes from;
	    | ``data.obs['n_genes']``, ``data.obs['n_counts']``, and ``data.obs['percent_mito']`` record the number of expressed genes, total UMI count, and mitochondrial rate for each cell respectively;
	    | ``data.obs['louvain_labels']``, ``data.obs['leiden_labels']``, ``data.obs['spectral_louvain_labels']``, and ``data.obs['spectral_leiden_labels']`` record each cell's cluster labels using different clustring algorithms;
	    | The ``var`` field contains gene related attributes.
	    | For example, ``data.var_names`` records gene symbols, ``data.var['gene_ids']`` records Ensembl gene IDs, and ``data.var['highly_variable_features']`` records selected variable genes.
	    | The ``obsm`` field records embedding coordinates.
	    | For example, ``data.obsm['X_pca']`` records PCA coordinates, ``data.obsm['X_tsne']`` records t-SNE coordinates,
	    | ``data.obsm['X_umap']`` records UMAP coordinates, ``data.obsm['X_diffmap']`` records diffusion map coordinates,
	    | ``data.obsm['X_diffmap_pca']`` records the first 3 PCs by projecting the diffusion components using PCA,
	    | and ``data.obsm['X_fle']`` records the force-directed layout coordinates from the diffusion components.
	    | The ``varm`` field records DE analysis results if performed: ``data.varm['de_res']``.
	    | The ``uns`` field stores other related information, such as reference genome (``data.uns['genome']``), kNN on PCA coordinates (``data.uns['pca_knn_indices']`` and ``data.uns['pca_knn_distances']``), etc.
	* - output_seurat_h5ad
	  - File
	  - h5ad file in seurat-compatible manner. This file can be loaded into R and converted into a Seurat object (see `here <./cumulus.html#load-cumulus-results-into-seurat>`_ for instructions)
	* - output_filt_xlsx
	  - File
	  - | Spreadsheet containing filtration results (output_name.filt.xlsx).
	    | This file has two sheets --- Cell filtration stats and Gene filtration stats.
	    | The first sheet records cell filtering results and it has 10 columns:
	    | Channel, channel name; kept, number of cells kept; median_n_genes, median number of expressed genes in kept cells; median_n_umis, median number of UMIs in kept cells;
	    | median_percent_mito, median mitochondrial rate as UMIs between mitochondrial genes and all genes in kept cells;
	    | filt, number of cells filtered out; total, total number of cells before filtration, if the input contain all barcodes, this number is the cells left after 'min_genes_on_raw' filtration;
	    | median_n_genes_before, median expressed genes per cell before filtration; median_n_umis_before, median UMIs per cell before filtration;
	    | median_percent_mito_before, median mitochondrial rate per cell before filtration.
	    | The channels are sorted in ascending order with respect to the number of kept cells per channel.
	    | The second sheet records genes that failed to pass the filtering.
	    | This sheet has 3 columns: gene, gene name; n_cells, number of cells this gene is expressed; percent_cells, the fraction of cells this gene is expressed.
	    | Genes are ranked in ascending order according to number of cells the gene is expressed.
	    | Note that only genes not expressed in any cell are removed from the data.
	    | Other filtered genes are marked as non-robust and not used for TPM-like normalization
	* - output_filt_plot
	  - Array[File]
	  - | If not empty, this array contains 3 PDF files.
	    | output_name.filt.gene.pdf, which contains violin plots contrasting gene count distributions before and after filtration per channel.
	    | output_name.filt.UMI.pdf, which contains violin plots contrasting UMI count distributions before and after filtration per channel.
	    | output_name.filt.mito.pdf, which contains violin plots contrasting mitochondrial rate distributions before and after filtration per channel
	* - output_loom_file
	  - File
	  - Outputted loom file (output_name.loom)
	* - output_parquet_file
	  - File
	  - Outputted PARQUET file that contains metadata and expression levels for every gene (output_name.parquet)

---------------------------------

de_analysis
^^^^^^^^^^^

de_analysis inputs
++++++++++++++++++

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - perform_de_analysis
	  - If perform differential expression (DE) analysis
	  - true
	  - true
	* - cluster_labels
	  - Specify the cluster label used for DE analysis
	  - "louvain_labels"
	  - "louvain_labels" 
	* - alpha
	  - Control false discovery rate at <alpha>
	  - 0.05
	  - 0.05
	* - auc
	  - Calculate area under ROC (AUROC)
	  - true
	  - true
	* - fisher
	  - Calculate Fisher’s exact test
	  - true
	  - true
	* - t_test
	  - Calculate Welch's t-test.
	  - true
	  - true
	* - mwu
	  - Calculate Mann-Whitney U test
	  - false
	  - false
	* - find_markers_lightgbm
	  - If also detect markers using LightGBM
	  - false
	  - false
	* - remove_ribo
	  - Remove ribosomal genes with either RPL or RPS as prefixes. Currently only works for human data
	  - false
	  - false
	* - min_gain
	  - Only report genes with a feature importance score (in gain) of at least <gain>
	  - 1.0
	  - 1.0 
	* - annotate_cluster
	  - If also annotate cell types for clusters based on DE results
	  - false
	  - false
	* - annotate_de_test
	  - Differential Expression test to use for inference on cell types. Options: "t", "fisher", or "mwu"
	  - "t"
	  - "t"
	* - organism
	  - Organism, could either be "human_immune", "mouse_immune", "human_brain", "mouse_brain" or a Google bucket link to a JSON file describing the markers
	  - "mouse_brain"
	  - "human_immune"
	* - minimum_report_score
	  - Minimum cell type score to report a potential cell type
	  - 0.5
	  - 0.5

de_analysis outputs
+++++++++++++++++++

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_de_h5ad
	  - File
	  - h5ad-formatted results with DE results updated (output_name.h5ad)
	* - output_de_xlsx
	  - File
	  - Spreadsheet reporting DE results (output_name.de.xlsx)
	* - output_markers_xlsx
	  - File
	  - An excel spreadsheet containing detected markers. Each cluster has one tab in the spreadsheet and each tab has three columns, listing markers that are strongly up-regulated, weakly up-regulated and down-regulated (output_name.markers.xlsx)
	* - output_anno_file
	  - File
	  - Annotation file (output_name.anno.txt)

How cell type annotation works
++++++++++++++++++++++++++++++

In this subsection, we will describe the format of input JSON cell type marker file, the *ad hoc* cell type inference algorithm, and the format of the output putative cell type file.

JSON file
*********

The top level of the JSON file is an object with two name/value pairs:

	- **title**: A string to describe what this JSON file is for (e.g. "Mouse brain cell markers").
	- **cell_types**: List of all cell types this JSON file defines. In this list, each cell type is described using a separate object with 2 to 3 name/value pairs:

		- **name**: Cell type name (e.g. "GABAergic neuron").
		- **markers**: List of gene-marker describing objects, each of which has 2 name/value pairs:

			- **genes**: List of positive and negative gene markers (e.g. ``["Rbfox3+", "Flt1-"]``).
			- **weight**: A real number between ``0.0`` and ``1.0`` to describe how much we trust the markers in **genes**. 
		
		All markers in **genes** share the weight evenly. For instance, if we have 4 markers and the weight is 0.1, each marker has a weight of ``0.1 / 4 = 0.025``. 

		The weights from all gene-marker describing objects of the same cell type should sum up to 1.0.

		- **subtypes**: Description on cell subtypes for the cell type. It has the same structure as the top level JSON object.

See below for an example JSON snippet::

	{
	  "title" : "Mouse brain cell markers",
	    "cell_types" : [
	      {
	        "name" : "Glutamatergic neuron",
	        "markers" : [
	          {
	            "genes" : ["Rbfox3+", "Reln+", "Slc17a6+", "Slc17a7+"],
	            "weight" : 1.0
	          }
	        ],
	        "subtypes" : {
	          "title" : "Glutamatergic neuron subtype markers",
	            "cell_types" : [
	              {
	                "name" : "Glutamatergic layer 4",
	                "markers" : [
	                  {
	                    "genes" : ["Rorb+", "Paqr8+"],
	                    "weight" : 1.0
	                  }
	                ]
	              }
	            ]
	        }
	      }
	    ]
	}

Inference Algorithm
*********************

We have already calculated the up-regulated and down-regulated genes for each cluster in the differential expression analysis step.

First, load gene markers for each cell type from the JSON file specified, and exclude marker genes, along with their associated weights, that are not expressed in the data. 

Then scan each cluster to determine its putative cell types. For each cluster and putative cell type, we calculate a score between ``0`` and ``1``, which describes how likely cells from the cluster are of this cell type. The higher the score is, the more likely cells are from the cell type. 

To calculate the score, each marker is initialized with a maximum impact value (which is ``2``). Then do case analysis as follows:

	- For a positive marker:

		- If it is not up-regulated, its impact value is set to ``0``. 

		- Otherwise, if it is up-regulated:

			- If it additionally has a fold change in percentage of cells expressing this marker (within cluster vs. out of cluster) no less than ``1.5``, it has an impact value of ``2`` and is recorded as a **strong supporting marker**. 

			- If its fold change (``fc``) is less than ``1.5``, this marker has an impact value of ``1 + (fc - 1) / 0.5`` and is recorded as a **weak supporting marker**. 

	- For a negative marker: 

		- If it is up-regulated, its impact value is set to ``0``. 

		- If it is neither up-regulated nor down-regulated, its impact value is set to ``1``. 

		- Otherwise, if it is down-regulated: 

			- If it additionally has ``1 / fc`` (where ``fc`` is its fold change) no less than ``1.5``, it has an impact value of ``2`` and is recorded as a **strong supporting marker**. 

			- If ``1 / fc`` is less than ``1.5``, it has an impact value of ``1 + (1 / fc - 1) / 0.5`` and is recorded as a **weak supporting marker**. 

The score is calculated as the weighted sum of impact values weighted over the sum of weights multiplied by 2 from all expressed markers. If the score is larger than 0.5 and the cell type has cell subtypes, each cell subtype will also be evaluated. 

Output annotation file
**********************

For each cluster, putative cell types with scores larger than ``minimum_report_score`` will be reported in descending order with respect to their scores. The report of each putative cell type contains the following fields:

	- **name**: Cell type name.
	- **score**: Score of cell type.
	- **average marker percentage**: Average percentage of cells expressing marker within the cluster between all positive supporting markers.
	- **strong support**: List of strong supporting markers. Each marker is represented by a tuple of its name and percentage of cells expressing it within the cluster.
	- **weak support**: List of week supporting markers. It has the same structure as **strong support**.

---------------------------------

plot
^^^^

The h5ad file contains a default cell attribute ``Channel``, which records which channel each that single cell comes from. If the input is a CSV format sample sheet, ``Channel`` attribute matches the ``Sample`` column in the sample sheet. Otherwise, it's specified in ``channel`` field of the cluster inputs. 

Other cell attributes used in plot must be added via ``attributes`` field in the ``aggregate_matrices`` inputs.


plot inputs
+++++++++++

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - plot_composition
	  - | Takes the format of "label:attr,label:attr,...,label:attr".
	    | If non-empty, generate composition plot for each "label:attr" pair. 
	    | "label" refers to cluster labels and "attr" refers to sample conditions
	  - "louvain_labels:Donor"
	  - None
	* - plot_fitsne
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored FIt-SNEs side by side
	  - "louvain_labels,Donor"
	  - None
	* - plot_tsne
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored t-SNEs side by side
	  - "louvain_labels,Channel"
	  - None
	* - plot_umap
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored UMAP side by side
	  - "louvain_labels,Donor"
	  - None
	* - plot_fle
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored FLE (force-directed layout embedding) side by side
	  - "louvain_labels,Donor"
	  - None
	* - plot_diffmap
	  - | Takes the format of "attr,attr,...,attr".
	    | If non-empty, generate attr colored 3D interactive plot. 
	    | The 3 coordinates are the first 3 PCs of all diffusion components
	  - "louvain_labels,Donor"
	  - None
	* - plot_citeseq_fitsne
	  - | plot cells based on FIt-SNE coordinates estimated from antibody expressions.
	    | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored FIt-SNEs side by side
	  - "louvain_labels,Donor"
	  - None
	* - plot_net_tsne
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored t-SNEs side by side based on net t-SNE result.
	  - "leiden_labels,Channel"
	  - None
	* - plot_net_umap
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored UMAP side by side based on net UMAP result.
	  - "leiden_labels,Donor"
	  - None
	* - plot_net_fle
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored FLE (force-directed layout embedding) side by side
	    | based on net FLE result.
	  - "leiden_labels,Donor"
	  - None

plot outputs
++++++++++++

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_pdfs
	  - Array[File]
	  - Outputted pdf files
	* - output_htmls
	  - Array[File]
	  - Outputted html files


---------------------------------

Generate SCP Output
^^^^^^^^^^^^^^^^^^^

Generate analysis result in `Single Cell Portal`_ (SCP) compatible format.

scp_output inputs
+++++++++++++++++


.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - generate_scp_outputs
	  - Whether to generate SCP format output or not.
	  - false
	  - false
	* - output_dense
	  - Output dense expression matrix, instead of the default sparse matrix format.
	  - false
	  - false


scp_output outputs
++++++++++++++++++

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_scp_files
	  - Array[File]
	  - Outputted SCP format files.

---------------------------------

Run CITE-Seq analysis
~~~~~~~~~~~~~~~~~~~~~~~~

To run CITE-Seq analysis, turn on ``cite_seq`` option in cluster inputs of cumulus workflow. 

An embedding of epitope expressions via FIt-SNE is available at basis ``X_citeseq_fitsne``. 

To plot this epitope embedding, specify attributes to plot in ``plot_citeseq_fitsne`` field of cluster inputs.

---------------------------------

Run subcluster analysis
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Once we have **cumulus** outputs, we could further analyze a subset of cells by running **cumulus_subcluster**. To run **cumulus_subcluster**, follow the following steps:

#. Import **cumulus_subcluster** method.

	See the Terra documentation for `adding a workflow`_. The cumulus workflow is under ``Broad Methods Repository`` with name "**cumulus/cumulus_subcluster**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export cumulus workflow in the drop-down menu.

#. In your workspace, open ``cumulus_subcluster`` in ``WORKFLOWS`` tab. Select ``Process single workflow from files`` as below

	.. image:: images/single_workflow.png

   and click the ``SAVE`` button.

cumulus_subcluster steps:
^^^^^^^^^^^^^^^^^^^^^^^^^^

*cumulus_subcluster* processes the subset of single cells in the following steps:

#. **subcluster**. In this step, **cumulus_subcluster** first select the subset of cells from **cumulus** outputs according to user-provided criteria. It then performs batch correction, dimension reduction, diffusion map calculation, graph-based clustering and 2D visualization calculation (e.g. t-SNE/UMAP/FLE).

#. **de_analysis** (optional). In this step, **cumulus_subcluster** calculates potential markers for each cluster by performing a variety of differential expression (DE) analysis. The available DE tests include Welch's t test, Fisher's exact test, and Mann-Whitney U test. **cumulus_subcluster** can also calculate the area under ROC curve (AUROC) values for putative markers. If the samples are human or mouse immune cells, **cumulus_subcluster** can optionally annotate putative cell types for each cluster based on known markers.

#. **plot** (optional). In this step, **cumulus_subcluster** can generate the following 5 types of figures based on the **subcluster** step results:

	- **composition** plots which are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects. 

	- **tsne**, **fitsne**, and **net_tsne**: t-SNE like plots based on different algorithms, respectively. Users can specify different cell attributes (e.g. cluster labels, conditions) for coloring side-by-side.

	- **umap** and **net_umap**: UMAP like plots based on different algorithms, respectively. Users can specify different cell attributes (e.g. cluster labels, conditions) for coloring side-by-side.

	- **fle** and **net_fle**: FLE (Force-directed Layout Embedding) like plots based on different algorithms, respectively. Users can specify different cell attributes (e.g. cluster labels, conditions) for coloring side-by-side.

	- **diffmap** plots which are 3D interactive plots showing the diffusion maps. The 3 coordinates are the first 3 PCs of all diffusion components.

cumulus_subcluster's inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**cumulus_subcluster** shares many inputs/outputs with **cumulus**, we will only cover inputs/outputs that are specific to **cumulus_subcluster** in this section.

Note that we will make the required inputs/outputs bold and all other inputs/outputs are optional.

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_h5ad**
	  - Google bucket URL of input h5ad file containing *cumulus* results
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results.h5ad"
	  - 
	* - **output_name**
	  - This is the prefix for all output files. It should contain the Google bucket URL, subdirectory name and output name prefix
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results_sub"
	  - 
	* - **subset_selections**
	  - | Specify which cells will be included in the subcluster analysis.
	    | This field contains one or more <subset_selection> strings separated by ';'. 
	    | Each <subset_selection> string takes the format of 'attr:value,…,value', which means select cells with attr in the values. 
	    | If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings
	  - "louvain_labels:3,6" or "louvain_labels:3,6;Donor:1,2"
	  - 
	* - calculate_pseudotime
	  - Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes
	  - "sample_1-ACCCGGGTTT-1,sample_1-TCCCGGGAAA-2"
	  - None
	* - num_cpu
	  - Number of cpus per cumulus job
	  - 32
	  - 64
	* - memory
	  - Memory size string
	  - "200G"
	  - "200G"
	* - disk_space
	  - Total disk space in GB
	  - 100
	  - 100
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

.. role:: red-bold

For other **cumulus_subcluster** inputs, please refer to `cumulus cluster inputs list`_ for details. Notice that some inputs (as listed below) in **cumulus cluster** inputs list are :red-bold:`DISABLED` for **cumulus_subcluster**:
	
	- cite_seq
	- cite_seq_capping
	- output_filtration_results
	- plot_filtration_results
	- plot_filtration_figsize
	- output_seurat_compatible
	- batch_group_by
	- min_genes
	- max_genes
	- min_umis
	- max_umis
	- mito_prefix
	- percent_mito
	- gene_percent_cells
	- min_genes_on_raw
	- counts_per_cell_after

.. _cumulus cluster inputs list: ./cumulus.html#cluster


cumulus_subcluster's outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - **output_h5ad**
	  - File
	  - | h5ad-formatted HDF5 file containing all results (output_name.h5ad). 
	    | If ``perform_de_analysis`` is on, this file should be the same as *output_de_h5ad*.
	    | To load this file in Python, it's similar as in `cumulus cluster outputs <./cumulus.html#cluster-outputs>`_ section.
	    | Besides, for subcluster results, there is a new cell attributes in ``data.obs['pseudo_time']``, which records the inferred pseudotime for each cell.
	* - output_loom_file
	  - File
	  - Generated loom file (output_name.loom)
	* - output_parquet_file
	  - File
	  - Generated PARQUET file that contains metadata and expression levels for every gene (output_name.parquet)
	* - output_de_h5ad
	  - File
	  - Generated h5ad-formatted results with DE results updated (output_name.h5ad)
	* - output_de_xlsx
	  - File
	  - Generated Spreadsheet reporting DE results (output_name.de.xlsx)
	* - output_pdfs
	  - Array[File]
	  - Generated pdf files
	* - output_htmls
	  - Array[File]
	  - Generated html files

---------------------------------

Load Cumulus results into Seurat  
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

First, you need to set ``output_seurat_compatible`` field to ``true`` in cumulus cluster inputs to generate a Seurat-compatible output file ``output_name.seurat.h5ad``, in addition to the normal result ``output_name.h5ad``.

Notice that python, the `anndata`_ python library with version at least ``0.6.22.post1``, and the ``reticulate`` R library are required to load the result into Seurat.

Execute the R code below to load the results into Seurat (working with both Seurat v2 and v3)::

	source("https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/workflows/cumulus/h5ad2seurat.R")
	ad <- import("anndata", convert = FALSE)
	test_ad <- ad$read_h5ad("output_name.seurat.h5ad")
	result <- convert_h5ad_to_seurat(test_ad)

The resulting Seurat object ``result`` has three data slots: 

	- **raw.data** records filtered raw count matrix. 
	- **data** records filtered and log-normalized expression matrix. 
	- **scale.data** records variable-gene-selected, standardized expression matrix that are ready to perform PCA.

---------------------------------

Visualize Cumulus results in Python
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Ensure you have `Pegasus <https://pegasus.readthedocs.io/en/latest/installation.html>`_ installed.

Download your analysis result data, say ``output_name.h5ad``, from Google bucket to your local machine.

Load the output::

	import pegasus as pg
	adata = pg.read_input("output_name.h5ad")

Violin plot of the computed quality measures::

	fig = pg.violin(adata, keys = ['n_genes', 'n_counts', 'percent_mito'], by = 'passed_qc')
	fig.savefig('output_file.qc.pdf', dpi = 500)

t-SNE plot colored by louvain cluster labels and channel::

	fig = pg.embedding(adata, basis = 'tsne', keys = ['louvain_labels', 'Channel'])
	fig.savefig('output_file.tsne.pdf', dpi = 500)

t-SNE plot colored by genes of interest::

	fig = pg.embedding(adata, basis = 'tsne', keys = ['CD4', 'CD8A'])
	fig.savefig('output_file.genes.tsne.pdf', dpi = 500)

For other embedding plots using FIt-SNE (``fitsne``), Net t-SNE (``net_tsne``), CITE-Seq FIt-SNE (``citeseq_fitsne``), UMAP (``umap``), Net UMAP (``net_umap``), FLE (``fle``), or Net FLE (``net_fle``) coordinates, simply substitute its basis name for ``tsne`` in the code above.

Composition plot on louvain cluster labels colored by channel::

	fig = pg.composition_plot(adata, by = 'louvain_labels', condition = 'Channel')
	fig.savefig('output_file.composition.pdf', dpi = 500)


.. _anndata: https://anndata.readthedocs.io/en/latest/
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/
.. _Single Cell Portal: https://portals.broadinstitute.org/single_cell
