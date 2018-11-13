Run single-cell cloud-based analysis module (scCloud) for scRNA-Seq data analysis
---------------------------------------------------------------------------------

Follow the steps below to run **scCloud** on FireCloud.

#. Create a sample sheet, **count_matrix.csv**, which describes the metadata for each 10x channel. The sample sheet should at least contain 2 columns --- *Sample* and *Location*. *Sample* refers to sample names and *Location* refers to the location of the channel-specific count matrix in either 10x format (e.g. ``gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5``) or dropseq format (e.g. ``gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_2/sample_2.umi.dge.txt.gz``). You are free to add any other columns and these columns will be used in selecting channels for futher analysis. In the example below, we have *Source*, which refers to the tissue of origin, *Platform*, which refers to the sequencing platform, *Donor*, which refers to the donor ID, and *Reference*, which refers to the reference genome.

	Example::

		Sample,Source,Platform,Donor,Reference,Location
		sample_1,bone_marrow,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5
		sample_2,bone_marrow,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_2/filtered_gene_bc_matrices_h5.h5
		sample_3,pbmc,NextSeq,1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_3/filtered_gene_bc_matrices_h5.h5
		sample_4,pbmc,NextSeq,2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_4/filtered_gene_bc_matrices_h5.h5

	If you ran **cellranger_mkfastq_count**, you should obtain a template **count_matrix.csv** file that you can modify from **cellranger_mkfastq_count**'s outputs. 

#. Upload your sample sheet to the workspace.  

	Example::
	
		gsutil cp /foo/bar/projects/my_count_matrix.csv gs://fc-e0000000-0000-0000-0000-000000000000/

#. Import **scCloud** method.

	In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Import From Method Repository". Type **scCloud/scCloud**.

#. Uncheck "Configure inputs/outputs using the Workspace Data Model"

---------------------------------

scCloud steps:
^^^^^^^^^^^^^^^

**scCloud** processes single cell data in the following steps:

#. **aggregate_matrix**. This step aggregates channel-specific count matrices into one big count matrix. Users could specify which channels they want to analyze and which sample attributes they want to import to the count matrix in this step.

#. **cluster**. This step is the main analysis step. In this step, **scCloud** performs low quality cell filtration, variable gene selection, batch correction, dimension reduction, diffusion map calculation, graph-based clustering and 2D visualization calculation (e.g. t-SNE/FLE).

#. **de_analysis**. This step is optional. In this step, **scCloud** could calculate potential markers for each cluster by performing a variety of differential expression (DE) analysis. The available DE tests include Welch's t test, Fisher's exact test, and Mann-Whitney U test. **scCloud** could also calculate the area under ROC curve values for putative markers. If the samples are human or mouse immune cells, **scCloud** could also optionally annotate putative cell types for each cluster based on known markers.

#. **plot**. This step is optional. In this step, **scCloud** could generate 3 types of figures based on the **cluster** step results. First, **composition** plots are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects. Second, **tsne** plot shows the same t-SNE colored by different attributes (e.g. cluster labels, conditions) side-by-side. Lastly, **diffmap** plots are 3D interactive plots showing the diffusion maps. The 3 coordinates are the first 3 PCs of all diffusion components.

In the following, we will first introduce global inputs and then introduce the WDL inputs and outputs for each step separately. But please note that you need to set inputs from all steps simultaneously in the FireCloud WDL. 

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
	* - **input_count_matrix_csv**
	  - Input CSV file describing metadata of each 10x channel
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_count_matrix.csv"
	  - 
	* - **output_name**
	  - This is the prefix for all output files. It should contain the google bucket url, subdirectory name and output name prefix
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results"
	  - 
	* - genome
	  - A string contains comma-separated genome names. scCloud will read all groups associated with genome names in the list from the hdf5 file. If genome is None, all groups will be considered.
	  - "GRCh38"
	  - 
	* - num_cpu
	  - Number of cpus per scCloud job
	  - 32
	  - 64
	* - memory
	  - Memory size in GB
	  - 200
	  - 200
	* - disk_space
	  - Total disk space
	  - 100
	  - 100
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

---------------------------------

aggregate_matrix
^^^^^^^^^^^^^^^^

aggregate_matrix inputs
+++++++++++++++++++++++

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
	* - select_only_singlets
	  - If we have demultiplexed data, turning on this option will make scCloud only include barcodes that are predicted as singlets
	  - true
	  - false
	* - minimum_number_of_genes
	  - Only keep barcodes with at least this number of expressed genes
	  - 100
	  - 100
	* - is_dropseq
	  - If inputs are dropseq data
	  - true
	  - false

aggregate_matrix output
+++++++++++++++++++++++

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - **output_10x_h5**
	  - File
	  - Aggregated count matrix in 10x format

---------------------------------

cluster
^^^^^^^

cluster inputs
++++++++++++++

Note that we will only list important inputs here. For other inputs, please refer to **scCloud** package documentation.

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - cite_seq
	  - | Data are CITE-Seq data. scCloud will perform analyses on RNA count matrix first. 
	    | Then it will attach the ADT matrix to the RNA matrix with all antibody names changing to 'AD-' + antibody_name. 
	    | Lastly, it will embed the antibody expression using t-SNE (the basis used for plotting is 'citeseq_tsne').
	  - true
	  - false
	* - output_filtration_results
	  - If output cell and gene filtration results to a spreadsheet
	  - true
	  - true
	* - output_seurat_compatible
	  - If output Seurat-compatible h5ad file
	  - true
	  - false
	* - output_loom
	  - If output loom-formatted file
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
	* - min_genes
	  - Only keep cells with at least <min_genes> of genes
	  - 500
	  - 500
	* - max_genes
	  - Only keep cells with less than <number> of genes
	  - 6000
	  - 6000
	* - mito_prefix
	  - Prefix for mitochondrial genes
	  - "mt-"
	  - "MT-"
	* - percent_mito
	  - Only keep cells with mitochondrial ratio less than <percent_mito>
	  - 0.1
	  - 0.1
	* - gene_percent_cells
	  - Only use genes that are expressed in at <ratio> * 100 percent of cells to select variable genes
	  - 0.0005
	  - 0.0005
	* - counts_per_cell_after
	  - Total counts per cell after normalization
	  - 1e5
	  - 1e5
	* - random_state
	  - Random number generator seed
	  - 0
	  - 0
	* - nPC
	  - Number of principal components
	  - 50
	  - 50
	* - nDC
	  - Number of diffusion components
	  - 50
	  - 50
	* - diffmap_K
	  - Number of neighbors used for constructing affinity matrix
	  - 100
	  - 100
	* - diffmap_alpha
	  - Power parameter for diffusion-based pseudotime
	  - 0.5
	  - 0.5
	* - run_louvain
	  - Run louvain clustering algorithm
	  - true
	  - true
	* - louvain_resolution
	  - Resolution parameter for the louvain clustering algorithm
	  - 1.3
	  - 1.3
	* - run_approximated_louvain
	  - Run approximated louvain clustering algorithm
	  - true
	  - false
	* - approx_louvain_ninit
	  - Number of Kmeans tries
	  - 30
	  - 20
	* - approx_louvain_nclusters
	  - Number of clusters for Kmeans initialization
	  - 40
	  - 30
	* - approx_louvain_resolution
	  - Resolution parameter for louvain
	  - 1.3
	  - 1.3
	* - run_tsne
	  - Run multi-core t-SNE for visualization
	  - true
	  - true
	* - tsne_perplexity
	  - t-SNE’s perplexity parameter
	  - 30
	  - 30
	* - run_fitsne
	  - Run FIt-SNE for visualization
	  - true
	  - false
	* - run_umap
	  - Run umap for visualization
	  - true
	  - false
	* - umap_on_diffmap
	  - Run umap on diffusion components
	  - ture
	  - false
	* - run_fle
	  - Run force-directed layout embedding
	  - true
	  - false
	* - fle_K
	  - K neighbors for building graph for FLE
	  - 50
	  - 50
	* - fle_n_steps
	  - Number of iterations for FLE
	  - 10000
	  - 10000

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
	  - h5ad-formatted HDF5 file containing all results (output_name.h5ad)
	* - output_filt_xlsx
	  - File
	  - Spreadsheet containing filtration results (output_name.filt.xlsx)
	* - output_seurat_h5ad
	  - File
	  - Seurat readable h5ad file (output_name.seurat.h5ad)
	* - output_loom_file
	  - File
	  - Outputted loom file (output_name.loom)

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
	  - If perform de analysis
	  - true
	  - true
	* - cluster_labels
	  - Specify the cluster labels used for differential expression analysis
	  - "louvain_labels"
	  - "louvain_labels" 
	* - alpha
	  - Control false discovery rate at <alpha>
	  - 0.05
	  - 0.05
	* - fisher
	  - Calculate Fisher’s exact test
	  - true
	  - true
	* - mwu
	  - Calculate Mann-Whitney U test
	  - true
	  - false
	* - roc
	  - Calculate area under cuver in ROC curve
	  - true
	  - false
	* - annotate_cluster
	  - If also annotate cell types for clusters based on DE results
	  - true
	  - false
	* - organism
	  - Organism, could either be "human_immune", "mouse_immune", or "mouse_brain"
	  - "mouse"
	  - "human"
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
	* - output_anno_file
	  - File
	  - Annotation file (output_name.anno.txt)

plot
^^^^

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
	* - plot_tsne
	  - | Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored t-SNEs side by side
	  - "louvain_labels,Donor"
	  - None
	* - plot_diffmap
	  - | Takes the format of "attr,attr,...,attr".
	    | If non-empty, generate attr colored 3D interactive plot. 
	    | The 3 coordinates are the first 3 PCs of all diffusion components
	  - "louvain_labels,Donor"
	  - None
	* - plot_citeseq_tsne
	  - | plot cells based on t-SNE coordinates estimated from antibody expressions.
		| Takes the format of "attr,attr,...,attr". 
	    | If non-empty, plot attr colored t-SNEs side by side
	  - "louvain_labels"
	  - None

plot outputs
++++++++++++

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_pngs
	  - Array[File]
	  - Outputted png files
	* - output_htmls
	  - Array[File]
	  - Outputted html files

---------------------------------

Run CITE-Seq analysis
---------------------

To run CITE-Seq analysis, turn on ``cite_seq`` option. 

An embedding of epitope expressions via t-SNE is available at basis ``X_citeseq_tsne``. 

To plot this epitope embedding, turn on ``plot_citeseq_tsne`` option.

---------------------------------

Run subcluster analysis
-----------------------

Once we have **scCloud** outputs, we could further analyze a subset of cells by running **scCloud_subcluster**. To run **scCloud_subcluster**, follow the following steps:

#. Import **scCloud_subcluster** method.

	In FireCloud, select the "Method Configurations" tab then click "Import Configuration". Click "Import From Method Repository". Type **scCloud/scCloud_subcluster**.

#. Uncheck "Configure inputs/outputs using the Workspace Data Model".

scCloud_subcluster steps:
^^^^^^^^^^^^^^^^^^^^^^^^^^

*scCloud_subcluster* processes the subset of single cells in the following steps:

#. **subcluster**. In this step, **scCloud_subcluster** first select the subset of cells from **scCloud** outputs according to user-provided criteria. It then performs batch correction, dimension reduction, diffusion map calculation, graph-based clustering and 2D visualization calculation (e.g. t-SNE/FLE).

#. **de_analysis**. This step is optional. In this step, **scCloud_subcluster** could calculate potential markers for each cluster by performing a variety of differential expression (DE) analysis. The available DE tests include Welch's t test, Fisher's exact test, and Mann-Whitney U test. **scCloud_subcluster** could also calculate the area under ROC curve values for putative markers. If the samples are human or mouse immune cells, **scCloud_subcluster** could also optionally annotate putative cell types for each cluster based on known markers.

#. **plot**. This step is optional. In this step, **scCloud_subcluster** could generate 3 types of figures based on the **subcluster** step results. First, **composition** plots are bar plots showing the cell compositions (from different conditions) for each cluster. This type of plots is useful to fast assess library quality and batch effects. Second, **tsne** plot shows the same t-SNE colored by different attributes (e.g. cluster labels, conditions) side-by-side. Lastly, **diffmap** plots are 3D interactive plots showing the diffusion maps. The 3 coordinates are the first 3 PCs of all diffusion components.

scCloud_subcluster's inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Since **scCloud_subcluster** shares many inputs/outputs with **scCloud**, we will only cover inputs/outputs that are specific to **scCloud_subcluster**.

Note that we will make the required inputs/outputs bold and all other inputs/outputs are optional.

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_h5ad**
	  - Input h5ad file containing *scCloud* results
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results.h5ad"
	  - 
	* - **output_name**
	  - This is the prefix for all output files. It should contain the google bucket url, subdirectory name and output name prefix
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir/my_results_sub"
	  - 
	* - **subset_selections**
	  - | Specify which cells will be included in the subcluster analysis.
	    | This field contains one or more <subset_selection> strings separated by ';'. 
	    | Each <subset_selection> string takes the format of 'attr:value,…,value', which means select cells with attr in the values. 
	    | If multiple <subset_selection> strings are specified, the subset of cells selected is the intersection of these strings
	  - "louvain_labels:3,6"
	  - 
	* - calculate_pseudotime
	  - Calculate diffusion-based pseudotimes based on <roots>. <roots> should be a comma-separated list of cell barcodes
	  - "sample_1-ACCCGGGTTT-1"
	  - None
	* - num_cpu
	  - Number of cpus per scCloud job
	  - 32
	  - 64
	* - memory
	  - Memory size in GB
	  - 200
	  - 200
	* - disk_space
	  - Total disk space
	  - 100
	  - 100
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

scCloud_subcluster's outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
	:widths: 5 5 20
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - **output_h5ad**
	  - File
	  - h5ad-formatted HDF5 file containing all results (output_name.h5ad)
	* - output_seurat_h5ad
	  - File
	  - Seurat readable h5ad file (output_name.seurat.h5ad)
	* - output_loom_file
	  - File
	  - Outputted loom file (output_name.loom)
	* - output_de_h5ad
	  - File
	  - h5ad-formatted results with DE results updated (output_name.h5ad)
	* - output_de_xlsx
	  - File
	  - Spreadsheet reporting DE results (output_name.de.xlsx)
	* - output_pngs
	  - Array[File]
	  - Outputted png files
	* - output_htmls
	  - Array[File]
	  - Outputted html files

---------------------------------


Load ``scCloud`` results into ``Seurat``  
-----------------------------------------

First, you need to set ``output_seurat_compatible`` to ``true`` in ``scCloud`` or ``scCloud_subcluster`` to obtain Seurat-compatible h5ad file ``output_name.seurat.h5ad``.

Then following the codes below to load the results into ``Seurat``::

	library(Seurat)
	library(reticulate)
	ad <- import("anndata", convert = FALSE)
	test_ad <- ad$read_h5ad("output_name.seurat.h5ad")
	test <- Convert(test_ad, to = "seurat")
