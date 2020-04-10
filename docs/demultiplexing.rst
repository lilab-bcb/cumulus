Demultiplex cell-hashing/nucleus-hashing/genetic-pooling sc/snRNA-Seq data
--------------------------------------------------------------------------

This ``demultiplexing`` workflow generates gene-count matrices from cell-hashing/nucleus-hashing/genetic-pooling data by demultiplexing.

demuxEM is used for analyzing cell-hashing/nucleus-hashing data, while souporcell and demuxlet are for genetic-pooling data.

Prepare input data and import workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Run ``cellranger_workflow``
++++++++++++++++++++++++++++++++

	To demultiplex, you'll need raw gene count and hashtag matrices for cell-hashing/nucleus-hashing data; raw gene count matrices and genome BAM files for genetic-pooling data. You can generate these data by running the ``cellranger_workflow``.

	Please refer to the `cellranger_workflow tutorial`_ for details.

	When finished, you should be able to find the raw gene count matrix (e.g. ``raw_gene_bc_matrices_h5.h5``), hashtag matrix (e.g. ``sample_1_ADT.csv``) / genome BAM file (e.g. ``possorted_genome_bam.bam``) for each sample.

2. Import ``demultiplexing`` 
++++++++++++++++++++++++++++++

Import *demultiplexing* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *demultiplexing* workflow is under ``Broad Methods Repository`` with name "**cumulus/demultiplexing**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *demultiplexing* workflow in the drop-down menu.

3. Prepare a sample sheet
++++++++++++++++++++++++++++

	**3.1 Sample sheet format:**

	Create a sample sheet, **sample_sheet_demux.csv**, which describes the metadata for each pair of RNA and hashtag data. A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **OUTNAME**
		  - Output name for one pair of RNA and hashtag data. Must be unique per pair.
		* - **RNA**
		  - Google bucket url to the raw gene count matrix generated in Step 1.
		* - **TagFile/ADT**
		  - Google bucket url to the hashtag file generated in Step 1. The column name can be either *TagFile* or *ADT*, where *ADT* is to be backward compatible with sample sheets working with *cumulus/cumulus_hashing_cite_seq* workflow.
		* - **TYPE**
		  - Assay type, which can be ``cell-hashing``, ``nucleus-hashing``, or ``genetic-pooling``.
		* - Genotype
		  - Google bucket url to the reference genotypes in ``vcf.gz`` format. This column is **not required** in the following cases:

		  	- *TYPE* is ``cell-hashing`` or ``nucleus-hashing``;

		  	- *TYPE* is ``genetic-pooling``, *demultiplexing_algorithm* input is ``souporcell``, user wish to run in *de novo* mode without reference genotypes, and don't need to rename cluster names by information from a known genotype vcf file.


	Example::

		OUTNAME,RNA,TagFile,TYPE,Genotype
		sample_1,gs://exp/data_1/raw_gene_bc_matrices_h5.h5,gs://exp/data_1/sample_1_ADT.csv,cell-hashing
		sample_2,gs://exp/data_2/raw_gene_bc_matrices_h5.h5,gs://exp/data_2/sample_2_ADT.csv,nucleus-hashing
		sample_3,gs://exp/data_3/raw_gene_bc_matrices_h5.h5,gs://exp/data_3/possorted_genome_bam.bam,genetic-pooling
		sample_4,gs://exp/data_4/raw_gene_bc_matrices_h5.h5,gs://exp/data_4/possorted_genome_bam.bam,genetic-pooling,gs://exp/variants/ref_genotypes.vcf.gz

	**3.2 Upload your sample sheet to the workspace bucket:**

	Use gsutil_ (you already have it if you've installed Google cloud SDK) in your unix terminal to upload your sample sheet to workspace bucket.

	Example::

			gsutil cp /foo/bar/projects/sample_sheet_demux.tsv gs://fc-e0000000-0000-0000-0000-000000000000/

---------------

Workflow inputs
^^^^^^^^^^^^^^^^

Below are inputs for *demultiplexing* workflow. Notice that required inputs are in bold.

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_sample_sheet**
	  - Input CSV file describing metadata of RNA and hashtag data pairing.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet_demux.csv"
	  -
	* - **output_directory**
	  - This is the output directory (gs url + path) for all results. There will be one folder per RNA-hashtag data pair under this directory.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/demux_output"
	  -
	* - **genome**
	  - Reference genome name. You should choose one from this `genome reference`_ list.
	  - "GRCh38"
	  -
	* - demultiplexing_algorithm
	  - demultiplexing algorithm to use for *genetic-pooling* data. Options:

	  	- "souporcell": Use souporcell_, a reference-genotypes-free algorithm for demultiplexing droplet scRNA-Seq data.

	  	- "demuxlet": Use demuxlet_, a canonical algorithm for demultiplexing droplet scRNA-Seq data.
	  - "souporcell"
	  - "souporcell"
	* - min_num_genes
	  - Only demultiplex cells/nuclei with at least <min_num_genes> expressed genes
	  - 500
	  - 500
	* - docker_registry
	  - Docker registry to use. Notice that docker image for Bustools is seperate.

	  	- "cumulusprod" for Docker Hub images; 

	  	- "quay.io/cumulus" for backup images on Red Hat registry.
	  - "cumulusprod"
	  - "cumulusprod"
	* - zones
	  - Google cloud zones to consider for execution.
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - preemptible
	  - Number of maximum preemptible tries allowed.
	  - 2
	  - 2
	* - demuxEM_alpha_on_samples
	  - demuxEM parameter. The Dirichlet prior concentration parameter (alpha) on samples. An alpha value < 1.0 will make the prior sparse.
	  - 0.0
	  - 0.0
	* - demuxEM_min_num_umis
	  - demuxEM parameter. Only demultiplex cells/nuclei with at least <demuxEM_min_num_umis> of UMIs.
	  - 100
	  - 100
	* - demuxEM_min_signal_hashtag
	  - demuxEM parameter. Any cell/nucleus with less than <demuxEM_min_signal_hashtag> hashtags from the signal will be marked as unknown.
	  - 10.0
	  - 10.0
	* - demuxEM_random_state
	  - demuxEM parameter. The random seed used in the KMeans algorithm to separate empty ADT droplets from others.
	  - 0
	  - 0
	* - demuxEM_generate_diagnostic_plots
	  - demuxEM parameter. If generate a series of diagnostic plots, including the background/signal between HTO counts, estimated background probabilities, HTO distributions of cells and non-cells, etc.
	  - true
	  - true
	* - demuxEM_generate_gender_plot
	  - demuxEM parameter. If generate violin plots using gender-specific genes (e.g. Xist). <demuxEM_generate_gender_plot> is a comma-separated list of gene names
	  - "XIST"
	  -
	* - demuxEM_version
	  - demuxEM version to use. Currently only support "0.1.1".
	  - "0.1.1"
	  - "0.1.1"
	* - demuxEM_num_cpu
	  - demuxEM parameter. Number of CPUs to request for demuxEM per pair.
	  - 8
	  - 8
	* - demuxEM_memory
	  - demuxEM parameter. Memory size (integer) in GB needed for demuxEM per pair.
	  - 10
	  - 10
	* - demuxEM_disk_space
	  - demuxEM parameter. Disk space (integer) in GB needed for demuxEM per pair.
	  - 20
	  - 20
	* - souporcell_version
	  - souporcell version to use. Currently only support "2020.03".
	  - "2020.03"
	  - "2020.03"
	* - souporcell_de_novo_mode
	  - souporcell parameter. If ``true``, run souporcell de novo mode without reference genotypes; otherwise, a reference genotype vcf file specified in sample sheet will be used.
	  - true
	  - true
	* - souporcell_num_clusters
	  - souporcell parameter. Number of expected clusters when doing clustering.
	  - 1
	  - 1
	* - souporcell_rename_donors
	  - | souporcell parameter. A comma-separated list of donor names for renaming clusters achieved by souporcell.
	    | By default, the resulting donors are *Donor1*, *Donor2*, ...
	  - "CB1,CB2,CB3,CB4"
	  - 
	* - souporcell_num_cpu
	  - souporcell parameter. Number of CPUs to request for souporcell per pair.
	  - 32
	  - 32
	* - souporcell_memory
	  - souporcell parameter. Memory size (integer) in GB needed for souporcell per pair.
	  - 120
	  - 120
	* - souporcell_disk_space
	  - souporcell parameter. Disk space (integer) in GB needed for souporcell per pair.
	  - 500
	  - 500
	* - demuxlet_version
	  - demuxlet version to use. Currently only support "0.1b".
	  - "0.1b"
	  - "0.1b"
	* - demuxlet_memory
	  - demuxlet parameter. Memory size (integer) in GB needed for demuxlet per pair.
	  - 10
	  - 10
	* - demuxlet_disk_space
	  - | demuxlet parameter. Disk space size (integer) in GB needed for demuxlet per pair.
	    | Notice that the overall disk space for demuxlet is this disk space plus the size of provided reference genotypes file in the sample sheet.
	  - 2
	  - 2



Workflow outputs
^^^^^^^^^^^^^^^^^^

See the table below for *demultiplexing* workflow outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_folder
	  - String
	  - Google Bucket URL of output directory. Within it, each folder is for one RNA-hashtag data pair in the input sample sheet.

In the output subfolder of each cell-hashing/nuclei-hashing RNA-hashtag data pair, you can find the following files:

.. list-table::
	:widths: 5 10
	:header-rows: 1

	* - Name
	  - Description
	* - output_name.out.demuxEM.zarr
	  - | RNA expression matrix with demultiplexed sample identities in zarr format.
	    | To load this file into Python, you need to first install `Pegasusio`_ on your local machine. Then use ``import pegasusio as io; data = io.read_input("output_name.out.demuxEM.zarr")`` in Python environment.
	    | It contains 2 UnimodalData objects: one with key ``hashing`` is the hashtag count matrix, the other one with genome name key is the demultiplexed RNA count matrix.
	    | To load the hashtag count matrix, type ``hash_data = data.get_data('hashing')``. The count matrix is ``hash_data.X``; cell barcode attributes are stored in ``hash_data.obs``; sample names are in ``hash_data.var_names``. Moreover, the estimated background probability regarding hashtags is in ``hash_data.uns['background_probs']``.
	    | To load the RNA matrix, type ``rna_data = data.get_data('<genome>')``, where ``<genome>`` is the genome name of the data. It only contains cells which have estimated sample assignments. The count matrix is ``rna_data.X``. Cell barcode attributes are stored in ``rna_data.obs``: ``rna_data.obs['demux_type']`` stores the estimated droplet types (singlet/doublet/unknown) of cells; ``rna_data.obs['assignment']`` stores the estimated hashtag(s) that each cell belongs to.
	* - output_name_demux.zarr
	  - Demultiplexed RNA count matrix in zarr format. Please refer to section `load demultiplexing results into Python and R`_ for its structure.
	* - output_name.ambient_hashtag.hist.png
	  - Optional output. A histogram plot depicting hashtag distributions of empty droplets and non-empty droplets.
	* - output_name.background_probabilities.bar.png
	  - Optional output. A bar plot visualizing the estimated hashtag background probability distribution.
	* - output_name.real_content.hist.png
	  - Optional output. A histogram plot depicting hashtag distributions of not-real-cells and real-cells as defined by total number of expressed genes in the RNA assay.
	* - output_name.rna_demux.hist.png
	  - Optional output. A histogram plot depicting RNA UMI distribution for singlets, doublets and unknown cells.
	* - output_name.gene_name.violin.png
	  - Optional outputs. Violin plots depicting gender-specific gene expression across samples. We can have multiple plots if a gene list is provided in ``demuxEM_generate_gender_plot`` field of cumulus_hashing_cite_seq inputs.

In the output subfolder of each genetic-pooling RNA-hashtag data pair generated by *souporcell*, you can find the following files:

.. list-table::
	:widths: 5 10
	:header-rows: 1

	* - Name
	  - Description
	* - output_name_demux.zarr
	  - Demultiplexed RNA count matrix in zarr format. Please refer to section `load demultiplexing results into Python and R`_ for its structure.
	* - clusters.tsv
	  - Inferred droplet type and cluster assignment for each cell barcode.
	* - cluster_genotypes.vcf
	  - Inferred genotypes for each cluster.
	* - match_donors.log
	  - Log of matching donors step, with information of donor matching included.

In the output subfolder of each genetic-pooling RNA-hashtag data pair generated by *demuxlet*, you can find the following files:

.. list-table::
	:widths: 5 10
	:header-rows: 1

	* - Name
	  - Description
	* - output_name_demux.zarr
	  - Demultiplexed RNA count matrix in zarr format. Please refer to section `load demultiplexing results into Python and R`_ for its structure.
	* - output_name.best
	  - Inferred droplet type and cluster assignment for each cell barcode.

---------------------------------

Load demultiplexing results into Python and R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load demultiplexed RNA count matrix into Python, you need to install Python package `pegasusio <https://pypi.org/project/pegasusio/>`_ first. Then follow the codes below::

	import pegasusio as io
	data = io.read_input('output_name_demux.zarr')

Once you load the data object, you can find estimated droplet types (singlet/doublet/unknown) in ``data.obs['demux_type']``. Notices that there are cell barcodes with no sample associated, and therefore have no droplet type.

You can also find estimated sample assignments in ``data.obs['assignment']``.

For cell-hashing/nucleus-hashing data, you can find estimated sample fractions (sample1, sample2, ..., samplen, background) for each droplet in ``data.obsm['raw_probs']``.

To load the results into R, you need to install R package ``reticulate`` in addition to Python package ``pegasusio``. Then follow the codes below::

	library(reticulate)
	ad <- import("pegasusio", convert = FALSE)
	data <- ad$read_input("output_name_demux.zarr")

Results are in ``data$obs['demux_type']``, ``data$obs['assignment']``, and ``data$obsm['raw_probs']`` (this exists in cell-hashing/nucleus-hashing results only).


.. _cellranger_workflow tutorial: ./cellranger.html
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _genome reference: ./cellranger.html#sample-sheet
.. _souporcell: https://github.com/wheaton5/souporcell
.. _demuxlet: https://github.com/statgen/demuxlet
.. _pegasusio: https://pypi.org/project/pegasusio/
.. _load demultiplexing results into Python and R: ./demultiplexing.html#load-demultiplexing-results-into-python-and-r