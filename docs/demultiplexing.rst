Demultiplex genetic-pooling/cell-hashing/nucleus-hashing sc/snRNA-Seq data
--------------------------------------------------------------------------

This ``demultiplexing`` workflow generates gene-count matrices from cell-hashing/nucleus-hashing/genetic-pooling data by demultiplexing.

In the workflow, ``demuxEM`` is used for analyzing cell-hashing/nucleus-hashing data, while ``souporcell`` and ``popscle`` (including *demuxlet* and *freemuxlet*) are for genetic-pooling data.

Prepare input data and import workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Run ``cellranger_workflow``
++++++++++++++++++++++++++++++++

	To demultiplex, you'll need raw gene count and hashtag matrices for cell-hashing/nucleus-hashing data, or raw gene count matrices and genome BAM files for genetic-pooling data. You can generate these data by running the ``cellranger_workflow``.

	Please refer to the `cellranger_workflow tutorial`_ for details.

	When finished, you should be able to find the raw gene count matrix (e.g. ``raw_gene_bc_matrices_h5.h5``), hashtag matrix (e.g. ``sample_1_ADT.csv``) / genome BAM file (e.g. ``possorted_genome_bam.bam``) for each sample.

2. Import ``demultiplexing``
++++++++++++++++++++++++++++++

	Import *demultiplexing* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose **github.com/klarman-cell-observatory/cumulus/Demultiplexing** to import.

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
		  - Google bucket url to the hashtag file generated in Step 1. The column name can be either *TagFile* or *ADT*, where *ADT* is for backward compatibility with older snapshots.
		* - **TYPE**
		  - Assay type, which can be ``cell-hashing``, ``nucleus-hashing``, or ``genetic-pooling``.
		* - Genotype
		  - Google bucket url to the reference genotypes in ``vcf.gz`` format. This column is **required** in the following cases:

		    - Run ``genetic-pooling`` assay with ``souporcell`` algorithm (i.e. *TYPE* is ``genetic-pooling``, *demultiplexing_algorithm* input is ``souporcell``):

		      - Run with reference genotypes, i.e. *souporcell_de_novo_mode* is ``false``.

		      - Run in *de novo* mode (i.e. *souporcell_de_novo_mode* is ``true``), but need to match the resulting cluster names by information from reference genotypes (see description of *souporcell_rename_donors* input below).

		    - Run ``genetic-pooling`` assay with ``popscle`` algorithm (i.e. *TYPE* is ``genetic-pooling``, *demultiplexing_algorithm* input is ``popscle``):

		      - *popscle_num_samples* input is ``0``. In this case, *demuxlet* will be run with reference genotypes.

		      - *popscle_num_samples* input is larger than ``0``. In this case, reference genotypes will be only used to generate pileups, then *freemuxlet* will be used for demultiplexing **without** reference genotypes.



	Example::

		OUTNAME,RNA,TagFile,TYPE,Genotype
		sample_1,gs://exp/data_1/raw_gene_bc_matrices_h5.h5,gs://exp/data_1/sample_1_ADT.csv,cell-hashing
		sample_2,gs://exp/data_2/raw_gene_bc_matrices_h5.h5,gs://exp/data_2/sample_2_ADT.csv,nucleus-hashing
		sample_3,gs://exp/data_3/raw_gene_bc_matrices_h5.h5,gs://exp/data_3/possorted_genome_bam.bam,genetic-pooling
		sample_4,gs://exp/data_4/raw_gene_bc_matrices_h5.h5,gs://exp/data_4/possorted_genome_bam.bam,genetic-pooling,gs://exp/variants/ref_genotypes.vcf.gz

	**3.2 Upload your sample sheet to the workspace bucket:**

	Use gsutil_ (you already have it if you've installed Google Cloud SDK) in your unix terminal to upload your sample sheet to workspace bucket.

	Example::

			gsutil cp /foo/bar/projects/sample_sheet_demux.csv gs://fc-e0000000-0000-0000-0000-000000000000/

---------------

Workflow inputs
^^^^^^^^^^^^^^^^

Below are inputs for *demultiplexing* workflow. We'll first introduce global inputs, and then inputs for each of the demultiplexing tools. Notice that required inputs are in bold.

global inputs
+++++++++++++++


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
	* - genome
	  - Reference genome name. Its usage depends on the assay type:

	    - For *cell-hashing* or *nucleus-hashing*, only write this name as an annotation into the resulting count matrix file.

	    - For *genetic-pooling*, if *demultiplexing_algorithm* input is ``souporcell``, you should choose one name from this `genome reference`_ list.

	    - For *genetic-pooling*, if *demultiplexing_algorithm* input is ``popscle``, reference genome name is not needed.
	  - "GRCh38"
	  -
	* - demultiplexing_algorithm
	  - demultiplexing algorithm to use for *genetic-pooling* data. Options:

	  	- "souporcell": Use souporcell_, a reference-genotypes-free algorithm for demultiplexing droplet scRNA-Seq data.

	  	- "popscle": Use popscle_, a canonical algorithm for demultiplexing droplet scRNA-Seq data, including *demuxlet* (with reference genotypes) and *freemuxlet* (reference-genotype-free) components.
	  - "souporcell"
	  - "souporcell"
	* - min_num_genes
	  - Only demultiplex cells/nuclei with at least <min_num_genes> expressed genes
	  - 100
	  - 100
	* - zones
	  - Google cloud zones to consider for execution.
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - docker_registry
	  - Docker registry to use.

	  	- "quay.io/cumulus" for images on Red Hat registry;

	  	- "cumulusprod" for backup images on Docker Hub.
	  - "quay.io/cumulus"
	  - "quay.io/cumulus"
	* - config_version
	  - Version of config docker image to use. This docker is used for parsing the input sample sheet for downstream execution. Available options: ``0.2``, ``0.1``.
	  - "0.2"
	  - "0.2"
	* - backend
	  - Cloud infrastructure backend to use. Available options:

	  	- "gcp" for Google Cloud;
	  	- "aws" for Amazon AWS;
	  	- "local" for local machine.
	  - "gcp"
	  - "gcp"
	* - ref_index_file
	  - The link/path of an index file in TSV format for fetching preset genome references, chemistry whitelists, etc. by their names. Set an GS URI if backend is ``gcp``; an S3 URI for ``aws`` backend; an absolute file path for ``local`` backend.
	  - "s3://xxxx/index.tsv"
	  - "gs://regev-lab/resources/cellranger/index.tsv"
	* - preemptible
	  - Number of maximum preemptible tries allowed. This works only when *backend* is ``gcp``.
	  - 2
	  - 2
	* - awsMaxRetries
	  - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
	  - 5
	  - 5

demuxEM inputs
++++++++++++++++

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1


	* - Name
	  - Description
	  - Example
	  - Default
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
	  - demuxEM version to use. Choose from "0.1.7", "0.1.6" and "0.1.5".
	  - "0.1.7"
	  - "0.1.7"
	* - demuxEM_num_cpu
	  - demuxEM parameter. Number of CPUs to request for demuxEM per pair.
	  - 8
	  - 8
	* - demuxEM_memory
	  - demuxEM parameter. Memory size string for demuxEM per pair.
	  - "10G"
	  - "10G"
	* - demuxEM_disk_space
	  - demuxEM parameter. Disk space (integer) in GB needed for demuxEM per pair.
	  - 20
	  - 20

souporcell inputs
++++++++++++++++++

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1


	* - Name
	  - Description
	  - Example
	  - Default
	* - souporcell_version
	  - souporcell version to use. Available versions:

	    - ``2021.03``: Based on commitment `1bd9f1 <https://github.com/wheaton5/souporcell/tree/1bd9f11d70eaee6ac14713de09c377c285ca2787>`_ on 2021/03/07.

	    - ``2020.07``: Based on commitment `0d09fb <https://github.com/wheaton5/souporcell/tree/0d09fbe26d878adb294b536c4f41a7718c0d0f9d>`_ on 2020/07/27.

	    - ``2020.03``: Based on commitment `eeddcd <https://github.com/wheaton5/souporcell/tree/eeddcde5892c5cbf8aba2149f0e77756f830a5ae>`_ on 2020/03/31.
	  - "2021.03"
	  - "2021.03"
	* - souporcell_num_clusters
	  - | souporcell parameter. Number of expected clusters when doing clustering.
	    | **This needs to be set when running souporcell.**
	  - 8
	  - 1
	* - souporcell_de_novo_mode
	  - souporcell parameter.

	    - If ``true``, run souporcell in de novo mode without reference genotypes:

		  - If input *souporcell_common_variants* is further provided, use this common variants list instead of calling SNPs de novo.

		  - If a reference genotype vcf file is provided in the sample sheet, use it **only** for matching the cluster labels computed by souporcell.

	    - If ``false``, run souporcell with ``--known_genotypes`` option using the reference genotype vcf file specified in sample sheet.
	  - true
	  - true
	* - souporcell_num_clusters
	  - | souporcell parameter. Number of expected clusters when doing clustering.
	    | **This needs to be set when running souporcell.**
	  - 8
	  - 1
	* - souporcell_common_variants
	  - | souporcell parameter. Users can provide a common variants list in VCF format for Souporcell to use, instead of calling SNPs de novo.
	    | **Notice:** This input is enabled only when *souporcell_de_novo_mode* is ``false``.
	  - "1000genome.common.variants.vcf.gz"
	  -
	* - souporcell_skip_remap
	  - souporcell parameter. Skip remap step. Only recommended in non denovo mode or common variants are provided.
	  - true
	  - false
	* - souporcell_rename_donors
	  - souporcell parameter. A comma-separated list of donor names for matching clusters achieved by souporcell. Must be consistent with *souporcell_num_clusters* input.

	    - If this input is empty, use cluster labels from the reference genotype vcf file if provided in the sample sheet; if this vcf file is not provided, simply name clusters as *Donor1*, *Donor2*, ...

	    - If this input is not empty, and a reference genotype vcf file is provided in the sample sheet, first match the cluster labels using those from this vcf file, then rename to donor names specified in this input.

	    - If this input is not empty, and **NO** reference genotype vcf file is provided in the sample sheet, simply match the cluster labels in one-to-one correspondence with donor names specified in this input.
	  - "CB1,CB2,CB3,CB4"
	  -
	* - souporcell_num_cpu
	  - souporcell parameter. Number of CPUs to request for souporcell per pair.
	  - 32
	  - 32
	* - souporcell_memory
	  - souporcell parameter. Memory size string for souporcell per pair.
	  - "120G"
	  - "120G"
	* - souporcell_disk_space
	  - souporcell parameter. Disk space (integer) in GB needed for souporcell per pair.
	  - 500
	  - 500

Popscle inputs
+++++++++++++++++

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1


	* - Name
	  - Description
	  - Example
	  - Default
	* - popscle_num_samples
	  - popscle parameter. Number of samples to be multiplexed together:

	    - If ``0``, run with *demuxlet* using reference genotypes.

	    - Otherwise, run with *freemuxlet* in de novo mode without reference genotypes.
	  - 4
	  - 0
	* - popscle_min_MQ
	  - popscle parameter. Minimum mapping quality to consider (lower MQ will be ignored).
	  - 20
	  - 20
	* - popscle_min_TD
	  - popscle parameter. Minimum distance to the tail (lower will be ignored).
	  - 0
	  - 0
	* - popscle_tag_group
	  - popscle parameter. Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use ``CB``.
	  - "CB"
	  - "CB"
	* - popscle_tag_UMI
	  - popscle parameter. Tag representing UMIs. For 10x genomics, use ``UB``.
	  - "UB"
	  - "UB"
	* - popscle_field
	  - popscle parameter. FORMAT field to extract from: genotype (``GT``), genotype likelihood (``GL``), or posterior probability (``GP``).
	  - "GT"
	  - "GT"
	* - popscle_alpha
	  - popscle parameter. Grid of alpha to search for, in a comma separated list format of all alpha values to be considered.
	  - "0.1,0.2,0.3,0.4,0.5"
	  - "0.1,0.2,0.3,0.4,0.5"
	* - popscle_rename_donors
	  - | popscle parameter. A comma-separated list of donor names for renaming clusters achieved by popscle. Must be consistent with *popscle_num_samples* input.
	    | By default, the resulting donors are *Donor1*, *Donor2*, ...
	  - "CB1,CB2,CB3,CB4"
	  -
	* - popscle_version
	  - popscle parameter. popscle version to use. Available options:

	    - ``2021.05``: Based on commitment `da70fc7 <https://github.com/statgen/popscle/tree/da70fc78da385ef049e0e890342acfd62842cae0>`_ on 2021/05/05.

	    - ``0.1b``: Based on version `0.1-beta <https://github.com/statgen/popscle/releases/tag/v0.1-beta>`_ released on 2019/10/03.
	  - "2021.05"
	  - "2021.05"
	* - popscle_num_cpu
	  - popscle parameter. Number of CPU used by popscle per pair.
	  - 1
	  - 1
	* - popscle_memory
	  - popscle parameter. Memory size string per pair.
	  - "120G"
	  - "120G"
	* - popscle_extra_disk_space
	  - popscle parameter. Extra disk space size (integer) in GB needed for popscle per pair, besides the disk size required to hold input files specified in the sample sheet.
	  - 100
	  - 100

---------------------

Workflow outputs
^^^^^^^^^^^^^^^^^^

See the table below for *demultiplexing* workflow outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_folders
	  - Array[String]
	  - A list of Google Bucket URLs of the output folders. Each folder is associated with one RNA-hashtag pair in the given sample sheet.
	* - output_zarr_files
	  - Array[File]
	  - A list of demultiplexed RNA count matrices in zarr format. Each zarr file is associated with one RNA-hashtag pair in the given sample sheet. Please refere to section `load demultiplexing results into Python and R`_ for its structure.

In the output subfolder of each cell-hashing/nuclei-hashing RNA-hashtag data pair, you can find the following files:

.. list-table::
	:widths: 5 10
	:header-rows: 1

	* - Name
	  - Description
	* - output_name_demux.zarr.zip
	  - Demultiplexed RNA raw count matrix in zarr format. Please refer to section `load demultiplexing results into Python and R`_ for its structure.
	* - output_name.out.demuxEM.zarr.zip
	  - | This file contains intermediate results for both RNA and hashing count matrices.
	    | To load this file into Python, you need to first install `Pegasusio`_ on your local machine. Then use ``import pegasusio as io; data = io.read_input("output_name.out.demuxEM.zarr.zip")`` in Python environment.
	    | It contains 2 UnimodalData objects: one with key name suffix ``-hashing`` is the hashtag count matrix, the other one with key name suffix ``-rna`` is the demultiplexed RNA count matrix.
	    | To load the hashtag count matrix, type ``hash_data = data.get_data('<genome>-hashing')``, where ``<genome>`` is the genome name of the data. The count matrix is ``hash_data.X``; cell barcode attributes are stored in ``hash_data.obs``; sample names are in ``hash_data.var_names``. Moreover, the estimated background probability regarding hashtags is in ``hash_data.uns['background_probs']``.
	    | To load the RNA matrix, type ``rna_data = data.get_data('<genome>-rna')``, where ``<genome>`` is the genome name of the data. It only contains cells which have estimated sample assignments. The count matrix is ``rna_data.X``. Cell barcode attributes are stored in ``rna_data.obs``: ``rna_data.obs['demux_type']`` stores the estimated droplet types (singlet/doublet/unknown) of cells; ``rna_data.obs['assignment']`` stores the estimated hashtag(s) that each cell belongs to. Moreover, for cell-hashing/nucleus-hashing data, you can find estimated sample fractions (sample1, sample2, ..., samplen, background) for each droplet in ``rna_data.obsm['raw_probs']``.
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
	* - output_name_demux.zarr.zip
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
	* - output_name_demux.zarr.zip
	  - Demultiplexed RNA count matrix in zarr format. Please refer to section `load demultiplexing results into Python and R`_ for its structure.
	* - output_name.best (demuxlet) or output_name.clust1.samples.gz (freemuxlet)
	  - Inferred droplet type and cluster assignment for each cell barcode.

---------------------------------

Load demultiplexing results into Python and R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load demultiplexed RNA count matrix into Python, you need to install Python package pegasusio_ first. Then follow the codes below::

	import pegasusio as io
	data = io.read_input('output_name_demux.zarr.zip')

Once you load the data object, you can find estimated droplet types (singlet/doublet/unknown) in ``data.obs['demux_type']``. Notices that there are cell barcodes with no sample associated, and therefore have no droplet type.

You can also find estimated sample assignments in ``data.obs['assignment']``.

For cell-hashing/nucleus-hashing data, if one sample name can correspond to multiple feature barcodes, each feature barcode is assigned to a unique sample name, and this deduplicated sample assignment results are in ``data.obs['assignment.dedup']``.

To load the results into R, you need to install R package ``reticulate`` in addition to Python package ``pegasusio``. Then follow the codes below::

	library(reticulate)
	ad <- import("pegasusio", convert = FALSE)
	data <- ad$read_input("output_name_demux.zarr.zip")

Results are in ``data$obs['demux_type']``, ``data$obs['assignment']``, and similarly as above, for cell-hashing/nucleus-hashing data, you'll find an additional field ``data$obs['assignment.dedup']`` for deduplicated sample assignment in the case that one sample name can correspond to multiple feature barcodes.


.. _cellranger_workflow tutorial: ./cellranger/index.html
.. _Import workflows to Terra: ./cumulus_import.html
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _genome reference: ./cellranger/index.html#sample-sheet
.. _souporcell: https://github.com/wheaton5/souporcell
.. _popscle: https://github.com/statgen/popscle
.. _pegasusio: https://pypi.org/project/pegasusio/
.. _load demultiplexing results into Python and R: ./demultiplexing.html#load-demultiplexing-results-into-python-and-r
