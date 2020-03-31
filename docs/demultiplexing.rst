Demultiplex cell-hashing/nucleus-hashing data
------------------------------------------------------

This ``demultiplexing`` workflow generates gene-count matrices from cell-hashing/nucleus-hashing 10X data by demultiplexing.

Prepare input data and import workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Run ``cellranger_workflow`` to generate data
++++++++++++++++++++++++++++++++++++++++++++++++

	You'll need raw gene count matrices and hashtag data for demultiplexing, which are achieved by Cell Ranger tool. You can skip this step if you already have them.

	Please refer to the `cellranger_workflow tutorial`_ for details.

	When finished, you should be able to find the raw gene count matrix (e.g. ``raw_gene_bc_matrices_h5.h5``) and hashtag data (e.g. ``sample_1_ADT.csv`` or ``possorted_genome_bam.bam``, depending on which demultiplexing algorithm you use) for each sample.

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
		  - Google bucket url to the reference genotypes in ``vcf.gz`` format. This column is required only when *TYPE* is ``genetic-pooling`` and *demultiplexing_algorithm* is ``demuxlet``.


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
	  - demuxEM version to use. Currently only support "0.1".
	  - "0.1"
	  - "0.1"
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
	  - Google Bucket URL of output directory. Within it, each folder is for one RNA and hashtag data pair in the input sample sheet.


.. _cellranger_workflow tutorial: ./cellranger.html
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _genome reference: ./cellranger.html#sample-sheet
.. _souporcell: https://github.com/wheaton5/souporcell
.. _demuxlet: https://github.com/statgen/demuxlet