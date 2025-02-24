To utilize ``cellranger arc`` for single-cell multiomics, follow the specific instructions below. Note that cumulus_feature_barcoding/demuxEM would not be triggered for hashing/citeseq in this setting.

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built Multiome ATAC + Gene Expression references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38-2020-A_arc_v2.0.0**
		  - Human GRCh38 sequences (GENCODE v32/Ensembl 98), cellranger arc reference 2.0.0
		* - **mm10-2020-A_arc_v2.0.0**
		  - Mouse GRCm38 sequences (GENCODE vM23/Ensembl 98), cellranger arc reference 2.0.0


#. *Chemistry* column.

	By default is **auto**, which will not specify a given chemistry.

#. *DataType* column.

	For each sample, choose a data type from the table below:

	.. list-table::
		:widths: 5 10
		:header-rows: 1

		* - DataType
		  - Description
		* - **rna**
		  - For scRNA-Seq modality of the data
		* - **atac**
		  - For scATAC-Seq modality of the data

#. *AuxFile* column.

	Leave it blank.

#. *Link* column.

	Put a unique link name for all modalities that are linked. **Notice:** The *Link* name must be different from all *Sample* column values.

#. Example::

	Link,Sample,Reference,Flowcell,DataType
	sample1,s1_rna,GRCh38-2020-A_arc_v2.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ/Fastq,rna
	sample1,s1_atac,GRCh38-2020-A_arc_v2.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ/Fastq,atac


In the above example, the linked samples will be processed by *cellranger arc* altogether.

Workflow input
++++++++++++++

For single-cell multiomics data, ``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files), and runs ``cellranger-arc ount``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

.. list-table::
	:widths: 5 30 30 20
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_csv_file**
	  - Sample Sheet (contains Sample, Reference, Flowcell as required and Chemistry, DataType, FeatureBarcodeFile, Link as optional)
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
	  -
	* - **output_directory**
	  - Output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
	  -
	* - include_introns
	  - Turn this option on to also count reads mapping to intronic regions. With this option, users do not need to use pre-mRNA references. Note that if this option is set, cellranger_version must be >= 5.0.0. This option is used by *cellranger multi* and *cellranger count*.
	  - true
	  - true
	* - no_bam
	  - Turn this option on to disable BAM file generation. This option is only available if cellranger_version >= 5.0.0. This option is used by *cellranger-arc count*, *cellranger multi* and *cellranger count*.
	  - false
	  - false
	* - expect_cells
	  - Expected number of recovered cells. Mutually exclusive with force_cells. This option is used by *cellranger multi* and *cellranger count*.
	  - 3000
	  -
	* - force_cells
	  - Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. This option is used by *cellranger multi* and *cellranger count*.
	  - 6000
	  -
	* - arc_gex_exclude_introns
	  - | Disable counting of intronic reads. In this mode, only reads that are exonic and compatible with annotated splice junctions in the reference are counted.
	    | **Note:** using this mode will reduce the UMI counts in the feature-barcode matrix.
	  - false
	  - false
	* - arc_min_atac_count
	  - | Cell caller override to define the minimum number of ATAC transposition events in peaks (ATAC counts) for a cell barcode.
	    | **Note:** this input must be specified in conjunction with ``arc_min_gex_count`` input.
	    | With both inputs set, a barcode is defined as a cell if it contains at least ``arc_min_atac_count`` ATAC counts AND at least ``arc_min_gex_count`` GEX UMI counts.
	  - 100
	  -
	* - arc_min_gex_count
	  - | Cell caller override to define the minimum number of GEX UMI counts for a cell barcode.
	    | **Note:** this input must be specified in conjunction with ``arc_min_atac_count``. See the description of ``arc_min_atac_count`` input for details.
	  - 200
	  -
	* - peaks
	  - A 3-column BED file of peaks to override cellranger arc peak caller. Peaks must be sorted by position and not contain overlapping peaks; comment lines beginning with ``#`` are allowed
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/common_peaks.bed"
	  -
	* - secondary
	  - Perform Cell Ranger secondary analysis (dimensionality reduction, clustering, etc.). This option is used by *cellranger multi* and *cellranger count*.
	  - false
	  - false
	* - cmo_set
	  - CMO set CSV file, delaring CMO constructs and associated barcodes. See `CMO reference`_ for details. Used only for *cellranger multi*.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cmo_set.csv"
	  -
	* - cellranger_arc_version
	  - cellranger-arc version, could be: ``2.0.2.strato`` (compatible with workflow v2.6.1+), ``2.0.2.custom-max-cell`` (with max_cell threshold set to 80,000), ``2.0.2`` (compatible with workflow v2.6.0 or earlier), ``2.0.1``, ``2.0.0``
	  - "2.0.2.strato"
	  - "2.0.2.strato"
	* - docker_registry
	  - Docker registry to use for cellranger_workflow. Options:

	  	- "quay.io/cumulus" for images on Red Hat registry;

	  	- "cumulusprod" for backup images on Docker Hub.
	  - "quay.io/cumulus"
	  - "quay.io/cumulus"
	* - acronym_file
	  - | The link/path of an index file in TSV format for fetching preset genome references, chemistry barcode inclusion lists, etc. by their names.
	    | Set an GS URI if running on GCP; an S3 URI for AWS; an absolute file path for HPC or local machines.
	  - "s3://xxxx/index.tsv"
	  - "gs://cumulus-ref/resources/cellranger/index.tsv"
	* - zones
	  - Google cloud zones. For GCP Batch backend, the zones are automatically restricted by the Batch settings.
	  - "us-central1-a us-west1-a"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - num_cpu
	  - Number of cpus to request for one node for cellranger vdj
	  - 32
	  - 32
	* - memory
	  - Memory size string for cellranger and cellranger vdj
	  - "120G"
	  - "120G"
	* - count_disk_space
	  - Disk space in GB needed for cellranger count
	  - 500
	  - 500
	* - arc_num_cpu
	  - Number of cpus to request for one node for cellranger-arc count
	  - 64
	  - 64
	* - arc_memory
	  - Memory size string for cellranger-arc count
	  - "160G"
	  - "160G"
	* - arc_disk_space
	  - Disk space in GB needed for cellranger-arc count
	  - 700
	  - 700
	* - preemptible
	  - Number of preemptible tries. Only works for GCP
	  - 2
	  - 2
	* - awsQueueArn
	  - The AWS ARN string of the job queue to be used. Only works for AWS
	  - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
	  - ""

Workflow output
+++++++++++++++

See the table below for important output:

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cellranger_arc_count.output_count_directory
	  - Array[String]
	  - Subworkflow output. A list of cloud URIs containing *cellranger-arc count* output, one URI per *Link* name.
	* - cellranger_arc_count.output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each *Link* name.
	* - collect_summaries_arc.metrics_summaries
	  - File
	  - An excel spreadsheet containing QCs for each *Link* name.


.. _Feature Barcode Reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/feature-bc-analysis#feature-ref
.. _CMO reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cmoreference
