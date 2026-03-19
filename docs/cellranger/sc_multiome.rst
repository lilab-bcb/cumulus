.. note::
	Cell Ranger ARC will send anonymized telemetry data to 10x Genomics starting from v2.1. Here is the details on `Cell Ranger ARC Pipeline Telemetry`_.

	This option has been turned off in this *cellranger_workflow*, thus **no data will be sent to 10x Genomics**.

To process `10x Multiome`_ (GEX + ATAC) data, follow the instructions below:

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built single-cell Multiome ATAC + Gene Expression references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38-2024-A_arc**
		  - Human GRCh38 (GENCODE v44/Ensembl 110) for cellranger arc
		* - **GRCm39-2024-A_arc**
		  - Mouse GRCm39 (GENCODE vM33/Ensembl 110) for cellranger arc
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


In the above example, the linked samples will be processed altogether. And the output will be one subfolder named ``sample1``.

Workflow input
++++++++++++++

For single-cell multiomics data, ``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files). Revalant workflow inputs are described below, with required inputs highlighted in bold.

.. list-table::
	:widths: 5 30 30 20
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_csv_file**
	  - Sample Sheet (contains Sample, Reference, Flowcell, Chemistry, DataType, and Link)
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
	  -
	* - **output_directory**
	  - Output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
	  -
	* - include_introns
	  - Turn this option on to also count reads mapping to intronic regions. With this option, users do not need to use pre-mRNA references.
	  - true
	  - true
	* - secondary
	  - | Perform Cell Ranger ARC secondary analysis (e.g. clustering).
	    | **Note:** This parameter works only for *cellranger_arc_version* ``2.1.0`` or later.
	  - false
	  - false
	* - no_bam
	  - Turn this option on to disable BAM file generation.
	  - false
	  - false
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
	* - cellranger_arc_version
	  - | cellranger-arc version, could be: ``2.1.0``, ``2.0.2.strato`` (compatible with workflow v2.6.1+), ``2.0.2.custom-max-cell`` (with max_cell threshold set to 80,000), ``2.0.2`` (compatible with workflow v2.6.0 or earlier), ``2.0.1``, ``2.0.0``
	  	| **Note:** The 20,000 total cell limit has been removed since version ``2.1.0``.
	  - "2.1.0"
	  - "2.1.0"
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
	* - arc_num_cpu
	  - Number of cpus to request for one link
	  - 64
	  - 64
	* - arc_memory
	  - Memory size string for one link
	  - "160G"
	  - "160G"
	* - arc_disk_space
	  - Disk space in GB needed for one link
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
	  - A list of cloud URIs to output, one URI per link
	* - cellranger_arc_count.output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each link
	* - collect_summaries_arc.metrics_summaries
	  - File
	  - An excel spreadsheet containing QCs for each link


.. _10x Multiome: https://www.10xgenomics.com/support/epi-multiome
.. _Cell Ranger ARC Pipeline Telemetry: https://www.10xgenomics.com/support/software/cell-ranger-arc/latest/tutorials/cr-arc-pipeline-telemetry
