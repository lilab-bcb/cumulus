To process scATAC-seq data, follow the specific instructions below.

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built scATAC-seq references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38-2020-A_arc_v2.0.0**
		  - Human GRCh38, cellranger-arc/atac reference 2.0.0
		* - **mm10-2020-A_arc_v2.0.0**
		  - Mouse mm10, cellranger-arc/atac reference 2.0.0
		* - **GRCh38_and_mm10-2020-A_atac_v2.0.0**
		  - Human GRCh38 and mouse mm10, cellranger-atac reference 2.0.0

#. *Chemistry* column.

	By default is **auto**, which will not specify a given chemistry. To analyze just the individual ATAC library from a 10x multiome assay using cellranger-atac count, use ``ARC-v1`` in the Chemistry column.

#. *DataType* column.

	Set it to **atac**.

#. *AuxFile* column.

	Leave it blank for scATAC-seq.

An example sample sheet is below::

	Sample,Reference,Flowcell,DataType
	sample_atac,GRCh38-2020-A_arc_v2.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9YB/Fastq,atac

Workflow input
++++++++++++++

``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files), and runs ``cellranger-atac count``. Please see the description of inputs below. Note that required inputs are shown in bold.

.. list-table::
	:widths: 5 30 30 20
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_csv_file**
	  - Sample Sheet (contains Sample, Reference, Flowcell as required and Chemistry, DataType, FeatureBarcodeFile as optional)
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
	  -
	* - **output_directory**
	  - Output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_atac_output"
	  -
	* - force_cells
	  - Force pipeline to use this number of cells, bypassing the cell detection algorithm
	  - 6000
	  -
	* - atac_dim_reduce
	  - Choose the algorithm for dimensionality reduction prior to clustering and tsne: "lsa", "plsa", or "pca"
	  - "lsa"
	  - "lsa"
	* - peaks
	  - A 3-column BED file of peaks to override cellranger atac peak caller. Peaks must be sorted by position and not contain overlapping peaks; comment lines beginning with ``#`` are allowed
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/common_peaks.bed"
	  -
	* - cellranger_atac_version
	  - cellranger-atac version. Available options: 2.1.0, 2.0.0
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
	* - atac_num_cpu
	  - Number of cpus for cellranger-atac count
	  - 64
	  - 64
	* - atac_memory
	  - Memory string for cellranger-atac count
	  - "57.6G"
	  - "57.6G"
	* - atac_disk_space
	  - Disk space in GB needed for cellranger-atac count
	  - 500
	  - 500
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

See the table below for important scATAC-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cellranger_atac_count.output_count_directory
	  - Array[String]
	  - Subworkflow output. A list of cloud URIs containing cellranger-atac count outputs, one URI per sample.
	* - cellranger_atac_count.output_web_summary
	  - Array[File]
	  - Subworkflow output. A list of htmls visualizing QCs for each sample (cellranger-atac count output).
	* - collect_summaries_atac.metrics_summaries
	  - File
	  - Task output. An Excel spreadsheet containing QCs for each sample.
