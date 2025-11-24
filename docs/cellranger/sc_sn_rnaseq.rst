.. note::
	Cell Ranger will send anonymized telemetry data to 10x Genomics starting from v9.0. Here is the details on `Cell Ranger Pipeline Telemetry`_.

	This option has been turned off in this *cellranger_workflow*, thus **no data will be sent to 10x Genomics**.

To process sc/snRNA-seq data, follow the specific instructions below.

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built scRNA-seq references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38-2024-A**
		  - Human GRCh38, comparable to cellranger reference 2024-A (GENCODE v44/Ensembl 110). *Notice: This reference only supports Cell Ranger v6.0.0+.*
		* - **GRCm39-2024-A**
		  - Mouse GRCm39, comparable to cellranger reference 2024-A (GENCODE vM33/Ensembl 110). *Notice: This reference only supports Cell Ranger v6.0.0+.*
		* - **GRCh38_and_GRCm39-2024-A**
		  - Human GRCh38 (v44/Ensembl 110) and mouse GRCm39 (GENCODE vM33/Ensembl 110). *Notice: This reference only supports Cell Ranger v6.0.0+.*
		* - **mRatBN7.2-2024-A**
		  - Rat mRatBN7.2 reference.
		* - **GRCh38-2020-A**
		  - Human GRCh38 (GENCODE v32/Ensembl 98)
		* - **mm10-2020-A**
		  - Mouse mm10 (GENCODE vM23/Ensembl 98)
		* - **GRCh38_and_mm10-2020-A**
		  - Human GRCh38 (GENCODE v32/Ensembl 98) and mouse mm10 (GENCODE vM23/Ensembl 98)

#. *Chemistry* column.

	The cellranger workflow fully supports all 10x assay configurations. The most widely used ones are listed below:

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Chemistry
		  - Explanation
		* - **auto**
		  - autodetection (default). If the index read has extra bases besides cell barcode and UMI, autodetection might fail. In this case, please specify the chemistry
		* - **threeprime**
		  - Single Cell 3′
		* - **fiveprime**
		  - Single Cell 5′
		* - **ARC-v1**
		  - Gene Expression portion of 10x Multiome data

	Please refer to the section of ``--chemistry`` option in `Cell Ranger Command Line Arguments`_ for all other valid chemistry keywords.

#. *Flowcell* column.

	See the `table <./index.html#prepare-a-sample-sheet>`_ in general steps section above.

	.. note::
		The workflow accepts input in TAR files which contain FASTQ files inside, and can automatically handle such cases.

#. *DataType* column.

	This column is optional with a default **rna**. If you want to put a value, put **rna** here.

#. Example::

	Sample,Reference,Flowcell,Chemistry,DataType
	sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,threeprime,rna
	sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/Fastq,threeprime,rna
	sample_2,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,fiveprime,rna
	sample_2,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/Fastq,fiveprime,rna
	sample_3,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,auto,rna

Workflow input
++++++++++++++

For sc/snRNA-seq data, ``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files), and runs ``cellranger count``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **input_csv_file**
		  - Sample Sheet (contains Sample, Reference, Flowcell, Chemistry, DataType) in CSV format
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
		  -
		* - **output_directory**
		  - Cloud URI of the output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
		  - Results are written under directory *output_directory* and will overwrite any existing files at this location.
		* - include_introns
		  - Turn this option on to also count reads mapping to intronic regions. With this option, users do not need to use pre-mRNA references. Note that if this option is set, cellranger_version must be >= 5.0.0.
		  - true
		  - true
		* - no_bam
		  - Turn this option on to disable BAM file generation. This option is only available if cellranger_version >= 5.0.0.
		  - false
		  - false
		* - expect_cells
		  - Expected number of recovered cells. Mutually exclusive with force_cells
		  - 3000
		  -
		* - force_cells
		  - Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells
		  - 6000
		  -
		* - secondary
		  - Perform Cell Ranger secondary analysis (dimensionality reduction, clustering, etc.)
		  - false
		  - false
		* - cellranger_version
		  - cellranger version, could be: 10.0.0, 9.0.1, 8.0.1, 7.2.0
		  - "10.0.0"
		  - "10.0.0"
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
		  - Number of cpus to request for one node for cellranger count
		  - 32
		  - 32
		* - memory
		  - Memory size string for cellranger count
		  - "120G"
		  - "120G"
		* - count_disk_space
		  - Disk space in GB needed for cellranger count
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

See the table below for important sc/snRNA-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cellranger_count.output_count_directory
	  - Array[String]
	  - Subworkflow output. A list of cloud URIs containing gene count matrices, one URI per sample.
	* - cellranger_count.output_web_summary
	  - Array[File]
	  - Subworkflow output. A list of htmls visualizing QCs for each sample (cellranger count output).
	* - collect_summaries.metrics_summaries
	  - File
	  - Task output. An excel spreadsheet containing QCs for each sample.


.. _Cell Ranger Pipeline Telemetry: https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-pipeline-telemetry
.. _Cell Ranger Command Line Arguments: https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments
