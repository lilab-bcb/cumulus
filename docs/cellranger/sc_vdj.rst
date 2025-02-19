To process single-cell immune profiling (scIR-seq) data, follow the specific instructions below.

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built scIR-seq references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCH38_vdj_v7.1.0**
		  - Human GRCh38 V(D)J sequences, cellranger reference 7.1.0, annotation built from Ensembl *Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf*
		* - **GRCh38_vdj_v7.0.0**
		  - Human GRCh38 V(D)J sequences, cellranger reference 7.0.0, annotation built from Ensembl *Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf*
		* - **GRCm38_vdj_v7.0.0**
		  - Mouse GRCm38 V(D)J sequences, cellranger reference 7.0.0, annotation built from Ensembl *Mus_musculus.GRCm38.94.gtf*

#. *Chemistry* column.

	This column is not used for scIR-seq data. Put **fiveprime** here as a placeholder if you decide to include the Chemistry column.

#. *DataType* column.

	Set it to **vdj**.

#. *FetureBarcodeFile* column.

	Leave it blank for scIR-seq.

#. Example::

	Sample,Reference,Flowcell,Chemistry,DataType
	sample_vdj,GRCh38_vdj_v7.1.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ/Fastq,fiveprime,vdj

Workflow input
++++++++++++++

For scIR-seq data, ``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files), and runs ``cellranger vdj``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

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
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
	  -
	* - vdj_denovo
	  - Do not align reads to reference V(D)J sequences before de novo assembly
	  - false
	  - false
	* - vdj_chain
	  - Force the analysis to be carried out for a particular chain type. The accepted values are:

		- "auto" for auto detection based on TR vs IG representation;

		- "TR" for T cell receptors;

		- "IG" for B cell receptors.
	  - "auto"
	  - "auto"
	* - vdj_inner_enrichment_primers
	  - | Cloud URI to a text file with custom inner enrichment primers. By default, use those provided in the 10x kit.
	    | If inner enrichment primers other than those provided in the 10x kits are used, they need to be specified here as a textfile with one primer per line. Disable secondary analysis, e.g. clustering
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/vdj_primers.txt"
	  -
	* - cellranger_version
	  - cellranger version, could be: 8.0.1, 8.0.0, 7.2.0, 7.1.0, 7.0.1, 7.0.0, 6.1.2, 6.1.1, 6.0.2, 6.0.1, 6.0.0, 5.0.1, 5.0.0
	  - "8.0.1"
	  - "8.0.1"
	* - docker_registry
	  - Docker registry to use for cellranger_workflow. Options:

	  	- "quay.io/cumulus" for images on Red Hat registry;

	  	- "cumulusprod" for backup images on Docker Hub.
	  - "quay.io/cumulus"
	  - "quay.io/cumulus"
	* - acronym_file
	  - | The link/path of an index file in TSV format for fetching preset genome references, chemistry barcode inclusion lists, etc. by their names.
	    | Set an GS URI if *backend* is ``gcp``; an S3 URI for ``aws`` backend; an absolute file path for ``local`` backend.
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
	  - Memory size string for cellranger vdj
	  - "120G"
	  - "120G"
	* - vdj_disk_space
	  - Disk space in GB needed for cellranger vdj
	  - 500
	  - 500
	* - backend
	  - Cloud backend for file transfer. Available options:

	  	- "gcp" for Google Cloud;
	  	- "aws" for Amazon AWS;
	  	- "local" for local machine.
	  - "gcp"
	  - "gcp"
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2
	* - awsQueueArn
	  - The AWS ARN string of the job queue to be used. This only works for ``aws`` backend.
	  - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
	  - ""


Workflow output
+++++++++++++++

See the table below for important scIR-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cellranger_vdj.output_count_directory
	  - Array[String]
	  - Subworkflow output. A list of cloud URIs containing vdj results, one URI per sample.
	* - cellranger_vdj.output_web_summary
	  - Array[File]
	  - Subworkflow output. A list of htmls visualizing QCs for each sample (cellranger vdj output).
	* - collect_summaries_vdj.metrics_summaries
	  - File
	  - Task output. An excel spreadsheet containing QCs for each sample.


.. _10x single cell V(D)J sample index set names: https://www.10xgenomics.com/support/single-cell-immune-profiling/documentation/steps/sequencing/sample-index-sets-for-single-cell-immune-profiling
