.. note::
	Cell Ranger will send anonymized telemetry data to 10x Genomics starting from v9.0. Here is the details on `Cell Ranger Pipeline Telemetry`_.

	This option has been turned off in this *cellranger_workflow*, thus **no data will be sent to 10x Genomics**.


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
		* - **GRCh38_vdj_v7.1.0**
		  - Human GRCh38 V(D)J sequences, cellranger reference 7.1.0, annotation built from Ensembl *Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf*
		* - **GRCh38_vdj_v7.0.0**
		  - Human GRCh38 V(D)J sequences, cellranger reference 7.0.0, annotation built from Ensembl *Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf*
		* - **GRCm38_vdj_v7.0.0**
		  - Mouse GRCm38 V(D)J sequences, cellranger reference 7.0.0, annotation built from Ensembl *Mus_musculus.GRCm38.94.gtf*

#. *Chemistry* column.

	This column is not used for scIR-seq data. Put **fiveprime** here as a placeholder if you decide to include the Chemistry column.

#. *DataType* column.

	Choose one from the availabe types below:

	* **vdj**: The VDJ library. Let the workflow auto-detect the chain type.
	* **vdj_t**: The VDJ-T library for T-cell receptor sequences.
	* **vdj_b**: The VDJ-B library for B-cell receptor sequences.
	* **vdj_t_gd**: The VDJ-T-GD library for T-cell receptor enriched for gamma (TRG) and delta (TRD) chains.

#. *AuxFile* column.

	Only need for **vdj_t_gd** type samples which use primer sequences to enrich cDNA for V(D)J sequences.
	In this case, provide a ``.txt`` file containing such sequences, one per line. Then this file would be given to ``--inner-enrichment-primers`` option in *cellranger vdj*.

.. note::
	The ``--chain`` option in ``cellranger vdj`` is automatically decided based on the *DataType* value specified:
		* For **vdj**: set to ``--chain auto``
		* For **vdj_t** and **vdj_t_gd**: set to ``--chain TR``
		* For **vdj_b**: set to ``--chain IG``

An example sample sheet is below::

	Sample,Reference,Flowcell,Chemistry,DataType,AuxFile
	sample1,GRCh38_vdj_v7.1.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ/Fastq,fiveprime,vdj,
	sample2,GRCh38_vdj_v7.1.0,gs://my-bucket/s2_fastqs,,vdj_t_gd,gs://my-bucket/s2_enrich_primers.txt

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
	  - Sample Sheet (contains Sample, Reference, Flowcell, DataType, Chemistry, and AuxFile)
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
