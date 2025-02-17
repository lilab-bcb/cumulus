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
		* - **GRCh38-2020-A**
		  - Human GRCh38 (GENCODE v32/Ensembl 98)
		* - **mm10-2020-A**
		  - Mouse mm10 (GENCODE vM23/Ensembl 98)
		* - **GRCh38_and_mm10-2020-A**
		  - Human GRCh38 (GENCODE v32/Ensembl 98) and mouse mm10 (GENCODE vM23/Ensembl 98)

#. *Chemistry* column.

	According to *cellranger count*'s documentation, chemistry can be

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

#. *Flowcell* column.


#. *DataType* column.

	This column is optional with a default **rna**. If you want to put a value, put **rna** here.

#. *FetureBarcodeFile* column.

	Put target panel CSV file here for targeted expressiond data. Note that if a target panel CSV is present, cell ranger version must be >= 4.0.0.

#. Example::

	Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
	sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,rna
	sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime,rna
	sample_2,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime,rna
	sample_2,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime,rna
	sample_3,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3,SI-TT-A1,auto,rna,gs://fc-e0000000-0000-0000-0000-000000000000/immunology_v1.0_GRCh38-2020-A.target_panel.csv

Workflow input
++++++++++++++

For sc/snRNA-seq data, ``cellranger_workflow`` takes Illumina outputs as input and runs ``cellranger mkfastq`` and ``cellranger count``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **input_csv_file**
		  - Sample Sheet (contains Sample, Reference, Flowcell, Lane, Index as required and Chemistry, DataType, FeatureBarcodeFile as optional)
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
		  - Results are written under directory *output_directory* and will overwrite any existing files at this location.
		* - force_cells
		  - Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells
		  - 6000
		  -
		* - expect_cells
		  - Expected number of recovered cells. Mutually exclusive with force_cells
		  - 3000
		  -
		* - include_introns
		  - Turn this option on to also count reads mapping to intronic regions. With this option, users do not need to use pre-mRNA references. Note that if this option is set, cellranger_version must be >= 5.0.0.
		  - true
		  - true
		* - no_bam
		  - Turn this option on to disable BAM file generation. This option is only available if cellranger_version >= 5.0.0.
		  - false
		  - false
		* - secondary
		  - Perform Cell Ranger secondary analysis (dimensionality reduction, clustering, etc.)
		  - false
		  - false
		* - cellranger_version
		  - cellranger version, could be: 9.0.0, 8.0.1, 8.0.0, 7.2.0, 7.1.0, 7.0.1, 7.0.0
		  - "9.0.0"
		  - "9.0.0"
		* - config_version
		  - config docker version used for processing sample sheets, could be 0.3, 0.2, 0.1
		  - "0.3"
		  - "0.3"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "quay.io/cumulus" for images on Red Hat registry;

		  	- "cumulusprod" for backup images on Docker Hub.
		  - "quay.io/cumulus"
		  - "quay.io/cumulus"
		* - mkfastq_docker_registry
		  - Docker registry to use for ``cellranger mkfastq``.
		    Default is the registry to which only Broad users have access.
		    See :ref:`bcl2fastq-docker` for making your own registry.
		  - "gcr.io/broad-cumulus"
		  - "gcr.io/broad-cumulus"
		* - acronym_file
		  - | The link/path of an index file in TSV format for fetching preset genome references, chemistry whitelists, etc. by their names.
		    | Set an GS URI if *backend* is ``gcp``; an S3 URI for ``aws`` backend; an absolute file path for ``local`` backend.
		  - "s3://xxxx/index.tsv"
		  - "gs://regev-lab/resources/cellranger/index.tsv"
		* - zones
		  - Google cloud zones
		  - "us-central1-a us-west1-a"
		  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - num_cpu
		  - Number of cpus to request for one node for cellranger mkfastq and cellranger count
		  - 32
		  - 32
		* - memory
		  - Memory size string for cellranger mkfastq and cellranger count
		  - "120G"
		  - "120G"
		* - mkfastq_disk_space
		  - Optional disk space in GB for mkfastq
		  - 1500
		  - 1500
		* - count_disk_space
		  - Disk space in GB needed for cellranger count
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

See the table below for important sc/snRNA-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cellranger_mkfastq.output_fastqs_directory
	  - Array[String]?
	  - Subworkflow output. A list of cloud urls containing FASTQ files, one url per flowcell.
	* - cellranger_count.output_count_directory
	  - Array[String]?
	  - Subworkflow output. A list of cloud urls containing gene count matrices, one url per sample.
	* - cellranger_count.output_web_summary
	  - Array[File]?
	  - Subworkflow output. A list of htmls visualizing QCs for each sample (cellranger count output).
	* - collect_summaries.metrics_summaries
	  - File?
	  - Task output. A excel spreadsheet containing QCs for each sample.
	* - count_matrix
	  - String
	  - Workflow output. Cloud url for a template count_matrix.csv to run Cumulus.


.. _Carly Ziegler: http://shaleklab.com/author/carly/
.. _[Kim et al. Cell 2020]: https://www.sciencedirect.com/science/article/pii/S0092867420304062
.. _10x single cell RNA-seq sample index set names: https://support.10xgenomics.com/single-cell-gene-expression/index/doc/specifications-sample-index-sets-for-single-cell-3
