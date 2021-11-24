We provide routines wrapping Cell Ranger tools to build references for sc/snRNA-seq, scATAC-seq and single-cell immune profiling data.

Build references for sc/snRNA-seq
+++++++++++++++++++++++++++++++++

We provide a wrapper of ``cellranger mkref`` to build sc/snRNA-seq references. Please follow the instructions below.

1. Import ``cellranger_create_reference``
==============================================

	Import *cellranger_create_reference* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose **github.com/kalarman-cell-observatory/cumulus/Cellranger_create_reference** to import.

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cellranger_create_reference* workflow in the drop-down menu.

2. Upload requred data to Google Bucket
=======================================

	Required data may include input sample sheet, genome FASTA files and gene annotation GTF files.

3. Input sample sheet
=====================

	If multiple species are specified, a sample sheet in CSV format is required. We describe the sample sheet format below, with required columns highlighted in bold:

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Genome**
		  - Genome name
		* - **Fasta**
		  - Location to the genome assembly in FASTA/FASTA.gz format
		* - **Genes**
		  - Location to the gene annotation file in GTF/GTF.gz format
		* - Attributes
		  - Optional, A list of ``key:value`` pairs separated by ``;``. If set, ``cellranger mkgtf`` will be called to filter the user-provided GTF file. See `10x filter with mkgtf`_ for more details

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	See below for an example for building
	Example::

		Genome,Fasta,Genes,Attributes
		GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/GRCh38.fa.gz,gs://fc-e0000000-0000-0000-0000-000000000000/GRCh38.gtf.gz,gene_biotype:protein_coding;gene_biotype:lincRNA;gene_biotype:antisense
		mm10,gs://fc-e0000000-0000-0000-0000-000000000000/mm10.fa.gz,gs://fc-e0000000-0000-0000-0000-000000000000/mm10.gtf.gz

	If multiple species are specified, the reference will built under **Genome** names concatenated by '_and_'s. In the above example, the reference is stored under 'GRCh38_and_mm10'.

4. Workflow input
=================

	Required inputs are highlighted in bold. Note that **input_sample_sheet** and **input_fasta**, **input_gtf** , **genome** and attributes are mutually exclusive.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **input_sample_sheet**
		  - A sample sheet in CSV format allows users to specify more than 1 genomes to build references (e.g. human and mouse). If a sample sheet is provided, **input_fasta**, **input_gtf**, and attributes will be ignored.
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/input_sample_sheet.csv"
		  -
		* - **input_fasta**
		  - Input genome reference in either FASTA or FASTA.gz format
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
		  -
		* - **input_gtf**
		  - Input gene annotation file in either GTF or GTF.gz format
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf.gz"
		  -
		* - **genome**
		  - Genome reference name. New reference will be stored in a folder named **genome**
		  - refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_reference"
		  -
		* - attributes
		  - A list of ``key:value`` pairs separated by ``;``. If this option is not None, ``cellranger mkgtf`` will be called to filter the user-provided GTF file. See `10x filter with mkgtf`_ for more details
		  - "gene_biotype:protein_coding;gene_biotype:lincRNA;gene_biotype:antisense"
		  -
		* - pre_mrna
		  - If we want to build pre-mRNA references, in which we use full length transcripts as exons in the annotation file. We follow `10x build Cell Ranger compatible pre-mRNA Reference Package`_ to build pre-mRNA references
		  - true
		  - false
		* - ref_version
		  - reference version string
		  - Ensembl v94
		  -
		* - cellranger_version
		  - cellranger version, could be 6.1.1, 6.0.2, 6.0.1, 6.0.0, 5.0.1, 5.0.0, 4.0.0, 3.1.0, 3.0.2, or 2.2.0
		  - "6.1.1"
		  - "6.1.1"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "quay.io/cumulus" for images on Red Hat registry;

		  	- "cumulusprod" for backup images on Docker Hub.
		  - "quay.io/cumulus"
		  - "quay.io/cumulus"
		* - zones
		  - Google cloud zones
		  - "us-central1-a us-west1-a"
		  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - num_cpu
		  - Number of cpus to request for one node for building indices
		  - 1
		  - 1
		* - memory
		  - Memory size string for cellranger mkref
		  - "32G"
		  - "32G"
		* - disk_space
		  - Optional disk space in GB
		  - 100
		  - 100
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
		* - awsMaxRetries
		  - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
		  - 5
		  - 5

5. Workflow output
==================

	.. list-table::
		:widths: 2 2 10
		:header-rows: 1

		* - Name
		  - Type
		  - Description
		* - output_reference
		  - File
		  - Gzipped reference folder with name *genome.tar.gz*. We will also store a copy of the gzipped tarball under **output_directory** specified in the input.

---------------------------------

Build references for scATAC-seq
+++++++++++++++++++++++++++++++

We provide a wrapper of ``cellranger-atac mkref`` to build scATAC-seq references. Please follow the instructions below.

1. Import ``cellranger_atac_create_reference``
==============================================

	Import *cellranger_atac_create_reference* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose **github.com/klarman-cell-observatory/cumulus/Cellranger_atac_create_reference** to import.

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cellranger_atac_create_reference* workflow in the drop-down menu.

2. Upload required data to Google Bucket
===========================================

	Required data include config JSON file, genome FASTA file, gene annotation file (GTF or GFF3 format) and motif input file (JASPAR format).

3. Workflow input
=================

	Required inputs are highlighted in bold.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **genome**
		  - Genome reference name. New reference will be stored in a folder named **genome**
		  - refdata-cellranger-atac-mm10-1.1.0
		  -
		* - **input_fasta**
		  - GSURL for input fasta file
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/GRCh38.fa"
		  -
		* - **input_gtf**
		  - GSURL for input GTF file
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/annotation.gtf"
		  -
		* - organism
		  - Name of the organism
		  - "human"
		  -
		* - non_nuclear_contigs
		  - A comma separated list of names of contigs that are not in nucleus
		  - "chrM"
		  - "chrM"
		* - input_motifs
		  - Optional file containing transcription factor motifs in JASPAR format
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/motifs.pfm"
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_atac_reference"
		  -
		* - cellranger_atac_version
		  - cellranger-atac version, could be: 2.0.0, 1.2.0, 1.1.0
		  - "2.0.0"
		  - "2.0.0"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "quay.io/cumulus" for images on Red Hat registry;

		  	- "cumulusprod" for backup images on Docker Hub.
		  - "quay.io/cumulus"
		  - "quay.io/cumulus"
		* - zones
		  - Google cloud zones
		  - "us-central1-a us-west1-a"
		  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - memory
		  - Memory size string for cellranger-atac mkref
		  - "32G"
		  - "32G"
		* - disk_space
		  - Optional disk space in GB
		  - 100
		  - 100
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
		* - awsMaxRetries
		  - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
		  - 5
		  - 5

4. Workflow output
==================

	.. list-table::
		:widths: 2 2 10
		:header-rows: 1

		* - Name
		  - Type
		  - Description
		* - output_reference
		  - File
		  - Gzipped reference folder with name *genome.tar.gz*. We will also store a copy of the gzipped tarball under **output_directory** specified in the input.

---------------------------------

Build references for single-cell immune profiling data
++++++++++++++++++++++++++++++++++++++++++++++++++++++

We provide a wrapper of ``cellranger mkvdjref`` to build single-cell immune profiling references. Please follow the instructions below.

1. Import ``cellranger_vdj_create_reference``
==============================================

	Import *cellranger_vdj_create_reference* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose **github.com/klarman-cell-observatory/cumulus/Cellranger_vdj_create_reference** to import.

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cellranger_vdj_create_reference* workflow in the drop-down menu.

2. Upload requred data to Google Bucket
=======================================

	Required data include genome FASTA file and gene annotation file (GTF format).

3. Workflow input
=================

	Required inputs are highlighted in bold.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **input_fasta**
		  - Input genome reference in either FASTA or FASTA.gz format
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/Homo_sapiens.GRCh38.dna.toplevel.fa.gz"
		  -
		* - **input_gtf**
		  - Input gene annotation file in either GTF or GTF.gz format
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf.gz"
		  -
		* - **genome**
		  - Genome reference name. New reference will be stored in a folder named **genome**
		  - refdata-cellranger-vdj-GRCh38-alts-ensembl-3.1.0
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_vdj_reference"
		  -
		* - ref_version
		  - reference version string
		  - Ensembl v94
		  -
		* - cellranger_version
		  - cellranger version, could be 6.1.1, 6.0.2, 6.0.1, 6.0.0, 5.0.1, 5.0.0, 4.0.0, 3.1.0, 3.0.2, or 2.2.0
		  - "6.1.1"
		  - "6.1.1"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "quay.io/cumulus" for images on Red Hat registry;

		  	- "cumulusprod" for backup images on Docker Hub.
		  - "quay.io/cumulus"
		  - "quay.io/cumulus"
		* - zones
		  - Google cloud zones
		  - "us-central1-a us-west1-a"
		  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - memory
		  - Memory size string for cellranger mkvdjref
		  - "32G"
		  - "32G"
		* - disk_space
		  - Optional disk space in GB
		  - 100
		  - 100
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
		* - awsMaxRetries
		  - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
		  - 5
		  - 5

4. Workflow output
==================

	.. list-table::
		:widths: 2 2 10
		:header-rows: 1

		* - Name
		  - Type
		  - Description
		* - output_reference
		  - File
		  - Gzipped reference folder with name *genome.tar.gz*. We will also store a copy of the gzipped tarball under **output_directory** specified in the input.


.. _Import workflows to Terra: ../cumulus_import.html
.. _10x filter with mkgtf: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf
.. _10x build Cell Ranger compatible pre-mRNA Reference Package: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#premrna
