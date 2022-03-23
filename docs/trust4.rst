Run TRUST4 for immune repertoire reconstruction from RNA-seq data
----------------------------------------------------------------------

This ``TRUST4`` workflow extracts reads, assembles reads into immune receptor sequences and performs annotation.

----------------------------

Workflow inputs
^^^^^^^^^^^^^^^^^^

Below are inputs for *TRUST4* workflow. Notice that required inputs are in bold. Options described with [trust4 option -...] refer directly to TRUST4 options and are mentioned in `TRUST4's manual`_.

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - trust4_version
      - TRUST4 version
      - "master"
      - "master"
    * - **sample_id**
      - Sample Id.
      - "trust4_vdj"
      -
    * - **input_fastqs_directories**
      - | A comma-separated list of input FASTQs directories (urls).
	| Provide input fastqs or BAM but not both
      - "/path/fastq_dir1,/path/fastq_dir2"
      -
    * - **output_directory**
      - GS URL of output directory.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/trust4_result"
      -	  	  
    * - **genome**
      - Genome reference. It can be either of the following two formats:

		- String. Pre-built `genome reference`_.

		- Google bucket URL of a custom reference, must be a ``.tar.gz`` file.
      - | "human_trust4",
	| or "gs://user-bucket/trust4.tar.gz"
      -
    * - **acronym_file**
      - | The link/path of an index file in TSV format for fetching preset genome references by their names.
	| Set an GS URI if *backend* is ``gcp``; an S3 URI for ``aws`` backend; an absolute file path for ``local`` backend.
      - "s3://xxxx/index.tsv" or "gs://xxxx/index.tsv"
      -
    * - se_fastq_pattern	
      - R2 Fastq pattern SE data.
      - "_S*_L*_R2_001.fastq.gz"
      - "_S*_L*_R2_001.fastq.gz"
    * - pe_read1_fastq_pattern
      - R1 Fastq pattern for PE data.
      - "_S*_L*_R1_001.fastq.gz"
      - 
    * - pe_read2_fastq_pattern
      - R2 Fastq pattern for PE data.
      - "_S*_L*_R2_001.fastq.gz"
      - 
    * - read1_range
      - start, end(-1 for length-1) in -1/-u files for genomic sequence
      - "0,-1"
      - "0,-1"
    * - read2_range
      - | start, end(-1 for length-1) in -2 files for genomic sequence
	| Only use in case of paired-end data
      - "0,-1"
      - "0,-1"
    * - barcode_fastq_pattern
      - Barcode Fastq pattern.
      - "_S*_L*_R1_001.fastq.gz"
      - "_S*_L*_R1_001.fastq.gz"
    * - barcode_range
      - start, end(-1 for lenght-1), strand in a barcode is the true barcode
      - "0,15,+"
      - "0,15,+"
    * - barcode_whitelist
      - | path to the barcode whitelist
	| Cell barcode whitelist file. This is supposed to be a txt file where each line is a whitelisted barcode.
      - 
      -
    * - umi_fastq_pattern
      - UMI Fastq pattern
      - "_S*_L*_R1_001.fastq.gz"
      - "_S*_L*_R1_001.fastq.gz"
    * - umi_range
      - start, end(-1 for length-1), strand in a UMI is the true UMI 
      - "16,-1,+"
      - "16,-1,+"
    * - umi_bam_field
      - If BAM file is provided as input; provide bam field for UMI
      - 
      - 
    * - **input_bam**
      - | [trust4 option -b] Path to bam file
        | Provide input fastqs or BAM but not both
      - 
      - 
    * - bam_barcode_field
      - | [trust4 option --barcode] BAM field for barcode
	| Only use when BAM file is used as input
      - 
      - 
    * - bam_abnormal_unmap_flag
      - | [trust4 option --abnormalUnmapFlag] The flag in BAM for the unmapped read-pair is nonconcordant.
	| Only use when BAM file is used as input
      - 
      - 
    * - skipMateExtension
      - [trust4 option --skipMateExtension] Do not extend assemblies with mate information, useful for SMART-seq
      -
      - 
    * - mateIdSuffixLen
      - [trust4 option --mateIdSuffixLen] The suffix length in read id for mate
      - 
      -
    * - noExtraction
      - [trust4 option --noExtraction] Directly use the files from provided -1 -2/-u to assemble
      - 
      -
    * - repseq
      - [trust4 option --repseq] The data is from TCR-seq or BCR-seq
      - 
      -
    * - outputReadAssignment
      - [trust4 option --outputReadAssignment] Output read assignment results to the prefix_assign.out file 
      - 
      -
    * - docker_registry
      - Docker registry to use:

	  	- "quay.io/cumulus" for images on Red Hat registry;

		- "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - zones
      - Google cloud zones to consider for execution.
      - "us-east1-d us-west1-a us-west1-b"
      - "us-central1-b"
    * - num_cpu
      - Number of CPUs to request for mapping, setting trust4 option -t.
      - 8
      - 8
    * - memory
      - Memory size string for count per sample.
      - "32G"
      - "32G"
    * - disk_space
      - Disk space in GB needed for count per sample.
      - 200
      - 200
    * - backend
      - Cloud infrastructure backend to use. Available options:

	    - "gcp" for Google Cloud;
	    - "aws" for Amazon AWS;
	    - "local" for local machine.
      - "gcp"
      - "gcp"
    * - preemptible
      - Number of maximum preemptible tries allowed. This works only when *backend* is ``gcp``.
      - 2
      - 2
    * - awsMaxRetries
      - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
      - 5
      - 5

Workflow outputs
^^^^^^^^^^^^^^^^^^^

See the table below for *trust4* workflow outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_vdj_directory
	  - String
	  - Google Bucket/S3 URI of output directory.

----------------------------

Prebuilt genome references
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We've built the following TRUST4 references for users' convenience:

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **trust4-human**
		  - GRCh38 reference
		* - **trsut4-mouse**
		  - mm10 reference

---------------------------

Build TRUST4 References
^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide a wrapper of TRUST4 to build custom V,J,C gene database. Please follow the instructions below.

1. Workflow input
+++++++++++++++++++

Required inputs are highlighted **in bold**.

.. list-table::
    :widths: 5 20 10 5
    :header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - **reference_fasta**
      - Input genome reference in FASTA format.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/mm-10/genome.fa"
      -
    * - **annotation_gtf**
      - Input gene annotation file in GTF format.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/mm-10/genes.gtf"
      -
    * - **gene_name_list**
      - Gene name list of interest
      - 
      -
    * - **species**
      - | Species name
	| The available species name can be found on `IMGT FTP`_.
      - "Homo sapien"
      -
    * - **ref_name**
      - Reference name
      - "trust4-human"
      -
    * - **output_directory**
      - Cloud bucket URI of the output directory.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/trust4-reference"
      -
    * - docker_registry
      - Docker registry to use:

        - ``quay.io/cumulus`` for images on Red Hat registry;

        - ``cumulusprod`` for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - trust4_version
      - TRUST4 version to use. Currently support: ``master``.
      - "master"
      - "master"
    * - memory
      - Memory size string for count per sample.
      - "8G"
      - "8G"
    * - disk_space
      - Disk space in GB needed for count per sample.
      - 50
      - 50
    * - zones
      - Google cloud zones to consider for execution.
      - "us-east1-d us-west1-a us-west1-b"
      - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    * - backend
      - Cloud infrastructure backend to use. Available options:

        - ``gcp`` for Google Cloud;
        - ``aws`` for Amazon AWS;
        - ``local`` for local machine.
      - "gcp"
      - "gcp"
    * - preemptible
      - Number of maximum preemptible tries allowed. This works only when *backend* is ``gcp``.
      - 2
      - 2
    * - awsMaxRetries
      - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
      - 5
      - 5

2. Workflow Output
+++++++++++++++++++

.. list-table::
    :widths: 2 2 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - output_reference
      - File
      - Gzipped reference folder with name **"<ref_name>.tar.gz"**, where *<ref_name>* is specified by workflow input **ref_name** above. The workflow will save a copy of it under **output_directory** specified in workflow input above.


.. _IMGT FTP: https://www.imgt.org/download/V-QUEST/IMGT_V-QUEST_reference_directory/
.. _TRUST4's manual: https://github.com/liulab-dfci/TRUST4#trust4
.. _genome reference: ./trust4.html#prebuilt-genome-references
