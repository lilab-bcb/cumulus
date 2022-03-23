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
	  - "shareseq_atac"
      -
	  -
	* - input_fastqs_directories
	  - | A comma-separated list of input FASTQs directories (urls).
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
	  - R1 Fastq pattern PE data.
	  - "_S*_L*_R1_001.fastq.gz"
	  - "_S*_L*_R1_001.fastq.gz"
	* - pe_read2_fastq_pattern
	  - R2 Fastq pattern PE data.
	  - "_S*_L*_R3_001.fastq.gz"
	  - "_S*_L*_R3_001.fastq.gz"
	* - read1_range
      - start, end(-1 for length-1) in -1/-u files for genomic sequence
	  - 0,-1
	  - 0,-1
	* - read2_range
      - | start, end(-1 for length-1) in -2 files for genomic sequence
	    | Only use in case of paired-end data
	  - 0,-1
	  - 0,-1	  
	* - barcode_fastq_pattern
	  - barcode index Fastq file identifier.
	  - "_S*_L*_R2_001.fastq.gz"
	  - "_S*_L*_R2_001.fastq.gz"
        * - read_format
          - [chromap option \-\-read-format] Format for read files and barcode files
          - "r1:0:-1,bc:0:-1"
          - "r1:0:-1,bc:0:-1"
        * - chromap_version
	  - Chromap version to use. Available options: 0.1.4, 0.1.5.
	  - "0.1.5"
          - "0.1.5"
	* - split_alignment
	  - [chromap option \-\-split-alignment] Allow split alignments. This option should be set only when mapping Hi-C reads.
	  - False
          -
	* - max_edit_dist_e
	  - [chromap option -e] Max edit distance allowed to map a read.
	  - 8
          - 8
	* - min_num_minimizer_s
	  - [chromap option -s] Min number of minimizers required to map a read.
	  - 2
          - 2
	* - ignore_minimizer_times_f
	  - [chromap option -f] Skip minimizers occuring > INT1 [500] times. INT2 [1000] is the threshold for a second round of seeding.
	  - "500,1000"
          - "500,1000"
	* - max_insert_size_l
	  - [chromap option -l] Max insert size, only for paired-end read mapping.
	  - 1000
          - 1000
	* - min_mapq_q
	  - [chromap option -q] Min MAPQ in range [0, 60] for mappings to be output.
	  - 30
          - 30
	* - min_read_length
	  - [chromap option \-\-min-read-length] Skip mapping the reads of length less than Min read length.
	  - 30
          - 30
	* - trim_adaptors
	  - | [chromap option \-\-trim-adapters]
            | Try to trim adapters on 3’. This only works for paired-end reads.
            | When the fragment length indicated by the read pair is less than the length of the reads,
            | the two mates are overlapped with each other. Then the regions outside the overlap are regarded as adapters and trimmed.
	  - True
          -
	* - remove_pcr_duplicates
	  - | [chromap option \-\-remove-pcr-duplicates]
            | Remove PCR duplicates.
	  - True
          -
	* - remove_pcr_duplicates_at_bulk_level
	  - | [chromap option \-\-remove-pcr-duplicates-at-bulk-level]
            | Remove PCR duplicates at bulk level for single cell data.
	  - False
          -
	* - remove_pcr_duplicates_at_cell_level
	  - | [chromap option \-\-remove-pcr-duplicates-at-cell-level]
            | Remove PCR duplicates at cell level for single cell data.
	  - False
          -
	* - tn5_shift
          - | [chromap option \-\-Tn5-shift]
	    | Perform Tn5 shift. When this option is turned on,
            | the forward mapping start positions are increased by 4bp and the reverse
            | mapping end positions are decreased by 5bp. Note that this works only when --SAM is NOT set.
	  - True
          -
	* - low_mem
          - | [chromap option \-\-low-mem]
	    | Use low memory mode. When this option is set,
            | multiple temporary intermediate mapping files might be
            | generated on disk and they are merged at the end of processing to reduce memory usage.
            | When this is NOT set, all the mapping results are kept in the memory before
            | they are saved on disk, which works more efficiently for datasets that are not too large.
	  - True
          -
	* - bc_error_threshold
          - | [chromap option \-\-bc-error-threshold]
	    | Max Hamming distance allowed to correct a barcode. Max allowed 2.
	  - 1
          - 1
	* - bc_probability_threshold
          - | [chromap option \-\-bc-probability-threshold]
	    | Min probability to correct a barcode.
	  - 0.9
          - 0.9
	* - output_mappings_not_in_whitelist
          - | [chromap option \-\-output-mappings-not-in-whitelist]
	    | Output mappings with barcode not in the whitelist.
	  -
          -
	* - output_format
	  - | Output format. The following formats are available:
            | bed, tagalign, sam, pairs
	  - "bed"
          -
	* - chr_order
          - | [chromap option \-\-chr-order]
	    | File with customized chromsome order.
	  -
          -
	* - pairs_natural_chr_order
	  - | [chromap option \-\-pairs-natural-chr-order]
            | File with natural chromosome order for pairs flipping.
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
	  - "64G"
	  - "64G"
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

See the table below for *chromap* workflow outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_aln_directory
	  - String
	  - Google Bucket URL of output directory. Within it, each folder is for one sample in the input sample sheet.

----------------------------

Prebuilt genome references
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We've built the following chromap references for users' convenience:

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38_and_mm10_chromap_v0.1.3**
		  - Human GRCh38 and Mouse mm10, comparable to cellranger reference GRCh38_and_mm10_atac_v1.2.0
		* - **GRCh38_chromap_v0.1.3**
		  - Mouse mm10, comparable to cellranger reference GRCh38-2020-A_arc_v2.0.0
		* - **mm10_chromap_v0.1.3**
		  - Human GRCh38, comparable to cellranger reference mm10-2020-A_arc_v2.0.0

.. _Chromap's manual: https://zhanghaowen.com/chromap/chromap.html
.. _genome reference: ./chromap.html#prebuilt-genome-references
