To process SHARE-Seq data, follow the specific instructions below.

Sample sheet
++++++++++++

#. **Sample** column.

#. *Reference* column.

	Pre-built Chromap references are summarized below.

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

	Pre-built STARsolo references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

                * - Keyword
                  - Description
                * - **GRCh38-2020-A**
                  - Human GRCh38, comparable to cellranger reference 2020-A (GENCODE v32/Ensembl 98)
                * - **mm10-2020-A**
                  - Mouse mm10, comparable to cellranger reference 2020-A (GENCODE vM23/Ensembl 98)
                * - **GRCh38-and-mm10-2020-A**
                  - Human GRCh38 (GENCODE v32/Ensembl 98) and mouse mm10 (GENCODE vM23/Ensembl 98)
                * - **GRCh38**
                  - Human GRCh38, comparable to cellranger reference 3.0.0, Ensembl v93 gene annotation
                * - **mm10**
                  - Mouse mm10, comparable to cellranger reference 3.0.0, Ensembl v93 gene annotation

#. **Flowcell** column.
        
        Put Flowcell URI/path here

#. **Lane** column.

        Put lane information here

#. **Index** column.

	Put index name here.

#. *Type* column.

	This column is optional when running only run_mkfastq. If you want to put a value, put **gex** or **atac** here.

#. Example::

	Sample,Reference,Flowcell,Lane,Index,Type
	sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,A1,gex
	sample_2,GRCh38_chromap_v0.1.3,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3,A2,atac

Workflow input
++++++++++++++

For SHARE-seq data, ``shareseq_workflow`` takes Illumina outputs as input and runs ``shareseq mkfastq``, ``Chromap`` and ``STARsolo``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **input_csv_file**
		  - Sample Sheet (contains Sample, Lane, Index, Flowcell, [optional for mkfastq: Reference, Type])
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/shareseq_output"
		  - Results are written under directory *output_directory* and will overwrite any existing files at this location.
		* - run_mkfastq
		  - If you want to run ``shareseq mkfastq``
		  - true
		  - true
		* - run_count
		  - If you want to run ``Chromap and STARsolo``
		  - true
		  - true
		* - delete_input_bcl_directory
		  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges
		  - false
		  - false
		* - dual_index_type
		  - Dual-index Paired-end flowcell workflow type
		  - "auto" or "forward" or "reverse"
		  - auto
		* - barcode_mismatches
		  - Number of allowed mismatches per index
		  - 0, 1 or 2
		  - 1
		* - use_bases_mask
		  - Override the read lengths as specified in RunInfo.xml
		  -
		  -
		* - read1_fastq_pattern
		  - read1 fastq pattern
		  - "_S*_L*_R1_001.fastq.gz"
		  - "_S*_L*_R1_001.fastq.gz"
		* - read2_fastq_pattern
		  - read2 fastq pattern
		  - "_S*_L*_R3_001.fastq.gz"
		  - "_S*_L*_R3_001.fastq.gz"
		* - index_fastq_pattern
		  - index fastq pattern
		  - "_S*_L*_R2_001.fastq.gz"
		  - "_S*_L*_R2_001.fastq.gz"
		* - docker_registry
		  - Docker registry for shareseq_workflow. Options:

                        - "quay.io/cumulus" for images on Red Hat registry;

                        - "cumulusprod" for backup images on Docker Hub.

		  - quay.io/cumulus
		  - quay.io/cumulus
		* - mkfastq_docker_registry
		  - docker registry for specifically mkfastq.
		  - gcr.io/broad-cumulus
		  - gcr.io/broad-cumulus
		* - zones
                  - Google cloud zones
                  - "us-central1-a us-west1-a"
                  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - backend
                  - Cloud backend for file transfer. Available options:

                    - "gcp" for Google Cloud;
                    - "aws" for Amazon AWS;
                    - "local" for local machine.
                  - "gcp"
                  - "gcp"
		* - acronym_file
		  - | The link/path of an index file in TSV format for fetching preset genome references, chemistry whitelists, etc. by their names.
		    | Set an GS URI if *backend* is ``gcp``; an S3 URI for ``aws`` backend; an absolute file path for ``local`` backend.
		  - "s3://xxxx/index.tsv"
		  - "gs://xxxx/ref-data/shareseq/index.tsv"
		* - config_version
		  - Config docker version used for processing sample sheets
		  - "0.2"
		  - "0.2"
		* - shareseq_mkfastq_version
		  - SHARE-seq mkfastq version
		  - "0.1.0"
		  - "0.1.0"
		* - shareseq_reorg_version
		  - SHARE-seq reorg version
		  - "0.1.0"
		  - "0.1.0"
		* - star_version
		  - STAR version
		  - "2.7.9a"
		  - "2.7.9a"
		* - chromap_version
		  - Chromap version 
		  - "0.1.4"
		  - "0.1.4"
		* - shareseq_mkfastq_num_cpu
		  - Number of CPUs for shareseq_mkfastq
		  - 32
		  - 32
		* - shareseq_mkfastq_memory
		  - Memory shareseq_mkfastq.
		  - "120G"
		  - "120G"
                * - sharseq_reorg_num_cpu
                  - Number of CPUs for shareseq_reorg
                  - 4
                  - 4
                * - sharseq_reorg_memory
                  - Memory sharseq_reorg.
                  - "8G"
                  - "8G"
                * - starsolo_num_cpu
                  - Number of CPUs for STARsolo.
                  - 32
                  - 32
                * - starsolo_memory
                  - Memory for STARsolo.
                  - "120G"
                  - "120G"
                * - chromap_num_cpu
                  - Number of CPUs for Chromap.
                  - 8
                  - 8
                * - chromap_memory
                  - Memory for Chromap.
                  - "64G"
                  - "64G"
                * - mkfastq_disk_space
                  - Disk space for shareseq_mkfastq
                  - 1500
                  - 1500
                * - shareseq_reorg_disk_space
                  - Disk space for shareseq_reorg
                  - 500
                  - 500
                * - starsolo_disk_space
                  - Disk space for STARsolo
                  - 500
                  - 500
                * - chromap_disk_space
                  - Disk space for Chromap
                  - 200
                  - 200
                * - preemptible
                  - Number of preemptible tries
                  - 2
                  - 2
                * - awsMaxRetries
                  - Max number of retries for AWS instance
                  - 5
                  - 5

Workflow output
+++++++++++++++

See the table below for important SHARE-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - demuxed_fastqs
	  - Array[String]
	  - 
            | A list of google bucket urls containing FASTQ files, one url per flowcell.
            | In SHARE-Seq: 

                           - R1 => read1 
                           - R2 => barcode
                           - R3 => read2

	* - reorg_gex_fastqs
	  - Array[String]
	  - 
            | A list of google bucket urls containing reorganized gene expression FASTQ files, one url per sample.
            | There are 2 FASTQs generated per sample: R1 (read1) and R2 (read2).
	* - reorg_atac_fastqs
	  - Array[String]
	  - 
            | A list of google bucket urls containing reorganized ATAC-seq FASTQ files, one url per sample.
            | There are 3 FASTQs generated per sample: I1 (index), R1 (read1) and R2 (read2).
        * - gex_outputs
	  - Array[String]
          -
            | Contains STARsolo (gene expression) output.
            | More detailed explaination on output files can be found in `STAR manual`_.  
        * - atac_outputs
	  - Array[String]
	  - 
            | Contains Chromap (ATAC-seq) output in BED format. 
            | Columns in BED file are:
                                 
                                  - chrom
                                  - chrom_start 
                                  - chrom_end
                                  - barcode
                                  - duplicate_count

	    | More detailed explaination on output files can be found in `Chromap README`_.

.. _STAR manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf
.. _Chromap README: https://github.com/haowenz/chromap/blob/master/README.md
