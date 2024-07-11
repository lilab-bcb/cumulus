``cellranger_workflow`` can extract feature-barcode count matrices in CSV format for feature barcoding assays such as *cell and nucleus hashing*, *CellPlex*, *CITE-seq*, and *Perturb-seq*. For cell and nucleus hashing as well as CITE-seq, the feature refers to antibody. For Perturb-seq, the feature refers to guide RNA. Please follow the instructions below to configure ``cellranger_workflow``.

Prepare feature barcode files
+++++++++++++++++++++++++++++

	Prepare a CSV file with the following format: feature_barcode,feature_name.
	See below for an example::

		TTCCTGCCATTACTA,sample_1
		CCGTACCTCATTGTT,sample_2
		GGTAGATGTCCTCAG,sample_3
		TGGTGTCATTCTTGA,sample_4

	The above file describes a cell hashing application with 4 samples.

	If cell hashing and CITE-seq data share a same sample index, you should concatenate hashing and CITE-seq barcodes together and add a third column indicating the feature type.
	See below for an example::

		TTCCTGCCATTACTA,sample_1,hashing
		CCGTACCTCATTGTT,sample_2,hashing
		GGTAGATGTCCTCAG,sample_3,hashing
		TGGTGTCATTCTTGA,sample_4,hashing
		CTCATTGTAACTCCT,CD3,citeseq
		GCGCAACTTGATGAT,CD8,citeseq

	Then upload it to your google bucket::

		gsutil antibody_index.csv gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv


Sample sheet
++++++++++++

#. **Reference** column.

	This column is not used for extracting feature-barcode count matrix. To be consistent, please put the reference for the associated scRNA-seq assay here.

#. **Index** column.

	The ADT/HTO index can be either Illumina index primer sequence (e.g. ``ATTACTCG``, also known as ``D701``), or `10x single cell RNA-seq sample index set names`_ (e.g. SI-GA-A12).

	**Note 1**: All ADT/HTO index sequences (including 10x's) should have the same length (8 bases). If one index sequence is shorter (e.g. ATCACG), pad it with P7 sequence (e.g. ATCACGAT).

	**Note 2**: It is users' responsibility to avoid index collision between 10x genomics' RNA indexes (e.g. SI-GA-A8) and Illumina index sequences for used here (e.g. ``ATTACTCG``).

	**Note 3**: For NextSeq runs, please reverse complement the ADT/HTO index primer sequence (e.g. use reverse complement ``CGAGTAAT`` instead of ``ATTACTCG``).

#. *Chemistry* column.

	The following keywords are accepted for *Chemistry* column:

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Chemistry
		  - Explanation
		* - **auto**
		  - Default. This is an alias for Single Cell 3' v3 (SC3Pv3)
		* - **threeprime**
		  - This is another alias for Single Cell 3' v3 (SC3Pv3)
		* - **SC3Pv4**
		  - Single Cell 3' v4. **Notice:** This is GEM-X chemistry, and only works for Cell Ranger v8.0.0+
		* - **SC3Pv3**
		  - Single Cell 3′ v3
		* - **SC3Pv2**
		  - Single Cell 3′ v2
		* - **fiveprime**
		  - Single Cell 5′
		* - **SC5P-PE**
		  - Single Cell 5′ paired-end (both R1 and R2 are used for alignment)
		* - **SC5P-PE-v3**
		  - Single Cell 5' paired-end v3 (both R1 and R2 are used for alignment). **Notice:** This is GEM-X chemistry, and only works for Cell Ranger v8.0.0+
		* - **SC5P-R2**
		  - Single Cell 5′ R2-only (where only R2 is used for alignment)
		* - **SC5P-R2-v3**
		  - Single Cell 5' R2-only v3 (where only R2 is used for alignment). **Notice:** This is GEM-X chemistry, and only works for Cell Rangrer v8.0.0+
		* - **multiome**
		  - 10x Multiome barcodes

#. *DataType* column.

	The following keywords are accepted for *DataType* column:

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - DataType
		  - Explanation
		* - **citeseq**
		  - CITE-seq
		* - **hashing**
		  - Cell or nucleus hashing
		* - **cmo**
		  - CellPlex
		* - **adt**
		  - Hashing and CITE-seq are in the same library
		* - **crispr**
		  - | Perturb-seq/CROP-seq
		    | If neither *crispr_barcode_pos* nor *scaffold_sequence* (see Workflow input) is set, **crispr** refers to 10x CRISPR assays. If in addition *Chemistry* is set to be **SC3Pv3** or its aliases, Cumulus automatically complement the middle two bases to convert 10x feature barcoding cell barcodes back to 10x RNA cell barcodes.
		    | Otherwise, **crispr** refers to non 10x CRISPR assays, such as CROP-Seq. In this case, we assume feature barcoding cell barcodes are the same as the RNA cell barcodes and no cell barcode convertion will be conducted.

#. *FetureBarcodeFile* column.

	Put Google Bucket URL of the feature barcode file here.

#. Example::

	Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
	sample_1_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,rna
	sample_1_adt,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,ATTACTCG,SC3Pv3,adt,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv
	sample_2_adt,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,TCCGGAGA,SC3Pv3,adt,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv
	sample_3_crispr,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,CGCTCATT,SC3Pv3,crispr,gs://fc-e0000000-0000-0000-0000-000000000000/crispr_index.csv

In the sample sheet above, despite the header row,

	- First row describes the normal 3' RNA assay;

	- Second row describes its associated antibody tag data, which can from either a CITE-seq, cell hashing, or nucleus hashing experiment.

	- Third row describes another tag data, which is in 10x genomics' V3 chemistry. For tag and crispr data, it is important to explicitly state the chemistry (e.g. ``SC3Pv3``).

	- Last row describes one gRNA guide data for Perturb-seq (see ``crispr`` in *DataType* field).

Workflow input
++++++++++++++

For feature barcoding data, ``cellranger_workflow`` takes Illumina outputs as input and runs ``cellranger mkfastq`` and ``cumulus adt``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

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
		  -
		* - run_mkfastq
		  - If you want to run ``cellranger mkfastq``
		  - true
		  - true
		* - run_count
		  - If you want to run ``cumulus adt``
		  - true
		  - true
		* - delete_input_bcl_directory
		  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges
		  - false
		  - false
		* - mkfastq_barcode_mismatches
		  - Number of mismatches allowed in matching barcode indices (bcl2fastq2 default is 1)
		  - 0
		  -
		* - mkfastq_force_single_index
		  - If 10x-supplied i7/i5 paired indices are specified, but the flowcell was run with only one sample index, allow the demultiplex to proceed using the i7 half of the sample index pair
		  - false
		  - false
		* - mkfastq_filter_single_index
		  - Only demultiplex samples identified by an i7-only sample index, ignoring dual-indexed samples. Dual-indexed samples will not be demultiplexed
		  - false
		  - false
		* - mkfastq_use_bases_mask
		  - Override the read lengths as specified in *RunInfo.xml*
		  - "Y28n*,I8n*,N10,Y90n*"
		  -
		* - mkfastq_delete_undetermined
		  - Delete undetermined FASTQ files generated by bcl2fastq2
		  - true
		  - false
		* - crispr_barcode_pos
		  - Barcode start position at Read 2 (0-based coordinate) for CRISPR
		  - 19
		  - 0
		* - scaffold_sequence
		  - Scaffold sequence in sgRNA for Purturb-seq, only used for crispr data type.
		  - "GTTTAAGAGCTAAGCTGGAA"
		  - ""
		* - max_mismatch
		  - Maximum hamming distance in feature barcodes for the adt task (changed to 2 as default)
		  - 2
		  - 2
		* - min_read_ratio
		  - Minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for the adt task and crispr data type
		  - 0.1
		  - 0.1
		* - cellranger_version
		  - cellranger version, could be: 8.0.0, 7.2.0, 7.1.0, 7.0.1, 7.0.0, 6.1.2, 6.1.1, 6.0.2, 6.0.1, 6.0.0, 5.0.1, 5.0.0
		  - "8.0.0"
		  - "8.0.0"
		* - cumulus_feature_barcoding_version
		  - Cumulus_feature_barcoding version for extracting feature barcode matrix. Version available: 0.11.3, 0.11.2, 0.11.1, 0.11.0, 0.10.0, 0.9.0, 0.8.0, 0.7.0, 0.6.0, 0.5.0, 0.4.0, 0.3.0, 0.2.0.
		  - "0.11.3"
		  - "0.11.3"
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
		  - Number of cpus to request for one node for cellranger mkfastq
		  - 32
		  - 32
		* - memory
		  - Memory size string for cellranger mkfastq
		  - "120G"
		  - "120G"
		* - feature_num_cpu
		  - Number of cpus for extracting feature count matrix
		  - 4
		  - 4
		* - feature_memory
		  - Optional memory string for extracting feature count matrix
		  - "32G"
		  - "32G"
		* - mkfastq_disk_space
		  - Optional disk space in GB for mkfastq
		  - 1500
		  - 1500
		* - feature_disk_space
		  - Disk space in GB needed for extracting feature count matrix
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
		* - awsQueueArn
		  - The AWS ARN string of the job queue to be used. This only works for ``aws`` backend.
		  - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
		  - ""

Parameters used for feature count matrix extraction
+++++++++++++++++++++++++++++++++++++++++++++++++++

If the chemistry is V2, `10x genomics v2 cell barcode white list`_ will be used, a hamming distance of 1 is allowed for matching cell barcodes, and the UMI length is 10.
If the chemistry is V3, `10x genomics v3 cell barcode white list`_ will be used, a hamming distance of 0 is allowed for matching cell barcodes, and the UMI length is 12.

For Perturb-seq data, a small number of sgRNA protospace sequences will be sequenced ultra-deeply and we may have PCR chimeric reads. Therefore, we generate filtered feature count matrices as well in a data driven manner:

#. First, plot the histogram of UMIs with certain number of read counts. The number of UMIs with ``x`` supporting reads decreases when ``x`` increases. We start from ``x = 1``, and a valley between two peaks is detected if we find ``count[x] < count[x + 1] < count[x + 2]``. We filter out all UMIs with ``< x`` supporting reads since they are likely formed due to chimeric reads.

#. In addition, we also filter out barcode-feature-UMI combinations that have their read count ratio, which is defined as total reads supporting barcode-feature-UMI over total reads supporting barcode-UMI, no larger than ``min_read_ratio`` parameter set above.

Workflow outputs
++++++++++++++++

See the table below for important outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cellranger_mkfastq.output_fastqs_directory
	  - Array[String]?
	  - Subworkflow output. A list of cloud urls containing FASTQ files, one url per flowcell.
	* - cumulus_adt.output_count_directory
	  - Array[String]?
	  - Subworkflow output. A list of cloud urls containing feature-barcode count matrices, one url per sample.

In addition, For each antibody tag or crispr tag sample, a folder with the sample ID is generated under ``output_directory``. In the folder, two files --- ``sample_id.csv`` and ``sample_id.stat.csv.gz`` --- are generated.

``sample_id.csv`` is the feature count matrix. It has the following format. The first line describes the column names: ``Antibody/CRISPR,cell_barcode_1,cell_barcode_2,...,cell_barcode_n``. The following lines describe UMI counts for each feature barcode, with the following format: ``feature_name,umi_count_1,umi_count_2,...,umi_count_n``.

``sample_id.stat.csv.gz`` stores the gzipped sufficient statistics. It has the following format. The first line describes the column names: ``Barcode,UMI,Feature,Count``. The following lines describe the read counts for every barcode-umi-feature combination.

If the feature barcode file has a third column, there will be two files for each feature type in the third column. For example, if ``hashing`` presents, ``sample_id.hashing.csv`` and ``sample_id.hashing.stat.csv.gz`` will be generated.

``sample_id.report.txt`` is a summary report in TXT format. The first lines describe the total number of reads parsed, the number of reads with valid cell barcodes (and percentage over all parsed reads), the number of reads with valid feature barcodes (and percentage over all parsed reads) and the number of reads with both valid cell and feature barcodes (and percentage over all parsed reads). It is then followed by sections describing each feature type. In each section, 7 lines are shown: section title, number of valid cell barcodes (with matching cell barcode and feature barcode) in this section, number of reads for these cell barcodes, mean number of reads per cell barcode, number of UMIs for these cell barcodes, mean number of UMIs per cell barcode and sequencing saturation.

If data type is ``crispr``, three additional files, ``sample_id.umi_count.pdf``, ``sample_id.filt.csv`` and ``sample_id.filt.stat.csv.gz``, are generated.

``sample_id.umi_count.pdf`` plots number of UMIs against UMI with certain number of reads and colors UMIs with high likelihood of being chimeric in blue and other UMIs in red. This plot is generated purely based on number of reads each UMI has. For better visualization, we do not show UMIs with > 50 read counts (rare in data).

``sample_id.filt.csv`` is the filtered feature count matrix. It has the same format as ``sample_id.csv``.

``sample_id.filt.stat.csv.gz`` is the filtered sufficient statistics. It has the same format as ``sample_id.stat.csv.gz``.


.. _10x genomics v2 cell barcode white list: gs://regev-lab/resources/cellranger/737K-august-2016.txt.gz
.. _10x genomics v3 cell barcode white list: gs://regev-lab/resources/cellranger/3M-february-2018.txt.gz
.. _10x single cell RNA-seq sample index set names: https://support.10xgenomics.com/single-cell-gene-expression/index/doc/specifications-sample-index-sets-for-single-cell-3
