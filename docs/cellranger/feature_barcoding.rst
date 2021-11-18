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
		  - This is another alias for Single Cell 3' v3.
		* - **SC3Pv3**
		  - Single Cell 3′ v3.
		* - **SC3Pv2**
		  - Single Cell 3′ v2
		* - **fiveprime**
		  - Single Cell 5′
		* - **SC5P-PE**
		  - Single Cell 5′ paired-end (both R1 and R2 are used for alignment)
		* - **SC5P-R2**
		  - Single Cell 5′ R2-only (where only R2 is used for alignment)

#. *DataType* column.

	Put **citeseq** if CITE-seq, **hashing** if cell or nucleus hashing, **cmo** if CellPlex, **adt** if a mix of hashing and Cite-seq, and **crispr** here if Perturb-seq.

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
		* - scaffold_sequence
		  - Scaffold sequence in sgRNA for Purturb-seq, only used for crispr data type. If it is "", we assume guide barcode starts at position 0 of read 2
		  - "GTTTAAGAGCTAAGCTGGAA"
		  - ""
		* - max_mismatch
		  - Maximum hamming distance in feature barcodes for the adt task
		  - 3
		  - 3
		* - min_read_ratio
		  - Minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for the adt task and crispr data type
		  - 0.1
		  - 0.1
		* - cellranger_version
		  - cellranger version, could be 6.1.1, 6.0.2, 6.0.1, 6.0.0, 5.0.1, 5.0.0, 4.0.0, 3.1.0, 3.0.2, 2.2.0
		  - "6.1.1"
		  - "6.1.1"
		* - cumulus_feature_barcoding_version
		  - Cumulus_feature_barcoding version for extracting feature barcode matrix. Version available: 0.6.0, 0.5.0, 0.4.0, 0.3.0, 0.2.0.
		  - "0.6.0"
		  - "0.6.0"
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
		* - awsMaxRetries
		  - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
		  - 5
		  - 5

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
	* - output_fastqs_directory
	  - Array[String]
	  - A list of google bucket urls containing FASTQ files, one url per flowcell.
	* - output_count_directory
	  - Array[String]
	  - A list of google bucket urls containing feature-barcode count matrices, one url per sample.
	* - count_matrix
	  - String
	  - gs url for a template count_matrix.csv to run cumulus.

In addition, For each antibody tag or crispr tag sample, a folder with the sample ID is generated under ``output_directory``. In the folder, two files --- ``sample_id.csv`` and ``sample_id.stat.csv.gz`` --- are generated.

``sample_id.csv`` is the feature count matrix. It has the following format. The first line describes the column names: ``Antibody/CRISPR,cell_barcode_1,cell_barcode_2,...,cell_barcode_n``. The following lines describe UMI counts for each feature barcode, with the following format: ``feature_name,umi_count_1,umi_count_2,...,umi_count_n``.

``sample_id.stat.csv.gz`` stores the gzipped sufficient statistics. It has the following format. The first line describes the column names: ``Barcode,UMI,Feature,Count``. The following lines describe the read counts for every barcode-umi-feature combination.

If the feature barcode file has a third column, there will be two files for each feature type in the third column. For example, if ``hashing`` presents, ``sample_id.hashing.csv`` and ``sample_id.hashing.stat.csv.gz`` will be generated.

If data type is ``crispr``, three additional files, ``sample_id.umi_count.pdf``, ``sample_id.filt.csv`` and ``sample_id.filt.stat.csv.gz``, are generated.

``sample_id.umi_count.pdf`` plots number of UMIs against UMI with certain number of reads and colors UMIs with high likelihood of being chimeric in blue and other UMIs in red. This plot is generated purely based on number of reads each UMI has.

``sample_id.filt.csv`` is the filtered feature count matrix. It has the same format as ``sample_id.csv``.

``sample_id.filt.stat.csv.gz`` is the filtered sufficient statistics. It has the same format as ``sample_id.stat.csv.gz``.


.. _10x genomics v2 cell barcode white list: gs://regev-lab/resources/cellranger/737K-august-2016.txt.gz
.. _10x genomics v3 cell barcode white list: gs://regev-lab/resources/cellranger/3M-february-2018.txt.gz
.. _10x single cell RNA-seq sample index set names: https://support.10xgenomics.com/single-cell-gene-expression/index/doc/specifications-sample-index-sets-for-single-cell-3
