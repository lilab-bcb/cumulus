``cellranger_workflow`` can extract feature-barcode count matrices in CSV format for feature barcoding assays such as *cell and nucleus hashing*, *CellPlex*, *CITE-seq*, and *Perturb-seq*.
For cell and nucleus hashing as well as CITE-seq, the feature refers to antibody. For Perturb-seq, the feature refers to guide RNA. Please follow the instructions below to configure ``cellranger_workflow``.

Tthe workflow uses `Cumulus Feature Barcoding`_ to process antibody and Perturb-Seq data.

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

		gcloud storage cp antibody_index.csv gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv


Sample sheet
++++++++++++

#. *Reference* column.

	This column is not used for extracting feature-barcode count matrix. To be consistent, you can put the reference for the associated scRNA-seq assay here.

#. *Chemistry* column.

	The following keywords are accepted for *Chemistry* column:

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Chemistry
		  - Explanation
		* - **auto**
		  - Default. Auto-detect the chemistry of your data from all possible 10x assay types.
		* - **threeprime**
		  - Auto-detect the chemistry of your data from all 3' assay types.
		* - **fiveprime**
		  - Auto-detect the chemistry of your data from all 5' assay types.
		* - **SC3Pv4**
		  - | Single Cell 3' v4. The workflow will auto-detect if Poly-A or CS1 capture method was applied to your data.
		    | **Notice:** This is a GEM-X chemistry, and only works for Cell Ranger v8.0.0+
		* - **SC3Pv3**
		  - Single Cell 3′ v3. This is a Next GEM chemistry. The workflow will auto-detect if Poly-A or CS1 capture method was applied to your data.
		* - **SC3Pv2**
		  - Single Cell 3′ v2
		* - **SC5Pv3**
		  - Single Cell 5' v3. **Notice:** This is a GEM-X chemistry, and only works for Cell Rangrer v8.0.0+
		* - **SC5Pv2**
		  - Single Cell 5′ v2
		* - **multiome**
		  - 10x Multiome barcodes

.. note::
	Not all 10x chemistry names are supported for feature barcoding, as the workflow uses *Cumulus Feature Barcoding* to process the data.

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

#. *AuxFile* column.

	Put cloud URI of the feature barcode file here.

Below is an example sample sheet::

	Sample,Reference,Flowcell,Chemistry,DataType,AuxFile
	sample_1_rna,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,auto,rna,
	sample_1_adt,,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,threeprime,hashing,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv
	sample_2_gex,GRCh38-2024-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,auto,rna
	sample_2_adt,GRCh38-2024-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,SC3Pv3,adt,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index2.csv
	sample_3_crispr,,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,fiveprime,crispr,gs://fc-e0000000-0000-0000-0000-000000000000/crispr_index.csv

In the sample sheet above, despite the header row,

	- Row 1 and 2 specify the GEX and Hashing libraries of the same sample.

	- Row 3 and 4 specify a sample which has GEX and **adt** (contains both Hashing and CITE-Seq data) libraries.

	- Row 5 describes one gRNA guide data for Perturb-seq (see ``crispr`` in *DataType* field).


Workflow input
++++++++++++++

For feature barcoding data, ``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files), and runs ``cumulus adt``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **input_csv_file**
		  - Sample Sheet (contains Sample, Reference, Flowcell, Chemistry, DataType, and AuxFile)
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
		  -
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
		* - cumulus_feature_barcoding_version
		  - Cumulus_feature_barcoding version for extracting feature barcode matrix.
		  - "1.0.0"
		  - "1.0.0"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "quay.io/cumulus" for images on Red Hat registry;

		  	- "cumulusprod" for backup images on Docker Hub.
		  - "quay.io/cumulus"
		  - "quay.io/cumulus"
		* - zones
		  - Google cloud zones. For GCP Batch backend, the zones are automatically restricted by the Batch settings.
		  - "us-central1-a us-west1-a"
		  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - feature_num_cpu
		  - Number of cpus for extracting feature count matrix
		  - 4
		  - 4
		* - feature_memory
		  - Optional memory string for extracting feature count matrix
		  - "32G"
		  - "32G"
		* - feature_disk_space
		  - Disk space in GB needed for extracting feature count matrix
		  - 100
		  - 100
		* - preemptible
		  - Number of preemptible tries. Only works for GCP
		  - 2
		  - 2
		* - awsQueueArn
		  - The AWS ARN string of the job queue to be used. Only works for AWS
		  - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
		  - ""

Parameters used for feature count matrix extraction
+++++++++++++++++++++++++++++++++++++++++++++++++++

Cell barcode inclusion lists (previously known as whitelists) are automatically decided based on the *Chemistry* specified in the sample sheet. The association table is `here <https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist>`_.

Cell barcode matching settings are also automatically decided based on the chemistry specified:

	* For 10x V3 and V4 chemistry: a hamming distance of ``0`` is allowed for matching cell barcodes, and the UMI length is ``12``;
	* For *multiome*: a hamming distance of ``1`` is allowed for matching cell barcodes, and the UMI length is ``12``;
	* For 10x V2 chemistry: a hamming distance of ``1`` is allowed for matching cell barcodes, and the UMI length is ``10``.

For Perturb-seq data, a small number of sgRNA protospace sequences will be sequenced ultra-deeply and we may have PCR chimeric reads. Therefore, we generate filtered feature count matrices as well in a data driven manner:

#. First, plot the histogram of UMIs with certain number of read counts. The number of UMIs with ``x`` supporting reads decreases when ``x`` increases. We start from ``x = 1``, and a valley between two peaks is detected if we find ``count[x] < count[x + 1] < count[x + 2]``. We filter out all UMIs with ``< x`` supporting reads since they are likely formed due to chimeric reads.

#. In addition, we also filter out barcode-feature-UMI combinations that have their read count ratio, which is defined as total reads supporting barcode-feature-UMI over total reads supporting barcode-UMI, no larger than ``min_read_ratio`` parameter set above.

Workflow outputs
++++++++++++++++

The table below lists important feature barcoding output when using Cumulus Feature Barcoding:

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - cumulus_adt.output_count_directory
	  - Array[String]
	  - Subworkflow output. A list of cloud URIs containing feature-barcode count matrices, one URI per sample.

In addition, For each antibody tag or crispr tag sample, a folder with the sample ID is generated under ``output_directory``. In the folder, two files --- ``sample_id.csv`` and ``sample_id.stat.csv.gz`` --- are generated.

``sample_id.csv`` is the feature count matrix. It has the following format. The first line describes the column names: ``Antibody/CRISPR,cell_barcode_1,cell_barcode_2,...,cell_barcode_n``. The following lines describe UMI counts for each feature barcode, with the following format: ``feature_name,umi_count_1,umi_count_2,...,umi_count_n``.

``sample_id.stat.csv.gz`` stores the gzipped sufficient statistics. It has the following format. The first line describes the column names: ``Barcode,UMI,Feature,Count``. The following lines describe the read counts for every barcode-umi-feature combination.

If the feature barcode file has a third column, there will be two files for each feature type in the third column. For example, if ``hashing`` presents, ``sample_id.hashing.csv`` and ``sample_id.hashing.stat.csv.gz`` will be generated.

``sample_id.report.txt`` is a summary report in TXT format. The first lines describe the total number of reads parsed, the number of reads with valid cell barcodes (and percentage over all parsed reads), the number of reads with valid feature barcodes (and percentage over all parsed reads) and the number of reads with both valid cell and feature barcodes (and percentage over all parsed reads). It is then followed by sections describing each feature type. In each section, 7 lines are shown: section title, number of valid cell barcodes (with matching cell barcode and feature barcode) in this section, number of reads for these cell barcodes, mean number of reads per cell barcode, number of UMIs for these cell barcodes, mean number of UMIs per cell barcode and sequencing saturation.

If data type is ``crispr``, three additional files, ``sample_id.umi_count.pdf``, ``sample_id.filt.csv`` and ``sample_id.filt.stat.csv.gz``, are generated.

``sample_id.umi_count.pdf`` plots number of UMIs against UMI with certain number of reads and colors UMIs with high likelihood of being chimeric in blue and other UMIs in red. This plot is generated purely based on number of reads each UMI has. For better visualization, we do not show UMIs with > 50 read counts (rare in data).

``sample_id.filt.csv`` is the filtered feature count matrix. It has the same format as ``sample_id.csv``.

``sample_id.filt.stat.csv.gz`` is the filtered sufficient statistics. It has the same format as ``sample_id.stat.csv.gz``.


.. _Cumulus Feature Barcoding: https://github.com/lilab-bcb/cumulus_feature_barcoding
