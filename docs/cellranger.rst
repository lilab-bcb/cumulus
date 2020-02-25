Run Cell Ranger tools using cellranger_workflow
-----------------------------------------------

``cellranger_workflow`` wraps Cell Ranger to process single-cell/nucleus RNA-seq, single-cell ATAC-seq and single-cell immune profiling data, and supports feature barcoding (cell/nucleus hashing, CITE-seq, Perturb-seq). It also provide routines to build cellranger references.

A general step-by-step instruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Import ``cellranger_workflow``
+++++++++++++++++++++++++++++++++

	Import *cellranger_workflow* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *cellranger_workflow* workflow is under ``Broad Methods Repository`` with name "**cumulus/cellranger_workflow**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cellranger_workflow* workflow in the drop-down menu.

2. Upload sequencing data to Google bucket
++++++++++++++++++++++++++++++++++++++++++

	Copy your sequencing output to your workspace bucket using gsutil_ (you already have it if you've installed Google cloud SDK) in your unix terminal.

	You can obtain your bucket URL in the dashboard tab of your Terra workspace under the information panel.

	.. image:: images/google_bucket_link.png
	
	Use ``gsutil cp [OPTION]... src_url dst_url`` to copy data to your workspace bucket. For example, the following command copies the directory at /foo/bar/nextseq/Data/VK18WBC6Z4 to a Google bucket::

		gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4
	
	``-m`` means copy in parallel, ``-r`` means copy the directory recursively, and ``gs://fc-e0000000-0000-0000-0000-000000000000`` should be replaced by your own workspace Google bucket URL.

	Note: Broad users need to be on an UGER node (not a login node) in order to use the ``-m`` flag

	Request an UGER node::

		reuse UGER
		qrsh -q interactive -l h_vmem=4g -pe smp 8 -binding linear:8 -P regevlab

	The above command requests an interactive node with 4G memory per thread and 8 threads. Feel free to change the memory, thread, and project parameters.

	Once you're connected to an UGER node, you can make gsutil_ available by running::

		reuse Google-Cloud-SDK

3. Prepare a sample sheet
+++++++++++++++++++++++++

	**3.1 Sample sheet format**:

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to demultiplex flowcells and generate channel-specific count matrices. Note that *Sample*, *Lane*, and *Index* columns are defined exactly the same as in 10x's simple CSV layout file.

	A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name.
		* - **Reference**
		  - 
		  	| Provides the reference genome used by Cell Ranger for each 10x channel. 
		  	| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as *GRCh38_v3.0.0*.
		  	| A full list of available keywords is included in each of the following data type sections (e.g. sc/snRNA-seq) below.
		* - **Flowcell**
		  - 
		    | Indicates the Google bucket URLs of uploaded BCL folders. 
		    | If starts with FASTQ files, this should be Google bucekt URLs of uploaded FASTQ folders.
		    | The FASTQ folders should contain one subfolder for each sample in the flowcell with the sample name as the subfolder name.
		    | Each subfolder contains FASTQ files for that sample. 
		* - **Lane**
		  - 
		    | Tells which lanes the sample was pooled into.
		    | Can be either single lane (e.g. 8) or a range (e.g. 7-8) or all (e.g. \*).
		* - **Index**
		  - Sample index (e.g. SI-GA-A12).
		* - Chemistry
		  - Describes the 10x chemistry used for the sample. This column is optional. 
		* - DataType
		  - 
			| Describes the data type of the sample --- *rna*, *vdj*, *adt*, or *crispr*. 
			| **rna** refers to gene expression data (*cellranger count*), 
			| **vdj** refers to V(D)J data (*cellranger vdj*), 
			| **adt** refers to antibody tag data, which can be either CITE-Seq, cell-hashing, or nucleus-hashing, 
			| **crispr** refers to Perturb-seq guide tag data,
			| **atac** refers to scATAC-Seq data (*cellranger-atac count*).
			| This column is optional and the default data type is *rna*.
		* - FeatureBarcodeFile
		  - 
		  	| Google bucket urls pointing to feature barcode files for *adt* and *crispr* data. 
		  	| Features can be either antibody for CITE-Seq, cell-hashing, nucleus-hashing or gRNA for Perburb-seq. 
		  	| This column is optional provided no *adt* or *crispr* data are in the sample sheet.

	The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 4 samples sequenced in two flowcells.

	Example::

		Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
		sample_1,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,rna
		sample_2,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,SC3Pv3,rna
		sample_3,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime,rna
		sample_4,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8,fiveprime,rna
		sample_1,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime,rna
		sample_2,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,SC3Pv3,rna
		sample_3,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime,rna
		sample_4,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8,fiveprime,rna

	**3.2 Upload your sample sheet to the workspace bucket:**

		Example::

			gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/

4. Launch analysis
++++++++++++++++++

	In your workspace, open ``cellranger_workflow`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Run workflow with inputs defined by file paths`` as below

		.. image:: images/single_workflow.png

	and click ``SAVE`` button. Select ``Use call caching`` and click ``INPUTS``. Then fill in appropriate values in the ``Attribute`` column. Alternative, you can upload a JSON file to configure input by clicking ``Drag or click to upload json``. 

	Once INPUTS are appropriated filled, click ``RUN ANALYSIS`` and then click ``LAUNCH``. 

5. Notice: run ``cellranger mkfastq`` if you are non Broad Institute users
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	Non Broad Institute users that wish to run ``cellranger mkfastq`` must create a custom docker image that contains ``bcl2fastq``.

		See :ref:`bcl2fastq-docker` instructions.

6. Do not run ``cellranger mkfastq``
++++++++++++++++++++++++++++++++++++

Sometimes, users might want to perform demultiplexing locally and only run the count part on the cloud. This section describes how to only run the count part via ``cellranger_workflow``.

#. Copy your FASTQ files to the workspace using gsutil_ in your unix terminal. 

	You should upload folders of FASTQ files. The uploaded folder (for one flowcell) should contain one subfolder for each sample belong to the this flowcell. In addition, the subfolder name should be the sample name. Each subfolder contains FASTQ files for that sample.

	Example::

		gsutil -m cp -r /foo/bar/fastq_path/K18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq

#. Create a sample sheet.

	**Flowcell** column should list Google bucket URLs of the FASTQ folders for flowcells.

	Example::

		Sample,Reference,Flowcell
		sample_1,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq

#. Set optional input ``run_mkfastq`` to ``false``.

---------------------------------

Single-cell and single-nucleus RNA-seq
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
		* - **GRCh38_v3.0.0**
		  - Human GRCh38, cellranger reference 3.0.0, Ensembl v93 gene annotation
		* - **hg19_v3.0.0**
		  - Human hg19, cellranger reference 3.0.0, Ensembl v87 gene annotation
		* - **mm10_v3.0.0**
		  - Mouse mm10, cellranger reference 3.0.0, Ensembl v93 gene annotation
		* - **GRCh38_and_mm10_v3.1.0**
		  - Human (GRCh38) and mouse (mm10), cellranger references 3.1.0, Ensembl v93 gene annotations for both human and mouse
		* - **GRCh38_v1.2.0** or **GRCh38**
		  - Human GRCh38, cellranger reference 1.2.0, Ensembl v84 gene annotation
		* - **hg19_v1.2.0** or **hg19**
		  - Human hg19, cellranger reference 1.2.0, Ensembl v82 gene annotation
		* - **mm10_v1.2.0** or **mm10**
		  - Mouse mm10, cellranger reference 1.2.0, Ensembl v84 gene annotation
		* - **GRCh38_and_mm10_v1.2.0** or **GRCh38_and_mm10**
		  - Human and mouse, built from GRCh38 and mm10 cellranger references, Ensembl v84 gene annotations are used

	Pre-built snRNA-seq references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38_premrna_v1.2.0** or **GRCh38_premrna**
		  - Human, introns included, built from GRCh38 cellranger reference 1.2.0, Ensembl v84 gene annotation, treating annotated transcripts as exons
		* - **mm10_premrna_v1.2.0** or **mm10_premrna**
		  - Mouse, introns included, built from mm10 cellranger reference 1.2.0, Ensembl v84 gene annotation, treating annotated transcripts as exons
		* - **GRCh38_premrna_and_mm10_premrna_v1.2.0** or **GRCh38_premrna_and_mm10_premrna**
		  - Human and mouse, introns included, built from GRCh38_premrna_v1.2.0 and mm10_premrna_v1.2.0

#. **Index** column.

	Put `10x single cell 3' sample index set names`_ (e.g. SI-GA-A12) here. 

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
		* - **SC3Pv1**
		  - Single Cell 3′ v1
		* - **SC3Pv2**
		  - Single Cell 3′ v2
		* - **SC3Pv3**
		  - Single Cell 3′ v3. You should set cellranger version input parameter to >= 3.0.2
		* - **SC5P-PE**
		  - Single Cell 5′ paired-end (both R1 and R2 are used for alignment)
		* - **SC5P-R2**
		  - Single Cell 5′ R2-only (where only R2 is used for alignment)

#. *DataType* column.
	
	This column is optional with a default **rna**. If you want to put a value, put **rna** here.

#. *FetureBarcodeFile* column.

	Leave it blank for scRNA-seq and snRNA-seq.

#. Example::

	Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
	sample_1,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,rna
	sample_2,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,SC3Pv3,rna
	sample_3,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime,rna
	sample_4,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8,fiveprime,rna
	sample_1,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime,rna
	sample_2,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,SC3Pv3,rna
	sample_3,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime,rna
	sample_4,mm10_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8,fiveprime,rna


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
		  - Results are written to $output_directory/$bcl_directory_fastqs/fastq_path/ and will overwrite any existing files at this location.
		* - run_mkfastq
		  - If you want to run ``cellranger mkfastq``
		  - true
		  - true
		* - run_count
		  - If you want to run ``cellranger count``
		  - true
		  - true
		* - delete_input_directory
		  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges 
		  - false
		  - false
		* - force_cells
		  - Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells
		  - 6000
		  - 
		* - expect_cells
		  - Expected number of recovered cells. Mutually exclusive with force_cells
		  - 3000
		  - 
		* - secondary
		  - Perform Cell Ranger secondary analysis (dimensionality reduction, clustering, etc.)
		  - false
		  - false
		* - cellranger_version
		  - cellranger version, could be 3.1.0, 3.0.2, or 2.2.0
		  - "3.1.0"
		  - "3.1.0"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "cumulusprod" for Docker Hub images; 

		  	- "quay.io/cumulus" for backup images on Red Hat registry.
		  - "cumulusprod"
		  - "cumulusprod"
		* - cellranger_mkfastq_docker_registry
		  - Docker registry to use for ``cellranger mkfastq``. 
		    Default is the registry to which only Broad users have access. 
		    See :ref:`bcl2fastq-docker` for making your own registry.
		  - "gcr.io/broad-cumulus"
		  - "gcr.io/broad-cumulus"
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
		* - preemptible
		  - Number of preemptible tries
		  - 2
		  - 2

Workflow output
+++++++++++++++

See the table below for important sc/snRNA-seq outputs.

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
	  - A list of google bucket urls containing count matrices, one url per sample.
	* - metrics_summaries
	  - File
	  - A excel spreadsheet containing QCs for each sample.
	* - output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each sample (cellranger count output).
	* - count_matrix
	  - String
	  - gs url for a template count_matrix.csv to run Cumulus.

---------------------------------

Feature barcoding assays (cell & nucleus hashing, CITE-seq and Perturb-seq)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cellranger_workflow`` can extract feature-barcode count matrices in CSV format for feature barcoding assays such as *cell and nucleus hashing*, *CITE-seq*, and *Perturb-seq*. For cell and nucleus hashing as well as CITE-seq, the feature refers to antibody. For Perturb-seq, the feature refers to guide RNA. Please follow the instructions below to configure ``cellranger_workflow``.

Prepare feature barcode files
+++++++++++++++++++++++++++++

	Prepare a CSV file with the following format: feature_barcode,feature_name.
	See below for an example::

		TTCCTGCCATTACTA,sample_1
		CCGTACCTCATTGTT,sample_2
		GGTAGATGTCCTCAG,sample_3
		TGGTGTCATTCTTGA,sample_4

	The above file describes a cell hashing application with 4 samples.

	Then upload it to your google bucket::

		gsutil antibody_index.csv gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv


Sample sheet
++++++++++++

#. **Reference** column.

	This column is not used for extracting feature-barcode count matrix. To be consistent, please put the reference for the associated scRNA-seq assay here.

#. **Index** column.

	The index can be either Illumina index primer sequence (e.g. ``ATTACTCG``, also known as ``D701``), or `10x single cell 3' sample index set names`_ (e.g. SI-GA-A12). 

	**Note 1**: All index sequences (including 10x's) should have the same length (8 bases). If one index sequence is shorter (e.g. ATCACG), pad it with P7 sequence (e.g. ATCACGAT).

	**Note 2**: It is users' responsibility to avoid index collision between 10x genomics' RNA indexes (e.g. SI-GA-A8) and Illumina index sequences for used here (e.g. ``ATTACTCG``).

#. *Chemistry* column.
	
	The following keywords are accepted for *Chemistry* column:

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Chemistry
		  - Explanation
		* - **SC3Pv3**
		  - Single Cell 3′ v3 (default).
		* - **SC3Pv2**
		  - Single Cell 3′ v2
		* - **fiveprime**
		  - Single Cell 5′
		* - **SC5P-PE**
		  - Single Cell 5′ paired-end (both R1 and R2 are used for alignment)
		* - **SC5P-R2**
		  - Single Cell 5′ R2-only (where only R2 is used for alignment)

#. *DataType* column.
	
	Put **adt** here if the assay is CITE-seq, cell or nucleus hashing. Put **crispr** here if Perturb-seq.

#. *FetureBarcodeFile* column.

	Put Google Bucket URL of the feature barcode file here.	

#. Example::

	Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
	sample_1_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,rna
	sample_1_adt,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,ATTACTCG,threeprime,adt,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv
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
		* - delete_input_directory
		  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges 
		  - false
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
		  - cellranger version, could be 3.1.0, 3.0.2, 2.2.0
		  - "3.1.0"
		  - "3.1.0"
		* - cumulus_feature_barcoding_version
		  - Cumulus_feature_barcoding version for extracting feature barcode matrix. Version available: 0.2.0.
		  - "0.2.0"
		  - "0.2.0"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "cumulusprod" for Docker Hub images; 

		  	- "quay.io/cumulus" for backup images on Red Hat registry.
		  - "cumulusprod"
		  - "cumulusprod"
		* - mkfastq_docker_registry
		  - Docker registry to use for ``cellranger mkfastq``. 
		    Default is the registry to which only Broad users have access. 
		    See :ref:`bcl2fastq-docker` for making your own registry.
		  - "gcr.io/broad-cumulus"
		  - "gcr.io/broad-cumulus"
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
		* - preemptible
		  - Number of preemptible tries
		  - 2
		  - 2

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

If data type is ``crispr``, three additional files, ``sample_id.umi_count.pdf``, ``sample_id.filt.csv`` and ``sample_id.filt.stat.csv.gz``, are generated.

``sample_id.umi_count.pdf`` plots number of UMIs against UMI with certain number of reads and colors UMIs with high likelihood of being chimeric in blue and other UMIs in red. This plot is generated purely based on number of reads each UMI has.

``sample_id.filt.csv`` is the filtered feature count matrix. It has the same format as ``sample_id.csv``.

``sample_id.filt.stat.csv.gz`` is the filtered sufficient statistics. It has the same format as ``sample_id.stat.csv.gz``.

---------------------------------

Single-cell ATAC-seq
^^^^^^^^^^^^^^^^^^^^

To process scATAC-seq data, follow the specific instructions below.

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built scATAC-seq references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38_atac_v1.2.0**
		  - Human GRCh38, cellranger-atac reference 1.2.0
		* - **mm10_atac_v1.2.0** 
		  - Mouse mm10, cellranger-atac reference 1.2.0
		* - **hg19_atac_v1.2.0**
		  - Human hg19, cellranger-atac reference 1.2.0
		* - **b37_atac_v1.2.0**
		  - Human b37 build, cellranger-atac reference 1.2.0
		* - **GRCh38_and_mm10_atac_v1.2.0**
		  - Human GRCh38 and mouse mm10, cellranger-atac reference 1.2.0
		* - **hg19_and_mm10_atac_v1.2.0**
		  - Human hg19 and mouse mm10, cellranger-atac reference 1.2.0
		* - **GRCh38_atac_v1.1.0**
		  - Human GRCh38, cellranger-atac reference 1.1.0
		* - **mm10_atac_v1.1.0** 
		  - Mouse mm10, cellranger-atac reference 1.1.0
		* - **hg19_atac_v1.1.0**
		  - Human hg19, cellranger-atac reference 1.1.0
		* - **b37_atac_v1.1.0**
		  - Human b37 build, cellranger-atac reference 1.1.0
		* - **GRCh38_and_mm10_atac_v1.1.0**
		  - Human GRCh38 and mouse mm10, cellranger-atac reference 1.1.0
		* - **hg19_and_mm10_atac_v1.1.0**
		  - Human hg19 and mouse mm10, cellranger-atac reference 1.1.0

#. **Index** column.

	Put `10x single cell ATAC sample index set names`_ (e.g. SI-NA-B1) here.

#. *Chemistry* column.
	
	This column is not used for scATAC-seq data. Put **auto** here as a placeholder if you decide to include the Chemistry column. 

#. *DataType* column.
	
	Set it to **atac**.

#. *FetureBarcodeFile* column.

	Leave it blank for scATAC-seq.

#. Example::

	Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType
	sample_atac,GRCh38_atac_v1.1.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9YB,*,SI-NA-A1,auto,atac

Workflow input
++++++++++++++

``cellranger_workflow`` takes Illumina outputs as input and runs ``cellranger-atac mkfastq`` and ``cellranger-atac count``. Please see the description of inputs below. Note that required inputs are shown in bold.

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
	  - If you want to run ``cellranger-atac mkfastq``
	  - true
	  - true
	* - run_count
	  - If you want to run ``cellranger-atac count``
	  - true
	  - true
	* - delete_input_directory
	  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges 
	  - false
	  - false
	* - force_cells
	  - Force pipeline to use this number of cells, bypassing the cell detection algorithm
	  - 6000
	  - 
	* - cellranger_atac_version
	  - cellranger-atac version, currently only 1.1.0
	  - "1.1.0"
	  - "1.1.0"
	* - docker_registry
	  - Docker registry to use for cellranger_workflow. Options:

	  	- "cumulusprod" for Docker Hub images; 

	  	- "quay.io/cumulus" for backup images on Red Hat registry.
	  - "cumulusprod"
	  - "cumulusprod"
	* - zones
	  - Google cloud zones
	  - "us-central1-a us-west1-a"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - atac_num_cpu
	  - Number of cpus for cellranger-atac count
	  - 64
	  - 64
	* - atac_memory
	  - Memory string for cellranger-atac count
	  - "57.6G"
	  - "57.6G"
	* - mkfastq_disk_space
	  - Optional disk space in GB for cellranger-atac mkfastq
	  - 1500
	  - 1500
	* - atac_disk_space
	  - Disk space in GB needed for cellranger-atac count
	  - 500
	  - 500
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

Workflow output
+++++++++++++++

See the table below for important scATAC-seq outputs.

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
	  - A list of google bucket urls containing cellranger-atac count outputs, one url per sample.
	* - metrics_summaries
	  - File
	  - A excel spreadsheet containing QCs for each sample.
	* - output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each sample (cellranger count output).
	* - count_matrix
	  - String
	  - gs url for a template count_matrix.csv to run cumulus.

Aggregate scATAC-Seq Samples
+++++++++++++++++++++++++++++

To aggregate multiple scATAC-Seq samples, follow the instructions below:

1. Import ``cellranger_atac_aggr`` workflow. Please see Step 1 `here <./cellranger.html#a-general-step-by-step-instruction>`_, and the name of workflow is "**cumulus/cellranger_atac_aggr**".

2. Set the inputs of workflow. Please see the description of inputs below. Notice that required inputs are shown in bold:

.. list-table::
	:widths: 5 30 30 20
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **aggr_id**
	  - Aggregate ID.
	  - "aggr_sample"
	  -
	* - **input_counts_directories**
	  - A string contains comma-separated URLs to directories of samples to be aggregated.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/data/sample1,gs://fc-e0000000-0000-0000-0000-000000000000/data/sample2"
	  -
	* - **output_directory**
	  - Output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/aggregate_result"
	  -
	* - **genome**
	  - The reference genome name used by Cell Ranger, can be either a keyword of pre-built genome, or a Google Bucket URL. See `this table <./cellranger.html#single-cell-and-single-nucleus-rna-seq>`_ for the list of keywords of pre-built genomes.
	  - "GRCh38_atac_v1.2.0"
	  -
	* - normalize
	  - Sample normalization mode. 
	    Options are: ``none``, ``depth``, or ``signal``.
	  - "none"
	  - "none"
	* - secondary
	  - Perform secondary analysis (dimensionality reduction, clustering and visualization).
	  - false
	  - false
	* - dim_reduce
	  - Chose the algorithm for dimensionality reduction prior to clustering and tsne. 
	    Options are: ``lsa``, ``plsa``, or ``pca``.
	  - "lsa"
	  - "lsa"
	* - cellranger_atac_version
	  - Cell Ranger ATAC version to use. 
	    Options: ``1.2.0``.
	  - "1.2.0"
	  - "1.2.0"
	* - zones
	  - Google cloud zones
	  - “us-central1-a us-west1-a”
	  - "us-central1-b"
	* - num_cpu
	  - Number of cpus to request for cellranger atac aggr.
	  - 64
	  - 64
	* - memory
	  - Memory size string for cellranger atac aggr.
	  - "57.6G"
	  - "57.6G"
	* - disk_space
	  - Disk space in GB needed for cellranger atac aggr.
	  - 500
	  - 500
	* - preemptible
	  - Number of preemptible tries.
	  - 2
	  - 2
	* - docker_registry
	  - Docker registry to use for cellranger_workflow. Options:

	  	- "cumulusprod" for Docker Hub images; 

	  	- "quay.io/cumulus" for backup images on Red Hat registry.
	  - "cumulusprod"
	  - "cumulusprod"

3. Check out the output in ``output_directory/aggr_id`` folder, where ``output_directory`` and ``aggr_id`` are the inputs you set in Step 2.

---------------------------------

Single-cell immune profiling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
		* - **GRCh38_vdj_v3.1.0**
		  - Human GRCh38 V(D)J sequences, cellranger reference 3.1.0, annotation built from Ensembl *Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf*
		* - **GRCm38_vdj_v3.1.0**
		  - Mouse GRCm38 V(D)J sequences, cellranger reference 3.1.0, annotation built from Ensembl *Mus_musculus.GRCm38.94.gtf*
		* - **GRCh38_vdj_v2.0.0** or **GRCh38_vdj**
		  - Human GRCh38 V(D)J sequences, cellranger reference 2.0.0, annotation built from Ensembl *Homo_sapiens.GRCh38.87.chr_patch_hapl_scaff.gtf* and *vdj_GRCh38_alts_ensembl_10x_genes-2.0.0.gtf*
		* - **GRCm38_vdj_v2.2.0** or **GRCm38_vdj**
		  - Mouse GRCm38 V(D)J sequences, cellranger reference 2.2.0, annotation built from Ensembl *Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf*

#. **Index** column.

	Put `10x single cell V(D)J sample index set names`_ (e.g. SI-GA-A3) here.

#. *Chemistry* column.
	
	This column is not used for scIR-seq data. Put **fiveprime** here as a placeholder if you decide to include the Chemistry column. 

#. *DataType* column.
	
	Set it to **vdj**.

#. *FetureBarcodeFile* column.

	Leave it blank for scIR-seq.

#. Example::

	Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType
	sample_vdj,GRCh38_vdj_v3.1.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ,1,SI-GA-A1,fiveprime,vdj

Workflow input
++++++++++++++

For scIR-seq data, ``cellranger_workflow`` takes Illumina outputs as input and runs ``cellranger mkfastq`` and ``cellranger vdj``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

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
	* - delete_input_directory
	  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges 
	  - false
	  - false
	* - force_cells
	  - Force pipeline to use this number of cells, bypassing the cell detection algorithm
	  - 6000
	  - 
	* - vdj_denovo
	  - Do not align reads to reference V(D)J sequences before de novo assembly
	  - false
	  - false
	* - cellranger_version
	  - cellranger version, could be 3.1.0, 3.0.2, 2.2.0 
	  - "3.1.0"
	  - "3.1.0"
	* - docker_registry
	  - Docker registry to use for cellranger_workflow. Options:

	  	- "cumulusprod" for Docker Hub images; 

	  	- "quay.io/cumulus" for backup images on Red Hat registry.
	  - "cumulusprod"
	  - "cumulusprod"
	* - cellranger_mkfastq_docker_registry
	  - Docker registry to use for ``cellranger mkfastq``. 
	    Default is the registry to which only Broad users have access. 
	    See :ref:`bcl2fastq-docker` for making your own registry.
	  - "gcr.io/broad-cumulus"
	  - "gcr.io/broad-cumulus"
	* - zones
	  - Google cloud zones
	  - "us-central1-a us-west1-a"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - num_cpu
	  - Number of cpus to request for one node for cellranger mkfastq and cellranger vdj
	  - 32
	  - 32
	* - memory
	  - Memory size string for cellranger mkfastq and cellranger vdj
	  - "120G"
	  - "120G"
	* - mkfastq_disk_space
	  - Optional disk space in GB for mkfastq
	  - 1500
	  - 1500
	* - vdj_disk_space
	  - Disk space in GB needed for cellranger vdj
	  - 500
	  - 500
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

Workflow output
+++++++++++++++

See the table below for important scIR-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_fastqs_directory
	  - Array[String]
	  - A list of google bucket urls containing FASTQ files, one url per flowcell.
	* - output_vdj_directory
	  - Array[String]
	  - A list of google bucket urls containing vdj results, one url per sample.
	* - metrics_summaries
	  - File
	  - A excel spreadsheet containing QCs for each sample.
	* - output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each sample (cellranger count output).
	* - count_matrix
	  - String
	  - gs url for a template count_matrix.csv to run cumulus.

---------------------------------

Build Cell Ranger References
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide routines wrapping Cell Ranger tools to build references for sc/snRNA-seq, scATAC-seq and single-cell immune profiling data.

Build references for sc/snRNA-seq
+++++++++++++++++++++++++++++++++

We provide a wrapper of ``cellranger mkref`` to build sc/snRNA-seq references. Please follow the instructions below.

1. Import ``cellranger_create_reference``
==============================================

	Import *cellranger_create_reference* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *cellranger_workflow* workflow is under ``Broad Methods Repository`` with name "**cumulus/cellranger_create_reference**".

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
		  - cellranger version, could be 3.1.0, 3.0.2, or 2.2.0
		  - "3.1.0"
		  - "3.1.0"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "cumulusprod" for Docker Hub images; 

		  	- "quay.io/cumulus" for backup images on Red Hat registry.
		  - "cumulusprod"
		  - "cumulusprod"
		* - zones
		  - Google cloud zones
		  - "us-central1-a us-west1-a"
		  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - num_cpu
		  - Number of cpus to request for one node for building indices
		  - 1
		  - 1
		* - memory
		  - Memory size string for cellranger-atac mkref
		  - "32G"
		  - "32G"
		* - disk_space
		  - Optional disk space in GB
		  - 100
		  - 100
		* - preemptible
		  - Number of preemptible tries
		  - 2
		  - 2

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

	Import *cellranger_atac_create_reference* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *cellranger_workflow* workflow is under ``Broad Methods Repository`` with name "**cumulus/cellranger_atac_create_reference**".

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
		* - **config_json**
		  - Configuration file defined in `10x genomics configuration file`_. Note that links to files in the JSON must be Google bucket URLs
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/config.json"
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_atac_reference"
		  -
		* - cellranger_atac_version
		  - cellranger-atac version, could be 1.1.0
		  - "1.1.0"
		  - "1.1.0"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "cumulusprod" for Docker Hub images; 

		  	- "quay.io/cumulus" for backup images on Red Hat registry.
		  - "cumulusprod"
		  - "cumulusprod"
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
		* - preemptible
		  - Number of preemptible tries
		  - 2
		  - 2

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

	Import *cellranger_vdj_create_reference* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *cellranger_workflow* workflow is under ``Broad Methods Repository`` with name "**cumulus/cellranger_vdj_create_reference**".

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
		  - cellranger version, could be 3.1.0, 3.0.2, or 2.2.0
		  - "3.1.0"
		  - "3.1.0"
		* - docker_registry
		  - Docker registry to use for cellranger_workflow. Options:

		  	- "cumulusprod" for Docker Hub images; 

		  	- "quay.io/cumulus" for backup images on Red Hat registry.
		  - "cumulusprod"
		  - "cumulusprod"
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
		* - preemptible
		  - Number of preemptible tries
		  - 2
		  - 2

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




.. _10x genomics v2 cell barcode white list: gs://regev-lab/resources/cellranger/737K-august-2016.txt.gz
.. _10x genomics v3 cell barcode white list: gs://regev-lab/resources/cellranger/3M-february-2018.txt.gz
.. _10x single cell 3' sample index set names: https://support.10xgenomics.com/single-cell-gene-expression/index/doc/specifications-sample-index-sets-for-single-cell-3
.. _10x single cell ATAC sample index set names: https://support.10xgenomics.com/single-cell-atac/sequencing/doc/specifications-sample-index-sets-for-single-cell-atac
.. _10x single cell V(D)J sample index set names: https://support.10xgenomics.com/single-cell-vdj/sequencing/doc/specifications-sample-index-sets-for-single-cell-vdj
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/
.. _10x genomics configuration file: https://support.10xgenomics.com/single-cell-atac/software/pipelines/latest/advanced/references#config
.. _10x filter with mkgtf: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#mkgtf
.. _10x build Cell Ranger compatible pre-mRNA Reference Package: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references#premrna

