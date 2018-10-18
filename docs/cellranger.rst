Run Cell Ranger mkfastq/count/vdj
---------------------------------

Follow the steps below to run CellRanger mkfastq/count/vdj on FireCloud.

#. Copy your sequencing output to your workspace bucket using gsutil in your unix terminal. You can obtain your bucket URL in the workspace summary tab in FireCloud under Google Bucket. You can also read `FireCloud instructions`_ on uploading data.
	
	Example of copying the directory at /foo/bar/nextseq/Data/VK18WBC6Z4 to a Google Cloud bucket::

		gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4
	
	``-m`` means copy in parallel, ``-r`` means copy the directory recursively.
	
	Note: Broad users need to be on an UGER node (not a login node) in order to use the ``-m`` flag

	Request an UGER server::

		reuse UGER
		qrsh -q interactive -l h_vmem=4g -pe smp 8 -binding linear:8 -P regevlab

	The above command requests an interactive shell with 4G memory per thread and 8 threads. Feel free to change the memory, thread, and project parameters.

	Once you've connected to an UGER node run::
		reuse Google-Cloud-SDK

	to make the Google Cloud tools available


#. Create a scRNA-Seq formatted sample sheet. 

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to demultiplex flowcells and generate channel-specific count matrices. Note that *Sample*, *Lane*, and *Index* columns are defined exactly the same as in 10x's simple CSV layout file.

	scRNA-Seq formatted sample sheet description (required column headers are shown in bold):

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name.
		* - **Reference**
		  - 
			| Provides the reference genome used by *cellranger count* for each 10x channel. 
			| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as
			| **GRCh38** for human GRCh38,
			| **hg19** for human hg19,
			| **mm10** for mouse, 
			| **GRCh38_and_mm10** for human and mouse,
			| **GRCh38_premrna** for human, introns included,
			| **mm10_premrna** for mouse, introns included, 
			| **GRCh38_premrna_and_mm10_premrna** for human and mouse, introns included,
			| **GRCh38_vdj** for human V(D)J sequences, and
			| **GRCm38_vdj** for mouse V(D)J sequences.
		* - **Flowcell**
		  - Indicates the Google bucket URL of uploaded BCL folders.
		* - **Lane**
		  - Tells which lanes the sample was pooled into.
		* - **Index**
		  - Contains 10x sample index set names (e.g. SI-GA-A12).
		* - Chemistry
		  - Describes the 10x chemistry used for the sample. If this column is omitted, *cellranger count* will try to determine the chemistry automatically. This column is optional.
		* - DataType
		  - Describes the data type of the sample --- *count*, *vdj*, or *adt*. *count* refers to gene expression data (*cellranger count*), *vdj* refers to V(D)J data (*cellranger vdj*), and *adt* refers to antibody tag data. This column is optional and the default data type is *count*.   

	The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 4 samples sequenced in two flowcells.

	Example::

		Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,count
		sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,threeprime,count
		sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime,count
		sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8,fiveprime,count
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime,count
		sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,threeprime,count
		sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime,count
		sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8,fiveprime,count
		sample_5,GRCh38_vdj,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ,1,SI-GA-A1,fiveprime,vdj
		sample_6,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ,2,AGATCCTT,threeprime,adt


#. Upload your sample sheet to the workspace bucket.

	Example::

		gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import cellranger_mkfastq_count method.

	In FireCloud, select the ``Method Configurations`` tab then click ``Import Configuration``. Click ``Import From Method Repository``. Type cellranger_mkfastq_count.

#. Uncheck ``Configure inputs/outputs using the Workspace Data Model``.


---------------------------------

Cell Ranger mkfastq/count/vdj inputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Cell Ranger mkfastq/count`` takes Illumina outputs as input and runs ``cellranger mkfastq`` and ``cellranger count``/``cellranger vdj``. Please see the description of inputs below. Note that required inputs are shown in bold.

.. list-table::
	:widths: 5 30 30 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_csv_file**
	  - Sample Sheet (contains Sample, Reference, Flowcell, Lane, Index)
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
	  - If you want to run ``cellranger count`` or ``cellranger vdj``
	  - true
	  - true
	* - delete_input_directory
	  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges 
	  - true
	  - true
	* - do_force_cells
	  - force cells
	  - true
	  - false
	* - force_cells
	  - Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells
	  - 3000
	  - 6000
	* - expect_cells
	  - Expected number of recovered cells. Mutually exclusive with force_cells
	  - 1000
	  - 3000
	* - secondary
	  - Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.)
	  - false
	  - false
	* - vdj_force_cells
	  - force pipeline to use this number of cells for the vdj task, bypassing the cell detection algorithm
	  - 2000
	  -
	* - vdj_denovo
	  - Do not align reads to reference V(D)J sequences before de novo assembly
	  - true
	  - false
	* - vdj_chain
	  - Force the web summary HTML and metrics summary CSV to only report on a particular chain type. The accepted values are: auto for autodetection based on TR vs IG representation, TR for T cell receptors, IG for B cell receptors, all for all chain types
	  - TR
	  - 
	* - antibody_barcode_file
	  - Antibody barcodes in csv format for the adt task
	  - antibody_barcodes.csv
	  -
	* - max_mismatch
	  - Maximum hamming distance in antibody barcodes for the adt task
	  - 3
	  - 3
	* - cellranger_version
	  - Cellranger version, could be 2.2.0 or 2.1.1
	  - "2.2.0"
	  - "2.2.0"
	* - num_cpu
	  - Number of cpus to request for one node
	  - 64
	  - 64
	* - memory
	  - Memory in GB
	  - 128
	  - 128
	* - adt_memory
	  - Optional memory in GB for extracting ADT count matrix
	  - 32
	  - 32
	* - mkfastq_disk_space
	  - Optional disk space in gigabytes for mkfastq
	  - 1500
	  - 1500
	* - count_disk_space
	  - Disk space in gigabytes needed for cellranger count
	  - 500
	  - 500
	* - vdj_disk_space
	  - Disk space in gigabytes needed for cellranger vdj
	  - 500
	  - 500
	* - adt_disk_space
	  - Disk space in gigabytes needed for extracting adt counts
	  - 100
	  - 100
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

---------------------------------

Cell Ranger mkfastq/count/vdj outputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the table below for important *Cell Ranger mkfastq/count* outputs.


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
	* - output_vdj_directory
	  - Array[String]
	  - A list of google bucket urls containing vdj results, one url per sample.
	* - output_adt_directory
	  - Array[String]
	  - A list of google bucket urls containing adt count matrices, one url per sample.	  
	* - metrics_summaries
	  - File
	  - A excel spreadsheet containing QCs for each sample.
	* - output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each sample (cellranger count output).
	* - count_matrix
	  - String
	  - gs url for a template count_matrix.csv to run scrtools.

---------------------------------

Only run ``cellranger count``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, people might want to perform demultiplexing locally and only run ``cellranger count`` on the cloud. This section describes how to only run ``cellranger count``  via ``cellranger_mkfastq_count``.

#. Copy your FASTQ files to the workspace using gsutil in your unix terminal. 

	You should upload folders of FASTQS. Each folder should contain all FASTQ files for one sample.

	Example::

		gsutil -m cp -r /foo/bar/fastq_path/K18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq

	``-m`` means copy in parallel, ``-r`` means copy the directory recursively.
	
	Note: Broad users need to be on an UGER node (not a login node) in order to use the ``-m`` flag
	
	You can also read `FireCloud instructions`_ on uploading data.

#. Create scRNA-Seq formatted sample sheet for cell ranger count only (required column headers are shown in bold):

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name.
		* - **Reference**
		  - 
			| Provides the reference genome used by *cellranger count*.
			| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as
			| **GRCh38** for human GRCh38,
			| **hg19** for human hg19,
			| **mm10** for mouse, 
			| **GRCh38_and_mm10** for human and mouse,
			| **GRCh38_premrna** for human, introns included,
			| **mm10_premrna** for mouse, introns included,
			| **GRCh38_premrna_and_mm10_premrna** for human and mouse, introns included,
			| **GRCh38_vdj** for human V(D)J sequences, and
			| **GRCm38_vdj** for mouse V(D)J sequences.
		* - **Flowcell**
		  - Indicates the Google bucket URL of the uploaded FASTQ folders. The full path to the FASTQ files is FlowCell/Sample
		* - Chemistry
		  - Describe the 10x chemistry used for the sample. This column is optional. If this column is omitted, *cellranger count* will try to determine the chemistry automatically.
		* - DataType
		  - Describes the data type of the sample --- *count*, *vdj*, or *adt*. *count* refers to gene expression data (*cellranger count*), *vdj* refers to V(D)J data (*cellranger vdj*), and *adt* refers to antibody tag data. This column is optional and the default data type is *count*.

	In the following example sample_1 is sequenced on 2 flowcells. The FASTQ files for flowcell_1 are located at gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_1/sample_1 while the FASTQ files for flowcell_2 are located at gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_2_sample1::

		Sample,Reference,Flowcell
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_1
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_2

#. Set optional input ``run_mkfastq`` to ``false``.

---------------------------------

Run CITE-Seq/Cell-hashing/Nuclei-hashing
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This WDL could extract ADT counts from *CITE-Seq/Cell-hashing/Nuclei-hashing* assays. Please follow the instructions below.

#. Add both RNA assay and ADT assay information into the sample sheet.

	See below for an example::

		Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType
		sample_1_rna,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,count
		sample_1_adt,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,ATTACTCG,threeprime,adt

	The second line in the above sheet describes the ADT part. The only difference between ADT and RNA parts is the *Index*. For the ADT part, the index is the Illumina index primer sequence (e.g. ATTACTCG).

#. Prepare an antibody barcode file and upload to the Google bucket.

	Prepare a CSV file with the following format: antibody_barcode,name.
	See below for an example::

		TTCCTGCCATTACTA,sample_1
		CCGTACCTCATTGTT,sample_2
		GGTAGATGTCCTCAG,sample_3
		TGGTGTCATTCTTGA,sample_4

	The above file describes a cell-hashing application with 4 samples.

#. Fill in the ADT-specific parameters:

	.. list-table::
		:widths: 5 30 30 5
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **antibody_barcode_file**
		  - Antibody barcode file in CSV format
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/antibody_barcode_file.csv"
		  -
		* - max_mismatch
		  - Maximum hamming distance in matching antibody barcodes
		  - 3
		  - 3
		* - adt_memory
		  - Optional memory in GB for extracting ADT count matrix
		  - 32
		  - 32
		* - adt_disk_space
		  - Optional disk space needed for extracting ADT count matrix
		  - 100
		  - 100


Extracted ADT output
++++++++++++++++++++

For each ADT sample, a folder with the sample ID is generated under ``cellranger_output_directory``. In the folder, two files --- ``sample_id.csv`` and ``sample_id.stat.csv`` are generated.

``sample_id.csv`` has the following format. The first line describes the column names: ``Antibody,cell_barcode_1,cell_barcode_2,...,cell_barcode_n``. The following lines describe UMI counts for each antibody barcode, with the following format: ``name,umi_count_1,umi_count_2,...,umi_count_n``.

``sample_id.stat.csv`` has the following format. The first line describes the column names: ``Barcode,Total_reads,Total_umis``. The following lines describe all cellular barcodes with at least UMI count, with the following format: ``cell_barcode,number_of_reads,number_of_umis``.


.. _FireCloud instructions: https://software.broadinstitute.org/firecloud/documentation/article?id=10574
