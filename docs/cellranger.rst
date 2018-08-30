Run Cell Ranger mkfastq/count
-----------------------------

Follow the steps below to run CellRanger mkfastq/count on FireCloud.

#. Copy your sequencing output to the workspace bucket using gsutil in your unix terminal. 

	It is highly recommended that you delete the **BCL files** after the pipeline is finished by turning on the **delete_input_directory** option.
    
	Example of copying a directory to a Google Cloud bucket::

		gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4
    
    ``-m`` means copy in parallel, ``-r`` means copy the directory recursively.
    
	Note: Broad users need to be on an UGER node (not a login node) in order to use the ``-m`` flag
    
	You can also read `FireCloud instructions`_ on uploading data.

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
		  - | Provides the reference genome used by *cellran ger count* for each 10x channel. 
		    | The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as
		    | **GRCh38** for human, 
		    | **mm10** for mouse, 
		    | **GRCh38_and_mm10** for human and mouse,
		    | **GRCh38_premrna** for human, introns included, and
		    | **mm10_premrna** for mouse, introns included.
		* - **Flowcell**
		  - Indicates the Google bucket URL of uploaded BCL folders.
		* - **Lane**
		  - Tells which lanes the sample was pooled into.
		* - **Index**
		  - Contains 10x sample index set names (e.g. SI-GA-A12).
		* - Chemistry
		  - Optionally describe the 10x chemistry used for the sample. If this column is omitted, *cellranger count* will try to determine the chemistry automatically.

	The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 4 samples sequenced in two flowcells.

	Example::
   
		Sample,Reference,Flowcell,Lane,Index,Chemistry
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime
		sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,threeprime
		sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime
		sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8,fiveprime
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime
		sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,threeprime
		sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime
		sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8,fiveprime


#. Upload your sample sheet to the workspace bucket.

	Example::

		gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import cellranger_mkfastq_count method.

	In FireCloud, select the ``Method Configurations`` tab then click ``Import Configuration``. Click ``Import From Method Repository``. Type cellranger_mkfastq_count.

#. Uncheck ``Configure inputs/outputs using the Workspace Data Model``.


---------------------------------

Cell Ranger mkfastq/count inputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``Cell Ranger mkfastq/count`` takes Illumina outputs as input and runs ``cellranger mkfastq`` and ``cellranger count``. Please see the description of inputs below. Note that required inputs are shown in bold.

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
	* - **cellranger_output_directory**
	  - Cellranger output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
	  -
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
	  - true
	  - false
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
	* - cellranger_version
	  - Cellranger version
	  - "2.1.1"
	  - "2.1.1"
	* - num_cpu
	  - Number of cpus to request for one node
	  - 64
	  - 64
	* - memory
	  - Memory in GB
	  - 128
	  - 128
	* - mkfastq_disk_space
	  - Optional disk space in gigabytes for mkfastq
	  - 1500
	  - 1500
	* - count_disk_space
	  - Disk space in gigabytes needed for cell ranger count
	  - 500
	  - 500
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

---------------------------------

Cell Ranger mkfastq/count outputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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
    	* - metrics_summaries
    	  - File
    	  - A excel spreadsheet containing QCs for each sample.
    	* - output_web_summary
    	  - Array[File]
    	  - A list of htmls visualizing QCs for each sample (cellranger count output).

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
	  - | Provides the reference genome used by *cellranger count*.
		| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as
		| **GRCh38** for human,
		| **mm10** for mouse,
		| **GRCh38_and_mm10** for human and mouse,
		| **GRCh38_premrna** for human, introns included, and
		| **mm10_premrna** for mouse, introns included.
	* - **Flowcell**
	  - Indicates the Google bucket URL of the uploaded FASTQ folders. The full path to the FASTQ files is FlowCell/Sample
	* - Chemistry
	  - Optionally describe the 10x chemistry used for the sample. If this column is omitted, *cellranger count* will try to determine the chemistry automatically.


	In the following example sample_1 is sequenced on 2 flowcells. The FASTQ files for flowcell_1 are located at gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_1/sample_1
	while the FASTQ files for flowcell_2 are located at gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_2_sample1 ::

		Sample,Reference,Flowcell
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_1
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_2

#. Set optional input ``run_mkfastq`` to ``false``.

.. _FireCloud instructions: https://software.broadinstitute.org/firecloud/documentation/article?id=10574
