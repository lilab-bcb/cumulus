Extract gene-count matrices from plated-based SMART-Seq2 data
-------------------------------------------------------------

Follow the steps below to extract gene-count matrices from SMART-Seq2 data on FireCloud. This WDL aligns reads using Bowtie 2 and estimates expression levels using RSEM.

#. Copy your sequencing output to your workspace bucket using gsutil in your unix terminal. You can obtain your bucket URL in the workspace summary tab in FireCloud under Google Bucket. You can also read `FireCloud instructions`_ on uploading data.
	
	Example of copying the directory at /foo/bar/nextseq/Data/VK18WBC6Z4 to a Google Cloud bucket::

		gsutil -m cp -r /foo/bar/Data/smartseq2 gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2
	
	``-m`` means copy in parallel, ``-r`` means copy the directory recursively.
	
	Note: Broad users need to be on an UGER node (not a login node) in order to use the ``-m`` flag

	Request an UGER server::

		reuse UGER
		qrsh -q interactive -l h_vmem=4g -pe smp 8 -binding linear:8 -P regevlab

	The above command requests an interactive shell with 4G memory per thread and 8 threads. Feel free to change the memory, thread, and project parameters.

	Once you've connected to an UGER node run::
		reuse Google-Cloud-SDK

	to make the Google Cloud tools available


#. Create a sample sheet. 

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet provides metadata for each cell:

Cell, Plate, Read1, and Read2

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - Cell
		  - Cell name.
		* - Plate
		  - Plate name. Cells with the same plate name are from the same plate.
		* - Read1
		  - Location of the FASTQ file for read1 in the cloud (gsurl).
		* - Read2
		  - Location of the FASTQ file for read1 in the cloud (gsurl).

	Example::

		Cell,Plate,Read1,Read2
		cell-1,plate-1,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-1_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-1_L001_R2_001.fastq.gz
		cell-2,plate-1,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-2_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-2_L001_R2_001.fastq.gz
		cell-3,plate-2,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-3_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-3_L001_R2_001.fastq.gz
		cell-4,plate-2,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-4_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-4_L001_R2_001.fastq.gz


#. Upload your sample sheet to the workspace bucket.

	Example::

		gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import smartseq2 method.

	In FireCloud, select the ``Method Configurations`` tab then click ``Import Configuration``. Click ``Import From Method Repository``. Type **scCloud/smartseq2**.

#. Uncheck ``Configure inputs/outputs using the Workspace Data Model``.


---------------------------------

Inputs:
^^^^^^^

Please see the description of inputs below. Note that required inputs are shown in bold.

.. list-table::
	:widths: 5 30 30 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_csv_file**
	  - Sample Sheet (contains Cell, Plate, Read1, Read2)
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
	  - 
	* - **output_directory**
	  - Output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2_output"
	  -
	* - **reference**
	  - Reference transcriptome to align reads to. Currently we only have ``GRCh38`` for human and ``GRCm38`` for mouse
	  - GRCh38
	  - 
	* - num_cpu
	  - Number of cpus to request for one node
	  - 4
	  - 4
	* - memory
	  - Memory in GB
	  - 10
	  - 10
	* - disk_space
	  - Disk space in gigabytes
	  - 10
	  - 10
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

---------------------------------

Outputs:
^^^^^^^^

See the table below for important outputs.


.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_count_matrix
	  - Array[String]
	  - A list of google bucket urls containing gene-count matrices, one per plate. Each gene-count matrix file has the suffix 'dge.txt.gz'.

This WDL generates one gene-count matrix per SMART-Seq2 plate. The gene-count matrix uses Drop-Seq format: The first line starts with 'Gene' and then gives cell barcodes separated by tabs. Starting from the second line, each line describes one gene. The first item in the line is the gene name and the rest items are TPM-normalized count values of this gene for each cell. The gene-count matrices can be fed directly into ``scCloud`` for downstream analysis.

TPM-normalized counts are calculated as follows: We first estimate the gene expression levels in TPM using RSEM. Suppose we have ``c`` reads for one cell, we then calculate TPM-normalized count for gene i as ``TPM_i / 1e6 * c``. 

TPM-normalized counts reflect both the relative expression levels and the cell sequencing depth.


.. _FireCloud instructions: https://software.broadinstitute.org/firecloud/documentation/article?id=10574
