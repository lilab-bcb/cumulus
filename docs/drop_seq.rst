Drop-seq pipeline
-------------------------------------------------------------

Follow the steps below to extract gene-count matrices from Drop-seq data.
This WDL follows the steps outlined in the `Drop-seq alignment cookbook`_ from the `McCarroll lab`_.

#. Copy your sequencing output to your workspace bucket using gsutil in your unix terminal. You can obtain your bucket URL in the workspace summary tab in FireCloud under Google Bucket. You can also read `FireCloud instructions`_ on uploading data.

	Example of copying the directory at /foo/bar/nextseq/Data/VK18WBC6Z4 to a Google Cloud bucket::

		gsutil -m cp -r /foo/bar/Data/dropseq gs://fc-e0000000-0000-0000-0000-000000000000/dropseq

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

	Please note that the columns in the CSV must be in the specified order and contains no header line.
	The sample sheet provides either the fastq files for each sample if you've already run bcl2fastq or the BCL directories if you're starting from BCL directories.
	Please note that BCL directories must contain a valid bcl2fastq sample sheet (SampleSheet.csv):


	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - Name
		  - Sample name.
		* - Read1
		  - Location of the FASTQ file for read1 in the cloud (gsurl).
		* - Read2
		  - Location of the FASTQ file for read2 in the cloud (gsurl).

	Example::


		sample-1,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample1-1_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample-1_L001_R2_001.fastq.gz
		sample-2,plate-1,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample-2_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample-2_L001_R2_001.fastq.gz
		sample-1,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-2/sample1-1_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-2/sample-1_L001_R2_001.fastq.gz


	Note that in this example, sample-1 was sequenced across two flowcells.

#. Upload your sample sheet to the workspace bucket.

	Example::

		gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import dropseq_workflow method.

	In FireCloud, select the ``Method Configurations`` tab then click ``Import Configuration``. Click ``Import From Method Repository``. Type **scCloud/dropseq_workflow**.

#. Uncheck ``Configure inputs/outputs using the Workspace Data Model``.


---------------------------------

Inputs:
^^^^^^^

Please see the description of important inputs below.

.. list-table::
	:widths: 5 30
	:header-rows: 1

	* - Name
	  - Description
	* - input_csv_file
	  - Sample Sheet (contains Name, Read1, Read2 or a list of BCL directories e.g. "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv")
	* - output_directory
	  - Pipeline output directory (gs URL e.g. "gs://fc-e0000000-0000-0000-0000-000000000000/dropseq_output")
	* - run_bcl2fastq
	  - Whether your sample sheet contains one BCL directory per line or one sample per line
	* - star_memory
	  - For large genomes (e.g. human and mouse combined genome), enter "120G" or else the STAR aligner will take forever ($$$)
	* - force_cells
	  - If supplied, bypass the cell detection algorithm and use this number of cells
	* - star_genome_file
	  - "gs://regev-lab/resources/DropSeq/human/STAR2_5_index_hg19.tar.gz" or "gs://regev-lab/resources/DropSeq/mouse/STAR2_5_index_mm10.tar.gz"
	* - refflat
	  - "gs://regev-lab/resources/DropSeq/human/hg19.refFlat" or "gs://regev-lab/resources/DropSeq/mouse/mm10.refFlat"
	* - gene_intervals
	  - :"gs://regev-lab/resources/DropSeq/human/hg19.genes.intervals" or "gs://regev-lab/resources/DropSeq/mouse/mm10.genes.intervals"
	* - genome_fasta
	  - "gs://regev-lab/resources/DropSeq/human/hg19.fasta" or "gs://regev-lab/resources/DropSeq/mouse/mm10.fasta"
	* - genome_dict
	  - "gs://regev-lab/resources/DropSeq/human/hg19.dict" or "gs://regev-lab/resources/DropSeq/mouse/mm10.dict"


---------------------------------

Outputs:
^^^^^^^^

The pipeline outputs a list of google bucket urls containing one gene-count matrix per sample. Each gene-count matrix file has the suffix 'dge.txt.gz'.

.. _FireCloud instructions: https://software.broadinstitute.org/firecloud/documentation/article?id=10574
.. _Drop-seq alignment cookbook: https://github.com/broadinstitute/Drop-seq/blob/master/doc/Drop-seq_Alignment_Cookbook.pdf
.. _McCarroll lab: http://mccarrolllab.org/dropseq-1/
