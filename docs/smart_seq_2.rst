Extract gene-count matrices from plated-based SMART-Seq2 data
-------------------------------------------------------------

Follow the steps below to extract gene-count matrices from SMART-Seq2 data on Terra_. This WDL aligns reads using Bowtie 2 and estimates expression levels using RSEM.

#. Copy your sequencing output to your workspace bucket using gsutil_ in your unix terminal.

	You can obtain your bucket URL in the dashboard tab of your Terra workspace under the information panel.

	.. image:: images/google_bucket_link.png

	Note: Broad users need to be on an UGER node (not a login node) in order to use the ``-m`` flag

	Request an UGER node::

		reuse UGER
		qrsh -q interactive -l h_vmem=4g -pe smp 8 -binding linear:8 -P regevlab

	The above command requests an interactive node with 4G memory per thread and 8 threads. Feel free to change the memory, thread, and project parameters.

	Once you're connected to an UGER node, you can make gsutil_ available by running::

		reuse Google-Cloud-SDK

	Use ``gsutil cp [OPTION]... src_url dst_url`` to copy data to your workspace bucket.
	For example, the following command copies the directory at /foo/bar/nextseq/Data/VK18WBC6Z4 to a Google bucket::

		gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4

	``-m`` means copy in parallel, ``-r`` means copy the directory recursively.


#. Create a sample sheet. 

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet provides metadata for each cell:

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


#. Import smartseq2 tool.

	In Terra, select the ``Tools`` tab, then click ``Find a Tool``. Click ``Broad Methods Repository``. Type **scCloud/smartseq2**.
 	You can also see the Terra documentation for `adding a tool`_.

#. Select ``Process single workflow from files``.

	.. image:: images/single_workflow.png


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
	* - smartseq2_version
	  - SMART-Seq2 docker version
	  - "1.0.0"
	  - "1.0.0"
	* - zones
	  - Google cloud zones
	  - "us-east1-b us-east1-c us-east1-d"
	  - "us-east1-b us-east1-c us-east1-d"
	* - num_cpu
	  - Number of cpus to request for one node
	  - 4
	  - 4
	* - memory
	  - Memory size string
	  - "3.60G"
	  - "3.60G"
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



---------------------------------

Custom Genome:
^^^^^^^^^^^^^^
#. Install RSEM 1.3.1 and Bowtie 2 v2.3.4.3.
#. Download a genome fasta file (e.g. Homo_sapiens.GRCh38.dna.primary_assembly.fa) and a GTF gene annotation file (e.g. Homo_sapiens.GRCh38.83.gtf).
#. Run the following, substituting your fasta file, GTF file, and the number of threads to use::

	mkdir rsem_ref
	rsem-prepare-reference --gtf Homo_sapiens.GRCh38.83.gtf --bowtie2 -p n_threads Homo_sapiens.GRCh38.dna.primary_assembly.fa rsem_ref/rsem_ref
	tar -czf rsem_ref.tar.gz rsem_ref

#. Upload the rsem_ref.tar.gz to your google bucket and set the URL to this file as the reference input field value.


.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a tool: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/
