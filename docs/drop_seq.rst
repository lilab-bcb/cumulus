Drop-seq pipeline
-------------------------------------------------------------

This workflow follows the steps outlined in the `Drop-seq alignment cookbook`_ from the `McCarroll lab`_ , except the default STAR aligner flags are *--limitOutSJcollapsed 1000000 --twopassMode Basic*.
Additionally the pipeline provides the option to generate count matrices using  `dropEst`_.

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


#. Non Broad Institute users that wish to run bcl2fastq must create a custom docker image.

    See :ref:`bcl2fastq-docker` instructions.

#. Create a sample sheet.

	Please note that the columns in the CSV must be in the order shown below and does not contain a header line.
	The sample sheet provides either the FASTQ files for each sample if you've already run bcl2fastq or a list of BCL directories if you're starting from BCL directories.
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

	Example using FASTQ input files::

		sample-1,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample1-1_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample-1_L001_R2_001.fastq.gz
		sample-2,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample-2_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-1/sample-2_L001_R2_001.fastq.gz
		sample-1,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-2/sample1-1_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/dropseq-2/sample-1_L001_R2_001.fastq.gz


	Note that in this example, sample-1 was sequenced across two flowcells.


	Example using BCL input directories::

		gs://fc-e0000000-0000-0000-0000-000000000000/flowcell-1
		gs://fc-e0000000-0000-0000-0000-000000000000/flowcell-2


	Note that the flow cell directory must contain a bcl2fastq sample sheet named SampleSheet.csv.

#. Upload your sample sheet to the workspace bucket.

	Example::

		gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import *dropseq_workflow* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *dropseq_workflow* is under ``Broad Methods Repository`` with name "**cumulus/dropseq_workflow**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace you want to export *dropseq_workflow* workflow in the drop-down menu.

#. In your workspace, open ``dropseq_workflow`` in ``WORKFLOWS`` tab. Select ``Process single workflow from files`` as below

	.. image:: images/single_workflow.png

   and click the ``SAVE`` button.

---------------------------------

Inputs
^^^^^^^

Please see the description of important inputs below.

.. list-table::
	:widths: 5 30
	:header-rows: 1

	* - Name
	  - Description
	* - input_csv_file
	  - CSV file containing sample name, read1, and read2 or a list of BCL directories.
	* - output_directory
	  - Pipeline output directory (gs URL e.g. "gs://fc-e0000000-0000-0000-0000-000000000000/dropseq_output")
	* - reference
	  - hg19, GRCh38, mm10, hg19_mm10, mmul_8.0.1 or a path to a custom reference JSON file
	* - run_bcl2fastq
	  - Whether your sample sheet contains one BCL directory per line or one sample per line (default false)
	* - run_dropseq_tools
	  - Whether to generate count matrixes using Drop-Seq tools from the `McCarroll lab`_ (default true)
	* - run_dropest
	  - Whether to generate count matrixes using `dropEst`_ (default false)
	* - cellular_barcode_whitelist
	  - Optional whitelist of known cellular barcodes
	* - drop_seq_tools_force_cells
	  - If supplied, bypass the cell detection algorithm (the elbow method) and use this number of cells.
	* - dropest_cells_max
	  - Maximal number of output cells
	* - dropest_genes_min
	  - Minimal number of genes for cells after the merge procedure (default 100)
	* - dropest_min_merge_fraction
	  - Threshold for the merge procedure (default 0.2)
	* - dropest_max_cb_merge_edit_distance
	  - Max edit distance between barcodes (default 2)
	* - dropest_max_umi_merge_edit_distance
	  - Max edit distance between UMIs (default 1)
	* - dropest_min_genes_before_merge
	  - Minimal number of genes for cells before the merge procedure. Used mostly for optimization. (default 10)
	* - dropest_merge_barcodes_precise
	  - Use precise merge strategy (can be slow), recommended to use when the list of real barcodes is not available (default true)
	* - dropest_velocyto
	  - Save separate count matrices for exons, introns and exon/intron spanning reads (default true)
	* - trim_sequence
	  - The sequence to look for at the start of reads for trimming (default "AAGCAGTGGTATCAACGCAGAGTGAATGGG")
	* - trim_num_bases
	  - How many bases at the beginning of the sequence must match before trimming occur (default 5)
	* - umi_base_range
	  - The base location of the molecular barcode (default 13-20)
	* - cellular_barcode_base_range
	  - The base location of the cell barcode (default 1-12)
	* - star_flags
	  - Additional options to pass to STAR aligner


Please note that run_bcl2fastq must be set to true if you're starting from BCL files instead of FASTQs.


Custom Genome JSON
===================

If you're reference is not one of the predefined choices, you can create a custom JSON file. Example::

	{
		"refflat":	  "gs://fc-e0000000-0000-0000-0000-000000000000/human_mouse/hg19_mm10_transgenes.refFlat",
		"genome_fasta":	   "gs://fc-e0000000-0000-0000-0000-000000000000/human_mouse/hg19_mm10_transgenes.fasta",
		"star_genome":	  "gs://fc-e0000000-0000-0000-0000-000000000000/human_mouse/STAR2_5_index_hg19_mm10.tar.gz",
		"gene_intervals":	 "gs://fc-e0000000-0000-0000-0000-000000000000/human_mouse/hg19_mm10_transgenes.genes.intervals",
		"genome_dict":	  "gs://fc-e0000000-0000-0000-0000-000000000000/human_mouse/hg19_mm10_transgenes.dict",
		"star_cpus": 32,
		"star_memory": "120G"
	}

The fields star_cpus and star_memory are optional and are used as the default cpus and memory for running STAR with your genome.


Outputs
^^^^^^^^

The pipeline outputs a list of google bucket urls containing one gene-count matrix per sample. Each gene-count matrix file produced by Drop-seq tools has the suffix 'dge.txt.gz', matrices produced by dropEst have the extension .rds.

.. _Drop-seq alignment cookbook: https://github.com/broadinstitute/Drop-seq/blob/master/doc/Drop-seq_Alignment_Cookbook.pdf
.. _McCarroll lab: http://mccarrolllab.org/dropseq-1/
.. _dropEst: https://github.com/hms-dbmi/dropEst
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/



Building a Custom Genome
==========================

The tool **dropseq_bundle** can be used to build a custom genome.
Please see the description of important inputs below.

.. list-table::
	:widths: 5 30
	:header-rows: 1

	* - Name
	  - Description
	* - fasta_file
	  - Array of fasta files. If more than one species, fasta and gtf files must be in the same order.
	* - gtf_file
	  - Array of gtf files. If more than one species, fasta and gtf files must be in the same order.
	* - genomeSAindexNbases
	  - Length (bases) of the SA pre-indexing string. Typically between 10 and 15. Longer strings will use much more memory, but allow faster searches. For small genomes, must be scaled down to min(14, log2(GenomeLength)/2 - 1)


dropseq_workflow Terra Release Notes
===================================================

**Version 4**

- Handle uncompressed fastq files as workflow input.
- Added optional prepare_fastq_disk_space_multiplier input.

**Version 3**

- Set default value for docker_registry input.

**Version 2**

- Added docker_registry input.

**Version 1**

- Renamed cumulus to cumulus
- Added use_bases_mask option when running bcl2fastq

**Version 18**

- Created a separate docker image for running bcl2fastq

**Version 17**

- Fixed bug that ignored WDL input star_flags (thanks to Carly Ziegler for reporting)
- Changed default value of star_flags to the empty string (Prior versions of the WDL incorrectly indicated that basic 2-pass mapping was done)

**Version 16**

- Use cumulus dockerhub organization
- Changed default dropEst version to 0.8.6

**Version 15**

- Added drop_deq_tools_prep_bam_memory and drop_deq_tools_dge_memory options

**Version 14**

- Fix for downloading files from user pays buckets

**Version 13**

- Set GCLOUD_PROJECT_ID for user pays buckets

**Version 12**

- Changed default dropEst memory from 52G to 104G

**Version 11**

- Updated formula for computing disk size for dropseq_count

**Version 10**

- Added option to specify merge_bam_alignment_memory and sort_bam_max_records_in_ram

**Version 9**

- Updated default drop_seq_tools_version from 2.2.0 to 2.3.0

**Version 8**

- Made additional options available for running dropEst

**Version 7**

- Changed default dropEst memory from 104G to 52G

**Version 6**

- Added option to run dropEst

**Version 5**

- Specify full version for bcl2fastq (2.20.0.422-2 instead of 2.20.0.422)

**Version 4**

- Fixed issue that prevented bcl2fastq from running

**Version 3**

- Set default run_bcl2fastq to false
- Create shortcuts for commonly used genomes

**Version 2**

- Updated QC report

**Version 1**

- Initial release
