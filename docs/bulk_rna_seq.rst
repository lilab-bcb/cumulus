Bulk RNA-Seq
-------------

Run Bulk RNA-Seq Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the steps below to generate TPM matrices from bulk RNA-Seq data on Terra_. This WDL aligns reads using *STAR*
and estimates expression levels using *RSEM*.

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

	Please note that the columns in the CSV must be in the order shown below and does not contain a header line.
	The sample sheet provides the FASTQ files for each sample.
	

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

		sample-1,gs://fc-e0000000-0000-0000-0000-000000000000/data-1/sample1-1_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/data-1/sample-1_L001_R2_001.fastq.gz
		sample-2,gs://fc-e0000000-0000-0000-0000-000000000000/data-1/sample-2_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/data-1/sample-2_L001_R2_001.fastq.gz
		sample-1,gs://fc-e0000000-0000-0000-0000-000000000000/data-2/sample1-1_L001_R1_001.fastq.gz,gs://fc-e0000000-0000-0000-0000-000000000000/data-2/sample-1_L001_R2_001.fastq.gz


	Note that in this example, sample-1 was sequenced across two flowcells.


#. Upload your sample sheet to the workspace bucket.

    Example::

        gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import *bulk_rna_seq* workflow to your workspace.

    See the Terra documentation for `adding a workflow`_. The *bulk_rna_seq* workflow is under ``Broad Methods Repository`` with name "**cumulus/bulk_rna_seq**".

    Moreover, in the workflow page, click ``Export to Workspace...`` button, and select the workspace to which you want to export *bulk_rna_seq* workflow in the drop-down menu.

#. In your workspace, open ``bulk_rna_seq`` in ``WORKFLOWS`` tab. Select ``Run workflow with inputs defined by file paths`` as below

    .. image:: images/single_workflow.png

   and click ``SAVE`` button.


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
      - Sample Sheet
      - CSV describing sample name, read1, and (optionally) read2
      - "gs://fc-e0000000/sample_sheet.csv"
      -
    * - **reference**
      - Reference genome to align reads to. Either a gs URL to a tar.gz file or one of the following pre-created references:
          - "GRCh38-2020-A" for human, GRCh38 (GENCODE v32/Ensembl 98), gene annotation is generated according to `cellranger mkgtf`_;
          - "mm10-2020-A" for mouse, mm10 (GENCODE vM23/Ensembl 98), gene annotation is generated according to `cellranger mkgtf`_;

        - Create a custom genome reference using `bulk_rna_seq_create_reference <./bulk_rna_seq.html#custom-genome>`_, and specify its Google bucket URL here.
      - | "GRCh38-2020-A", or
        | "gs://fc-e0000000-0000-0000-0000-000000000000/rsem_ref.tar.gz"
      -
    * - output_genome_bam
      - Whether to output bam file with alignments mapped to genomic coordinates and annotated with their posterior probabilities.
      -
      - false
    * - qc_vars
      - Maps qc name to comma separated list of genes
      -
      - {"mitochrondrial":"MT-", "ribosome":"RPL,RPS"}
    * - extra_disk_space
      - Extra disk space for RSEM
      -
      - 5
    * - num_cpu
      - Number of cpus to request for one sample
      - 4
      - 4
    * - memory
      - Memory size string for RSEM
      - "32G"
      -
    * - disk_space_multiplier
      - Factor to multiply size of R1 and R2 by for RSEM
      - Float
      - 11

---------------------------------

Outputs:
^^^^^^^^

.. list-table::
    :widths: 5 5 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - output_count_matrix
      - File
      - Array containing aggregated TPM matrices.
    * - qc_report
      - File
      - File containing quality control statistics. Each file contains one line per sample and each line has three columns: Total reads, Alignment rate and Unique rate.
    * - rsem_gene
      - Array[Array[File]]
      - A 2D array of RSEM gene expression estimation files.
    * - rsem_isoform
      - Array[Array[File]]
      - A 2D array of RSEM isoform expression estimation files.
    * - rsem_trans_bam
      - Array[Array[File]]
      - A 2D array of RSEM transcriptomic BAM files.
    * - rsem_genome_bam
      - Array[Array[File]]
      - A 2D array of RSEM genomic BAM files if ``output_genome_bam`` is ``true``.
    * - rsem_time
      - Array[Array[File]]
      - A 2D array of RSEM execution time log files.
    * - aligner_log
      - Array[Array[File]]
      - A 2D array of aligner log files.
    * - rsem_cnt
      - Array[Array[File]]
      - A 2D array of RSEM count files.
    * - rsem_model
      - Array[Array[File]]
      - A 2D array of RSEM model files.
    * - rsem_theta
      - Array[Array[File]]
      - A 2D array of RSEM generated theta files.


The aggregated TPM matrix format is:

- The first line starts with ``"Gene"`` and then gives samples separated by tabs.
- Starting from the second line, each line describes one gene. 
  The first item in the line is the gene name and the rest items are TPM values of this gene for each sample.



Building a Custom Genome
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The tool **bulk_rna_seq_create_reference** can be used to build a custom genome.
Please see the description of important inputs below.

.. list-table::
	:widths: 5 30
	:header-rows: 1

	* - Name
	  - Description
	* - fasta
	  - fasta file.
	* - gtf
	  - gtf file.

.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _cellranger mkgtf: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
.. _Terra: https://app.terra.bio/

