Bulk RNA-Seq
-------------

Run Bulk RNA-Seq Workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Follow the steps below to generate count matrices from bulk RNA-Seq data on Terra_. This WDL estimates expression levels using *RSEM*.

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


#. Create a Terra `data table`_

    Example::

        entity:sample_id  read1 read2
        sample-1  gs://fc-e0000000/data-1/sample1-1_L001_R1_001.fastq.gz    gs://fc-e0000000/data-1/sample-1_L001_R2_001.fastq.gz
        sample-2 gs://fc-e0000000/data-1/sample-2_L001_R1_001.fastq.gz  gs://fc-e0000000/data-1/sample-2_L001_R2_001.fastq.gz

    You are free to add more columns, but sample ids and URLs to fastq files are required.

#. Upload your TSV file to your workspace. Open the ``DATA`` tab on your workspace. Then click the upload button on left ``TABLE`` panel, and select the TSV file above. When uploading is done, you'll see a new data table with name "sample":


#. Import *bulk_rna_seq* workflow to your workspace. Then open ``bulk_rna_seq`` in the ``WORKFLOW`` tab. Select ``Run workflow(s) with inputs defined by data table``, and choose *sample* from the drop-down menu.



Inputs:
^^^^^^^

Please see the description of important inputs below. Note that required inputs are in bold.

.. list-table::
    :header-rows: 1
    :widths: 5 20 5

    * - Name
      - Description
      - Default
    * - **sample_name**
      - Sample name
      -
    * - **read1**
      - Array of URLs to read 1
      -
    * - **read2**
      - Array of URLs to read 2
      -
    * - **reference**
      - Reference to align reads to
         - Pre-created genome references:
            - "GRCh38_ens93filt" for human, genome version is GRCh38, gene annotation is generated using human Ensembl 93 GTF according to `cellranger mkgtf`_;
            - "GRCm38_ens93filt" for mouse, genome version is GRCm38, gene annotation is generated using mouse Ensembl 93 GTF according to `cellranger mkgtf`_;
         - Create a custom genome reference using `smartseq2_create_reference workflow <./smart_seq_2.html#custom-genome>`_, and specify its Google bucket URL here.
      -
    * - aligner
      - Which aligner to use for read alignment. Options are "hisat2-hca", "star" and "bowtie"
      - "star"
    * - output_genome_bam
      - Whether to output bam file with alignments mapped to genomic coordinates and annotated with their posterior probabilities.
      - false


---------------------------------

Outputs:
^^^^^^^^

.. list-table::
    :header-rows: 1
    :widths: 5 20

    * - Name
      - Description
    * - rsem_gene
      - RSEM gene expression estimation.
    * - rsem_isoform
      - RSEM isoform expression estimation.
    * - rsem_trans_bam
      - RSEM transcriptomic BAM.
    * - rsem_genome_bam
      - RSEM genomic BAM files if ``output_genome_bam`` is ``true``.
    * - rsem_time
      - RSEM execution time log.
    * - aligner_log
      - Aligner log.
    * - rsem_cnt
      - RSEM count.
    * - rsem_model
      - RSEM model.
    * - rsem_theta
      - RSEM theta.


.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/
.. _data table: https://support.terra.bio/hc/en-us/articles/360025758392
.. _cellranger mkgtf: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
