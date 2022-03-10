Extract gene-count matrices from plated-based SMART-Seq2 data
-------------------------------------------------------------

Run SMART-Seq2 Workflow
~~~~~~~~~~~~~~~~~~~~~~~~

Follow the steps below to extract gene-count matrices from SMART-Seq2 data on Terra_. This WDL aligns reads using *STAR*, *HISAT2*, or *Bowtie 2* and estimates expression levels using *RSEM*.

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

    Please note that the columns in the TSV can be in any order, but that the column names must match the recognized headings.

    The sample sheet provides metadata for each cell:

    .. list-table::
        :widths: 5 30
        :header-rows: 1

        * - Column
          - Description
        * - entity:sample
          - Cell name.
        * - plate
          - Plate name. Cells with the same plate name are from the same plate.
        * - read1
          - Location of the FASTQ file for read1 in the cloud (gsurl).
        * - read2
          - (Optional). Location of the FASTQ file for read2 in the cloud (gsurl). This field can be skipped for single-end reads.

    Example::

        entity:sample	plate	read1	read2
        cell-1	plate-1	gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-1_L001_R1_001.fastq.gz	gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-1_L001_R2_001.fastq.gz
        cell-2	plate-1	gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-2_L001_R1_001.fastq.gz	gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-2_L001_R2_001.fastq.gz
        cell-3	plate-2	gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-3_L001_R1_001.fastq.gz
        cell-4	plate-2	gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2/cell-4_L001_R1_001.fastq.gz


#. Upload your sample sheet to the workspace bucket.

    Example::

        gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import *smartseq2* workflow to your workspace.

    Import by following instructions in `Import workflows to Terra`_. You should choose **github.com/lilab-bcb/cumulus/Smart-Seq2** to import.

    Moreover, in the workflow page, click ``Export to Workspace...`` button, and select the workspace to which you want to export *smartseq2* workflow in the drop-down menu.

#. In your workspace, open ``smartseq2`` in ``WORKFLOWS`` tab. Select ``Run workflow with inputs defined by file paths`` as below

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
    * - **input_tsv_file**
      - Sample Sheet (contains entity:sample, plate, read1, read2)
      - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.tsv"
      -
    * - **output_directory**
      - Output directory
      - "gs://fc-e0000000-0000-0000-0000-000000000000/smartseq2_output"
      -
    * - **reference**
      - Reference transcriptome to align reads to. Acceptable values:

        - Pre-created genome references:

          - "GRCh38_ens93filt" for human, genome version is GRCh38, gene annotation is generated using human Ensembl 93 GTF according to `cellranger mkgtf`_;

          - "GRCm38_ens93filt" for mouse, genome version is GRCm38, gene annotation is generated using mouse Ensembl 93 GTF according to `cellranger mkgtf`_;

        - Create a custom genome reference using `smartseq2_create_reference workflow <./smart_seq_2.html#custom-genome>`_, and specify its Google bucket URL here.
      - | "GRCh38_ens93filt", or
        | "gs://fc-e0000000-0000-0000-0000-000000000000/rsem_ref.tar.gz"
      -
    * - aligner
      - Which aligner to use for read alignment. Options are "hisat2-hca", "star" and "bowtie"
      - "star"
      - "hisat2-hca"
    * - output_genome_bam
      - Whether to output bam file with alignments mapped to genomic coordinates and annotated with their posterior probabilities.
      - false
      - false
    * - smartseq2_version
      - SMART-Seq2 version to use. Versions available: 1.3.0.
      - "1.3.0"
      - "1.3.0"
    * - docker_registry
      - Docker registry to use. Options:

        - "quay.io/cumulus" for images on Red Hat registry;

        - "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - zones
      - Google cloud zones
      - "us-east1-d us-west1-a us-west1-b"
      - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    * - num_cpu
      - Number of cpus to request for one node
      - 4
      - 4
    * - memory
      - Memory size string
      - "3.60G"
      - If aligner is bowtie2 or hisat2-hca, "3.6G"; otherwise "32G"
    * - disk_space_multiplier
      - Factor to multiply size of R1 and R2 by for RSEM
      - Float
      - 11
    * - generate_count_matrix_disk_space
      - Disk space for count matrix generation task in GB
      - Integer
      - 10
    * - backend
      - Cloud infrastructure backend to use. Available options:

        - "gcp" for Google Cloud;
        - "aws" for Amazon AWS;
        - "local" for local machine.
      - "gcp"
      - "gcp"
    * - preemptible
      - Number of preemptible tries. This works only when *backend* is ``gcp``.
      - 2
      - 2
    * - awsMaxRetries
      - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
      - 5
      - 5

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
      - String
      - Point to a Google bucket URL for count matrix in matrix market format.
    * - rsem_trans_bam
      - Array[String?]
      - An array of Google bucket URLs for RSEM transcriptomic BAM files
    * - rsem_genome_bam
      - Array[String?]
      - An array of Google bucket URLs for RSEM genomic BAM files if ``output_genome_bam`` is ``true``.
    * - rsem_gene
      - Array[File?]
      - An array of RSEM gene expression estimation files.
    * - rsem_isoform
      - Array[File?]
      - An array of RSEM isoform expression estimation files.
    * - rsem_time
      - Array[File?]
      - An array of RSEM execution time log files.
    * - aligner_log
      - Array[File?]
      - An array of Aligner log files.
    * - rsem_cnt
      - Array[File?]
      - An array of RSEM count files.
    * - rsem_model
      - Array[File?]
      - An array of RSEM model files.
    * - rsem_theta
      - Array[File?]
      - An array of RSEM generated theta files.


This WDL generates one gene-count matrix in matrix market format:

- output_count_matrix is a folder containing three files: matrix.mtx.gz, barcodes.tsv.gz, and features.tsv.gz.
- matrix.mtx.gz is a gzipped matrix in matrix market format.
- barcodes.tsv.gz is a gzipped TSV file, containing 5 columns. 'barcodekey' is cell name. 'plate' is the plate name, which can be used for batch correction. 'total_reads' is the total number of reads. 'alignment_rate' is the alignment rate obtained from the aligner. 'unique_rate' is the percentage of reads aligned uniquely to a gene. Cells sequenced with single-end reads appear first in 'barcodekey'.
- features.tsv.gz is a gzipped TSV file, containing 2 columns. 'featurekey' is gene symbol. 'featureid' is Ensembl ID.

The gene-count matrix can be fed directly into **cumulus** for downstream analysis.

TPM-normalized counts are calculated as follows:

#. Estimate the gene expression levels in TPM using *RSEM*.

#. Suppose ``c`` reads are achieved for one cell, then calculate TPM-normalized count for gene ``i`` as ``TPM_i / 1e6 * c``.

TPM-normalized counts reflect both the relative expression levels and the cell sequencing depth.


---------------------------------

Custom Genome
~~~~~~~~~~~~~~~~

We also provide a way of generating user-customized Genome references for SMART-Seq2 workflow.

#. Import smartseq2_create_reference workflow to your workspace.

    Import by following instructions in `Import workflows to Terra`_. You should choose **github.com/lilab-bcb/cumulus/Smart-Seq2_create_reference** to import.

    Moreover, in the workflow page, click ``Export to Workflow...`` button, and select the workspace to which you want to export ``smartseq2_create_reference`` in the drop-down menu.

#. In your workspace, open ``smartseq2_create_reference`` in ``WORKFLOWS`` tab. Select ``Run workflow with inputs defined by file paths`` as below

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
      - Type or Example
      - Default
    * - **fasta**
      - Genome fasta file
      - | File.
        | For example, "gs://fc-e0000000-0000-0000-0000-000000000000/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
      -
    * - **gtf**
      - GTF gene annotation file (e.g. Homo_sapiens.GRCh38.83.gtf)
      - | File.
        | For example, "gs://fc-e0000000-0000-0000-0000-000000000000/Homo_sapiens.GRCh38.83.gtf"
      -
    * - **output_directory**
      - Google bucket url for the output folder
      - "gs://fc-e0000000-0000-0000-0000-000000000000/output_refs"
      -
    * - **genome**
      - Output reference genome name. Output reference is a gzipped tarball with name genome_aligner.tar.gz
      - "GRCm38_ens97filt"
      -
    * - aligner
      - Build indices for which aligner, choices are hisat2-hca, star, or bowtie2.
      - "hisat2-hca"
      - "hisat2-hca"
    * - smartseq2_version
      - | SMART-Seq2 version to use.
        | Versions available: 1.3.0.
      - "1.3.0"
      - "1.3.0"
    * - docker_registry
      - Docker registry to use. Options:

        - "quay.io/cumulus" for images on Red Hat registry;

        - "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - zones
      - Google cloud zones
      - "us-central1-c"
      - "us-central1-b"
    * - cpu
      - Number of CPUs
      - Integer
      - If aligner is bowtie2 or hisat2-hca, 8; otherwise 32
    * - memory
      - Memory size string
      - String
      - If aligner is bowtie2 or hisat2-hca, "7.2G"; otherwise "120G"
    * - disk_space
      - Disk space in GB
      - Integer
      - If aligner is bowtie2 or hisat2-hca, 40; otherwise 120
    * - backend
      - Cloud infrastructure backend to use. Available options:

        - "gcp" for Google Cloud;
        - "aws" for Amazon AWS;
        - "local" for local machine.
      - "gcp"
      - "gcp"
    * - preemptible
      - Number of preemptible tries. This works only when *backend* is ``gcp``.
      - 2
      - 2
    * - awsMaxRetries
      - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
      - 5
      - 5

Outputs
^^^^^^^^

.. list-table::
    :widths: 5 5 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - output_reference
      - File
      - The custom Genome reference generated. Its default file name is ``genome_aligner.tar.gz``.
    * - monitoring_log
      - File
      - CPU and memory profiling log.

---------------------------------


.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _Import workflows to Terra: ./cumulus_import.html
.. _cellranger mkgtf: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/advanced/references
.. _Terra: https://app.terra.bio/
