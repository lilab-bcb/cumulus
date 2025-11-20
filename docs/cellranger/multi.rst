.. note::
  Cell Ranger will send anonymized telemetry data to 10x Genomics starting from v9.0. Here is the details on `Cell Ranger Pipeline Telemetry`_.

  This option has been turned off in this *cellranger_workflow*, thus **no data will be sent to 10x Genomics**.


The cellranger workflow supports processing data of 10x Flex and Sample Multiplexing type, as well as multiomics data.
Follow the corresponding sections below based on your data type:

Flex Gene Expression
++++++++++++++++++++++

This section covers preparing the sample sheet for Flex_ (previously named *Fixed RNA Profiling*) data.

1. *Sample* and *Link* column.

  *Sample* column is for specifying the name of each sample in your data. They must be unique to each other in the sample sheet.

  *Link* column is for specifying the name of your whole data, so that the workflow knows which samples should be put together to run ``cellranger multi``.
    * **Notice 1:** You should use a unique *Link* name for all samples belonging to the same data/experiment. Moreover, the *Link* name must be different from all *Sample* names.
    * **Notice 2:** If there is only a scRNA-Seq sample in the data, you don't need to specify *Link* name. Then the workflow would use its *Sample* name for the whole data.

2. *DataType*, *Reference*, and *AuxFile* column.

  For each sample, choose a data type from the table below, and prepare its corresponding auxiliary file if needed:

  .. list-table::
    :widths: 5 5 10 10
    :header-rows: 1

    * - DataType
      - Reference
      - AuxFile
      - Description
    * - **frp**
      - | Select one from prebuilt genome references in `scRNA-seq section`_,
        | or provide a cloud URI of a custom reference in ``.tar.gz`` format.
      - Path to a text file including the sample name to Flex probe barcode association (see an example below this table).
      - For RNA-Seq samples
    * - Choose one from: **citeseq**, **crispr**
      - No need to specify a reference
      - Path to its feature reference file of `10x Feature Reference`_ format. **Notice:** If multiple antibody capture samples, you need to combine feature barcodes used in all of them in one reference file.
      - For antibody capture samples:

        - ``citeseq``: For CITE-Seq samples.

        - ``crispr``: For Perturb-Seq samples. **Notice:** This data type used in Flex is supported only in Cell Ranger v8.0+.

  An example sample name to Flex probe barcode association file is the following (see `here <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-flex-multi-frp#example-configs>`_ for examples of different Flex experiment settings)::

    sample_id,probe_barcode_ids,description
    sample1,BC001,Control
    sample2,BC002,Treated

  The *description* column is optional, which specifies the description of the samples.

  .. note::
    In the sample name to Flex probe barcode file, the header line is optional. But if users don't specify this header line, the order of columns must be fixed as *sample_id*, *probe_barcode_ids*, and *description* (optional).

  Below is an example sample sheet for Flex data::

    Sample,Reference,Flowcell,DataType,AuxFile
    s1,GRCh38-2020-A,gs://my-bucket/s1_fastqs,frp,gs://my-bucket/s1_flex.csv

  Notice that *Link* column is not required for this case.

  An example sample sheet for a more complex Flex data::

    Link,Sample,Reference,Flowcell,DataType,AuxFile
    s2,s2_gex,GRCh38-2020-A,gs://my-bucket/s2_fastqs,frp,gs://my-bucket/s2_flex.csv
    s2,s2_citeseq,,gs://my-bucket/s2_fastqs,citeseq,gs://my-bucket/s2_fbc.csv
    s2,s2_crispr,,gs://my-bucket/s2_fastqs,crispr,gs://my-bucket/s2_fbc.csv

3. Flex Probe Set.

  Flex uses probes that target protein-coding genes in the human or mouse transcriptome. It's automatically determined by the genome reference specified by users for the scRNA-Seq sample by following the table below:

    .. list-table::
        :widths: 5 5 5
        :header-rows: 1

        * - Genome Reference
          - Probe Set
          - Cell Ranger version
        * - GRCh38-2024-A
          - `Flex_human_probe_v1.1`_
          - v9.0+
        * - GRCh38-2020-A
          - `Flex_human_probe_v1.0.1`_
          - v7.1+
        * - GRCm39-2024-A
          - `Flex_mouse_probe_v1.1`_
          - v9.0+
        * - mm10-2020-A
          - `Flex_mouse_probe_v1.0.1`_
          - v7.1+

  See `Flex probe sets overview`_ for details on these probe sets.

On Chip Multiplexing
+++++++++++++++++++++

This section covers preparing the sample sheet for `On-Chip Multiplexing`_ (OCM) data.

1. *Sample* and *Link* column.

  *Sample* column is for specifying the name of each sample in your data. They must be unique to each other in the sample sheet.

  *Link* column is for specifying the name of your whole data, so that the workflow knows which samples should be put together to run ``cellranger multi``.
  **Notice:** You should use a unique *Link* name for all samples belonging to the same data/experiment. Moreover, the *Link* name must be different from all *Sample* names.

2. *DataType*, *Reference*, and *AuxFile* column.

  For each sample, choose a data type from the table below, and prepare its corresponding auxiliary file if needed:

  .. list-table::
    :widths: 5 5 5 10
    :header-rows: 1

    * - DataType
      - Reference
      - AuxFile
      - Description
    * - **rna**
      - Select one from prebuilt genome references in `scRNA-seq section`_, or provide a cloud URI of a custom reference in ``.tar.gz`` format.
      - Path to a text file including the sample name to OCM barcode association (see an example below this table).
      - For RNA-Seq samples
    * - Choose one from: **vdj**, **vdj_t**, **vdj_b**, **vdj_t_gd**
      - Select one from prebuilt VDJ references in `Single-cell immune profiling`_ section.
      - *Optional*. For ``vdj_t_gd`` type samples only: path to a text file containing inner enrichment primers info. This is the ``inner-enrichment-primers`` option in `VDJ section`_ of Cell Ranger multi config CSV.
      - For each VDJ sample, choose one from the 4 provided VDJ data types:

        - ``vdj``: Leave the workflow to auto-detect.

        - ``vdj_t``: VDJ-T library for T-cell receptor sequences.

        - ``vdj_b``: VDJ-B library for B-cell receptor sequences.

        - ``vdj_t_gd``: VDJ-T-GD library for T-cell receptor enriched for gamma (TRG) and delta (TRD) chains. **Notice:** For such sample, A text file containing inner enrichment primers info must provided in *AuxFile* column.
    * - Choose one from: **citeseq**, **adt**
      - No need to specify a reference
      - Path to its feature reference file of `10x Feature Reference`_ format. **Notice:** If ``adt`` type, you need to combine feature barcodes of both CITE-Seq and Hashing modalities in one file.
      - For antibody capture samples:

        - ``citeseq``: For samples only containing CITE-Seq modality.

        - ``adt``: For samples containing both CITE-Seq and Hashing modalities.

  An example sample name to OCM barcode association file is the following::

    sample_id,ocm_barcode_ids,description
    sample1,OB1,Control
    sample2,OB2,Treated

  where *description* column is optional, which specifies the description of the samples.

.. note::
  In the sample name to OCM barcode file, the header line is optional. But if users don't specify this header line, the order of columns must be fixed as *sample_id*, *ocm_barcode_ids*, and *description* (optional).

Below is an example sample sheet for OCM::

  Sample,Reference,Flowcell,DataType,AuxFile,Link
  s1_gex,GRCh38-2020-A,gs://my-bucket/s1_fastqs,rna,gs://my-bucket/s1_ocm.csv,s1
  s1_vdj,GRCh38_vdj_v7.1.0,gs://my-bucket/s1_fastqs,vdj,,s1
  s1_adt,,gs://my-bucket/s1_fastqs,citeseq,gs://my-bucket/s1_fbc.csv,s1

In the case where there is only scRNA-Seq library in your data, the *Link* column is optional::

  Sample,Reference,Flowcell,DataType,AuxFile
  s2,GRCh38-2020-A,gs://my-bucket/s2_fastqs,rna,gs://my-bucket/s2_ocm.csv

----------

Hashing with Antibody Capture
+++++++++++++++++++++++++++++++

This section covers preparing the sample sheet for `non-OCM hashtag oligo`_ (HTO) data.

1. *Sample* and *Link* column.

  *Sample* column is for specifying the name of each sample in your data. They must be unique to each other in the sample sheet.

  *Link* column is for specifying the name of your whole data, so that the workflow knows which samples should be put together to run ``cellranger multi``.
  **Notice:** You should use a unique *Link* name for all samples belonging to the same data/experiment. Moreover, the *Link* name must be different from all *Sample* names.

2. *DataType*, *Reference*, and *AuxFile* column.

  For each sample, choose a data type from the table below, and prepare its corresponding auxiliary file if needed:

  .. list-table::
    :widths: 5 5 5 10
    :header-rows: 1

    * - DataType
      - Reference
      - AuxFile
      - Description
    * - **rna**
      - Select one from prebuilt genome references in `scRNA-seq section`_, or provide a cloud URI of a custom reference in ``.tar.gz`` format.
      - Path to a text file including the sample name to HTO barcode association (see an example below this table).
      - For RNA-Seq samples
    * - Choose one from: **vdj**, **vdj_t**, **vdj_b**, **vdj_t_gd**
      - Select one from prebuilt VDJ references in `Single-cell immune profiling`_ section.
      - *Optional*. For ``vdj_t_gd`` type samples only: path to a text file containing inner enrichment primers info. This is the ``inner-enrichment-primers`` option in `VDJ section`_ of Cell Ranger multi config CSV.
      - For each VDJ sample, choose one from the 4 provided VDJ data types:

        - ``vdj``: Leave the workflow to auto-detect.

        - ``vdj_t``: VDJ-T library for T-cell receptor sequences.

        - ``vdj_b``: VDJ-B library for B-cell receptor sequences.

        - ``vdj_t_gd``: VDJ-T-GD library for T-cell receptor enriched for gamma (TRG) and delta (TRD) chains. **Notice:** For such sample, A text file containing inner enrichment primers info must provided in *AuxFile* column.
    * - **hashing**
      - No need to specify a reference
      - Path to its feature reference file of `10x Feature Reference`_ format, which specifies the oligonucleotide sequences used in the data.
      - For antibody capture samples

An example sample name to HTO barcode association file is the following::

    sample_id,hashtag_ids,description
    sample1,TotalSeqB_Hashtag_1,Control
    sample2,CD3_TotalSeqB,Treated

where names in *hashtag_ids* column must be consistent with ``id`` column in the feature reference file. The *description* column is optional, which specifies the description of the samples.

.. note::
  In the sample name to HTO barcode file, the header line is optional. But if users don't specify this header line, the order of columns must be fixed as *sample_id*, *hashtag_ids*, and *description* (optional).

Below is an example sample sheet for HTO::

  Link,Sample,Reference,Flowcell,DataType,AuxFile
  s1,s1_gex,GRCh38-2020-A,gs://my-bucket/s1_fastqs,rna,gs://my-bucket/s1_hto.csv
  s1,s1_vdj,GRCh38_vdj_v7.1.0,gs://my-bucket/s1_fastqs,vdj,
  s1,s1_hto,,gs://my-bucket/s1_fastqs,hashing,gs://my-bucket/s1_fbc_ref.csv

Or if your data contain only scRNA-Seq and antibody capture libraries::

  Link,Sample,Reference,Flowcell,DataType,AuxFile
  s2,s2_gex,GRCh38-2020-A,gs://my-bucket/s2_fastqs,rna,gs://my-bucket/s2_hto.csv
  s2,s2_hto,,gs://my-bucket/s2_fastqs,hashing,gs://my-bucket/s2_fbc_ref.csv

------------

Cell Multiplexing with CMO (CellPlex)
+++++++++++++++++++++++++++++++++++++++

This section covers preparing the sample sheet for CellPlex data using `Cell Multiplexing Oligos`_ (CMO).

1. *Sample* and *Link* column.

  *Sample* column is for specifying the name of each sample in your data. They must be unique to each other in the sample sheet.

  *Link* column is for specifying the name of your whole data, so that the workflow knows which samples should be put together to run ``cellranger multi``.
  **Notice:** You should use a unique *Link* name for all samples belonging to the same data/experiment. Moreover, the *Link* name must be different from all *Sample* names.

2. *DataType*, *Reference*, and *AuxFile* column.

  For each sample, choose a data type from the table below, and prepare its corresponding auxiliary file if needed:

  .. list-table::
    :widths: 5 5 5 10
    :header-rows: 1

    * - DataType
      - Reference
      - AuxFile
      - Description
    * - **rna**
      - Select one from prebuilt genome references in `scRNA-seq section`_, or provide a cloud URI of a custom reference in ``.tar.gz`` format.
      - Path to a text file including the sample name to CMO barcode association (see an example below this table).
      - For RNA-Seq samples
    * - **cmo**
      - No need to specify a reference
      - *Optional*. If using custom CMOs, provide the path to their ``cmo-set`` reference file of `10x Feature Reference`_ format. See `here <https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi#cmo-ref>`_ for an example.
      - For CMO samples.
    * - **citeseq**
      - No need to specify a reference
      - Path to its feature reference file of `10x Feature Reference`_ format.
      - For CITE-Seq samples.

An example sample name to CMO barcode association file is the following::

    sample_id,cmo_ids,description
    sample1,CMO301,Control
    sample2,CMO302,Treated

If using a ``cmo-set`` reference file, the names in *cmo_ids* must be consistent with ``id`` column in the CMO reference file. The *description* column is optional, which specifies the description of the samples.

.. note::
  In the sample name to CMO barcode file, the header line is optional. But if users don't specify this header line, the order of columns must be fixed as *sample_id*, *cmo_ids*, and *description* (optional).

Below is an example sample sheet for CellPlex::

  Link,Sample,Reference,Flowcell,DataType,AuxFile
  s1,s1_gex,GRCh38-2020-A,gs://my-bucket/s1_fastqs,rna,gs://my-bucket/s1_cmo.csv
  s1,s1_cellplex,,gs://my-bucket/s1_fastqs,cmo,

Or if a CITE-Seq sample/library is also included in the data::

  Link,Sample,Reference,Flowcell,DataType,AuxFile
  s2,s2_gex,GRCh38-2020-A,gs://my-bucket/s2_fastqs,rna,gs://my-bucket/s2_cmo.csv
  s2,s2_cellplex,,gs://my-bucket/s2_fastqs,cmo,
  s2,s2_citeseq,,gs://my-bucket/s2_fastqs,citeseq,gs://my-bucket/s2_fbc.csv

--------------

Multiomics
++++++++++++

To analyze multiomics (GEX + CITE-Seq/CRISPR) data, prepare your sample sheet as follows:

1. *Link* column.

  A unique link name for all modalities of the same data

2. *Chemistry* column.

  The workflow supports all 10x assay configurations. The most widely used ones are listed below:

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
    * - **ARC-v1**
      - Gene Expression portion of 10x Multiome data

  Please refer to the section of ``--chemistry`` option in `Cell Ranger Command Line Arguments`_ for all other valid chemistry keywords.

3. *DataType* column.

  The following keywords are accepted for *DataType* column:

  .. list-table::
    :widths: 5 20
    :header-rows: 1

    * - DataType
      - Explanation
    * - **rna**
      - For scRNA-seq samples
    * - **citeseq**
      - For CITE-seq samples
    * - **crispr**
      - For 10x CRISPR samples

4. *AuxFile* column.

  Prepare your feature reference file in `10x Feature Reference`_ format.

  **Notice:** If multiple antibody samples are used, you need to merge them into one feature reference file, and assign it for each of the samples.

Below is an example sample sheet::

  Link,Sample,Reference,DataType,Flowcell,Chemistry,AuxFile
  sample_4,s4_gex,GRCh38-2020-A,rna,gs://my-bucket/s4_fastqs,auto,
  sample_4,s4_citeseq,,citeseq,gs://my-bucket/s4_fastqs,SC3Pv4,gs://my-bucket/s4_feature_ref.csv

Here, by specifying ``sample_4`` in *Link* column, the two modalities will be processed together. The output will be one subfolder named ``sample_4``.


Workflow Input
++++++++++++++++

All the sample multiplexing assays share the same workflow input settings. ``cellranger_workflow`` takes sequencing reads as input (FASTQ files, or TAR files containing FASTQ files), and runs ``cellranger multi``. Revalant workflow inputs are described below, with required inputs highlighted in **bold**:

.. list-table::
    :widths: 5 30 30 20
    :header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - **input_csv_file**
      - Sample Sheet (contains Link, Sample, Reference, DataType, Flowcell, and AuxFile columns)
      - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
      -
    * - **output_directory**
      - Output directory
      - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
      -
    * - include_introns
      - Turn this option on to also count reads mapping to intronic regions. With this option, users do not need to use pre-mRNA references
      - true
      - true
    * - no_bam
      - Turn this option on to disable BAM file generation
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
      - Cell Ranger version to use. Available versions: 10.0.0, 9.0.1, 8.0.1, 7.2.0.
      - "10.0.0"
      - "10.0.0"
    * - docker_registry
      - Docker registry to use for cellranger_workflow. Options:

        - "quay.io/cumulus" for images on Red Hat registry;

        - "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - acronym_file
      - | The link/path of an index file in TSV format for fetching preset genome references, probe set references, chemistry whitelists, etc. by their names.
        | Set an GS URI if running on GCP; an S3 URI if running on AWS; an absolute file path if running on HPC or local machines.
      - "s3://xxxx/index.tsv"
      - "gs://cumulus-ref/resources/cellranger/index.tsv"
    * - zones
      - Google cloud zones. For GCP Batch backend, the zones are automatically restricted by the Batch settings.
      - "us-central1-a us-west1-a"
      - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    * - num_cpu
      - Number of cpus to request per link
      - 32
      - 32
    * - memory
      - Memory size string to request per link
      - "120G"
      - "120G"
    * - multi_disk_space
      - Used by Flex and Sample Multiplexing data. Disk space in GB to request per link.
      - 1500
      - 1500
    * - count_disk_space
      - Only used by Multiomics data. Disk space in GB to request per link
      - 500
      - 500
    * - preemptible
      - Number of preemptible tries. This only works for GCP.
      - 2
      - 2
    * - awsQueueArn
      - The AWS ARN string of the job queue to be used. This only works for AWS.
      - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
      - ""

Workflow Output
+++++++++++++++++

All the sample multiplexing assays share the same workflow output structure. See the table below for important outputs:

.. list-table::
    :widths: 5 5 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - cellranger_multi.output_multi_directory
      - Array[String]
      - Flex and Sample Multiplexing output. A list of cloud URIs to output folders, one URI per link.
    * - cellranger_count_fbc.output_count_directory
      - Array[String]
      - Multiomics output. A list of cloud URIs to output folders, one URI per link.



.. _scRNA-seq section: ./index.html#single-cell-and-single-nucleus-rna-seq
.. _Single-cell immune profiling: ./index.html#single-cell-immune-profiling
.. _VDJ section: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-multi-config-csv-opts#vdj
.. _On-Chip Multiplexing: https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-3p-what-is-cellplex#on-chip
.. _non-OCM hashtag oligo: https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-3p-what-is-cellplex#antibody-capture
.. _Cell Multiplexing Oligos: https://www.10xgenomics.com/support/software/cell-ranger/latest/getting-started/cr-3p-what-is-cellplex#antibody-capture
.. _10x Feature Reference: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-feature-ref-csv
.. _Flex: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-flex-multi-frp
.. _Flex probe sets overview: https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-probe-sets-overview
.. _Cell Ranger Command Line Arguments: https://www.10xgenomics.com/support/software/cell-ranger/latest/resources/cr-command-line-arguments
.. _Flex_human_probe_v1.1: https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-human-transcriptome-probe-set-1-1
.. _Flex_human_probe_v1.0.1: https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-human-transcriptome-probe-set
.. _Flex_mouse_probe_v1.1: https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-mouse-transcriptome-probe-set-1-1
.. _Flex_mouse_probe_v1.0.1: https://www.10xgenomics.com/support/flex-gene-expression/documentation/steps/probe-sets/chromium-frp-mouse-transcriptome-probe-set
