Cellranger multi supports `Flex`_ (previous "Fixed RNA Profiling") since version 7.0.0.

Sample Sheet
++++++++++++++

#. **Reference** column.

    Prebuilt scRNA-seq references for Flex data processing are summarized below.

    .. list-table::
        :widths: 5 20
        :header-rows: 1

        * - Keyword
          - Description
        * - **GRCh38-2024-A**
          - Human GRCh38 (GENCODE v44/Ensembl 110)
        * - **GRCh38-2020-A**
          - Human GRCh38 (GENCODE v32/Ensembl 98)
        * - **GRCm39-2024-A**
          - Mouse GRCm39 (GENCODE vM33/Ensembl 110)
        * - **mm10-2020-A**
          - Mouse mm10 (GENCODE vM23/Ensembl 98)

#. *DataType* column.

    Set ``frp`` for RNA-Seq modalities of your Flex samples. For other modalities (e.g. citeseq or antibody), set to their corresponding data types.

#. **ProbeSet** column.

    Preset Flex probe set references have their own compatible scRNA-seq references and Cell Ranger versions, respectively. This is summarized as follows:

    .. list-table::
        :widths: 5 5 5
        :header-rows: 1

        * - Keyword
          - Genome Reference
          - Cell Ranger version
        * - **Flex_human_probe_v1.1**
          - GRCh38-2024-A
          - v9.0+
        * - **Flex_human_probe_v1.0.1**
          - GRCh38-2020-A
          - v7.1+
        * - **Flex_human_probe_v1**
          - GRCh38-2020-A
          - v7.0+
        * - **Flex_mouse_probe_v1.1**
          - GRCm39-2024-A
          - v9.0+
        * - **Flex_mouse_probe_v1.0.1**
          - mm10-2020-A
          - v7.1+
        * - **Flex_mouse_probe_v1**
          - mm10-2020-A
          - v7.0+

    Custom probe set references are also accepted. Simply put the GS or S3 URI of the custom probe set CSV file for this column.

#. *FeatureBarcodeFile* column.

    Provide sample name - Probe Barcode association as follows::

        sample1,BC001|BC002,Control
        sample2,BC003|BC004,Treated

    where the third column (i.e. ``Control`` and ``Treated`` above) is optional, which specifies the description of the samples.

.. note::
  In the case of Singleplex Flex with Antibody Capture, for ``citeseq`` sample, the *FeatureBarcodeFile* you prepare should be in 10x format (see `here <https://cf.10xgenomics.com/samples/cell-exp/7.0.0/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex/10k_Human_PBMC_TotalSeqB_singleplex_Multiplex_count_feature_reference.csv>`_ for an example).

#. *Link* column.

    Put a sample unique link name for all modalities that are linked.

    If *Link* column is not set, only consider RNA-seq modalities (i.e. samples of *DataType* ``frp``) and use their *Sample* names as the *Link* names.

#. Example::

    Sample,Reference,ProbeSet,Flowcell,DataType,FeatureBarcodeFile,Link
    sample1,GRCh38-2020-A,Flex_human_probe_v1.0.1,/path/to/sample1/fastq/folder,frp,/path/to/sample1/fbf/file,
    sample2_rna,GRCh38-2020-A,Flex_human_probe_v1.0.1,/path/to/sample2/rna/fastq/folder,frp,/path/to/sample2/rna/fbf/file,sample2
    sample2_citeseq,GRCh38-2020-A,,/path/to/sample2/citeseq/fastq/folder,citeseq,/path/to/sample2/citeseq/fbf/file,sample2

In the example above, two linked samples are provided.


Workflow Input
++++++++++++++++

For Flex data, ``cellranger_workflow`` takes Illumina outputs as input and runs ``cellranger mkfastq`` and ``cellranger multi``. Revalant workflow inputs are described below, with required inputs highlighted in **bold**:

.. list-table::
    :widths: 5 30 30 20
    :header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - **input_csv_file**
      - Sample Sheet (contains Sample, Reference, DataType, Flowcell as required; Lane and Index are required if *run_mkfastq* is ``true``; ProbeSet, FeatureBarcodeFile and Link are optional)
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
      - If you want to run ``cellranger multi``
      - true
      - true
    * - delete_input_bcl_directory
      - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges
      - false
      - false
    * - mkfastq_barcode_mismatches
      - Number of mismatches allowed in matching barcode indices (bcl2fastq2 default is 1)
      - 0
      - 1
    * - mkfastq_force_single_index
      - If 10x-supplied i7/i5 paired indices are specified, but the flowcell was run with only one sample index, allow the demultiplex to proceed using the i7 half of the sample index pair
      - false
      - false
    * - mkfastq_filter_single_index
      - Only demultiplex samples identified by an i7-only sample index, ignoring dual-indexed samples. Dual-indexed samples will not be demultiplexed
      - false
      - false
    * - mkfastq_use_bases_mask
      - Override the read lengths as specified in RunInfo.xml
      - "“Y28n*,I8n*,N10,Y90n*”"
      -
    * - mkfastq_delete_undetermined
      - Delete undetermined FASTQ files generated by bcl2fastq2
      - false
      - false
    * - force_cells
      - Force pipeline to use this number of cells, bypassing the cell detection algorithm, mutually exclusive with expect_cells. This option is used by ``cellranger multi``.
      - 6000
      -
    * - expect_cells
      - Expected number of recovered cells. Mutually exclusive with force_cells. This option is used by ``cellranger multi``.
      - 3000
      -
    * - include_introns
      - Turn this option on to also count reads mapping to intronic regions. With this option, users do not need to use pre-mRNA references. Note that if this option is set, cellranger_version must be >= 5.0.0. This option is used by ``cellranger multi``.
      - true
      - true
    * - no_bam
      - Turn this option on to disable BAM file generation. This option is only available if cellranger_version >= 5.0.0. This option is used by ``cellranger multi``.
      - false
      - false
    * - secondary
      - Perform Cell Ranger secondary analysis (dimensionality reduction, clustering, etc.). This option is used by ``cellranger multi``.
      - false
      - false
    * - cellranger_version
      - Cell Ranger version to use. Available versions working for Flex data: 9.0.0, 8.0.1, 8.0.0, 7.2.0, 7.1.0, 7.0.1, 7.0.0.
      - "9.0.0"
      - "9.0.0"
    * - docker_registry
      - Docker registry to use for cellranger_workflow. Options:

        - "quay.io/cumulus" for images on Red Hat registry;

        - "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - mkfastq_docker_registry
      - Docker registry to use for ``cellranger mkfastq``. Default is the registry to which only Broad users have access. See :ref:`bcl2fastq-docker` for making your own registry.
      - "gcr.io/broad-cumulus"
      - "gcr.io/broad-cumulus"
    * - acronym_file
      - | The link/path of an index file in TSV format for fetching preset genome references, probe set references, chemistry whitelists, etc. by their names.
        | Set an GS URI if *backend* is ``gcp``; an S3 URI for ``aws`` backend; an absolute file path for ``local`` backend.
      - "s3://xxxx/index.tsv"
      - "gs://regev-lab/resources/cellranger/index.tsv"
    * - zones
      - Google cloud zones
      - "us-central1-a us-west1-a"
      - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    * - num_cpu
      - Number of cpus to request for one node for cellranger mkfastq and cellranger multi
      - 32
      - 32
    * - memory
      - Memory size string for cellranger mkfastq and cellranger multi
      - "120G"
      - "120G"
    * - mkfastq_disk_space
      - Optional disk space in GB for mkfastq
      - 1500
      - 1500
    * - multi_disk_space
      - Disk space in GB needed for cellranger multi
      - 1500
      - 1500
    * - backend
      - Cloud backend for file transfer and computation. Available options:

        - "gcp" for Google Cloud;
        - "aws" for Amazon AWS;
        - "local" for local machines.
      - "gcp"
      - "gcp"
    * - preemptible
      - Number of preemptible tries
      - 2
      - 2
    * - awsQueueArn
      - The AWS ARN string of the job queue to be used. This only works for ``aws`` backend.
      - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
      - ""

Workflow Output
+++++++++++++++++

See the table below for important outputs:

.. list-table::
    :widths: 5 5 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - fastq_outputs
      - Array[Array[String]]
      - ``fastq_outputs[0]`` gives the list of cloud urls containing FASTQ files for RNA-Seq modalities of Flex data, one url per flowcell.
    * - count_outputs
      - Map[String, Array[String]]
      - ``count_outputs["multi"]`` gives the list of cloud urls containing *cellranger multi* outputs, one url per sample.


.. _Flex: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-flex-multi-frp
