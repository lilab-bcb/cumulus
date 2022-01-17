Run STARsolo to generate gene-count matrices from FASTQ files
----------------------------------------------------------------------

This ``starsolo_workflow`` workflow generates gene-count matrices from FASTQ data using STARsolo.

Prepare input data and import workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Run ``starsolo_workflow`` to generate FASTQ data
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

You can skip this step if your data are already in FASTQ format.

Otherwise, for 10X data, you need to first run *cellranger_workflow* to generate FASTQ files from BCL raw data for each sample. Please follow `cellranger_workflow manual <./cellranger/index.html>`_.

Notice that you should set **run_mkfastq** to ``true`` to get FASTQ output. You can also set **run_count** to ``false`` to skip Cell Ranger count step.

For Non-Broad users, you'll need to build your own docker for *bcl2fastq* step. Instructions are `here <bcl2fastq.html>`_.

2. Import ``starsolo_workflow``
++++++++++++++++++++++++++++++++++

Import *starsolo_workflow* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose workflow **github.com/klarman-cell-observatory/cumulus/STARsolo** to import.

Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *starsolo_workflow* in the drop-down menu.

3. Prepare a sample sheet
++++++++++++++++++++++++++++

**3.1 Sample sheet format:**

Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

The sample sheet describes how to identify flowcells and generate sample/channel-specific count matrices.

A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

.. list-table::
    :widths: 5 30
    :header-rows: 1

    * - Column
      - Description
    * - **Sample**
      - Contains the sample name. Each sample should have a unique sample name.
    * - **Reference**
      - | Provides the reference genome used by STARSolo for each sample.
        | The elements in this column can be either Cloud bucket URIs to reference tarballs or keywords such as *GRCh38-2020-A*.
        | A full list of available keywords is included in `genome reference`_ section below.
    * - **Location**
      - Indicates the Cloud bucket URI of the folder holding FASTQ files of each sample.
    * - Assay
      - | Indicates the assay type of each sample.
        | Available options: ``tenX_v3`` for 10X v3, ``tenX_v2`` for 10X v2, ``DropSeq``, ``SeqWell``, ``SlideSeq``, ``ShareSeq`` and ``None``.
        | If not specified, use the default ``tenX_v3``.
        | If ``tenX_v3``, the following STARsolo options would be applied (could be overwritten by user-specified options): ``--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 12 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR``
        | ``--soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB``
        | If ``tenX_v2``, the following STARsolo options would be applied (could be overwritten by user-specified options): ``--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR``
        | ``--soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB``
        | If ``ShareSeq``, the following STARsolo options would be applied (could be overwritten by user-specific options): ``--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 24 --soloUMIstart 25 --soloUMIlen 10 --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts --soloUMIfiltering MultiGeneUMI_CR``
        | ``--soloUMIdedup 1MM_CR --clipAdapterType CellRanger4 --outFilterScoreMin 30 --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB``
        | If ``SeqWell`` or ``DropSeq``, the following STARsolo options would be applied (could be overwritten by user-specified options): ``--soloType CB_UMI_Simple --soloCBstart 1 --soloCBlen 12 --soloUMIstart 13 --soloUMIlen 8 --outSAMtype BAM SortedByCoordinate --outSAMattributes CR UR CY UY CB UB``
        | If ``None``, not preset options would be applied.

The sample sheet supports sequencing the same sample across multiple flowcells. In case of multiple flowcells, you should specify one line for each flowcell using the same sample name. In the following example, we have 2 samples and ``sample_1`` is sequenced in two flowcells.

Example::

    Sample,Location
    sample_1,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/sample_1_fastqs
    sample_1,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/sample_1_fastqs
    sample_2,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/sample_2_fastqs


**3.2 Upload your sample sheet to the workspace bucket:**

Example::

    gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/

4. Launch analysis
+++++++++++++++++++

In your workspace, open ``starsolo_workflow`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Process single workflow from files`` as below

    .. image:: images/single_workflow.png

and click ``SAVE`` button. Select ``Use call caching`` and click ``INPUTS``. Then fill in appropriate values in the ``Attribute`` column. Alternative, you can upload a JSON file to configure input by clicking ``Drag or click to upload json``.

Once INPUTS are appropriated filled, click ``RUN ANALYSIS`` and then click ``LAUNCH``.

----------------------------

Workflow inputs
^^^^^^^^^^^^^^^^^^

Below are inputs for *count* workflow. Notice that required inputs are in bold.

.. list-table::
    :widths: 5 20 10 5
    :header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - **input_csv_file**
      - Input CSV sample sheet describing metadata of each sample.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.tsv"
      -
    * - **output_directory**
      - Cloud bucket URI of output directory.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/count_result"
      -
    * - read1_fastq_pattern
      - | Filename suffix pattern in wildcards for Read 1. This is used for looking for Read 1 fastq files.
        | If fastq files are generated by CellRanger count, use ``_S*_L*_R1_001.fastq.gz``, which means Read 1 files must have names such as "<Sample>_S1_L1_R1_001.fastq.gz", where *<Sample>* is specified in **input_csv_file**.
        | If fastq files are Sequence Read Archive (SRA) data, use something like ``_1.fastq.gz``, where ``_1`` refers to the first reads, so that Read 1 files must have names such as "<Sample>_1.fastq.gz" where *<Sample>* is specified in **input_csv_file**.
        | If fastq files are not zipped, substitute ``.fastq`` for ``.fastq.gz`` in the corresponding pattern above.
      - "_S*_L*_R1_001.fastq.gz"
      - "_S*_L*_R1_001.fastq.gz"
    * - read2_fastq_pattern
      - | Filename suffix pattern in wildcards for Read 2. This is used for looking for Read 2 fastq files.
        | If fastq files are generated by CellRanger count, use ``_S*_L*_R2_001.fastq.gz``, which means Read 2 files must have names such as "<Sample>_S1_L1_R2_001.fastq.gz", where *<Sample>* is specified in **input_csv_file**.
        | If fastq files are Sequence Read Archive (SRA) data, use something like ``_2.fastq.gz``, where ``_2`` refers to the second reads, so that Read 2 files must have names such as "<Sample>_2.fastq.gz" where *<Sample>* is specified in **input_csv_file**.
        | If fastq files are not zipped, substitute ``.fastq`` for ``.fastq.gz`` in the corresponding pattern above.
      - "_S*_L*_R2_001.fastq.gz"
      - "_S*_L*_R2_001.fastq.gz"
    * - barcode_read
      - | Specify which read contains cell barcodes and UMIs: either ``read1`` or ``read2``. This only applies to samples with *Assay* ``None`` in **input_csv_file**.
        | Otherwise, samples with *Assay* type ``ShareSeq`` automatically specify ``read2`` for cell barcodes and UMIs, while ``read1`` for cDNAs;
        | samples of all the other know *Assay* types automatically specify ``read1`` for cell barcodes and UMIs, while ``read2`` for cDNAs.
      - "read1"
      - "read1"
    * - soloType
      - [STARsolo option] Type of single-cell RNA-seq, choosing from *CB_UMI_Simple*, *CB_UMI_Complex*, *CB_samTagOut*, *SmartSeq*.
      - "CB_UMI_Simple"
      -    None
    * - soloCBwhitelist
      - [STARsolo option] Cell barcode white list in either plain text or gzipped format.
      - gs://my_bucket/my_white_list.txt
      - None
    * - soloFeatures
      - [STARsolo option] Genomic features for which the UMI counts per Cell Barcode are collected (can choose multiple items):

        - *Gene*: reads match the gene transcript
        - *SJ*: splice junctions reported in SJ.out.tab
        - *GeneFull*: count all reads overlapping genes' exons and introns
        - *Velocyto*: calculate Spliced, Unspliced, and Ambiguous counts per cell per gene similar to the velocyto.py tool developed by LaManno et al. Note that *Velocyto* requires *Gene*.
      - "Gene GeneFull SJ Velocyto"
      - "Gene"
    * - soloMultiMappers
      - [STARsolo option] Counting method for reads mapping to multiple genes (can choose multiple items):

        - *Unique*: count only reads that map to unique genes
        - *Uniform*: uniformly distribute multi-genic UMIs to all genes
        - *Rescue*: distribute UMIs proportionally to unique+uniform counts (first iteartion of EM)
        - *PropUnique*: distribute UMIs proportionally to unique mappers, if present, and uniformly if not
        - *EM*: use Maximum Likelihood Estimation (MLE) to distribute multi-gene UMIs among their genes
      - "Unique"
      - "Unique"
    * - soloCBstart
      - [STARsolo option] Cell barcode start position (1-based coordinate).
      - 1
      - 1
    * - soloCBlen
      - [STARsolo option] Cell barcode length.
      - 16
      - 16
    * - soloUMIstart
      - [STARsolo option] UMI start position (1-based coordinate).
      - 17
      - 17
    * -    soloUMIlen
      - [STARsolo option] UMI length.
      - 10
      - 10
    * - soloBarcodeReadLength
      - [STARsolo option] Length of the barcode read

        - 1: equals to sum of *soloCBlen* and *soloUMIlen*.
        - 0: not defined, do not check.
      - 1
      - 1
    * - soloBarcodeMate
      - [STARsolo option] Identifies which read mate contains the barcode (CB+UMI) sequence:

        - 0: barcode sequence is on separate read, which should always be the last file in the input Read1 file list
        - 1: barcode sequence is a part of mate 1
        - 2: barcode sequence is a part of mate 2
      - 0
      - 0
    * - soloCBposition
      - | [STARsolo option] Position of Cell Barcode(s) on the barcode read.
        | Presently only works when *solo_type* is ``CB_UMI_Complex``, and barcodes are assumed to be on Read2.
        | Format for each barcode: "startAnchor_startPosition_endAnchor_endPosition"
        | start(end)Anchor defines the Anchor Base for the CB: 0: read start; 1: read end; 2: adapter start; 3: adapter end
        | start(end)Position is the 0-based position with of the CB start(end) with respect to the Anchor Base
        | String for different barcodes are separated by space.
      - "0\_0\_2\_-1 3\_1\_3\_8"
      -
    * - soloUMIposition
      - [STARsolo option] Position of the UMI on the barcode read, same as soloCBposition
      - "3\_9\_3\_14"
      -
    * - soloAdapterSequence
      - [STARsolo option] Adapter sequence to anchor barcodes.
      -
      -
    * - soloAdapterMismatchesNmax
      - [STARsolo option] Maximum number of mismatches allowed in adapter sequence.
      - 1
      - 1
    * - soloCBmatchWLtype
      - [STARsolo option] Matching the Cell Barcodes to the WhiteList, choosing from

        - *Exact*: only exact matches allowed
        - *1MM*: only one match in whitelist with 1 mismatched base allowed. Allowed CBs have to have at least one read with exact match
        - *1MM_multi*: multiple matches in whitelist with 1 mismatched base allowed, posterior probability calculation is used choose one of the matches. Allowed CBs have to have at least one read with exact match. This option matches best with CellRanger 2.2.0
        - *1MM_multi_pseudocounts*: same as *1MM_multi*, but pseudocounts of 1 are added to all whitelist barcodes
        - *1MM_multi_Nbase_pseudocounts*: same as *1MM_multi_pseudocounts*, multimatching to WL is allowed for CBs with N-bases. This option matches best with CellRanger >= 3.0.0
      - "1MM_multi"
      - "1MM_multi"
    * - soloInputSAMattrBarcodeSeq
      - [STARsolo option] When inputting reads from a SAM file (``--readsFileType SAM SE/PE``), these SAM attributes mark the barcode qualities (in proper order). For instance, for 10X CellRanger or STARsolo BAMs, use ``--soloInputSAMattrBarcodeSeq CR UR``. This parameter is required when running STARsolo with input from SAM.
      - "CR UR"
      -
    * - soloInputSAMattrBarcodeQual
      - [STARsolo option] When inputting reads from a SAM file (``--readsFileType SAM SE/PE``), these SAM attributes mark the barcode sequence (in proper order). For instance, for 10X CellRanger or STARsolo BAMs, use ``--soloInputSAMattrBarcodeQual CY UY``. If this parameter is ``-`` (default), the quality 'H' will be assigned to all bases.
      - "CY UY"
      -
    * - soloStrand
      - [STARsolo option] Strandedness of the solo libraries:

        - *Unstranded*: no strand information
        - *Forward*: read strand same as the original RNA molecule
        - *Reverse*: read strand opposite to the original RNA molecule
      - "Forward"
      - "Forward"
    * - soloUMIdedup
      - [STARsolo option] Type of UMI deduplication (collapsing) algorithm:

        - *1MM_All*: all UMIs with 1 mismatch distance to each other are collapsed (i.e. counted once)
        - *1MM Directional UMItools*: follows the "directional" method from the UMI-tools by Smith, Heger and Sudbery (Genome Research 2017)
        - *1MM Directional*: same as 1MM Directional UMItools, but with more stringent criteria for duplicate UMIs
        - *Exact*: only exactly matching UMIs are collapsed
        - *NoDedup*: no deduplication of UMIs, count all reads
        - *1MM CR*: CellRanger2-4 algorithm for 1MM UMI collapsing
      - "1MM_All"
      - "1MM_All"
    * - soloUMIfiltering
      - [STARsolo option] Type of UMI filtering (for reads uniquely mapping to genes):

        - *-*: basic filtering: remove UMIs with N and homopolymers (similar to CellRanger 2.2.0)
        - *MultiGeneUMI*: basic + remove lower-count UMIs that map to more than one gene
        - *MultiGeneUMI_All*: basic + remove all UMIs that map to more than one gene
        - *MultiGeneUMI_CR*: basic + remove lower-count UMIs that map to more than one gene, matching CellRanger > 3.0.0. Only works with ``--soloUMIdedup 1MM CR``
      - "MultiGeneUMI"
      - "-"
    * - soloCellFilter
      - [STARsolo option] Cell filtering type and parameters:

        - *None*: do not output filtered cells
        - *TopCells*: only report top cells by UMI count, followed by the exact number of cells
        - *CellRanger2.2*: simple filtering of CellRanger 2.2. Can be followed by numbers: number of expected cells, robust maximum percentile for UMI count, maximum to minimum ratio for UMI count. The harcoded values are from CellRanger: nExpectedCells=3000; maxPercentile=0.99; maxMinRatio=10
        - *EmptyDrops CR*: EmptyDrops filtering in CellRanger flavor. Please cite the original EmptyDrops paper: A.T.L Lun et al, Genome Biology, 20, 63 (2019): https://genomebiology.biomedcentral.com/articles/10.1186/s13059-019-1662-y. Can be followed by 10 numeric parameters: nExpectedCells maxPercentile maxMinRatio indMin indMax umiMin umiMinFracMedian candMaxN FDR simN. The harcoded values are from CellRanger: 3000 0.99 10 45000 90000 500 0.01 20000 0.01 10000
      - "CellRanger2.2 3000 0.99 10"
      - "CellRanger2.2 3000 0.99 10"
    * - soloOutFormatFeaturesGeneField3
      - [STARsolo option] Field 3 in the Gene features.tsv file. If "-", then no 3rd field is output.
      - "Gene Expression"
      - "Gene Expression"
    * - outSAMtype
      - [STAR option] Type of SAM/BAM output.
      - "BAM SortedByCoordinate"
      - | "BAM SortedByCoordinate" for *tenX_v3*, *tenX_v2*, *SeqWell* and *DropSeq* assay types,
        | "BAM Unsorted" otherwise.
    * - star_version
      - STAR version to use. Currently support: ``2.7.9a``.
      - "2.7.9a"
      - "2.7.9a"
    * - docker_registry
      - Docker registry to use:

        - ``quay.io/cumulus`` for images on Red Hat registry;

        - ``cumulusprod`` for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - zones
      - Google cloud zones to consider for execution.
      - "us-east1-d us-west1-a us-west1-b"
      - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    * - num_cpu
      - Number of CPUs to request for count per sample.
      - 32
      - 32
    * - memory
      - Memory size string for count per sample.
      - "120G"
      - "120G"
    * - disk_space
      - Disk space in GB needed for count per sample.
      - 500
      - 500
    * - backend
      - Cloud infrastructure backend to use. Available options:

        - ``gcp`` for Google Cloud;
        - ``aws`` for Amazon AWS;
        - ``local`` for local machine.
      - "gcp"
      - "gcp"
    * - preemptible
      - Number of maximum preemptible tries allowed. This works only when *backend* is ``gcp``.
      - 2
      - 2
    * - awsMaxRetries
      - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
      - 5
      - 5

Workflow outputs
^^^^^^^^^^^^^^^^^^^

See the table below for *star_solo* workflow outputs.

.. list-table::
    :widths: 5 5 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - output_folder
      - String
      - Google Bucket URL of output directory. Within it, each folder is for one sample in the input sample sheet.

----------------------------

Prebuilt genome references
^^^^^^^^^^^^^^^^^^^^^^^^^^^

We've built the following scRNA-seq references for users' convenience:

.. list-table::
    :widths: 5 20
    :header-rows: 1

    * - Keyword
      - Description
    * - **GRCh38-2020-A**
      - Human GRCh38, comparable to cellranger reference 2020-A (GENCODE v32/Ensembl 98)
    * - **mm10-2020-A**
      - Mouse mm10, comparable to cellranger reference 2020-A (GENCODE vM23/Ensembl 98)
    * - **GRCh38-and-mm10-2020-A**
      - Human GRCh38 (GENCODE v32/Ensembl 98) and mouse mm10 (GENCODE vM23/Ensembl 98)
    * - **GRCh38**
      - Human GRCh38, comparable to cellranger reference 3.0.0, Ensembl v93 gene annotation
    * - **mm10**
      - Mouse mm10, comparable to cellranger reference 3.0.0, Ensembl v93 gene annotation

We've built the following snRNA-seq references for users' convenience:

.. list-table::
    :widths: 5 20
    :header-rows: 1

    * - Keyword
      - Description
    * - **GRCh38-2020-A-premrna**
      - Human, introns included, built from GRCh38 cellranger reference 2020-A, GENCODE v32/Ensembl 98 gene annotation, treating annotated transcripts as exons
    * - **mm10-2020-A-premrna**
      - Mouse, introns included, built from mm10 cellranger reference 2020-A, GENCODE vM23/Ensembl 98 gene annotation, treating annotated transcripts as exons

---------------------------

Build STARSolo References
^^^^^^^^^^^^^^^^^^^^^^^^^^

We provide a wrapper of STAR to build sc/snRNA-seq references. Please follow the instructions below.

1. Import ``starsolo_create_reference``
+++++++++++++++++++++++++++++++++++++++++

Import *starsolo_create_reference* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose **github.com/klarman-cell-observatory/STARsolo_create_reference** to import.

Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *starsolo_create_reference* workflow in the drop-down menu.

2. Upload required data to Cloud bucket
++++++++++++++++++++++++++++++++++++++++++

Required data include the genome FASTA file and gene annotation GTF file of the target genome reference.

3. Workflow input
+++++++++++++++++++

Required inputs are highlighted **in bold**.

.. list-table::
    :widths: 5 20 10 5
    :header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - **input_fasta**
      - Input genome reference in FASTA format.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/mm-10/genome.fa"
      -
    * - **input_gtf**
      - Input gene annotation file in GTF format.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/mm-10/genes.gtf"
      -
    * - **genome**
      - Genome reference name. This is used for specifying the name of the genome index generated.
      - "mm-10"
      -
    * - **output_directory**
      - Cloud bucket URI of the output directory.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/starsolo-reference"
      -
    * - docker_registry
      - Docker registry to use:

        - ``quay.io/cumulus`` for images on Red Hat registry;

        - ``cumulusprod`` for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - star_version
      - STAR version to use. Currently support: ``2.7.9a``.
      - "2.7.9a"
      - "2.7.9a"
    * - num_cpu
      - Number of CPUs to request for count per sample.
      - 32
      - 32
    * - memory
      - Memory size string for count per sample.
      - "80G"
      - "80G"
    * - disk_space
      - Disk space in GB needed for count per sample.
      - 100
      - 100
    * - zones
      - Google cloud zones to consider for execution.
      - "us-east1-d us-west1-a us-west1-b"
      - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    * - backend
      - Cloud infrastructure backend to use. Available options:

        - ``gcp`` for Google Cloud;
        - ``aws`` for Amazon AWS;
        - ``local`` for local machine.
      - "gcp"
      - "gcp"
    * - preemptible
      - Number of maximum preemptible tries allowed. This works only when *backend* is ``gcp``.
      - 2
      - 2
    * - awsMaxRetries
      - Number of maximum retries when running on AWS. This works only when *backend* is ``aws``.
      - 5
      - 5

4. Workflow Output
+++++++++++++++++++

.. list-table::
    :widths: 2 2 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - output_reference
      - File
      - Gzipped reference folder with name **"<genome>-starsolo.tar.gz"**, where *<genome>* is specified by workflow input **genome** above. The workflow will save a copy of it under **output_directory** specified in workflow input above.

.. _Import workflows to Terra: ./cumulus_import.html
.. _genome reference: ./starsolo.html#prebuilt-genome-references
