Cumulus workflows overview
-----------------------------

`Cumulus workflows`_ are written in WDL_ language, and published on Dockstore_. Below is an overview of them:

.. list-table::
    :widths: 8 5 30
    :header-rows: 1

    * - Workflow
      - Introduced in version
      - Function
    * - `Cellranger <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger>`_
      -
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generating count matrix using cellranger count or cellranger-atac count, running cellranger vdj or feature-barcode extraction.
    * - `Spaceranger <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Spaceranger>`_
      -
      - Run Space Ranger tools to process spatial transcriptomics data, which includes extracting sequence reads using spaceranger mkfastq, and generating count matrix using spaceranger count.
    * - `STARsolo <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/STARsolo>`_
      -
      - Run STARsolo to generate gene-count matrices fro FASTQ files.
    * - `STARsolo_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/STARsolo_create_reference>`_
      -
      -
    * - `Demultiplexing <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Demultiplexing>`_
      -
      - Run tools (demuxEM, souporcell, or popscle) for cell-hashing/nucleus-hashing/genetic-pooling analysis.
    * - `Cellranger_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_create_reference>`_
      -
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - `Cellranger_atac_aggr <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_atac_aggr>`_
      -
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - `Cellranger_atac_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_atac_create_reference>`_
      -
      - Run Cell Ranger tools to build scATAC-seq references.
    * - `cellranger_vdj_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_vdj_create_reference>`_
      -
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - `Smart-Seq2 <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Smart-Seq2>`_
      -
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files.
    * - `Smart-Seq2_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Smart-Seq2_create_reference>`_
      -
      - Generate user-customized genome references for SMART-Seq2 data.
    * - `Cumulus <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cumulus>`_
      -
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.

.. toctree::
   :maxdepth: 1
   :hidden:

   broad_method_registry


.. _Cumulus workflows: https://dockstore.org/organizations/lilab/collections/Cumulus
.. _WDL: https://openwdl.org/
.. _Dockstore: https://dockstore.org
