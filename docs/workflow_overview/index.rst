Cumulus workflows overview
-----------------------------

Cumulus workflows are written in WDL_ language, and published on Dockstore_. Below is an overview of them:

.. list-table::
    :widths: 8 5 5 25
    :header-rows: 1

    * - Workflow
      - First Version
      - Date Added
      - Function
    * - `Cellranger <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger>`_
      - 0.1.0
      - 2018-07-27
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generating count matrix using cellranger count or cellranger-atac count, running cellranger vdj or feature-barcode extraction.
    * - `Spaceranger <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Spaceranger>`_
      - 1.2.0
      - 2021-01-19
      - Run Space Ranger tools to process spatial transcriptomics data, which includes extracting sequence reads using spaceranger mkfastq, and generating count matrix using spaceranger count.
    * - `STARsolo <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/STARsolo>`_
      - 1.2.0
      - 2021-01-19
      - Run STARsolo to generate gene-count matrices fro FASTQ files.
    * - `GeoMx_fastq_to_dcc <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/GeoMx_fastq_to_dcc>`_
      - 2.2.0
      - 2022-10-04
      - Run Nanostring GeoMx Digital Spatial NGS Pipeline, and convert FASTQ files into DCC files.
    * - `GeoMx_dcc_to_count_matrix <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/GeoMx_dcc_to_count_matrix>`_
      - 2.2.0
      - 2022-10-04
      - Take the DCC zip file from *GeoMxFastqToDCC* workflow, as well as other output of GeoMx DSP machine as the input, and generate an Area Of Interest (AOI) by probe count matrix with pathologists' annotation.
    * - `Demultiplexing <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Demultiplexing>`_
      - 0.3.0
      - 2018-10-24
      - Run tools (demuxEM, souporcell, or popscle) for cell-hashing/nucleus-hashing/genetic-pooling analysis.
    * - `Cumulus <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cumulus>`_
      - 0.1.0
      - 2018-07-27
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - `Cellbender <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/CellBender>`_
      - 2.1.0
      - 2022-07-13
      - Run CellBender tool to remove technical artifacts from high-throughput single-cell/single-nucleus RNA sequencing data.
    * - `Cellranger_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_create_reference>`_
      - 0.12.0
      - 2019-12-14
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - `Cellranger_atac_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_atac_create_reference>`_
      - 0.12.0
      - 2019-12-14
      - Run Cell Ranger tools to build scATAC-seq references.
    * - `cellranger_vdj_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_vdj_create_reference>`_
      - 0.12.0
      - 2019-12-14
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - `STARsolo_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/STARsolo_create_reference>`_
      - 2.0.0
      - 2022-03-14
      - Run STAR to build sc/snRNA-seq references for STARsolo count.
    * - `Cellranger_atac_aggr <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Cellranger_atac_aggr>`_
      - 0.13.0
      - 2020-02-07
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - `Smart-Seq2 <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Smart-Seq2>`_
      - 0.5.0
      - 2018-11-18
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files.
    * - `Smart-Seq2_create_reference <https://dockstore.org/workflows/github.com/lilab-bcb/cumulus/Smart-Seq2_create_reference>`_
      - 0.12.0
      - 2019-12-14
      - Generate user-customized genome references for SMART-Seq2 data.


.. toctree::
   :maxdepth: 1
   :hidden:

   broad_method_registry


.. _Dockstore: https://dockstore.org/organizations/lilab/collections/Cumulus
.. _WDL: https://openwdl.org/
