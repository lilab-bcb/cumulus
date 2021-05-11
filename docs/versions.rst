Latest and stable versions on Terra_
------------------------------------

Cumulus is a fast growing project. As a result, we frequently update WDL snapshot versions on Terra_.
See below for latest and stable WDL versions you can use.

Latest version
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `26 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/23>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generating count matrix using cellranger count or cellranger-atac count, running cellranger vdj or feature-barcode extraction.
    * - cumulus/spaceranger_workflow
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/spaceranger_workflow/3>`_
      - Run Space Ranger tools to process spatial transcriptomics data, which includes extracting sequence reads using spaceranger mkfastq, and generating count matrix using spaceranger count.
    * - cumulus/star_solo
      - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/star_solo/6>`_
      - Run STARsolo to generate gene-count matrices fro FASTQ files.
    * - cumulus/count
      - `18 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/18>`__
      - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
    * - cumulus/demultiplexing
      - `28 <https://portal.firecloud.org/?return=terra#methods/cumulus/demultiplexing/26>`_
      - Run tools (demuxEM, souporcell, or popscle) for cell-hashing/nucleus-hashing/genetic-pooling analysis.
    * - cumulus/cellranger_create_reference
      - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/10>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/5>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/3>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/4>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/10>`__
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files.
    * - cumulus/smartseq2_create_reference
      - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/10>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `38 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/38>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.


Stable version - v1.3.0
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `15 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/15>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generating count matrix using cellranger count or cellranger-atac count, running cellranger vdj or feature-barcode extraction.
    * - cumulus/spaceranger_workflow
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/spaceranger_workflow/1>`_
      - Run Space Ranger tools to process spatial transcriptomics data, which includes extracting sequence reads using spaceranger mkfastq, and generating count matrix using spaceranger count.
    * - cumulus/star_solo
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/star_solo/3>`_
      - Run STARsolo to generate gene-count matrices fro FASTQ files.
    * - cumulus/count
      - `18 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/18>`__
      - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
    * - cumulus/demultiplexing
      - `22 <https://portal.firecloud.org/?return=terra#methods/cumulus/demultiplexing/22>`_
      - Run tools (demuxEM, souporcell, or demuxlet) for cell-hashing/nucleus-hashing/genetic-pooling analysis.
    * - cumulus/cellranger_create_reference
      - `9 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/9>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/2>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/2>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/3>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`__
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files.
    * - cumulus/smartseq2_create_reference
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `36 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/36>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.


Stable version - v1.2.0
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `15 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/15>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generating count matrix using cellranger count or cellranger-atac count, running cellranger vdj or feature-barcode extraction.
    * - cumulus/spaceranger_workflow
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/spaceranger_workflow/1>`_
      - Run Space Ranger tools to process spatial transcriptomics data, which includes extracting sequence reads using spaceranger mkfastq, and generating count matrix using spaceranger count.
    * - cumulus/star_solo
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/star_solo/3>`_
      - Run STARsolo to generate gene-count matrices fro FASTQ files.
    * - cumulus/count
      - `18 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/18>`__
      - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
    * - cumulus/demultiplexing
      - `22 <https://portal.firecloud.org/?return=terra#methods/cumulus/demultiplexing/22>`_
      - Run tools (demuxEM, souporcell, or demuxlet) for cell-hashing/nucleus-hashing/genetic-pooling analysis.
    * - cumulus/cellranger_create_reference
      - `9 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/9>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/2>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/2>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/3>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`__
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files.
    * - cumulus/smartseq2_create_reference
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `35 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/34>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.


Stable version - v1.1.0
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `14 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/14>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/star_solo
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/star_solo/3>`_
      - Run STARsolo to generate gene-count matrices fro FASTQ files.
    * - cumulus/count
      - `16 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/16>`__
      - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
    * - cumulus/demultiplexing
      - `21 <https://portal.firecloud.org/?return=terra#methods/cumulus/demultiplexing/21>`_
      - Run tools (demuxEM, souporcell, or demuxlet) for cell-hashing/nucleus-hashing/genetic-pooling analysis.
    * - cumulus/cellranger_create_reference
      - `9 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/9>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/2>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/2>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/3>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`__
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/smartseq2_create_reference
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `34 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/34>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.

Stable version - v1.0.0
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `12 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/12>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/count
      - `14 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/14>`__
      - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
    * - cumulus/demultiplexing
      - `20 <https://portal.firecloud.org/?return=terra#methods/cumulus/demultiplexing/20>`_
      - Run tools (demuxEM, souporcell, or demuxlet) for cell-hashing/nucleus-hashing/genetic-pooling analysis.
    * - cumulus/cellranger_create_reference
      - `9 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/9>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/2>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/2>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/3>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`__
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/smartseq2_create_reference
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `31 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/31>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - cumulus/cumulus_hashing_cite_seq
      - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/10>`__
      - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis


Stable version - v0.15.0
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/10>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/count
      - `14 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/14>`__
      - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
    * - cumulus/cellranger_create_reference
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/8>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/2>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/2>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/2>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`__
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/smartseq2_create_reference
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `24 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/24>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - cumulus/cumulus_subcluster
      - `16 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/16>`__
      - Run subcluster analysis using cumulus
    * - cumulus/cumulus_hashing_cite_seq
      - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/10>`__
      - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis


Stable version - v0.14.0
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/8>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/count
      - `11 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/11>`__
      - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
    * - cumulus/cellranger_create_reference
      - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/6>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/1>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/1>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/1>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`__
      - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/smartseq2_create_reference
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `16 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/16>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - cumulus/cumulus_subcluster
      - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/10>`__
      - Run subcluster analysis using cumulus
    * - cumulus/cumulus_hashing_cite_seq
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/8>`__
      - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis


Stable version - v0.13.0
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/7>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/cellranger_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/1>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_aggr
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/1>`__
      - Run Cell Ranger tools to aggregate scATAC-seq samples.
    * - cumulus/cellranger_atac_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/1>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/1>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/5>`__
      - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/smartseq2_create_reference
      - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/4>`__
      - Generate user-customized genome references for SMART-Seq2 data.
    * - cumulus/cumulus
      - `14 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/14>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - cumulus/cumulus_subcluster
      - `9 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/9>`__
      - Run subcluster analysis using cumulus
    * - cumulus/cumulus_hashing_cite_seq
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/7>`__
      - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis


Stable version - v0.12.0
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/6>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/cellranger_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/1>`__
      - Run Cell Ranger tools to build sc/snRNA-seq references.
    * - cumulus/cellranger_atac_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/1>`__
      - Run Cell Ranger tools to build scATAC-seq references.
    * - cumulus/cellranger_vdj_create_reference
      - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/1>`__
      - Run Cell Ranger tools to build single-cell immune profiling references.
    * - cumulus/smartseq2
      - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/5>`__
      - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/smartseq2_create_reference
      - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/4>`__
      - Generate user-customized genome references for SMART-Seq2 workflow.
    * - cumulus/cumulus
      - `11 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/11>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - cumulus/cumulus_subcluster
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/8>`__
      - Run subcluster analysis using cumulus
    * - cumulus/cumulus_hashing_cite_seq
      - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/6>`__
      - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis


Stable version - v0.11.0
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/4>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/smartseq2
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/3>`__
      - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/cumulus
      - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/8>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - cumulus/cumulus_subcluster
      - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/5>`__
      - Run subcluster analysis using cumulus
    * - cumulus/cumulus_hashing_cite_seq
      - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/5>`__
      - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis

Stable version - v0.10.0
^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - cumulus/cellranger_workflow
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/3>`__
      - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
    * - cumulus/smartseq2
      - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/3>`__
      - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - cumulus/cumulus
      - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/7>`__
      - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
    * - cumulus/cumulus_subcluster
      - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/4>`__
      - Run subcluster analysis using cumulus
    * - cumulus/cumulus_hashing_cite_seq
      - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/4>`__
      - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis

Stable version - HTAPP v2
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - regev/cellranger_mkfastq_count
      - 45
      - Run Cell Ranger to extract FASTQ files and generate gene-count matrices for 10x genomics data
    * - scCloud/smartseq2
      - `5 <https://portal.firecloud.org/?return=terra#methods/scCloud/smartseq2/5>`__
      - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
    * - scCloud/scCloud
      - `14 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud/14>`__
      - Run scCloud analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering and more
    * - scCloud/scCloud_subcluster
      - `9 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud_subcluster/9>`__
      - Run subcluster analysis using scCloud
    * - scCloud/scCloud_hashing_cite_seq
      - `9 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud_hashing_cite_seq/9>`__
      - Run scCloud for cell-hashing/nucleus-hashing/CITE-Seq analysis

Stable version - HTAPP v1
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
    :widths: 15 5 30
    :header-rows: 1

    * - WDL
      - Snapshot
      - Function
    * - regev/cellranger_mkfastq_count
      - 39
      - Run Cell Ranger to extract FASTQ files and generate gene-count matrices for 10x genomics data
    * - scCloud/scCloud
      - `3 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud/3>`__
      - Run scCloud analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering and more

.. _Terra: https://app.terra.bio
