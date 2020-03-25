Latest and stable versions on Terra_
---------------------------------------

Cumulus is a fast growing project. As a result, we frequently update WDL snapshot versions on Terra_.
See below for latest and stable WDL versions you can use.

Latest version
^^^^^^^^^^^^^^^^

.. list-table::
	:widths: 15 5 30
	:header-rows: 1

	* - WDL
	  - Snapshot
	  - Function
	* - cumulus/cellranger_workflow
	  - `9 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/9>`_
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/count
	  - `13 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/13>`_
	  - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
	* - cumulus/cellranger_create_reference
	  - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/6>`_
	  - Run Cell Ranger tools to build sc/snRNA-seq references.
	* - cumulus/cellranger_atac_aggr
	  - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/2>`_
	  - Run Cell Ranger tools to aggregate scATAC-seq samples.
	* - cumulus/cellranger_atac_create_reference
	  - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/2>`_
	  - Run Cell Ranger tools to build scATAC-seq references.
	* - cumulus/cellranger_vdj_create_reference
	  - `2 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/2>`_
	  - Run Cell Ranger tools to build single-cell immune profiling references.
	* - cumulus/smartseq2
	  - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`_
	  - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/smartseq2_create_reference
	  - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`_
	  - Generate user-customized genome references for SMART-Seq2 data.
	* - cumulus/cumulus
	  - `18 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/18>`_
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - `12 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/12>`_
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/10>`_
	  - Run cumulus for cell-hashing/nucleus-hashing/CITE-Seq analysis

Stable version - v0.14.0
^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
	:widths: 15 5 30
	:header-rows: 1

	* - WDL
	  - Snapshot
	  - Function
	* - cumulus/cellranger_workflow
	  - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/8>`_
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/count
	  - `11 <https://portal.firecloud.org/?return=terra#methods/cumulus/count/11>`_
	  - Run alternative tools (STARsolo, Optimus, Salmon alevin, or Kallisto BUStools) to generate gene-count matrices from FASTQ files.
	* - cumulus/cellranger_create_reference
	  - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/6>`_
	  - Run Cell Ranger tools to build sc/snRNA-seq references.
	* - cumulus/cellranger_atac_aggr
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/1>`_
	  - Run Cell Ranger tools to aggregate scATAC-seq samples.
	* - cumulus/cellranger_atac_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/1>`_
	  - Run Cell Ranger tools to build scATAC-seq references.
	* - cumulus/cellranger_vdj_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/1>`_
	  - Run Cell Ranger tools to build single-cell immune profiling references.
	* - cumulus/smartseq2
	  - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/7>`_
	  - Run HISAT2/STAR/Bowtie2-RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/smartseq2_create_reference
	  - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/8>`_
	  - Generate user-customized genome references for SMART-Seq2 data.
	* - cumulus/cumulus
	  - `16 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/16>`_
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - `10 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/10>`_
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/8>`_
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
	  - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/7>`_
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/cellranger_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/1>`_
	  - Run Cell Ranger tools to build sc/snRNA-seq references.
	* - cumulus/cellranger_atac_aggr
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_aggr/1>`_
	  - Run Cell Ranger tools to aggregate scATAC-seq samples.
	* - cumulus/cellranger_atac_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/1>`_
	  - Run Cell Ranger tools to build scATAC-seq references.
	* - cumulus/cellranger_vdj_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/1>`_
	  - Run Cell Ranger tools to build single-cell immune profiling references.
	* - cumulus/smartseq2
	  - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/5>`_
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/smartseq2_create_reference
	  - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/4>`_
	  - Generate user-customized genome references for SMART-Seq2 data.
	* - cumulus/cumulus
	  - `14 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/14>`_
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - `9 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/9>`_
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/7>`_
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
	  - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/6>`_
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/cellranger_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_create_reference/1>`_
	  - Run Cell Ranger tools to build sc/snRNA-seq references.
	* - cumulus/cellranger_atac_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_atac_create_reference/1>`_
	  - Run Cell Ranger tools to build scATAC-seq references.
	* - cumulus/cellranger_vdj_create_reference
	  - `1 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_vdj_create_reference/1>`_
	  - Run Cell Ranger tools to build single-cell immune profiling references.
	* - cumulus/smartseq2
	  - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/5>`_
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/smartseq2_create_reference
	  - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2_create_reference/4>`_
	  - Generate user-customized genome references for SMART-Seq2 workflow.
	* - cumulus/cumulus
	  - `11 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/11>`_
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/8>`_
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - `6 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/6>`_
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
	  - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/4>`_
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/smartseq2
	  - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/3>`_
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/cumulus
	  - `8 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/8>`_
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/5>`_
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - `5 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/5>`_
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
	  - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/cellranger_workflow/3>`_
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/smartseq2
	  - `3 <https://portal.firecloud.org/?return=terra#methods/cumulus/smartseq2/3>`_
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/cumulus
	  - `7 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus/7>`_
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_subcluster/4>`_
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - `4 <https://portal.firecloud.org/?return=terra#methods/cumulus/cumulus_hashing_cite_seq/4>`_
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
	  - `5 <https://portal.firecloud.org/?return=terra#methods/scCloud/smartseq2/5>`_
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - scCloud/scCloud
	  - `14 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud/14>`_
	  - Run scCloud analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering and more
	* - scCloud/scCloud_subcluster
	  - `9 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud_subcluster/9>`_
	  - Run subcluster analysis using scCloud
	* - scCloud/scCloud_hashing_cite_seq
	  - `9 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud_hashing_cite_seq/9>`_
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
	  - `3 <https://portal.firecloud.org/?return=terra#methods/scCloud/scCloud/3>`_
	  - Run scCloud analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering and more

.. _Terra: https://app.terra.bio
