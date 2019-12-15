Latest and stable versions on Terra_
---------------------------------------

Cumulus is a fast growing project. As a result, we frequently update WDL snapshot versions on Terra_.
See below for latest and stable WDL versions you can use.

Latest version
^^^^^^^^^^^^^^^

.. list-table::
	:widths: 15 5 30
	:header-rows: 1

	* - WDL
	  - Snapshot
	  - Function
	* - cumulus/cellranger_workflow
	  - 5
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/cellranger_create_reference
	  - 1
	  - Run Cell Ranger tools to build sc/snRNA-seq references.
	* - cumulus/cellranger_atac_create_reference
	  - 1
	  - Run Cell Ranger tools to build scATAC-seq references.
	* - cumulus/cellranger_vdj_create_reference
	  - 1
	  - Run Cell Ranger tools to build single-cell immune profiling references.
	* - cumulus/smartseq2
	  - 4
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/smartseq2_create_reference
	  - 4
	  - Generate user-customized genome references for SMART-Seq2 workflow.
	* - cumulus/cumulus
	  - 11
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - 8
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - 6
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
	  - 5
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/cellranger_create_reference
	  - 1
	  - Run Cell Ranger tools to build sc/snRNA-seq references.
	* - cumulus/cellranger_atac_create_reference
	  - 1
	  - Run Cell Ranger tools to build scATAC-seq references.
	* - cumulus/cellranger_vdj_create_reference
	  - 1
	  - Run Cell Ranger tools to build single-cell immune profiling references.
	* - cumulus/smartseq2
	  - 4
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/smartseq2_create_reference
	  - 4
	  - Generate user-customized genome references for SMART-Seq2 workflow.
	* - cumulus/cumulus
	  - 11
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - 8
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - 6
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
	  - 4
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/smartseq2
	  - 3
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/cumulus
	  - 8
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - 5
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - 5
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
	  - 3
	  - Run Cell Ranger tools, which include extracting sequence reads using cellranger mkfastq or cellranger-atac mkfastq, generate count matrix using cellranger count or cellranger-atac count, run cellranger vdj or feature-barcode extraction
	* - cumulus/smartseq2
	  - 3
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - cumulus/cumulus
	  - 7
	  - Run cumulus analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering, visualization, differential expression analysis, cell type annotation, etc.
	* - cumulus/cumulus_subcluster
	  - 4
	  - Run subcluster analysis using cumulus
	* - cumulus/cumulus_hashing_cite_seq
	  - 4
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
	  - 5
	  - Run Bowtie2 and RSEM to generate gene-count matrices for SMART-Seq2 data from FASTQ files
	* - scCloud/scCloud
	  - 14
	  - Run scCloud analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering and more
	* - scCloud/scCloud_subcluster
	  - 9
	  - Run subcluster analysis using scCloud
	* - scCloud/scCloud_hashing_cite_seq
	  - 9
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
	  - 3
	  - Run scCloud analysis module for variable gene selection, batch correction, PCA, diffusion map, clustering and more

.. _Terra: https://app.terra.bio
