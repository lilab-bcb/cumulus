Run Cell Ranger tools using cellranger_workflow
------------------------------------------------


Follow the steps below to run CellRanger mkfastq/count/vdj on Terra_.

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
	
#. Non Broad Institute users that wish to run cellrankger mkfastq must create a custom docker image that contains bcl2fastq.

    See :ref:`bcl2fastq-docker` instructions.

#. Create a scRNA-Seq formatted sample sheet. 

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to demultiplex flowcells and generate channel-specific count matrices. Note that *Sample*, *Lane*, and *Index* columns are defined exactly the same as in 10x's simple CSV layout file.

	scRNA-Seq formatted sample sheet description (required column headers are shown in bold):

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name.
		* - **Reference**
		  - 
			| Provides the reference genome used by *cellranger count* for each 10x channel. 
			| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as
			| **GRCh38_v3.0.0** for human GRCh38, cellranger reference 3.0.0, Ensembl v93 gene annotation
			| **hg19_v3.0.0** for human hg19, cellranger reference 3.0.0, Ensembl v87 gene annotation
			| **mm10_v3.0.0** for mouse mm10, cellranger reference 3.0.0, Ensembl v93 gene annotation
			| **GRCh38_and_mm10_v3.1.0** for human and mouse, built from GRCh38 and mm10 cellranger references (3.1.0), Ensembl v93 gene annotations for both human and mouse
			| **GRCh38_premrna_v1.2.0** (GRCh38_premrna, obsoleted) for human, introns included, built from GRCh38 cellranger reference 1.2.0, Ensembl v84 gene annotation, treating annotated transcripts as exons
			| **mm10_premrna_v1.2.0** (mm10_premrna, obsoleted) for mouse, introns included, built from mm10 cellranger reference 1.2.0, Ensembl v84 gene annotation, treating annotated transcripts as exons
			| **GRCh38_premrna_and_mm10_premrna_v1.2.0** (GRCh38_premrna_and_mm10_premrna, obsoleted) for human and mouse, introns included, built from GRCh38_premrna and mm10_premrna
			| **GRCh38_vdj_v3.1.0** for human V(D)J sequences, cellranger reference 3.1.0, annotation built from Ensembl *Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf*
			| **GRCm38_vdj_v3.1.0** for mouse V(D)J sequences, cellranger reference 3.1.0, annotation built from Ensembl *Mus_musculus.GRCm38.94.gtf*
			| **GRCh38_atac_v1.1.0** for scATAC-Seq, human GRCh38, cellranger-atac reference 1.1.0
			| **mm10_atac_v1.1.0** for scATAC-Seq, mouse mm10, cellranger-atac reference 1.1.0
			| **GRCh38_v1.2.0** (GRCh38, obsoleted) for human GRCh38, cellranger reference 1.2.0, Ensembl v84 gene annotation
			| **hg19_v1.2.0** (hg19, obsoleted) for human hg19, cellranger reference 1.2.0, Ensembl v82 gene annotation
			| **mm10_v1.2.0** (mm10, obsoleted) for mouse mm10, cellranger reference 1.2.0, Ensembl v84 gene annotation
			| **GRCh38_and_mm10_v1.2.0** (GRCh38_and_mm10, obsoleted) for human and mouse, built from GRCh38 and mm10 cellranger references (1.2.0), Ensembl v84 gene annotations for both human and mouse
			| **GRCh38_vdj_v2.0.0** (GRCh38_vdj, obsoleted) for human V(D)J sequences, cellranger reference 2.0.0, annotation built from Ensembl *Homo_sapiens.GRCh38.87.chr_patch_hapl_scaff.gtf* and *vdj_GRCh38_alts_ensembl_10x_genes-2.0.0.gtf*
			| **GRCm38_vdj_v2.0.0** (GRCm38_vdj, obsoleted) for mouse V(D)J sequences, cellranger reference 2.0.0, annotation built from Ensembl *Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf*.
		* - **Flowcell**
		  - Indicates the Google bucket URL of uploaded BCL folders.
		* - **Lane**
		  - Tells which lanes the sample was pooled into.
		* - **Index**
		  - 
			| Contains 10x sample index set names (e.g. SI-GA-A12). 
			| For cell-hashing/nucleus-hashing/CITE-Seq, this field should be the index sequence. 
			| Note that all index sequences (including 10x's) should have the same length (8 bases). If one index sequence is shorter (e.g. ATCACG), pad it with P7 sequence (e.g. ATCACGAT).
		* - Chemistry
		  - 
			| Describes the 10x chemistry used for the sample. 
			| This column is optional. The default chemistry for **rna** data (see DataType column below) is **auto**, which will try to detect the chemistry automatically. 
			| The default chemistry for **adt** and **crispr** is **sC3Pv3**, which stands for single cell 3' v3. 
			| Note that if the index read has extra bases besides cell barcode and UMI, autodetection might fail. In this case, please specify the chemistry.
			| According to *cellranger count*'s documentation, chemistry can be
			| **auto** for autodetection,
			| **threeprime** for Single Cell 3′,
			| **fiveprime** for Single Cell 5′,
			| **SC3Pv1** for Single Cell 3′ v1,
			| **SC3Pv2** for Single Cell 3′ v2,
			| **SC3Pv3** for Single Cell 3′ v3 (set cellranger version input parameter to >= 3.0.2),
			| **SC5P-PE** for Single Cell 5′ paired-end (both R1 and R2 are used for alignment),
			| **SC5P-R2** for Single Cell 5′ R2-only (where only R2 is used for alignment).
		* - DataType
		  - 
			| Describes the data type of the sample --- *count*, *vdj*, *adt*, or *crispr*. 
			| **rna** refers to gene expression data (*cellranger count*), 
			| **vdj** refers to V(D)J data (*cellranger vdj*), 
			| **adt** refers to antibody tag data, which can be either CITE-Seq, cell-hashing, or nucleus-hashing, 
			| **crispr** refers to Perturb-seq guide tag data,
			| **atac** refers to scATAC-Seq data (*cellranger-atac count*).
			| This column is optional and the default data type is *rna*.
		* - FeatureBarcodeFile
		  - Google bucket urls pointing to feature barcode files for *adt* and *crispr* data. Features can be either antibody for CITE-Seq, cell-hashing, nucleus-hashing or gRNA for Perburb-seq. This column is optional provided no *adt* or *crispr* data are in the sample sheet.

	The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 4 samples sequenced in two flowcells.

	Example::

		Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,rna
		sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,SC3Pv3,rna
		sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime,rna
		sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8,fiveprime,rna
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime,rna
		sample_2,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,SC3Pv3,rna
		sample_3,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime,rna
		sample_4,mm10,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8,fiveprime,rna
		sample_5,GRCh38_vdj,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ,1,SI-GA-A1,fiveprime,vdj
		sample_6,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ,2,AGATCCTT,SC3Pv3,adt,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv
		sample_7,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9ZZ,3,TCCGGAGA,threeprime,crispr,gs://fc-e0000000-0000-0000-0000-000000000000/crispr_index.csv
		sample_8,GRCh38_atac_v1.1.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9YB,*,SI-NA-A1,auto,atac



#. Upload your sample sheet to the workspace bucket.

	Example::

		gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/


#. Import cellranger_workflow tool.

	In Terra, select the ``Tools`` tab, then click ``Find a Tool``. Click ``Broad Methods Repository``. Type **cumulus/cellranger_workflow**.
 	You can also see the Terra documentation for `adding a tool`_.

#. Select ``Process single workflow from files``.

	.. image:: images/single_workflow.png

---------------------------------

cellranger_workflow inputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cellranger_workflow`` takes Illumina outputs as input and runs ``cellranger mkfastq``/``cellranger-atac mkfastq`` and ``cellranger count``/``cellranger vdj``/``cellranger-atac count``/cumulus feature extraction. Please see the description of inputs below. Note that required inputs are shown in bold.

.. list-table::
	:widths: 5 30 30 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_csv_file**
	  - Sample Sheet (contains Sample, Reference, Flowcell, Lane, Index as required and Chemistry, DataType, FeatureBarcodeFile as optional)
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
	  - 
	* - **output_directory**
	  - Output directory
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/cellranger_output"
	  -
	* - run_mkfastq
	  - If you want to run ``cellranger mkfastq`` or ``cellranger-atac mkfastq``
	  - true
	  - true
	* - run_count
	  - If you want to run steps after ``mkfastq``, such as ``cellranger count``, ``cellranger vdj``, ``cellranger-atac count`` or ``cumulus feature_extraction``
	  - true
	  - true
	* - delete_input_directory
	  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges 
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
	  - Perform cell ranger secondary analysis (dimensionality reduction, clustering, etc.)
	  - false
	  - false
	* - vdj_denovo
	  - Do not align reads to reference V(D)J sequences before de novo assembly
	  - true
	  - false
	* - vdj_chain
	  - Force the web summary HTML and metrics summary CSV to only report on a particular chain type. The accepted values are: auto for autodetection based on TR vs IG representation, TR for T cell receptors, IG for B cell receptors, all for all chain types
	  - TR
	  - 
	* - scaffold_sequence
	  - Scaffold sequence in sgRNA for Purturb-seq, only used for crispr data type
	  - "GTTTAAGAGCTAAGCTGGAA"
	  - 
	* - max_mismatch
	  - Maximum hamming distance in feature barcodes for the adt task
	  - 3
	  - 3
	* - min_read_ratio
	  - Minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for the adt task and crispr data type
	  - 0.1
	  - 0.1
	* - cellranger_version
	  - cellranger version, could be 2.2.0, 3.0.2, 3.1.0
	  - "3.0.2"
	  - "3.0.2"
	* - cellranger_atac_version
	  - cellranger-atac version, currently only 1.1.0
	  - "1.1.0"
	  - "1.1.0"
	* - cumulus_version
	  - Cumulus version for extracting feature barcode matrix, currently only "0.10.0" is available.
	  - "0.10.0"
	  - "0.10.0"
	* - zones
	  - Google cloud zones
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-east1-d us-west1-a us-west1-b"
	* - num_cpu
	  - Number of cpus to request for one node
	  - 64
	  - 64
	* - atac_num_cpu
	  - Number of cpus for cellranger-atac count
	  - 64
	  - 64
	* - memory
	  - Memory size string
	  - "128G"
	  - "128G"
	* - feature_memory
	  - Optional memory string for extracting feature count matrix
	  - "32G"
	  - "32G"
	* - atac_memory
	  - Memory string for cellranger-atac count
	  - "57.6G"
	  - "57.6G"
	* - mkfastq_disk_space
	  - Optional disk space in gigabytes for mkfastq
	  - 1500
	  - 1500
	* - count_disk_space
	  - Disk space in gigabytes needed for cellranger count
	  - 500
	  - 500
	* - vdj_disk_space
	  - Disk space in gigabytes needed for cellranger vdj
	  - 500
	  - 500
	* - feature_disk_space
	  - Disk space in gigabytes needed for extracting feature count matrix
	  - 100
	  - 100
	* - atac_disk_space
	  - Disk space in gigabytes needed for cellranger-atac count
	  - 500
	  - 500
	* - preemptible
	  - Number of preemptible tries
	  - 2
	  - 2

---------------------------------

cellranger_workflow outputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the table below for important *Cell Ranger mkfastq/count* outputs.


.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_fastqs_directory
	  - Array[String]
	  - A list of google bucket urls containing FASTQ files, one url per flowcell.
	* - output_count_directory
	  - Array[String]
	  - A list of google bucket urls containing count matrices, one url per sample.
	* - output_vdj_directory
	  - Array[String]
	  - A list of google bucket urls containing vdj results, one url per sample.
	* - output_adt_directory
	  - Array[String]
	  - A list of google bucket urls containing adt count matrices, one url per sample.
	* - output_atac_count_directory
	  - Array[String]
	  - A list of google bucket urls containing cellranger-atac count results, one url per sample.
	* - metrics_summaries
	  - File
	  - A excel spreadsheet containing QCs for each sample.
	* - output_web_summary
	  - Array[File]
	  - A list of htmls visualizing QCs for each sample (cellranger count output).
	* - count_matrix
	  - String
	  - gs url for a template count_matrix.csv to run cumulus.

---------------------------------

Only run the count part
^^^^^^^^^^^^^^^^^^^^^^^

Sometimes, users might want to perform demultiplexing locally and only run the count part on the cloud. This section describes how to only run the count part via ``cellranger_workflow``.


#. Copy your FASTQ files to the workspace using gsutil in your unix terminal. 

	You should upload folders of FASTQS. Each folder should contain all FASTQ files for one sample.

	Example::

		gsutil -m cp -r /foo/bar/fastq_path/K18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq


#. Create a scRNA-Seq formatted sample sheet.


	This sample sheet is the same as the sample sheet described for running mkfastq and count except that the full path to the FASTQ files is ``FlowCell/Sample``.

	scRNA-Seq formatted sample sheet description (required column headers are shown in bold):

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name.
		* - **Reference**
		  -
			| Provides the reference genome used by *cellranger count* for each 10x channel.
			| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as
			| **GRCh38_v3.0.0** for human GRCh38, cellranger reference 3.0.0, Ensembl v93 gene annotation
			| **hg19_v3.0.0** for human hg19, cellranger reference 3.0.0, Ensembl v87 gene annotation
			| **mm10_v3.0.0** for mouse mm10, cellranger reference 3.0.0, Ensembl v93 gene annotation
			| **GRCh38_and_mm10_v3.1.0** for human and mouse, built from GRCh38 and mm10 cellranger references (3.1.0), Ensembl v93 gene annotations for both human and mouse
			| **GRCh38_premrna_v1.2.0** (GRCh38_premrna, obsoleted) for human, introns included, built from GRCh38 cellranger reference 1.2.0, Ensembl v84 gene annotation, treating annotated transcripts as exons
			| **mm10_premrna_v1.2.0** (mm10_premrna, obsoleted) for mouse, introns included, built from mm10 cellranger reference 1.2.0, Ensembl v84 gene annotation, treating annotated transcripts as exons
			| **GRCh38_premrna_and_mm10_premrna_v1.2.0** (GRCh38_premrna_and_mm10_premrna, obsoleted) for human and mouse, introns included, built from GRCh38_premrna and mm10_premrna
			| **GRCh38_vdj_v3.1.0** for human V(D)J sequences, cellranger reference 3.1.0, annotation built from Ensembl *Homo_sapiens.GRCh38.94.chr_patch_hapl_scaff.gtf*
			| **GRCm38_vdj_v3.1.0** for mouse V(D)J sequences, cellranger reference 3.1.0, annotation built from Ensembl *Mus_musculus.GRCm38.94.gtf*
			| **GRCh38_atac_v1.1.0** for scATAC-Seq, human GRCh38, cellranger-atac reference 1.1.0
			| **mm10_atac_v1.1.0** for scATAC-Seq, mouse mm10, cellranger-atac reference 1.1.0
			| **GRCh38_v1.2.0** (GRCh38, obsoleted) for human GRCh38, cellranger reference 1.2.0, Ensembl v84 gene annotation
			| **hg19_v1.2.0** (hg19, obsoleted) for human hg19, cellranger reference 1.2.0, Ensembl v82 gene annotation
			| **mm10_v1.2.0** (mm10, obsoleted) for mouse mm10, cellranger reference 1.2.0, Ensembl v84 gene annotation
			| **GRCh38_and_mm10_v1.2.0** (GRCh38_and_mm10, obsoleted) for human and mouse, built from GRCh38 and mm10 cellranger references (1.2.0), Ensembl v84 gene annotations for both human and mouse
			| **GRCh38_vdj_v2.0.0** (GRCh38_vdj, obsoleted) for human V(D)J sequences, cellranger reference 2.0.0, annotation built from Ensembl *Homo_sapiens.GRCh38.87.chr_patch_hapl_scaff.gtf* and *vdj_GRCh38_alts_ensembl_10x_genes-2.0.0.gtf*
			| **GRCm38_vdj_v2.0.0** (GRCm38_vdj, obsoleted) for mouse V(D)J sequences, cellranger reference 2.0.0, annotation built from Ensembl *Mus_musculus.GRCm38.90.chr_patch_hapl_scaff.gtf*.
		* - **Flowcell**
		  - Indicates the Google bucket URL of uploaded FASTQ folders.
		* - Chemistry
		  -
			| Describes the 10x chemistry used for the sample.
			| This column is optional. The default chemistry for **rna** data (see DataType column below) is **auto**, which will try to detect the chemistry automatically.
			| The default chemistry for **adt** and **crispr** is **sC3Pv3**, which stands for single cell 3' v3.
			| Note that if the index read has extra bases besides cell barcode and UMI, autodetection might fail. In this case, please specify the chemistry.
			| According to *cellranger count*'s documentation, chemistry can be
			| **auto** for autodetection,
			| **threeprime** for Single Cell 3′,
			| **fiveprime** for Single Cell 5′,
			| **SC3Pv1** for Single Cell 3′ v1,
			| **SC3Pv2** for Single Cell 3′ v2,
			| **SC3Pv3** for Single Cell 3′ v3 (set cellranger version input parameter to >= 3.0.2),
			| **SC5P-PE** for Single Cell 5′ paired-end (both R1 and R2 are used for alignment),
			| **SC5P-R2** for Single Cell 5′ R2-only (where only R2 is used for alignment).
		* - DataType
		  -
			| Describes the data type of the sample --- *count*, *vdj*, *adt*, or *crispr*.
			| **rna** refers to gene expression data (*cellranger count*),
			| **vdj** refers to V(D)J data (*cellranger vdj*),
			| **adt** refers to antibody tag data, which can be either CITE-Seq, cell-hashing, or nucleus-hashing,
			| **crispr** refers to Perturb-seq guide tag data,
			| **atac** refers to scATAC-Seq data (*cellranger-atac count*).
			| This column is optional and the default data type is *rna*.
		* - FeatureBarcodeFile
		  - Google bucket urls pointing to feature barcode files for *adt* and *crispr* data. Features can be either antibody for CITE-Seq, cell-hashing, nucleus-hashing or gRNA for Perburb-seq. This column is optional provided no *adt* or *crispr* data are in the sample sheet.



	In the following example sample_1 is sequenced on 2 flowcells. The FASTQ files for flowcell_1 are located at gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_1/sample_1 while the FASTQ files for flowcell_2 are located at gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_2_sample1::

		Sample,Reference,Flowcell
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_1
		sample_1,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/flowcell_2

#. Set optional input ``run_mkfastq`` to ``false``.

---------------------------------

Extract feature count matrices from CITE-Seq/Cell-hashing/Nucleus-hashing/Perturb-seq assays
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``cellranger_workflow`` can optionally extract feature count matrices from *CITE-Seq/Cell-hashing/Nucleus-hashing/Perturb-seq* assays. For *CITE-Seq/Cell-hashing/Nucleus-hashing*, the feature refers to antibody. Note that for *CITE-Seq/Cell-hashing*, only Biolegend TotalSeq-A is supported. For *Perturb-seq*, the feature refers to guide RNA. To extract feature count matrices, please follow the instructions below.

Instructions to configure ``cellranger_workflow``
++++++++++++++++++++++++++++++++++++++++++++++++++++++

#. Prepare one feature barcode file per assay and upload the files to the Google bucket.

	Prepare a CSV file with the following format: feature_barcode,feature_name.
	See below for an example::

		TTCCTGCCATTACTA,sample_1
		CCGTACCTCATTGTT,sample_2
		GGTAGATGTCCTCAG,sample_3
		TGGTGTCATTCTTGA,sample_4

	The above file describes a cell-hashing application with 4 samples.

#. Add assay information into the sample sheet.

	See below for an example::

		Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
		sample_1_rna,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,count
		sample_1_adt,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,ATTACTCG,threeprime,adt,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv
		sample_2_adt,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,TCCGGAGA,SC3Pv3,adt,gs://fc-e0000000-0000-0000-0000-000000000000/antibody_index.csv
		sample_3_crispr,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,CGCTCATT,SC3Pv3,crispr,gs://fc-e0000000-0000-0000-0000-000000000000/crispr_index.csv

	In the above sample sheet, the first line describes the normal 3' RNA assay. The second line describes its associated antibody tag data, which can from either a CITE-Seq, cell-hashing, or nucleus-hashing experiment. Note that for the tag data, the *Index* field is different. The index for tag and crispr data should be Illumina index primer sequence (e.g. D701 in line two). In addition, the *DataType* field is changed to *adt*. The third line describes another tag data, which is in 10x genomics' V3 chemistry. For tag and crispr data, it is important to explicitly state the chemistry (e.g. *SC3Pv3*). The last line describes one gRNA guide data for Perturb-seq (see the *crispr* in *DataType* field). Note that it is users' responsibility to avoid index collision between 10x genomics' RNA indexes (e.g. SI-GA-A8) and Illumina index sequences for tag and crispr data (e.g. ATTACTCG).


#. Fill in the ADT-specific parameters:

	.. list-table::
		:widths: 5 30 30 5
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - scaffold_sequence
		  - Scaffold sequence in sgRNA for Purturb-seq, only used for crispr data type
		  - "GTTTAAGAGCTAAGCTGGAA"
		  - 
		* - max_mismatch
		  - Maximum hamming distance in feature barcodes for the adt task
		  - 3
		  - 3
		* - min_read_ratio
		  - Minimum read count ratio (non-inclusive) to justify a feature given a cell barcode and feature combination, only used for the adt task and crispr data type
		  - 0.1
		  - 0.1
		* - adt_memory
		  - Optional memory in GB for extracting ADT count matrix
		  - 32
		  - 32
		* - adt_disk_space
		  - Optional disk space needed for extracting ADT count matrix
		  - 100
		  - 100

Parameters used for feature count matrix extraction
+++++++++++++++++++++++++++++++++++++++++++++++++++

If the chemistry is V2, `10x genomics v2 cell barcode white list`_ will be used, a hamming distance of 1 is allowed for matching cell barcodes, and the UMI length is 10. 
If the chemistry is V3, `10x genomics v3 cell barcode white list`_ will be used, a hamming distance of 0 is allowed for matching cell barcodes, and the UMI length is 12.

For Perturb-seq data, a small number of sgRNA protospace sequences will be sequenced ultra-deeply and we may have PCR chimeric reads. Therefore, we generate filtered feature count matrices as well in a data driven manner: 

We first plot the histogram of UMIs with certain number of read counts. The number of UMIs with x supporting reads decreases when x increases. We start from x = 1 and if we find that count[x] < count[x + 1] < count[x + 2], we detect a valley between two peaks. We filter out all UMIs with < x supporting reads since they are likely formed due to chimeric reads. 

In addition, we also filter out barcode-feature-UMI combinations that have their read count ratio, which is defined as total reads supporting barcode-feature-UMI over total reads supporting barcode-UMI, no larger than *min_read_ratio*.

Extracted feature count matrix output
+++++++++++++++++++++++++++++++++++++

For each antibody tag or crispr tag sample, a folder with the sample ID is generated under ``cellranger_output_directory``. In the folder, two files --- ``sample_id.csv`` and ``sample_id.stat.csv.gz`` --- are generated.

``sample_id.csv`` is the feature count matrix. It has the following format. The first line describes the column names: ``Antibody/CRISPR,cell_barcode_1,cell_barcode_2,...,cell_barcode_n``. The following lines describe UMI counts for each feature barcode, with the following format: ``feature_name,umi_count_1,umi_count_2,...,umi_count_n``.

``sample_id.stat.csv.gz`` stores the gzipped sufficient statistics. It has the following format. The first line describes the column names: ``Barcode,UMI,Feature,Count``. The following lines describe the read counts for every barcode-umi-feature combination.

If data type is ``crispr``, three additional files, ``sample_id.umi_count.pdf``, ``sample_id.filt.csv`` and ``sample_id.filt.stat.csv.gz``, are generated.

``sample_id.umi_count.pdf`` plots number of UMIs against UMI with certain number of reads and colors UMIs with high likelihood of being chimeric in blue and other UMIs in red. This plot is generated purely based on number of reads each UMI has.

``sample_id.filt.csv`` is the filtered feature count matrix. It has the same format as ``sample_id.csv``.

``sample_id.filt.stat.csv.gz`` is the filtered sufficient statistics. It has the same format as ``sample_id.stat.csv.gz``.



.. _10x genomics v2 cell barcode white list: gs://regev-lab/resources/cellranger/737K-august-2016.txt.gz
.. _10x genomics v3 cell barcode white list: gs://regev-lab/resources/cellranger/3M-february-2018.txt.gz
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a tool: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/
