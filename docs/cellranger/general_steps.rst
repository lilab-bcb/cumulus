This section mainly considers jobs starting from BCL files. If your job starts with FASTQ files, and only need to run ``cellranger count`` part, please refer to `this subsection <./index.html#run-cellranger-count-only>`_.

1. Import ``cellranger_workflow``
+++++++++++++++++++++++++++++++++

	Import *cellranger_workflow* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose workflow **github.com/lilab-bcb/cumulus/CellRanger** to import.

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cellranger_workflow* workflow in the drop-down menu.

2. Upload sequencing data to Google bucket
++++++++++++++++++++++++++++++++++++++++++

	Copy your sequencing output to your workspace bucket using gsutil_ (you already have it if you've installed Google cloud SDK) in your unix terminal.

	You can obtain your bucket URL in the dashboard tab of your Terra workspace under the information panel.

	.. image:: ../images/google_bucket_link.png

	Use ``gsutil cp [OPTION]... src_url dst_url`` to copy data to your workspace bucket. For example, the following command copies the directory at /foo/bar/nextseq/Data/VK18WBC6Z4 to a Google bucket::

		gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4

	``-m`` means copy in parallel, ``-r`` means copy the directory recursively, and ``gs://fc-e0000000-0000-0000-0000-000000000000`` should be replaced by your own workspace Google bucket URL.

.. note::
	If input is a folder of BCL files, users do not need to upload the whole folder to the Google bucket. Instead, they only need to upload the following files::

		RunInfo.xml
		RTAComplete.txt
		runParameters.xml
		Data/Intensities/s.locs
		Data/Intensities/BaseCalls

	If data are generated using MiSeq or NextSeq, the location files are inside lane subfloders ``L001`` under ``Data/Intensities/``. In addition, if users' data only come from a subset of lanes (e.g. ``L001`` and ``L002``), users only need to upload lane subfolders from the subset (e.g. ``Data/Intensities/BaseCalls/L001, Data/Intensities/BaseCalls/L002`` and ``Data/Intensities/L001, Data/Intensities/L002`` if sequencer is MiSeq or NextSeq).

Alternatively, users can submit jobs through command line interface (CLI) using `altocumulus <./command_line.html>`_, which will smartly upload BCL folders according to the above rules.

.. note:: Broad users need to be on an UGER node (not a login node) in order to use the ``-m`` flag

	Request an UGER node::

		reuse UGER
		qrsh -q interactive -l h_vmem=4g -pe smp 8 -binding linear:8 -P regevlab

	The above command requests an interactive node with 4G memory per thread and 8 threads. Feel free to change the memory, thread, and project parameters.

	Once you're connected to an UGER node, you can make gsutil_ available by running::

		reuse Google-Cloud-SDK

3. Prepare a sample sheet
+++++++++++++++++++++++++

	**3.1 Sample sheet format**:

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to demultiplex flowcells and generate channel-specific count matrices. Note that *Sample*, *Lane*, and *Index* columns are defined exactly the same as in 10x's simple CSV layout file.

	A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name. Sample name can only contain characters from [a-zA-Z0-9\_-].
		* - **Reference**
		  -
		  	| Provides the reference genome used by Cell Ranger for each 10x channel.
		  	| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as *GRCh38-2020-A*.
		  	| A full list of available keywords is included in each of the following data type sections (e.g. sc/snRNA-seq) below.
		* - **Flowcell**
		  -
		    | Indicates the Google bucket URLs of uploaded BCL folders.
		    | If starts with FASTQ files, this should be Google bucket URLs of uploaded FASTQ folders.
		    | The FASTQ folders should contain one subfolder for each sample in the flowcell with the sample name as the subfolder name.
		    | Each subfolder contains FASTQ files for that sample.
		* - **Lane**
		  -
		    | Tells which lanes the sample was pooled into.
		    | Can be either single lane (e.g. 8) or a range (e.g. 7-8) or all (e.g. \*).
		* - **Index**
		  - Sample index (e.g. SI-GA-A12).
		* - Chemistry
		  - Describes the 10x chemistry used for the sample. This column is optional.
		* - DataType
		  -
			| Describes the data type of the sample --- *rna*, *vdj*, *citeseq*, *hashing*, *cmo*, *crispr*, *atac*.
			| **rna** refers to gene expression data (*cellranger count*),
			| **vdj** refers to V(D)J data (*cellranger vdj*),
			| **citeseq** refers to CITE-Seq tag data,
			| **hashing** refers to cell-hashing or nucleus-hashing tag data,
			| **adt**, which refers to the case where *hashing* and *citeseq* reads are in a sample library.
			| **cmo** refers to cell multiplexing oligos used in 10x Genomics' CellPlex assay,
			| **crispr** refers to Perturb-seq guide tag data,
			| **atac** refers to scATAC-Seq data (*cellranger-atac count*),
			| This column is optional and the default data type is *rna*.
		* - FeatureBarcodeFile
		  -
		  	| Google bucket urls pointing to feature barcode files for *rna*, *citeseq*, *hashing*, *cmo* and *crispr* data.
		  	| Features can be either targeted genes for targeted gene expression analysis, antibody for CITE-Seq, cell-hashing, nucleus-hashing or gRNA for Perburb-seq.
		  	| If *cmo* data is analyzed separately using *cumulus_feature_barcoding*, file format should follow the guide in Feature barcoding assays section, otherwise follow the guide in Single-cell multiomics section.
		  	| This column is only required for targeted gene expression analysis (*rna*), CITE-Seq (*citeseq*), cell-hashing or nucleus-hashing (*hashing*), CellPlex (*cmo*) and Perturb-seq (*crispr*).
		* - Link
		  -
			| Designed for Single Cell Multiome	ATAC + Gene Expression, Feature Barcoding, or CellPlex.
			| Link multiple modalities together using a single link name.
			| cellranger-arc count, cellranger count, or cellranger multi will be triggered automatically depending on the modalities.
			| If empty string is provided, no link is assumed.
			| Link name can only contain characters from [a-zA-Z0-9\_-].


	The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 4 samples sequenced in two flowcells.

	Example::

		Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,threeprime,rna
		sample_2,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,SC3Pv3,rna
		sample_3,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,5-6,SI-GA-C8,fiveprime,rna
		sample_4,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,7-8,SI-GA-D8,fiveprime,rna
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,threeprime,rna
		sample_2,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,SC3Pv3,rna
		sample_3,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,5-6,SI-GA-C8,fiveprime,rna
		sample_4,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,7-8,SI-GA-D8,fiveprime,rna

	**3.2 Upload your sample sheet to the workspace bucket:**

		Example::

			gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/

4. Launch analysis
++++++++++++++++++

	In your workspace, open ``cellranger_workflow`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Run workflow with inputs defined by file paths`` as below

		.. image:: ../images/single_workflow.png

	and click ``SAVE`` button. Select ``Use call caching`` and click ``INPUTS``. Then fill in appropriate values in the ``Attribute`` column. Alternative, you can upload a JSON file to configure input by clicking ``Drag or click to upload json``.

	Once INPUTS are appropriated filled, click ``RUN ANALYSIS`` and then click ``LAUNCH``.

5. Notice: run ``cellranger mkfastq`` if you are non Broad Institute users
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	Non Broad Institute users that wish to run ``cellranger mkfastq`` must create a custom docker image that contains ``bcl2fastq``.

		See :ref:`bcl2fastq-docker` instructions.

6. Run ``cellranger count`` only
++++++++++++++++++++++++++++++++++++

	Sometimes, users might want to perform demultiplexing locally and only run the count part on the cloud. This section describes how to only run the count part via ``cellranger_workflow``.

	#. Copy your FASTQ files to the workspace using gsutil_ in your unix terminal. There are two cases:

		- **Case 1**: All the FASTQ files are in one top-level folder. Then you can simply upload this folder to Cloud, and in your sample sheet, make sure **Sample** names are consistent with the filename prefix of their corresponding FASTQ files.
		- **Case 2**: In the top-level folder, each sample has a dedicated subfolder containing its FASTQ files. In this case, you need to upload the whole top-level folder, and in your sample sheet, make sure **Sample** names and their corresponding subfolder names are identical.

		Notice that if your FASTQ files are downloaded from the Sequence Read Archive (SRA) from NCBI, you must rename your FASTQs to follow the bcl2fastq `file naming conventions`_.

		Example::

			gsutil -m cp -r /foo/bar/fastq_path/K18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq

	#. Create a sample sheet following the similar structure as `above <./index.html#prepare-a-sample-sheet>`_, except the following differences:

		- **Flowcell** column should list Google bucket URLs of the FASTQ folders for flowcells.
		- **Lane** and **Index** columns are NOT required in this case.

		Example::

			Sample,Reference,Flowcell
			sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq

	#. Set optional input ``run_mkfastq`` to ``false``.


7. Workflow outputs
+++++++++++++++++++

	See the table below for workflow level outputs.

	.. list-table::
		:widths: 5 5 10
		:header-rows: 1

		* - Name
		  - Type
		  - Description
		* - fastq_outputs
		  - Array[Array[String]?]
		  - The top-level array contains results (as arrays) for different data modalities. The inner-level array contains cloud locations of FASTQ files, one url per flowcell.
		* - count_outputs
		  - Array[Array[String]?]
		  - The top-level array contains results (as arrays) for different data modalities. The inner-level array contains cloud locations of count matrices, one url per sample.
		* - count_matrix
		  - String
		  - Cloud url for a template count_matrix.csv to run Cumulus. It only contains sc/snRNA-Seq samples.


.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _Import workflows to Terra: ../cumulus_import.html
.. _file naming conventions: https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-
