The workflow starts with FASTQ files.

.. note::
	Starting from v3.0.0, Cumulus cellranger_workflow drops support for ``mkfastq``. If your data start from BCL files, please first run `BCL Convert`_ to demultiplex flowcells to generate FASTQ files.

1. Import ``cellranger_workflow``
+++++++++++++++++++++++++++++++++

	Import *cellranger_workflow* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose workflow **github.com/lilab-bcb/cumulus/CellRanger** to import.

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cellranger_workflow* workflow in the drop-down menu.

2. Upload sequencing data to Google bucket
++++++++++++++++++++++++++++++++++++++++++

	Copy your FASTQ files to your workspace bucket using `gcloud storage`_ command (you already have it if you've installed Google cloud SDK) in your unix terminal.

	You can obtain your bucket URL in the dashboard tab of your Terra workspace under the information panel.

	.. image:: ../images/google_bucket_link.png

	There are three cases:

		- **Case 1**: All the FASTQ files are in one top-level folder. Then you can simply upload this folder to Cloud, and in your sample sheet, make sure **Sample** names are consistent with the filename prefix of their corresponding FASTQ files.
		- **Case 2**: In the top-level folder, each sample has a dedicated subfolder containing its FASTQ files. In this case, you need to upload the whole top-level folder, and in your sample sheet, make sure **Sample** names and their corresponding subfolder names are identical.
		- **Case 3**: Each sample's FASTQ files are wrapped in a TAR file. In this case, upload the folder which contains this TAR file. Also, make sure **Sample** names are consistent with the filename prefix of their corresponding FASTQ files inside the TAR files.

		Notice that if your FASTQ files are downloaded from the Sequence Read Archive (SRA) from NCBI, you must rename your FASTQs to follow the Illumina `file naming conventions`_.

		Example::

			gcloud storage cp -r /foo/bar/K18WBC6Z4/Fastq gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq

		where ``-r`` means copy the directory recursively, and ``fc-e0000000-0000-0000-0000-000000000000`` should be replaced by your own workspace Google bucket name.


Alternatively, users can submit jobs through command line interface (CLI) using `altocumulus <../command_line.html>`_, which will smartly upload FASTQ files to cloud.


3. Prepare a sample sheet
+++++++++++++++++++++++++

	**3.1 Sample sheet format**:

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to generate count matrices from sequencing reads. A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  -
		  	| Sample name. This name must be consistent with its corresponding FASTQ filename prefix in the folder specified in **Flowcell** column. Sample names can only contain characters from ``[a-zA-Z0-9\_-]`` to be recognized by Cell Ranger.
		  	| Notice that if a sample has multiple sequencing runs, each of which has FASTQ files stored in dedicated location, you can specify multiple entries in the sample sheet with the same name in **Sample** column, and each entry accounts for one FASTQ folder location.
		* - **Reference**
		  -
		  	| Provides the reference genome used by Cell Ranger for processing the sample.
		  	| The reference can be a keyword of prebuilt references (e.g. ``GRCh38-2020-A``) that stored in Cumulus bucket, or a user specified cloud URI to a custom reference (in tarball ``.tar.gz`` format).
		  	| A full list of available keywords is included in each of the following data type sections (e.g. sc/snRNA-seq) below.
		* - **Flowcell**
		  - Indicates the cloud URI of the uploaded folder containing FASTQ files for each sample.
		* - Chemistry
		  - Keywords to describe the 10x chemistry used for the sample. This column is optional. Check data type sections (e.g. sc/snRNA-seq) below for the corresponding list of available keywords.
		* - DataType
		  - Describes the data type of each sample, with keywords chosen from the list below. This column is optional, and the default is **rna**.

		  	- **rna**: Gene expression (GEX) data

		  	- **vdj**: V(D)J data

			- **citeseq**: CITE-Seq tag data

			- **hashing**: Cell-hashing or nucleus-hashing tag data

			- **adt**: For the case where *hashing* and *citeseq* reads are in the same sample library

			- **cmo**: Cell multiplexing oligos used in 10x Genomics' CellPlex assay

			- **crispr**: Perturb-seq guide tag data

			- **atac**: scATAC-Seq data

			- **frp**: 10x Flex gene expression (old name is Fixed RNA Profiling) data
		* - AuxFile
		  - The Cloud URI pointing to auxiliary files of the corresponding samples, with different usage depending on *DataType* values:

		  	- For *rna*: It's used by Sample Multiplexing methods, which specifies the sample name to multiplexing barcode mapping.

			- For *frp*: It's used by Flex data, which specifies the sample name to Flex probe barcode mapping.

			- For *citeseq*, *hashing*, *adt*, and *crispr*: It's the feature barcode file, which contains the information of antibody for CITE-Seq, cell-hashing, nucleus-hashing, or gNRA for Perturb-Seq.

				- If analyzing using *cumulus_feature_barcoding*, the feature barcode file should be in format specified in `Feature barcoding assays`_ section below;

				- If analyzing as part of the Sample Multiplexing data using ``cellranger multi``, the feature barcode file should be in `10x Feature Reference`_ format.

			- For *cmo*: It's the CMO reference file (``cmo-set`` option) when using custom CMOs in CellPlex data.

			- For *vdj_t_gd*: It's the inner enrichment primer file (``inner-enrichment-primers`` option) for VDJ-T-GD data.
		* - Link
		  -
			| Designed for Single Cell Multiome	ATAC + Gene Expression, Feature Barcoding, Sample Multiplexing, or Flex.
			| Link multiple modalities together using a single link name.
			| ``cellranger-arc count``, ``cellranger count``, or ``cellranger multi`` will be triggered automatically depending on the modalities.
			| If empty string is provided, no link is assumed.
			| Link name can only contain characters from ``[a-zA-Z0-9\_-]`` for Cell Ranger to recognize.
			| **Notice:** The Link names must be unique to *Sample* values to avoid overwriting each other's settings.



	The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 4 samples sequenced in two flowcells.

	Example::

		Sample,Reference,Flowcell,Chemistry,DataType
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,threeprime,rna
		sample_2,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,SC3Pv3,rna
		sample_3,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,fiveprime,rna
		sample_4,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/Fastq,fiveprime,rna
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/Fastq,threeprime,rna
		sample_2,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/Fastq,SC3Pv3,rna
		sample_3,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/Fastq,fiveprime,rna
		sample_4,mm10-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/Fastq,fiveprime,rna

	**3.2 Upload your sample sheet to the workspace bucket:**

		Example::

			gcloud storage cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/

Alternatively, users can submit jobs through command line interface (CLI) using `altocumulus <../command_line.html>`_, which will smartly upload FASTQ files to cloud.

4. Launch analysis
++++++++++++++++++

	In your workspace, open ``cellranger_workflow`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Run workflow with inputs defined by file paths`` as below

		.. image:: ../images/single_workflow.png

	and click ``SAVE`` button. Select ``Use call caching`` and click ``INPUTS``. Then fill in appropriate values in the ``Attribute`` column. Alternative, you can upload a JSON file to configure input by clicking ``Drag or click to upload json``.

	Once INPUTS are appropriated filled, click ``RUN ANALYSIS`` and then click ``LAUNCH``.



5. Workflow outputs
+++++++++++++++++++

	See the table below for workflow level outputs.

	.. list-table::
		:widths: 5 5 10
		:header-rows: 1

		* - Name
		  - Type
		  - Description
		* - count_outputs
		  - Map[String, Array[String]?]
		  - A modality-to-output map showing output URIs for all samples, organized by modality and one URI per sample.
		* - count_matrix
		  - String
		  - Cloud URI for a template count_matrix.csv to run Cumulus. It only contains sc/snRNA-Seq samples (i.e. with ``rna`` value in **DataType** column).


.. _BCL Convert: https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html
.. _gcloud storage: https://cloud.google.com/sdk/gcloud/reference/storage#COMMAND
.. _Import workflows to Terra: ../cumulus_import.html
.. _file naming conventions: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-specifying-fastqs#file-naming-convention
.. _Feature barcoding assays: ./index.html#feature-barcoding-assays-cell-nucleus-hashing-cite-seq-and-perturb-
.. _10x Feature Reference: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/inputs/cr-feature-ref-csv
