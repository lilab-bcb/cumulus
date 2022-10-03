10x Visium
----------

Run Space Ranger tools using spaceranger_workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

``spaceranger_workflow`` wraps Space Ranger to process 10x Visium data.

A general step-by-step instruction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section mainly considers jobs starting from BCL files. If your job starts with FASTQ files, and only need to run ``spaceranger count`` part, please refer to `this subsection <./spaceranger.html#run-spaceranger-count-only>`_.

1. Import ``spaceranger_workflow``
++++++++++++++++++++++++++++++++++

	Import *spaceranger_workflow* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose workflow **github.com/lilab-bcb/cumulus/Spaceranger** to import.

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *spaceranger_workflow* workflow in the drop-down menu.

2. Upload sequencing and image data to Google bucket
++++++++++++++++++++++++++++++++++++++++++++++++++++

	Copy your sequencing output to your workspace bucket using gsutil_ (you already have it if you've installed Google cloud SDK) in your unix terminal.

	You can obtain your bucket URL in the dashboard tab of your Terra workspace under the information panel.

	.. image:: images/google_bucket_link.png

	Use ``gsutil cp [OPTION]... src_url dst_url`` to copy data to your workspace bucket. For example, the following command copies the directory at /foo/bar/nextseq/Data/VK18WBC6Z4 to a Google bucket::

		gsutil -m cp -r /foo/bar/nextseq/Data/VK18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4

	``-m`` means copy in parallel, ``-r`` means copy the directory recursively, and ``gs://fc-e0000000-0000-0000-0000-000000000000`` should be replaced by your own workspace Google bucket URL.

	Similarly, copy all images for spatial data to the same google bucket.

.. note::
	If input is a folder of BCL files, users do not need to upload the whole folder to the Google bucket. Instead, they only need to upload the following files::

		RunInfo.xml
		RTAComplete.txt
		runParameters.xml
		Data/Intensities/s.locs
		Data/Intensities/BaseCalls

	If data are generated using MiSeq or NextSeq, the location files are inside lane subfloders ``L001`` under ``Data/Intensities/``. In addition, if users' data only come from a subset of lanes (e.g. ``L001`` and ``L002``), users only need to upload lane subfolders from the subset (e.g. ``Data/Intensities/BaseCalls/L001, Data/Intensities/BaseCalls/L002`` and ``Data/Intensities/L001, Data/Intensities/L002`` if sequencer is MiSeq or NextSeq).

Alternatively, users can submit jobs through command line interface (CLI) using `altocumulus <./command_line.html>`_, which will smartly upload BCL folders according to the above rules.


3. Prepare a sample sheet
+++++++++++++++++++++++++

	**3.1 Sample sheet format**:

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	For **FFPE** data, **ProbeSet** column is mandatory.

	The sample sheet describes how to demultiplex flowcells and generate channel-specific count matrices. Note that *Sample*, *Lane*, and *Index* columns are defined exactly the same as in 10x's simple CSV layout file.

	A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name.
		* - **Reference**
		  -
		  	| Provides the reference genome used by Space Ranger for each 10x channel.
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
		* - ProbeSet
		  - Probe set for FFPE samples. **Choosing** from ``human_probe_v1`` (10x human probe set, CytoAssist-incompatible), ``human_probe_v2`` (10x human probe set, CytoAssist-compatible) and ``mouse_probe_v1`` (10x mouse probe set). Alternatively, a CSV file describing the probe set can be directly used. Setting ProbeSet to ```` for a sample implies the sample is not FFPE.
		* - Image
		  - Cloud bucket url for a brightfield tissue H&E image in .jpg or .tiff format. This column is mutually exclusive with DarkImage and ColorizedImage columns.
		* - DarkImage
		  - Cloud bucket urls for Multi-channel, dark-background fluorescence image as either a single, multi-layer .tiff file, multiple .tiff or .jpg files, or a pre-combined color .tiff or .jpg file. If multiple files are provided, please separate them by ';'. This column is mutually exclusive with Image and ColorizedImage columns.
		* - ColorizedImage
		  - Cloud bucket url for a color composite of one or more fluorescence image channels saved as a single-page, single-file color .tiff or .jpg. This column is mutually exclusive with Image and DarkImage columns.
		* - CytaImage
		  - Cloud bucket url for a brightfield image generated by the CytAssist instrument.
		* - Slide
		  - Visium slide serial number. If both Slide and Area are empty, the --unknown-slide option would be set.
		* - Area
		  - Visium capture area identifier. Options for Visium are A1, B1, C1, D1. If both Slide and Area are empty, the --unknown-slide option would be set.
		* - SlideFile
		  - Slide layout file indicating capture spot and fiducial spot positions. Only required if internet access is not available.
		* - LoupeAlignment
		  - Alignment file produced by the manual Loupe alignment step.
		* - TargetPanel
		  - Cloud bucket url for a target panel CSV for targeted gene expression analysis.

	The sample sheet supports sequencing the same 10x channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 2 samples sequenced in two flowcells.

	Example::

		Sample,Reference,Flowcell,Lane,Index,Image,Slide,Area
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,SI-GA-A8,gs://image/image1.tif,V19J25-123,A1
		sample_2,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,3-4,SI-GA-B8,gs://image/image2.tif,V19J25-123,B1
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,SI-GA-A8,gs://image/image1.tif,V19J25-123,A1
		sample_2,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,SI-GA-B8,gs://image/image2.tif,V19J25-123,B1

	**3.2 Upload your sample sheet to the workspace bucket:**

		Example::

			gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/

4. Launch analysis
++++++++++++++++++

	In your workspace, open ``spaceranger_workflow`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Run workflow with inputs defined by file paths`` as below

		.. image:: images/single_workflow.png

	and click ``SAVE`` button. Select ``Use call caching`` and click ``INPUTS``. Then fill in appropriate values in the ``Attribute`` column. Alternative, you can upload a JSON file to configure input by clicking ``Drag or click to upload json``.

	Once INPUTS are appropriated filled, click ``RUN ANALYSIS`` and then click ``LAUNCH``.

5. Notice: run ``spaceranger mkfastq`` if you are non Broad Institute users
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	Non Broad Institute users that wish to run ``spaceranger mkfastq`` must create a custom docker image that contains ``bcl2fastq``.

		See :ref:`bcl2fastq-docker` instructions.

6. Run ``spaceranger count`` only
++++++++++++++++++++++++++++++++++++

Sometimes, users might want to perform demultiplexing locally and only run the count part on the cloud. This section describes how to only run the count part via ``spaceranger_workflow``.

#. Copy your FASTQ files to the workspace using gsutil_ in your unix terminal. There are two cases:

	- **Case 1**: All the FASTQ files are in one top-level folder. Then you can simply upload this folder to Cloud, and in your sample sheet, make sure **Sample** names are consistent with the filename prefix of their corresponding FASTQ files.
	- **Case 2**: In the top-level folder, each sample has a dedicated subfolder containing its FASTQ files. In this case, you need to upload the whole top-level folder, and in your sample sheet, make sure **Sample** names and their corresponding subfolder names are identical.

	Notice that if your FASTQ files are downloaded from the Sequence Read Archive (SRA) from NCBI, you must rename your FASTQs to follow the bcl2fastq `file naming conventions`_.

	Example::

		gsutil -m cp -r /foo/bar/fastq_path/K18WBC6Z4 gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq

#. Create a sample sheet following the similar structure as `above <./spaceranger.html#prepare-a-sample-sheet>`_, except the following differences:

	- **Flowcell** column should list Google bucket URLs of the FASTQ folders for flowcells.
	- **Lane** and **Index** columns are NOT required in this case.

	Example::

		Sample,Reference,Flowcell,Image,Slide,Area
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq,gs://image/image1.tif,V19J25-123,A1

#. Set optional input ``run_mkfastq`` to ``false``.

---------------------------------

Visium spatial transcriptomics data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To process spatial transcriptomics data, follow the specific instructions below.

Sample sheet
++++++++++++

#. **Reference** column.

	Pre-built scRNA-seq references are summarized below.

	.. list-table::
		:widths: 5 20
		:header-rows: 1

		* - Keyword
		  - Description
		* - **GRCh38-2020-A**
		  - Human GRCh38 (GENCODE v32/Ensembl 98)
		* - **mm10-2020-A**
		  - Mouse mm10 (GENCODE vM23/Ensembl 98)

Workflow input
++++++++++++++

For spatial data, ``spaceranger_workflow`` takes Illumina outputs and related images as input and runs ``spaceranger mkfastq`` and ``spaceranger count``. Revalant workflow inputs are described below, with required inputs highlighted in bold.

	.. list-table::
		:widths: 5 30 30 20
		:header-rows: 1

		* - Name
		  - Description
		  - Example
		  - Default
		* - **input_csv_file**
		  - Sample Sheet (contains Sample, Reference, Flowcell, Lane, Index as required and ProbeSet, Image, DarkImage, ColorizedImage, CytaImage, Slide, Area, SlideFile, LoupeAlignment, TargetPanel as optional)
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
		  -
		* - **output_directory**
		  - Output directory
		  - "gs://fc-e0000000-0000-0000-0000-000000000000/spaceranger_output"
		  - Results are written under directory *output_directory* and will overwrite any existing files at this location.
		* - run_mkfastq
		  - If you want to run ``spaceranger mkfastq``
		  - true
		  - true
		* - run_count
		  - If you want to run ``spaceranger count``
		  - true
		  - true
		* - delete_input_bcl_directory
		  - If delete BCL directories after demux. If false, you should delete this folder yourself so as to not incur storage charges
		  - false
		  - false
		* - mkfastq_barcode_mismatches
		  - Number of mismatches allowed in matching barcode indices (bcl2fastq2 default is 1)
		  - 0
		  -
		* - reorient_images
		  - For use with automatic fiducial alignment. This option will apply to all samples in the sample sheet. Spaceranger will attempt to find the best alignment of the fiducial markers by any rotation or mirroring of the image.
		  - true
		  - true
		* - filter_probes
		  - Whether to filter the probe set using the "included" column of the probe set CSV.
		  - true
		  - true
		* - dapi_index
		  - Index of DAPI channel (1-indexed) of fluorescence image, only used in the CytaAssist case, with dark background image.
		  - 2
		  -
		* - unknown_slide
		  - Use this option if the slide serial number and area identifier have been lost. Choose from visium-1, visium-2 and visium-2-large.
		  - visium-2
		  -
		* - no_bam
		  - Turn this option on to disable BAM file generation.
		  - false
		  - false
		* - secondary
		  - Perform Space Ranger secondary analysis (dimensionality reduction, clustering, etc.)
		  - false
		  - false
		* - r1_length
		  - Hard trim the input Read 1 to this length before analysis
		  - 28
		  -
		* - r2_length
		  - Hard trim the input Read 1 to this length before analysis. This value will be set to 50 automatically for FFPE samples if spaceranger version < 2.0.0.
		  - 50
		  -
		* - spaceranger_version
		  - spaceranger version, could be: ``2.0.0``, ``1.3.1``, ``1.3.0``
		  - "2.0.0"
		  - "2.0.0"
		* - config_version
		  - config docker version used for processing sample sheets, could be 0.3.
		  - "0.3"
		  - "0.3"
		* - docker_registry
		  - Docker registry to use for spaceranger_workflow. Options:

		  	- "quay.io/cumulus" for images on Red Hat registry;

		  	- "cumulusprod" for backup images on Docker Hub.
		  - "quay.io/cumulus"
		  - "quay.io/cumulus"
		* - spaceranger_mkfastq_docker_registry
		  - Docker registry to use for ``spaceranger mkfastq``.
		    Default is the registry to which only Broad users have access.
		    See :ref:`bcl2fastq-docker` for making your own registry.
		  - "gcr.io/broad-cumulus"
		  - "gcr.io/broad-cumulus"
		* - zones
		  - Google cloud zones
		  - "us-central1-a us-west1-a"
		  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
		* - num_cpu
		  - Number of cpus to request for one node for spaceranger mkfastq and spaceranger count
		  - 32
		  - 32
		* - memory
		  - Memory size string for spaceranger mkfastq and spaceranger count
		  - "120G"
		  - "120G"
		* - mkfastq_disk_space
		  - Optional disk space in GB for mkfastq
		  - 1500
		  - 1500
		* - count_disk_space
		  - Disk space in GB needed for spaceranger count
		  - 500
		  - 500
		* - backend
		  - Cloud infrastructure backend to use. Available options:

		    - "gcp" for Google Cloud;
		    - "aws" for Amazon AWS;
		    - "local" for local machine.
		  - "gcp"
		  - "gcp"
		* - preemptible
		  - Number of preemptible tries. This works only when *backend* is ``gcp``.
		  - 2
		  - 2
		* - awsQueueArn
		  - The AWS ARN string of the job queue to be used. This only works for ``aws`` backend.
		  - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
		  - ""

Workflow output
+++++++++++++++

See the table below for important sc/snRNA-seq outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - fastq_outputs
	  - Array[String]?
	  - A list of cloud urls containing FASTQ files, one url per flowcell.
	* - count_outputs
	  - Array[String]?
	  - A list of cloud urls containing spaceranger count outputs, one url per sample.
	* - metrics_summaries
	  - File?
	  - A excel spreadsheet containing QCs for each sample.
	* - spaceranger_count.output_web_summary
	  - Array[File]?
	  - A list of htmls visualizing QCs for each sample (spaceranger count output).

---------------------------------

Build Space Ranger References
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Reference built by Cell Ranger for sc/snRNA-seq should be compatible with Space Ranger. For more details on building references uing Cell Ranger, please refer to `here <./cellranger/index.html#build-references-for-sc-snrna-seq>`_.


.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _Import workflows to Terra: ./cumulus_import.html
.. _`file naming conventions`: https://kb.10xgenomics.com/hc/en-us/articles/115003802691-How-do-I-prepare-Sequence-Read-Archive-SRA-data-from-NCBI-for-Cell-Ranger-
