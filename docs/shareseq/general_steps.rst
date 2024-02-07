This section mainly considers jobs starting from BCL files. If your job starts with FASTQ files, and only need to run ``count`` part, please refer to `this subsection <./index.html#run-count-only>`_.

1. Import ``shareseq_workflow``
+++++++++++++++++++++++++++++++++

	Import *shareseq_workflow* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose workflow **github.com/kalarman-cell-observatory/cumulus/SHARE-seq** to import.

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *shareseq_workflow* workflow in the drop-down menu.

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


3. Prepare a sample sheet
+++++++++++++++++++++++++

	**3.1 Sample sheet format**:

	Please note that the columns in the CSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to demultiplex flowcells and generate channel-specific count matrices.

	A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each channel should have a unique sample name. Sample name can only contain characters from [a-zA-Z0-9\_-].
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
                    | Required for mkfastq.
		* - **Index**
		  -
                    | Sample index (e.g. A1).
                    | Required for mkfastq.
		* - Reference
		  -
		  	| Provides the reference genome used by tools (Chromap or STARsolo) for each channel.
		  	| The elements in the *reference* column can be either Google bucket URLs to reference tarballs or keywords such as *GRCh38-2020-A*.
		  	| A full list of available keywords are included in `this sample sheet subsection <./index.html#sample-sheet>`_.
                        | Required for count step.
		* - Type
		  -
			| Describes the data type of the sample --- *gex*, *atac*.
			| **gex** refers to gene expression data (*STARsolo*),
			| **atac** refers to ATAC-seq data (*Chromap*),
                        | Required for count step.

	The sample sheet supports sequencing the same channels across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list it in multiple rows, with one flowcell per row. In the following example, we have 2 samples sequenced in two flowcells.

	Example::

		Sample,Reference,Flowcell,Lane,Index,Type
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,A1,gex
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,A2,gex
		sample_2,GRCh38_chromap_v0.1.3,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4,1-2,A1,atac
		sample_2,GRCh38_chromap_v0.1.3,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,3-4,A2,atac

	**3.2 Upload your sample sheet to the workspace bucket:**

		Example::

			gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/

4. Launch analysis
++++++++++++++++++

	In your workspace, open ``shareseq_workflow`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Run workflow with inputs defined by file paths`` as below

		.. image:: ../images/single_workflow.png

	and click ``SAVE`` button. Select ``Use call caching`` and click ``INPUTS``. Then fill in appropriate values in the ``Attribute`` column. Alternative, you can upload a JSON file to configure input by clicking ``Drag or click to upload json``.

	Once INPUTS are appropriated filled, click ``RUN ANALYSIS`` and then click ``LAUNCH``.

5. Notice: run ``shareseq mkfastq`` if you are non Broad Institute users
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	Non Broad Institute users that wish to run ``shareseq mkfastq`` must create a custom docker image that contains ``bcl2fastq``.

		See :ref:`bcl2fastq-docker` instructions.

6. Run ``count`` only
++++++++++++++++++++++++++++++++++++

Sometimes, users might want to perform demultiplexing locally and only run the count part on the cloud. This section describes how to only run the count part via ``shareseq_workflow``.

#. Create a sample sheet following the similar structure as `above <./index.html#prepare-a-sample-sheet>`_, except the following differences:

	- **Flowcell** column should list Google bucket URLs of the FASTQ folders for flowcells.
	- **Lane** and **Index** columns are NOT required in this case.

	Example::

		Sample,Reference,Flowcell,Type
		sample_1,GRCh38-2020-A,gs://fc-e0000000-0000-0000-0000-000000000000/K18WBC6Z4_fastq,gex
                sample_2,GRCh38_chromap_v0.1.3,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2,1-2,A1,atac

#. Set optional input ``run_mkfastq`` to ``false``.


.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _Import workflows to Terra: ../cumulus_import.html
