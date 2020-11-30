Run STARsolo to generate gene-count matrices from FASTQ files
----------------------------------------------------------------------

This ``star_solo`` workflow generates gene-count matrices from FASTQ data using STARsolo.

Prepare input data and import workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Run ``cellranger_workflow`` to generate FASTQ data
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	You can skip this step if your data are already in FASTQ format.

	Otherwise, for 10X data, you need to first run *cellranger_workflow* to generate FASTQ files from BCL raw data for each sample. Please follow `cellranger_workflow manual <cellranger.html>`_.

	Notice that you should set **run_mkfastq** to ``true`` to get FASTQ output. You can also set **run_count** to ``false`` to skip Cell Ranger count step.

	For Non-Broad users, you'll need to build your own docker for *bcl2fastq* step. Instructions are `here <bcl2fastq.html>`_.

2. Import ``star_solo``
+++++++++++++++++++++++

	Import *star_solo* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *star_solo* workflow is under ``Broad Methods Repository`` with name "**cumulus/star_solo**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *star_solo* workflow in the drop-down menu.

3. Prepare a sample sheet
++++++++++++++++++++++++++++

	**3.1 Sample sheet format:**

	The sample sheet for *star_solo* workflow should be in TSV format, i.e. columns are separated by tabs (NOT commas). Please note that the columns in the TSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to identify flowcells and generate sample/channel-specific count matrices.

	A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each sample or 10X channel should have a unique sample name.
		* - **Flowcells**
		  - Indicates the Google bucket URLs of folder(s) holding FASTQ files of this sample.

	For 10X data, the sample sheet supports sequencing the same 10X channel across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list all of its flowcells in a comma-seperated way. In the following example, we have 2 samples sequenced in two flowcells.

	Example::

		Sample	Flowcells
		sample_1	gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/sample_1_fastqs,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/sample_1_fastqs
		sample_2	gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/sample_2_fastqs

	**3.2 Upload your sample sheet to the workspace bucket:**

	Use gsutil_ (you already have it if you've installed Google cloud SDK) in your unix terminal to upload your sample sheet to workspace bucket.

	Example::

			gsutil cp /foo/bar/projects/sample_sheet.tsv gs://fc-e0000000-0000-0000-0000-000000000000/

4. Launch analysis
+++++++++++++++++++

	In your workspace, open ``star_solo`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Process single workflow from files`` as below

		.. image:: images/single_workflow.png

	and click ``SAVE`` button. Select ``Use call caching`` and click ``INPUTS``. Then fill in appropriate values in the ``Attribute`` column. Alternative, you can upload a JSON file to configure input by clicking ``Drag or click to upload json``.

	Once INPUTS are appropriated filled, click ``RUN ANALYSIS`` and then click ``LAUNCH``.

----------------------------

Workflow inputs
^^^^^^^^^^^^^^^^^^

Below are inputs for *count* workflow. Notice that required inputs are in bold.

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_tsv_file**
	  - Input TSV sample sheet describing metadata of each sample.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.tsv"
	  -
	* - **genome**
	  - Genome reference. It can be either of the following two formats:

		- String. Pre-built genome reference. Currently support: ``GRCh38``, ``mm10``.

		- Google bucket URL of a custom reference, must be a ``.tar.gz`` file.
	  - | "GRCh38",
	    | or "gs://user-bucket/starsolo.tar.gz"
	  -
	* - **chemistry**
	  - Chemistry name. Available options: "tenX_v3" (for 10X V3 chemistry), "tenX_v2" (for 10X V2 chemistry), "DropSeq", and "SeqWell".
	  - "tenX_v3"
	  -
	* - **output_directory**
	  - GS URL of output directory.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/count_result"
	  -
	* - docker_registry
	  - Docker registry to use:

	  	- "quay.io/cumulus" for docker images on Red Hat registry;

		- "cumulusprod" for backup images on Docker Hub.
	  -
	  - "quay.io/cumulus"
	* - zones
	  - Google cloud zones to consider for execution.
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - num_cpu
	  - Number of CPUs to request for count per sample.
	  - 32
	  - 32
	* - disk_space
	  - Disk space in GB needed for count per sample.
	  - 500
	  - 500
	* - memory
	  - Memory size in GB needed for count per sample.
	  - 120
	  - 120
	* - preemptible
	  - Number of maximum preemptible tries allowed.
	  - 2
	  - 2
	* - star_version
	  - STAR version to use. Currently only support ``2.7.6a``.
	  - "2.7.6a"
	  - "2.7.6a"
	* - config_version
	  - Version of docker image to run configuration on the sample sheet. Currently only has version "0.1".
	  - "0.1"
	  - "0.1"


Workflow outputs
^^^^^^^^^^^^^^^^^^^

See the table below for *star_solo* workflow outputs.

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_folder
	  - String
	  - Google Bucket URL of output directory. Within it, each folder is for one sample in the input sample sheet.

.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
