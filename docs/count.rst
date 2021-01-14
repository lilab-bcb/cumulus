Cell Ranger alternatives to generate gene-count matrices for 10X data
----------------------------------------------------------------------

This ``count`` workflow generates gene-count matrices from 10X FASTQ data using alternative methods other than Cell Ranger.

Prepare input data and import workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Run ``cellranger_workflow`` to generate FASTQ data
++++++++++++++++++++++++++++++++++++++++++++++++++++++++

	You can skip this step if your data are already in FASTQ format.

	Otherwise, you need to first run *cellranger_workflow* to generate FASTQ files from BCL raw data for each sample. Please follow `cellranger_workflow manual <cellranger.html>`_.

	Notice that you should set **run_mkfastq** to ``true`` to get FASTQ output. You can also set **run_count** to ``false`` if you want to skip Cell Ranger count, and only use the result from *count* workflow.

	For Non-Broad users, you'll need to build your own docker for *bcl2fastq* step. Instructions are `here <bcl2fastq.html>`_.

2. Import ``count``
+++++++++++++++++++++++

	Import *count* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The *count* workflow is under ``Broad Methods Repository`` with name "**cumulus/count**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *count* workflow in the drop-down menu.

3. Prepare a sample sheet
++++++++++++++++++++++++++++

	**3.1 Sample sheet format:**

	The sample sheet for *count* workflow should be in TSV format, i.e. columns are seperated by tabs not commas. Please note that the columns in the TSV can be in any order, but that the column names must match the recognized headings.

	The sample sheet describes how to identify flowcells and generate channel-specific count matrices.

	A brief description of the sample sheet format is listed below **(required column headers are shown in bold)**.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - **Sample**
		  - Contains sample names. Each 10x channel should have a unique sample name.
		* - **Flowcells**
		  - Indicates the Google bucket URLs of folder(s) holding FASTQ files of this sample.

	The sample sheet supports sequencing the same 10x channel across multiple flowcells. If a sample is sequenced across multiple flowcells, simply list all of its flowcells in a comma-seperated way. In the following example, we have 2 samples sequenced in two flowcells.

	Example::

		Sample	Flowcells
		sample_1	gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/sample_1_fastqs,gs://fc-e0000000-0000-0000-0000-000000000000/VK10WBC9Z2/sample_1_fastqs
		sample_2	gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/sample_2_fastqs

	Moreover, if one flowcell of a sample contains multiple FASTQ files for each read, i.e. sequences from multiple lanes, you should keep your sample sheet as the same, and *count* workflow will automatically merge lanes altogether for the sample before performing counting.

	**3.2 Upload your sample sheet to the workspace bucket:**

	Use gsutil_ (you already have it if you've installed Google cloud SDK) in your unix terminal to upload your sample sheet to workspace bucket.

	Example::

			gsutil cp /foo/bar/projects/sample_sheet.tsv gs://fc-e0000000-0000-0000-0000-000000000000/

4. Launch analysis
+++++++++++++++++++

	In your workspace, open ``count`` in ``WORKFLOWS`` tab. Select the desired snapshot version (e.g. latest). Select ``Process single workflow from files`` as below

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
	  - Genome reference name. Current support: GRCh38, mm10.
	  - "GRCh38"
	  -
	* - **chemistry**
	  - 10X genomics' chemistry name. Current support: "tenX_v3" (for V3 chemistry), "tenX_v2" (for V2 chemistry), "dropseq" (for Drop-Seq).
	  - "tenX_v3"
	  -
	* - **output_directory**
	  - GS URL of output directory.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/count_result"
	  -
	* - run_count
	  - If you want to run count tools to generate gene-count matrices.
	  - true
	  - true
	* - count_tool
	  - Count tool to generate result. Options:

	  	- "StarSolo": Use STARsolo_.

	  	- "Optimus": Use Optimus_ pipeline, developed by the Data Coordination Platform team of the Human Cell Atlas.

	  	- "Bustools": Use `Kallisto BUSTools`_.

	  	- "Alevin": Use `Salmon Alevin`_.
	  - "StarSolo"
	  - "StarSolo"
	* - docker_registry
	  - Docker registry to use. Notice that docker image for Bustools is seperate.

	  	- "quay.io/cumulus" for images on Red Hat registry;

	  	- "cumulusprod" for backup images on Docker Hub.
	  - "quay.io/cumulus"
	  - "quay.io/cumulus"
	* - config_version
	  - Version of config docker image to use. This docker is used for parsing the input sample sheet for downstream execution. Currently only one version is available: ``0.1``.
	  - "0.1"
	  - "0.1"
	* - zones
	  - Google cloud zones to consider for execution.
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - num_cpu
	  - | Number of CPUs to request for count per channel.
	    | Notice that when use Optimus for count, this input only affects steps of copying files. Optimus uses CPUs due to its own strategy.
	  - 32
	  - 32
	* - disk_space
	  - | Disk space in GB needed for count per channel.
	    | Notice that when use Optimus for count, this input only affects steps of copying files. Optimus uses disk space due to its own strategy.
	  - 500
	  - 500
	* - memory
	  - | Memory size in GB needed for count per channel.
	    | Notice that when use Optimus for count, this input only affects steps of copying files. Optimus uses memory size due to its own strategy.
	  - 120
	  - 120
	* - preemptible
	  - | Number of maximum preemptible tries allowed.
	    | Notice that when use Optimus for count, this input only affects steps of copying files. Optimus uses preemptible tries due to its own strategy.
	  - 2
	  - 2
	* - merge_fastq_memory
	  - Memory size in GB needed for merge fastq per channel.
	  - 32
	  - 32
	* - starsolo_star_version
	  - | STAR version to use. Currently only support "2.7.3a".
	    | This input only works when setting *count_tool* to ``StarSolo``.
	  - "2.7.3a"
	  - "2.7.3a"
	* - alevin_version
	  - | Salmon version to use. Currently only support "1.1".
	    | This input only works when setting *count_tool* to ``Alevin``.
	  - "1.1"
	  - "1.1"
	* - bustools_output_loom
	  - | If BUSTools generates gene-count matrices in ``loom`` format.
	    | This input only works when setting *count_tool* to ``Bustools``.
	  - false
	  - false
	* - bustools_output_h5ad
	  - | If BUSTools generates gene-count matrices in ``h5ad`` format.
	    | This input only works when setting *count_tool* to ``Bustools``.
	  - false
	  - false
	* - bustools_docker
	  - | Docker image used for Kallisto BUSTools count.
	    | This input only works when setting *count_tool* to ``Bustools``.
	  - "shaleklab/kallisto-bustools"
	  - "shaleklab/kallisto-bustools"
	* - bustools_version
	  - | kb version to use. Currently only support "0.24.4".
	    | This input only works when setting *count_tool* to ``Bustools``.
	  - "0.24.4"
	  - "0.24.4"
	* - optimus_output_loom
	  - | If Optimus generates gene-count matrices in ``loom`` format.
	    | This input only works when setting *count_tool* to ``Optimus``.
	  - true
	  - true


Workflow outputs
^^^^^^^^^^^^^^^^^^^

See the table below for *count* workflow outputs.

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
.. _STARsolo: https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md
.. _Optimus: https://github.com/HumanCellAtlas/skylab/tree/master/pipelines/optimus
.. _Kallisto BUSTools: https://www.kallistobus.tools/introduction
.. _Salmon Alevin: https://salmon.readthedocs.io/en/latest/alevin.html
