:orphan:

sccutil sample-sheet
----------------

*sccutil sample-sheet* creates a sample sheet.

Synopsis
^^^^^^^^

.. code-block:: shell

	sccutil sample-sheet [OPTION]... dir

where *dir* can be one or more directories to look for fastq files. Directories are searched recursively for fastq files.

Options
^^^^^^^

.. list-table::
	:widths: 5 20
	:header-rows: 1

	* - Name
	  - Description
	* - | -f
	    | --format
	  - Sample sheet format. Can be directory, r1_r2, r1_r2_i1. For cellranger, use directory, for drop-seq, use r1_r2
	* - | **-o**
	    | **--output**
	  - Output sample sheet

sccutil fc-upload
-------------

*sccutil fc-upload* uploads local files/directories to your Terra workspace Google Cloud bucket.
 
Your sample sheet can point to local file paths and sccutil upload will take care of uploading directories (e.g. fastq directories) and modifying the sample sheet to point to a Google Cloud bucket.

Synopsis
^^^^^^^^

.. code-block:: shell

	sccutil fc-upload [OPTION]... file

where file can be one or more files, such as a sample sheet or an input JSON.

Options
^^^^^^^

.. list-table::
	:widths: 5 20
	:header-rows: 1

	* - Name
	  - Description
	* - | **-w**
	    | **--workspace**
	  - Workspace name (e.g. foo/bar). The workspace will be created if it does not exist
	* - --dry-run
	  - Causes upload to run in "dry run" mode, i.e., just outputting what would be uploaded without actually doing any uploading.

Example: Upload fastq directories and sample sheet, convert sample sheet paths to gs:// URLs.

example_sample_sheet.txt:

.. list-table::
	:widths: 5 20
	:header-rows: 0

	* - s1
	  - /mylocalpath/flowcell1/s1
	* - s2
	  - /mylocalpath/flowcell1/s2
	* - s1
	  - /mylocalpath/flowcell2/s1


command::

	sccutil upload -w myworkspace_namespace/myworkspace_name example_sample_sheet.txt 


sccutil fc-inputs
-------------

*sccutil fc-inputs* generates a stub JSON input file that can be used as input to *sccutil run*. The JSON file can optionally be based on a published method config.

Options
^^^^^^^

.. list-table::
	:widths: 5 20
	:header-rows: 1

	* - Name
	  - Description

	* - | **-m**
	    | **--method**
	  - Method namespace/name (e.g. scCloud/cellranger_workflow). A version can optionally be specified (e.g. scCloud/cellranger_workflow/4), otherwise the latest version of the method is used.
	* - | -c
	    | --config
	  - Repository config to use for generating input.json stub (e.g. regev/drop-seq-MMUL_8_0_1
	* - | -o
	    | --out
	  - JSON output file

Example: Generate a stub JSON file based on the published config for running Drop-Seq using the MMUL_8_0_1 genome::

	sccutil fc-inputs -m regev/drop-seq -c regev/drop-seq-MMUL_8_0_1
