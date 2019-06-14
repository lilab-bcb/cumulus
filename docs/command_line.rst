Run Terra pipelines via command line
----------------------------------------------

You can run Terra pipelines via the command line by installing the **kco** tools.

Install ``kco`` tools for Broad users
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Request an UGER node::

	reuse UGER
	qrsh -q interactive -l h_vmem=4g -pe smp 8 -binding linear:8 -P regevlab

The above command requests an interactive shell using the regevlab project with 4G memory per thread, 8 threads. Feel free to change the memory, thread, and project parameters.

Add conda to your path::

	reuse Anaconda3

Activate the kco virtual environment::

	source activate /seq/regev_genome_portal/conda_env/kco_tools

Install ``kco`` tools for non-Broad users
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Make sure you have ``conda`` installed. If you haven't installed ``conda``, use the following commands to install::

	wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh .
	bash Miniconda3-latest-Linux-x86_64.sh -p /users/foo/miniconda3
	mv Miniconda3-latest-Linux-x86_64.sh /users/foo/miniconda3

Next, create an virtual environment ``kco`` and install ``kco`` tools::

	conda create -n kco -y pip
	source activate kco
	git clone https://github.com/klarman-cell-observatory/KCO.git kco_tools
	cd kco_tools
	pip install -e .

---------------------------------

Run Terra pipelines via ``kco fc_run``
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**kco fc_run** runs a Terra method. Features:

#. Uploads local files/directories in your inputs to a Google Cloud bucket updates the file paths to point to the Google Cloud bucket. Your sample sheet can point to local file paths and kco run will take care of uploading directories (e.g. fastq directories) and modifying the sample sheet to point to a Google Cloud bucket.

#. Creates or uses an existing workspace.

#. Uses the latest version of a method unless the method version is specified.

Options
+++++++

.. list-table::
	:widths: 5 20
	:header-rows: 1

	* - Name
	  - Description
	* - | **-m**
		| **--method**
	  - | Method namespace/name (e.g. scCloud/cellranger_workflow).
		| A version can optionally be specified (e.g. scCloud/cellranger_workflow/4), otherwise the latest version of the method is used
	* - | **-w**
		| **--workspace**
	  - Workspace name (e.g. foo/bar). The workspace will be created if it does not exist
	* - | **-i**
		| **--inputs**
	  - WDL input JSON. Can be specified as file or string
	* - --bucket-folder
	  - Store inputs to <folder> under workspace's google bucket
	* - | -o
		| --upload
	  - Upload files/directories to the workspace Google Cloud bucket and output updated input json (with local path replaced by google bucket urls) to <updated_json>

Examples
++++++++

Upload BCL directories and sample sheet, convert sample sheet paths to gs:// URLs, and run cell ranger mkfastq and count.

example_sample_sheet.csv::

	Sample,Reference,Flowcell,Lane,Index,Chemistry
	sample_1,GRCh38,/mylocalpath/flowcell1,1-2,SI-GA-A8,threeprime
	sample_2,GRCh38,/mylocalpath/flowcell1,3-4,SI-GA-B8,threeprime
	sample_3,mm10,/mylocalpath/flowcell1,5-6,SI-GA-C8,fiveprime
	sample_4,mm10,/mylocalpath/flowcell1,7-8,SI-GA-D8,fiveprime
	sample_1,GRCh38,/mylocalpath/flowcell2,1-2,SI-GA-A8,threeprime
	sample_2,GRCh38,/mylocalpath/flowcell2,3-4,SI-GA-B8,threeprime
	sample_3,mm10,/mylocalpath/flowcell2,5-6,SI-GA-C8,fiveprime
	sample_4,mm10,/mylocalpath/flowcell2,7-8,SI-GA-D8,fiveprime

Note that sample_1, sample_2, sample_3, and sample_4 are sequenced on 2 flowcells and for each sample, all of its FASTQ files will be passed to cell ranger count in one command by the pipeline.

inputs.json:

.. code-block:: json

	{
		"cellranger_workflow.input_csv_file" : "mylocalpath/sample_sheet.csv",
		"cellranger_workflow.output_directory" : "gs://url/outputs",
		"cellranger_workflow.delete_input_bcl_directory": true
	}

Run the following command to kick off your Terra pipeline::

	kco fc_run -m scCloud/cellranger_workflow -i inputs.json -w myworkspace_namespace/myworkspace_name --bucket-folder inputs -o inputs_updated.json

Upon success, **kco fc_run** returns a url pointing the the submitted Terra job.

If for any reason, your job failed. You could rerun it without uploading files again via the following command::

	kco fc_run -m scCloud/cellranger_workflow -i inputs_updated.json -w myworkspace_namespace/myworkspace_name
