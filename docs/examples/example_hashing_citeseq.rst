Example of Cell-Hashing and CITE-Seq Analysis on Cloud
++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this example, you'll learn how to perform Cell-Hashing and CITE-Seq analysis using **cumulus** on Terra_.

-----------------------------

0. Workspace and Data Preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After registering on Terra and creating a workspace there, you'll need the following two information:

	* **Terra workspace name**. This is shown on your Terra workspace webpage, with format *"<workspace-namespace>/<workspace-name>"*. Let it be ``ws-lab/ws-01`` in this example, which means that your workspace has namespace **ws-lab** and name **ws-01**.
	* The corresponding **Google Cloud Bucket location** of your workspace. You can check it by clicking the link under **"Google Bucket"** title on your Terra workspace webpage. Let it be ``gs://fc-e0000000-0000-0000-0000-000000000000`` in this example.

Then upload your BCL directories to Google bucket of your workspace using gsutil_::

	gsutil -m cp -r /my-local-path/BCL/* gs://fc-e0000000-0000-0000-0000-000000000000/data-source

where ``/my-local-path/BCL`` is the path to the top-level directory of your BCL files on your local machine, and ``data-source`` is the folder on Google bucket to hold the uploaded data.

------------------------

1. Extract Gene-Count Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First step is to extract gene-count matrices from sequencing output. 


You need two original files from your dataset to start:

	* Cell-Hashing Index CSV file, say its filename is ``cell_hashing_index.csv``, of format *feature_barcode,feature_name*. See an example below::

		AATCATCACAAGAAA,CB1
		GGTCACTGTTACGTA,CB2
		... ...

	  where each line is a pair of feature barcode and feature name of a sample.

	* CITE-Seq Index CSV file, say its filename is ``cite_seq_index.csv``, of the same format as above. See an example below::

		TTACATGCATTACGA,CD19
		GCATTAGCATGCAGC,HLA-ABC
		... ...

	  where each line is a pair of Barcode and Specificity of an Antibody.

Then upload them to your Google Bucket using gsutil_. Assuming both files are in folder ``/Users/foo/data-source`` on your local machine, type the following command to upload::

	gsutil -m cp -r /Users/foo/data-source gs://fc-e0000000-0000-0000-0000-000000000000/data-source

where option ``-m`` means copy in parallel, ``-r`` means copy the directory recursively, and ``gs://fc-e0000000-0000-0000-0000-000000000000/data-source`` is your working directory at cloud side, which can be changed at your will.

Next, create a sample sheet, ``cellranger_sample_sheet.csv``, on your local machine with content below::

	Sample,Reference,Flowcell,Lane,Index,DataType,FeatureBarcodeFile
	sample_control,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,2,SI-GA-F1,rna
	sample_cc,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,3,SI-GA-A1,rna
	sample_cell_hashing,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,3,ATTACTCG,adt,cell_hashing_index.csv
	sample_cite_seq,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,3,CGTGAT,adt,cite_seq_index.csv

For the details on how to prepare this sample sheet, please refer to Step 3 of `Cell Ranger sample sheet instruction`_.

	.. _Cell Ranger sample sheet instruction: ../cellranger.html

When you are done with the sample sheet, upload it to Google bucket::

	gsutil cp cellranger_sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/

Now we are ready to set up **cellranger_workflow** workflow for this phase. If your workspace doesn't have this workflow, import it to your workspace by following Step 5 and 6 of `cellranger_workflow documentation <../cellranger.html>`_. 

Then prepare a JSON file, ``cellranger_inputs.json``, which is used to set up the workflow inputs::

	{
		"cellranger_workflow.input_csv_file" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/cellranger_sample_sheet.csv",
		"cellranger_workflow.output_directory" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir"
	}

where ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir`` is the remote directory in which the output of cellranger_workflow will be generated. For the details on the options above, please refer to `Cell Ranger workflow inputs`_.

	.. _Cell Ranger workflow inputs: ../cellranger.html#cellranger-workflow-inputs

When you are done with the JSON file, on cellranger_workflow workflow page, upload ``cellranger_inputs.json`` by clicking ``upload json`` link as below:

	.. image:: ../images/upload_json.png 
	   :scale: 70%

Then Click ``SAVE`` button to save the inputs, and click ``RUN ANALYSIS`` button as below to start the job:

	.. image:: ../images/run_analysis.png
	   :scale: 70%

When the execution is done, all the output results will be in folder ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir``. 

You'll need 4 files for the next phases. 3 are from the output:

	* RNA count matrix of the sample group of interest: ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cc/raw_feature_bc_matrix.h5``;
	* Cell-Hashing Antibody count matrix: ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cell_hashing/sample_cell_hashing.csv``;
	* CITE-Seq Antibody count matrix: ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cite_seq/sample_cite_seq.csv``.

Besides, create a sample sheet, ``citeseq_antibody_control.csv``, with content as the following example::

	Antibody,Control
	CD3-0034,Mouse_IgG1
	CD4-0045,Mouse_IgG1
	... ...

where each line is a pair of Antibody name and the Control group name to which it is assigned. You should be able to get this information from your experiment setting or the original dataset.

Copy or upload them to ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir``.

-------------------------------------

2. Demultiplex Cell-Hashing Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	#. Prepare a sample sheet, ``cell_hashing_sample_sheet.csv``, with the following content::

		OUTNAME,RNA,ADT,TYPE
		exp,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/raw_feature_bc_matrix.h5,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cell_hashing.csv,cell-hashing

	   where **OUTNAME** specifies the subfolder and file names of output, which is free to change, **RNA** and **ADT** columns specify the RNA and ADT meta-data of samples, and **TYPE** is ``cell-hashing`` for this phase.

	   Then upload it to Google bucket::

	   	gsutil cp cell_hashing_sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/

	#. If your workspace doesn't have **cumulus_hashing_cite_seq** workflow, import it to your workspace by following Step 5 and 6 of `cumulus_hashing_cite_seq documentation <../hashing_cite_seq.html>`_.
	
	#. Prepare an input JSON file, ``cell_hashing_inputs.json`` with the following content to set up cumulus_hashing_cite_seq workflow inputs::

		{
			"cumulus_hashing_cite_seq.input_sample_sheet" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/cell_hashing_sample_sheet.csv",
			"cumulus_hashing_cite_seq.output_directory" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/",
			"cumulus_hashing_cite_seq.demuxEM_min_num_genes" : 500,
			"cumulus_hashing_cite_seq.demuxEM_generate_diagnostic_plots" : true
		}

	   For the details on these options, please refer to `cell-hashing/nuclei-hashing inputs`_.

	   .. _cell-hashing/nuclei-hashing inputs: ../hashing_cite_seq.html#cumulus-hashing-cite-seq-inputs

	#. On the page of cumulus_hashing_cite_seq workflow, upload ``cell_hashing_inputs.json`` by clicking ``upload json`` link. Save the inputs, and click ``RUN ANALYSIS`` button to start the job.

When the execution is done, you'll get a processed file, ``exp_demux.h5sc``, stored on cloud ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp/``.


----------------------------------------------------

3. Merge RNA and ADT Matrices for CITE-Seq Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	#. Prepare a sample sheet, ``cite_seq_sample_sheet.csv``, with the following content::

		OUTNAME,RNA,ADT,TYPE
		exp_raw,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp/exp_demux.h5sc,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cite_seq.csv,cite-seq

	   The structure of sample sheet here is the same as Phase 2. The difference is that you are now using the demultiplexed output ``h5sc`` file from Phase 2 as **RNA** here, and the sample **TYPE** is now ``cite-seq``.

	   Then upload it to Google bucket::

	   	gsutil cp cite_seq_sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/

	#. Prepare an input JSON file, ``cite_seq_inputs.json``, in the same directory as above, with the following content::

		{
			"cumulus_hashing_cite_seq.input_sample_sheet" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/cite_seq_sample_sheet.csv",
			"cumulus_hashing_cite_seq.output_directory" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/",
			"cumulus_hashing_cite_seq.antibody_control_csv" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/citeseq_antibody_control.csv"
		}

	   For the details on these options, please refer to `cell-hashing/nuclei-hashing inputs`_.

	#. On **cumulus_hashing_cite_seq** workflow page, clear all previous inputs, and then upload ``cite_seq_inputs.json`` by clicking ``upload json`` link. Save the new inputs, and click ``RUN ANALYSIS`` button to start the job.

When the execution is done, you'll get a merged raw matrices file, ``exp_raw.h5sc``, stored on cloud ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp_raw``.


-------------------

4. Data Analysis
^^^^^^^^^^^^^^^^^^^

	#. Prepare a sample sheet, ``cumulus_count_matrix.csv``, with the following content::

		Sample,Location
		exp,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp_raw/exp_raw.h5sc

	   This sample sheet describes the metadata for each 10x channel (as one row in the sheet). **Sample** specifies the name for each channel, which can be renamed; **Location** specifies the file location, which is the output of Phase 3.

	   Then upload it to Google bucket::

	   	gsutil cp cumulus_count_matrix.csv gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/

	   **Alternative**, if you have only one count matrix for analysis, which is the case here, you can skip this step. See `this manual`_ for input file formats that cumulus currently supports.

	#. If your workspace doesn't have **cumulus** workflow, import it to your workspace by following Step 2 and 3 of `cumulus documentation <../cumulus.html>`_.

	#. Prepare a JSON file, ``cumulus_inputs.json`` with the following content to set up **cumulus** workflow inputs::

		{
			"cumulus.input_file" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/cumulus_count_matrix.csv",
			"cumulus.output_name" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/results/exp_merged_out",
			"cumulus.num_cpu" : 8,
			"cumulus.select_only_singlets" : true,
			"cumulus.cite_seq" : true,
			"cumulus.run_louvain" : true,
			"cumulus.find_markers_lightgbm" : true,
			"cumulus.remove_ribo" : true,
			"cumulus.mwu" : true,
			"cumulus.annotate_cluster" : true,
			"cumulus.plot_fitsne" : "louvain_labels,assignment",
			"cumulus.plot_citeseq_fitsne" : "louvain_labels,assignment",
			"cumulus.plot_composition" : "louvain_labels:assignment"
		}

	   Alternatively, if you have only one count matrix for analysis and has skipped Step 1, directly set its location in ``cumulus.input_file`` parameter above. For this example, it is::

		{
			"cumulus.input_file" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp_raw/exp_raw.h5sc",
			... ...
		}

	   All the rest parameters remain the same.

	   Notice that for some file formats, ``cumulus.genome`` is required.

	   A typical cumulus pipeline consists of 4 steps, which is given here_. For the details of options above, please refer to `cumulus inputs`_.

	   .. _this manual: ../cumulus.html#prepare-input-data
	   .. _here: ../cumulus.html#cumulus-steps
	   .. _cumulus inputs: ../cumulus.html#global-inputs

	#. On the page of cumulus workflow, upload ``cumulus_inputs.json`` by clicking ``upload json`` link. Save the inputs, and click ``RUN ANALYSIS`` button to start the job.

When the execution is done, you'll get the following results stored on cloud ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/results/`` to check:
	
	* ``exp_merged_out.h5sc``: The aggregated count matrix data. This file doesn't exist if your ``cumulus.input_file`` parameter is not a sample sheet.
	* ``exp_merged_out.h5ad``: The processed RNA matrix data.
	* ``exp_merged_out.filt.xlsx``: The Quality-Control (QC) summary of the raw data.
	* ``exp_merged_out.filt.{UMI, gene, mito}.pdf``: The QC plots of the raw data.
	* ``exp_merged_out.de.xlsx``: Differential Expression analysis result.
	* ``exp_merged_out.markers.xlsx``: Result on cluster-specific markers predicted by gradient boosting machine.
	* ``exp_merged_out.anno.txt``: Cell type annotation output.
	* ``exp_merged_out.fitsne.pdf``: FIt-SNE plot.
	* ``exp_merged_out.citeseq.fitsne.pdf``: CITE-Seq FIt-SNE plot.
	* ``exp_merged_out.louvain_labels.assignment.composition.pdf``: Composition plot.

You can directly go to your Google Bucket to view or download these results.

----------------------

(optional) Run Terra Workflows in Command Line
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

For Phase 1, 2, and 3, besides uploading sample sheets and setting-up workflow inputs on workflow pages, you can also start the workflow execution via command line using **altocumulus** tool.

First, install *altocumulus* by following `altocumulus installation instruction <../command_line.html#install-altocumulus-for-non-broad-users>`_.

#. For Phase 1 above, when you are done with creating a sample sheet ``cellranger_sample_sheet.csv`` on your local machine, in the same directory, prepare JSON file ``cellranger_inputs.json`` as below::

	{
		"cellranger_workflow.input_csv_file" : "cellranger_sample_sheet.csv",
		... ...
	}

   where all the rest parameters remain the same as in Phase 1. Import **cellranger_workflow** workflow to your workspace as usual.

   Now run the following command in the same directory on your local machine::

   	alto fc_run -m cumulus/cellranger_workflow -w ws-lab/ws-01 --bucket-folder my-dir -i cellranger_input.json -o cellranger_input_updated.json

   Notice that if the execution failed, you could rerun the execution by setting ``cellranger_input_updated.json`` for ``-i`` option to use the sample sheet already uploaded to Google bucket. Similarly below.

#. For Phase 2 above, similarly, in the same directory of your ``cell_hashing_sample_sheet.csv`` file, prepare JSON file ``cell_hashing_inputs.json`` as below::

	{
		"cumulus_hashing_cite_seq.input_sample_sheet" : "cell_hashing_sample_sheet.csv",
		... ...
	}

   where all the rest parameters remain the same as in Phase 2. Import **cumulus_hashing_cite_seq** workflow to your workspace as usual.

   Run the following command in the same directory on your local machine::

	alto fc_run -m cumulus/cumulus_hashing_cite_seq -w ws-lab/ws-01 --bucket-folder my-dir -i cell_hashing_inputs.json -o cell_hashing_inputs_updated.json

#. For Phase 3 above, similarly, in the same directory of your ``cite_seq_sample_sheet.csv`` file, prepare JSON file ``cite_seq_inputs.json`` as below::

	{
		"cumulus_hashing_cite_seq.input_sample_sheet" : "cite_seq_sample_sheet.csv",
		... ...
	}

   where all the rest parameters remain the same as in Phase 3.

   Run the following command in the same directory on your local machine::

	alto fc_run -m cumulus/cumulus_hashing_cite_seq -w ws-lab/ws-01 --bucket-folder my-dir -i cite_seq_inputs.json -o cite_seq_inputs_updated.json

#. For Phase 4 above, similarly, in the same directory of your ``cumulus_count_matrix.csv`` file, prepare JSON file ``cumulus_inputs.json`` as below::

	{
		"cumulus.input_file" : "cumulus_count_matrix.csv",
		... ...
	}

   where all the rest parameters remain the same as in Phase 4. 

   **Alternatively**, if your input is not a sample sheet, simply set your ``cumulus_inputs.json`` as::

	{
		"cumulus.input_file" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp_raw/exp_raw.h5sc",
		... ...
	}

   where all the rest parameters remain the same.

   Run the following command in the same directory of your ``cumulus_inputs.json`` file::

	alto fc_run -m cumulus/cumulus -w ws-lab/ws-01 --bucket-folder my-dir/results -i cumulus_inputs.json -o cumulus_inputs_updated.json


.. _Terra: https://app.terra.bio/
.. _gsutil: https://cloud.google.com/storage/docs/gsutil