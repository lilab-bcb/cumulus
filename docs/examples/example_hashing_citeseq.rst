Example of Cell-Hashing and CITE-Seq Analysis on Cloud
++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this example, you'll learn how to perform Cell-Hashing and CITE-Seq analysis using **cumulus** on Terra_.

-----------------------------

0. Workspace and Data Preparation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After registering on Terra and creating a workspace there, you'll need the following two information:

	* **Terra workspace name**. This is shown on your Terra workspace webpage, with format *"<workspace-namespace>/<workspace-name>"*. Let it be ``ws-lab/ws-01`` in this example, which means that your workspace has namespace **ws-lab** and name **ws-01**.
	* The corresponding **Google Cloud Bucket** of your workspace. You can check it under **"Google Bucket"** title on the right panel of your Terra workspace's *Dashboard* tab. The bucket name associated with your workspace starts with ``fc-`` followed by a sequence of heximal numbers. In this example, let it be: ``gs://fc-e0000000``, where *"gs://"* is the head of Google bucket URL.

Then upload your BCL directories to Google bucket of your workspace using gsutil_::

	gsutil -m cp -r /my-local-path/BCL/* gs://fc-e0000000/data-source

where option ``-m`` means copy in parallel, ``-r`` means copy the directory recursively, ``/my-local-path/BCL`` is the path to the top-level directory of your BCL files on your local machine, and ``data-source`` is the folder on Google bucket to hold the uploaded data.

------------------------

1. Extract Gene-Count Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First step is to extract gene-count matrices from sequencing output. 


You need two original files from your dataset to start:

	* Cell-Hashing Index CSV file, say its filename is ``cell_hashing_index.csv``, of format "*feature_barcode,feature_name*". See an example below::

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

	gsutil -m cp -r /Users/foo/data-source gs://fc-e0000000/data-source

where ``gs://fc-e0000000/data-source`` is your working directory at cloud side, which can be changed at your will.

Next, create a sample sheet, ``cellranger_sample_sheet.csv``, for Cell Ranger processing. Below is an example::

	Sample,Reference,Flowcell,Lane,Index,DataType,FeatureBarcodeFile
	sample_control,GRCh38,gs://fc-e0000000/data-source,2,SI-GA-F1,rna
	sample_cc,GRCh38,gs://fc-e0000000/data-source,3,SI-GA-A1,rna
	sample_cell_hashing,GRCh38,gs://fc-e0000000/data-source,3,ATTACTCG,adt,cell_hashing_index.csv
	sample_cite_seq,GRCh38,gs://fc-e0000000/data-source,3,CGTGAT,adt,cite_seq_index.csv

For the details on how to prepare this sample sheet, please refer to Step 3 of `Cell Ranger sample sheet instruction`_.

	.. _Cell Ranger sample sheet instruction: ../cellranger.html#prepare-a-sample-sheet

When you are done with the sample sheet, upload it to Google bucket::

	gsutil cp cellranger_sample_sheet.csv gs://fc-e0000000/my-dir/

Now we are ready to set up **cellranger_workflow** workflow for this phase. If your workspace doesn't have this workflow, import it to your workspace by following `cellranger_workflow import instructions <../cellranger.html#import-cellranger-workflow>`_. 

Then prepare a JSON file, ``cellranger_inputs.json``, which is used to set up the workflow inputs::

	{
		"cellranger_workflow.input_csv_file" : "gs://fc-e0000000/my-dir/cellranger_sample_sheet.csv",
		"cellranger_workflow.output_directory" : "gs://fc-e0000000/my-dir"
	}

where ``gs://fc-e0000000/my-dir`` is the remote directory in which the output of cellranger_workflow will be generated. For the details on the options above, please refer to `Cell Ranger workflow inputs`_.

	.. _Cell Ranger workflow inputs: ../cellranger.html#workflow-input

When you are done with the JSON file, on cellranger_workflow workflow page, upload ``cellranger_inputs.json`` by clicking ``upload json`` link as below:

	.. image:: ../images/upload_json.png 
	   :scale: 70%

Then Click ``SAVE`` button to save the inputs, and click ``RUN ANALYSIS`` button as below to start the job:

	.. image:: ../images/run_analysis.png
	   :scale: 70%

When the execution is done, all the output results will be in folder ``gs://fc-e0000000/my-dir``. 

For the next phases, you'll need 3 files from the output:

	* RNA count matrix of the sample group of interest: ``gs://fc-e0000000/my-dir/sample_cc/raw_feature_bc_matrix.h5``;
	* Cell-Hashing Antibody count matrix: ``gs://fc-e0000000/my-dir/sample_cell_hashing/sample_cell_hashing.csv``;
	* CITE-Seq Antibody count matrix: ``gs://fc-e0000000/my-dir/sample_cite_seq/sample_cite_seq.csv``.

-------------------------------------

2. Demultiplex Cell-Hashing Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

	#. Prepare a sample sheet, ``demux_sample_sheet.csv``, with the following content::

		OUTNAME,RNA,TagFile,TYPE
		exp,gs://fc-e0000000/my-dir/raw_feature_bc_matrix.h5,gs://fc-e0000000/my-dir/sample_cell_hashing.csv,cell-hashing

	   where **OUTNAME** specifies the subfolder and file names of output, which is free to change, **RNA** and **TagFile** columns specify the RNA and hashing tag meta-data of samples, and **TYPE** is ``cell-hashing`` for this phase.

	   Then upload it to Google bucket::

	   	gsutil cp demux_sample_sheet.csv gs://fc-e0000000/my-dir/

	#. If your workspace doesn't have **demultiplexing** workflow, import it to your workspace by following Step 2 of `demultiplexing workflow preparation instructions <../demultiplexing.html#prepare-input-data-and-import-workflow>`_.
	
	#. Prepare an input JSON file, ``demux_inputs.json`` with the following content to set up cumulus_hashing_cite_seq workflow inputs::

		{
			"demultiplexing.input_sample_sheet" : "gs://fc-e0000000/my-dir/demultiplex_sample_sheet.csv",
			"demultiplexing.output_directory" : "gs://fc-e0000000/my-dir/"
		}

	   For the details on these options, please refer to `demultiplexing workflow inputs <../demultiplexing.html#workflow-inputs>`_.

	#. On the page of *demultiplexing* workflow, upload ``demux_inputs.json`` by clicking ``upload json`` link. Save the inputs, and click ``RUN ANALYSIS`` button to start the job.

When the execution is done, you'll get a processed file, ``exp_demux.zarr.zip``, stored on cloud in directory ``gs://fc-e0000000/my-dir/exp/``.


-----------------------------------

3. Data Analysis on CITE-Seq Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

In this step, we need to merge RNA and ADT matrices for CITE-Seq data, and perform the downstream analysis.

	#. Prepare a sample sheet, ``cumulus_count_matrix.csv``, with the following content::

		Sample,Location,Modality
		exp,gs://fc-e0000000/my-dir/exp/exp_demux.zarr.zip,rna
		exp,gs://fc-e0000000/my-dir/sample_cite_seq/sample_cite_seq.csv,citeseq

	   This sample sheet describes the metadata for each modality (as one row in the sheet): 
	   
	   	* **Sample** specifies the name of the modality, and all modalities must have *the same name*, as otherwise their count matrices won't be aggregated together;
		* **Location** specifies the file location. For RNA data, it's the output of Phase 2; for CITE-Seq antibody data, it's the output of Phase 1.
		* **Modality** specifies the modality type, which is either ``rna`` for RNA matrix, or ``citeseq`` for CITE-Seq antibody matrix.

	   Then upload it to Google bucket::

	   	gsutil cp cumulus_count_matrix.csv gs://fc-e0000000/my-dir/

	#. If your workspace doesn't have **cumulus** workflow, import it to your workspace by following Step 2 and 3 of `cumulus documentation <../cumulus.html>`_.

	#. Prepare a JSON file, ``cumulus_inputs.json`` with the following content to set up **cumulus** workflow inputs::

		{
			"cumulus.input_file" : "gs://fc-e0000000/my-dir/cumulus_count_matrix.csv",
			"cumulus.output_directory" : "gs://fc-e0000000/my-dir/results",
			"cumulus.output_name" : "exp_merged_out",
			"cumulus.select_only_singlets" : true,
			"cumulus.run_louvain" : true,
			"cumulus.run_umap" : true,
			"cumulus.citeseq" : true,
			"cumulus.citeseq_umap" : true,
			"cumulus.citeseq_umap_exclude" : "Mouse_IgG1,Mouse_IgG2a,Mouse_IgG2b,Rat_IgG2b",
			"cumulus.plot_composition" : "louvain_labels:assignment",
			"cumulus.plot_umap" : "louvain_labels,assignment",
			"cumulus.plot_citeseq_umap" : "louvain_labels,assignment",
			"cumulus.cluster_labels" : "louvain_labels",
			"cumulus.annotate_cluster" : true
		}

	   A typical cumulus pipeline consists of 4 steps, which is given here_. For the details of options above, please refer to `cumulus inputs`_.

	   .. _this manual: ../cumulus.html#prepare-input-data
	   .. _here: ../cumulus.html#cumulus-steps
	   .. _cumulus inputs: ../cumulus.html#global-inputs

	#. On the page of cumulus workflow, upload ``cumulus_inputs.json`` by clicking ``upload json`` link. Save the inputs, and click ``RUN ANALYSIS`` button to start the job.

When the execution is done, you'll get the following results stored on cloud ``gs://fc-e0000000/my-dir/results/exp_merged_out/`` to check:
	
	* ``exp_merged_out.aggr.zarr.zip``: The *ZARR* format file containing both the aggregated count matrix in ``<genome>-rna`` modality, as well as CITE-Seq antibody count matrix in ``<genome>-citeseq`` modality, where ``<genome>`` is the genome reference name of your count matrices, e.g. GRCh38.
	* ``exp_merged_out.zarr.zip``: The *ZARR* format file containing the analysis results in ``<genome>-rna`` modality, and CITE-Seq antibody count matrix in ``<genome>-citeseq`` modality.
	* ``exp_merged_out.<genome>-rna.h5ad``: The processed RNA matrix data in *H5AD* format.
	* ``exp_merged_out.<genome>-rna.filt.xlsx``: The Quality-Control (QC) summary of the raw data.
	* ``exp_merged_out.<genome>-rna.filt.{UMI, gene, mito}.pdf``: The QC plots of the raw data.
	* ``exp_merged_out.<genome>-rna.de.xlsx``: Differential Expression analysis result.
	* ``exp_merged_out.<genome>-rna.anno.txt``: Cell type annotation output.
	* ``exp_merged_out.<genome>-rna.umap.pdf``: UMAP plot.
	* ``exp_merged_out.<genome>-rna.citeseq.umap.pdf``: CITE-Seq UMAP plot.
	* ``exp_merged_out.<genome>-rna.louvain_labels.assignment.composition.pdf``: Composition plot.

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

   	alto run -m cumulus/cellranger_workflow -w ws-lab/ws-01 --bucket-folder my-dir -i cellranger_input.json

   Notice that if the execution failed, you could rerun the execution by setting ``cellranger_input_updated.json`` for ``-i`` option to use the sample sheet already uploaded to Google bucket. Similarly below.

#. For Phase 2 above, similarly, in the same directory of your ``demux_sample_sheet.csv`` file, prepare JSON file ``demux_inputs.json`` as below::

	{
		"demultiplexing.input_sample_sheet" : "demux_sample_sheet.csv",
		... ...
	}

   where all the rest parameters remain the same as in Phase 2. Import **demultiplexing** workflow to your workspace as usual.

   Run the following command in the same directory on your local machine::

	alto run -m cumulus/demultiplexing -w ws-lab/ws-01 --bucket-folder my-dir -i demux_inputs.json

#. For Phase 3 above, similarly, in the same directory of your ``cumulus_count_matrix.csv`` file, prepare JSON file ``cumulus_inputs.json`` as below::

	{
		"cumulus.input_file" : "cumulus_count_matrix.csv",
		... ...
	}

   where all the rest parameters remain the same as in Phase 3. 

   Run the following command in the same directory of your ``cumulus_inputs.json`` file::

	alto run -m cumulus/cumulus -w ws-lab/ws-01 --bucket-folder my-dir/results -i cumulus_inputs.json


.. _Terra: https://app.terra.bio/
.. _gsutil: https://cloud.google.com/storage/docs/gsutil