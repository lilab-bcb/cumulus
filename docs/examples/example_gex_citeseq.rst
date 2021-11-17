Example of Gene expression, CITE-Seq & Hashing Analysis on Cloud
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this example, you’ll learn how to perform Gene expression, CITE-Seq & Hashing data analysis on Cloud using Cromwell.

-------------------------------------------------------------------------------------------------------------------------

0. Prerequisite
^^^^^^^^^^^^^^^^^

You need to install the corresponding Cloud SDK tool on your local machine if not:

* `Google Cloud SDK <https://cloud.google.com/sdk/docs/install>`_ for Google Cloud.
* `AWS CLI v2 <https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html>`_ for Amazon AWS Cloud.

And then install `Altocumulus <https://github.com/lilab-bcb/altocumulus>`_ in your Python environment. This is the tool for data transfer between local machine and Cloud VM instance.

In this example, we assume that your Cromwell server is already deployed on Cloud at IP address ``10.0.0.0`` with port ``8000``, and also assume using Google Cloud with bucket ``gs://my-bucket``.


1. Extract Gene-Count Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First step is to extract gene-count matrices from sequencing output.


You need to prepare a feature index file from your dataset to start:

	* Feature barcode CSV file, say its filename is ``antibody_index.csv``, of format "*feature_barcode,feature_name,feature_type*". See an example below::

		TTCCTGCCATTACTA,HTO_1,hashing
                CCGTACCTCATTGTT,HTO_2,hashing
                GGTAGATGTCCTCAG,HTO_3,hashing
                TGGTGTCATTCTTGA,Ab1,hashing
                CTCATTGTAACTCCT,Ab2,citeseq
                GCGCAACTTGATGAT,Ab3,citeseq
		... ...

	  where each line is a pair of feature barcode, feature name and feature type of a sample.

Then upload antibody_index.csv to your Google Bucket using gsutil_. Assuming both files are in folder ``/Users/foo/data-source`` on your local machine, type the following command to upload::

	gsutil -m cp -r /Users/foo/data-source gs://fc-e0000000/data-source

where ``gs://fc-e0000000/data-source`` is your working directory at cloud side, which can be changed at your will.

Next, create a sample sheet, ``cellranger_sample_sheet.csv``, for Cell Ranger processing. Below is an example::

       Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
       sample_1_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/gex,fiveprime,rna
       sample_2_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/gex,fiveprime,rna
       sample_3_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/gex,fiveprime,rna
       sample_4_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/hashing_citeseq,fiveprime,citeseq,antibody_index.csv
       sample_5_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/hashing_citeseq,fiveprime,citeseq,antibody_index.csv
       sample_6_rna,GRCh38_v3.0.0,gs://fc-e0000000-0000-0000-0000-000000000000/VK18WBC6Z4/hashing_citeseq,fiveprime,citeseq,antibody_index.csv

For the details on how to prepare this sample sheet, please refer to Step 3 of `Cell Ranger sample sheet instruction`_.

	.. _Cell Ranger sample sheet instruction: ../cellranger/index.html#prepare-a-sample-sheet

When you are done with the sample sheet, upload it to Google bucket::

	gsutil cp cellranger_sample_sheet.csv gs://fc-e0000000/my-dir/

Now we are ready to set up **cellranger_workflow** workflow for this phase. If your workspace doesn't have this workflow, import it to your workspace by following `cellranger_workflow import instructions <../cellranger/index.html#import-cellranger-workflow>`_.

Then prepare a JSON file, ``cellranger_inputs.json``, which is used to set up the workflow inputs::

	{
		"cellranger_workflow.input_csv_file" : "gs://fc-e0000000/my-dir/cellranger_sample_sheet.csv",
		"cellranger_workflow.output_directory" : "gs://fc-e0000000/my-dir"
	}

where ``gs://fc-e0000000/my-dir`` is the remote directory in which the output of cellranger_workflow will be generated. For the details on the options above, please refer to `Cell Ranger workflow inputs`_.

	.. _Cell Ranger workflow inputs: ../cellranger/index.html#workflow-input

When the execution is done, all the output results will be in folder ``gs://fc-e0000000/my-dir``.

For the next phases, you'll need 2 files from the output:

	* RNA count matrix of the sample group of interest: ``gs://fc-e0000000/my-dir/sample_cc/raw_feature_bc_matrix.h5``;
	* Cell-Hashing Antibody count matrix: ``gs://fc-e0000000/my-dir/sample_cell_hashing/sample_cell_hashing.csv``;


2. Demultiplex Cell-Hashing Data using DemuxEM
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

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


3. Extract Demultiplexing results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


      #. Load Libraries::

              import numpy as np
              import pandas as pd
              import pegasus as pg
              import matplotlib.pyplot as plt
              import seaborn as sns

      #. Load demuxEM output. For demuxEM, load RNA expression matrix with demultiplexed sample identities in Zarr format. These can be found in Google cloud console. QC 500 <= # of genes < 6000, % mito <= 10%::

              data = pg.read_input('exp_demux.zarr.zip')
              pg.qc_metrics(data, min_genes=500, max_genes=6000, mito_prefix='Mt-', percent_mito=10)
              pg.filter_data(data)

      #. Demultiplexing results showing singlets, doublets and unknown::

              data.obs['demux_type'].value_counts()

      #. Show assignments in singlets::

              idx = data.obs['demux_type'] == 'singlet'
              data.obs.loc[idx, 'assignment'].value_counts()[0:10]

      #. Write assignment outputs to CSV::

              data.obs[['demux_type', 'assignment']].to_csv('demux_exp.csv')


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
