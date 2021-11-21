Example of Gene expression, Hashing and CITE-Seq Analysis on Cloud
--------------------------------------------------------------------

In this example, you’ll learn how to perform Gene expression, Hashing and CITE-Seq data analysis on Cloud.

This example covers the cases of both Terra_ platform and a custom cloud server running Cromwell_. When reading through the tutorial, you may check out the corresponding part based on your working situation.

-------------------------------------------------------------------------------------------------------------------------

0. Prerequisite
~~~~~~~~~~~~~~~~~

0-a. Cromwell server
^^^^^^^^^^^^^^^^^^^^^

If you use a Cromwell server on Cloud, on your local machine, you need to install the corresponding Cloud SDK tool if not:

	* `Google Cloud SDK`_ if your Cloud bucket is on Google Cloud.
	* `AWS CLI v2`_ if your Cloud bucket is on Amazon AWS Cloud.

And then install Altocumulus_ in your Python environment. This is the tool for data transfer between local machine and cloud bucket, as well as communication with the Cromwell server on cloud.


0-b. Terra Platform
^^^^^^^^^^^^^^^^^^^^^

If you use Terra, after registering on Terra and creating a workspace there, you'll need the following information:

	* **Terra workspace name**. This is shown on your Terra workspace webpage, with format "*<workspace-namespace>/<workspace-name>*". For example, if your Terra workspace has full name ``ws-lab/ws-01``, then **ws-lab** is the namespace and **ws-01** is the workspace name winthin that namespace.
	* The corresponding **Google Cloud Bucket** of your Terra workspace. You can check it under "*Google Bucket*" title on the right panel of your Terra workspace’s *Dashboard* tab. The bucket name associated with your workspace starts with ``fc-`` followed by a sequence of heximal numbers. For example, ``gs://fc-e0000000``, where "*gs://*" is the header of GS URI.

Besides, install `Google Cloud SDK`_ and Altocumulus_ on your local machine for data uploading. These tools will be used for data transfer between local machine and Cloud bucket.

Alternatively, you can also use Terra web UI for job submission instead of command-line submission. This will be discussed in Section `Run Analysis with Terra Web UI`_ below.

------------------------------------

1. Extract Gene-Count Matrices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This phase is to extract gene-count matrices from sequencing output.

There are two cases: (1) from BCL data, which includes *mkfastq* step to generate FASTQ files and *count* step to generate gene-count matrices; (2) from FASTQ files, which only runs the *count* step.

1-a. Extract Genen-Count Matrices from BCL data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section covers the case starting from BCL data.

Step 1. Sample Sheet Preparation
++++++++++++++++++++++++++++++++++

First, prepare a feature index file for your dataset. Say its filename is ``antibody_index.csv``, which has format "*feature_barcode, feature_name,feature_type*". See an example below::

	TTCCTGCCATTACTA,HTO_1,hashing
	CCGTACCTCATTGTT,HTO_2,hashing
	GGTAGATGTCCTCAG,HTO_3,hashing
	TGGTGTCATTCTTGA,Ab1,citeseq
	CTCATTGTAACTCCT,Ab2,citeseq
	GCGCAACTTGATGAT,Ab3,citeseq
	... ...

where each line contains the barcode and the name of a Hashing/CITE-Seq index: ``hashing`` indicates a Cell/Nucleus-Hashing index, while ``citeseq`` indicates a CITE-Seq index.

Next, create a sample sheet ``cellranger_sample_sheet.csv`` for Cell Ranger processing on your local machine. Below is an example::

	Sample,Reference,Flowcell,Lane,Index,Chemistry,DataType,FeatureBarcodeFile
	sample_control,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-TT-A1,fiveprime,rna
	sample_gex,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-TT-A2,fiveprime,rna
	sample_cell_hashing,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-NN-A1,fiveprime,hashing,/path/to/antibody_index.csv
	sample_cite_seq,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-NN-A2,fiveprime,citeseq,/path/to/antibody_index.csv

where

	* ``GRCh38-2020-A`` is the is the Human GRCh38 (GENCODE v32/Ensembl 98) genome reference prebuilt by Cumulus. See `Cumulus single-cell genome reference list`_ for a complete list of genome references.
	* ``/path/to/flowcell/folder`` should be replaced by the actual local path to the BCL folder of your sequencing data.
	* ``/path/to/antibody_index.csv`` should be replaced by the actual local path to ``antibody_index.csv`` file we just created above.
	* ``rna``, ``hashing`` and ``citeseq`` refer to gene expression data, cell/nucleus-hashing data, and CITE-Seq data, respectively.
	* Samples of type ``rna`` do not need any feature barcode file for indexing.

For the details on how to prepare this sample sheet, please refer to Step 3 of `Cell Ranger sample sheet instruction`_.

Step 2. Workflow Input Preparation
++++++++++++++++++++++++++++++++++++

Now prepare a JSON file for **cellranger_workflow** WDL workflow input on your local machine (say named ``cellranger_inputs.json``)::

	{
		"cellranger_workflow.input_csv_file": "/path/to/cellranger_sample_sheet.csv",
		"cellranger_workflow.output_directory": "gs://my-bucket/cellranger_output"
	}

where

	* ``/path/to/cellranger_sample_sheet.csv`` should be replaced by the actual local path to your sample sheet created above.
	* ``gs://my-bucket/cellranger_output`` is the target folder on Google bucket to store your result when the workflow job is finished, where ``my-bucket`` should be replaced by your own Google bucket name.

For details on the all the workflow inputs of *cellranger_workflow*, please refer to `Cell Ranger workflow inputs`_.

Step 3. Job Submission
++++++++++++++++++++++++

Now we are ready to submit a job to cloud for computing:

* If you use a Cromwell server on cloud, run the following Altocumulus command::

	alto cromwell run -s <server-address> -p <port-number> -m broadinstitute:cumulus:cellranger -i /path/to/cellranger_inputs.json -o cellranger_inputs_updated.json -b gs://my-bucket/data_source

where

	* ``-s`` specifies the server's IP address (or hostname), where ``<server-address>`` should be replaced by the actual IP address (or hostname).
	* ``-p`` specifies the server's port number for Cromwell, where ``<port-number>`` should be replaced by the actual port number.
	* ``-m`` specifies which WDL workflow to use. You should use the Dockstore name of Cumulus cellranger_workflow_. Here, the version is omitted, so that the default version will be used. Alternatively, you can explicitly specify which version to use, e.g. ``broadinstitute:cumulus:cellranger:master`` to use its development version in *master* branch.
	* ``-i`` specifies the workflow input JSON file.
	* ``-o`` and ``-b`` are used when the input data (which are specified in the workflow input JSON file and sample sheet CSV file) are local and need to be uploaded to Cloud bucket first.
	* ``-o`` specifies the updated workflow input JSON file after uploading the input data, with all the local paths updated to Cloud bucket URIs. This is useful when resubmitting jobs running the same input data, without uploading the same input data again.
	* ``-b`` specifies which folder on Cloud bucket to upload the local input data, where ``my-bucket`` should be replaced by your own Google bucket name. Feel free to choose the folder name other than ``data_source``.

Notice that ``-o`` and ``-b`` options can be dropped if all of your input data are already on Cloud bucket.

After submission, you'll get the job's ID for tracking its status::

	alto cromwell check_status -s <server-address> -p <port-number> --id <your-job-ID>

where ``<your-job-ID>`` should be replaced by the actual Cromwell job ID.

* If you use Terra, run the following Altocumulus command::

	alto terra run -m broadinstitute:cumulus:cellranger -w ws-lab/ws-01 --bucket-folder data_source -i /path/to/cellranger_inputs.json -o cellranger_inputs_updated.json

where

	* ``-m`` specifies which WDL workflow to use. You should use the Dockstore name of Cumulus cellranger_workflow_. Here, the version is omitted, so that the default version will be used. Alternatively, you can explicitly specify which version to use, e.g. ``broadinstitute:cumulus:cellranger:master`` to use its development version in *master* branch.
	* ``-w`` specifies the Terra workspace full name to use, where ``ws-lab/ws-01`` should be replaced by your own Terra workspace full name.
	* ``--bucket-folder`` specifies the folder name on the Google bucket associated with the Terra workspace to store the uploaded data. Feel free to choose folder name other than ``data_source``.
	* ``-i`` specifies the workflow input JSON file, where ``/path/to/cellranger_inputs.json`` should be replaced by the actual local path to ``cellranger_inputs.json`` file.
	* ``-o`` specifies the updated workflow input JSON file after uploading the input data, with all the local paths updated to Cloud bucket URIs. This is useful when resubmitting jobs running the same input data, without uploading the same input data again.

Notice that ``--bucket-folder`` and ``-o`` options can be dropped if all of your input data are already on Cloud bucket.

After submission, you can check the job's status in the *Job History* tab of your Terra workspace page.

When the job is done, you'll get results in ``gs://my-bucket/cellranger_output``, which is specified in ``cellranger_inputs.json`` above. It should contain 4 subfolders, each of which is associated with one sample specified in ``cellranger_sample_sheet.csv`` above.

For the next phases, you'll need 3 files from the output:

	* RNA count matrix of the sample group of interest: ``gs://my-bucket/cellranger_output/sample_gex/raw_feature_bc_matrix.h5``;
	* Cell-Hashing Antibody count matrix: ``gs://my-bucket/cellranger_output/sample_cell_hashing/sample_cell_hashing.csv``;
	* CITE-Seq Antibody count matrix: ``gs://my-bucket/cellranger_output/sample_cite_seq/sample_cite_seq.csv``.


1-b. Extract Gene-Cound Matrices from FASTQ files
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This section covers the case starting from FASTQ files.

Similarly as above, First, prepare a feature index file for your dataset. Say its filename is ``antibody_index.csv``, which has format "*feature_barcode, feature_name,feature_type*". See an example below::

	TTCCTGCCATTACTA,HTO_1,hashing
	CCGTACCTCATTGTT,HTO_2,hashing
	GGTAGATGTCCTCAG,HTO_3,hashing
	TGGTGTCATTCTTGA,Ab1,citeseq
	CTCATTGTAACTCCT,Ab2,citeseq
	GCGCAACTTGATGAT,Ab3,citeseq
	... ...

where each line contains the barcode and the name of a Hashing/CITE-Seq index: ``hashing`` indicates a Cell/Nucleus-Hashing index, while ``citeseq`` indicates a CITE-Seq index.

Next, create a sample sheet ``cellranger_sample_sheet.csv`` for Cell Ranger processing on your local machine. Below is an example::

	Sample,Reference,Flowcell,Chemistry,DataType,FeatureBarcodeFile
	sample_1_rna,GRCh38-2020-A,/path/to/fastq/gex,fiveprime,rna
	sample_2_rna,GRCh38-2020-A,/path/to/fastq/gex,fiveprime,rna
	sample_3_rna,GRCh38-2020-A,/path/to/fastq/gex,fiveprime,rna
	sample_1_adt,GRCh38-2020-A,/path/to/fastq/hashing_citeseq,fiveprime,adt,/path/to/antibody_index.csv
	sample_2_adt,GRCh38-2020-A,/path/to/fastq/hashing_citeseq,fiveprime,adt,/path/to/antibody_index.csv
	sample_3_adt,GRCh38-2020-A,/path/to/fastq/hashing_citeseq,fiveprime,adt,/path/to/antibody_index.csv

where

	* ``GRCh38-2020-A`` is the is the Human GRCh38 (GENCODE v32/Ensembl 98) genome reference prebuilt by Cumulus. See `Cumulus single-cell genome reference list`_ for a complete list of genome references.
	* ``/path/to/fastq/gex`` should be replaced by the actual local path to the folder containing FASTQ files of RNA samples.
	* ``/path/to/fastq/hashing_citeseq`` should be replaced by the actual local path to the folder containing FASTQ files of Cell/Nucleus-Hashing and CITE-Seq samples.
	* ``/path/to/antibody_index.csv`` should be replaced by the actual local path to ``antibody_index.csv`` file we just created above.
	* ``rna`` and ``adt`` refer to gene expression data and antibody data, respectively. In specific, ``adt`` covers both ``citeseq`` and ``hashing`` types, i.e. it includes both Hashing and CITE-Seq data types.
	* Samples of type ``rna`` do not need any feature barcode file for indexing.
	* Columns *Lane* and *Index* are not needed if starting from FASTQ files, as *mkfastq* step will be skipped.

For the details on how to prepare this sample sheet, please refer to Step 3 of `Cell Ranger sample sheet instruction`_.

Now prepare a JSON file for **cellranger_workflow** WDL workflow input on your local machine (say named ``cellranger_inputs.json``)::

	{
		"cellranger_workflow.input_csv_file": "/path/to/cellranger_sample_sheet.csv",
		"cellranger_workflow.output_directory": "gs://my-bucket/cellranger_output",
		"cellranger_workflow.run_mkfastq": false
	}

where

	* ``/path/to/cellranger_sample_sheet.csv`` should be replaced by the actual local path to your sample sheet created above.
	* ``gs://my-bucket/cellranger_output`` is the target folder on Google bucket to store your result when the workflow job is finished, where ``my-bucket`` should be replaced by your own Google bucket name.
	* Set *run_mkfastq* to ``false`` to skip the *mkfastq* step, as we start from FASTQ files.

For details on the all the workflow inputs of *cellranger_workflow*, please refer to `Cell Ranger workflow inputs`_.

Now we are ready to submit a job to cloud for computing. Follow instructions in `Section 1-a <./example_hashing_citeseq.html#step-3-job-submission>`_ above.

When finished, you'll get results in ``gs://my-bucket/cellranger_output``, which is specified in ``cellranger_inputs.json`` above. It should contain 6 subfolders, each of which is associated with one sample specified in ``cellranger_sample_sheet.csv`` above.

In specific, for each ``adt`` type sample, there are both count matrix of Hashing data and that of CITE-Seq data generated inside its corresponding subfolder, with filename suffix ``.hashing.csv`` and ``.citeseq.csv``, respectively.

-------------------------------------------------------

2. Demultiplex Cell-Hashing Data using DemuxEM
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Run Workflow on Cloud
^^^^^^^^^^^^^^^^^^^^^^^^

Next, we need to demultiplex the resulting RNA gene-count matrices. We use `DemuxEM <https://demuxem.readthedocs.io/>`_ method in this example.

To be brief, we use the output of Section 1-a for illustration:

	#. On your local machine, prepare a CSV-format sample sheet ``demux_sample_sheet.csv`` with the following content::

		OUTNAME,RNA,TagFile,TYPE
		exp,gs://my-bucket/cellranger_output/sample_gex/raw_feature_bc_matrix.h5,gs://my-bucket/cellranger_output/sample_cell_hashing/sample_cell_hashing.csv,cell-hashing

	   where **OUTNAME** specifies the subfolder and file names of output, which is free to be changed, **RNA** and **TagFile** columns specify the RNA and hashing tag meta-data of samples, and **TYPE** is ``cell-hashing`` for this phase.

	#. On your local machine, also prepare an input JSON file ``demux_inputs.json`` for **demultiplexing** WDL workflow, ``demux_inputs.json`` with the following content::

		{
			"demultiplexing.input_sample_sheet" : "/path/to/demux_sample_sheet.csv",
			"demultiplexing.output_directory" : "gs://my-bucket/demux_output"
		}

	   where ``/path/to/demux_sample_sheet.csv`` should be replaced by the actual local path to ``demux_sample_sheet.csv`` created above.

	   For the details on these options, please refer to `demultiplexing workflow inputs <../demultiplexing.html#workflow-inputs>`_.

	#. Submit a *demultiplexing* job with ``demux_inputs.json`` input above to cloud for execution.

For job submission:

* If you use a Cromwell server on cloud, run the following Altocumulus command on your local machine::

	alto cromwell run -s <server-address> -p <port-number> -m broadinstitute:cumulus:demultiplexing -i /path/to/demux_inputs.json -o demux_inputs_updated.json -b gs://my-bucket/data_source

where

	* ``broadinstitute:cumulus:demultiplexing`` refers to demultiplexing_ WDL workflow published on Dockstore. Here, the version is omitted, so that the default version will be used. Alternatively, you can explicitly specify which version to use, e.g. ``broadinstitute:cumulus:demultiplexing:master`` to use its development version in *master* branch.
	* ``/path/to/demux_inputs.json`` should be replaced by the actual local path to ``demux_inputs.json`` created above.
	* Replace ``my-bucket`` in ``-b`` option by your own Google bucket name, and feel free to choose folder name other than ``data_source`` for uploading.
	* We still need ``-o`` and ``-b`` options because ``demux_sample_sheet.csv`` is on the local machine.

Similarly, when the submission succeeds, you'll get another job ID for demultiplexing. You can use it to track the job status.

* If you use Terra, run the following Altocumulus command::

	alto terra run -m broadinstitute:cumulus:demultiplexing -w ws-lab/ws-01 --bucket-folder data_source -i /path/to/demux_inputs.json -o demux_inputs_updated.json

where

	* ``broadinstitute:cumulus:demultiplexing`` refers to demultiplexing_ WDL workflow published on Dockstore. Here, the version is omitted, so that the default version will be used. Alternatively, you can explicitly specify which version to use, e.g. ``broadinstitute:cumulus:demultiplexing:master`` to use its development version in *master* branch.
	* ``/path/to/demux_inputs.json`` should be replaced by the actual local path to ``demux_inputs.json`` created above.
	* ``ws-lab/ws-01`` should be replaced by your own Terra workspace full name.
	* ``--bucket-folder``: Feel free to choose folder name other than ``data_source`` for uploading.
	* We still need ``-o`` and ``--bucket-folder`` options because ``demux_sample_sheet.csv`` is on the local machine.

After submission, you can check the job's status in the *Job History* tab of your Terra workspace page.

When finished, demultiplexing results are in ``gs://my-bucket/demux_output/exp`` folder, with the following important output files:

	* ``exp_demux.zarr.zip``: Demultiplexed RNA raw count matrix. This will be used for downstram analysis.
	* ``exp.out.demuxEM.zarr.zip``: This file contains intermediate results for both RNA and hashing count matrices, which is useful for compare with other demultiplexing methods.
	* DemuxEM plots in PDF format. They are used for evaluating the performance of DemuxEM on the data.

(Optional) Extract Demultiplexing results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

This is performed on your local machine with demultiplexing results downloaded from cloud to your machine.

To download the demultiplexed count matrix ``exp_demux.zarr.zip``, you can either do it in Google cloud console, or using gsutil_ in command line::

	gsutil -m cp gs://my-bucket/demux_output/exp/exp_demux.zarr.zip .

After that, in your Python environment, install Pegasus_ package, and follow the steps below to extract the demultiplexing results:

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

----------------------------------------

3. Data Analysis on CITE-Seq Data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this phase, we merge RNA and ADT matrices for CITE-Seq data, and perform the downstream analysis.

To be brief, we use the CITE-Seq count matrix generated from Section 1-a and demultiplexing results in Section 2 for illustraion here:

	1. On your local machine, prepare a CSV-format sample sheet ``count_matrix.csv`` with the following content::

		Sample,Location,Modality
		exp,gs://my-bucket/demux_output/exp/exp_demux.zarr.zip,rna
		exp,gs://my-bucket/cellranger_output/sample_cite_seq/sample_cite_seq.csv,citeseq

	   This sample sheet describes the metadata for each modality (as one row in the sheet):

	   	* **Sample** specifies the name of the modality, and all the modalities of the same sample should have one common name, as otherwise their count matrices won't be aggregated together;
	   	* **Location** specifies the file location. For RNA data, this is the output of Phase 2; for CITE-Seq antibody data, it's the output of Phase 1.
	   	* **Modality** specifies the modality type, which is either ``rna`` for RNA matrix, or ``citeseq`` for CITE-Seq antibody matrix.

	2. On your local machine, also prepare a JSON file ``cumulus_inputs.json`` for **cumulus** WDL workflow, with the following content::

		{
			"cumulus.input_file": "/path/to/count_matrix.csv",
			"cumulus.output_directory": "gs://my-bucket/cumulus_output",
			"cumulus.output_name": "exp_merged_out",
			"cumulus.select_only_singlets": true,
			"cumulus.run_louvain": true,
			"cumulus.run_umap": true,
			"cumulus.citeseq": true,
			"cumulus.citeseq_umap": true,
			"cumulus.citeseq_umap_exclude": "Mouse_IgG1,Mouse_IgG2a,Mouse_IgG2b,Rat_IgG2b",
			"cumulus.plot_composition": "louvain_labels:assignment",
			"cumulus.plot_umap": "louvain_labels,assignment",
			"cumulus.plot_citeseq_umap": "louvain_labels,assignment",
			"cumulus.cluster_labels": "louvain_labels",
			"cumulus.annotate_cluster": true,
			"cumulus.organism": "human_immune"
		}

	   where ``/path/to/count_matrix.csv`` should be replaced by the actual local path to ``count_matrix.csv`` created above.

	   A typical Cumulus WDL pipeline consists of 4 steps, which is given `here <../cumulus.html#cumulus-steps>`_. For details on Cumulus workflow inputs above, please refer to `cumulus inputs`_.

	3. Submit a *demultiplexing* job with ``cumulus_inputs.json`` input above to cloud for execution.

For job submission:

* If you use a Cromwell server on cloud, run the following Altocumulus command to submit the job::

	alto cromwell run -s <server-address> -p <port-number> -m broadinstitute:cumulus:cumulus -i /path/to/cumulus_inputs.json -o cumulus_inputs_updated.json -b gs://my-bucket/data_source

where

	* ``broadinstitute:cumulus:cumulus`` refers to cumulus_ WDL workflow published on Dockstore. Here, the version is omitted, so that the default version will be used. Alternatively, you can explicitly specify which version to use, e.g. ``broadinstitute:cumulus:cumulus:master`` to use its development version in *master* branch.
	* ``/path/to/cumulus_inputs.json`` should be replaced by the actual local path to ``cumulus_inputs.json`` created above.
	* ``my-bucket`` in ``-b`` option should be replaced by your own Google bucket name, and feel free to choose folder name other than ``data_source`` for uploading data.
	* We still need ``-o`` and ``-b`` options because ``count_matrix.csv`` is on the local machine.

Similarly, when the submission succeeds, you'll get another job ID for demultiplexing. You can use it to track the job status.

* If you use Terra, run the following Altocumulus command::

	alto terra run -m broadinstitute:cumulus:cumulus -w ws-lab/ws-01 --bucket-folder data_source -i /path/to/cumulus_inputs.json -o cumulus_inputs_updated.json

where

	* ``broadinstitute:cumulus:cumulus`` refers to cumulus_ WDL workflow published on Dockstore. Here, the version is omitted, so that the default version will be used. Alternatively, you can explicitly specify which version to use, e.g. ``broadinstitute:cumulus:cumulus:master`` to use its development version in *master* branch.
	* ``ws-lab/ws-01`` should be replaced by your own Terra workspace full name.
	* ``--bucket-folder``: Feel free to choose folder name other than ``data_source`` for uploading data.
	* ``/path/to/cumulus_inputs.json`` should be replaced by the actual local path to ``cumulus_inputs.json`` created above.
	* We still need ``-o`` and ``--bucket-folder`` options because ``count_matrix.csv`` is on the local machine.

After submission, you can check the job's status in the *Job History* tab of your Terra workspace page.

When finished, all the output files are in ``gs://my-bucket/cumulus_output`` folder, with the following important files:

	* ``exp_merged_out.aggr.zarr.zip``: The *ZARR* format file containing both the aggregated count matrix in ``<genome>-rna`` modality, as well as CITE-Seq antibody count matrix in ``<genome>-citeseq`` modality, where ``<genome>`` is the genome reference name of your count matrices, e.g. *GRCh38-2020-A*.
	* ``exp_merged_out.zarr.zip``: The *ZARR* format file containing the analysis results in ``<genome>-rna`` modality, and CITE-Seq antibody count matrix in ``<genome>-citeseq`` modality.
	* ``exp_merged_out.<genome>-rna.h5ad``: The processed RNA matrix data in *H5AD* format.
	* ``exp_merged_out.<genome>-rna.filt.xlsx``: The Quality-Control (QC) summary of the raw data.
	* ``exp_merged_out.<genome>-rna.filt.{UMI, gene, mito}.pdf``: The QC plots of the raw data.
	* ``exp_merged_out.<genome>-rna.de.xlsx``: Differential Expression analysis result.
	* ``exp_merged_out.<genome>-rna.anno.txt``: The putative cell type annotation output.
	* ``exp_merged_out.<genome>-rna.umap.pdf``: UMAP plot.
	* ``exp_merged_out.<genome>-rna.citeseq.umap.pdf``: CITE-Seq UMAP plot.
	* ``exp_merged_out.<genome>-rna.louvain_labels.assignment.composition.pdf``: Composition plot.

------------------------------

Run Analysis with Terra Web UI
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

For Terra users, instead of using Altocumulus to submit jobs in command line, they can also use the Terra web UI.

First, upload the local BCL data or FASTQ files to the Google bucket associated with your Terra workspace (say ``gs://fc-e0000000``) using gsutil_::

	gsutil -m cp -r /path/to/your/data/folder gs://fc-e0000000/data_source/

where ``/path/to/your/data/folder`` should be replaced by the actual local path to your data folder, and ``data_source`` is the folder on Google bucket to store the uploaded data.

Then for each of the 3 phases above:

	1. When preparing the sample sheet, remember to replace all the local paths by the GS URIs of the corresponding folders/files that you uploaded to Google bucket. Then upload it to Google bucket as well::

		gsutil cp /path/to/sample/sheet gs://fc-e0000000/data_source/

	   where ``/path/to/sample/sheet`` should be replaced by the actual local path to your sample sheet.
	   Notice that for Phase 1, ``antibody_index.csv`` file should also be uploaded to Google bucket, and its references in the sample sheet must be replaced by its GS URI.
	2. When preparing the workflow input JSON file, change the field of sample sheet to its GS URI on cloud.
	3. Import the corresponding WDL workflow to your Terra workspace by following steps in `import workflows`_ tutorial.
	4. In the workflow page (Workspace -> Workflows -> your WDL workflow), upload your input JSON file by clicking the "*upload json*" button:

		.. image:: ../images/upload_json.png
		   :scale: 70%

	5. Click "*SAVE*" button to save the configuration, and click "*RUN ANALYSIS*" button to submit the job:

		.. image:: ../images/run_analysis.png
		   :scale: 70%

You can check the job's status in the *Job History* tab of your Terra workspace page.


.. _Terra: https://app.terra.bio/
.. _Cromwell: https://cromwell.readthedocs.io
.. _Google Cloud SDK: https://cloud.google.com/sdk/docs/install
.. _AWS CLI v2: https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html
.. _Altocumulus: https://github.com/lilab-bcb/altocumulus
.. _Cell Ranger sample sheet instruction: ../cellranger/index.html#prepare-a-sample-sheet
.. _Cell Ranger workflow inputs: ../cellranger/index.html#workflow-input
.. _cellranger_workflow: https://dockstore.org/workflows/github.com/klarman-cell-observatory/cumulus/Cellranger
.. _demultiplexing: https://dockstore.org/workflows/github.com/klarman-cell-observatory/cumulus/Demultiplexing
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _Cumulus single-cell genome reference list: ../cellranger/index.html#sample-sheet
.. _cumulus inputs: ../cumulus.html#global-inputs
.. _cumulus: https://dockstore.org/workflows/github.com/klarman-cell-observatory/cumulus/Cumulus
.. _import workflows: ../cumulus_import.html
.. _Run Analysis with Terra Web UI: ./example_hashing_citeseq.html#run-analysis-with-terra-web-ui
.. _Pegasus: https://pegasus.readthedocs.io
