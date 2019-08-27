Example of Cell-Hashing and CITE-Seq Analysis on Cloud
++++++++++++++++++++++++++++++++++++++++++++++++++++++

In this example, you'll learn how to perform Cell-Hashing and CITE-Seq analysis using **scCloud** on Terra_.

-----------------------------

0. Create Terra Workspace
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

After registering on Terra and creating a workspace there, you'll need the following two information:

	* **Terra workspace name**. This is shown on your Terra workspace webpage, with format *"<workspace-namespace>/<workspace-name>"*. Let it be ``ws-lab/ws-01`` in this example, which means that your workspace has namespace **ws-lab** and name **ws-01**.
	* The corresponding **Google Cloud Bucket location**. You can check it by clicking the link under **"Google Bucket"** title on your Terra workspace webpage. Let it be ``gs://fc-e0000000-0000-0000-0000-000000000000`` in this example.



------------------------

1. Run Cell Ranger Tools
^^^^^^^^^^^^^^^^^^^^^^^^

You'll need two original files from your dataset to start the process:

	* Cell-Hashing Index CSV file, say its filename is ``cell_hashing_index.csv``, with content as the example below::

		AATCATCACAAGAAA,CB1
		GGTCACTGTTACGTA,CB2
		... ...

	  where each line is a pair of Barcode and [TODO]

	* CITE-Seq Index CSV file, say its filename is ``cite_seq_index.csv``, with content as the example below::

		TTACATGCATTACGA,CD19
		GCATTAGCATGCAGC,HLA-ABC
		... ...

	  where each line is a pair of Barcode and Specificity of an Antibody.

Then upload them to your Google Bucket. Assuming both files are in folder ``/Users/foo/data-source`` on your local machine, then type the following command to upload::

	gsutil -m cp -r /Users/foo/data-source gs://fc-e0000000-0000-0000-0000-000000000000/data-source

where option ``-m`` means copy in parallel, ``-r`` means copy the directory recursively, and ``gs://fc-e0000000-0000-0000-0000-000000000000/data-source`` is your working directory at cloud side, which can be changed at your will.

Next, create a sample sheet, ``ranger_sample_sheet.csv``, on your local machine with content below::

	Sample,Reference,Flowcell,Lane,Index,DataType,FeatureBarcodeFile
	sample_control,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,2,SI-GA-F1,rna
	sample_cc,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,3,SI-GA-A1,rna
	sample_cell_hashing,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,3,ATTACTCG,adt,cell_hashing_index.csv
	sample_cite_seq,GRCh38,gs://fc-e0000000-0000-0000-0000-000000000000/data-source,3,CGTGAT,adt,cite_seq_index.csv

For the details on how to prepare this sample sheet, please refer to Step 2 of `Cell-Ranger sample sheet instruction`_.

	.. _Cell-Ranger sample sheet instruction: ../cellranger.html

Moreover, prepare a JSON file, ``ranger_input.json``, in the same directory to set up the workflow::

	{
		"cellranger_workflow.input_csv_file" : "ranger_sample_sheet.csv",
		"cellranger_workflow.output_directory" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir"
	}

where ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir`` is the remote directory in which the output will be generated. For the detaills on the options above, please refer to `Cell-Ranger workflow inputs`_.

	.. _Cell-Ranger workflow inputs: ../cellranger.html#cellranger-workflow-inputs

Now run the following sccutil command in the same directory on your local machine::

	sccutil fc_run -m scCloud/cellranger_workflow -w ws-lab/ws-01 --bucket-folder my-dir -i ranger_input.json -o ranger_input_updated.json

where options

	   	* ``-m`` specifies the method to be used, 
	   	* ``-w`` specifies the Terra workspace name,
	   	* ``--bucket-folder`` specifies the working directory on your Google Bucket, 
	   	* ``-i`` specifies which input JSON file to be used, 
	   	* ``-o`` specifies the filename of updated JSON file after execution, which can be renamed at your will.

When the execution is done, all the output results will be in folder ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir``. 

You'll need 4 files for the next phases. 3 are from the output:

	* RNA information on the sample group of interest: ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cc/raw_gene_bc_matrices_h5.h5``;
	* Antibody Cell-Hashing matrix: ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cell_hashing/sample_cell_hashing.csv``;
	* Antibody CITE-Seq matrix: ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cite_seq/sample_cite_seq.csv``.

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

On your local machine:

	#. Prepare a sample sheet, ``sample_sheet_01.csv``, with the following content::

		OUTNAME,RNA,ADT,TYPE
		exp,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/raw_gene_bc_matrices_h5.h5,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cell_hashing.csv,cell-hashing

	   where **OUTNAME** specifies the subfolder name and prefix of output files, which can be renamed, **RNA** and **ADT** columns specify the RNA and ADT meta-data of samples, and **TYPE** is ``cell-hashing`` for this phase.
	
	#. Prepare an input JSON file, ``input_01.json``, in the same directory as above, with the following content::

		{
			"scCloud_hashing_cite_seq.input_sample_sheet" : "sample_sheet_01.csv",
			"scCloud_hashing_cite_seq.output_directory" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/",
			"scCloud_hashing_cite_seq.demuxEM_min_num_genes" : 500,
			"scCloud_hashing_cite_seq.demuxEM_generate_diagnostic_plots" : true
		}

	   For the details on these options, please refer to `cell-hashing/nuclei-hashing inputs`_.

	   .. _cell-hashing/nuclei-hashing inputs: ../hashing_cite_seq.html#sccloud-hashing-cite-seq-inputs

	#. In the same directory on your local machine, type the following command::

		sccutil fc_run -m scCloud/scCloud_hashing_cite_seq -w ws-lab/ws-01 --bucket-folder my-dir -i input_01.json -o input_01_updated.json

	   Notice that the method here is changed to ``scCloud/scCloud_hashing_cite_seq``, with new JSON file ``input_01.json``.

When the execution is done, you'll get a processed file, ``exp_demux_10x.h5``, stored on cloud ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp/``.


----------------------------------------------------

3. Merge RNA and ADT Matrices for CITE-Seq Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The following steps are done on your local machine:

	#. Prepare a sample sheet, ``sample_sheet_02.csv``, with the following content::

		OUTNAME,RNA,ADT,TYPE
		exp_raw,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp/exp_demux_10x.h5,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/sample_cite_seq.csv,cite-seq

	   The structure of sample sheet here is the same as Phase 2. The difference is that you are now using the output ``h5`` file from Phase 2 as **RNA** here, and the sample **TYPE** is now ``cite-seq``.

	#. Prepare an input JSON file, ``input_02.json``, in the same directory as above, with the following content::

		{
			"scCloud_hashing_cite_seq.input_sample_sheet" : "sample_sheet_02.csv",
			"scCloud_hashing_cite_seq.output_directory" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/",
			"scCloud_hashing_cite_seq.antibody_control_csv" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/citeseq_antibody_control.csv"
		}

	   For the details on these options, please refer to `cell-hashing/nuclei-hashing inputs`_.

	#. In the same directory on your local machine, type the following command::

		sccutil fc_run -m scCloud/scCloud_hashing_cite_seq -w ws-lab/ws-01 --bucket-folder my-dir -i input_02.json -o input_02_updated.json

	   Notice that the input JSON file after ``-i`` option is now ``input_02.json``.

When the execution is done, you'll get a merged raw matrices file, ``exp_raw_merged_10x.h5``, stored on cloud ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp_raw``.


-------------------

4. scCloud Pipeline
^^^^^^^^^^^^^^^^^^^

The following steps are done on your local machine:

	#. Prepare a sample sheet, ``count_matrix_03.csv``, with the following content::

		Sample,Location
		exp,gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/exp_raw/exp_raw_merged_10x.h5

	   This sample sheet describes the metadata for each 10x channel. **Sample** specifies the name for each channel, which can be renamed; **Location** specifies the file location, which is the output of Phase 3.

	#. Prepare an input JSON file, ``input_03.json``, in the same directory as above, with the following content::

		{
			"scCloud.input_count_matrix_csv" : "count_matrix_03.csv",
			"scCloud.output_name" : "gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/results/exp_merged_out",
			"scCloud.num_cpu" : 8,
			"scCloud.select_only_singlets" : true,
			"scCloud.cite_seq" : true,
			"scCloud.run_louvain" : true,
			"scCloud.find_markers_lightgbm" : true,
			"scCloud.remove_ribo" : true,
			"scCloud.mwu" : true,
			"scCloud.annotate_cluster" : true,
			"scCloud.plot_fitsne" : "louvain_labels,assignment",
			"scCloud.plot_citeseq_fitsne" : "louvain_labels,assignment",
			"scCloud.plot_composition" : "louvain_labels:assignment"
		}

	   A typical scCloud pipeline consists of 4 steps, which is given here_. For the details of options above, please refer to `scCloud inputs`_.

	   .. _here: ../scCloud.html#sccloud-steps
	   .. _scCloud inputs: ../scCloud.html#global-inputs

	#. In the same directory on your local machine, type the following command::

		sccutil fc_run -m scCloud/scCloud -w ws-lab/ws-01 --bucket-folder my-dir/results -i input_03.json -o input_03_updated.json

	   Notice that the method is changed to ``scCloud/scCloud``, and the input file after ``-i`` is now ``input_03.json``.

When the execution is done, you'll get the following results stored on cloud ``gs://fc-e0000000-0000-0000-0000-000000000000/my-dir/results/`` to check:
	
	* ``exp_merged_out.h5ad``: The processed RNA matrix data;
	* ``exp_merged_out.de.xlsx``: de_analysis result;
	* ``exp_merged_out.markers.xlsx``: Markers result;
	* ``exp_merged_out.anno.txt``: Annotation output;
	* ``exp_merged_out.fitsne.pdf``: FIt-SNE plot;
	* ``exp_merged_out.citeseq.fitsne.pdf``: CITE-Seq FIt-SNE plot;
	* ``exp_merged_out.louvain_labels+assignment.composition.pdf``: Composition plot.

You can directly go to your Google Bucket to view or download these results.

.. _Terra: https://app.terra.bio/
