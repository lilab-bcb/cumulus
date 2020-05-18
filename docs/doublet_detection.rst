Doublet Detection
--------------------

Prepare input data
~~~~~~~~~~~~~~~~~~~

Case One: Sample Sheet
^^^^^^^^^^^^^^^^^^^^^^^^

Follow the steps below to run **doublet_detection** on Terra_.

1. Create a sample sheet, **count_matrix.csv**, which describes the metadata for each sample count matrix. The sample sheet should at least contain 2 columns --- *Sample* and *Location*. *Sample* refers to sample names and *Location* refers to the location of the channel-specific count matrix in either of

  - 10x format with v2 chemistry. For example, ``gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_gene_bc_matrices_h5.h5``.
  - 10x format with v3 chemistry. For example, ``gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/filtered_feature_bc_matrices.h5``.

Example::

	Sample,Location
	sample_1,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/raw_gene_bc_matrices_h5.h5
	sample_2,gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_2/raw_feature_bc_matrices.h5

If you ran **cellranger_workflow** ahead, you should already obtain a template **count_matrix.csv** file that you can modify from **generate_count_config**'s outputs.

2. Upload your sample sheet to the workspace.  

    Example::
    
        gsutil cp /foo/bar/projects/sample_sheet.csv gs://fc-e0000000-0000-0000-0000-000000000000/

    where ``/foo/bar/projects/sample_sheet.csv`` is the path to your sample sheet in local machine, and ``gs://fc-e0000000-0000-0000-0000-000000000000/`` is the location on Google bucket to hold it.

3. Import *doublet_detection* workflow to your workspace.

    See the Terra documentation for `adding a workflow`_. The *cumulus* workflow is under ``Broad Methods Repository`` with name "**cumulus/doublet_detection**".

    Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *doublet_detection* workflow in the drop-down menu.

4. In your workspace, open ``doublet_detection`` in ``WORKFLOWS`` tab. Select ``Run workflow with inputs defined by file paths`` as below

    .. image:: images/single_workflow.png

   and click the ``SAVE`` button.

Case Two: Terra Data Table
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Alternatively, instead of sample sheet, you can also set up sample information in `Terra data table`_. In brief,

1. You need to first `create a data table`_ on Terra. An example TSV file is the following::

    entity:test_sample_id  input_h5
    sample_1  gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_1/raw_gene_bc_matrix_h5.h5
    sample_2  gs://fc-e0000000-0000-0000-0000-000000000000/my_dir/sample_2/raw_feature_bc_matrix.h5

You are free to add more columns, but sample ids and URLs to RNA count matrix files are required. I'll use this example TSV file for the rest of steps in this case.

2. Upload your TSV file to your workspace. Open the ``DATA`` tab on your workspace. Then click the upload button on left ``TABLE`` panel, and select the TSV file above. When uploading is done, you'll see a new data table with name "test_sample".

3. Import *doublet_detection* workflow to your workspace as in Case one. Then open ``doublet_detection`` in ``WORKFLOW`` tab. Select ``Run workflow(s) with inputs defined by data table``, and choose *test* from the drop-down menu.

4. In the input field, specify:

  - ``input_h5_file``: Type ``this.input_h5``, where ``this`` refers to the data table selected, and ``input_h5`` is the column name in this data table for RNA count matrices.
  - ``output_directory``: Type Google bucket URL for the main output folder. For example, ``gs://fc-e0000000-0000-0000-0000-000000000000/doublet_detection_results``.
  - ``sample_id``: Type ``this.test_sample_id``, where ``test_sample_id`` is the column name in data table for sample IDs.

Notice that you should make ``input_sample_sheet`` field empty; otherwise, your job will ignore the above settings, and only use the sample sheet as input data. 

Finish setting up other inputs following the description in sections below. When you are done, click ``SAVE``, and then ``RUN ANALYSIS``.

When all the jobs are done, you'll find output for the 2 samples in subfolders ``gs://fc-e0000000-0000-0000-0000-000000000000/doublet_detection_results/sample_1`` and ``gs://fc-e0000000-0000-0000-0000-000000000000/doublet_detection_results/sample_2``, respectively.

Case Three: Single File
^^^^^^^^^^^^^^^^^^^^^^^^^

If you have only one sample for doublet detection, you can directly use its URL to RNA count matrix file without sample sheet or data table.

Import *doublet_detection* workflow to your workspace as described in cases above. In workflow page, select ``Run workflow with inputs defined by file paths``. Then in the input field, specify:

  - ``input_h5_file``: Click the folder button to select the sample's hdf5 file.
  - ``sample_id``: Give a sample name for the result files.

Notice that you should make ``input_sample_sheet`` field empty; otherwise, your job will ignore the above settings, and only use the sample sheet as input data.

Workflow input
~~~~~~~~~~~~~~~~

Below are inputs for *doublet_detection* workflow. Notice that required inputs are in bold.

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **output_directory**
	  - Google bucket URL of the output directory.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_results_dir"
	  -
	* - input_sample_sheet
	  - Input CSV sample sheet describing metadata of each 10x channel.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet.csv"
	  -
	* - sample_id
	  - | This is the name of subdirectory for the current sample; and all output files within the subdirectory will have this string as the common filename prefix.
	    | Notice that this input will be ignored if ``input_sample_sheet`` is set.
	  - "sample_1"
	  -
	* - input_file
	  - Input count matrix file. Notice that this input will be ignored if ``input_sample_sheet`` is set.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_1/raw_feature_bc_matrix.h5"
	  -
	* - select_singlets
	  - Whether select singlets only or not.
	  - true
	  - true
	* - mito_prefix
	  - Prefix of mitochondrial gene names. This is to identify mitochondrial genes.
	  - "mt-"
	  - "MT-"
	* - min_genes
	  - Only keep cells with at least <min_genes> of genes.
	  - 500
	  - 500
	* - max_genes
	  - Only keep cells with less than <max_genes> of genes.
	  - 6000
	  - 6000
	* - min_umis
	  - Only keep cells with at least <min_umis> of UMIs.
	  - 100
	  - 100
	* - max_umis
	  - Only keep cells with less than <max_umis> of UMIs.
	  - 600000
	  - 600000
	* - percent_mito
	  - Only keep cells with mitochondrial ratio less than <percent_mito>% of total counts.
	  - 50
	  - 10.0
	* - gene_percent_cells
	  - Only use genes that are expressed in at <gene_percent_cells>% of cells to select variable genes.
	  - 50
	  - 0.05
	* - expected_doublet_rate
	  - The expected doublet rate in the experiment.
	  - 0.1
	  - 0.1
	* - random_state
	  - Random state for doublet simulation, approximate nearest neighbor search, and PCA/TruncatedSVD.
	  - 0
	  - 0
	* - nPC
	  - Number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction.
	  - 30
	  - 30
	* - docker_registry
	  - Docker registry to use. Options:

	    - "cumulusprod" for Docker Hub images; 

	    - "quay.io/cumulus" for backup images on Red Hat registry.
	  - "cumulusprod"
	  - "cumulusprod"
	* - config_version
	  - Version of config docker image for processing sample sheet. Currently only 0.1 is available.
	  - "0.1"
	  - "0.1"
	* - scrublet_version
	  - Scrublet version for doublet detection. Currently only 0.2.1 is available.
	  - "0.2.1"
	  - "0.2.1"
	* - zones
	  - Google cloud zones to consider for execution.
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - memory
	  - Memory size in GB needed per sample.
	  - 10
	  - 10
	* - disk_space
	  - Disk space in GB per sample.
	  - 10
	  - 10
	* - preemptible
	  - Number of maximum preemptible tries allowed.
	  - 2
	  - 2

----------------------

Workflow output
~~~~~~~~~~~~~~~~

.. list-table::
    :widths: 5 5 20
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - output_zarr
      - File
      - The count matrix in ``zarr`` format, with doublet scores for cells added. This output is for the cases of data table and single file input data.
    * - output_histogram_pdf
      - File
      - Histogram of doublet scores for observed transcriptomes and simulated doublets. This output is for the cases of data table and single file input data.
    * - output_zarr_list
      - Array[File]
      - A list of count matrix files in ``zarr`` format, each of which has doublet scores for cells added per sample. This output is for the case of sample sheet input data.
    * - output_histogram_pdf_list
      - Array[File]
      - A list of histograms of doublet scores for observed transcriptomes and simulated doublets, each of which is associated with one sample. This output is for the case of sample sheet input data.


.. _Terra: https://app.terra.bio/
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra data table: https://support.terra.bio/hc/en-us/articles/360025758392-Managing-data-with-tables-
.. _create a data table: https://support.terra.bio/hc/en-us/articles/360025758392