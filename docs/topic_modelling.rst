Topic Modelling
--------------------

Prepare input data
~~~~~~~~~~~~~~~~~~~


Follow the steps below to run **topic_modelling** on Terra_.

1. Prepare your count matrix. **Cumulus** currently supports the following formats:
- anndata (h5ad);
- 10x genomics v2/v3 format (hdf5);
- Drop-seq dge format;
- csv (no HCA DCP format), tsv or loom formats.


2. Upload your count matrix to the workspace.

    Example::
    
        gsutil cp /foo/bar/projects/dataset.h5ad gs://fc-e0000000-0000-0000-0000-000000000000/

    where ``/foo/bar/projects/dataset.h5ad`` is the path to your dataset on your local machine, and
    ``gs://fc-e0000000-0000-0000-0000-000000000000/`` is the Google bucket destination.

3. Import *topic_modelling* workflow to your workspace.

    See the Terra documentation for `adding a workflow`_. The *cumulus* workflow is under ``Broad Methods Repository`` with name "**cumulus/topic_modelling**".

    Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *topic_modelling* workflow in the drop-down menu.

4. In your workspace, open ``topic_modelling`` in ``WORKFLOWS`` tab. Select ``Run workflow with inputs defined by file paths`` as below

    .. image:: images/single_workflow.png

   and click the ``SAVE`` button.


Workflow input
~~~~~~~~~~~~~~~~

Below are inputs for *topic_modelling* workflow. Required inputs are in bold.

.. list-table::
    :widths: 5 20 10 5
    :header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - **input_file**
      - Google bucket URL of the input count matrix.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/my_dataset.h5ad"
      -
    * - **number_of_topics**
      - Array of number of topics.
      - [10,15,20]
      -
    * - prefix_exclude
      - Comma separated list of features to exclude that start with prefix.
      - "mt-,Rpl,Rps"
      - "mt-,Rpl,Rps"
    * - min_percent_expressed
      - Exclude features expressed below min_percent.
      - 2
      -
    * - max_percent_expressed
      - Exclude features expressed below min_percent.
      - 98
      -
    * - random_number_seed
      - Random number seed for reproducibility.
      - 0
      - 0

----------------------

Workflow output
~~~~~~~~~~~~~~~~

.. list-table::
    :widths: 5 5 20
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - coherence_plot
      - File
      - Plot of coherence scores vs. number of topics
    * - perplexity_plot
      - File
      - Plot of perplexity values vs. number of topics
    * - cell_scores
      - Array[File]
      - Topic by cells (one file for each topic number)
    * - feature_topics
      - Array[File]
      - Topic by features (one file for each topic number)
    * - report
      - Array[File]
      - HTML visualization report (one file for each topic number)
    * - stats
      - Array[File]
      - Computed coherence and perplexity (one file for each topic number)
    * - model
      - Array[File]
      - Serialized LDA model (one file for each topic number)
    * - corpus
      - File
      - Serialized corpus
    * - dictionary
      - File
      - Serialized dictionary
.. _Terra: https://app.terra.bio/
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository






