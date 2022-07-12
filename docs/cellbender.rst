Run CellBender for ambient RNA removal
------------------------------------------

**cellbender** workflow wraps CellBender_ tool for removing technical artifacts from high-throughput single-cell/single-nucleus RNA sequencing data.

This workflow is modified from the official BSD-3-Clause licensed `CellBender WDL workflow`_, with adding the support on scattering over multiple samples simultaneously.

Prepare input data and import workflow
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

1. Run ``cellranger_workflow``
++++++++++++++++++++++++++++++++

    To demultiplex, you'll need raw gene count and hashtag matrices for cell-hashing/nucleus-hashing data, or raw gene count matrices and genome BAM files for genetic-pooling data. You can generate these data by running the ``cellranger_workflow``.

    Please refer to the `cellranger_workflow tutorial`_ for details.

    When finished, you should be able to find the raw gene count matrix (e.g. ``raw_feature_bc_matrix.h5`` or ``raw_gene_bc_matrices_h5.h5``) for each sample.

2. Import ``cellbender``
++++++++++++++++++++++++++++++

    Import *cellbender* workflow to your workspace by following instructions in `Import workflows to Terra`_. You should choose **github.com/lilab-bcb/cumulus/CellBender** to import.

    Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cellbender* workflow in the drop-down menu.

3. Prepare a sample sheet
++++++++++++++++++++++++++++

    Create a TSV-format sample sheet (say **cellbender_sheet.tsv**), which describes the metadata for each scRNA-Seq sample. Notice that the first column specifies sample names, and the second column specifies the Cloud URI of the **raw** gene-count matrices generated in Step 1.

    An example sample sheet is the follow::

        sample_1	gs://exp/data_1/raw_feature_bc_matrix.h5
        sample_2	gs://exp/data_2/raw_feature_bc_matrix.h5
        sample_3	gs://exp/data_3/raw_feature_bc_matrix.h5

    Then upload your sample sheet to your Terra workspace's bucket using gsutil_. For example::

        gsutil cp /foo/bar/projects/cellbender_sheet.tsv gs://fc-e0000000-0000-0000-0000-000000000000/my-project/

-------------

Workflow inputs
^^^^^^^^^^^^^^^^^^

Below are inputs for *CellBender* workflow. Notice that required inputs are **in bold**:

.. list-table::
    :widths: 5 20 10 5
    :header-rows: 1

    * - Name
      - Description
      - Example
      - Default
    * - **input_tsv_file**
      - Input TSV file describing metadata of scRNA-Seq samples.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/my-project/cellbender_sheet.tsv"
      -
    * - **output_directory**
      - This is the output directory URI for all results. There will be one subfolder per sample under this directory.
      - "gs://fc-e0000000-0000-0000-0000-000000000000/my-project/cellbender_output"
      -
    * - expected_cells
      - Number of cells expected in the dataset (a rough estimate within a factor of 2 is sufficient).
      - 2048
      - ``None``
    * - total_droplets_included
      - The number of droplets from the rank-ordered UMI plot that will be analyzed. The largest *total_droplets_included* droplets will have their cell probabilities inferred as an output.
      - 25000
      - 25000
    * - model
      - Which model is being used for count data:

        - "simple" does not model either ambient RNA or random barcode swapping (for debugging purposes -- not recommended).

        - "ambient" assumes background RNA is incorporated into droplets.

        - "swapping" assumes background RNA comes from random barcode swapping.

        - "full" uses a combined ambient and swapping model.
      - "full"
      - "full"
    * - low_count_threshold
      - Droplets with UMI counts below this number are completely excluded from the analysis. This can help identify the correct prior for empty droplet counts in the rare case where empty counts are extremely high (over 200).
      - 15
      - 15
    * - fpr
      - Target false positive rate in (0, 1). A false positive is a true signal count that is erroneously removed. More background removal is accompanied by more signal removal at high values of FPR. You can specify multiple values by giving a space-separated string, which will create multiple output files.
      - "0.01 0.05 0.1"
      - "0.01"
    * - epochs
      - Number of epochs to train.
      - 150
      - 150
    * - z_dim
      - Dimension of latent variable *z*.
      - 100
      - 100
    * - z_layers
      - Dimension of hidden layers in the encoder for *z*. For multiple layers, specify them in space-separated string format.
      - "500 100 300"
      - "500"
    * - empty_drop_training_fraction
      - Training detail: the fraction of the training data each epoch that is drawn (randomly sampled) from surely empty droplets.
      - 0.5
      - 0.5
    * - blacklist_genes
      - Integer indices of genes to ignore entirely. In the output count matrix, the counts for these genes will be set to zero. For multiple genes, specify them in space-separated string format.
      - "0 1 2"
      - ""
    * - learning_rate
      - Training detail: lower learning rate for inference. A OneCycle learning rate schedule is used, where the upper learning rate is ten times this value. (For this value, probably do not exceed 1e-3).
      - 1e-4
      - 1e-4
    * - exclude_antibody_capture
      - Enalbe this flag will cause remove-background to operate on gene counts only, ignoring other features.
      - false
      - false
    * - docker_registry
      - Docker registry to use.

        - "quay.io/cumulus" for images on Red Hat registry;

        - "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - cellbender_version
      - CellBender version to use. Currently available: ``0.2.0``.
      - "0.2.0"
      - "0.2.0"
    * - zones
      - Google cloud zones to consider for execution. Only works if *backend* is ``gcp``.
      - "us-east1-d us-west1-a us-west1-b"
      - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
    * - num_cpu
      - Number of CPUs used for each sample.
      - 4
      - 4
    * - memory
      - Memory size in string used for each sample.
      - "15G"
      - "15G"
    * - backend
      - Cloud infrastructure backend to use. Available options:

        - "gcp" for Google Cloud;
        - "aws" for Amazon AWS;
        - "local" for local machine.
      - "gcp"
      - "gcp"
    * - gpu_type
      - The GPU type to be used. Only works for ``gcp`` backend. See `here <https://cromwell.readthedocs.io/en/stable/RuntimeAttributes/#gpucount-gputype-and-nvidiadriverversion>`_ for a complete list of available GPU types.
      - "nvidia-tesla-t4"
      - "nvidia-tesla-t4"
    * - disk_space
      - Disk space (integer) in GB needed for each sample.
      - 50
      - 50
    * - preemptible
      - Number of maximum preemptible tries allowed. Only works when *backend* is ``gcp``.
      - 2
      - 2
    * - awsMaxRetries
      - Number of maximum retries when running on AWS. Only works when *backend* is ``aws``.
      - 5
      - 5
    * - awsQueueArn
      - The Arn URI of the AWS job queue to be used. Only works when *backend* is ``aws``.
      - "arn:aws:batch:us-east-1:xxxxxx"
      - ""

----------

Workflow outputs
^^^^^^^^^^^^^^^^^^^

See the table below for *cellbender* workflow outputs:

.. list-table::
    :widths: 5 5 10
    :header-rows: 1

    * - Name
      - Type
      - Description
    * - cellbender_outputs
      - Array[String]
      - A list of Cloud URIs of the output folders. Each folder is associated with one scRNA-seq sample in the given sample sheet.


.. _CellBender: https://github.com/broadinstitute/CellBender
.. _CellBender WDL workflow: https://portal.firecloud.org/#methods/cellbender/remove-background/11/wdl
.. _cellranger_workflow tutorial: ./cellranger/index.html
.. _Import workflows to Terra: ./cumulus_import.html
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
