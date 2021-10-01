Example of 10X Genomics CellPlex Analysis on Cloud
+++++++++++++++++++++++++++++++++++++++++++++++++++

In this example, you'll learn how to perform Cellplex analysis on Cloud using `Cromwell <https://cromwell.readthedocs.io>`_.

---------------

0. Prerequisite
^^^^^^^^^^^^^^^^^

You need to install the corresponding Cloud SDK tool on your local machine if not:

* `Google Cloud SDK <https://cloud.google.com/sdk/docs/install>`_ for Google Cloud.
* `AWS CLI v2 <https://docs.aws.amazon.com/cli/latest/userguide/install-cliv2.html>`_ for Amazon AWS Cloud.

And then install `Altocumulus <https://github.com/lilab-bcb/altocumulus>`_ in your Python environment. This is the tool for data transfer between local machine and Cloud VM instance.

In this example, we assume that your Cromwell server is already deployed on Cloud at IP address ``10.0.0.0`` with port ``8000``, and also assume using Google Cloud with bucket ``gs://my-bucket``.


1. Extract Genen-Count Matrices
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

First step is to extract gene-count matrices from sequencer output.

In this example, we have the following experiment setting:

* A sample named ``cellplex_gex`` by pooling all RNA data together for sequencing, with index ``SI-TT-A1``;
* A sample named ``cellplex_barcode`` for hashing data, with index ``SI-NN-A1``;
* Three samples to perform individual control, respectively (named ``A``, ``B`` and ``C``), with indexes ``SI-TT-A2``, ``SI-TT-A3``, and ``SI-TT-A4``, respectively.

Besides, you should also have a CellRanger Multi config file like the following::

    [samples],,,
    sample_id,cmo_ids,description,
    CP_1_I,CMO_301,CP_1_I
    CP_1_C,CMO_302,CP_1_C
    CP_2_I,CMO_303,CP_2_I
    CP_2_C,CMO_304,CP_2_C
    CP_3_I,CMO_305,CP_3_I
    CP_3_C,CMO_306,CP_3_C

which means we have 3 donors for the experiment, with hashing IDs shown by ``cmo_ids``.

We need to create a feature barcoding file (say named ``feature_barcode.csv``) for the gene-count matrix extraction on the hashing data. Please refer to `10X Multi CMO Reference`_ for the sequence information of these CMO IDs::

    ATGAGGAATTCCTGC,CP_1_I
    CATGCCAATAGAGCG,CP_1_C
    CCGTCGTCCAAGCAT,CP_2_I
    AACGTTAATCACTCA,CP_2_C
    CGCGATATGGTCGGA,CP_3_I
    AAGATGAGGTCTGTG,CP_3_C

After that, create a sample sheet in CSV format (say named ``cellranger_sample_sheet.csv``) as the following::

    Sample,Reference,Flowcell,DataType,FeatureBarcodeFile
    cellplex_gex,GRCh38-2020-A,/path/to/flowcell/folder,rna
    cellplex_barcode,GRCh38-2020-A,/path/to/flowcell/folder,cmo,/path/to/feature_barcode.csv
    A,GRCh38-2020-A,/path/to/flowcell/folder,rna
    B,GRCh38-2020-A,/path/to/flowcell/folder,rna
    C,GRCh38-2020-A,/path/to/flowcell/folder,rna

where

* ``GRCh38-2020-A`` is the Human GRCh38 (GENCODE v32/Ensembl 98) genome reference prebuilt by Cumulus. See `Cumulus single-cell genome reference list`_ for a complete list of genome references.
* ``/path/to/flowcell/folder`` should be replaced by the local path to the BCL folder of your sequencer output.
* ``/path/to/feature_barcode.csv`` should be replaced by the local path to ``feature_barcode.csv`` file we just created above.
* ``rna`` and ``cmo`` refer to gene expression data and cell multiplexing oligos used in 10X Genomics CellPlex assay, respectively.
* Only the sample of ``cmo`` type needs a feature barcode file for indexing.

For details on preparing this sample sheet, please refer to `CellRanger workflow sample sheet format`_.

Now let's prepare an input JSON file for **cellranger_workflow** WDL workflow to execute (say named ``cellranger_inputs.json``)::

    {
        "cellranger_workflow.input_csv_file": "/path/to/cellranger_sample_sheet.csv",
        "cellranger_workflow.output_directory": "gs://my-bucket/cellplex/cellranger_output"
    }

where

* ``/path/to/cellranger_sample_sheet.csv`` should be replaced by the local path to your sample sheet created above.
* ``gs://my-bucket/cellplex/cellranger_output`` is the target folder on Google bucket to store your result when the workflow job is finished.

For details on these workflow inputs, please refer to `CellRanger workflow inputs`_.

Now we are ready to submit a job to the Cromwell server on Cloud for computing. On your local machine, run the following command::

    alto cromwell run -s 10.0.0.0 -p 8000 -m broadinstitute:cumulus:cellranger:master -i /path/to/cellranger_inputs.json -o cellranger_inputs_updated.json -b gs://my-bucket/cellplex

where

* ``-s`` to specify the server's IP address (or hostname), ``-p`` to specify the server's port number.
* ``-m`` to specify which WDL workflow to use. You should use the Dockstore name of Cumulus `cellranger_workflow`_. Here, the latest version ``master`` is used. If omit the version info, i.e. ``broadinstitute:cumulus:cellranger``, the default version will be used.
* ``-i`` to specify the workflow input JSON file.
* ``-o`` and ``-b`` are used when the input data are local and need to be uploaded to Cloud bucket first. This can be inferred from the workflow input JSON file and sample sheet CSV file.
* ``-o`` to specify the updated workflow input JSON file after uploading the input data, with all the local paths updated to Cloud bucket URLs.
* ``-b`` to specify which folder on Cloud bucket to upload the local input data.

Notice that ``-o`` and ``-b`` options can be dropped if all of your input data are already on Cloud bucket.

After submission, you'll get the job's ID for tracking its status::

    alto cromwell check_status -s 10.0.0.0 -p 8000 --id <your-job-ID>

where ``<your-job-ID>`` should be replaced by the actual Cromwell job ID.

When the job is done, you'll get results in ``gs://my-bucket/cellplex/cellranger_output``. It should contain 6 subfolders, each of which is associated with one sample in ``cellranger_sample_sheet.csv``.

2. Demultiplexing
^^^^^^^^^^^^^^^^^^^

Next, we need to demultiplex the resulting gene-count matrices. In this example, we perform both DemuxEM_ and Souporcell_ methods, respectively.

For **DemuxEM**, we'll need the RNA raw count matrix in HDF5 format (``gs://my-bucket/cellplex/cellranger_output/cellplex_gex/raw_feature_bc_matrix.h5``) and the hashing count matrix in CSV format (``gs://my-buckjet/cellplex/cellranger_output/cellplex_barcode/cellplex_barcode.csv``).

For **Souporcell**, both the RNA raw count matrix above and its corresponding BAM file (``gs://my-bucket/cellplex/cellranger_output/cellplex_gex/possorted_genome_bam.bam``) are needed.

Prepare a sample sheet in CSV format (say named ``demux_sample_sheet.csv``) for demultiplexing, one line for DemuxEM, one for Souporcell::

    OUTNAME,RNA,TagFile,TYPE
    cellplex_demux,gs://my-bucket/cellplex/cellranger_output/cellplex_gex/raw_feature_bc_matrix.h5,gs://my-buckjet/cellplex/cellranger_output/cellplex_barcode/cellplex_barcode.csv,cell-hashing
    cellplex_souporcell,gs://gs://my-bucket/cellplex/cellranger_output/cellplex_gex/raw_feature_bc_matrix.h5,gs://my-bucket/cellplex/cellranger_output/cellplex_gex/possorted_genome_bam.bam,genetic-pooling

where

* ``cell-hashing`` indicates using DemuxEM for demultiplexing, while ``genetic-pooling`` indicates using genetic pooling methods for demultiplexing, with Souporcell being the default.

For details on this sample sheet, please refer to `Demultiplexing workflow sample sheet format`_.

Then prepare a workflow input JSON file (say named ``demux_inputs.json``) for demultiplexing::

    {
        "demultiplexing.input_sample_sheet": "/path/to/demux_sample_sheet.csv",
        "demultiplexing.output_directory": "gs://my-bucket/cellplex/demux_output",
        "demultiplexing.genome": "GRCh38-2020-A",
        "demultiplexing.souporcell_num_clusters": 3
    }

where

* ``/path/to/demux_sample_sheet.csv`` should be replaced by the local path to your ``demux_sample_sheet.csv`` created above.
* ``gs://my-bucket/cellplex/demux_output`` is the Bucket folder to write the results when the job is finished.
* ``GRCh38-2020-A`` is the genome reference used by Souporcell, which should be consistent with your settings in Step 1.
* ``souporcell_num_clusters`` is to set the number of clusters you expect to see for Souporcell clustering. Since we have 3 donors, so set it to 3.

For details, please refer to `Demultiplexing workflow inputs`.

Now submit the demultiplexing job to Cromwell server on Cloud::

    alto cromwell run -s 10.0.0.0 -p 8000 -m broadinstitute:cumulus:demultiplexing:master -i demux_inputs.json -o demux_inputs_updated.json -b gs://my-bucket/cellplex

where

* ``broadinstitute:cumulus:demultiplexing`` refers to demultiplexing_ workflow published on Dockstore.
* We still need ``-o`` and ``-b`` options because ``demux_sample_sheet.csv`` is on the local machine.

Similarly, when the submission succeeds, you'll get another job ID for demultiplexing. You can use it to track the job status.

3. Interactive Data Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^


.. _10X Multi CMO Reference: https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/latest/using/multi#cmoreference
.. _CellRanger workflow sample sheet format: ../cellranger/index.html#prepare-a-sample-sheet
.. _Cumulus single-cell genome reference list: ../cellranger/index.html#sample-sheet
.. _CellRanger workflow inputs: ../cellranger/index.html#workflow-input
.. _cellranger_workflow: https://dockstore.org/workflows/github.com/klarman-cell-observatory/cumulus/Cellranger
.. _DemuxEM: https://demuxem.readthedocs.io
.. _Souporcell: https://github.com/wheaton5/souporcell
.. _Demultiplexing workflow sample sheet format: ../demultiplexing.html#prepare-a-sample-sheet
.. _Demultiplexing workflow inputs: ../demultiplexing.html#workflow-inputs
.. _demultiplexing: https://dockstore.org/workflows/github.com/klarman-cell-observatory/cumulus/Demultiplexing
