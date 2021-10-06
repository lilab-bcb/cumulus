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
* Three samples to perform individual control:

  * Sample ``A`` with index ``SI-TT-A2`` and CMO ID ``CMO_301``,
  * Sample ``B`` with index ``SI-TT-A3`` and CMO ID ``CMO_302``,
  * Sample ``C`` with index ``SI-TT-A4`` and CMO ID ``CMO_303``

To extract feature barcodes for the hashing data, we need to create a feature barcoding file (say named ``feature_barcode.csv``). Please refer to `10X Multi CMO Reference`_ for the sequence information of these CMO IDs::

    ATGAGGAATTCCTGC,A
    CATGCCAATAGAGCG,B
    CCGTCGTCCAAGCAT,C

After that, create a sample sheet in CSV format (say named ``cellranger_sample_sheet.csv``) as the following::

    Sample,Reference,Flowcell,Lane,Index,DataType,FeatureBarcodeFile
    cellplex_gex,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-TT-A1,rna
    cellplex_barcode,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-NN-A1,cmo,/path/to/feature_barcode.csv
    A,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-TT-A2,rna
    B,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-TT-A3,rna
    C,GRCh38-2020-A,/path/to/flowcell/folder,*,SI-TT-A4,rna

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
    cellplex_souporcell,gs://my-bucket/cellplex/cellranger_output/cellplex_gex/raw_feature_bc_matrix.h5,gs://my-bucket/cellplex/cellranger_output/cellplex_gex/possorted_genome_bam.bam,genetic-pooling

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

For details, please refer to `Demultiplexing workflow inputs`_.

Now submit the demultiplexing job to Cromwell server on Cloud::

    alto cromwell run -s 10.0.0.0 -p 8000 -m broadinstitute:cumulus:demultiplexing:master -i demux_inputs.json -o demux_inputs_updated.json -b gs://my-bucket/cellplex

where

* ``broadinstitute:cumulus:demultiplexing`` refers to demultiplexing_ workflow published on Dockstore.
* We still need ``-o`` and ``-b`` options because ``demux_sample_sheet.csv`` is on the local machine.

Similarly, when the submission succeeds, you'll get another job ID for demultiplexing. You can use it to track the job status.

When finished, below are the important output files:

* DemuxEM output: In folder ``gs://my-bucket/cellplex/demux_output/cellplex_demux``,

  * ``cellplex_demux_demux.zarr.zip``: Demultiplexed RNA raw count matrix. This will be used for downstream analysis.
  * ``cellplex_demux.out.demuxEM.zarr.zip``: This file contains intermediate results for both RNA and hashing count matrices, which is useful for compare with other demultiplexing methods.
  * DemuxEM plots in PDF format. They are used for estimating the performance of DemuxEM on the data.

* Souporcell output: In folder ``gs://my-bucket/cellplex/demux_output/cellplex_souporcell``,

  * ``cellplex_souporcell_demux.zarr.zip``: Demultiplexed RNA raw count matrix. This will be used for downstream analysis.
  * ``clusters.tsv``: Inferred droplet type and cluster assignment for each cell barcode.
  * ``cluster_genotypes.vcf``: Inferred genotypes for each cluster.

3. Interactive Data Analysis
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

You may use `Cumulus workflow`_ to perform the downstream analysis in a batch way.
Alternatively, you can also download the demultiplexing results from the Cloud bucket to your local machine, and perform the analysis interactively.
This section introduces how to use Cumulus' analysis module Pegasus to load demultiplexing results, perform quality control (QC), and compare the performance of the two methods.

You'll need to first install `Pegasus`_ in your local Python environment. Also, download the demultiplexed raw counts in ``.zarr.zip`` format mentioned above to your local machine.

3.1. Extract Singlet/Doublet Type and Assignment
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can load the DemuxEM result, and perform QC by::

    import pegasus as pg
    data_demuxEM = pg.read_input("cellplex_demux_demux.zarr.zip")
    pg.qc_metrics(data_demuxEM, min_genes=500, max_genes=6000, mito_prefix='MT-', percent_mito=20)
    pg.filter_data(data_demuxEM)

where ``qc_metrics`` and ``filter_data`` are Pegasus functions to filter out low quality cells, and keep those with number of genes within range ``[500, 6000)``
and having expression of mitochondrial genes ``<= 20%``. Please see `Pegasus preprocess tools`_ for details.

There are two columns in `data_demuxEM.obs` field related to demultiplexing results:

* **demux_type**: This column stores the singlet/doublet type of each cell: ``singlet``, ``doublet``, or ``unknown``.
* **assignment**: This column stores the more detailed assignment of cells regarding samples/donors.

To get the distribution regarding these columns, e.g. *demux_type*::

    data_demuxEM.obs['demux_type'].value_counts()

Besides, you can export the cell barcodes along with their singlet/doublet type and assignment as a CSV file by::

    data_demuxEM.obs[['demux_type', 'assignment']].to_csv("demuxEM_assignment.csv")

We can also do it similarly for the Souporcell result as above, by reading ``cellplex_souporcell_demux.zarr.zip`` instead.

3.2. Compare the Two Demultiplexing Methods
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

We can compare the performance of DemuxEM and Souporcell by plotting a heatmap showing their singlet/doublet assignment results.

Assume we've already loaded the two results (``data_demuxEM`` for DemuxEM result, ``data_souporcell`` for Souporcell result), and performed QC as in 3.1.
The following Python code will generate this heatmap in an interactive Python environment (e.g. in a Jupyter notebook)::

    import pandas as pd
    import matplotlib.pyplot as plt
    import seaborn as sns

    def extract_assignment(data):
        assign = data.obs['demux_type'].values.astype('object')
        idx_singlet = (data.obs['demux_type'] == 'singlet').values
        assign[idx_singlet] = data.obs.loc[idx_singlet, 'assignment'].values.astype(object)
        return assign

    assign_demuxEM = extract_assignment(data_demuxEM)
    assign_souporcell = extract_assignment(data_souporcell)

    df = pd.crosstab(assign_demuxEM, assign_souporcell)
    df.columns.name = df.index.name = ""
    ax = plt.gca()
    ax.xaxis.tick_top()
    ax = sns.heatmap(df, annot=True, fmt='d', cmap='inferno', ax=ax)
    plt.tight_layout()
    plt.gcf().dpi=500

3.3. Downstream Analysis
~~~~~~~~~~~~~~~~~~~~~~~~~

To perform further downstream analysis on the singlets, please refer to `Pegasus tutorials`_.


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
.. _Cumulus workflow: ../cumulus.html
.. _Pegasus: https://pegasus.readthedocs.io/en/stable/installation.html
.. _Pegasus preprocess tools: https://pegasus.readthedocs.io/en/stable/api/index.html#preprocess
.. _Pegasus tutorials: https://pegasus.readthedocs.io/en/stable/tutorials.html
