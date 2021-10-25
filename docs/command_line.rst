Run Terra pipelines via command line
----------------------------------------------

You can run Terra pipelines via the command line by installing the **altocumulus** package.

Install ``altocumulus`` for Broad users
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Request an UGER node::

    reuse UGER
    qrsh -q interactive -l h_vmem=4g -pe smp 8 -binding linear:8 -P regevlab

The above command requests an interactive shell using the regevlab project with 4G memory per thread, 8 threads. Feel free to change the memory, thread, and project parameters.

Add conda to your path::

    reuse Anaconda3

Activate the alto virtual environment::

    source activate /seq/regev_genome_portal/conda_env/cumulus

Install ``altocumulus`` for non-Broad users
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

#. Make sure you have ``conda`` installed. If you haven't installed conda_, use the following commands to install it on Linux::

    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh .
    bash Miniconda3-latest-Linux-x86_64.sh -p /home/foo/miniconda3
    mv Miniconda3-latest-Linux-x86_64.sh /home/foo/miniconda3

   where ``/home/foo/miniconda3`` should be replaced by your own folder holding Miniconda3.

Or use the following commdands for MacOS installation::

    curl -O curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-x86_64.sh
    bash Miniconda3-latest-MacOSX-x86_64.sh -p /Users/foo/miniconda3
    mv Miniconda3-latest-MacOSX-x86_64.sh /Users/foo/miniconda3

    where ``/Users/foo/miniconda3`` should be replaced by your own folder holding Miniconda3.

#. Create a conda environment named "alto" and install ``altocumulus``::

    conda create -n alto -y pip
    source activate alto
    pip install altocumulus

When the installation is done, type ``alto -h`` in terminal to see if you can see the help information.

Set up Google Cloud Account
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Install `Google Cloud SDK <https://cloud.google.com/cloud-sdk>`_ on your local machine.

Then type the following command in your terminal

.. code-block:: bash

    gcloud auth application-default login

and follow the pop-up instructions to set up your Google cloud account.

Run Terra workflows
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**alto run** runs a Terra method. Features:

- Uploads local files/directories in your inputs to a Google Cloud bucket updates the file paths to point to the Google Cloud bucket.

   Your sample sheet can point to local file paths. In this case, ``alto run`` will take care of uploading directories smartly (e.g. only upload necessary files in BCL folders) and modifying the sample sheet to point to a Google Cloud bucket.

- Creates or uses an existing workspace.

- Uses the latest version of a method unless the method version is specified.

**alto cromwell run** runs on a chosen backend. Features:

- Uploads local files/directories in your inputs to an appropriate location depending on backend chosen and updates the file paths to point to the bucket information.

- Uses the method parameter to pull in appropriate worflow to import and run.   

Options
+++++++

Required options are in bold.

.. list-table::
    :widths: 5 20
    :header-rows: 1

    * - Name
      - Description
    * - | **-m <METHOD>**
        | **\\-\\-method <METHOD>**
      - | Specify a Terra workflow *<METHOD>* to use.
        | *<METHOD>* is of format *Namespace/Name* (e.g. ``cumulus/cellranger_workflow``).
        | A snapshot version number can optionally be specified (e.g. ``cumulus/cellranger_workflow/4``); otherwise the latest snapshot of the method is used.
    * - | **-w <WORKSPACE>**
        | **\\-\\-workspace <WORKSPACE>**
      - | Specify which Terra workspace *<WORKSPACE>* to use.
        | *<WORKSPACE>* is also of format *Namespace/Name* (e.g. ``foo/bar``). The workspace will be created if it does not exist.
    * - | **-i <WDL_INPUTS>**
        | **\\-\\-inputs <WDL_INPUTS>**
      - | Specify the WDL input JSON file to use.
        | It can be a local file, a JSON string, or a Google bucket URL directing to a remote JSON file.
    * - | \\-\\-bucket-folder <folder>
      - | Store inputs to <folder> under workspace's google bucket.
    * - | -o <updated_json>
        | \\-\\-upload <updated_json>
      - | Upload files/directories to Google bucket of the workspace, and generate an updated input JSON file (with local paths replaced by Google bucket URLs) to <updated_json> on local machine.
    * - | \\-\\-no-cache
      - | Disable Terra cache calling


#. Example run on Terra
+++++++++++++++++++++++++

This example shows how to use ``alto run`` to run cellranger_workflow to extract gene-count matrices from sequencing output.

#. Prepare your sample sheet ``example_sample_sheet.csv`` as the following::

    Sample,Reference,Flowcell,Lane,Index,Chemistry
    sample_1,GRCh38,/my-local-path/flowcell1,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,/my-local-path/flowcell1,3-4,SI-GA-B8,threeprime
    sample_3,mm10,/my-local-path/flowcell1,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,/my-local-path/flowcell1,7-8,SI-GA-D8,fiveprime
    sample_1,GRCh38,/my-local-path/flowcell2,1-2,SI-GA-A8,threeprime
    sample_2,GRCh38,/my-local-path/flowcell2,3-4,SI-GA-B8,threeprime
    sample_3,mm10,/my-local-path/flowcell2,5-6,SI-GA-C8,fiveprime
    sample_4,mm10,/my-local-path/flowcell2,7-8,SI-GA-D8,fiveprime

   where ``/my-local-path`` is the top-level directory of your BCL files on your local machine.

   Note that ``sample_1``, ``sample_2``, ``sample_3``, and ``sample_4`` are sequenced on 2 flowcells.


#. Prepare your JSON input file ``inputs.json`` for cellranger_workflow::

    {
        "cellranger_workflow.input_csv_file" : "/my-local-path/sample_sheet.csv",
        "cellranger_workflow.output_directory" : "gs://url/outputs",
        "cellranger_workflow.delete_input_bcl_directory": true
    }

   where ``gs://url/outputs`` is the folder on Google bucket of your workspace to hold output.

#. Run the following command to kick off your Terra workflow::

    alto run -m cumulus/cellranger_workflow -i inputs.json -w myworkspace_namespace/myworkspace_name -o inputs_updated.json

   where ``myworkspace_namespace/myworkspace_name`` should be replaced by your workspace namespace and name.


Upon success, ``alto run`` returns a URL pointing to the submitted Terra job for you to monitor.

If for any reason, your job failed. You could rerun it without uploading files again via the following command::

    alto run -m cumulus/cellranger_workflow -i inputs_updated.json -w myworkspace_namespace/myworkspace_name

because ``inputs_updated.json`` is the updated version of ``inputs.json`` with all local paths being replaced by their corresponding Google bucket URLs after uploading.

#. Example import of any Cumulus workflow 
++++++++++++++++++++++++++++++++++++++++++

This example shows how to use ``alto cromwell run`` to run demultiplexing workflow on any backend.

#. Prepare your sample sheet ``demux_sample_sheet.csv`` as the following::

     OUTNAME,RNA,TagFile,TYPE
     sample_1,gs://exp/data_1/raw_feature_bc_matrix.h5,gs://exp/data_1/sample_1_ADT.csv,cell-hashing
     sample_2,gs://exp/data_2/raw_feature_bc_matrix.h5,gs://exp/data_3/possorted_genome_bam.bam,genetic-pooling

#. Prepare your JSON input file ``cumulus_inputs.json`` for cellranger_workflow::

     {
        "demultiplexing.input_sample_sheet" : "demux_sample_sheet.csv",
        "demultiplexing.output_directory" : "gs://url/outputs",
        "demultiplexing.zones" : "us-west1-a us-west1-b us-west1-c",
        "demultiplexing.genome" : "GRCh38-2020-A"
     }

   where ``gs://url/outputs`` is the folder on Google bucket of your workspace to hold output.

#. Run the following command to kick off your run on a chosen backend::

    alto cromwell run -s 10.10.10.10 -p 3000 -m broadinstitute:cumulus:Demultiplexing:master \
                      -i cumulus_inputs.json --no-ssl-verify

.. _conda: https://docs.conda.io/en/latest/miniconda.html
