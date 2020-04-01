Prepare for CITE-Seq analysis
------------------------------

Follow the steps below to run **cumulus** for CITE-Seq data on Terra_.

#. Run Cell Ranger tool to generate raw gene count matrices and antibody hashtag data.

	Please refer to the `cellranger_workflow tutorial`_ for details.

	When finished, you should be able to find the raw gene count matrix (e.g. ``raw_gene_be_matrices_h5.h5``) and ADT count matrix (e.g. ``sample_1_ADT.csv``) for each sample.

#. Create a sample sheet, **sample_sheet_cite_seq.csv**, which describes the metadata for each pair of RNA and antibody hashtag data. The sample sheet should contain 3 columns --- *OUTNAME*, *RNA*, and *ADT*. *OUTNAME* is the output name for one pair of RNA and ADT data. *RNA* and *ADT* are the raw gene count matrix and the ADT count matrix generated in Step 1, respectively.

	There is one optional column *TYPE*. It is kept for backward compatibility with sample sheets working with old versions of Cumulus WDL for CITE-Seq analysis. ``cite-seq`` is the only value accepted for this column.

	Example::

		OUTNAME,RNA,ADT
		sample_1,gs://exp/sample_1/raw_gene_bc_matrices_h5.h5,gs://exp/sample_1_ADT/sample_1_ADT.csv
		sample_2,gs://exp/sample_2/raw_feature_bc_matrices.h5,gs://exp/sample_2_ADT/sample_2_ADT.csv

	Note that in the example above, sample_2 is 10x genomics' v3 chemistry. Cumulus can automatically detect v2/v3 chemistry when loading hdf5 files.

#. (Optional) Create an additional antibody-control sheet **antibody_control.csv** if you have IgG control for each antibody. This sheet contains 2 columns --- *Antibody* and *Control*. 

	Example::

		Antibody,Control
		CD8,Mouse-IgG1
		HLA-ABC,Mouse-IgG2a
		CD45RA,Mouse-IgG2b

#. Upload your sample sheets to the Google bucket of your workspace.  

	Example::
	
		gsutil cp /foo/bar/projects/sample_sheet_cite_seq.csv gs://fc-e0000000-0000-0000-0000-000000000000/
		gsutil cp /foo/bar/projects/antibody_control.csv gs://fc-e0000000-0000-0000-0000-000000000000/

#. Import *cumulus_cite_seq* to your workspace.

	See the Terra documentation for `adding a workflow`_. The *cumulus_hashing_cite_seq* workflow is under ``Broad Methods Repository`` with name "**cumulus/cumulus_cite_seq**".

	Moreover, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace to which you want to export *cumulus_cite_seq* workflow in the drop-down menu.

#. In your workspace, open ``cumulus_cite_seq`` in ``WORKFLOWS`` tab. Select ``Run workflow with inputs defined by file paths`` as below

	.. image:: images/single_workflow.png

   and click the ``SAVE`` button.

---------------------------------

cumulus_cite_seq inputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. list-table::
	:widths: 5 20 10 5
	:header-rows: 1

	* - Name
	  - Description
	  - Example
	  - Default
	* - **input_sample_sheet**
	  - Input CSV file describing metadata of RNA and ADT data pairing.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/sample_sheet_cite_seq.csv"
	  - 
	* - **output_directory**
	  - This is the output directory (gs url + path) for all results. There will be one folder per RNA-ADT data pair under this directory.
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/my_cite_seq_dir"
	  - 
	* - genome
	  - Reference genome name. If not provided, Cumulus will infer the genome name from data.
	  - "GRCh38"
	  - 
	* - antibody_control_csv
	  - Optional parameter. If there is no IgG control information, leave this option blank. Otherwise, specify a CSV file containing the IgG control information for each antibody. 
	  - "gs://fc-e0000000-0000-0000-0000-000000000000/antibody_control.csv"
	  - 
	* - docker_registry
	  - Docker registry to use. Options:

	  	- "cumulusprod" for Docker Hub images; 

	  	- "quay.io/cumulus" for backup images on Red Hat registry.
	  - "cumulusprod"
	  - "cumulusprod"
	* - cumulus_version
	  - cumulus version to use. Versions available: 0.15.0, 0.14.0, 0.13.0, 0.12.0, 0.11.0, 0.10.0.
	  - "0.15.0"
	  - "0.15.0"
	* - zones
	  - Google cloud zones.
	  - "us-east1-d us-west1-a us-west1-b"
	  - "us-central1-a us-central1-b us-central1-c us-central1-f us-east1-b us-east1-c us-east1-d us-west1-a us-west1-b us-west1-c"
	* - memory
	  - Memory size (integer) in GB.
	  - 10
	  - 10
	* - disk_space
	  - Total disk space (integer) in GB.
	  - 20
	  - 20
	* - preemptible
	  - Number of preemptible tries.
	  - 2
	  - 2

---------------------------------

cumulus_cite_seq outputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

See the table below for important *cumulus_cite_seq* outputs:

.. list-table::
	:widths: 5 5 10
	:header-rows: 1

	* - Name
	  - Type
	  - Description
	* - output_folder
	  - String
	  - Google Bucket URL of output directory. Within it, each subfolder is for one RNA-ADT data pair in the input sample sheet.

In the output subfolder of each CITE-Seq RNA-ADT data pair, you can find the following file:

.. list-table::
	:widths: 5 10
	:header-rows: 1

	* - Name
	  - Description
	* - output_name.h5sc
	  - A Cumulus hdf5 format (h5sc) file containing both RNA and ADT count matrices.

---------------------------------

Load CITE-Seq assay into Python and R
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

To load CITE-Seq assay into Python, you need to install Python package `pegasusio <https://pypi.org/project/pegasusio/>`_ first. Then follow the codes below::

	import pegasusio as io
	data = io.read_input("output_name.h5sc")

Once you load the data object, you can get its CITE-Seq count matrix by ``data.get_data('CITE_Seq_<ref-genome>')``, and RNA count matrix by ``data.get_data('<ref-genome>')``, where ``<ref-genome>`` is the reference genome name of the assay.

To load the assay into R, you need to install R package ``reticulate`` in addition to Python package ``pegasusio``. Then follow the codes below::

	library(reticulate)
	ad <- import("pegasusio", convert = FALSE)
	data <- ad$read_input("output_name.h5sc")

And similarly, its CITE-Seq count matrix is achieved by ``data$get_data('CITE_Seq_<ref-genome>')``, and its RNA count matrix by ``data$get_data('<ref-genome>')``, where ``<ref-genome>`` is the reference genome name of the assay.


.. _cellranger_workflow tutorial: ./cellranger.html
.. _gsutil: https://cloud.google.com/storage/docs/gsutil
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/
