Demuxlet
-------------------------------------------------------------

This workflow runs demuxlet_ to deconvolute sample identity when multiple samples are pooled by barcoded single-cell sequencing.

#. Align your single-cell sequencing data (for example using the cellranger_ or drop_seq_ workflows).

#. Create a sample sheet.

	Please note that the columns in the tab separated file must be in the order shown below and does not contain a header line.

	.. list-table::
		:widths: 5 30
		:header-rows: 1

		* - Column
		  - Description
		* - Name
		  - Sample name.
		* - BAM
		  - Location of the BAM file in the cloud (gs:// URL).
		* - Barcodes
		  - Location of the valid cellular barcodes file in the cloud (gs:// URL).
		* - VCF
		  - Location of the VCF file to use for this sample in the cloud (gs:// URL).

	Example::

		sample-1,gs://fc-e0000000/sample-1/out/possorted_genome_bam.bam,gs://fc-e0000000/sample-1/out/filtered_feature_bc_matrix/barcodes.tsv.gz,gs://fc-e0000000/sample-1.vcf
		sample-2,gs://fc-e0000000/sample-2/out/possorted_genome_bam.bam,gs://fc-e0000000/sample-2/out/filtered_feature_bc_matrix/barcodes.tsv.gz,gs://fc-e0000000/sample-2.vcf




#. Upload your sample sheet to the workspace bucket.

	Example::

		gsutil cp /foo/bar/projects/sample_sheet.tsv gs://fc-e0000000/


#. Import *demuxlet* workflow to your workspace.

	See the Terra documentation for `adding a workflow`_. The workflow is under ``Broad Methods Repository`` with the name "**cumulus/demuxlet**".

	Next, in the workflow page, click the ``Export to Workspace...`` button, and select the workspace you want to export to in the drop-down menu.

#. In your workspace, open ``demuxlet`` in ``WORKFLOWS`` tab. Select ``Process single workflow from files`` as below

	.. image:: images/single_workflow.png

   and click the ``Save`` button.

---------------------------------

Inputs
^^^^^^^

Please see the description of important inputs below.

.. list-table::
    :widths: 5 30
    :header-rows: 1

    * - Column
      - Description
    * - tsv_file
      - Four column tab-separated file without a header with name, coordinate sorted bam, barcodes, and vcf
    * - min_MQ
      - Minimum mapping quality to consider (default 20)
    * - alpha
      - Grid of alpha to search for (default [0.1, 0.2, 0.3, 0.4, 0.5]).
    * - min_TD
      - Minimum distance to the tail (default 0)
    * - tag_group
      - Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups (default "CB")
    * - tag_UMI
      - Tag representing UMIs (default "UB"")
    * - field
      - FORMAT field to extract the genotype, likelihood, or posterior from (default "GT")
    * - geno_error
      - Offset of genotype error rate (default 0.1)

Outputs
^^^^^^^^

The demuxlet output file contains the best guess of the sample identity, with detailed statistics to reach to the best guess.

.. _demuxlet: https://github.com/statgen/popscle
.. _adding a workflow: https://support.terra.bio/hc/en-us/articles/360025674392-Finding-the-tool-method-you-need-in-the-Methods-Repository
.. _Terra: https://app.terra.bio/
.. _cellranger: cellranger.html
.. _drop_seq: drop_seq.html
