Version 1.5.1 :small:`September 15, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Fix the issue of WDLs after Terra platform updates the Cromwell engine.

Version 1.5.0 :small:`July 20, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* On *demultiplexing* workflow
    * Update *demuxEM* to v0.1.6.
* On *cumulus* workflow
    * Add Nonnegative Matrix Factorization (NMF) feature: ``run_nmf`` and ``nmf_n`` inputs.
    * Add integrative NMF (iNMF) data integration method: ``inmf`` option in ``correction_method`` input; the number of expected factors is also specified by ``nmf_n`` input.
    * When NMF or iNMF is enabled, word cloud plots and gene program UMAP plots of NMF/iNMF results will be generated.
    * Update *Pegasus* to v1.4.2.

Version 1.4.0 :small:`May 17, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* On *cellranger* workflow
    * Add support for multiomics analysis using linked samples, *cellranger-arc count*, *cellranger multi* and *cellranger count* will be automatically triggered based on the sample sheet
    * Add support for cellranger version 6.0.1 and 6.0.0
    * Add support for cellranger-arc version 2.0.0, 1.0.1, 1.0.0
    * Add support for cellranger-atac version 2.0.0
    * Add support for cumulus_feature_barcoding version 0.6.0, which handles CellPlex CMO tags
    * Add *GRCh38-2020-A_arc_v2.0.0*, *mm10-2020-A_arc_v2.0.0*, *GRCh38-2020-A_arc_v1.0.0* and *mm10-2020-A_arc_v1.0.0* references for *cellranger-arc*.
    * Fixed bugs in cellranger_atac_create_reference
    * Add delete undetermined FASTQs option for mkfastq
* On *demultiplexing* workflow
    * Replace *demuxlet* with *popscle*, which includes both *demuxlet* and *freemuxlet*
* On *cumulus* workflow
    * Fixed bug that ``remap_singlets`` and ``subset_singlets`` don't work when input is in sample sheet format.
* Modified workflows to remove trailing spaces and support spaces within output_directory

Version 1.3.0 :small:`February 2, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* On *cumulus* workflow:
    * Change ``cumulus_version`` to ``pegasus_version`` to avoid confusion.
    * Update to use Pegasus v1.3.0 for analysis.

Version 1.2.0 :small:`January 19, 2021`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Add *spaceranger* workflow:
    * Wrap up spaceranger version 1.2.1
* On *cellranger* workflow:
    * Fix workflow WDL to support both single index and dual index
    * Add support for cellranger version 5.0.1 and 5.0.0
    * Add support for targeted gene expression analysis
    * Add support for ``--include-introns`` and ``--no-bam`` options for cellranger count
    * Remove ``--force-cells`` option for cellranger vdj as noted in cellranger 5.0.0 release note
    * Add *GRCh38_vdj_v5.0.0* and *GRCm38_vdj_v5.0.0* references
* Bug fix on *cumulus* workflow.
* Reorganize the sidebar of Cumulus documentation website.

Version 1.1.0 :small:`December 28, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* On *cumulus* workflow:
    * Add CITE-Seq data analysis back. (See section `Run CITE-Seq analysis <./cumulus.html#run-cite-seq-analysis>`_ for details)
    * Add doublet detection. (See ``infer_doublets``, ``expected_doublet_rate``, and ``doublet_cluster_attribute`` input fields)
    * For tSNE visualization, only support FIt-SNE algorithm. (see ``run_tsne`` and ``plot_tsne`` input fields)
    * Improve efficiency on log-normalization and DE tests.
    * Support multiple marker JSON files used in cell type annotation. (see ``organism`` input field)
    * More preset gene sets provided in gene score calculation. (see ``calc_signature_scores`` input field)
* Add *star_solo* workflow (see `STARsolo section <./starsolo.html>`_ for details):
    * Use `STARsolo <https://github.com/alexdobin/STAR/blob/master/docs/STARsolo.md>`_ to generate count matrices from FASTQ files.
    * Support chemistry protocols such as 10X-V3, 10X-V2, DropSeq, and SeqWell.
* Update the example of analyzing hashing and CITE-Seq data (see `Example section <./examples/example_hashing_citeseq.html>`_) with the new workflows.
* Bug fix.

Version 1.0.0 :small:`September 23, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Add *demultiplexing* workflow for cell-hashing/nucleus-hashing/genetic-pooling analysis.
* Add support on CellRanger version ``4.0.0``.
* Update *cumulus* workflow with Pegasus version ``1.0.0``:
    * Use ``zarr`` file format to handle data, which has a better I/O performance in general.
    * Support focus analysis on Unimodal data, and appending other Unimodal data to it. (``focus`` and ``append`` inputs in *cluster* step).
    * Quality-Control: Change ``percent_mito`` default from ``10.0`` to ``20.0``; by default remove bounds on UMIs (``min_umis`` and ``max_umis`` inputs in *cluster* step).
    * Quality-Control: Automatically figure out name prefix of mitochondrial genes for ``GRCh38`` and ``mm10`` genome reference data.
    * Support signature / gene module score calculation. (``calc_signature_scores`` input in *cluster* step)
    * Add *Scanorama* method to batch correction. (``correction_method`` input in *cluster* step).
    * Calculate UMAP embedding by default, instead of FIt-SNE.
    * Differential Expression (DE) analysis: remove inputs ``mwu`` and ``auc`` as they are calculated by default. And cell-type annotation uses MWU test result by default.
* Remove *cumulus_subcluster* workflow.
