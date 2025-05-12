3.1.0 :small:`May 12, 2025`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * Upgrade *cumulus_feature_barcoding_version* default to ``2.0.0``.

        * Apply UMI correction to raw count matrices; for ``crispr`` type samples, further apply PCR chimeric filtering.
    * Rename ``min_read_ratio`` to ``read_ratio_cutoff``, and default is changed from ``0.1`` to ``0.5``. Still only applied to ``crispr`` type samples.
    * Samples of ``crispr``, ``hashing``, ``citeseq``, ``adt`` and ``cmo`` types now accept *Reference* column values.

3.0.0 :small:`March 7, 2025`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Overall highlights:

    * For data localization/delocalization on GCP, now use *gcloud storage* commands instead of *gsutil*, as *gsutil* is deprecated, and *gcloud storage* achieves some speed improvement.
    * Users no longer need to specify ``backend`` input for AWS or local backend, and it is now automatically figured out from ``output_directory`` location.
    * Remove support on *mkfastq* in Cellranger and Spaceranger workflows, as it will soon be removed from Cell Ranger and Space Ranger. Users need to run `BCL Convert`_ themselves first to generate FASTQ data.
    * For Cellranger and Spaceranger workflows, provide better support for shared computing environments like AWS Batch, GCP Batch and sHPC.
    * For resources like prebuilt references, they are now held in a single-region bucket ``gs://cumulus-ref`` in **US-CENTRAL1** to reduce potentially higher network cost at users' side from the previous multi-region bucket.

* Cellranger workflow:

    * Upgrade *cellranger_version* default to ``9.0.1``.
    * Remove mkfastq related workflow inputs.
    * Remove *run_count* input as it's always ``true``.
    * Support input data format as a TAR file containing FASTQ files.
    * Change *FeatureBarcodeFile* column header to **AuxFile** (for backward compatibility, *FeatureBarcodeFile* is still accepted but not recommended).
    * Support all the 3 `Sample Multiplexing`_ methods provided since Cell Ranger v9.0. See `details <./cellranger/index.html#flex-sample-multiplexing-and-multiomics>`_
    * For single-cell and single-nucleus RNA-seq:

        * Add new genome reference ``mRatBN7.2-2024-A``
        * Remove Target Gene Expression-related inputs, as it's no longer supported since Cell Ranger v7.2.0.
    * For feature barcoding:

        * Upgrade *cumulus_feature_barcoding_version* to `1.0.0 <https://github.com/lilab-bcb/cumulus_feature_barcoding/releases/tag/1.0.0>`_.
        * The workflow now can automatically detect the chemistry type of the data:

            * ``auto`` by default: The workflow checks all possible assay types to decide the correct one.
            * ``threeprime``: The workflow checks all 3' assay types to decide the correct one.
            * ``fiveprime``: The workflow checks all 5' assay types to decide the correct one.
        * Reorg the keywords in *Chemistry* column of sample sheet: 5' now only has 2 types, ``SC5Pv2`` and ``SC5Pv3`` for 5' v2 and v3 chemistries where only R2 is used for alignment.
        * Support the new format of `10x cell barcode inclusion lists`_ provided in Cell Ranger v9.0+.
        * Fix issue in processing UTF-encoded feature barcode files

    * For immune profiling:

        * Add types ``vdj_t``, ``vdj_b`` and ``vdj_t_gd`` for *DataType* column. See `10x 5' Immune Profiling Kit`_ for details.
        * Remove ``chain`` input, as it is now automatically decided by user-specified *DataType* types.
        * For ``vdj_t_gd`` type samples, support the feature of specifying primer sequences used to enrich cDNA for V(D)J sequences. To enable it, provide a ``.txt`` file in *AuxFile* column of the sample, and it will be passed to ``--inner-enrichment-primers`` option of *cellranger vdj* in execution.

    * For Flex Gene Expression:

        * Remove *ProbeSet* column. The probe set is now automatically decided based on user-specified *Reference* name.
        * Support Flex probe sets v1.1 which are associated with 2024-A genome references.

    * For CellPlex using CMO:

        * Remove ``cmo_set`` input. If using custom CMOs in your experiment, just provide the custom feature reference file in *AuxFile* column of the ``cmo`` type sample.


* Spaceranger workflow:

    * Upgrade *spaceranger_version* default to ``3.1.3``.
    * Remove mkfastq related workflow inputs.
    * Remove *run_count* input as it's always ``true``.
    * Remove support on Targeted Gene Expression analysis, as it's no longer supported since Space Ranger v2.1.1.

.. _BCL Convert: https://emea.support.illumina.com/sequencing/sequencing_software/bcl-convert.html
.. _10x cell barcode inclusion lists: https://kb.10xgenomics.com/hc/en-us/articles/115004506263-What-is-a-barcode-inclusion-list-formerly-barcode-whitelist
.. _Sample Multiplexing: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-3p-multi
.. _10x 5' Immune Profiling Kit: https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-5p-multi
