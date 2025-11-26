4.0.0 :small:`Nov 26, 2025`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * Support library-level chemistry for Flex data processing. This is needed when Flex and CITE-Seq libraries have probe barcodes at different read pairs, which causes the auto-detection to fail (e.g. ``MFRP-Ab-R1`` for CITE-Seq library, ``MFRP-RNA`` for Flex library). (PR `446 <https://github.com/lilab-bcb/cumulus/pull/446>`_)
    * Upgrade *cellranger_version* default to ``10.0.0``.
    * Upgrade *cellranger_arc_version** default to ``2.1.0``, with associated important feature changes:

        * The ``20,000`` total cell limit is removed.
        * Add *secondary* input for running secondary analysis or not, with default ``false``.
    * Upgrade *cellranger_atatc_version* default to ``2.2.0``, with associated imporant feature changes:

        * The ``20,000`` total cell limit is removed.
        * The *force_cells* input can be any positive integer, which is no longer restricted to be smaller than ``20,000``.
        * Add *secondary* input for running secondary analysis or not, with default ``false``.
    * Add 2024-A transcriptome references for Cell Ranger ARC and ATAC: ``GRCh38-2024-A_arc`` and ``GRCm39-2024-A_arc``.
* Spaceranger workflow:
    * Upgrade *spaceranger_version* default to ``4.0.1``, with associated important feature changes:

        * Support Visium HD 3'. To enable it, in the sample sheet, leave the *ProbeSet* column blank for such samples, or skip this column if all samples are Visium HD 3'.
        * Add Visium probesets v2.1 which are bundled with 2024-A transcriptome references: ``human_probe_v2.1`` and ``mouse_probe_v2.1``.
