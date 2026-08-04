4.1.0 :small:`Aug 4, 2026`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * Upgrade *cumulus_feature_barcoding* default to ``2.1.0`` to improve memory efficiency.
    * Upgrade *cellranger_version* default to ``10.1.0``.
    * Upgrade *cellranger_arc_version* default to ``2.2.0``.
    * Fix a bug when processing Flex+VDJ modality. (PR `469 <https://github.com/lilab-bcb/cumulus/pull/469>`_)
* Spaceranger workflow:
    * Upgrade *spaceranger_version* default to ``4.1.0``, which supports Visium HD 11mm slides.

4.0.4 :small:`Apr 14, 2026`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * Support Flex v2 data processing. Users need to specify ``flex-v2`` in ``DataType`` column of the input sample sheet. ``flex-v1`` and ``frp`` are for Flex v1.
* SmartSeq2 workflow:
    * Default resources are stored in ``gs://cumulus-ref`` bucket in GCP region ``us-central1``.

4.0.3 :small:`Apr 6, 2026`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * Fix a bug in processing 10x On-Chip Multiplexing (OCM) data.
    * Support input FASTQ files in unzipped ``.fastq`` format.
* Demultiplexing workflow:
    * For Souporcell, add *souporcell_umi_tag* to allow users specify a custom UMI tag other than the default ``UB`` for their data.

4.0.2 :small:`Jan 6, 2026`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * Bug fix on processing VDJ data.

4.0.1 :small:`Dec 14, 2025`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * For Flex data, avoid using genome reference when ``no_bam = true`` and ``cellranger_version >= "8.0.0"``.


4.0.0 :small:`Nov 26, 2025`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Cellranger workflow:
    * Support library-level chemistry for Flex data processing. This is needed when Flex and CITE-Seq libraries have probe barcodes at different read pairs, which causes the auto-detection to fail (e.g. ``MFRP-Ab-R1`` for CITE-Seq library, ``MFRP-RNA`` for Flex library). (PR `446 <https://github.com/lilab-bcb/cumulus/pull/446>`_)
    * Upgrade *cellranger_version* default to ``10.0.0``.
    * Upgrade *cellranger_arc_version** default to ``2.1.0``, with associated important feature changes:

        * The ``20,000`` total cell limit is removed.
        * Add *secondary* input for running secondary analysis or not, with default ``false``.
    * Upgrade *cellranger_atatc_version* default to ``2.2.0``, with associated important feature changes:

        * The ``20,000`` total cell limit is removed.
        * The *force_cells* input can be any positive integer, which is no longer restricted to be smaller than ``20,000``.
        * Add *secondary* input for running secondary analysis or not, with default ``false``.
    * Add 2024-A transcriptome references for Cell Ranger ARC and ATAC: ``GRCh38-2024-A_arc`` and ``GRCm39-2024-A_arc``.
* Spaceranger workflow:
    * Upgrade *spaceranger_version* default to ``4.0.1``, with associated important feature changes:

        * Support Visium HD 3'. To enable it, in the sample sheet, leave the *ProbeSet* column blank for such samples, or skip this column if all samples are Visium HD 3'.
        * Involve a new nucleus and cell segmentation algorithm for Visium HD data to infer cells instead of using bins. Check the corresponding section for how to enable this feature.
        * Add Visium probesets v2.1 which are bundled with 2024-A transcriptome references: ``human_probe_v2.1`` and ``mouse_probe_v2.1``.
