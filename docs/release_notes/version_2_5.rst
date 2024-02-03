2.5.0 :small:`February 3, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Improvements:**

* In Cellranger workflow:

    * Add *multi_disk_space* input for specifying disk size used by Fixed RNA Profiling or 10x multiome jobs, which usually take large amount.

* In Cellbender workflow:

    * With ``0.3.0`` default Cell Bender version, users don't need to always specify *expected_cells* or *total_droplets_included* input.

**Updates:**

* In Cellranger workflow:

    * Upgrade *cellranger_version*  default to ``7.2.0``.
    * Upgrade *cumulus_feature_barcoding* version default to ``0.11.2``.
    * Add GRCh38 VDJ v7.1.0 reference: ``GRCh38_vdj_v7.1.0``.

* In Spaceranger workflow:

    * Upgrade *spaceranger_version* default to ``2.1.1``.

* In Demultiplexing workflow:

    * Upgrade *souporcell_version* default to ``2.5``.

* In Cellbender workflow:

    * Upgrade *cellbender_version* default to ``0.3.0``.
