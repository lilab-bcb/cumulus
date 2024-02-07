2.4.1 :small:`May 30, 2023`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Improve:**

* In Cellranger workflow:

    * Fixed RNA Profiling now accepts custom probe set references.

* In STARsolo workflow:

    * Add *limitBAMsortRAM* and *outBAMsortingBinsN* inputs to handle out-of-memory error in BAM sorting phase.

* In GeoMx_fastq_to_dcc workflow:

    * Support multiple FASTQ folders as input.

**Updates:**

* In Demultiplexing workflow:

    * Upgrade *souporcell_version* to ``2022.12``, which is based on commit `9fb527 <https://github.com/wheaton5/souporcell/tree/9fb5271ae9f2257ea9a8552dfda3d4b7080be194>`_ on 2022/12/13.

* In STARsolo workflow:

    * Upgrade *star_version* to ``2.7.10b``.

**Bug Fixes:**

* In Spaceranger workflow:

    * Fix the image localization issue for CytAssist samples.

2.4.0 :small:`January 28, 2023`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Updates:**

* In Cellranger workflow:

    * Upgrade *cellranger_version* default to ``7.1.0``.
    * Add Mouse probe set v1.0 for `Fixed RNA Profiling`_ analysis.
    * Add probe set v1.0.1 of both Human and Mouse for `Fixed RNA Profiling`_ analysis.
    * Upgrade *cumulus_feature_barcoding* version default to ``0.11.1``.

* In Cellranger_create_reference workflow: upgrade *cellranger_version* default to ``7.1.0``.

* In Cellranger_vdj_create_reference workflow: upgrade *cellranger_version* default to ``7.1.0``.

* In Spaceranger workflow: upgrade *spaceranger_version* default to ``2.0.1``.


.. _Fixed RNA Profiling: ./cellranger/index.html#fixed-rna-profiling
