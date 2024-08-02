2.6.3 :small:`August 2, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Update Demultiplexing workflow to work with Stratocumulus v0.2.4.

2.6.2 :small:`June 19, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Solve the issue of Cellranger workflow with cellranger-arc. ``cellranger_arc_version`` default is now "2.0.2.strato" which is compatible with workflow v2.6.1 or later.

2.6.1 :small:`May 8, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^

* In Cellranger workflow: Add new genome references for single-cell/nucleus RNA-Seq: **GRCh38-2024-A** for human, **GRCm39-2024-A** for mouse, and **GRCh38_and_GRCm39-2024-A** for human and mouse.
* In Spaceranger workflow: Add a new probe set **mouse_probe_v2** for mouse.
* Some underlying workflow improvement.


2.6.0 :small:`April 22, 2024`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**Update:**

* In Cellranger workflow:

    * Upgrade *cellranger_version* default to ``8.0.0``.
    * Upgrade *cumulus_feature_barcoding* version default to ``0.11.3``.

* In Spaceranger workflow:

    * Upgrade *spaceranger_version* default to ``3.0.0``.
    * Support Visium HD data.
