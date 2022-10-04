2.2.0 :small:`October 4, 2022`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Features:**

* Add `Nanostring GeoMx DSP`_ workflows. It consists of two steps:

    * **GeoMxFastqToDCC** workflow to convert FASTQ files into DCC files by wrapping Nanostring GeoMx Digital Spatial NGS Pipeline.
    * **GeoMxDCCToCountMatrix** workflow to generate probe count matrix from DCC files with pathologists' annotation.

**Updates:**

* Spaceranger workflow:

    * Add support on 10x Space Ranger v2.0.0.
    * Add ``human_probe_v2`` Probe Set for FFPE samples, which is compatible with CytAssist FFPE samples.

* Upgrade ``cumulus_feature_barcoding_version`` default to `v0.11.0`_ for Feature Barcoding in Cellranger workflow.

**API Changes:**

* Across all workflows, for AWS backend:

    * All workflows now have ``awsQueueArn`` input, which is used for explicitly specifying the Arn string of an AWS Compute Environment.
    * Remove ``awsMaxRetries`` input for all workflows. Namely, use Cromwell's default value ``0``.

**Bug Fix:**

* Fix the issue on localizing GCP folders in Cellranger workflow for ATAC-Seq and 10x Multiome data.


.. _Nanostring GeoMx DSP: ./geomxngs/index.html
.. _v0.11.0: https://github.com/lilab-bcb/cumulus_feature_barcoding/releases/tag/0.11.0
