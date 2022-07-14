2.1.0 :small:`July 13, 2022`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**New Features:**

* Add CellBender_ workflow for ambient RNA removal.
* *CellRanger*:

    * For **ATAC-Seq** data, add ``ARC-v1`` chemistry keyword for analyzing only the ATAC part of 10x multiome data. See `CellRanger scATAC-seq sample sheet`_ section for details.
    * For **antibody/hashing/citeseq/crispr** data, add ``multiome`` chemistry keyword for the feature barcoding on 10x multiome data.
* *STARsolo*:

    * In workflow output, besides ``mtx`` format gene-count matrices, the workflow also generates matrices in 10x-compatible ``hdf5`` format.


**Improvements:**

* *CellRanger*: For **antibody/hashing/citeseq/crispr** data,

    * *cumulus_feature_barcoding* v0.9.0+ now supports multi-threading and faster gzip file I/O.

**Updates:**

* Genome Reference:

    * Add Cellranger VDJ v7.0.0 genome references: ``GRCh38_vdj_v7.0.0`` and ``GRCm38_vdj_v7.0.0`` in `CellRanger scIR-seq sample sheet`_ section.

* Default version upgrade:

    * Update *cellranger* default to `v7.0.0 <https://support.10xgenomics.com/single-cell-gene-expression/software/pipelines/7.0/release-notes>`_.
    * Update *cellranger-atac* default to `v2.1.0 <https://support.10xgenomics.com/single-cell-atac/software/pipelines/2.1/release-notes>`_.
    * Update *cellranger-arc* default to `v2.0.1 <https://support.10xgenomics.com/single-cell-multiome-atac-gex/software/pipelines/2.0/release-notes>`_.
    * Update *cumulus_feature_barcoding* default to `v0.10.0 <https://github.com/lilab-bcb/cumulus_feature_barcoding/releases/tag/0.10.0>`_.
    * *STARsolo* workflow uses STAR v2.7.10a by default.

.. _CellBender: ./cellbender.html
.. _CellRanger scATAC-seq sample sheet: ./cellranger/index.html#id5
.. _CellRanger scIR-seq sample sheet: ./cellranger/index.html#id8
