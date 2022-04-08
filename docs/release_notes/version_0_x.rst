Version 0.15.0 :small:`May 6, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Update all workflows to OpenWDL version 1.0.
* Cumulus now supports multi-job execution from Terra data table input.
* Cumulus generates Cirrocumulus input in ``.cirro`` folder, instead of a huge ``.parquet`` file.

Version 0.14.0 :small:`February 28, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added support for gene-count matrices generation using alternative tools (STARsolo, Optimus, Salmon alevin, Kallisto BUStools).
* Cumulus can process demultiplexed data with remapped singlets names and subset of singlets.
* Update VDJ related inputs in Cellranger workflow.
* SMART-Seq2 and Count workflows are in OpenWDL version 1.0.

Version 0.13.0 :small:`February 7, 2020`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added support for aggregating scATAC-seq samples.
* Cumulus now accepts mtx format input.

Version 0.12.0 :small:`December 14, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added support for building references for sc/snRNA-seq, scATAC-seq, single-cell immune profiling, and SMART-Seq2 data.

Version 0.11.0 :small:`December 4, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Reorganized Cumulus documentation.

Version 0.10.0 :small:`October 2, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* scCloud is renamed to Cumulus.
* Cumulus can accept either a sample sheet or a single file.

Version 0.7.0 :small:`Feburary 14, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added support for 10x genomics scATAC assays.
* scCloud runs FIt-SNE as default.

Version 0.6.0 :small:`January 31, 2019`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added support for 10x genomics V3 chemistry.
* Added support for extracting feature matrix for Perturb-Seq data.
* Added R script to convert output_name.seurat.h5ad to Seurat object. Now the raw.data slot stores filtered raw counts.
* Added min_umis and max_umis to filter cells based on UMI counts.
* Added QC plots and improved filtration spreadsheet.
* Added support for plotting UMAP and FLE.
* Now users can upload their JSON file to annotate cell types.
* Improved documentation.
* Added lightGBM based marker detection.

Version 0.5.0 :small:`November 18, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added support for plated-based SMART-Seq2 scRNA-Seq data.

Version 0.4.0 :small:`October 26, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added CITE-Seq module for analyzing CITE-Seq data.

Version 0.3.0 :small:`October 24, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added the demuxEM module for demultiplexing cell-hashing/nuclei-hashing data.

Version 0.2.0 :small:`October 19, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* Added support for V(D)J and CITE-Seq/cell-hashing/nuclei-hashing.

Version 0.1.0 :small:`July 27, 2018`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

* KCO tools released!
