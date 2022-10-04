Nanostring GeoMx DSP
--------------------

This section contains two workflows: **geomxngs_fastq_to_dcc** and **geomxngs_dcc_to_count_matrix**.

**geomxngs_fastq_to_dcc** workflow wraps `Nanostring GeoMx Digital Spatial NGS Pipeline`_ and can convert FASTQ files into DCC files.

**geomxngs_dcc_to_count_matrix** workflow takes the DCC zip file from **geomxngs_fastq_to_dcc** and other files produced by the GeoMx DSP machine as inputs, and outputs an area of illumination (AOI) by probe count matrix with pathologists' annotation.

---------------------------------

.. include:: geomxngs_fastq_to_dcc.rst

---------------------------------

.. include:: geomxngs_dcc_to_count_matrix.rst
