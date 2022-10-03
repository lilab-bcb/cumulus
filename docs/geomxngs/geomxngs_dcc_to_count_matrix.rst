Generate probe count matrix with pathologists' annotation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **geomxngs_dcc_to_count_matrix** workflow generates an area of illumination (AOI) by probe count matrix with patholgoists' annotation from the output of the **geomxngs_fastq_to_dcc** workflow and user inputs.

Workflow Input
++++++++++++++

Workflow inputs are described below (required inputs in bold).


.. list-table::
    :header-rows: 1
    :widths: 5 30 15 15

    * - Name
      - Description
      - Example
      - Default
    * - **dcc_zip**
      - DCC zip file from **geomxngs_fastq_to_dcc** workflow output
      - "gs://foo/bar/out/DCC-20221001.zip"
      -
    * - **ini**
      - Configuration file in INI format, containing pipeline processing parameters
      - "gs://foo/bar/config.ini"
      - 
    * - **lab_worksheet**
      - A text file containing library setups
      - "gs://foo/bar/LabWorksheet.txt"
      -
    * - **dataset**
      - Data QC and annotation file (Excel) downloaded from instrument after uploading DCC zip file; we only use the first tab (SegmentProperties)
      - "gs://foo/bar/BioprobeQC.xlsx"
      -
    * - **pkc**
      - GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers. Options:
        - CTA_v1.0-4 for Cancer Transcriptome Atlas
        - COVID-19_v1.0 for COVID-19 Immune Response Atlas
        - Human_WTA_v1.0 for Human Whole Transcriptome Atlas
        - Mouse_WTA_v1.0 for Mouse Whole Transcriptome Atlas
        If your configuration file is not listed, you can provide a URL to a PKC zip file or PKC file instead.
      - "Human_WTA_v1.0"
      - 
    * - **output_directory**
      - URL to write results
      - "gs://foo/bar/out" or "s3://foo/bar/out"
      - 
    * - backend
      - Backend for computation. Available options:
        - "gcp" for Google Cloud
        - "aws" for Amazon AWS
        - "local" for local machine
      - "aws"
      - "gcp"
    * - docker_registry
      - Docker registry to use for this workflow. Options:

        - "quay.io/cumulus" for images on Red Hat registry;

        - "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - docker_version
      - Docker image version.
      - "1.0.0"
      - "1.0.0"
    * - preemptible
      - Number of preemptible tries
      - 2
      - 2
    * - memory
      - Memory string
      - "8GB"
      - "8GB"
    * - cpu
      - Number of CPUs
      - 1
      - 1
    * - extra_disk_space
      - Extra disk space in GB.
      - 5
      - 5
    * - aws_queue_arn
      - The arn URI of the AWS job queue to be used. Only works when backend is aws
      - "arn:aws:batch:us-east-1:xxx:job-queue/priority-gwf"
      - ""

Workflow Output
+++++++++++++++

.. list-table::
    :header-rows: 1
    :widths: 5 20 5

    * - Name
      - Description
      - Type
    * - count_matrix_h5ad
      - URL to a count matrix in h5ad format. X contains the count matrix, obs contains AOI information, and .var contains probe metadata
      - String
    * - count_matrix_text
      - URL to a count matrix in text format.  Each row is one probe and each column is one AOI. First column is RTS_ID (Readout Tag Sequence-ID (RTS-ID)). Second column is Gene (if multiple probes map to the same gene, their values are the same). Third columns is Probe (if multiple probes map to the same gene, values are different control_1, control_2). Starting from column 4, we have counts.
      - String
    * - count_matrix_metadata
      - URL to a count matrix metadata in text format. All columns from dataset file are included; each row describes one AOI (area of illumination)
      - String
