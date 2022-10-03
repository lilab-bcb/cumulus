Convert FASTQ files into DCC files by the Nanostring GeoMx Digital Spatial NGS Pipeline
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The **geomxngs_fastq_to_dcc** workflow converts FASTQ files to DCC files by wrapping the `Nanostring GeoMx Digital Spatial NGS Pipeline`_.
After generating DCC files, use the **geomxngs_dcc_to_count_matrix** workflow to generate an area of interest by probe count matrix.

Workflow Input
++++++++++++++

Relevant workflow inputs are described below (required inputs in bold)


.. list-table::
    :header-rows: 1
    :widths: 5 30 15 15

    * - Name
      - Description
      - Example
      - Default
    * - **fastq_directory**
      - FASTQ directory URL
      - "gs://foo/bar/fastqs" or "s3://foo/bar/fastqs"
      - 
    * - **ini**
      - Configuration file in INI format, containing pipeline processing parameters
      - "gs://foo/bar/config.ini"
      -
    * - **output_directory**
      - URL to write results
      - "gs://foo/bar/out" or "s3://foo/bar/out"
      - 
    * - fastq_rename
      - Optional 2 column TSV file with no header used to map original FASTQ names to FASTQ names that GeoMX recognizes.
      - "gs://foo/bar/fastq_rename.tsv"
      - 
    * - delete_fastq_directory
      - Whether to delete the input fastqs upon successful completion
      - true
      - false
    * - geomxngs_version
      - Version of the geomx software, currently only "2.3.3.10".
      - "2.3.3.10"
      - "2.3.3.10"
    * - docker_registry
      - Docker registry to use for this workflow. Options:

        - "quay.io/cumulus" for images on Red Hat registry;

        - "cumulusprod" for backup images on Docker Hub.
      - "quay.io/cumulus"
      - "quay.io/cumulus"
    * - backend
      - Backend for computation. Available options:
		    - "gcp" for Google Cloud
		    - "aws" for Amazon AWS
		    - "local" for local machine
      - "aws"
      - "gcp"
    * - zones
      - Google cloud zones
      - "us-central1-a"
      - "us-central1-a us-central1-b us-central1-c us-central1-f"
    * - preemptible
      - Number of preemptible tries
      - 2
      - 2
    * - memory
      - Memory string
      - "64GB"
      - "64GB"
    * - cpu
      - Number of CPUs
      - 4
      - 4
    * - disk_space
      - Disk space in GB
      - 500
      - 500
    * - aws_queue_arn
      - The arn URI of the AWS job queue to be used. Only works when backend is aws.
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
    * - dcc_zip
      - URL to the output DCC zip file
      - String
    * - geomxngs_output
      - URL to the output of geomxngspipeline; the DCC zip file is part of the output here
      - String

.. _`Nanostring GeoMx Digital Spatial NGS Pipeline`: https://nanostring.com/products/geomx-digital-spatial-profiler/software-updates/
