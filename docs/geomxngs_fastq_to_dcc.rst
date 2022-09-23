geomxngs_fastq_to_dcc
-----------------------------------------------------------------------------------------
Processes FASTQ sequencing files into digital count conversion (DCC) files using NanoString's GeoMx NGS Pipeline

Non Broad Institute users that wish to run geomxngs_fastq_to_dcc must create a custom docker image containing the
NanoString's GeoMx NGS software.

Docker Instructions (for Non-Broad Institute Users)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Please note that if you're a Broad Institute member and are not able to pull the docker image, please check
https://app.terra.bio/#groups to see that you're a member of the all_broad_users group. If not, please contact
Terra support and ask to be added to the all_broad_users@firecloud.org group.

Download NanoString's GeoMx NGS Pipeline software from `the NanoString website <https://nanostring.com/products/geomx-digital-spatial-profiler/software-updates/>`_,
which requires registration. Download the Linux zip file (for example GeoMxNGSPipeline_Linux_2.3.3.10.zip)
to the same directory as your Dockerfile.

You can host your private docker images in the `Google Artifact Registry`_.

Docker Example
^^^^^^^^^^^^^^^
In this example we create a docker image for running ``geomxngs_fastq_to_dcc`` using GeoMx NGS Pipeline software 2.3.3.10.

#. Create a GCP project or reuse an existing project.
#. Enable the Google Artifact Registry
#. Clone the cumulus repository::

    git clone https://github.com/lilab-bcb/cumulus.git

#. Ensure you have `Docker installed`_
#. Download the GeoMx NGS Pipeline software from `the NanoString website <https://nanostring.com/products/geomx-digital-spatial-profiler/software-updates/>`_,
#. Build, tag, and push the docker. Remember to replace us-central1-docker.pkg.dev/my-lab-1234 with your docker registry URL::

    cd cumulus/docker/geomxngs_fastq_to_dcc/
    docker build --build-arg GEOMX_ZIP=GeoMxNGSPipeline_Linux_2.3.3.10.zip --platform linux/x86_64 -t geomxngs_fastq_to_dcc-2.3.3.10 .
    docker tag geomxngs_fastq_to_dcc-2.3.3.10 us-central1-docker.pkg.dev/my-lab-1234/geomxngs_fastq_to_dcc:2.3.3.10
    docker push us-central1-docker.pkg.dev/my-lab-1234/geomxngs_fastq_to_dcc:2.3.3.10



Workflow Input
^^^^^^^^^^^^^^
Workflow inputs are described below (required inputs in bold).


.. list-table::
    :header-rows: 1
    :widths: 5 20 5

    * - Name
      - Description
      - Default
    * - **dcc_zip**
      - DCC zip file from **geomxngs_fastq_to_dcc** workflow output
      -
    * - **ini**
      - Configuration file in INI format, containing pipeline processing parameters
      -
    * - **dataset**
      - Data QC and annotation file (Excel) downloaded from instrument after uploading DCC zip file; we only use the first tab (SegmentProperties)
      -
    * - **pkc**
      - GeoMx DSP configuration file to associate assay targets with GeoMx HybCode barcodes and Seq Code primers
         - Options:
            - CTA_v1.0-4 for Cancer Transcriptome Atlas
            - COVID-19_v1.0 for COVID-19 Immune Response Atlas
            - Human_WTA_v1.0 for Human Whole Transcriptome Atlas
            - Mouse_WTA_v1.0 for Mouse Whole Transcriptome Atlas
        If your configuration file is not listed, you can provide a URL to a PKC zip file or PKC file instead.
      -
    * - **lab_worksheet**
      - A text file containing library setups
      -
    * - **output_directory**
      - URL to write results (e.g. s3://foo/bar/out or gs://foo/bar/out)
      -
    * - backend
      - Backend for computation. Available options:
        - "gcp" for Google Cloud
        - "aws" for Amazon AWS
        - "local" for local machine
      - "gcp"
    * - docker_registry
      - Docker registry
      - quay.io/cumulus
    * - docker_version
      - Docker image version.
      - "1.0.0"
    * - preemptible
      - Number of preemptible tries.
      - 2
    * - memory
      - Memory string.
      - "8GB"
    * - cpu
      - Number of CPUs.
      - 1
    * - extra_disk_space
      - Extra disk space in GB.
      - 5
    * - aws_queue_arn
      - The arn URI of the AWS job queue to be used (e.g. arn:aws:batch:us-east-1:xxxxx). Only works when backend is aws.
      -


Workflow Output
^^^^^^^^^^^^^^^^

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


.. _`Google Artifact Registry`: https://cloud.google.com/artifact-registry
.. _`Docker installed`: https://www.docker.com/products/docker-desktop