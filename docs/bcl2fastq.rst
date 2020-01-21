.. _bcl2fastq-docker:

bcl2fastq
-----------

License
^^^^^^^^^
`bcl2fastq license`_

Docker
^^^^^^^^^

Read `this tutorial <https://docs.docker.com/get-started/>`_ if you are new to Docker and don't know how to write your Dockerfile.

First, you need to download ``bcl2fastq`` software from `its official website <https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html>`_, which requires your registration. After registration, choose its ``Linux rpm`` format file for downloading. Unzip the file and move the ``bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm`` to the same directory where your Dockerfile will be saved.

Then for a Debian-based docker (e.g. `continuumio/miniconda3 <https://hub.docker.com/r/continuumio/miniconda3>`_), add the following chunks below into its Dockerfile to install ``bcl2fastq``

.. code-block:: Docker

    FROM continuumio/miniconda3:4.7.12

.. code-block:: Docker

    # Install Google Cloud SDK
    RUN apt-get update \
        && apt-get install --no-install-recommends -y build-essential dpkg-dev gnupg lsb-release procps curl \
        && export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" \
        && echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
        && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
        && apt-get update && apt-get install -y google-cloud-sdk

.. code-block:: Docker

    # Add Bcl2Fastq: Make sure the downloaded Bcl2Fastq RPM file is in the same directory as your Dockerfile 
    ADD bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm /software/
    RUN apt-get update && apt-get install --no-install-recommends -y alien \
        && alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm

Make sure ``bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm`` is located in the same directory as your ``Dockerfile`` file.

You can host your private docker images in the `Google Container Registry`_.

Workflows
^^^^^^^^^^^^
Workflows such as **cellranger_workflow** and **dropseq_workflow** provide the option of running ``bcl2fastq``. We provide dockers
containing ``bcl2fastq`` that are accessible only by members of the Broad Institute. Non-Broad Institute members will have to provide
their own docker images.


Example
^^^^^^^^^
In this example we create a docker image for running ``cellranger mkfastq`` version 3.0.2.

#. Create a GCP project or reuse an existing project.
#. Enable the Google Container Registry
#. Clone the cumulus repository::

    git clone https://github.com/klarman-cell-observatory/cumulus.git

#. Add the lines to cumulus/docker/cellranger/3.0.2/Dockerfile to include bcl2fastq (see Docker_).
#. Ensure you have `Docker installed`_
#. Download cellranger from https://support.10xgenomics.com/single-cell-gene-expression/software/downloads/3.0
#. Build, tag, and push the docker. Remember to replace PROJECT_ID with your GCP project id::

    cd cumulus/docker/cellranger/3.0.2/
    docker build -t cellranger-3.0.2 .
    docker tag cellranger-3.0.2 gcr.io/PROJECT_ID/cellranger:3.0.2
    gcr.io/PROJECT_ID/cellranger:3.0.2

#. Import **cellranger_workflow** workflow to your workspace (see `cellranger_workflow steps <./cellranger.html>`_), and enter your docker registry URL (in this example, ``"gcr.io/PROJECT_ID/"``) in ``cellranger_mkfastq_docker_registry`` field of `cellranger_workflow inputs <./cellranger.html#cellranger-workflow-inputs>`_.

.. _`Google Container Registry`: https://cloud.google.com/container-registry/docs/
.. _`bcl2fastq license`: https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-eula.pdf
.. _`Docker installed`: https://www.docker.com/products/docker-desktop

