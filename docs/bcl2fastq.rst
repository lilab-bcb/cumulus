.. _bcl2fastq-docker:

bcl2fastq
-----------

License
^^^^^^^^^
`bcl2fastq license`_

Docker
^^^^^^^^^
Add the lines below into your Dockerfile to install bcl2fastq into a debian based docker::

    RUN apt-get install --no-install-recommends -y alien unzip
    ADD https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-0-linux-x86-64.zip /software
    RUN unzip -d /software/ /software/bcl2fastq2-v2-20-0-linux-x86-64.zip && alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && rm /software/bcl2fastq2-v2*

You can host your private docker images in the `Google Container Registry`_.

Workflows
^^^^^^^^^^^^
Workflows such as cellranger_workflow and dropseq_workflow provide the option of running bcl2fastq. We provide dockers
containing bcl2fastq that are accessible only by members of the Broad Institute. Non-Broad Institute members will have to provide
their own docker images.


Example
^^^^^^^^^
In this example we create a docker image for running cellranger mkfastq version 3.0.2.

- Create a GCP project or reuse an existing project.
- Enable the Google Container Registry
- Clone the cumulus repository::

    git clone https://github.com/klarman-cell-observatory/cumulus.git

- Add the lines to cumulus/docker/cellranger/3.0.2/Dockerfile to include bcl2fastq (see Docker_).
- Ensure you have `Docker installed`_
- Build, tag, and push the docker. Remember to replace PROJECT_ID with your GCP project id::

    cd cumulus/docker/cellranger/3.0.2/
    docker build -t cellranger-3.0.2 .
    docker tag cellranger-3.0.2 gcr.io/PROJECT_ID/cellranger:3.0.2
    gcr.io/PROJECT_ID/cellranger:3.0.2

- Run cell_ranger_workflow, entering your docker registry URL for the input ``cellranger_mkfastq_docker_registry``

.. _`Google Container Registry`: https://cloud.google.com/container-registry/docs/
.. _`bcl2fastq license`: https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-eula.pdf
.. _`Docker installed`: https://www.docker.com/products/docker-desktop

