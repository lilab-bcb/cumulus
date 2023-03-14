.. _bcl2fastq-docker:

bcl2fastq
-----------

License
^^^^^^^^^
`bcl2fastq license`_

Workflows
^^^^^^^^^^^^
Workflows such as **cellranger_workflow**, **spaceranger_workflow** and **dropseq_workflow** provide the option of running ``bcl2fastq``. We provide dockers
containing ``bcl2fastq`` that are accessible only by members of the Broad Institute. Non-Broad Institute members will have to provide
their own docker images. Please note that if you're a Broad Institute member and are not able to pull the docker image, please check
https://app.terra.bio/#groups to see that you're a member of the all_broad_users group. If not, please contact
Terra support and ask to be added to the all_broad_users@firecloud.org group.

Example
^^^^^^^^^
In this example, we create a docker image on Google Cloud (GCP) for running ``cellranger mkfastq`` version 7.1.0. For AWS users, you should use AWS Elastic Container Registry (ECR) for this purpose.

#. Create a GCP project or reuse an existing project.
#. Enable the `Google Container Registry`_.
#. Ensure you have `Docker installed`_
#. Prepare your ``Dockerfile`` with the following content::

    FROM cumulusprod/cellranger:7.1.0
    SHELL ["/bin/bash", "-c"]

    RUN apt-get update && \
        apt-get install -y --no-install-recommends zlib1g-dev

    ADD bcl2fastq-v2-20-0-tar.zip /software/
    RUN cd /software && \
        unzip -d /software/ /software/bcl2fastq2-v2-20-0-tar.zip && \
        tar -zxf /software/bcl2fastq2-v2.20.0.422-Source.tar.gz

    ENV C_INCLUDE_PATH=/usr/include/x86_64-linux-gnu
    ENV INSTALL_DIR=/usr/local/bcl2fastq
    ENV SOURCE=/software/bcl2fastq
    ENV BUILD=/software/bcl2fastq-build

    RUN mkdir ${BUILD} && \
    cd ${BUILD} && \
    chmod ugo+x ${SOURCE}/src/configure && \
    chmod ugo+x ${SOURCE}/src/cmake/bootstrap/installCmake.sh && \
    ${SOURCE}/src/configure --prefix=${INSTALL_DIR} && \
    make && \
    make install && \
    rm -rf /software/bcl2fastq-build

    ENV PATH=$INSTALL_DIR/bin:$PATH

#. From `Illumina website`_, download *bcl2fastq* **Linux tarball** format source code to the same folder where your ``Dockerfile`` lives.
#. In the same folder where your ``Dockerfile`` lives, build, tag, and push the docker image. Remember to replace ``PROJECT_ID`` with your GCP project id::

    docker build -t cellranger:7.1.0 .
    docker tag cellranger:7.1.0 gcr.io/PROJECT_ID/cellranger:7.1.0
    docker push gcr.io/PROJECT_ID/cellranger:7.1.0

#. Import **cellranger_workflow** workflow to your workspace (see `cellranger_workflow steps <./cellranger/index.html>`_), and enter your docker registry URL (in this example, ``"gcr.io/PROJECT_ID"``) in ``mkfastq_docker_registry`` field of `cellranger_workflow inputs <./cellranger/index.html#workflow-input>`_.

Similarly for other workflows, just change ``FROM`` part in your *Dockerfile*. We provide a list of images containing only open-source softwares on Docker Hub under `cumulusprod <https://hub.docker.com/u/cumulusprod>`_ organization.

AWS users simply need to push the docker images to their AWS ECR registry, i.e. replace ``gcr.io/PROJECT_ID`` by ECR project ID (e.g. ``ACCOUNT_ID.dkr.ecr.REGION.amazonaws.com``, where ``ACCOUNT_ID`` and ``REGION`` should be replaced by the actual AWS account ID and region).

Build with bcl2fastq rpm package
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
On `Illumina website`_, there is also a **Linux rpm** format package of *bcl2fastq* which has a much smaller size.

You may switch to a RedHat/Fedora image base to build your bcl2fastq docker with this rpm package, as all Cumulus docker images are based on Debian.

If using a Debian based image, however, this way no longer works for Debian 11 (Bullseye) or later, equivalently Ubuntu 20.04 or later. Thus you need to switch to an earlier image base, as all Cumulus docker images built since 2021 are all based on Debian 11.

Below shows an example code on installing bcl2fastq from its Linux rpm package in ``Dockerfile``::

    RUN apt-get update && apt-get install --no-install-recommends -y alien unzip
    ADD bcl2fastq2-v2-20-0-linux-x86-64.zip /software/
    RUN unzip -d /software/ /software/bcl2fastq2-v2-20-0-linux-x86-64.zip && alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && rm /software/bcl2fastq2-v2*

Besides, you also need to install Google Cloud CLI and 10x Cell Ranger (or Space Ranger for **spaceranger_workflow**) in your ``Dockerfile`` by yourself.

.. _Google Container Registry: https://cloud.google.com/container-registry/docs/
.. _bcl2fastq license: https://support.illumina.com/content/dam/illumina-support/documents/downloads/software/bcl2fastq/bcl2fastq2-v2-20-eula.pdf
.. _Docker installed: https://www.docker.com/products/docker-desktop
.. _Illumina website: https://support.illumina.com/downloads/bcl2fastq-conversion-software-v2-20.html
