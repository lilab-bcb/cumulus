# Start from nvidia-docker image with drivers pre-installed to use a GPU
FROM nvcr.io/nvidia/cuda:11.7.1-base-ubuntu22.04

LABEL description="This docker is modified from the official BSD-3-Clause licensed CellBender dockerfile (https://github.com/broadinstitute/CellBender/blob/master/docker/Dockerfile)."

# Install curl and sudo and git and miniconda and pytorch, cudatoolkit, pytables, and cellbender
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        curl \
        wget \
        rsync \
        ca-certificates \
        unzip

RUN curl -so miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-py37_23.1.0-1-Linux-x86_64.sh && \
    bash miniconda.sh -b -p /opt/miniconda3 && \
    rm miniconda.sh

ENV PATH=/opt/miniconda3/bin:$PATH \
    DOCKER=true \
    CONDA_AUTO_UPDATE_CONDA=false

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=483.0.0-0

ENV CLOUDSDK_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python
ENV CLOUDSDK_GSUTIL_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.15.10.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN pip install stratocumulus==0.2.4

RUN wget https://github.com/broadinstitute/CellBender/archive/refs/tags/v0.3.0.tar.gz && \
    tar -xzf v0.3.0.tar.gz && \
    rm v0.3.0.tar.gz && \
    mkdir /software && \
    mv CellBender-0.3.0 /software/CellBender && \
    cd /software/CellBender && \
    pip install -e .

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

# Add cellbender command to PATH
ENV PATH=/software:$PATH
