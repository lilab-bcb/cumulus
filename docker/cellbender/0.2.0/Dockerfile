# Start from nvidia-docker image with drivers pre-installed to use a GPU
FROM nvcr.io/nvidia/cuda:11.7.0-base-ubuntu22.04

# Install curl and sudo and git and miniconda and pytorch, cudatoolkit, pytables, and cellbender
RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    wget \
    ca-certificates \
    unzip \
    sudo \
&& sudo apt-get install -y --no-install-recommends \
    git \
    bzip2 \
    python3 \
    python3-pip \
&& sudo rm -rf /var/lib/apt/lists/* \
&& curl -so ~/miniconda.sh https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
&& chmod +x ~/miniconda.sh \
&& ~/miniconda.sh -b -p /home/user/miniconda \
&& rm ~/miniconda.sh

RUN pip3 install --upgrade pip --no-cache-dir && \
    pip3 install stratocumulus==0.1.7 --no-cache-dir

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=385.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.6.4.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

ENV PATH=/home/user/miniconda/bin:$PATH
ENV CONDA_AUTO_UPDATE_CONDA=false

RUN conda install -y pytorch==1.11.0 cudatoolkit=10.2 -c pytorch \
&& conda clean -ya

ENV DOCKER='true'

RUN wget https://github.com/broadinstitute/CellBender/archive/refs/tags/v0.2.0.tar.gz && \
    tar -xzvf v0.2.0.tar.gz && \
    rm v0.2.0.tar.gz && \
    mkdir /software && \
    mv CellBender-0.2.0 /software/CellBender && \
    cd /software/CellBender && \
    pip install -e .

RUN apt-get -qq -y remove curl gnupg wget curl && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    ln -s /usr/bin/python3 /usr/bin/python

ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

# Add cellbender command to PATH
ENV PATH=/software:~/miniconda/bin:$PATH