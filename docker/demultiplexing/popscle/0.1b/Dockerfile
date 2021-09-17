FROM continuumio/miniconda3:4.7.12
SHELL ["/bin/bash", "-c"]

RUN apt-get update --allow-releaseinfo-change && \
    apt-get -qq -y install --no-install-recommends && \
        build-essential \
        curl \
        wget \
        gnupg \
        unzip \
        python3 \
        python3-pip \
        cmake \
        liblzma-dev \
        zlib1g-dev \
        libbz2-dev \
        libcurl4-openssl-dev \
        libssl-dev && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=354.0.0-0

RUN conda install -c bioconda popscle

RUN pip install --upgrade pip && \
    pip install Cython==0.29.23 && \
    pip install setuptools==56.0.0 && \
    pip install importlib-metadata==4.0.1 && \
    pip install docopt==0.6.2 && \
    pip install pandas==1.2.4 && \
    pip install loompy==3.0.6 && \
    pip install anndata==0.7.6 && \
    pip install scipy==1.6.3 && \
    pip install zarr==2.8.1 && \
    pip install natsort==7.1.1 && \
    pip install h5py==3.2.1 && \
    pip install numba==0.53.1 && \
    pip install argparse==1.4.0 && \
    pip install networkx==2.5.1 && \
    pip install pegasusio==0.3.1.post2 && \
    pip install stratocumulus==0.1.1

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/demultiplexing/souporcell/extract_barcodes_from_rna.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/demultiplexing/popscle/generate_zarr.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH
