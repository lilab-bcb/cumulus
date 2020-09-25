FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y build-essential lsb-release gnupg && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=311.0.0-0 && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN conda install -c anaconda pandas==1.1.1

RUN conda install -c bioconda popscle

RUN pip install --upgrade pip && \
    pip install Cython==0.29.21 && \
    pip install importlib-metadata==1.7.0 && \
    pip install docopt==0.6.2 && \
    pip install loompy==3.0.6 && \
    pip install anndata==0.7.4 && \
    pip install scipy==1.5.2 && \
    pip install zarr==2.4.0 && \
    pip install natsort==7.0.1 && \
    pip install h5py==2.10.0 && \
    pip install numba==0.51.2 && \
    pip install argparse==1.4.0 && \
    pip install pegasusio==0.2.5.post2

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/demultiplexing/souporcell/extract_barcodes_from_rna.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/demultiplexing/demuxlet/generate_zarr.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH
