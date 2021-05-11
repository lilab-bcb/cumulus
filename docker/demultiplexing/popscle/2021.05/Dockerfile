FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get -qq -y install --no-install-recommends \
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
        libssl-dev

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=339.0.0-0

RUN mkdir /software && \
    wget https://github.com/samtools/htslib/releases/download/1.12/htslib-1.12.tar.bz2 && \
    tar -xvjf htslib-1.12.tar.bz2 && \
    mv htslib-1.12 /software/ && \
    rm htslib-1.12.tar.bz2 && \
    cd /software/htslib-1.12 && \
    make && \
    make install

RUN wget https://github.com/statgen/popscle/archive/da70fc78da385ef049e0e890342acfd62842cae0.zip -O popscle.zip && \
    unzip popscle.zip && \
    rm popscle.zip && \
    mv popscle-da70fc78da385ef049e0e890342acfd62842cae0 /software/popscle && \
    cd /software/popscle && \
    mkdir build && \
    cd build && \
    cmake .. && \
    make

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install Cython==0.29.23 --no-cache-dir && \
    python -m pip install setuptools==56.0.0 --no-cache-dir && \
    python -m pip install importlib-metadata==4.0.1 --no-cache-dir && \
    python -m pip install docopt==0.6.2 --no-cache-dir && \
    python -m pip install pandas==1.2.4 --no-cache-dir && \
    python -m pip install loompy==3.0.6 --no-cache-dir && \
    python -m pip install anndata==0.7.6 --no-cache-dir && \
    python -m pip install scipy==1.6.3 --no-cache-dir && \
    python -m pip install zarr==2.8.1 --no-cache-dir && \
    python -m pip install natsort==7.1.1 --no-cache-dir && \
    python -m pip install h5py==3.2.1 --no-cache-dir && \
    python -m pip install numba==0.53.1 --no-cache-dir && \
    python -m pip install argparse==1.4.0 --no-cache-dir && \
    python -m pip install networkx==2.5.1 --no-cache-dir && \
    python -m pip install pegasusio==0.2.10 --no-cache-dir

RUN apt-get -qq -y remove curl gnupg wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/demultiplexing/souporcell/extract_barcodes_from_rna.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/demultiplexing/popscle/generate_zarr.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software/popscle/bin:/software:$PATH
