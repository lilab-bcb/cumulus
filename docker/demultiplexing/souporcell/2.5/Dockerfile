FROM debian:bullseye-slim

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        gnupg \
        ca-certificates \
        automake \
        build-essential \
        bzip2 \
        wget \
        rsync \
        git \
        default-jre \
        curl \
        python3 \
        python3-pip \
        python3-dev \
        libncurses5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        zip \
        unzip \
        pkg-config

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=452.0.1-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.13.30.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN cd /opt && \
    wget https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download && \
    mv download hisat2.zip && \
    unzip hisat2.zip && \
    wget https://github.com/lh3/minimap2/releases/download/v2.26/minimap2-2.26_x64-linux.tar.bz2 && \
    tar -xf minimap2-2.26_x64-linux.tar.bz2 && \
    wget https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    tar -xf bedtools-2.28.0.tar.gz && \
    cd bedtools2 && \
    make && \
    cd .. && \
    rm minimap2-2.26_x64-linux.tar.bz2 hisat2.zip bedtools-2.28.0.tar.gz && \
    cd ..

ENV PATH=/opt/hisat2-2.2.1:/opt/minimap2-2.26_x64-linux:/opt/bedtools2/bin:$PATH
ENV CARGO_HOME=/opt/.cargo
ENV RUSTUP_HOME=/opt/.cargo

RUN bash -c 'curl https://sh.rustup.rs -sSf | sh -s -- -y'

ENV PATH=/opt/.cargo/bin:$PATH

RUN cd /opt && \
    wget https://github.com/wheaton5/souporcell/archive/refs/tags/v2.5.tar.gz && \
    tar -xzf v2.5.tar.gz && \
    rm v2.5.tar.gz && \
    mv souporcell-2.5 souporcell && \
    cd souporcell/troublet && \
    cargo build --release && \
    cd /opt/souporcell/souporcell && \
    cargo build --release

ENV PATH=/opt/souporcell:/opt/souporcell/troublet/target/release:$PATH

RUN python -m pip install numpy==1.26.4 && \
    python -m pip install pandas==2.2.2 && \
    python -m pip install pysam==0.22.0 && \
    python -m pip install git+https://github.com/stan-dev/pystan2.git@master && \
    python -m pip install pyfaidx==0.7.2.2 && \
    python -m pip install "setuptools<58" --upgrade && \
    python -m pip install pyvcf==0.6.8 && \
    python -m pip install scipy==1.11.3 && \
    python -m pip install anndata==0.10.3 && \
    python -m pip install zarr==2.16.1 && \
    python -m pip install networkx==3.2.1 && \
    python -m pip install pegasusio==0.8.1 && \
    python -m pip install stratocumulus==0.2.4

RUN cd /opt && \
    wget -O htslib.tar.bz2 https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 && \
    mkdir htslib && \
    tar -xvjf htslib.tar.bz2 -C htslib --strip-components 1 && \
    rm htslib.tar.bz2 && \
    cd htslib && \
    ./configure && \
    make && \
    make install

RUN cd /opt && \
    wget -O samtools.tar.bz2 https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 && \
    mkdir samtools && \
    tar -xvjf samtools.tar.bz2 -C samtools --strip-components 1 && \
    rm samtools.tar.bz2 && \
    cd samtools && \
    ./configure && \
    make && \
    make install

RUN cd /opt && \
    wget -O bcftools.tar.bz2 https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 && \
    mkdir bcftools && \
    tar -xvjf bcftools.tar.bz2 -C bcftools --strip-components 1 && \
    rm bcftools.tar.bz2 && \
    cd bcftools && \
    ./configure && \
    make && \
    make install

RUN cd /opt && \
    wget https://github.com/freebayes/freebayes/releases/download/v1.3.6/freebayes-1.3.6-linux-amd64-static.gz && \
    gunzip freebayes-1.3.6-linux-amd64-static.gz && \
    mv freebayes-1.3.6-linux-amd64-static freebayes && \
    chmod 777 freebayes

RUN cd /opt && \
    wget https://github.com/10XGenomics/vartrix/releases/download/v1.1.22/vartrix_linux && \
    mv vartrix_linux vartrix && \
    chmod 777 vartrix

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/demultiplexing/souporcell/extract_barcodes_from_rna.py /opt
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/demultiplexing/souporcell/match_donors.py /opt

ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /opt
RUN chmod a+rx /opt/monitor_script.sh

ENV PATH=/opt/freebayes/build:/opt/vartrix:/opt:/usr/local/aws-cli/v2/current/bin:$PATH
ENV TMPDIR=/tmp
