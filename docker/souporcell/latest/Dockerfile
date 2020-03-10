FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y \
        wget \ 
        build-essential \
        curl \
        libncurses5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        git

RUN cd /opt && \
    wget -O minimap2.tar.gz https://github.com/lh3/minimap2/archive/v2.7.tar.gz && \
    mkdir minimap2 && \
    tar -zxvf minimap2.tar.gz -C minimap2 --strip-components 1 && \
    rm minimap2.tar.gz && \
    cd minimap2 && \
    make

ENV PATH=/opt/minimap2/:$PATH

RUN cd /opt && \
    wget -O bedtools.tar.gz https://github.com/arq5x/bedtools2/releases/download/v2.28.0/bedtools-2.28.0.tar.gz && \
    mkdir bedtools2 && \
    tar -zxvf bedtools.tar.gz -C bedtools2 --strip-components 1 && \
    rm bedtools.tar.gz && \
    cd bedtools2 && \
    make

ENV PATH=/opt/bedtools2/bin:$PATH

ENV CARGO_HOME=/opt/.cargo
ENV RUSTUP_HOME=/opt/.cargo

RUN bash -c 'curl https://sh.rustup.rs -sSf | sh -s -- -y'

ENV PATH=/opt/.cargo/bin:$PATH

RUN cd /opt && \
    git clone https://github.com/wheaton5/souporcell.git && \
    cd souporcell/troublet && \
    cargo build --release && \
    cd /opt/souporcell/souporcell && \
    cargo build --release

ENV PATH=/opt/souporcell:/opt/souporcell/troublet/target/release:$PATH

RUN conda install -y python==3.6.10 && \
    pip install pysam && \
    pip install pyvcf && \
    pip install numpy && \
    pip install pandas && \
    pip install scipy && \
    pip install pystan==2.17.1.0 && \
    pip install pyfaidx

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
    wget -O freebayes https://github.com/ekg/freebayes/releases/download/v1.3.1/freebayes-v1.3.1 && \
    chmod 777 freebayes

RUN cd /opt && \
    wget -O vartrix.tar.gz https://github.com/10XGenomics/vartrix/releases/download/v1.1.3/vartrix-v1.1.3-x86_64-linux.tar.gz && \
    mkdir vartrix && \
    tar -zxvf vartrix.tar.gz -C vartrix --strip-components 1 && \
    rm vartrix.tar.gz

ENV PATH=/opt/vartrix:/opt:$PATH