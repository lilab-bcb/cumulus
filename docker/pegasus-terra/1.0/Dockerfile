FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:0.0.11
USER root
#this makes it so pip runs as root, not the user
ENV PIP_USER=false

RUN apt-get update && apt-get install -yq --no-install-recommends \
        build-essential \
        automake \
        zlib1g-dev \
        libxml2-dev \
        cmake \
        python3.7-dev \
        gnupg \
        lsb-release \
        libfftw3-dev \
        default-jdk && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 -V && \
    pip3 install --upgrade pip && \
    pip3 install numpy==1.19.0 && \
    pip3 install Cython==0.29.21 && \
    pip3 install pybind11==2.5.0 && \
    pip3 install h5py==2.10.0 && \
    pip3 install fitsne==1.1.1 && \
    pip3 install joblib==0.16.0 && \
    pip3 install leidenalg==0.8.1 && \
    pip3 install lightgbm==2.2.1 && \
    pip3 install llvmlite==0.33.0 && \
    pip3 install loompy==3.0.6 && \
    pip3 install natsort==7.0.1 && \
    pip3 install numba==0.50.1 && \
    pip3 install pandas==1.0.5 && \
    pip3 install scikit-learn==0.23.1 && \
    pip3 install scikit-misc==0.1.3 && \
    pip3 install scipy==1.5.1 && \
    pip3 install statsmodels==0.11.1 && \
    pip3 install torch==1.5.1 && \
    pip3 install anndata==0.7.4 && \
    pip3 install harmony-pytorch==0.1.4 && \
    pip3 install hnswlib==0.4.0 && \
    pip3 install fisher==0.1.9 && \
    pip3 install python-igraph==0.8.2 && \
    pip3 install louvain==0.7.0 && \
    pip3 install MulticoreTSNE-modified==0.1 && \
    pip3 install umap-learn==0.4.6 && \
    pip3 install importlib-metadata==1.7.0 && \
    pip3 install zarr==2.4.0 && \
    pip3 install pegasusio==0.2.4.post1 && \
    pip3 install pegasuspy==1.0.0rc3 && \
    pip3 install cirrocumulus==1.1.1

RUN wget https://bin.equinox.io/c/4VmDzA7iaHb/ngrok-stable-linux-amd64.zip && \
    unzip ngrok-stable-linux-amd64.zip && \
    rm ngrok-stable-linux-amd64.zip && \
    mkdir -p /software && \
    mv ngrok /software/

ENV PATH=/software:$PATH

ENV USER jupyter-user
USER $USER
#we want pip to install into the user's dir when the notebook is running
ENV PIP_USER=true

ENTRYPOINT ["/usr/local/bin/jupyter", "notebook"]
