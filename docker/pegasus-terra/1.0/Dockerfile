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
    python3 -m pip install numpy==1.19.0 && \
    python3 -m pip install Cython==0.29.21 && \
    python3 -m pip install pybind11==2.5.0 && \
    python3 -m pip install h5py==2.10.0 && \
    python3 -m pip install fitsne==1.1.1 && \
    python3 -m pip install joblib==0.16.0 && \
    python3 -m pip install leidenalg==0.8.1 && \
    python3 -m pip install lightgbm==2.2.1 && \
    python3 -m pip install llvmlite==0.33.0 && \
    python3 -m pip install loompy==3.0.6 && \
    python3 -m pip install natsort==7.0.1 && \
    python3 -m pip install numba==0.50.1 && \
    python3 -m pip install pandas==1.0.5 && \
    python3 -m pip install scikit-learn==0.23.1 && \
    python3 -m pip install scikit-misc==0.1.3 && \
    python3 -m pip install scipy==1.5.1 && \
    python3 -m pip install statsmodels==0.11.1 && \
    python3 -m pip install torch==1.5.1 && \
    python3 -m pip install anndata==0.7.4 && \
    python3 -m pip install harmony-pytorch==0.1.3 && \
    python3 -m pip install hnswlib==0.4.0 && \
    python3 -m pip install fisher==0.1.9 && \
    python3 -m pip install python-igraph==0.8.2 && \
    python3 -m pip install louvain==0.7.0 && \
    python3 -m pip install MulticoreTSNE-modified==0.1 && \
    python3 -m pip install umap-learn==0.4.6 && \
    python3 -m pip install bokeh==2.1.1 && \
    python3 -m pip install scplot==0.0.16 && \
    python3 -m pip install importlib-metadata==1.7.0 && \
    python3 -m pip install zarr==2.4.0 && \
    python3 -m pip install pegasusio==0.2.4 && \
    python3 -m pip install pegasuspy==1.0.0rc1 && \
    python3 -m pip install cirrocumulus==1.1.1

ENV USER jupyter-user
USER $USER
#we want pip to install into the user's dir when the notebook is running
ENV PIP_USER=true

ENTRYPOINT ["/usr/local/bin/jupyter", "notebook"]
