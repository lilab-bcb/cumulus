FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:0.0.8
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
    pip3 install numpy==1.16.4 && \
    pip3 install Cython==0.29.13 && \
    pip3 install pybind11==2.3.0 && \
    pip3 install h5py==2.9.0 && \
    pip3 install fitsne==1.0.1 && \
    pip3 install joblib==0.13.2 && \
    pip3 install leidenalg==0.7.0 && \
    pip3 install lightgbm==2.2.1 && \
    pip3 install llvmlite==0.31.0 && \
    pip3 install loompy==2.0.17 && \
    pip3 install natsort==6.0.0 && \
    pip3 install numba==0.47.0 && \
    pip3 install pandas==0.24.2 && \
    pip3 install scikit-learn==0.21.3 && \
    pip3 install scikit-misc==0.1.1 && \
    pip3 install scipy==1.2.2 && \
    pip3 install statsmodels==0.10.1 && \
    pip3 install tables==3.5.2 && \
    pip3 install torch==1.4.0 && \
    pip3 install anndata==0.6.22.post1 && \
    pip3 install harmony-pytorch==0.1.1 && \
    pip3 install hnswlib==0.3.2.0 && \
    pip3 install fisher==0.1.9 && \
    pip3 install python-igraph==0.7.1.post6 && \
    pip3 install louvain-github==0.6.1.post1 && \
    pip3 install MulticoreTSNE-modified==0.1 && \
    pip3 install umap-learn==0.3.10 && \
    pip3 install scplot==0.0.16 && \
    pip3 install pegasuspy==0.16.11

ENV USER jupyter-user
USER $USER
#we want pip to install into the user's dir when the notebook is running
ENV PIP_USER=true

ENTRYPOINT ["/usr/local/bin/jupyter", "notebook"]	