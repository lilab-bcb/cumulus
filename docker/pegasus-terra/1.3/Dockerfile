FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:0.0.19
USER root
#this makes it so pip runs as root, not the user
ENV PIP_USER=false

RUN apt-get update && apt-get install -yq --no-install-recommends \
        build-essential \
        python3-dev \
        libfftw3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN pip3 -V && \
    pip3 install --upgrade pip && \
    pip3 install setuptools==53.0.0 && \
    pip3 install numpy==1.19.5 && \
    pip3 install pandas==1.2.1 && \
    pip3 install scipy==1.5.4 && \
    pip3 install Cython==0.29.21 && \
    pip3 install pybind11==2.6.1 && \
    pip3 install scikit-image==0.18.1 && \
    pip3 install scikit-learn==0.24.1 && \
    pip3 install h5py==3.1.0 && \
    pip3 install fitsne==1.1.1 && \
    pip3 install importlib-metadata==3.4.0 && \
    pip3 install joblib==1.0.0 && \
    pip3 install psutil==5.8.0 && \
    pip3 install threadpoolctl==2.1.0 && \
    pip3 install python-igraph==0.8.3 && \
    pip3 install leidenalg==0.8.3 && \
    pip3 install lightgbm==3.1.1 && \
    pip3 install loompy==3.0.6 && \
    pip3 install matplotlib==3.3.4 && \
    pip3 install natsort==7.1.1 && \
    pip3 install numba==0.52.0 && \
    pip3 install scanorama==1.7 && \
    pip3 install scikit-misc==0.1.3 && \
    pip3 install seaborn==0.11.1 && \
    pip3 install statsmodels==0.12.2 && \
    pip3 install numcodecs==0.7.3 && \
    pip3 install networkx==2.5 && \
    pip3 install zarr==2.6.1 && \
    pip3 install anndata==0.7.5 && \
    pip3 install hnswlib==0.5.0 && \
    pip3 install louvain==0.7.0 && \
    pip3 install umap-learn==0.4.6 && \
    pip3 install torch==1.7.1 && \
    pip3 install harmony-pytorch==0.1.6 && \
    pip3 install cirrocumulus==1.1.15.post4 && \
    pip3 install annoy==1.17.0 && \
    pip3 install pegasusio==0.2.10 && \
    pip3 install demuxEM==0.1.5.post1 && \
    pip3 install forceatlas2-python==1.1 && \
    pip3 install pegasuspy==1.3.0

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
