FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:1.0.4
USER root
#this makes it so pip runs as root, not the user
ENV PIP_USER=false

RUN apt-get update && apt-get install -yq --no-install-recommends \
        python3-dev && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN conda install -y -n base -c conda-forge mamba && \
    conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    mamba install -y -c conda-forge \
        pyfit-sne=1.2.1 \
        hnswlib=0.6.2 \
        loompy=3.0.6 \
        louvain=0.7.1 && \
    mamba install -y -c bioconda \
        scanorama=1.7.1 \
        bioconductor-deseq2=1.34.0 \
        bioconductor-fgsea=1.20.0

# Uninstall TensorFlow packages
RUN pip3 uninstall -y \
        tfx-bsl \
        tensorboard \
        tensorboard-data-server \
        tensorboard-plugin-wit \
        tensorflow \
        tensorflow-cloud \
        tensorflow-datasets \
        tensorflow-estimator \
        tensorflow-hub \
        tensorflow-io-gcs-filesystem \
        tensorflow-metadata \
        tensorflow-probability \
        tensorflow-serving-api \
        tensorflow-transform \
        tensorflow-io

RUN pip3 -V && \
    pip3 install numpy==1.21.6 && \
    pip3 install pandas==1.3.5 && \
    pip3 install scipy==1.7.3 && \
    pip3 install Cython==0.29.28 && \
    pip3 install pybind11==2.9.2 && \
    pip3 install scikit-image==0.19.2 && \
    pip3 install scikit-learn==1.0.2 && \
    pip3 install h5py==3.6.0 && \
    pip3 install importlib-metadata==4.11.3 && \
    pip3 install joblib==1.1.0 && \
    pip3 install psutil==5.9.0 && \
    pip3 install threadpoolctl==3.1.0 && \
    pip3 install python-igraph==0.9.10 && \
    pip3 install leidenalg==0.8.9 && \
    pip3 install lightgbm==3.3.2 && \
    pip3 install matplotlib==3.5.1 && \
    pip3 install natsort==8.1.0 && \
    pip3 install llvmlite --ignore-installed && \
    pip3 install numba==0.55.1 && \
    pip3 install scikit-misc==0.1.4 && \
    pip3 install seaborn==0.11.2 && \
    pip3 install statsmodels==0.13.2 && \
    pip3 install networkx==2.6.3 && \
    pip3 install zarr==2.11.3 && \
    pip3 install anndata==0.8.0 && \
    pip3 install umap-learn==0.5.3 && \
    pip3 install torch==1.11.0 && \
    pip3 install harmony-pytorch==0.1.7 && \
    pip3 install cirrocumulus==1.1.29 && \
    pip3 install pegasusio==0.5.1.post1 && \
    pip3 install demuxEM==0.1.7 && \
    pip3 install forceatlas2-python==1.1 && \
    pip3 install nmf-torch==0.1.1 && \
    pip3 install pegasuspy==1.6.0 && \
    pip3 install rpy2==3.5.1 && \
    pip3 install scvi-tools==0.16.0

RUN conda clean -y -a

RUN wget https://bin.equinox.io/c/4VmDzA7iaHb/ngrok-stable-linux-amd64.zip && \
    unzip ngrok-stable-linux-amd64.zip && \
    rm ngrok-stable-linux-amd64.zip && \
    mkdir -p /software && \
    mv ngrok /software/

ENV PATH=/software:$PATH

ENV USER jupyter
USER $USER
#we want pip to install into the user's dir when the notebook is running
ENV PIP_USER=true

ENTRYPOINT ["/opt/conda/bin/jupyter", "notebook"]
