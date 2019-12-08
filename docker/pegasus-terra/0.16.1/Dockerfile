FROM us.gcr.io/broad-dsp-gcr-public/terra-jupyter-base:0.0.1
USER root
#this makes it so pip runs as root, not the user
ENV PIP_USER=false

RUN apt-get update && apt-get upgrade -y && \
    apt-get install --no-install-recommends -y \
    	build-essential \
    	automake \
    	zlib1g-dev \
    	libxml2-dev \
    	cmake \
    	python3.7-dev \
    	gnupg \
    	lsb-release \
    	libfftw3-dev \
    	libboost-iostreams-dev \
    	default-jdk && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk

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
	pip3 install loompy==2.0.17 && \
	pip3 install natsort==6.0.0 && \
	pip3 install numba==0.45.1 && \
	pip3 install pandas==0.24.2 && \
	pip3 install scikit-learn==0.21.3 && \
	pip3 install scikit-misc==0.1.1 && \
	pip3 install scipy==1.2.2 && \
	pip3 install statsmodels==0.10.1 && \
	pip3 install tables==3.5.2 && \
	pip3 install anndata==0.6.22.post1 && \
	pip3 install hnswlib==0.3.2.0 && \
	pip3 install fisher==0.1.9 && \
	pip3 install python-igraph==0.7.1.post6 && \
	pip3 install louvain-github==0.6.1.post1 && \
	pip3 install MulticoreTSNE-modified==0.1 && \
	pip3 install umap-learn==0.3.10 && \
	pip3 install pegasuspy==0.16.1

ENV USER jupyter-user
USER $USER
#we want pip to install into the user's dir when the notebook is running
ENV PIP_USER=true

ENTRYPOINT ["/usr/local/bin/jupyter", "notebook"]	