FROM continuumio/miniconda3:4.6.14	
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install --no-install-recommends -y \
		build-essential \
		automake \
		zlib1g-dev \
		libxml2-dev \
		cmake \
		gnupg \
		lsb-release \
		libfftw3-dev \
		libboost-iostreams-dev \
		default-jdk && \
	export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
	echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
	curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
	apt-get update && apt-get install -y google-cloud-sdk

RUN pip install --upgrade pip && \
	pip install numpy==1.16.4 && \
	pip install Cython==0.29.13 && \
	pip install pybind11==2.3.0 && \
	pip install h5py==2.9.0 && \
	pip install joblib==0.13.2 && \
	pip install leidenalg==0.7.0 && \
	pip install lightgbm==2.2.1 && \
	pip install loompy==2.0.17 && \
	pip install natsort==6.0.0 && \
	pip install numba==0.45.1 && \
	pip install pandas==0.24.2 && \
	pip install scikit-learn==0.21.3 && \
	pip install scikit-misc==0.1.1 && \
	pip install scipy==1.2.2 && \
	pip install statsmodels==0.10.1 && \
	pip install tables==3.5.2 && \
	pip install anndata==0.6.22.post1 && \
	pip install hnswlib==0.3.2.0 && \
	pip install fisher-modified==0.1.5 && \
	pip install python-igraph==0.7.1.post6 && \
	pip install louvain-github==0.6.1.post1 && \
	pip install MulticoreTSNE-modified==0.1 && \
	pip install umap-learn==0.3.10 && \
	pip install sccloud[fitsne]==0.14.0rc3

COPY extract_ADT_UMIs /software/extract_ADT_UMIs
ADD https://raw.githubusercontent.com/klarman-cell-observatory/scCloud/master/docker/monitor_script.sh /software

RUN cd /software/extract_ADT_UMIs && make clean && make all && cd ../.. && \
	chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/extract_ADT_UMIs:$PATH
