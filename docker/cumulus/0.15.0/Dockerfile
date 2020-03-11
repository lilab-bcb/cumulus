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
		default-jdk && \
	export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
	echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
	curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
	apt-get update && apt-get install -y google-cloud-sdk


RUN pip install --upgrade pip && \
	pip install numpy==1.16.4 && \
	pip install Cython==0.29.13 && \
	pip install pybind11==2.3.0 && \
	pip install h5py==2.10.0 && \
	pip install fitsne==1.0.1 && \
	pip install joblib==0.13.2 && \
	pip install leidenalg==0.7.0 && \
	pip install lightgbm==2.2.1 && \
	pip install llvmlite==0.31.0 && \
	pip install loompy==3.0.6 && \
	pip install natsort==6.0.0 && \
	pip install numba==0.47.0 && \
	pip install pandas==1.0.1 && \
	pip install scikit-learn==0.21.3 && \
	pip install scikit-misc==0.1.1 && \
	pip install scipy==1.2.2 && \
	pip install statsmodels==0.11.1 && \
	pip install tables==3.5.2 && \
	pip install anndata==0.7.1 && \
	pip install hnswlib==0.3.4 && \
	pip install fisher==0.1.9 && \
	pip install python-igraph==0.7.1.post6 && \
	pip install louvain-github==0.6.1.post1 && \
	pip install MulticoreTSNE-modified==0.1 && \
	pip install umap-learn==0.3.10 && \
	pip install torch==1.4.0 && \
	pip install harmony-pytorch==0.1.2 && \
	pip install pegasuspy==0.17.0


RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH

# Install gcloud 264.0.0
RUN gsutil -m cp gs://cloud-sdk-release/google-cloud-sdk-264.0.0-linux-x86_64.tar.gz /software && \
	cd /software && \
	tar -zxvf google-cloud-sdk-264.0.0-linux-x86_64.tar.gz && \
	rm google-cloud-sdk-264.0.0-linux-x86_64.tar.gz && \
	cd google-cloud-sdk && \
	./install.sh

ENV PATH=/software/google-cloud-sdk/bin:$PATH
