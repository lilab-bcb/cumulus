FROM continuumio/miniconda3:latest	
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install --no-install-recommends -y \
		build-essential \
		automake \
		zlib1g-dev \
		python3-igraph \
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
	pip install numpy && \
	pip install Cython && \
	pip install pybind11 && \
	mkdir /software && \
	git clone https://github.com/nmslib/hnsw.git /software/hnswlib && \
	cd /software/hnswlib/python_bindings && python setup.py install && cd ../../..

COPY scCloud-0.6.0.tar.gz /software
COPY extract_ADT_UMIs /software/extract_ADT_UMIs
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software

RUN tar -xzf /software/scCloud-0.6.0.tar.gz -C /software && \
	cd /software/scCloud-0.6.0 && pip install . && cd ../.. && \
	cd /software/extract_ADT_UMIs && make all && cd ../.. && \
	chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/extract_ADT_UMIs:$PATH
