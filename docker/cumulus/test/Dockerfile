FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN mkdir -p /usr/share/man/man1 && \
    apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \ 
        build-essential \
        automake \
        zlib1g-dev \
        libxml2-dev \
        cmake \
        gnupg \
        libfftw3-dev \
        default-jdk \
        curl \
        python3 \
        python3-dev \
        python3-pip 

RUN pip3 install setuptools==46.1.3 --no-cache-dir && \
    pip3 install numpy==1.18.2 --no-cache-dir && \
    pip3 install pandas==1.0.3 --no-cache-dir && \
    pip3 install scipy==1.4.1 --no-cache-dir && \
    pip3 install Cython==0.29.16 --no-cache-dir && \
    pip3 install pybind11==2.4.3 --no-cache-dir && \
    pip3 install scikit-learn==0.22.2.post1 --no-cache-dir && \
    pip3 install h5py==2.10.0 --no-cache-dir && \
    pip3 install fitsne==1.1.1 --no-cache-dir && \
    pip3 install joblib==0.14.1 --no-cache-dir && \
    pip3 install python-igraph==0.7.1.post6 --no-cache-dir && \
    pip3 install leidenalg==0.7.0 --no-cache-dir && \
    pip3 install lightgbm==2.2.1 --no-cache-dir && \
    pip3 install llvmlite==0.31.0 --no-cache-dir && \
    pip3 install loompy==3.0.6 --no-cache-dir && \
    pip3 install natsort==7.0.1 --no-cache-dir && \
    pip3 install numba==0.48.0 --no-cache-dir && \
    pip3 install scikit-misc==0.1.1 --no-cache-dir && \
    pip3 install statsmodels==0.11.1 --no-cache-dir && \
    pip3 install tables==3.6.1 --no-cache-dir && \
    pip3 install anndata==0.7.1 --no-cache-dir && \
    pip3 install hnswlib==0.3.4 --no-cache-dir && \
    pip3 install fisher==0.1.9 --no-cache-dir && \
    pip3 install louvain-github==0.6.1.post1 --no-cache-dir && \
    pip3 install MulticoreTSNE-modified==0.1 --no-cache-dir && \
    pip3 install umap-learn==0.3.10 --no-cache-dir && \
    pip3 install torch==1.4.0 --no-cache-dir && \
    pip3 install harmony-pytorch==0.1.3 --no-cache-dir && \
    pip3 install bokeh==2.0.0 --no-cache-dir && \
    pip3 install scplot==0.0.16 --no-cache-dir && \
    pip3 install pegasuspy==0.17.1 --no-cache-dir

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \ 
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=290.0.0-0

RUN apt-get -qq -y remove curl gnupg cmake python3-pip ca-certificates && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    rm /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH
