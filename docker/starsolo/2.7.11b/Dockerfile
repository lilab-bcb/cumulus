FROM debian:bookworm-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        wget \
        build-essential \
        curl \
        zlib1g-dev \
        unzip \
        gnupg \
        python3 \
        python3-pip \
        python3-venv \
        xxd

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=502.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.22.4.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN python3 -m venv /software/python
ENV PATH=/software/python/bin:$PATH

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install numpy==2.1.3 --no-cache-dir && \
    python -m pip install pandas==2.2.3 --no-cache-dir && \
    python -m pip install h5py==3.12.1 --no-cache-dir && \
    python -m pip install pegasusio==0.9.1 --no-cache-dir && \
    python -m pip install stratocumulus==0.2.4 --no-cache-dir

ENV CLOUDSDK_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python
ENV CLOUDSDK_GSUTIL_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python

RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/2.7.11b.tar.gz && \
    tar -xzvf 2.7.11b.tar.gz && \
    rm 2.7.11b.tar.gz && \
    mv STAR-2.7.11b /software/STAR && \
    cd /software/STAR/source && \
    make STAR && \
    mkdir -p /software/STAR/bin && \
    ln -s /software/STAR/source/STAR /software/STAR/bin/STAR

RUN apt-get -qq -y remove curl gnupg wget curl && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    ln -s /usr/bin/python3 /usr/bin/python

ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/STAR/bin:$PATH
ENV TMPDIR=/tmp
