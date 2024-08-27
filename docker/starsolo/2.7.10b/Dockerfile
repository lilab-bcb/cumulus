FROM debian:bullseye-slim
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
        python3-pip

RUN pip3 install --upgrade pip --no-cache-dir && \
    pip3 install numpy==1.24.2 --no-cache-dir && \
    pip3 install pandas==2.0.0 --no-cache-dir && \
    pip3 install h5py==3.8.0 --no-cache-dir && \
    pip3 install pegasusio==0.8.0 --no-cache-dir && \
    pip3 install stratocumulus==0.2.4 --no-cache-dir

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=483.0.0-0

ENV CLOUDSDK_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python
ENV CLOUDSDK_GSUTIL_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.11.13.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN wget https://github.com/alexdobin/STAR/archive/refs/tags/STAR_2.7.10b_alpha_230301.tar.gz && \
    tar -xzvf STAR_2.7.10b_alpha_230301.tar.gz && \
    rm STAR_2.7.10b_alpha_230301.tar.gz && \
    mkdir /software && \
    mv STAR-STAR_2.7.10b_alpha_230301 /software/STAR && \
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
