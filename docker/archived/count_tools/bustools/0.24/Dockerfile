FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y build-essential cmake lsb-release gnupg && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    apt-get clean

RUN pip install --upgrade pip && \
    pip install numpy==1.18.1 && \
    pip install pandas==0.25.3 && \
    pip install anndata==0.7.1 && \
    pip install scipy==1.4.1 && \
    pip install numba==0.47.0 && \
    pip install loompy==3.0.6 && \
    pip install h5py==2.10.0 && \
    pip install natsort==7.0.0 && \
    pip install kb-python==0.24.4

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