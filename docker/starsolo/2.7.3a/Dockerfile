FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y wget build-essential lsb-release curl gnupg && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    apt-get clean

RUN pip install --upgrade pip && \
    pip install numpy==1.18.1 && \
    pip install pandas==0.25.3

RUN wget https://github.com/alexdobin/STAR/archive/2.7.3a.tar.gz && \
    tar -xzvf 2.7.3a.tar.gz && \
    rm 2.7.3a.tar.gz && \
    mkdir /software && \
    mv STAR-2.7.3a/ /software/STAR

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software

RUN chmod a+rx /software/monitor_script.sh

# Install gcloud 264.0.0
RUN gsutil -m cp gs://cloud-sdk-release/google-cloud-sdk-264.0.0-linux-x86_64.tar.gz /software && \
    cd /software && \
    tar -zxvf google-cloud-sdk-264.0.0-linux-x86_64.tar.gz && \
    rm google-cloud-sdk-264.0.0-linux-x86_64.tar.gz && \
    cd google-cloud-sdk && \
    ./install.sh

ENV PATH=/software/google-cloud-sdk/bin:/software:/software/STAR/bin/Linux_x86_64:$PATH