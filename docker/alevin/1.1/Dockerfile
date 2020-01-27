FROM ubuntu:18.04
SHELL ["/bin/bash", "-c"]

RUN apt-get update -qq && \
    apt-get install -y wget git cmake && \
    apt-get install --no-install-recommends -y build-essential cmake gnupg lsb-release curl bison zlib1g-dev && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk && \
    apt-get clean

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz && \
    tar -zxvf salmon-1.1.0_linux_x86_64.tar.gz && \
    rm salmon-1.1.0_linux_x86_64.tar.gz && \
    mkdir /software && \
    mv salmon-latest_linux_x86_64/ /software/salmon

RUN git clone https://github.com/lh3/bioawk.git /software/bioawk && \
    cd /software/bioawk && make

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
ADD make_transcripts_fa.sh /software
ADD make_tgmap.sh /software

RUN chmod a+rx /software/monitor_script.sh /software/make_transcripts_fa.sh /software/make_tgmap.sh

# Install gcloud 264.0.0
RUN gsutil -m cp gs://cloud-sdk-release/google-cloud-sdk-264.0.0-linux-x86_64.tar.gz /software && \
    cd /software && \
    tar -zxvf google-cloud-sdk-264.0.0-linux-x86_64.tar.gz && \
    rm google-cloud-sdk-264.0.0-linux-x86_64.tar.gz && \
    cd google-cloud-sdk && \
    ./install.sh

ENV PATH=/software/google-cloud-sdk/bin:/software:/software/bioawk:/software/salmon/bin:$PATH