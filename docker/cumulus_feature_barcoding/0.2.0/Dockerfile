FROM debian:stretch-slim
SHELL ["/bin/bash", "-c"]

ADD https://github.com/klarman-cell-observatory/cumulus_feature_barcoding/archive/0.2.0.tar.gz /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl tar\ 
        apt-transport-https gnupg lsb-release python3 python3-pip \
        build-essential zlib1g-dev libboost-iostreams-dev \
    && pip3 --no-cache-dir install numpy pandas matplotlib \
    && export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" \
    && echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    && apt-get update && apt-get -qq -y install --no-install-recommends google-cloud-sdk \
    && tar -xzf /software/0.2.0.tar.gz -C /software \
    && cd /software/cumulus_feature_barcoding-0.2.0 && make clean && make all && cd ../.. \
    && rm -f /software/0.2.0.tar.gz \ 
    && apt-get -qq -y remove ca-certificates curl apt-transport-https gnupg lsb-release python3-pip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && chmod a+rx /software/*

ENV PATH=/software:/software/cumulus_feature_barcoding-0.2.0:$PATH
