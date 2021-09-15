FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

ADD https://github.com/klarman-cell-observatory/cumulus_feature_barcoding/archive/0.6.0.tar.gz /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        curl \
        gnupg \
        python3 \
        python3-dev \
        python3-pip \
        build-essential \
        libboost-iostreams-dev && \
    pip3 --no-cache-dir install numpy pandas matplotlib && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=340.0.0-0

RUN tar -xzf /software/0.6.0.tar.gz -C /software && \
    cd /software/cumulus_feature_barcoding-0.6.0 && make clean && make all && cd ../.. && \
    rm -f /software/0.6.0.tar.gz && \
    apt-get -qq -y remove gnupg && \
    apt-get -qq -y autoremove && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    chmod a+rx /software/*

RUN python -m pip install stratocumulus==0.1 --no-cache-dir 

ENV PATH=/software:/software/cumulus_feature_barcoding-0.6.0:$PATH
