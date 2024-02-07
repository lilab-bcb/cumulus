FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

ADD https://github.com/klarman-cell-observatory/cumulus_feature_barcoding/archive/0.6.0.tar.gz /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        curl \
        unzip \
        gnupg \
        python3 \
        python3-dev \
        python3-pip \
        build-essential \
        libboost-iostreams-dev && \
    pip3 install --upgrade pip && \
    pip3 install numpy==1.21.0 --no-cache-dir && \
    pip3 install pandas==1.2.5 --no-cache-dir && \
    pip3 install matplotlib==3.4.2 --no-cache-dir && \
    pip3 install stratocumulus==0.1.2 --no-cache-dir && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=357.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.2.39.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN tar -xzf /software/0.6.0.tar.gz -C /software && \
    cd /software/cumulus_feature_barcoding-0.6.0 && make clean && make all && cd ../.. && \
    rm -f /software/0.6.0.tar.gz && \
    apt-get -qq -y remove gnupg && \
    apt-get -qq -y autoremove && \
    apt-get autoclean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    chmod a+rx /software/*

ENV PATH=/software:/software/cumulus_feature_barcoding-0.6.0:$PATH
