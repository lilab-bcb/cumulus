FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

ADD https://github.com/lilab-bcb/cumulus_feature_barcoding/archive/refs/tags/0.9.0.tar.gz /software/
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software/

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        curl \
        unzip \
        gnupg \
        python3 \
        python3-dev \
        python3-pip \
        build-essential \
        libisal2 libisal-dev \
        libdeflate0 libdeflate-dev

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=389.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.7.7.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install numpy==1.22.4 --no-cache-dir && \
    python -m pip install pandas==1.4.2 --no-cache-dir && \
    python -m pip install matplotlib==3.5.2 --no-cache-dir && \
    python -m pip install stratocumulus==0.1.7 --no-cache-dir

RUN tar -xzf /software/0.9.0.tar.gz -C /software && \
    cd /software/cumulus_feature_barcoding-0.9.0 && make clean && make all && cd ../.. && \
    rm -f /software/0.9.0.tar.gz

RUN apt-get -qq -y remove gnupg && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN chmod a+rx /software/*

ENV PATH=/software:/software/cumulus_feature_barcoding-0.9.0:$PATH
