FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y wget build-essential curl unzip wget zlib1g-dev python3 python3-pip gnupg && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update && apt-get install -y google-cloud-cli=392.0.0-0 && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.7.11.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install stratocumulus==0.1.7 --no-cache-dir

RUN mkdir /software && \
    wget https://github.com/haowenz/chromap/archive/refs/tags/v0.1.5.zip && \
    unzip v0.1.5.zip && \
    rm v0.1.5.zip && \
    mv chromap-0.1.5 /software/ && \
    cd /software/chromap-0.1.5 && \
    make

RUN apt-get -qq -y remove curl gnupg unzip wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/chromap-0.1.5:$PATH
ENV TMPDIR=/tmp
