FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y wget build-essential curl unzip wget zlib1g-dev python3 python3-pip && \
    apt-get install --no-install-recommends -y lsb-release gnupg && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=357.0.0-0 && \
    apt-get -qq -y autoremove && \
    apt-get clean

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.2.39.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install stratocumulus==0.1.3 --no-cache-dir

RUN mkdir /software && \
    wget https://github.com/haowenz/chromap/archive/refs/tags/v0.1.3.zip && \
    unzip v0.1.3.zip && \
    rm v0.1.3.zip && \
    mv chromap-0.1.3 /software/ && \
    cd /software/chromap-0.1.3 && \
    make

RUN apt-get -qq -y remove curl gnupg unzip wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log
 
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/chromap-0.1.3:$PATH
ENV TMPDIR=/tmp
