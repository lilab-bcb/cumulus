FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install -y build-essential curl unzip zlib1g-dev python3 python3-pip perl git wget && \
    apt-get install --no-install-recommends -y lsb-release gnupg && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-cli=375.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.4.20.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install stratocumulus==0.1.5 --no-cache-dir

RUN mkdir /software && \
    git clone https://github.com/liulab-dfci/TRUST4.git && \
    mv TRUST4 /software/ && \
    cd /software/TRUST4 && \
    make

RUN apt-get -qq -y remove curl unzip python3-pip gnupg && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python
 
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh /software/TRUST4/BuildDatabaseFa.pl /software/TRUST4/BuildImgtAnnot.pl

RUN ln -s /software/TRUST4/BuildDatabaseFa.pl .
RUN ln -s /software/TRUST4/BuildImgtAnnot.pl .

ENV PATH=/software:/software/TRUST4:$PATH
ENV TMPDIR=/tmp
