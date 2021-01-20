FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        build-essential \
        apt-transport-https \
        ca-certificates \
        gnupg \
        curl \
        wget \
        bzip2 \
        libncurses5-dev \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        python3

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \      
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=324.0.0-0

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN wget https://github.com/samtools/samtools/releases/download/1.11/samtools-1.11.tar.bz2 -O /opt/samtools-1.11.tar.bz2 && \
    cd /opt && tar -xvf samtools-1.11.tar.bz2 && rm samtools-1.11.tar.bz2 && \
    cd /opt/samtools-1.11 && ./configure && make && make install

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/samtools/docker/samtools/process_fastq.py /software

ENV PATH=/software:$PATH
