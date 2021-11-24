FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y alien unzip build-essential dpkg-dev curl gnupg lsb-release procps python3 python3-pip build-essential libboost-iostreams-dev && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=357.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.2.39.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install pandas==1.2.5 --no-cache-dir && \
    python -m pip install packaging==21.0 --no-cache-dir && \
    python -m pip install stratocumulus==0.1.3 --no-cache-dir

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/scripts/shareseq2bcl /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/scripts/gzip_utils.hpp /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/scripts/barcode_utils.hpp /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/scripts/shareseq_reorg_barcodes.cpp /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/scripts/Makefile /software
RUN cd /software && make clean && make all && cd ..

ADD bcl2fastq2-v2-20-0-linux-x86-64.zip /software/
RUN unzip -d /software/ /software/bcl2fastq2-v2-20-0-linux-x86-64.zip && alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && rm -f /software/bcl2fastq2-v2*

RUN mkdir /indices
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/indices/shareseq_sample_index.csv /indices/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/indices/shareseq_barcode_index.csv /indices/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/shareseq/shareseqdemux/indices/shareseq_flanking_sequence.csv /indices/

RUN apt-get -qq -y remove alien curl gnupg python3-pip unzip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN chmod a+rx /software/monitor_script.sh /software/shareseq2bcl /software/shareseq_reorg_barcodes

ENV PATH=/software:$PATH
ENV TMPDIR=/tmp
