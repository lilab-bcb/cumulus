FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y alien unzip dpkg-dev curl gnupg lsb-release procps python3 python3-pip build-essential libboost-iostreams-dev pigz && \
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
    python -m pip install stratocumulus==0.1.7 --no-cache-dir

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/shareseq/shareseq_reorg/scripts/gzip_utils.hpp /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/shareseq/shareseq_reorg/scripts/barcode_utils.hpp /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/shareseq/shareseq_reorg/scripts/shareseq_reorg_barcodes.cpp /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/shareseq/shareseq_reorg/scripts/Makefile /software
RUN cd /software && make clean && make all && cd ..

RUN mkdir /indices
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/shareseq/shareseq_reorg/indices/shareseq_barcode_index.csv /indices/
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/shareseq/shareseq_reorg/indices/shareseq_flanking_sequence.csv /indices/

RUN apt-get -qq -y remove alien curl gnupg python3-pip unzip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN chmod a+rx /software/monitor_script.sh /software/shareseq_reorg_barcodes

ENV PATH=/software:$PATH
ENV TMPDIR=/tmp
