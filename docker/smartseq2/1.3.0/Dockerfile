FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

ADD https://github.com/BenLangmead/bowtie2/releases/download/v2.4.2/bowtie2-2.4.2-linux-x86_64.zip /software/
ADD https://github.com/alexdobin/STAR/archive/2.7.7a.zip /software/
ADD https://github.com/DaehwanKimLab/hisat2/archive/v2.2.1.zip /software/
ADD https://github.com/deweylab/RSEM/archive/v1.3.3.zip /software/
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software/
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/smartseq2/generate_matrix_ss2.py /software/

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl apt-transport-https gnupg lsb-release unzip build-essential zlib1g-dev libncurses5-dev python3 python3-pip \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && python -m pip install --upgrade pip --no-cache-dir \
    && python -m pip install setuptools==53.0.0 --no-cache-dir \
    && python -m pip install numpy==1.20.1 --no-cache-dir \
    && python -m pip install pandas==1.2.2 --no-cache-dir \
    && python -m pip install pegasusio==0.2.10 --no-cache-dir \
    && python -m pip install stratocumulus==0.1.7 --no-cache-dir

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - \
    && apt-get update -y && apt-get install -y google-cloud-cli=478.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.15.62.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN unzip /software/bowtie2-2.4.2-linux-x86_64.zip -d /software/ \
    && unzip /software/2.7.7a.zip -d /software/ \
    && unzip /software/v2.2.1.zip -d /software/ && cd /software/hisat2-2.2.1/ && make && cd / \
    && unzip /software/v1.3.3.zip -d /software/ && cd /software/RSEM-1.3.3/ && make && cd / \
    && rm -f /software/*.zip \
    && apt-get -qq -y remove ca-certificates curl apt-transport-https gnupg lsb-release python3-pip unzip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && chmod a+rx /software/monitor_script.sh \
    && chmod a+rx /software/generate_matrix_ss2.py

ENV PATH=/software:/software/bowtie2-2.4.2-linux-x86_64:/software/STAR-2.7.7a/bin/Linux_x86_64_static:/software/hisat2-2.2.1:/software/RSEM-1.3.3:$PATH
ENV TMPDIR=/tmp

ENV CLOUDSDK_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python
ENV CLOUDSDK_GSUTIL_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python
