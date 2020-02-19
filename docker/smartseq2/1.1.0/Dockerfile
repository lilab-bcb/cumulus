FROM debian:stretch-slim
SHELL ["/bin/bash", "-c"]

ADD https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5.1/bowtie2-2.3.5.1-linux-x86_64.zip /software/
ADD https://github.com/alexdobin/STAR/archive/2.7.3a.zip /software/
ADD https://github.com/DaehwanKimLab/hisat2/archive/v2.1.0.zip /software/
ADD https://github.com/deweylab/RSEM/archive/v1.3.3.zip /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl apt-transport-https gnupg lsb-release unzip build-essential zlib1g-dev libncurses5-dev python3 python3-pip \ 
    && pip3 --no-cache-dir install numpy pandas \
    && export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" \
    && echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
    && curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    && apt-get update && apt-get -qq -y install --no-install-recommends google-cloud-sdk \
    && unzip /software/bowtie2-2.3.5.1-linux-x86_64.zip -d /software/ \
    && unzip /software/2.7.3a.zip -d /software/ \
    && unzip /software/v2.1.0.zip -d /software/ && cd /software/hisat2-2.1.0/ && make && cd / \
    && unzip /software/v1.3.3.zip -d /software/ && cd /software/RSEM-1.3.3/ && make && cd / \
    && rm -f /software/*.zip \
    && apt-get -qq -y remove ca-certificates curl apt-transport-https gnupg lsb-release python3-pip unzip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/bowtie2-2.3.5.1-linux-x86_64:/software/STAR-2.7.3a/bin/Linux_x86_64_static:/software/hisat2-2.1.0:/software/RSEM-1.3.3:$PATH
