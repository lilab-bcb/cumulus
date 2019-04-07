FROM continuumio/miniconda3:latest
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
apt-get install --no-install-recommends -y build-essential dpkg-dev gnupg lsb-release procps unzip zlib1g-dev libncurses5-dev && \
export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
apt-get update && apt-get install -y google-cloud-sdk

RUN pip install --upgrade pip && \
pip install pandas && \
pip install xlsxwriter

RUN mkdir /software
ADD https://github.com/deweylab/RSEM/archive/v1.3.2.tar.gz /software/RSEM_v1.3.2.tar.gz
ADD https://github.com/BenLangmead/bowtie2/releases/download/v2.3.5/bowtie2-2.3.5-linux-x86_64.zip /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software

RUN chmod a+rx /software/monitor_script.sh && \
	unzip /software/bowtie2-2.3.5-linux-x86_64.zip -d /software && \
	tar -xzf /software/RSEM_v1.3.2.tar.gz -C /software && cd /software/RSEM-1.3.2 && make all && cd ../.. 

ENV PATH=/software:/software/bowtie2-2.3.5-linux-x86_64:/software/RSEM-1.3.2:$PATH
ENV TMPDIR=/tmp
