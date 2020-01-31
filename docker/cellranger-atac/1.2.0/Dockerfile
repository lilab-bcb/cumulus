FROM debian:stretch-slim
SHELL ["/bin/bash", "-c"]

#ADD bcl2fastq2-v2-20-0-linux-x86-64.zip /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/
ADD cellranger-atac-1.2.0.tar.gz /software/

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl apt-transport-https gnupg lsb-release python3 python3-pip \
  	&& pip3 --no-cache-dir install pandas \
	&& export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" \
	&& echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
	&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    && apt-get update && apt-get -qq -y install --no-install-recommends google-cloud-sdk \
#    && apt-get install --no-install-recommends -y alien unzip \
#    && unzip -d /software/ /software/bcl2fastq2-v2-20-0-linux-x86-64.zip && alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && rm /software/bcl2fastq2-v2* \
#    && apt-get -qq -y remove alien unzip \
    && apt-get -qq -y remove ca-certificates curl apt-transport-https gnupg lsb-release python3-pip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && ln -s /usr/bin/python3 /usr/bin/python \
	&& chmod a+rx /software/*

ENV PATH=/software:/software/cellranger-atac-1.2.0/:$PATH
