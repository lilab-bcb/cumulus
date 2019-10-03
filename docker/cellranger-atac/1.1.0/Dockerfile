FROM debian:stretch-slim
SHELL ["/bin/bash", "-c"]

#ADD bcl2fastq2_0v2.20.0.422-2_amd64.deb /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/
ADD cellranger-atac-1.1.0.tar.gz /software/

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl apt-transport-https gnupg lsb-release python3 python3-pip \
  	&& pip3 --no-cache-dir install pandas \
	&& export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" \
	&& echo "deb https://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list \
	&& curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - \
    && apt-get update && apt-get -qq -y install --no-install-recommends google-cloud-sdk \
    && apt-get -qq -y remove ca-certificates curl apt-transport-https gnupg lsb-release python3-pip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && ln -s /usr/bin/python3 /usr/bin/python \
    #&& dpkg -i /software/bcl2fastq2_0v2.20.0.422-2_amd64.deb && rm /software/bcl2fastq2_0v2.20.0.422-2_amd64.deb \
	&& chmod a+rx /software/*


ENV PATH=/software:/software/cellranger-atac-1.1.0/:$PATH
