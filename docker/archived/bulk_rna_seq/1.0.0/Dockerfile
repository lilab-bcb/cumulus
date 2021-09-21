FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

ADD https://github.com/alexdobin/STAR/archive/2.7.3a.zip /software/
ADD https://github.com/deweylab/RSEM/archive/v1.3.3.zip /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl apt-transport-https gnupg lsb-release unzip build-essential zlib1g-dev libncurses5-dev \
    && unzip /software/2.7.3a.zip -d /software/ \
    && unzip /software/v1.3.3.zip -d /software/ && cd /software/RSEM-1.3.3/ && make && cd / \
    && rm -f /software/*.zip \
    && apt-get -qq -y remove ca-certificates curl apt-transport-https gnupg lsb-release python3-pip unzip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/STAR-2.7.3a/bin/Linux_x86_64_static:/software/RSEM-1.3.3:$PATH
