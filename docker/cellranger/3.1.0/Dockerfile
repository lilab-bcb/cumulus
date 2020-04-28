FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y build-essential dpkg-dev curl gnupg lsb-release procps python3 python3-pip && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=290.0.0-0

RUN pip3 install pandas && \
    pip3 install xlsxwriter

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
ADD cellranger-3.1.0.tar.gz /software

# ADD bcl2fastq2_0v2.20.0.422-2_amd64.deb /software
# RUN dpkg -i /software/bcl2fastq2_0v2.20.0.422-2_amd64.deb && rm /software/bcl2fastq2_0v2.20.0.422-2_amd64.deb

RUN apt-get -qq -y remove curl gnupg python3-pip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    rm /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:/software/cellranger-3.1.0:$PATH
ENV TMPDIR=/tmp