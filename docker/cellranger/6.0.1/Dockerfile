FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y alien unzip build-essential dpkg-dev curl gnupg lsb-release procps python3 python3-pip python3-dev && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-sdk=338.0.0-0

RUN pip3 install numpy==1.20.2 && \
    pip3 install pandas==1.2.4 && \
    pip3 install packaging==20.9 && \
    pip3 install stratocumulus==0.2.4

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
ADD cellranger-6.0.1.tar.gz /software

#ADD bcl2fastq2-v2-20-0-linux-x86-64.zip /software/
#RUN unzip -d /software/ /software/bcl2fastq2-v2-20-0-linux-x86-64.zip && alien -i /software/bcl2fastq2-v2.20.0.422-Linux-x86_64.rpm && rm -f /software/bcl2fastq2-v2*

RUN apt-get -qq -y remove alien curl gnupg python3-pip unzip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:/software/cellranger-6.0.1:$PATH
ENV TMPDIR=/tmp
