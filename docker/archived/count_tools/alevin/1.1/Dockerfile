FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        build-essential \
        cmake \
        gnupg \
        curl \
        bison \
        git \
        cmake \
        wget \
        zlib1g-dev \
        python3 \
        python3-pip

RUN pip3 install numpy==1.18.1 --no-cache-dir && \
    pip3 install pandas==1.0.3 --no-cache-dir

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \ 
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=290.0.0-0

RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v1.1.0/salmon-1.1.0_linux_x86_64.tar.gz && \
    tar -zxvf salmon-1.1.0_linux_x86_64.tar.gz && \
    rm salmon-1.1.0_linux_x86_64.tar.gz && \
    mkdir /software && \
    mv salmon-latest_linux_x86_64/ /software/salmon

RUN git clone https://github.com/lh3/bioawk.git /software/bioawk && \
    cd /software/bioawk && make

RUN apt-get -qq -y remove curl gnupg cmake python3-pip git wget && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    rm /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:/software/bioawk:/software/salmon/bin:$PATH