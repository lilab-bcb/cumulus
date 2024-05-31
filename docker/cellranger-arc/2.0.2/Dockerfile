FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y unzip rsync build-essential dpkg-dev curl gnupg procps python3 python3-pip && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=392.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.7.11.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN pip3 install --upgrade pip && \
    pip3 install pandas==1.4.3 && \
    pip3 install packaging==21.3 && \
    pip3 install stratocumulus==0.2.4

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
ADD cellranger-arc-2.0.2.tar.gz /software

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    rm -f /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:/software/cellranger-arc-2.0.2:/usr/local/aws-cli/v2/current/bin:$PATH
ENV TMPDIR=/tmp
