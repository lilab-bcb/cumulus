FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y alien unzip build-essential dpkg-dev curl gnupg procps python3 python3-pip

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=415.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.9.18.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install pandas==1.5.3 --no-cache-dir && \
    python -m pip install packaging==23.0 --no-cache-dir && \
    python -m pip install stratocumulus==0.1.7 --no-cache-dir

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
ADD spaceranger-2.0.1.tar.gz /software

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:/software/spaceranger-2.0.1:$PATH
ENV TMPDIR=/tmp
