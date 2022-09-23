FROM debian:bullseye-slim

SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y curl unzip ca-certificates python3 python3-pip rsync

RUN curl -o google-cloud-cli.tar.gz https://dl.google.com/dl/cloudsdk/channels/rapid/downloads/google-cloud-cli-403.0.0-linux-x86_64.tar.gz && \
    tar xf google-cloud-cli.tar.gz && ./google-cloud-sdk/install.sh -q --usage-reporting=false

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.2.39.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install pandas==1.2.5 --no-cache-dir && \
    python -m pip install packaging==21.0 --no-cache-dir && \
    python -m pip install stratocumulus==0.1.3 --no-cache-dir


RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
ARG GEOMX_ZIP
ADD $GEOMX_ZIP /software
RUN unzip -d /software/GeoMxNGSPipeline /software/$GEOMX_ZIP && \
    cd /software/GeoMxNGSPipeline/ && tar xf GeoMxNGSPipeline.tgz && mv /software/GeoMxNGSPipeline /var/GeoMxNGSPipeline && rm /software/$GEOMX_ZIP
ADD src /software/scripts
RUN chmod a+rx /software/monitor_script.sh && chmod a+rx /var/GeoMxNGSPipeline/geomxngspipeline
RUN touch /var/GeoMxNGSPipeline/agreement_cli_yes
RUN apt-get -qq -y remove git python3-pip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/*
ENV PATH=/software:/google-cloud-sdk/bin/:/var/GeoMxNGSPipeline:$PATH
