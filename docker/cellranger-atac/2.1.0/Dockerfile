FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y unzip build-essential dpkg-dev curl gnupg procps python3 python3-pip && \
    echo "deb https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get install -y google-cloud-cli=510.0.0-0

ENV CLOUDSDK_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python
ENV CLOUDSDK_GSUTIL_PYTHON=/usr/lib/google-cloud-sdk/platform/bundledpythonunix/bin/python

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.24.6.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip && \
    python -m pip install pandas==2.2.3 && \
    python -m pip install packaging==24.2 && \
    python -m pip install stratocumulus==0.3.0

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
ADD cellranger-atac-2.1.0.tar.gz /software

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:/software/cellranger-atac-2.1.0:/usr/local/aws-cli/v2/current/bin:$PATH
ENV TMPDIR=/tmp
