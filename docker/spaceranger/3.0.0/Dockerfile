FROM debian:bookworm-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y unzip build-essential dpkg-dev curl gnupg procps python3 python3-pip python3-venv

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=469.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.15.32.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN python3 -m venv /software/python
ENV PATH=/software/python/bin:$PATH

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install pandas==2.2.1 --no-cache-dir && \
    python -m pip install packaging==24.0 --no-cache-dir && \
    python -m pip install stratocumulus==0.2.4 --no-cache-dir

ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software
ADD spaceranger-3.0.0.tar.gz /software

RUN apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:/software/spaceranger-3.0.0:/usr/local/aws-cli/v2/current/bin:$PATH
ENV TMPDIR=/tmp
