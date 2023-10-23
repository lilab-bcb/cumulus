FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        ca-certificates \
        rsync \
        curl \
        unzip \
        gnupg \
        python3 \
        python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install numpy==1.26.1 && \
    python -m pip install pandas==2.1.1 && \
    python -m pip install XlsxWriter==3.1.9 && \
    python -m pip install boto3==1.28.68 && \
    python -m pip install stratocumulus==0.2.1

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
    apt-get update -y && apt-get install -y google-cloud-cli=450.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.13.26.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/config/check_uri.py /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:$PATH
