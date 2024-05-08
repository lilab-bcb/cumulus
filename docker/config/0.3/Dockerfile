FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        python3 \
        python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install numpy==1.22.4 --no-cache-dir && \
    python -m pip install pandas==1.4.2 --no-cache-dir && \
    python -m pip install XlsxWriter==3.0.3 --no-cache-dir && \
    apt-get -qq -y remove python3-pip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/yiming/docker/config/check_uri.py /software
RUN chmod a+rx /software/*
