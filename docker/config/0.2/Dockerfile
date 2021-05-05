FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        python3 \
        python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install numpy==1.19.5 --no-cache-dir && \
    python -m pip install pandas==1.2.0 --no-cache-dir && \
    python -m pip install XlsxWriter==1.3.7 --no-cache-dir && \
    apt-get -qq -y remove python3-pip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log
