FROM debian:bullseye-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        python3 \
        python3-pip && \
    ln -s /usr/bin/python3 /usr/bin/python && \
    python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install numpy==1.26.1 && \
    python -m pip install pandas==2.1.1 && \
    python -m pip install XlsxWriter==3.1.9 && \
    apt-get -qq -y remove python3-pip && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN mkdir /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/config/check_uri.py /software
ADD https://raw.githubusercontent.com/lilab-bcb/cumulus/master/docker/monitor_script.sh /software

RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software:$PATH
