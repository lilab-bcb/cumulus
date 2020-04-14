FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends curl lsb-release curl gnupg python3 python3-pip && \
    pip3 install numpy==1.18.2 --no-cache-dir && \
    pip3 install pandas==1.0.3 --no-cache-dir && \
    export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
    echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
    apt-get update && apt-get -qq -y install --no-install-recommends google-cloud-sdk=264.0.0-0 && \
    apt-get -qq -y remove curl gnupg lsb-release python3-pip ca-certificates && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log && \
    ln -s /usr/bin/python3 /usr/bin/python 

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software

RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH