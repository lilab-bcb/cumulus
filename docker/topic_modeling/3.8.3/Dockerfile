FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/
ADD ./lda.py /software

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl apt-transport-https gnupg lsb-release build-essential python3 python3-dev python3-pip \
    && pip3 --no-cache-dir install setuptools Cython numba \
    && pip3 --no-cache-dir install gensim pyLDAvis matplotlib pegasusio \
    && apt-get -qq -y remove ca-certificates curl apt-transport-https gnupg lsb-release python3-pip unzip \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* /var/log/dpkg.log \
    && ln -s /usr/bin/python3 /usr/bin/python \
    && chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH
