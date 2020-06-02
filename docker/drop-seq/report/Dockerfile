FROM continuumio/miniconda3:4.6.14
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/drop-seq/filter_histogram.py /software/

RUN apt-get update && apt-get install -y curl apt-transport-https ca-certificates gnupg && \
echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] https://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg add - && \
apt-get update && apt-get install -y google-cloud-sdk

RUN conda install -y pandas numpy \
&& pip install --upgrade pip \
&& pip install cufflinks \
&& apt-get -qq -y autoremove \
&& apt-get autoclean \
&& rm -rf /var/lib/apt/lists/* \
&& chmod a+rx /software/*
ENV PATH=/software/:${PATH}

ADD ./report.py /software/