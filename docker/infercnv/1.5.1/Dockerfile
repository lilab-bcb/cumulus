FROM trinityctat/infercnv:1.5.1

RUN R -e "install.packages(c('optparse', 'logging'))"

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        build-essential \
        automake \
        zlib1g-dev \
        cmake \
        gnupg \
        curl \
        python3 \
        python3-pip \
        python3-dev && \
    echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=312.0.0-0

RUN rm /usr/bin/python && ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install Cython --no-cache-dir && \
    python -m pip install pegasuspy==1.0.0 --no-cache-dir

RUN apt-get -qq -y remove gnupg curl && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
ADD gen_raw_count.py /software
ADD gen_anno.py /software
ADD run_inferCNV.R /software
ADD calc_cnv_metrics.R /software
ADD insert_cnv_scores.py /software

ENV PATH=${PATH}:/software
