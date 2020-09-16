FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        curl \
        gnupg \
        python3 \
        build-essential \
        python3-dev \
        python3-pip

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=310.0.0-0

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install setuptools==50.3.0 --no-cache-dir && \
    python -m pip install numpy==1.19.2 --no-cache-dir && \
    python -m pip install scipy==1.5.2 --no-cache-dir && \
    python -m pip install pandas==1.1.2 --no-cache-dir && \
    python -m pip install Cython==0.29.21 --no-cache-dir && \
    python -m pip install docopt==0.6.2 --no-cache-dir && \
    python -m pip install h5py==2.10.0 --no-cache-dir && \
    python -m pip install joblib==0.16.0 --no-cache-dir && \
    python -m pip install loompy==3.0.6 --no-cache-dir && \
    python -m pip install matplotlib==3.3.2 --no-cache-dir && \
    python -m pip install natsort==7.0.1 --no-cache-dir && \
    python -m pip install numba==0.51.2  --no-cache-dir && \
    python -m pip install seaborn==0.11.0 --no-cache-dir && \
    python -m pip install anndata==0.7.4 --no-cache-dir && \
    python -m pip install zarr==2.4.0 --no-cache-dir && \
    python -m pip install importlib-metadata==1.7.0 --no-cache-dir && \
    python -m pip install pegasusio==0.2.5.post2 --no-cache-dir && \
    python -m pip install demuxEM==0.1.5 --no-cache-dir

RUN apt-get -qq -y remove curl gnupg && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH
