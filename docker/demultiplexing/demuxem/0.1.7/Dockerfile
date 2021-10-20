FROM debian:buster-slim
SHELL ["/bin/bash", "-c"]

RUN apt-get -qq update && \
    apt-get -qq -y install --no-install-recommends \
        curl \
        unzip \
        gnupg \
        python3 \
        build-essential \
        python3-dev \
        python3-pip

RUN echo "deb [signed-by=/usr/share/keyrings/cloud.google.gpg] http://packages.cloud.google.com/apt cloud-sdk main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
    curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key --keyring /usr/share/keyrings/cloud.google.gpg  add - && \
    apt-get update -y && apt-get install -y google-cloud-sdk=357.0.0-0

RUN curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64-2.2.39.zip" -o "awscliv2.zip" && \
    unzip awscliv2.zip && \
    ./aws/install && \
    rm awscliv2.zip

RUN ln -s /usr/bin/python3 /usr/bin/python

RUN python -m pip install --upgrade pip --no-cache-dir && \
    python -m pip install setuptools==57.0.0 --no-cache-dir && \
    python -m pip install numpy==1.21.0 --no-cache-dir && \
    python -m pip install scipy==1.7.0 --no-cache-dir && \
    python -m pip install pandas==1.2.5 --no-cache-dir && \
    python -m pip install scikit-learn==0.24.2 --no-cache-dir && \
    python -m pip install Cython==0.29.23 --no-cache-dir && \
    python -m pip install h5py==3.3.0 --no-cache-dir && \
    python -m pip install joblib==1.0.1 --no-cache-dir && \
    python -m pip install loompy==3.0.6 --no-cache-dir && \
    python -m pip install matplotlib==3.4.2 --no-cache-dir && \
    python -m pip install natsort==7.1.1 --no-cache-dir && \
    python -m pip install numba==0.53.1 --no-cache-dir && \
    python -m pip install seaborn==0.11.1 --no-cache-dir && \
    python -m pip install anndata==0.7.6 --no-cache-dir && \
    python -m pip install zarr==2.8.3 --no-cache-dir && \
    python -m pip install pegasusio==0.3.1.post2 --no-cache-dir && \
    python -m pip install demuxEM==0.1.7 --no-cache-dir && \
    python -m pip install stratocumulus==0.1.2 --no-cache-dir

RUN apt-get -qq -y remove curl gnupg && \
    apt-get -qq -y autoremove && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH
