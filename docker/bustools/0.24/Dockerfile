FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y build-essential cmake && \
    apt-get clean

RUN pip install --upgrade pip && \
    pip install numpy==1.18.1 && \
    pip install pandas==0.25.3 && \
    pip install anndata==0.7.1 && \
    pip install scipy==1.4.1 && \
    pip install numba==0.47.0 && \
    pip install loompy==3.0.6 && \
    pip install h5py==2.10.0 && \
    pip install natsort==7.0.0 && \
    pip install kb-python==0.24.4
