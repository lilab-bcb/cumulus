FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && \
    apt-get install --no-install-recommends -y build-essential && \
    apt-get clean

RUN conda install pandas

RUN conda install -c bioconda popscle

RUN pip install Cython==0.29.16 && \
    pip install pegasusio==0.1.1

RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/souporcell/extract_barcodes_from_rna.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software
RUN chmod a+rx /software/monitor_script.sh

ENV PATH=/software:$PATH