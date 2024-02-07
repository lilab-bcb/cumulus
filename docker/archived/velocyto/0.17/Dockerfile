FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install --no-install-recommends -y \
		build-essential zlib1g-dev libbz2-dev
RUN apt-get install --no-install-recommends -y libcurl4-openssl-dev libssl-dev
RUN conda install -c bioconda samtools
RUN conda install -y numpy scipy cython numba matplotlib scikit-learn h5py click
RUN pip install velocyto
RUN mkdir /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/
RUN chmod a+rx /software/monitor_script.sh
ENV PATH=/software/:$PATH
