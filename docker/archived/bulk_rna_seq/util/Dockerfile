FROM continuumio/miniconda3:4.6.14
ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/monitor_script.sh /software/

ADD https://raw.githubusercontent.com/klarman-cell-observatory/cumulus/master/docker/drop-seq/filter_histogram.py /software/


RUN conda install -y pandas \
&& conda install -c plotly plotly \
&& apt-get -qq -y autoremove \
&& apt-get autoclean \
&& rm -rf /var/lib/apt/lists/* \
&& chmod a+rx /software/*
ENV PATH=/software/:${PATH}

ADD ./rsem_agg.py /software/
ADD ./sample_sheet.py /software