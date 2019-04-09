FROM continuumio/miniconda3:latest
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/drop-seq/summary.py /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software

RUN apt-get update && \
apt-get install --no-install-recommends -y curl dpkg-dev gnupg lsb-release procps && \
export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
apt-get update && apt-get install -y google-cloud-sdk

RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends default-jre ca-certificates curl r-base libgomp1 build-essential zlib1g-dev \
&& conda install -y pandas numpy \
&& pip install --upgrade pip \
&& pip install cufflinks plotly \
&& curl -L -s -o "/Drop-seq_tools-2.2.0.zip" "https://github.com/broadinstitute/Drop-seq/releases/download/v2.2.0/Drop-seq_tools-2.2.0.zip" \
&& unzip /Drop-seq_tools-2.2.0.zip \
&& rm /Drop-seq_tools-2.2.0.zip \
&& mv /Drop-seq_tools-2.2.0 /software/Drop-seq_tools \
&& echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile \
&& Rscript -e "install.packages(c('ineq', 'optparse'))" \
&& curl -L -s -o /software/v1.0.4.zip https://github.com/broadinstitute/DropSeq_utility/archive/v1.0.4.zip \
&& unzip -d /software/ /software/v1.0.4.zip \
&& rm /software/v1.0.4.zip \
&& curl -L -s -o "/2.7.0c.tar.gz" "https://github.com/alexdobin/STAR/archive/2.7.0c.tar.gz" \
&& tar xf /2.7.0c.tar.gz \
&& rm -rf /STAR-2.7.0c/bin/Linux_x86_64_static /STAR-2.7.0c/bin/MacOSX_x86_64 \
&& mv /STAR-2.7.0c /software/ \
&& rm -f /2.7.0c.tar.gz \
&& curl -L -s -o "/software/picard.jar" "https://github.com/broadinstitute/picard/releases/download/2.19.0/picard.jar" \
&& curl -L -s -o "/software/filter_histogram.py" "https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/drop-seq/filter_histogram.py" \
&& apt-get -qq -y remove curl build-essential \
&& apt-get -qq -y autoremove \
&& apt-get autoclean \
&& rm -rf /var/lib/apt/lists/* \
&& chmod a+rx /software/* \
&& chmod a+rx /software/DropSeq_utility-1.0.4/* \
&& chmod a+rx /software/Drop-seq_tools/*
ENV PATH=/software/bin/:/software/:/software/STAR-2.7.0c/bin/Linux_x86_64/:/software/DropSeq_utility-1.0.4/:/software/Drop-seq_tools/:${PATH}



