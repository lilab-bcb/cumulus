FROM rocker/r-ver:3.4

RUN apt-get update && \
apt-get install --no-install-recommends -y curl dpkg-dev gnupg lsb-release procps && \
export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
apt-get update && apt-get install -y google-cloud-sdk


RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  less \
  libbamtools-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libz-dev \
  nano \
  pandoc

RUN \
  R -e 'chooseCRANmirror(ind=52); install.packages(c("devtools", "Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown", "withr"))'

RUN mkdir -p /software
ENTRYPOINT ["/bin/bash"]
WORKDIR "/software"

RUN \
  git clone https://github.com/hms-dbmi/dropEst.git && \
  mkdir -p ~/R/x86_64-redhat-linux-gnu-library/3.4

RUN \
  echo '.libPaths(c("~/R/x86_64-redhat-linux-gnu-library/3.4", .libPaths()))' > .Rprofile && \
  R -e 'chooseCRANmirror(ind=52); install.packages("ks", dependencies=c("Depends", "Imports", "LinkingTo"))' && \
  R -e 'devtools::install_local("/software/dropEst/dropestr/", dependencies=T)' \
  R -e 'install.packages(c("optparse","rmarkdown"))'
RUN \
  mkdir -p dropEst/build && \
  cd dropEst/build && \
  cmake ../ && \
  make



ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/dropest/dropest_umi_correct.R /software/
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/dropest/dropseq_config.xml /software/
RUN chmod +x /software/*
RUN apt-get -qq -y remove git cmake build-essential \
&& apt-get -qq -y autoremove \
&& apt-get autoclean \
&& rm -rf /var/lib/apt/lists/*

ENV PATH=/software/dropEst/build:$PATH \
  LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH \
  R_PROFILE=~/.Rprofile
