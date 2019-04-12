FROM rocker/r-ver:3.5

RUN apt-get update && \
apt-get install --no-install-recommends -y ca-certificates curl dpkg-dev lsb-release procps gnupg2 && \
export CLOUD_SDK_REPO="cloud-sdk-$(lsb_release -c -s)" && \
echo "deb http://packages.cloud.google.com/apt $CLOUD_SDK_REPO main" | tee -a /etc/apt/sources.list.d/google-cloud-sdk.list && \
curl https://packages.cloud.google.com/apt/doc/apt-key.gpg | apt-key add - && \
apt-get update && apt-get install -y google-cloud-sdk


RUN apt-get update --yes && apt-get install --no-install-recommends --yes \
  build-essential \
  cmake \
  git \
  libbamtools-dev \
  libboost-dev \
  libboost-iostreams-dev \
  libboost-log-dev \
  libboost-system-dev \
  libboost-test-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  libz-dev \
  libatlas3-base \
  pandoc

RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile && \
  R -e 'install.packages(c("devtools", "Rcpp","RcppEigen", "RInside", "Matrix", "optparse", "rmarkdown", "withr"), quiet=TRUE)' && \
  R -e 'install.packages("ks", dependencies=c("Depends", "Imports", "LinkingTo"), quiet=TRUE)'

RUN mkdir -p /software
WORKDIR "/software"

RUN \
  git clone https://github.com/hms-dbmi/dropEst.git

RUN \
  R -e 'devtools::install_local("/software/dropEst/dropestr/", dependencies=T)'

RUN \
  mkdir -p dropEst/build && \
  cd dropEst/build && \
  cmake ../ && \
  make


ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/dropest/dropest_umi_correct.R /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/dropest/create_dropest_config.py /software
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software

RUN rm ~/.Rprofile
RUN chmod +x /software/*
RUN apt-get -qq -y remove curl git cmake build-essential \
&& apt-get -qq -y autoremove \
&& apt-get autoclean \
&& rm -rf /var/lib/apt/lists/*

ENV PATH=/software:/software/dropEst/build:$PATH \
  LD_LIBRARY_PATH=/usr/local/lib/R/lib/:$LD_LIBRARY_PATH \
  R_PROFILE=~/.Rprofile
