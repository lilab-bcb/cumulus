FROM java:openjdk-8-jre
MAINTAINER <kco-help@broadinstitute.org>
ENV TERM=xterm-256color
ARG SAMTOOLS_VER=1.3.1
ARG R_VER=3.4.3
ARG PYTHON_VER=2.7
ARG PICARD_VER=1.141
WORKDIR /home
# Install python
RUN apt-get update --fix-missing && apt-get upgrade -y && \
  apt-get install -y python2.7 && \
  apt-get install -y r-base && \
  apt-get install -y zlib1g-dev && \
  apt-get install -y git-all
  
RUN ln -s /usr/bin/python2.7 /usr/bin/python

RUN wget --quiet -O Drop-seq_tools-1.13-3.zip http://mccarrolllab.com/download/1276/ && \
  unzip Drop-seq_tools-1.13-3.zip && \
  rm Drop-seq_tools-1.13-3.zip
  
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('ineq')"
RUN Rscript -e "install.packages('optparse')"
RUN Rscript -e "install.packages('MASS')"
RUN Rscript -e "install.packages('data.table')"
RUN Rscript -e "install.packages('RColorBrewer')"
RUN Rscript -e "install.packages('caroline')"
RUN Rscript -e "install.packages('gplots')"
RUN Rscript -e "install.packages('stringdist')"
RUN Rscript -e "install.packages('IRanges')"
RUN Rscript -e "install.packages('gridExtra')"
RUN Rscript -e "install.packages('gtable')"
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("Biostrings")'
RUN Rscript -e 'source("http://bioconductor.org/biocLite.R")' -e 'biocLite("IRanges")'
RUN R CMD INSTALL Drop-seq_tools-1.13/DropSeq.barnyard_1.13.tar.gz
ENV R_LIBS="/home/Drop-seq_tools-1.13/DropSeq.barnyard_1.13"
# Install STAR
RUN git clone https://github.com/alexdobin/STAR.git && \
    cd STAR && \
    git checkout 2.5.1b && \
    cd source && \
    make STAR && \
    cd ../..
# Install Utility for dropseq
RUN git clone https://github.com/broadinstitute/DropSeq_utility.git
# Install Picard    
RUN wget https://github.com/broadinstitute/picard/releases/download/1.141/picard-tools-1.141.zip -O picard-tools-1.141.zip && \
  unzip -a picard-tools-1.141.zip && \
  rm picard-tools-1.141.zip && \
  ln picard-tools-1.141/picard.jar picard.jar
# Install Samtools
RUN wget https://github.com/samtools/samtools/releases/download/1.3.1/samtools-${SAMTOOLS_VER}.tar.bz2 -O samtools-1.3.1.tar.bz2 && \
  tar xf samtools-1.3.1.tar.bz2 && \
  rm samtools-1.3.1.tar.bz2 && \
  cd samtools-1.3.1 && \
  make && \
  make install && \
  cd ../ && \
  rm -r samtools-1.3.1
# Install Star
ENV PATH "/home/Drop-seq_tools-1.13:/home/STAR/source:/home/DropSeq_utility:/usr/bin/python:${PATH}"
ENTRYPOINT []
