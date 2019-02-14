FROM openjdk:8-jre-slim
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software/
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl \
&& curl -L -s -o "/Drop-seq_tools-2.1.0.zip" "https://github.com/broadinstitute/Drop-seq/releases/download/v2.1.0/Drop-seq_tools-2.1.0.zip" \
&& unzip /Drop-seq_tools-2.1.0.zip \
&& rm /Drop-seq_tools-2.1.0.zip \
&& mv /Drop-seq_tools-2.1.0 /software/Drop-seq_tools-2.1.0 \
&& apt-get -qq -y remove curl \
&& apt-get -qq -y autoremove \
&& apt-get autoclean \
&& rm -rf /var/lib/apt/lists/* \
&& chmod a+rx /software/*

ENV PATH=/software/Drop-seq_tools-2.1.0/:${PATH}

