FROM debian:stretch-slim
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software/
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl libgomp1 \
	&& curl -L -s -o "/2.7.0c.tar.gz" "https://github.com/alexdobin/STAR/archive/2.7.0c.tar.gz" \
	&& tar xf /2.7.0c.tar.gz \
	&& rm -rf /STAR-2.7.0c/bin/Linux_x86_64_static /STAR-2.7.0c/bin/MacOSX_x86_64 \
	&& mv /STAR-2.7.0c /software/ \
	&& rm -f /STAR-2.7.0c.tgz \
    && apt-get -qq -y remove curl \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* \
    && chmod a+rx /software/*

ENV PATH=/software:/software/STAR-2.7.0c/bin/Linux_x86_64/:${PATH}


