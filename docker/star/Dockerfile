FROM debian:stretch-slim
ADD https://raw.githubusercontent.com/klarman-cell-observatory/KCO/master/docker/monitor_script.sh /software/
RUN apt-get -qq update && apt-get -qq -y install --no-install-recommends ca-certificates curl libgomp1 \
	&& curl -L -s -o "/2.5.4b.tar.gz" "https://github.com/alexdobin/STAR/archive/2.5.4b.tar.gz" \
	&& tar xf /2.5.4b.tar.gz \
	&& rm -rf /STAR-2.5.4b/bin/Linux_x86_64_static /STAR-2.5.4b/bin/MacOSX_x86_64 \
	&& mv /STAR-2.5.4b /software/ \
	&& rm -f /STAR-2.5.4b.tgz \
    && apt-get -qq -y remove curl \
    && apt-get -qq -y autoremove \
    && apt-get autoclean \
    && rm -rf /var/lib/apt/lists/* \
    && chmod a+rx /software/*

ENV PATH=/software:/software/STAR-2.5.4b/bin/Linux_x86_64/:${PATH}


