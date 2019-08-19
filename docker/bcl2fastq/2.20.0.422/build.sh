#!/usr/bin/env bash

docker build -t bcl2fastq-2.20.0.422 .
docker tag bcl2fastq-2.20.0.422 sccloud/bcl2fastq:2.20.0.422
docker push sccloud/bcl2fastq:2.20.0.422

