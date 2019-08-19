#!/usr/bin/env bash

docker build -t dropseq-2.3.0 .
docker tag dropseq-2.3.0 sccloud/dropseq:2.3.0
docker push sccloud/dropseq
