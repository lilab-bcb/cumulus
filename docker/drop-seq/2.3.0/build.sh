#!/usr/bin/env bash

docker build -t dropseq-2.3.0 .
docker tag dropseq-2.3.0 cumulusprod/dropseq:2.3.0
docker push cumulusprod/dropseq
