#!/usr/bin/env bash

docker build -t dropseq-2.3.0 .
docker tag dropseq-2.3.0 cumulus/dropseq:2.3.0
docker push cumulus/dropseq
