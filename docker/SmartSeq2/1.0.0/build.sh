#!/usr/bin/env bash

docker build -t smartseq2-1.0.0 .
docker tag smartseq2-1.0.0 cumulus/smartseq2:1.0.0
docker push cumulus/smartseq2
