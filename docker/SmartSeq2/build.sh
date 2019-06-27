#!/usr/bin/env bash

docker build -t smartseq2-0.2.0 .
docker tag smartseq2-0.2.0 regevlab/smartseq2:0.2.0
docker push regevlab/smartseq2:0.2.0
