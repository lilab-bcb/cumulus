#!/usr/bin/env bash

docker build -t smartseq2 .
docker tag smartseq2 regevlab/smartseq2
docker push regevlab/smartseq2
