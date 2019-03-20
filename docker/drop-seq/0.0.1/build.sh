#!/usr/bin/env bash

docker build -t dropseq-0.0.1 .
docker tag dropseq-0.0.1 regevlab/dropseq-0.0.1
docker push regevlab/dropseq-0.0.1

