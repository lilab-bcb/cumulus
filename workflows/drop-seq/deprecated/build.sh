#!/usr/bin/env bash

docker build -t dropseq_v5 .
docker tag dropseq_v5 regevlab/dropseq_v5
docker push regevlab/dropseq_v5
