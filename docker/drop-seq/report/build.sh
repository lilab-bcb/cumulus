#!/usr/bin/env bash

docker build -t dropseq_report-1.0.1 .
docker tag dropseq_report-1.0.1 "$DOCKER"/dropseq_report:1.0.1
docker push "$DOCKER"/dropseq_report:1.0.1
