#!/usr/bin/env bash

docker build -t dropseq_report-1.0.0 .
docker tag dropseq_report-1.0.0 "$DOCKER"/dropseq_report:1.0.0
docker push "$DOCKER"/dropseq_report:1.0.0
