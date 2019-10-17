#!/usr/bin/env bash

docker build -t dropseq-2.3.0 .
docker tag dropseq-2.3.0 "$DOCKER"/dropseq:2.3.0
docker push "$DOCKER"/dropseq:2.3.0
