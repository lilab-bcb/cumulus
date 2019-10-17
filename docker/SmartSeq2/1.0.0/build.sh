#!/usr/bin/env bash

docker build -t smartseq2-1.0.0 .
docker tag smartseq2-1.0.0 "$DOCKER"/smartseq2:1.0.0
docker push "$DOCKER"/smartseq2:1.0.0
