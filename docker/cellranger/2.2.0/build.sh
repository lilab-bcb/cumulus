#!/usr/bin/env bash

docker build -t cellranger-2.2.0 .
docker tag cellranger-2.2.0 "$DOCKER"/cellranger:2.2.0
docker push "$DOCKER"/cellranger:2.2.0
