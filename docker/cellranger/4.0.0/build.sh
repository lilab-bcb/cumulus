#!/usr/bin/env bash

# Download cellranger-4.0.0.tar.gz into this directory from cell ranger website before building
docker build -t cellranger-4.0.0 .
docker tag cellranger-4.0.0 "$DOCKER"/cellranger:4.0.0
docker push "$DOCKER"/cellranger:4.0.0
