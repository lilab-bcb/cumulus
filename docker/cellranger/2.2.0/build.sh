#!/usr/bin/env bash

docker build -t cellranger-2.2.0 .
docker tag cellranger-2.2.0 regevlab/cellranger-2.2.0
docker push regevlab/cellranger-2.2.0
