#!/usr/bin/env bash

# Download cellranger-2.0.2.tar into this directory from cell ranger website before building
docker build -t cellranger-2.0.2 .
docker tag cellranger2.0.2 regevlab/cellranger:2.0.2
docker push regevlab/cellranger:2.0.2
