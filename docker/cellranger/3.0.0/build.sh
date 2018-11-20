#!/usr/bin/env bash

# Download cellranger-3.0.0.tar into this directory from cell ranger website before building
docker build -t cellranger-3.0.0 .
docker tag cellranger-3.0.0 regevlab/cellranger-3.0.0
docker push regevlab/cellranger-3.0.0
