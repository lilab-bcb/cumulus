#!/usr/bin/env bash

# Download cellranger-3.0.2.tar.gz into this directory from cell ranger website before building
docker build -t cellranger-3.0.2 .
docker tag cellranger-3.0.2 regevlab/cellranger-3.0.2
docker push regevlab/cellranger-3.0.2
