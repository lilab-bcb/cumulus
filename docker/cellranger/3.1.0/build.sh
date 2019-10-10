#!/usr/bin/env bash

# Download cellranger-3.1.0.tar.gz into this directory from cell ranger website before building
docker build -t cellranger-3.1.0 .
docker tag cellranger-3.1.0 cumulusprod/cellranger:3.1.0
docker push cumulusprod/cellranger
