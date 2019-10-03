#!/usr/bin/env bash

docker build -t cellranger-2.2.0 .
docker tag cellranger-2.2.0 cumulus/cellranger:2.2.0
docker push cumulus/cellranger
