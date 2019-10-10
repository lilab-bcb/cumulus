#!/usr/bin/env bash

docker build -t cellranger-2.2.0 .
docker tag cellranger-2.2.0 cumulusprod/cellranger:2.2.0
docker push cumulusprod/cellranger
