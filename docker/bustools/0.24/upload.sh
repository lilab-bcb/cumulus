#!/usr/bin/env bash

docker tag bustools-0.24 cumulusprod/bustools:0.24
docker login
docker push cumulusprod/bustools:0.24
docker tag bustools-0.24 quay.io/cumulus/bustools:0.24
docker login quay.io
docker push quay.io/cumulus/bustools:0.24