#!/usr/bin/env bash

docker tag demuxlet cumulusprod/demuxlet:0.1b
docker login
docker push cumulusprod/demuxlet:0.1b
docker tag demuxlet quay.io/cumulus/demuxlet:0.1b
docker login quay.io
docker push quay.io/cumulus/demuxlet:0.1b