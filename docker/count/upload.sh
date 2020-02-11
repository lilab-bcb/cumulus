#!/usr/bin/env bash

docker tag count cumulusprod/count
docker login
docker push cumulusprod/count
docker tag count quay.io/cumulus/count
docker login quay.io
docker push quay.io/cumulus/count