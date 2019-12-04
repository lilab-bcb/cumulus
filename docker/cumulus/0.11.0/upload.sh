#!/usr/bin/env bash

docker tag cumulus-0.11.0 cumulusprod/cumulus:0.11.0
docker login
docker push cumulusprod/cumulus:0.11.0
docker tag cumulus-0.11.0 quay.io/cumulus/cumulus:0.11.0
docker login quay.io
docker push quay.io/cumulus/cumulus:0.11.0