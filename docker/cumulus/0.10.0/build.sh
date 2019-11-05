#!/usr/bin/env bash

docker build -t cumulus-0.10.0 .
docker tag cumulus-test cumulusprod/cumulus:0.10.0
docker push cumulusprod/cumulus:0.10.0
#docker tag cumulus-test quay.io/cumulus:0.10.0
#docker push quay.io/cumulus:0.10.0