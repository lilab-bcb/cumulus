#!/usr/bin/env bash

docker build -t cumulus-test .
docker tag cumulus-test cumulusprod/cumulus:test
docker push cumulusprod/cumulus:test
#docker tag cumulus-test quay.io/cumulus:test
#docker push quay.io/cumulus:test