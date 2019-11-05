#!/usr/bin/env bash

docker build -t cumulus-test .
docker tag cumulus-test cumulusprod/cumulus:test
docker login
docker push cumulusprod/cumulus:test
docker tag cumulus-test quay.io/cumulus/cumulus:test
docker login quay.io
docker push quay.io/cumulus/cumulus:test