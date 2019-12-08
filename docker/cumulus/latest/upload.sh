#!/usr/bin/env bash

docker tag cumulus-latest cumulusprod/cumulus:latest
docker login
docker push cumulusprod/cumulus:latest
docker tag cumulus-latest quay.io/cumulus/cumulus:latest
docker login quay.io
docker push quay.io/cumulus/cumulus:latest