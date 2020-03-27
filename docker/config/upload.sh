#!/usr/bin/env bash

docker tag config cumulusprod/config
docker login
docker push cumulusprod/config
docker tag config quay.io/cumulus/config
docker login quay.io
docker push quay.io/cumulus/config