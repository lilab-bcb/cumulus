#!/usr/bin/env bash

docker tag alevin-1.1 cumulusprod/alevin:1.1
docker login
docker push cumulusprod/alevin:1.1
docker tag alevin-1.1 quay.io/cumulus/alevin:1.1
docker login quay.io
docker push quay.io/cumulus/alevin:1.1