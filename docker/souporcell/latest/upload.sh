#!/usr/bin/env bash

docker tag souporcell:latest cumulusprod/souporcell:latest
docker login
docker push cumulusprod/souporcell:latest
docker tag souporcell:latest quay.io/cumulus/souporcell:latest
docker login quay.io
docker push quay.io/cumulus/souporcell:latest