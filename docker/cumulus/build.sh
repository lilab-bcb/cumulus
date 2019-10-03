#!/usr/bin/env bash

docker build --no-cache -t cumulus-0.10.0 .
docker tag cumulus-0.10.0 cumulusprod/cumulus:0.10.0
docker push cumulusprod/cumulus
