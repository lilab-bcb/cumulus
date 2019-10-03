#!/usr/bin/env bash

docker build -t cumulus-0.9.1 .
docker tag cumulus-0.9.1 cumulus/cumulus:0.9.1
docker push cumulus/cumulus
