#!/usr/bin/env bash

docker build --no-cache -t cumulus-0.10.0 .
docker tag cumulus-0.10.0 "$DOCKER"/cumulus:0.10.0
docker push "$DOCKER"/cumulus:0.10.0
