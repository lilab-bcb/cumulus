#!/usr/bin/env bash

docker build --no-cache -t cumulusprod-0.10.0 .
docker tag cumulusprod-0.10.0  "$DOCKER"/cumulusprod:0.10.0
docker push "$DOCKER"/cumulusprod:0.10.0
