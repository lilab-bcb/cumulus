#!/usr/bin/env bash

docker build -t dropest-0.8.6 .
docker tag dropest-0.8.6 "$DOCKER"/dropest:0.8.6
docker push "$DOCKER"/dropest:0.8.6
