#!/usr/bin/env bash

docker build -t demuxlet .
docker tag demuxlet "$DOCKER"/demuxlet:0.1-beta
docker push "$DOCKER"/demuxlet:0.1-beta
