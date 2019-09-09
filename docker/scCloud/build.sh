#!/usr/bin/env bash

docker build -t sccloud-0.9.0 .
docker tag sccloud-0.9.0 sccloud/sccloud:0.9.0
docker push sccloud/sccloud
