#!/usr/bin/env bash

docker build -t sccloud-0.14.0 .
docker tag sccloud-0.14.0 sccloud/sccloud:0.14.0
docker push sccloud/sccloud
