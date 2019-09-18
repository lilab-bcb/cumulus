#!/usr/bin/env bash

docker build -t sccloud-0.9.1 .
docker tag sccloud-0.9.1 sccloud/sccloud:0.9.1
docker push sccloud/sccloud
