#!/usr/bin/env bash

docker build -t dropest-0.8.5 .
docker tag dropest-0.8.5 sccloud/dropest:0.8.5
docker push sccloud/dropest:0.8.5

