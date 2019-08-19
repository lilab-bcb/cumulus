#!/usr/bin/env bash

docker build -t dropest-0.8.6 .
docker tag dropest-0.8.6 sccloud/dropest:0.8.6
docker push sccloud/dropest:0.8.6

