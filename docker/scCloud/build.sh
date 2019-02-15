#!/usr/bin/env bash

docker build -t sccloud-0.7.0 .
docker tag sccloud-0.7.0 regevlab/sccloud-0.7.0
docker push regevlab/sccloud-0.7.0
