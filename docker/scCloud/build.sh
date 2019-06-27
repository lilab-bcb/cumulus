#!/usr/bin/env bash

docker build -t sccloud-0.8.0 .
docker tag sccloud-0.8.0 regevlab/sccloud-0.8.0:latest
docker push regevlab/sccloud-0.8.0:latest
