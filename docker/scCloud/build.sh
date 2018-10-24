#!/usr/bin/env bash

docker build -t sccloud .
docker tag sccloud regevlab/sccloud
docker push regevlab/sccloud
