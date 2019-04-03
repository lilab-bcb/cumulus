#!/usr/bin/env bash

docker build -t dropest-0.8.5 .
docker tag dropest-0.8.5 regevlab/dropest-0.8.5
docker push regevlab/dropest-0.8.5

