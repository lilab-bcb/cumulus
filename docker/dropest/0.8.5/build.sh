#!/usr/bin/env bash

docker build -t dropest-0.8.5 .
docker tag dropest-0.8.5 cumulus/dropest:0.8.5
docker push cumulus/dropest:0.8.5

