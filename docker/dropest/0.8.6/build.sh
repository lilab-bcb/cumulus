#!/usr/bin/env bash

docker build -t dropest-0.8.6 .
docker tag dropest-0.8.6 cumulus/dropest:0.8.6
docker push cumulus/dropest:0.8.6

