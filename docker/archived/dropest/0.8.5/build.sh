#!/usr/bin/env bash

docker build -t dropest-0.8.5 .
docker tag dropest-0.8.5 cumulusprod/dropest:0.8.5
docker push cumulusprod/dropest:0.8.5

