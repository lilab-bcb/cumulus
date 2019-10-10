#!/usr/bin/env bash

docker build --no-cache -t cumulusprod-0.10.0 .
docker tag cumulusprod-0.10.0 cumulusprod/cumulusprod:0.10.0
docker push cumulusprod/cumulusprod
