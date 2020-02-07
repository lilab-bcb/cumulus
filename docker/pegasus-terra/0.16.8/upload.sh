#!/usr/bin/env bash

docker tag pegasus-terra-0.16.8 cumulusprod/pegasus-terra:0.16.8
docker login
docker push cumulusprod/pegasus-terra:0.16.8
docker tag pegasus-terra-0.16.8 quay.io/cumulus/pegasus-terra:0.16.8
docker login quay.io
docker push quay.io/cumulus/pegasus-terra:0.16.8