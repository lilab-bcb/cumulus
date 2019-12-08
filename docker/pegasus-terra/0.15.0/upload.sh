#!/usr/bin/env bash

docker tag pegasus-terra-0.15.0 cumulusprod/pegasus-terra:0.15.0
docker login
docker push cumulusprod/pegasus-terra:0.15.0
docker tag pegasus-terra-0.15.0 quay.io/cumulus/pegasus-terra:0.15.0
docker login quay.io
docker push quay.io/cumulus/pegasus-terra:0.15.0