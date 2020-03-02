#!/usr/bin/env bash

docker tag pegasus-terra-0.16.11 cumulusprod/pegasus-terra:0.16.11
docker login
docker push cumulusprod/pegasus-terra:0.16.11
docker tag pegasus-terra-0.16.11 quay.io/cumulus/pegasus-terra:0.16.11
docker login quay.io
docker push quay.io/cumulus/pegasus-terra:0.16.11