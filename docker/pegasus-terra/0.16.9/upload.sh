#!/usr/bin/env bash

docker tag pegasus-terra-0.16.9 cumulusprod/pegasus-terra:0.16.9
docker login
docker push cumulusprod/pegasus-terra:0.16.9
docker tag pegasus-terra-0.16.9 quay.io/cumulus/pegasus-terra:0.16.9
docker login quay.io
docker push quay.io/cumulus/pegasus-terra:0.16.9