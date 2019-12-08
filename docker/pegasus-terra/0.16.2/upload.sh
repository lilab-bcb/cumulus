#!/usr/bin/env bash

docker tag pegasus-terra-0.16.2 cumulusprod/pegasus-terra:0.16.2
docker login
docker push cumulusprod/pegasus-terra:0.16.2
docker tag pegasus-terra-0.16.2 quay.io/cumulus/pegasus-terra:0.16.2
docker login quay.io
docker push quay.io/cumulus/pegasus-terra:0.16.2