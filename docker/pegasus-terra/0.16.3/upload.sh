#!/usr/bin/env bash

docker tag pegasus-terra-0.16.3 cumulusprod/pegasus-terra:0.16.3
docker login
docker push cumulusprod/pegasus-terra:0.16.3
docker tag pegasus-terra-0.16.3 quay.io/cumulus/pegasus-terra:0.16.3
docker login quay.io
docker push quay.io/cumulus/pegasus-terra:0.16.3