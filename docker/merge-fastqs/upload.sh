#!/usr/bin/env bash

docker tag merge-fastqs cumulusprod/merge-fastqs
docker login
docker push cumulusprod/merge-fastqs
docker tag merge-fastqs quay.io/cumulus/merge-fastqs
docker login quay.io
docker push quay.io/cumulus/merge-fastqs