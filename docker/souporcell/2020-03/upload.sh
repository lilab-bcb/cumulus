#!/usr/bin/env bash

docker tag souporcell:2020.03 cumulusprod/souporcell:2020.03
docker login
docker push cumulusprod/souporcell:2020.03
docker tag souporcell:2020.03 quay.io/cumulus/souporcell:2020.03
docker login quay.io
docker push quay.io/cumulus/souporcell:2020.03