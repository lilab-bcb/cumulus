#!/usr/bin/env bash

docker tag starsolo-2.7.3a cumulusprod/starsolo:2.7.3a
docker login
docker push cumulusprod/starsolo:2.7.3a
docker tag starsolo-2.7.3a quay.io/cumulus/starsolo:2.7.3a
docker login quay.io
docker push quay.io/cumulus/starsolo:2.7.3a