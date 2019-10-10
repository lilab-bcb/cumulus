#!/usr/bin/env bash

docker build -t smartseq2-1.0.0 .
docker tag smartseq2-1.0.0 cumulusprod/smartseq2:1.0.0
docker push cumulusprod/smartseq2
