#!/usr/bin/env bash

docker build -t star_2.7.0c .
docker tag star_2.7.0c regevlab/star_2.7.0c
docker push regevlab/star_2.7.0c
