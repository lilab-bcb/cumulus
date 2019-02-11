#!/usr/bin/env bash

docker build -t star .
docker tag star regevlab/star
docker push regevlab/star
