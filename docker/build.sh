#!/usr/bin/env bash

# Download cellranger-atac-1.2.0.tar.gz into this directory from cell ranger website before building

usage () {
    echo "Usage: build.sh [-n] docker_registry"
    echo "  -n: do not build docker, just push"
    echo "  docker_registry: cumulusprod, quay.io/cumulus or gcr.io/broad-cumulus"
    echo "Note: run this script under where the dockerfile locate. We assume the Dockerfile locates under SOFTWARE/VERSION. For example, cellranger/3.1.0/Dockerfile. In this case, we can automatically determine the DOCKER name and tag VERSION based on the working directory."
    exit 1
}

IFS='/'
read -ra TOKENS <<< "$PWD"
TLEN=${#TOKENS[@]}
SOFTWARE="${TOKENS[$TLEN-2]}"
VERSION="${TOKENS[$TLEN-1]}"
DOCKER="$SOFTWARE:$VERSION"


build=true

while getopts ":n" opt; do
    case ${opt} in
        n ) build=false
            ;;
        \? ) usage 
            ;;
    esac
    shift $((OPTIND -1))
done

if [ $# -ne 1 ]
then
    usage
fi

REGISTRY=$1

if [ "$build" = true ]
then
    echo "docker build -t ${DOCKER} ."
    docker build -t "${DOCKER}" .
fi

echo "docker tag $DOCKER $REGISTRY/$DOCKER"
docker tag "$DOCKER" "$REGISTRY/$DOCKER"

echo "docker push $REGISTRY/$DOCKER"
docker push "$REGISTRY/$DOCKER"
