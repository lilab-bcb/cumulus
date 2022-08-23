#!/usr/bin/env bash

# Download cellranger-atac-1.2.0.tar.gz into this directory from cell ranger website before building

usage () {
    echo "Usage: build.sh [-b] [-p] docker_registry"
    echo "  -b: build docker but not push"
    echo "  -p: do not build docker, just push"
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
push=true

while getopts ":bp" opt; do
    case ${opt} in
        b ) push=false
            ;;
        p ) build=false
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
PLATFORM="linux/x86_64"

if [ "$build" = true ]
then
    echo "docker build -t ${DOCKER} --pull --platform ${PLATFORM} ."
    docker build -t "${DOCKER}" --pull --platform "${PLATFORM}" .
fi

if [ "$push" = true ]
then
    echo "docker tag $DOCKER $REGISTRY/$DOCKER"
    docker tag "$DOCKER" "$REGISTRY/$DOCKER"

    echo "docker push $REGISTRY/$DOCKER"
    docker push "$REGISTRY/$DOCKER"
fi
