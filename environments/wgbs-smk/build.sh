#!/bin/bash

# # Test build
# docker build -t clarity001/wgbs-smk:test -f Dockerfile --progress=plain . 2>&1 | tee build.log

# Build and push to dockerhub
docker build -t clarity001/wgbs-smk:latest -f Dockerfile --progress=plain --push . 2>&1 | tee build.log
