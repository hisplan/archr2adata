#!/bin/bash -e

source config.sh

docker build -t archr2adata:${version} .
