#!/bin/bash -e

source config.sh

hub="hisplan"

echo "Packaging ${hub}/archr2adata:${version}..."

#
# tag it and push it to docker hub
#

docker login
docker tag archr2adata:${version} ${hub}/archr2adata:${version}
docker push ${hub}/archr2adata:${version}
