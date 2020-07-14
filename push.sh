#!/bin/bash -e

hub="hisplan"
version="0.0.1"

echo "Packaging ${hub}/archr2adata:${version}..."

#
# tag it and push it to docker hub
#

docker login
docker tag archr2adata:${version} ${hub}/archr2adata:${version}
docker push ${hub}/archr2adata:${version}