#!/usr/bin/env bash
# Build the sv-imputation Docker image and push to Docker Hub (ndatth/sv-imputation).
#
# The same image is also compatible with HPC environments via Singularity/Apptainer:
#   singularity pull sv-imputation.sif docker://ndatth/sv-imputation:latest
#   apptainer  pull sv-imputation.sif docker://ndatth/sv-imputation:latest
#
# Usage:
#   bash docker/build_and_upload.sh [TAG]
#
# TAG defaults to "latest".
# Example with a versioned tag:
#   bash docker/build_and_upload.sh 1.0.0
set -euo pipefail

DOCKERHUB_USER="ndatth"
IMAGE_NAME="sv-imputation"
TAG="${1:-latest}"
FULL_TAG="${DOCKERHUB_USER}/${IMAGE_NAME}:${TAG}"

# Resolve the directory containing this script so it works from any cwd
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"


echo "==> Building ${FULL_TAG} (linux/amd64) from ${SCRIPT_DIR}/Dockerfile"
docker build \
    --platform linux/amd64 \
    --file "${SCRIPT_DIR}/Dockerfile" \
    --tag "${FULL_TAG}" \
    "${SCRIPT_DIR}"

# Also tag as latest when a specific version is provided
if [[ "${TAG}" != "latest" ]]; then
    docker tag "${FULL_TAG}" "${DOCKERHUB_USER}/${IMAGE_NAME}:latest"
    echo "==> Also tagged as ${DOCKERHUB_USER}/${IMAGE_NAME}:latest"
fi

echo "==> Pushing ${FULL_TAG}"
docker push "${FULL_TAG}"

if [[ "${TAG}" != "latest" ]]; then
    echo "==> Pushing ${DOCKERHUB_USER}/${IMAGE_NAME}:latest"
    docker push "${DOCKERHUB_USER}/${IMAGE_NAME}:latest"
fi

echo "==> Done. Image available at: https://hub.docker.com/r/${DOCKERHUB_USER}/${IMAGE_NAME}"
