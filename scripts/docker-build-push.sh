#!/bin/bash
# ==========================================
# Script de build et push des images Docker v3
# Avec les corrections Unicode et R
# ==========================================

set -e

REGISTRY="docker.io/serignebamba398"
VERSION="v3"
BACKEND_IMAGE="${REGISTRY}/annotation-nmr-backend:${VERSION}"
FRONTEND_IMAGE="${REGISTRY}/annotation-nmr-frontend:${VERSION}"

echo "=========================================="
echo "Building Docker images version ${VERSION}"
echo "=========================================="

# Build Backend
echo ""
echo "[1/4] Building backend image..."
docker build \
    -f docker/Dockerfile.backend \
    -t ${BACKEND_IMAGE} \
    --platform linux/amd64 \
    .

echo ""
echo "[2/4] Building frontend image..."
docker build \
    -f docker/Dockerfile.frontend \
    -t ${FRONTEND_IMAGE} \
    --platform linux/amd64 \
    .

echo ""
echo "[3/4] Pushing backend image..."
docker push ${BACKEND_IMAGE}

echo ""
echo "[4/4] Pushing frontend image..."
docker push ${FRONTEND_IMAGE}

echo ""
echo "=========================================="
echo "Done! Images pushed:"
echo "  - ${BACKEND_IMAGE}"
echo "  - ${FRONTEND_IMAGE}"
echo "=========================================="
echo ""
echo "Update docker-compose.shareable.yml with:"
echo "  backend:  image: ${BACKEND_IMAGE}"
echo "  frontend: image: ${FRONTEND_IMAGE}"
