#!/bin/bash
# ==========================================
# Script de build et push des images Docker
# Multi-platform: linux/amd64 + linux/arm64
# Compatible Windows (x64), Mac Intel, Mac ARM, Linux
# ==========================================

set -e

REGISTRY="docker.io/serignebamba398"
VERSION="${1:-v1}"
PLATFORMS="linux/amd64,linux/arm64"
BACKEND_IMAGE="${REGISTRY}/annotation-nmr-backend:${VERSION}"
FRONTEND_IMAGE="${REGISTRY}/annotation-nmr-frontend:${VERSION}"
BUILDER_NAME="nmr-multiarch"

echo "=========================================="
echo "Multi-platform Docker build (${PLATFORMS})"
echo "Version: ${VERSION}"
echo "=========================================="

# ---- Ensure buildx builder exists ----
if ! docker buildx inspect "${BUILDER_NAME}" > /dev/null 2>&1; then
    echo "[setup] Creating buildx builder '${BUILDER_NAME}'..."
    docker buildx create --name "${BUILDER_NAME}" --driver docker-container --use
else
    docker buildx use "${BUILDER_NAME}"
fi

# ---- Build & push backend (multi-arch) ----
echo ""
echo "[1/2] Building & pushing backend (multi-arch)..."
docker buildx build \
    --platform "${PLATFORMS}" \
    -f docker/Dockerfile.backend \
    -t "${BACKEND_IMAGE}" \
    --push \
    .

# ---- Build & push frontend (multi-arch) ----
echo ""
echo "[2/2] Building & pushing frontend (multi-arch)..."
docker buildx build \
    --platform "${PLATFORMS}" \
    -f docker/Dockerfile.frontend \
    -t "${FRONTEND_IMAGE}" \
    --push \
    .

echo ""
echo "=========================================="
echo "Done! Multi-arch images pushed:"
echo "  - ${BACKEND_IMAGE}"
echo "  - ${FRONTEND_IMAGE}"
echo ""
echo "Supported platforms: ${PLATFORMS}"
echo "=========================================="
echo ""
echo "Users can now run on any OS with:"
echo "  docker compose -f docker-compose.shareable.yml up -d"
