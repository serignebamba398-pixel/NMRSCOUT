# ==========================================
# Script de build et push des images Docker
# Multi-platform: linux/amd64 + linux/arm64
# Compatible Windows (x64), Mac Intel, Mac ARM, Linux
# ==========================================

$ErrorActionPreference = "Stop"

$REGISTRY = "docker.io/serignebamba398"
$VERSION = if ($args[0]) { $args[0] } else { "v1" }
$PLATFORMS = "linux/amd64,linux/arm64"
$BACKEND_IMAGE = "${REGISTRY}/annotation-nmr-backend:${VERSION}"
$FRONTEND_IMAGE = "${REGISTRY}/annotation-nmr-frontend:${VERSION}"
$BUILDER_NAME = "nmr-multiarch"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Multi-platform Docker build ($PLATFORMS)" -ForegroundColor Cyan
Write-Host "Version: $VERSION" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan

# ---- Ensure buildx builder exists ----
$ErrorActionPreference = "Continue"
docker buildx inspect $BUILDER_NAME 2>$null | Out-Null
$builderFound = $LASTEXITCODE -eq 0
$ErrorActionPreference = "Stop"

if (-not $builderFound) {
    Write-Host "[setup] Creating buildx builder '$BUILDER_NAME'..." -ForegroundColor Yellow
    docker buildx create --name $BUILDER_NAME --driver docker-container --use
    if ($LASTEXITCODE -ne 0) { throw "Failed to create buildx builder" }
} else {
    docker buildx use $BUILDER_NAME
}

# ---- Build & push backend (multi-arch) ----
Write-Host "`n[1/2] Building & pushing backend (multi-arch)..." -ForegroundColor Yellow
docker buildx build `
    --platform $PLATFORMS `
    -f docker/Dockerfile.backend `
    -t $BACKEND_IMAGE `
    --push `
    .
if ($LASTEXITCODE -ne 0) { throw "Backend build failed" }

# ---- Build & push frontend (multi-arch) ----
Write-Host "`n[2/2] Building & pushing frontend (multi-arch)..." -ForegroundColor Yellow
docker buildx build `
    --platform $PLATFORMS `
    -f docker/Dockerfile.frontend `
    -t $FRONTEND_IMAGE `
    --push `
    .
if ($LASTEXITCODE -ne 0) { throw "Frontend build failed" }

Write-Host "`n==========================================" -ForegroundColor Green
Write-Host "Done! Multi-arch images pushed:" -ForegroundColor Green
Write-Host "  - $BACKEND_IMAGE" -ForegroundColor Green
Write-Host "  - $FRONTEND_IMAGE" -ForegroundColor Green
Write-Host ""
Write-Host "Supported platforms: $PLATFORMS" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host "`nUsers can now run on any OS with:"
Write-Host "  docker compose -f docker-compose.shareable.yml up -d"
