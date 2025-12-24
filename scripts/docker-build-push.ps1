# ==========================================
# Script de build et push des images Docker v3
# Pour Windows PowerShell
# ==========================================

$ErrorActionPreference = "Stop"

$REGISTRY = "docker.io/serignebamba398"
$VERSION = "v3"
$BACKEND_IMAGE = "${REGISTRY}/annotation-nmr-backend:${VERSION}"
$FRONTEND_IMAGE = "${REGISTRY}/annotation-nmr-frontend:${VERSION}"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "Building Docker images version $VERSION" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan

# Build Backend
Write-Host "`n[1/4] Building backend image..." -ForegroundColor Yellow
docker build -f docker/Dockerfile.backend -t $BACKEND_IMAGE .
if ($LASTEXITCODE -ne 0) { throw "Backend build failed" }

Write-Host "`n[2/4] Building frontend image..." -ForegroundColor Yellow
docker build -f docker/Dockerfile.frontend -t $FRONTEND_IMAGE .
if ($LASTEXITCODE -ne 0) { throw "Frontend build failed" }

Write-Host "`n[3/4] Pushing backend image..." -ForegroundColor Yellow
docker push $BACKEND_IMAGE
if ($LASTEXITCODE -ne 0) { throw "Backend push failed" }

Write-Host "`n[4/4] Pushing frontend image..." -ForegroundColor Yellow
docker push $FRONTEND_IMAGE
if ($LASTEXITCODE -ne 0) { throw "Frontend push failed" }

Write-Host "`n==========================================" -ForegroundColor Green
Write-Host "Done! Images pushed:" -ForegroundColor Green
Write-Host "  - $BACKEND_IMAGE" -ForegroundColor Green
Write-Host "  - $FRONTEND_IMAGE" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host "`nUpdate docker-compose.shareable.yml with:"
Write-Host "  backend:  image: $BACKEND_IMAGE"
Write-Host "  frontend: image: $FRONTEND_IMAGE"
