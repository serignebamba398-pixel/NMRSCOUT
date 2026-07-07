# ==========================================
# NMRScout — Installation rapide (Windows)
# ==========================================
# Usage: powershell -ExecutionPolicy Bypass -File install.ps1
# ==========================================

$ErrorActionPreference = "Stop"

$APP_NAME = "NMRScout"
$COMPOSE_FILE = "docker-compose.shareable.yml"
$ENV_FILE = ".env"
$ENV_EXAMPLE = ".env.example"

Write-Host "==========================================" -ForegroundColor Cyan
Write-Host "  $APP_NAME - Docker Installer (Windows)" -ForegroundColor Cyan
Write-Host "==========================================" -ForegroundColor Cyan

# ---- Check Docker ----
try {
    $null = docker --version 2>&1
} catch {
    Write-Host "[ERROR] Docker is not installed." -ForegroundColor Red
    Write-Host "Install Docker Desktop: https://docs.docker.com/desktop/install/windows-install/"
    exit 1
}

try {
    $null = docker info 2>&1
    if ($LASTEXITCODE -ne 0) { throw }
} catch {
    Write-Host "[ERROR] Docker daemon is not running. Start Docker Desktop and retry." -ForegroundColor Red
    exit 1
}

Write-Host "[OK] Docker is running" -ForegroundColor Green

# ---- Check Docker Compose ----
$COMPOSE_CMD = $null
try {
    $null = docker compose version 2>&1
    if ($LASTEXITCODE -eq 0) { $COMPOSE_CMD = "docker compose" }
} catch {}

if (-not $COMPOSE_CMD) {
    try {
        $null = docker-compose version 2>&1
        if ($LASTEXITCODE -eq 0) { $COMPOSE_CMD = "docker-compose" }
    } catch {}
}

if (-not $COMPOSE_CMD) {
    Write-Host "[ERROR] Docker Compose not found." -ForegroundColor Red
    Write-Host "Install Docker Desktop (includes Compose)."
    exit 1
}

Write-Host "[OK] Docker Compose available ($COMPOSE_CMD)" -ForegroundColor Green

# ---- Check required files ----
if (-not (Test-Path $COMPOSE_FILE)) {
    Write-Host "[ERROR] $COMPOSE_FILE not found." -ForegroundColor Red
    Write-Host "Make sure you run this script from the project root directory."
    exit 1
}

# ---- Create .env if missing ----
if (-not (Test-Path $ENV_FILE)) {
    if (Test-Path $ENV_EXAMPLE) {
        Copy-Item $ENV_EXAMPLE $ENV_FILE
        Write-Host "[INFO] Created .env from .env.example" -ForegroundColor Yellow
    } else {
        @"
PORT=4000
NODE_ENV=production
FRONTEND_URL=http://localhost:3000
CORS_ORIGINS=http://localhost:3000,http://127.0.0.1:3000
FRONTEND_PORT=3000
BACKEND_PORT=4000
R_PATH=/usr/bin/Rscript
PRED_PYTHON=/opt/conda/envs/nmr/bin/python
MODEL_1H_PATH=/app/backend/models/1H_sparsified_proposed_proposed_set2set_1.pt
MODEL_13C_PATH=/app/backend/models/13C_sparsified_proposed_proposed_set2set_1.pt
NMR_DATA_DIR=/app/backend/ml_dependencies/data
PYTHONIOENCODING=utf-8
PYTHONUTF8=1
PYTHONUNBUFFERED=1
DGLBACKEND=pytorch
BODY_LIMIT_MB=500
MAX_UPLOAD_MB=500
TASK_TTL_HOURS=48
PROCESS_TIMEOUT_HOURS=8
MEMORY_LIMIT_MB=6144
OMP_NUM_THREADS=4
MKL_NUM_THREADS=4
"@ | Out-File -FilePath $ENV_FILE -Encoding utf8
        Write-Host "[INFO] Created default .env file" -ForegroundColor Yellow
    }
} else {
    Write-Host "[OK] .env file exists" -ForegroundColor Green
}

# ---- Pull images ----
Write-Host ""
Write-Host "[1/2] Pulling Docker images (this may take a while on first run)..." -ForegroundColor Cyan
Invoke-Expression "$COMPOSE_CMD -f $COMPOSE_FILE pull"
if ($LASTEXITCODE -ne 0) { throw "Failed to pull images" }

# ---- Start services ----
Write-Host ""
Write-Host "[2/2] Starting services..." -ForegroundColor Cyan
Invoke-Expression "$COMPOSE_CMD -f $COMPOSE_FILE up -d"
if ($LASTEXITCODE -ne 0) { throw "Failed to start services" }

# ---- Wait for health ----
Write-Host ""
Write-Host "Waiting for backend to be healthy..." -ForegroundColor Yellow
$retries = 0
$maxRetries = 30
$healthy = $false

while ($retries -lt $maxRetries) {
    try {
        $response = Invoke-WebRequest -Uri "http://localhost:4000/health" -UseBasicParsing -TimeoutSec 3 -ErrorAction Stop
        if ($response.StatusCode -eq 200) {
            $healthy = $true
            break
        }
    } catch {}
    $retries++
    Start-Sleep -Seconds 3
}

if (-not $healthy) {
    Write-Host "[WARN] Backend not healthy yet. It may still be starting (ML models take time to load)." -ForegroundColor Yellow
    Write-Host "  Check logs: $COMPOSE_CMD -f $COMPOSE_FILE logs -f backend"
} else {
    Write-Host "[OK] Backend is healthy" -ForegroundColor Green
}

# ---- Done ----
Write-Host ""
Write-Host "==========================================" -ForegroundColor Green
Write-Host "  $APP_NAME is running!" -ForegroundColor Green
Write-Host "==========================================" -ForegroundColor Green
Write-Host ""
Write-Host "  Frontend : http://localhost:3000" -ForegroundColor White
Write-Host "  Backend  : http://localhost:4000" -ForegroundColor White
Write-Host ""
Write-Host "  Commands:" -ForegroundColor White
Write-Host "    Stop:    $COMPOSE_CMD -f $COMPOSE_FILE down"
Write-Host "    Logs:    $COMPOSE_CMD -f $COMPOSE_FILE logs -f"
Write-Host "    Restart: $COMPOSE_CMD -f $COMPOSE_FILE restart"
Write-Host "==========================================" -ForegroundColor Green
