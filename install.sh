#!/bin/bash
# ==========================================
# NMRScout — Installation rapide (macOS / Linux)
# ==========================================
# Usage: curl -fsSL <url>/install.sh | bash
#   ou : bash install.sh
# ==========================================

set -e

APP_NAME="NMRScout"
COMPOSE_FILE="docker-compose.shareable.yml"
VOLUMES_FILE="docker-compose.volumes.yml"
ENV_FILE=".env"
ENV_EXAMPLE=".env.example"
REPO_URL="https://github.com/serignebamba398/annotation-nmr"

RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
NC='\033[0m'

echo -e "${CYAN}=========================================="
echo "  ${APP_NAME} — Docker Installer"
echo -e "==========================================${NC}"

# ---- Check Docker ----
if ! command -v docker &> /dev/null; then
    echo -e "${RED}[ERROR] Docker is not installed.${NC}"
    echo "Install Docker Desktop:"
    echo "  macOS:  https://docs.docker.com/desktop/install/mac-install/"
    echo "  Linux:  https://docs.docker.com/engine/install/"
    exit 1
fi

if ! docker info &> /dev/null; then
    echo -e "${RED}[ERROR] Docker daemon is not running. Start Docker Desktop and retry.${NC}"
    exit 1
fi

echo -e "${GREEN}[OK] Docker is running${NC}"

# ---- Check Docker Compose ----
if docker compose version &> /dev/null; then
    COMPOSE_CMD="docker compose"
elif command -v docker-compose &> /dev/null; then
    COMPOSE_CMD="docker-compose"
else
    echo -e "${RED}[ERROR] Docker Compose not found.${NC}"
    echo "Install Docker Desktop (includes Compose) or install docker-compose separately."
    exit 1
fi

echo -e "${GREEN}[OK] Docker Compose available ($COMPOSE_CMD)${NC}"

# ---- Check required files ----
if [ ! -f "$COMPOSE_FILE" ]; then
    echo -e "${RED}[ERROR] $COMPOSE_FILE not found.${NC}"
    echo "Make sure you run this script from the project root directory."
    echo "  git clone $REPO_URL && cd annotation-nmr"
    exit 1
fi

# ---- Create .env if missing ----
if [ ! -f "$ENV_FILE" ]; then
    if [ -f "$ENV_EXAMPLE" ]; then
        cp "$ENV_EXAMPLE" "$ENV_FILE"
        echo -e "${YELLOW}[INFO] Created .env from .env.example${NC}"
    else
        cat > "$ENV_FILE" << 'ENVEOF'
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
ENVEOF
        echo -e "${YELLOW}[INFO] Created default .env file${NC}"
    fi
else
    echo -e "${GREEN}[OK] .env file exists${NC}"
fi

# ---- Generate docker-compose.volumes.yml (detect mount points) ----
echo ""
echo -e "${YELLOW}[INFO] Detecting mount points for spectrum access...${NC}"

VOLUMES_LINES=""

if [[ "$(uname)" == "Darwin" ]]; then
    # macOS: mount /Users and any /Volumes/* (external drives, USB)
    VOLUMES_LINES="      - /Users:/host/Users:ro"
    echo -e "  /Users -> /host/Users"

    if [ -d "/Volumes" ]; then
        for vol in /Volumes/*/; do
            vol="${vol%/}"  # remove trailing slash
            volname="$(basename "$vol")"
            # Skip "Macintosh HD" (already covered by /Users)
            if [ "$volname" != "Macintosh HD" ] && [ -d "$vol" ]; then
                VOLUMES_LINES="${VOLUMES_LINES}
      - ${vol}:/host/Volumes/${volname}:ro"
                echo -e "  ${vol} -> /host/Volumes/${volname}"
            fi
        done
    fi
else
    # Linux: mount /home and any /media, /mnt, /run/media (external drives, USB)
    VOLUMES_LINES="      - /home:/host/home:ro"
    echo -e "  /home -> /host/home"

    # /media/<user>/* (Ubuntu/Debian auto-mount)
    if [ -d "/media" ]; then
        for mnt in /media/*/*/; do
            mnt="${mnt%/}"
            if [ -d "$mnt" ]; then
                VOLUMES_LINES="${VOLUMES_LINES}
      - ${mnt}:/host${mnt}:ro"
                echo -e "  ${mnt} -> /host${mnt}"
            fi
        done
    fi

    # /mnt/* (manual mounts)
    if [ -d "/mnt" ]; then
        for mnt in /mnt/*/; do
            mnt="${mnt%/}"
            if [ -d "$mnt" ]; then
                VOLUMES_LINES="${VOLUMES_LINES}
      - ${mnt}:/host${mnt}:ro"
                echo -e "  ${mnt} -> /host${mnt}"
            fi
        done
    fi

    # /run/media/<user>/* (Fedora/Arch auto-mount)
    if [ -d "/run/media" ]; then
        for mnt in /run/media/*/*/; do
            mnt="${mnt%/}"
            if [ -d "$mnt" ]; then
                VOLUMES_LINES="${VOLUMES_LINES}
      - ${mnt}:/host${mnt}:ro"
                echo -e "  ${mnt} -> /host${mnt}"
            fi
        done
    fi
fi

cat > "$VOLUMES_FILE" << VYEOF
# Auto-generated by install.sh - DO NOT EDIT
# Mounts host directories so Docker can read Bruker spectra
services:
  backend:
    volumes:
${VOLUMES_LINES}
VYEOF

echo -e "${GREEN}[OK] Generated $VOLUMES_FILE${NC}"

# ---- Pull images ----
echo ""
echo -e "${CYAN}[1/2] Pulling Docker images (this may take a while on first run)...${NC}"
$COMPOSE_CMD -f "$COMPOSE_FILE" -f "$VOLUMES_FILE" pull

# ---- Start services ----
echo ""
echo -e "${CYAN}[2/2] Starting services...${NC}"
$COMPOSE_CMD -f "$COMPOSE_FILE" -f "$VOLUMES_FILE" up -d

# ---- Wait for health ----
echo ""
echo -e "${YELLOW}Waiting for backend to be healthy...${NC}"
RETRIES=0
MAX_RETRIES=30
while [ $RETRIES -lt $MAX_RETRIES ]; do
    if curl -fsS http://localhost:4000/health > /dev/null 2>&1; then
        break
    fi
    RETRIES=$((RETRIES + 1))
    sleep 3
done

if [ $RETRIES -eq $MAX_RETRIES ]; then
    echo -e "${YELLOW}[WARN] Backend not healthy yet. It may still be starting (ML models take time to load).${NC}"
    echo "  Check logs: $COMPOSE_CMD -f $COMPOSE_FILE -f $VOLUMES_FILE logs -f backend"
else
    echo -e "${GREEN}[OK] Backend is healthy${NC}"
fi

# ---- Done ----
echo ""
echo -e "${GREEN}=========================================="
echo "  ${APP_NAME} is running!"
echo "=========================================="
echo ""
echo "  Frontend : http://localhost:3000"
echo "  Backend  : http://localhost:4000"
echo ""
echo "  Commands:"
echo "    Start:   ./start.sh"
echo "    Stop:    ./stop.sh"
echo "    Logs:    $COMPOSE_CMD -f $COMPOSE_FILE -f $VOLUMES_FILE logs -f"
echo ""
echo "  Si vous branchez un nouveau disque externe,"
echo "  relancez ./install.sh pour le detecter."
echo -e "==========================================${NC}"
