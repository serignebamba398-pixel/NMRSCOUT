#!/usr/bin/env bash
set -euo pipefail

echo ""
echo "============================================"
echo "  ANNOTATION-NMR - install_linux_mac.sh"
echo "============================================"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
BACKEND="${ROOT}/backend"
FRONTEND="${ROOT}/frontend"

have_cmd() { command -v "$1" >/dev/null 2>&1; }

# --- Node ---
if have_cmd node; then
  echo "[OK] Node $(node -v)"
else
  echo "[WARN] Node.js non détecté. Installe Node 18+."
fi

# --- Python env: conda si dispo sinon venv ---
PY_BIN=""
if have_cmd conda; then
  echo "[INFO] Utilisation de conda"
  # Crée l'env nmr s'il n'existe pas
  if ! conda env list | grep -E '^[[:space:]]*nmr[[:space:]]' >/dev/null; then
    conda create -y -n nmr python=3.9
  fi
  # shellcheck disable=SC1091
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate nmr
  echo "[INFO] Installation deps Python (CPU) via conda..."
  conda install -y -c conda-forge rdkit pandas numpy scipy
  conda install -y -c pytorch pytorch cpuonly
  conda install -y -c dglteam dgl
  pip install --upgrade pip
  pip install dgllife tqdm
  PY_BIN="$(which python)"
else
  echo "[INFO] conda non disponible, utilisation venv Python"
  PYTHON3=$(command -v python3 || true)
  if [[ -z "${PYTHON3}" ]]; then
    echo "[ERR] python3 introuvable."; exit 1
  fi
  VENV_DIR="${BACKEND}/.venv"
  "${PYTHON3}" -m venv "${VENV_DIR}"
  # shellcheck disable=SC1091
  source "${VENV_DIR}/bin/activate"
  pip install --upgrade pip wheel
  # CPU only: torch via index CPU officiel
  pip install torch --index-url https://download.pytorch.org/whl/cpu
  # rdkit via pip (rdkit-pypi); si souci, préférer conda-forge
  pip install rdkit-pypi pandas numpy scipy dgl dgllife tqdm
  PY_BIN="$(which python)"
fi

# --- Node deps backend ---
if [[ -f "${BACKEND}/package.json" ]]; then
  echo "[INFO] npm install (backend)..."
  pushd "${BACKEND}" >/dev/null
  if [[ -f package-lock.json ]]; then
    npm ci
  else
    npm install
  fi
  popd >/dev/null
else
  echo "[WARN] package.json introuvable dans backend."
fi

# --- .env backend ---
if [[ ! -f "${BACKEND}/.env" ]]; then
  echo "[INFO] Génération backend/.env"
  cat > "${BACKEND}/.env" <<EOF
PORT=4000
CORS_ORIGINS=http://localhost:3000,http://127.0.0.1:3000
FRONTEND_URL=http://localhost:3000
BODY_LIMIT_MB=500
MAX_UPLOAD_MB=500
TIMEOUT_MS=900000
TASK_TTL_HOURS=24
PYTHON=${PY_BIN}
# R_PATH=/usr/bin/Rscript   # décommente si Rscript n'est pas dans le PATH
EOF
else
  echo "[OK] backend/.env déjà présent (non écrasé)."
fi

echo ""
echo "[OK] Installation terminée."
echo "  - Place les poids .pt dans backend/models/"
echo "  - Lance le serveur :"
echo "      cd \"${BACKEND}\" && node server.js"
echo "  - Vérifie: http://localhost:4000/health"
echo ""
