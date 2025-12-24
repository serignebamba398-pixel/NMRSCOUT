# 🐳 NMR Annotation - Docker Deployment

## Quick Start (Utilisateurs)

### Prérequis
- Docker Desktop installé
- ~4 GB d'espace disque

### Lancement rapide

```bash
# 1. Télécharger les fichiers
curl -O https://raw.githubusercontent.com/YOUR_REPO/docker-compose.shareable.yml
curl -O https://raw.githubusercontent.com/YOUR_REPO/.env

# 2. Lancer l'application
docker compose -f docker-compose.shareable.yml up -d

# 3. Ouvrir dans le navigateur
# http://localhost:3000
```

### Arrêter l'application

```bash
docker compose -f docker-compose.shareable.yml down
```

---

## 🔧 Pour les Développeurs

### Build des images Docker v3

**Windows (PowerShell):**
```powershell
cd annotation-nmr
.\scripts\docker-build-push.ps1
```

**Linux/Mac:**
```bash
cd annotation-nmr
chmod +x scripts/docker-build-push.sh
./scripts/docker-build-push.sh
```

### Build manuel

```bash
# Backend
docker build -f docker/Dockerfile.backend -t serignebamba398/annotation-nmr-backend:v3 .

# Frontend
docker build -f docker/Dockerfile.frontend -t serignebamba398/annotation-nmr-frontend:v3 .

# Push
docker push serignebamba398/annotation-nmr-backend:v3
docker push serignebamba398/annotation-nmr-frontend:v3
```

### Mode développement

```bash
docker compose -f docker-compose.dev.yml up --build
```

---

## 📦 Images Docker

| Image | Description | Taille |
|-------|-------------|--------|
| `serignebamba398/annotation-nmr-backend:v3` | Backend Node.js + Python ML + R | ~4 GB |
| `serignebamba398/annotation-nmr-frontend:v3` | Frontend React + Nginx | ~50 MB |

### Contenu du Backend
- Node.js 18 LTS
- Python 3.10 (Conda) avec:
  - PyTorch 2.2.2 CPU
  - DGL 2.2.1
  - RDKit 2022.09.5
  - nmrglue, matplotlib
- R 4.x avec packages statistiques
- Modèles ML (1H et 13C)

---

## ⚙️ Configuration (.env)

| Variable | Description | Défaut |
|----------|-------------|--------|
| `FRONTEND_PORT` | Port frontend | 3000 |
| `BACKEND_PORT` | Port backend API | 4000 |
| `MAX_UPLOAD_MB` | Taille max upload | 500 |
| `MEMORY_LIMIT_MB` | Limite mémoire ML | 6144 |
| `OMP_NUM_THREADS` | Threads OpenMP | 4 |

---

## 🔍 Healthcheck

- **Backend:** `http://localhost:4000/health`
- **Frontend:** `http://localhost:3000`

---

## 📁 Volumes Docker

| Volume | Chemin container | Description |
|--------|------------------|-------------|
| `nmr_uploads` | `/app/backend/uploads` | Fichiers uploadés |
| `nmr_outputs` | `/app/backend/outputs` | Résultats |
| `nmr_cache` | `/app/backend/cache` | Cache |

---

## 🚀 Changelog

### v3 (2024-12-20)
- ✅ Fix encodage Unicode (caractères ✓ → [OK])
- ✅ Fix syntaxe R script annotation_13C.R
- ✅ Mise à jour architecture projet

### v2
- Initial Docker release
