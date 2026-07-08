# NMRScout

Application d'annotation et de prédiction RMN (NMR) pour les produits naturels.

Disponible en deux versions :
- **Docker** : application web accessible via navigateur (recommandé)
- **Electron** : application desktop Windows

## Installation Docker (recommandé)

Fonctionne sur **Windows, macOS (Intel + Apple Silicon), et Linux**.

### Prérequis

- [Docker Desktop](https://docs.docker.com/desktop/install/) installé et démarré

### Installation

**Windows** : double-cliquez sur `install.bat`

**macOS / Linux** :
```bash
chmod +x install.sh
./install.sh
```

L'application sera accessible sur **http://localhost:3000**

### Commandes utiles

```bash
# Arrêter
docker compose -f docker-compose.shareable.yml down

# Relancer
docker compose -f docker-compose.shareable.yml up -d

# Voir les logs
docker compose -f docker-compose.shareable.yml logs -f
```

## Installation Electron (Windows)

Téléchargez l'installateur depuis la page [Releases](https://github.com/serignebamba398-pixel/annotation-nmr/releases) :
- `NMRScout-vX.X.X-Setup.exe` : installateur classique
- `NMRScout-vX.X.X-portable.zip` : version portable (pas d'installation requise)

## Structure du projet

```
frontend/          React (Create React App)
backend/           Node.js API + Python ML/NMR engine
docker/            Dockerfiles + nginx.conf
scripts/           Scripts de build et déploiement
```

## Développement

### Prérequis dev

- Node.js 18+
- Python 3.10+ (avec les dépendances de `backend/requirements.cpu.txt`)
- R (pour les scripts d'annotation)

### Lancer en dev

```bash
# Tout (backend + frontend + Electron)
npm install
npm start

# Frontend seul
cd frontend && npm install && npm start

# Backend seul
cd backend && npm install && npm start
```

### Build Electron

```bash
npm run build:production
```

## CI/CD (GitHub Actions)

| Workflow | Déclencheur | Action |
|----------|-------------|--------|
| **CI** (`ci.yml`) | Push / PR | Tests, lint, build frontend |
| **Release** (`release.yml`) | Tag `v*` | Build Electron Windows + publication GitHub Release |
| **Docker Publish** (`docker-publish.yml`) | Manuel / Tag `v*` | Build multi-arch + push sur Docker Hub |

### Secrets GitHub requis

| Secret | Description |
|--------|-------------|
| `DOCKERHUB_USERNAME` | Nom d'utilisateur Docker Hub |
| `DOCKERHUB_TOKEN` | Token d'accès Docker Hub |
| `MODELS_ZIP_URL` | URL du zip contenant les modèles ML et données |
| `R_PORTABLE_ZIP_URL` | URL du zip R-Portable pour le build Electron |

### Fichiers volumineux

Les modèles ML et données (`.pt`, `.npz`, `.pkl`, `.sd`, `.jar`) ne sont pas dans le repo Git. Ils sont :
- Intégrés dans les **images Docker** sur Docker Hub
- Téléchargés depuis les **GitHub Release assets** lors des builds CI

## Licence

Projet privé - Tous droits réservés.
