# annotation-nmr (NMRScout)

Desktop app (Electron) for NMR annotation.

## Repo structure
- `frontend/`: React (Create React App)
- `backend/`: Node.js API + Python ML/NMR engine entrypoints
- `docker/` + `docker-compose*.yml`: container workflows
- `R-Portable/`: bundled R runtime (used by annotation scripts)

## Prerequisites
- Git LFS: `git lfs install` (required for large models/data)
- Node.js (for Electron + frontend/backend)
- Python (for ML/NMR scripts; see `requirements.txt` and `backend/requirements.cpu.txt`)

## Development
Run everything (backend + frontend + Electron):

```powershell
npm install
npm start
```

Frontend only:

```powershell
cd frontend
npm install
npm start
```

Backend only:

```powershell
cd backend
npm install
npm start
```

## Build / packaging
Build the Python bundle + frontend build, then package Electron:

```powershell
npm run build:production
```

Python bundle only:

```powershell
npm run build:python
```

## Notes
- Large files are tracked via Git LFS (`.gitattributes`).
- Outputs/caches/uploads are ignored by Git (`.gitignore`).
