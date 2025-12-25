# annotation-nmr (NMRScout)

Desktop app (Electron) for NMR annotation.

## Repo structure

- `frontend/`: React (Create React App)
- `backend/`: Node.js API + Python ML/NMR engine entrypoints
- `docker/` + `docker-compose*.yml`: container workflows
- `R-Portable/`: bundled R runtime (used by annotation scripts; required at runtime)

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

## Releases (GitHub Actions)

- Pushing a tag like `v1.2.3` triggers `.github/workflows/release.yml`.
- The Release publishes installers for **Windows** (Setup + Portable `.exe`).
- The workflow requires **bundled R** which is not in Git; provide it via a repo secret:

| Secret name          | Contents                                               |
| -------------------- | ------------------------------------------------------ |
| `R_PORTABLE_ZIP_URL` | Archive of `R-Portable/` folder with `bin/Rscript.exe` |

> The Python engine is built in CI via PyInstaller. Alternatively, set `PYTHON_ENGINE_ZIP_URL` secret to skip PyInstaller build.
