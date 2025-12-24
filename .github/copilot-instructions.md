# Copilot Instructions for annotation-nmr

## Project Architecture
- **Monorepo** with `frontend` (React, Create React App) and `backend` (Python, Node.js, ML, API) in top-level folders.
- **Frontend**: React app in `frontend/` (see `frontend/README.md` for scripts and build/test commands).
- **Backend**: Python and Node.js services in `backend/`. Entry points: `backend/server.js` (Node), `backend/py_entry.py` (Python ML logic), `backend/nmr_engine_entry.py` (NMR engine).
- **Build/Dist**: Build outputs in `build/`, `dist_final/`, and `recovered_dist/`.
- **Docker**: Multiple Dockerfiles and Compose files for local/dev/prod (`docker/`, `docker-compose*.yml`).
- **R Integration**: R-Portable in `R-Portable/` for R-based workflows.

## Developer Workflows
- **Frontend**: Use `npm start`, `npm run build`, `npm test` in `frontend/`.
- **Backend**:
  - Run Python ML: `python backend/py_entry.py` or via `run_python.bat`.
  - Node server: `node backend/server.js`.
  - Use `backend/build_bundle_complete.ps1` for Windows build automation.
- **Docker**: Use `docker-compose.dev.yml` for dev, `docker-compose.yml` for prod. Backend/Frontend containers defined in `docker/`.
- **Install/Setup**: Use scripts in `installation/` for environment setup (`install_windows.bat`, `setup_environment.py`).

## Conventions & Patterns
- **API**: Frontend talks to backend via REST endpoints (see `frontend/src/api.js`).
- **ML/Engine**: Python ML code is invoked from Node.js backend (see `backend/py_entry.py`, `backend/nmr_engine_entry.py`).
- **Config**: Central configs in `backend/config/` and `backend/constants/`.
- **Testing**: Test files prefixed with `test_` in root and backend (e.g., `test_ml_system.py`, `test_imports.py`).
- **Data**: Test data in CSV/SDF/JSON files in root and `test_output_*/`.

## Integration Points
- **Node ↔ Python**: Backend uses child process or similar to invoke Python scripts for ML/NMR tasks.
- **R**: R-Portable is bundled for reproducible R execution (see `R-Portable/`).
- **Docker**: All major services can be containerized; see Compose files for service boundaries.

## Examples
- To run the full stack locally: `docker-compose -f docker-compose.dev.yml up --build`
- To test ML system: `python test_ml_system.py` (from root)
- To build frontend: `cd frontend && npm run build`

## Key Files/Dirs
- `frontend/README.md`: Frontend dev scripts
- `backend/py_entry.py`, `backend/server.js`: Backend entry points
- `docker/`, `docker-compose*.yml`: Containerization
- `installation/`: Setup scripts
- `test_*.py`, `test_output_*/`: Tests and sample data

---
For more, see per-folder READMEs and comments in entry-point scripts.