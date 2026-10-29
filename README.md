# NMRScout

NMRScout is an NMR prediction and annotation application for natural products.

Two versions are available:

- **Docker**: a browser-based web application. Recommended for most users.
- **Electron**: a Windows desktop application.
## Tutorials

The following video tutorials guide users through the main NMRScout features, from installation to result visualization.

1. **[Part 1 – Download and run NMRScout](https://mediaserver.univ-nantes.fr/lti/v126d637166d4ooeprg7/)**

2. **[Part 2 – Prediction](https://mediaserver.univ-nantes.fr/lti/v126d637166dax0tgggb/)**

3. **[Part 3 – Annotation](https://mediaserver.univ-nantes.fr/lti/v126d637166e0ew7wvke/)**

4. **[Part 4 – Results visualization](https://mediaserver.univ-nantes.fr/edit/v126d637166cfzi37gw2/)**
## Docker Installation

The Docker version runs on Windows, macOS (Intel and Apple Silicon), and Linux.

### Requirements

- [Docker Desktop](https://docs.docker.com/desktop/install/) installed and running

### Install

1. Download the archive for your operating system:

| Operating system | Download                                                                                                                                    |
| ---------------- | ------------------------------------------------------------------------------------------------------------------------------------------- |
| Windows          | [NMRScout-Docker-Windows.zip](https://github.com/serignebamba398-pixel/annotation-nmr/releases/download/assets/NMRScout-Docker-Windows.zip) |
| macOS or Linux   | [NMRScout-Docker-Unix.zip](https://github.com/serignebamba398-pixel/annotation-nmr/releases/download/assets/NMRScout-Docker-Unix.zip)       |

2. Extract the archive to a folder of your choice.

3. Run the installer:
   - **Windows**: double-click `install.bat`.
   - **macOS or Linux**:
     ```bash
     chmod +x install.sh
     ./install.sh
     ```

The installer detects available internal and external drives, including USB drives, so NMRScout can access Bruker spectra.

After installation, NMRScout starts automatically and opens [http://localhost:3000](http://localhost:3000) in your default browser. No separate start or stop script is needed.

Projects are stored in `Documents/NMRScout/` on your computer.

### Archive contents

- **Windows**: `install.bat`, `install.ps1`, and Docker Compose configuration files
- **macOS and Linux**: `install.sh` and Docker Compose configuration files

### View logs

```bash
docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml logs -f
```

To detect a newly connected external drive, run the installer again (`install.bat` on Windows or `./install.sh` on macOS/Linux).

## Electron Installation for Windows

The Windows desktop installers are distributed from the [GitHub Releases page](https://github.com/serignebamba398-pixel/annotation-nmr/releases):

- `NMRScout-vX.X.X-Setup.exe`: standard Windows installer
- `NMRScout-vX.X.X-portable.zip`: portable version that does not require installation

The GitHub Actions release workflow builds and publishes the Electron application when a version tag matching `v*` is pushed.

## Project structure

```text
frontend/          React application (Create React App)
backend/           Node.js API and Python NMR engine
docker/            Dockerfiles and nginx configuration
scripts/           Build and deployment scripts
```

## Development

### Requirements

- Node.js 18 or later
- Python 3.10 or later, with dependencies from `backend/requirements.cpu.txt`
- R for the annotation scripts

### Run locally

Start the backend, frontend, and Electron application from the repository root:

```bash
npm install
npm start
```

To run one component, open a terminal in its directory:

```bash
cd frontend
npm install
npm start
```

```bash
cd backend
npm install
npm start
```

### Build Electron

```bash
npm run build:production
```

## GitHub Actions

| Workflow                              | Trigger                         | Action                                                      |
| ------------------------------------- | ------------------------------- | ----------------------------------------------------------- |
| CI (`ci.yml`)                         | Push or pull request            | Run tests, lint, and build the frontend                     |
| Release (`release.yml`)               | Tag matching `v*`               | Build the Windows Electron app and publish a GitHub Release |
| Docker Publish (`docker-publish.yml`) | Manual run or tag matching `v*` | Build multi-architecture images and push them to Docker Hub |

### Required GitHub secrets

| Secret               | Description                                            |
| -------------------- | ------------------------------------------------------ |
| `DOCKERHUB_USERNAME` | Docker Hub username                                    |
| `DOCKERHUB_TOKEN`    | Docker Hub access token                                |
| `MODELS_ZIP_URL`     | URL of the archive containing ML models and data       |
| `R_PORTABLE_ZIP_URL` | URL of the R-Portable archive used for Electron builds |

### Large files

ML models and data files (`.pt`, `.npz`, `.pkl`, `.sd`, and `.jar`) are not stored in the Git repository. They are included in the Docker images on Docker Hub or downloaded from GitHub Release assets during CI builds.

## License

Private project. All rights reserved.
