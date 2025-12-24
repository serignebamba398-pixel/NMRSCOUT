@echo off
setlocal EnableExtensions EnableDelayedExpansion
chcp 65001 >nul

echo.
echo ============================================
echo   ANNOTATION-NMR - install_windows.bat
echo ============================================

REM -- Localisation des dossiers
set "SCRIPT_DIR=%~dp0"
for %%I in ("%SCRIPT_DIR%..") do set "ROOT=%%~fI"
set "BACKEND=%ROOT%\backend"
set "FRONTEND=%ROOT%\frontend"

REM -- Vérif Node
where node >nul 2>nul
if errorlevel 1 (
  echo [WARN] Node.js non detecte. Installe Node 18+ puis relance ce script.
) else (
  for /f "tokens=2 delims=v" %%v in ('node -v') do echo [OK] Node v%%v
)

REM -- Trouver conda
set "ACTIVATE_BAT="
if exist "%USERPROFILE%\anaconda3\Scripts\activate.bat" set "ACTIVATE_BAT=%USERPROFILE%\anaconda3\Scripts\activate.bat"
if exist "C:\ProgramData\anaconda3\Scripts\activate.bat" set "ACTIVATE_BAT=C:\ProgramData\anaconda3\Scripts\activate.bat"
if exist "C:\Users\seck-s\anaconda3\Scripts\activate.bat" set "ACTIVATE_BAT=C:\Users\seck-s\anaconda3\Scripts\activate.bat"

if not defined ACTIVATE_BAT (
  echo [ERR] Anaconda/Miniconda non trouve. Installe conda puis relance.
  goto END
)

REM -- Creer/MAJ l'env nmr
call "%ACTIVATE_BAT%"
for /f "delims=" %%e in ('conda env list ^| findstr /i "^nmr " 2^>nul') do set ENV_EXISTS=1

if not defined ENV_EXISTS (
  echo [INFO] Creation de l'environnement conda 'nmr' (python=3.9)...
  conda create -y -n nmr python=3.9
  if errorlevel 1 ( echo [ERR] Echec creation env. & goto END )
)

call "%ACTIVATE_BAT%" nmr
if errorlevel 1 ( echo [ERR] Echec activation env nmr. & goto END )

echo [INFO] Installation des dependances Python (CPU)...
conda install -y -c conda-forge rdkit pandas numpy scipy
conda install -y -c pytorch pytorch cpuonly
conda install -y -c dglteam dgl
pip install --upgrade pip
pip install dgllife tqdm
REM Si besoin d'autres libs: pip install networkx==3.1

REM -- Install deps Node backend
if exist "%BACKEND%\package.json" (
  echo [INFO] npm install (backend)...
  pushd "%BACKEND%"
  if exist package-lock.json (
    npm ci
  ) else (
    npm install
  )
  popd
) else (
  echo [WARN] package.json introuvable dans backend.
)

REM -- Creer .env backend (s'il n'existe pas)
if not exist "%BACKEND%\.env" (
  echo [INFO] Generation backend\.env ...
  >"%BACKEND%\.env" (
    echo PORT=4000
    echo CORS_ORIGINS=http://localhost:3000,http://127.0.0.1:3000
    echo FRONTEND_URL=http://localhost:3000
    echo BODY_LIMIT_MB=500
    echo MAX_UPLOAD_MB=500
    echo TIMEOUT_MS=900000
    echo TASK_TTL_HOURS=24
    echo PY_WRAPPER=%BACKEND:\=\\%\\run_python.bat
    REM Si tu veux forcer Python direct au lieu du wrapper:
    REM echo PYTHON=C:\Users\seck-s\anaconda3\envs\nmr\python.exe
    REM R_PATH optionnel si R n'est pas dans le PATH:
    REM echo R_PATH=C:\Program Files\R\R-4.3.1\bin\Rscript.exe
  )
) else (
  echo [OK] backend\.env deja present (non ecrase).
)

echo.
echo [OK] Installation terminee.
echo   - Place les poids .pt dans backend\models\
echo   - Lance le serveur :
echo       cd "%BACKEND%" ^&^& node server.js
echo   - Verifie: http://localhost:4000/health
echo.

:END
exit /b 0
