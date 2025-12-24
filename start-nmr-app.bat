@echo off
setlocal enabledelayedexpansion

echo ===========================================
echo LANCEUR NMR APP v4 - Mode Debug Actif
echo ===========================================
echo.

REM Vérification Docker
echo  Verification de Docker...
docker --version >nul 2>&1
if errorlevel 1 (
    echo  Docker n'est pas installe ou non accessible dans le PATH
    pause & exit /b 1
)

docker info >nul 2>&1
if errorlevel 1 (
    echo  Docker n'est pas demarre. Veuillez lancer Docker Desktop.
    pause & exit /b 1
)
echo  Docker OK

REM Detection du disque en MINUSCULE
echo.
echo  Detection du disque...
set "DISK_LETTER=%CD:~0,1%"
set "DISK_LOWER=%DISK_LETTER%"

REM Conversion en minuscule (A->a, B->b, ..., Z->z)
for %%I in (a b c d e f g h i j k l m n o p q r s t u v w x y z) do (
    set "DISK_LOWER=!DISK_LOWER:%%~I=%%~I!"
)

echo  Disque detecte : %DISK_LETTER%: (montage: %DISK_LOWER%:)

REM Gestion du fichier .env avec VERIFICATION EXPLICITE
echo.
echo  Gestion du fichier .env...
echo Repertoire actuel : %CD%
echo Fichiers presents :
dir /b | findstr /i ".env"

if exist .env (
    echo  .env existe deja
) else if exist .env.example (
    echo   .env manquant, copie depuis .env.example...
    copy .env.example .env >nul
    if errorlevel 1 (
        echo  ERREUR lors de la copie
        pause & exit /b 1
    ) else (
        echo  Copie reussie
    )
) else (
    echo  ERREUR CRITIQUE : Ni .env ni .env.example n'existent !
    echo Assurez-vous que .env.example est bien present dans : %CD%
    pause & exit /b 1
)

REM Generation docker-compose.yml
echo.
echo  Generation du docker-compose.yml...

(
echo # Genere pour Windows le %date% %time%
echo services:
echo   backend:
echo     image: docker.io/serignebamba398/annotation-nmr-backend:v3
echo     container_name: nmr-backend
echo     restart: unless-stopped
echo     env_file:
echo       - .env
echo     volumes:
echo       - backend_uploads:/app/backend/uploads
echo       - backend_outputs:/app/backend/outputs
echo       - backend_cache:/app/backend/cache
echo       - %DISK_LOWER%:/:/host/%DISK_LOWER%:ro
echo     networks:
echo       - nmr-network
echo     ports:
echo       - "4000:4000"
echo     healthcheck:
echo       test: [ "CMD", "curl", "-fsS", "http://localhost:4000/health" ]
echo       interval: 15s
echo       timeout: 10s
echo       retries: 5
echo       start_period: 90s
echo.
echo   frontend:
echo     image: docker.io/serignebamba398/annotation-nmr-frontend:v3
echo     container_name: nmr-frontend
echo     restart: unless-stopped
echo     env_file:
echo       - .env
echo     ports:
echo       - "3000:80"
echo     depends_on:
echo       backend:
echo         condition: service_healthy
echo     networks:
echo       - nmr-network
echo     healthcheck:
echo       test: [ "CMD", "curl", "-fsS", "http://localhost/health" ]
echo       interval: 30s
echo       timeout: 5s
echo       retries: 3
echo.
echo volumes:
echo   backend_uploads:
echo     driver: local
echo     name: nmr_uploads
echo   backend_outputs:
echo     driver: local
echo     name: nmr_outputs
echo   backend_cache:
echo     driver: local
echo     name: nmr_cache
echo.
echo networks:
echo   nmr-network:
echo     driver: bridge
echo     name: nmr-network
) > docker-compose.yml

echo  docker-compose.yml cree

REM Pull et lancement
echo.
echo  Pull des images Docker...
docker compose pull

echo  Lancement des conteneurs...
docker compose up -d

echo.
echo  Attente du demarrage...
timeout /t 8 /nobreak >nul

REM Ouverture navigateur
echo  Ouverture du navigateur...
start http://localhost:3000

echo.
echo ===========================================
echo  APPLICATION NMR LANCEE !
echo ===========================================
echo.
echo  Chemins d'acces : /host/%DISK_LOWER%/
echo   IMPORTANT : Evitez les parentheses, espaces et caracteres speciaux
echo.
pause
