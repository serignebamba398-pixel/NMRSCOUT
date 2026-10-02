@echo off
title NMRScout - Start
echo.
echo  Starting NMRScout...
echo.

if not exist docker-compose.volumes.yml (
    echo  First run detected, launching installer...
    powershell -ExecutionPolicy Bypass -File "%~dp0install.ps1"
    exit /b
)

echo  Pulling latest Docker images...
docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml pull
if errorlevel 1 (
    echo  Failed to pull Docker images. Check Docker Desktop and your internet connection.
    pause
    exit /b 1
)

docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml up -d
if errorlevel 1 (
    echo  Failed to start NMRScout.
    pause
    exit /b 1
)
echo.
echo  NMRScout: http://localhost:3000
echo.
pause
