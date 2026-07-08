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

docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml up -d
echo.
echo  NMRScout: http://localhost:3000
echo.
pause
