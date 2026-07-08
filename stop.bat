@echo off
title NMRScout - Stop
echo.
echo  Stopping NMRScout...
echo.
docker compose -f docker-compose.shareable.yml -f docker-compose.volumes.yml down
echo.
echo  NMRScout stopped.
echo.
pause
