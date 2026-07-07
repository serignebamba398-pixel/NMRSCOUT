@echo off
title NMRScout - Installation
echo.
echo  Starting NMRScout installer...
echo.
powershell -ExecutionPolicy Bypass -File "%~dp0install.ps1"
echo.
pause
