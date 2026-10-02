@echo off
title NMRScout - Installation
echo.
echo  Starting NMRScout installer...
echo.
powershell -ExecutionPolicy Bypass -File "%~dp0install.ps1"
if errorlevel 1 (
	echo.
	echo  Installation failed. Review the errors above.
	pause
	exit /b 1
)
