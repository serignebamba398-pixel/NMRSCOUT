@echo off
setlocal EnableExtensions EnableDelayedExpansion
chcp 65001 >nul

REM === Directories ===
set "SCRIPT_DIR=%~dp0"
set "ENV_NAME=nmr"

REM === Search for conda activate.bat ===
set "ACTIVATE_BAT="
if exist "%USERPROFILE%\anaconda3\Scripts\activate.bat" set "ACTIVATE_BAT=%USERPROFILE%\anaconda3\Scripts\activate.bat"
if exist "C:\ProgramData\anaconda3\Scripts\activate.bat" set "ACTIVATE_BAT=C:\ProgramData\anaconda3\Scripts\activate.bat"
if exist "C:\Users\seck-s\anaconda3\Scripts\activate.bat" set "ACTIVATE_BAT=C:\Users\seck-s\anaconda3\Scripts\activate.bat"

if defined ACTIVATE_BAT (
  call "%ACTIVATE_BAT%" %ENV_NAME%
  if errorlevel 1 (
    echo [WARN] Failed conda env "%ENV_NAME%" activation, trying system python...
    goto RUN_DIRECT
  )
) else (
  echo [WARN] activate.bat introuvable, on tente python systeme...
  goto RUN_DIRECT
)

REM === Env CPU + Pythonpath ===
set "PYTHONUNBUFFERED=1"
set "CUDA_VISIBLE_DEVICES=-1"
set "DGLBACKEND=pytorch"
set "PYTHONPATH=%SCRIPT_DIR%ml_dependencies"

python %*
exit /b %ERRORLEVEL%

:RUN_DIRECT
set "PYTHONUNBUFFERED=1"
set "CUDA_VISIBLE_DEVICES=-1"
set "DGLBACKEND=pytorch"
set "PYTHONPATH=%SCRIPT_DIR%ml_dependencies"
python %*
exit /b %ERRORLEVEL%
