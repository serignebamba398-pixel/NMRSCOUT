# Build script for NMR Engine Bundle
# 1. Cleans previous build
# 2. Runs PyInstaller (which now handles most dependencies via spec)
# 3. Verifies the build

$ErrorActionPreference = "Stop"
$AnacondaPath = "C:\Users\seck-s\anaconda3"
$ProjectDir = "C:\Users\seck-s\Documents\_UNC\VERSIONDOCKER00Last\annotation-nmr\backend"
$AppName = "nmr_engine_production"
$DistDir = "$ProjectDir\dist\$AppName"

Set-Location $ProjectDir

Write-Host "=== Step 1: Cleaning previous build ===" -ForegroundColor Cyan
Remove-Item -Recurse -Force dist, build -ErrorAction SilentlyContinue

Write-Host "=== Step 2: Running PyInstaller ===" -ForegroundColor Cyan
& "$AnacondaPath\python.exe" -m PyInstaller --clean --noconfirm nmr_engine_production.spec
if ($LASTEXITCODE -ne 0) {
    Write-Error "PyInstaller failed with code $LASTEXITCODE"
}

Write-Host "=== Step 3: Post-Build Cleanup (Ensuring no bad DLLs) ===" -ForegroundColor Cyan
# Ensure bad DLLs are gone (redundant if spec works, but safe)
$BadPatterns = @("mkl_blacs_*.dll", "mkl_pgi_*.dll", "omptarget.*.dll")
foreach ($pattern in $BadPatterns) {
    Get-ChildItem -Path "$DistDir\_internal" -Filter $pattern -Recurse -ErrorAction SilentlyContinue | Remove-Item -Force -Verbose
}

Write-Host "=== Step 4: Verifying Build ===" -ForegroundColor Cyan
$ExePath = "$DistDir\$AppName.exe"
if (-not (Test-Path $ExePath)) {
    Write-Error "Executable not found at $ExePath"
}

Write-Host "Running test_imports..."
& $ExePath test_imports
if ($LASTEXITCODE -ne 0) {
    Write-Error "Verification failed! See output above."
}

Write-Host "=== Build Complete & Verified! ===" -ForegroundColor Green
Write-Host "Executable is at: $ExePath"
