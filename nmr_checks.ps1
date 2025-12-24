# ================================
# nmr_checks.ps1  (PowerShell)
# Smoke-tests Nginx + Backend /api + rewrites + endpoints clefs
# ================================

# --- PARAMS A AJUSTER SI BESOIN ---
$BaseUrl = "http://localhost"   # ex: http://localhost, http://127.0.0.1, ou https://mon-domaine
$TIMEOUT = 30

# (facultatif) chemins Bruker pour tester /spectral-bg/interactive avec un vrai dossier
$PDATA_1H = ""   # ex: "/host/C/Users/xxx/mon-projet/data/1H/pdata/1"
$PDATA_13C = ""
$PDATA_COSY = ""
$PDATA_HSQC = ""
# ----------------------------------

function Get-HeaderValue($resp, $name) {
  try { return ($resp.Headers[$name] | Select-Object -First 1) } catch { return $null }
}

function Show-Result($title, $resp) {
  $code = $null
  try { $code = $resp.StatusCode.value__ } catch { $code = $resp.StatusCode }

  $len = Get-HeaderValue $resp "Content-Length"
  $type = Get-HeaderValue $resp "Content-Type"
  $xBackend = Get-HeaderValue $resp "X-Backend-Route"
  $xCache = Get-HeaderValue $resp "X-Cache-Status"

  Write-Host "=== $title ===" -ForegroundColor Cyan
  if ($code) { Write-Host ("Status: {0}" -f $code) }
  if ($type) { Write-Host ("Content-Type: {0}" -f $type) }
  if ($len) { Write-Host ("Content-Length: {0}" -f $len) }
  if ($xBackend) { Write-Host ("X-Backend-Route: {0}" -f $xBackend) }
  if ($xCache) { Write-Host ("X-Cache-Status: {0}" -f $xCache) }
  try {
    $txt = $resp.Content
    if ($txt) {
      if ($txt.Length -gt 600) { $txt = $txt.Substring(0, 600) + "...(tronc.)" }
      Write-Host "Body:" -ForegroundColor DarkGray
      Write-Host $txt
    }
  }
  catch {}
  Write-Host ""
}

function Invoke-GET($url) {
  try {
    return Invoke-WebRequest -Uri $url -Method GET -TimeoutSec $TIMEOUT -UseBasicParsing
  }
  catch {
    if ($_.Exception.Response) { return $_.Exception.Response }
    throw
  }
}

function Invoke-POSTJSON($url, $obj) {
  $json = $obj | ConvertTo-Json -Depth 10
  try {
    return Invoke-WebRequest -Uri $url -Method POST -Body $json -ContentType "application/json" -TimeoutSec $TIMEOUT -UseBasicParsing
  }
  catch {
    if ($_.Exception.Response) { return $_.Exception.Response }
    throw
  }
}

Write-Host ("Hôte ciblé: {0}" -f $BaseUrl) -ForegroundColor Yellow
Write-Host ""

# 1) FRONT OK ?  /health (servi par Nginx / React)
$resp1 = Invoke-GET "$BaseUrl/health"
Show-Result "1) Frontend /health" $resp1

# 2) BACKEND OK ?  /api/health (proxy direct vers backend, pas de rewrite)
$resp2 = Invoke-GET "$BaseUrl/api/health"
Show-Result "2) Backend /api/health" $resp2

# 3) REWRITE /api/outputs/* -> /outputs/* (test côté Nginx)
#    On cible un chemin absent pour juste valider le passage et les entêtes
$testOutput = "$BaseUrl/api/outputs/__ping__/__not_found__.txt"
$resp3 = Invoke-GET $testOutput
Show-Result "3) Rewrite /api/outputs -> /outputs" $resp3

# 4) Endpoints clefs atteignables (sans tout lancer pour de vrai)
# 4a) /api/spectral-bg/interactive (JSON POST). Si pdata vide -> attendu: 4xx mais route OK.
$payloadInteractive = @{
  method    = "1H"
  pdata_dir = if ([string]::IsNullOrWhiteSpace($PDATA_1H)) { "__dummy__" } else { $PDATA_1H }
  size      = @{ width_px = 800; height_px = 300; dpi = 96 }
  window    = @{ ppm_min = 0; ppm_max = 10 }
}
$resp4a = Invoke-POSTJSON "$BaseUrl/api/spectral-bg/interactive" $payloadInteractive
Show-Result "4a) POST /api/spectral-bg/interactive (reachability)" $resp4a

# 4b) /api/predict/start (SSE côté backend).
#     On envoie un body minimal pour obtenir une réponse 2xx/4xx (mais pas 404), confirmant que la route est joignable.
$payloadPredict = @{
  mode    = "SMILES"
  smiles  = "CCO"
  methods = @("1H")
  options = @{ nb_cores = 1 }
}
$resp4b = Invoke-POSTJSON "$BaseUrl/api/predict/start" $payloadPredict
Show-Result "4b) POST /api/predict/start (reachability)" $resp4b

Write-Host "✅ Terminé. Si les statuts sont 200/204/4xx (mais PAS 404) sur 3) et 4), la chaîne Nginx→Backend est OK." -ForegroundColor Green
