const path = require('path');
const fs = require('fs');

// Determine scripts directory
let scriptsDir = path.join(__dirname, '..', 'scripts');

// In Electron production, scripts are in resources or unpacked asar
if (process.env.RESOURCES_PATH) {
  // Prioritize bundled Python scripts first (most likely in production)
  const candidates = [
    // Bundled Python: scripts are inside _internal folder (PRIORITY)
    path.join(process.env.RESOURCES_PATH, 'python', '_internal', 'scripts'),
    path.join(process.env.RESOURCES_PATH, 'app.asar.unpacked', 'backend', 'scripts'),
    path.join(process.env.RESOURCES_PATH, 'scripts'),
  ];
  
  // Check that the candidate has the main prediction script
  const mainScript = 'integrated_nmr_prediction.py';
  
  for (const candidate of candidates) {
    const mainScriptPath = path.join(candidate, mainScript);
    console.log('[paths] Checking candidate:', candidate, '-> exists:', fs.existsSync(mainScriptPath));
    if (fs.existsSync(mainScriptPath)) {
      scriptsDir = candidate;
      console.log('[paths] [OK] Using scriptsDir:', scriptsDir);
      break;
    }
  }
  
  // Fallback: if no candidate found, log warning
  if (!fs.existsSync(path.join(scriptsDir, mainScript))) {
    console.error('[paths] [WARN] WARNING: No valid scriptsDir found with', mainScript);
    console.error('[paths] Candidates checked:', candidates);
  }
}

let baseDir = path.join(__dirname, '..');
if (process.env.UPLOADS_DIR) {
  baseDir = path.dirname(process.env.UPLOADS_DIR);
}

const uploadsDir = process.env.UPLOADS_DIR || path.join(baseDir, 'uploads');
const outputsDir = process.env.OUTPUTS_DIR || path.join(baseDir, 'outputs');
const cacheDir = process.env.CACHE_DIR || path.join(baseDir, 'cache');
const projectsFile = process.env.PROJECTS_FILE || path.join(baseDir, 'projects.json');

// Folder creation
[uploadsDir, outputsDir, cacheDir].forEach(d => {
  if (d && !fs.existsSync(d)) {
    try {
      fs.mkdirSync(d, { recursive: true });
    } catch (e) {
      console.error(`[paths] Error creating dir ${d}:`, e.message);
    }
  }
});

// Init projects file if not exists
if (projectsFile && !fs.existsSync(projectsFile)) {
  try {
    fs.writeFileSync(projectsFile, JSON.stringify([], null, 2));
  } catch (e) {
    console.error(`[paths] Error creating projects file ${projectsFile}:`, e.message);
  }
}

// Scripts paths
const PRED_SCRIPT = path.join(scriptsDir, 'integrated_nmr_prediction.py');
const BG_SCRIPT = path.join(scriptsDir, 'nmr_export_bg.py');
const INTERACT_SCRIPT = path.join(scriptsDir, 'nmr_interact.py');

const EXTRACTORS = {
  '1H': path.join(scriptsDir, 'predict_1H_real.py'),
  '13C': path.join(scriptsDir, 'predict_13C_real.py'),
  COSY: path.join(scriptsDir, 'predict_COSY_real.py'),
  HSQC: path.join(scriptsDir, 'predict_HSQC_real.py'),
};

const scriptMap = {
  '1H': 'annotation_1H.R',
  '13C': 'annotation_13C.R',
  COSY: 'annotation_COSY.R',
  HSQC: 'annotation_HSQC.R',
};

const legacyOut = {
  '1H': '1H_With_Ransy_Julien.csv',
  '13C': '13C_With_Ransy_and_ECFP4.csv',
  COSY: 'COSY_With_Ransy_and_ECFP4.csv',
  HSQC: 'HSQC_With_Ransy_and_ECFP4.csv',
};

module.exports = {
  scriptsDir,
  uploadsDir,
  outputsDir,
  cacheDir,
  PRED_SCRIPT,
  BG_SCRIPT,
  INTERACT_SCRIPT,
  EXTRACTORS,
  scriptMap,
  legacyOut,
  projectsFile
};