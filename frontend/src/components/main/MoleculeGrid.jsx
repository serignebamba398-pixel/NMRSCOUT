import React, { useState, useEffect, useCallback, useRef } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import './MoleculeGrid.css';

// Chargeur RDKit
let rdkitPromise = null;
function loadRDKit() {
  if (rdkitPromise) return rdkitPromise;
  if (window.RDKitModule) return Promise.resolve(window.RDKitModule);

  rdkitPromise = new Promise((resolve, reject) => {
    const ex = document.querySelector('script[data-rdkit]');
    const onReady = () => {
      if (typeof window.initRDKitModule !== 'function') {
        reject(new Error('initRDKitModule non défini'));
        return;
      }
      window.initRDKitModule()
        .then((mod) => { window.RDKitModule = mod; resolve(mod); })
        .catch(reject);
    };

    if (ex) {
      ex.addEventListener('load', onReady, { once: true });
      return;
    }

    const s = document.createElement('script');
    s.async = true;
    s.dataset.rdkit = 'RDKit_minimal';
    s.src = 'https://unpkg.com/@rdkit/rdkit/dist/RDKit_minimal.js';
    s.onload = onReady;
    s.onerror = () => reject(new Error('Échec RDKit.js'));
    document.head.appendChild(s);
  });

  return rdkitPromise;
}

// API base URL
const API_BASE_URL = window.location.origin.includes('localhost') 
  ? 'http://localhost:3001' 
  : window.location.origin;

// Fetch helper
async function tryFetchText(url) {
  try {
    const r = await fetch(url);
    if (r.ok) return await r.text();
  } catch {}
  return null;
}

// Extraire le premier molblock d'un fichier SDF
function extractFirstMolBlock(text) {
  if (!text) return null;
  const endMarker = text.indexOf('M  END');
  if (endMarker < 0) return text;
  return text.substring(0, endMarker + 6);
}

/**
 * Composant MoleculeGrid - Affichage en grille style SIRIUS/MS-FINDER
 */
export const MoleculeGrid = ({ 
  molecules = [], 
  onMoleculeSelect, 
  selectedMolecule,
  onClose 
}) => {
  const { viewerData } = useAppState();
  const [rdkit, setRdkit] = useState(null);
  const [moleculesSvg, setMoleculesSvg] = useState({});
  const [loading, setLoading] = useState(true);
  const [downloadCount, setDownloadCount] = useState(50);
  const [columnsPerRow, setColumnsPerRow] = useState(5);
  const [showDownloadDialog, setShowDownloadDialog] = useState(false);
  const [moleculesData, setMoleculesData] = useState([]);
  const gridRef = useRef(null);

  // Charger RDKit
  useEffect(() => {
    loadRDKit()
      .then(mod => setRdkit(mod))
      .catch(err => console.error('Erreur chargement RDKit:', err))
      .finally(() => setLoading(false));
  }, []);

  // Charger les molblocks pour chaque molécule
  useEffect(() => {
    const loadMolblocks = async () => {
      if (!viewerData?.taskId || molecules.length === 0) {
        // Utiliser les molécules directement sans molblock
        setMoleculesData(molecules.map((mol, idx) => ({
          ...mol,
          rank: idx + 1,
          mw: mol.mw || 'N/A',
          score: mol.avgRansy?.toFixed(2) || 'N/A',
          deviation: mol.deviation || 'N/A',
          molBlock: null
        })));
        return;
      }

      const taskId = viewerData.taskId;
      const baseUrl = `${API_BASE_URL}/outputs/${taskId}/`;
      
      const enrichedMols = await Promise.all(
        molecules.map(async (mol, idx) => {
          const baseName = mol.source
            .replace(/_predicted$/i, '')
            .replace(/\.sdf$/i, '');
          
          // Essayer différentes variantes de noms de fichiers
          const candidates = [
            `${baseUrl}${baseName}.sdf`,
            `${baseUrl}${baseName}_predicted.sdf`,
            `${baseUrl}${encodeURIComponent(baseName)}.sdf`,
          ];
          
          let molBlock = null;
          for (const url of candidates) {
            const text = await tryFetchText(url);
            if (text) {
              molBlock = extractFirstMolBlock(text);
              break;
            }
          }
          
          return {
            ...mol,
            rank: idx + 1,
            mw: mol.mw || 'N/A',
            score: mol.avgRansy?.toFixed(2) || 'N/A',
            deviation: mol.avgDeviation?.toFixed(2) || 'N/A',
            molBlock
          };
        })
      );
      
      setMoleculesData(enrichedMols);
    };

    loadMolblocks();
  }, [viewerData, molecules]);

  // Générer les SVG des molécules
  useEffect(() => {
    if (!rdkit || moleculesData.length === 0) return;

    const generateSvgs = async () => {
      const svgs = {};
      
      for (const mol of moleculesData) {
        try {
          let rdkMol = null;
          
          if (mol.molBlock) {
            rdkMol = rdkit.get_mol(mol.molBlock);
          } else if (mol.smiles) {
            rdkMol = rdkit.get_mol(mol.smiles);
          }
          
          if (rdkMol) {
            const svg = rdkMol.get_svg(JSON.stringify({
              width: 150,
              height: 120,
              bondLineWidth: 1.2,
              atomLabelFontSize: 9,
              addAtomIndices: true,
              addStereoAnnotation: true,
              explicitMethyl: false
            }));
            svgs[mol.source] = svg;
            rdkMol.delete();
          }
        } catch (err) {
          console.warn(`Erreur génération SVG pour ${mol.source}:`, err);
        }
      }
      
      setMoleculesSvg(svgs);
    };

    generateSvgs();
  }, [rdkit, moleculesData]);

  // Téléchargement de la grille
  const handleDownload = useCallback(async (format = 'png') => {
    if (!gridRef.current) return;

    const count = Math.min(downloadCount, moleculesData.length);
    const modsToExport = moleculesData.slice(0, count);
    
    // Créer un canvas pour le rendu
    const cols = columnsPerRow;
    const rows = Math.ceil(count / cols);
    const cellWidth = 180;
    const cellHeight = 200;
    const padding = 10;
    
    const canvas = document.createElement('canvas');
    canvas.width = cols * cellWidth + padding * 2;
    canvas.height = rows * cellHeight + padding * 2;
    const ctx = canvas.getContext('2d');
    
    // Fond blanc
    ctx.fillStyle = 'white';
    ctx.fillRect(0, 0, canvas.width, canvas.height);
    
    // Dessiner chaque molécule
    for (let i = 0; i < modsToExport.length; i++) {
      const mol = modsToExport[i];
      const col = i % cols;
      const row = Math.floor(i / cols);
      const x = padding + col * cellWidth;
      const y = padding + row * cellHeight;
      
      // Cadre
      ctx.strokeStyle = '#e0e0e0';
      ctx.strokeRect(x, y, cellWidth - 4, cellHeight - 4);
      
      // SVG de la molécule
      const svg = moleculesSvg[mol.source];
      if (svg) {
        const img = new Image();
        const svgBlob = new Blob([svg], { type: 'image/svg+xml' });
        const url = URL.createObjectURL(svgBlob);
        
        await new Promise((resolve) => {
          img.onload = () => {
            ctx.drawImage(img, x + 10, y + 5, cellWidth - 24, 120);
            URL.revokeObjectURL(url);
            resolve();
          };
          img.onerror = resolve;
          img.src = url;
        });
      }
      
      // Informations texte
      ctx.fillStyle = '#333';
      ctx.font = '10px Arial';
      ctx.fillText(`Rank: ${mol.rank}  MW: ${mol.mw}`, x + 5, y + 135);
      ctx.fillText(mol.npaid || mol.source, x + 5, y + 148);
      ctx.fillText(`Score: ${mol.score}`, x + 5, y + 161);
      ctx.fillText(`Deviation: ${mol.deviation}`, x + 5, y + 174);
    }
    
    // Télécharger
    const link = document.createElement('a');
    link.download = `molecules_grid_${count}.${format}`;
    link.href = canvas.toDataURL(`image/${format}`, 0.95);
    link.click();
    
    setShowDownloadDialog(false);
  }, [downloadCount, columnsPerRow, moleculesData, moleculesSvg]);

  if (loading) {
    return (
      <div className="molecule-grid-loading">
        <div className="spinner"></div>
        <p>Chargement de RDKit...</p>
      </div>
    );
  }

  return (
    <div className="molecule-grid-container">
      {/* Header with controls */}
      <div className="molecule-grid-header">
        <div className="grid-header-left">
          <h3>Molecules Grid</h3>
          <span className="molecule-count-badge">{moleculesData.length} molecules</span>
        </div>
        <div className="grid-header-right">
          <label className="grid-columns-label">
            Columns:
            <select 
              value={columnsPerRow} 
              onChange={(e) => setColumnsPerRow(Number(e.target.value))}
              className="grid-columns-select"
            >
              <option value={3}>3</option>
              <option value={4}>4</option>
              <option value={5}>5</option>
              <option value={6}>6</option>
              <option value={8}>8</option>
            </select>
          </label>
          <button 
            className="grid-download-btn"
            onClick={() => setShowDownloadDialog(true)}
            title="Télécharger la grille"
          >
            <svg viewBox="0 0 24 24" width="16" height="16">
              <path fill="currentColor" d="M5 20h14v-2H5v2zm7-18v12l4-4h3L12 18 5 10h3l4 4V2z"/>
            </svg>
            Télécharger
          </button>
          {onClose && (
            <button className="grid-close-btn" onClick={onClose} title="Fermer">
              ✕
            </button>
          )}
        </div>
      </div>

      {/* Dialogue de téléchargement */}
      {showDownloadDialog && (
        <div className="download-dialog-overlay">
          <div className="download-dialog">
            <h4>Download Grid</h4>
            <div className="download-option">
              <label>Number of molecules:</label>
              <input 
                type="number" 
                min={1} 
                max={moleculesData.length}
                value={downloadCount}
                onChange={(e) => setDownloadCount(Math.min(Math.max(1, Number(e.target.value)), moleculesData.length))}
              />
              <span className="download-hint">/ {moleculesData.length} max</span>
            </div>
            <div className="download-option">
              <label>Colonnes par ligne:</label>
              <input 
                type="number" 
                min={2} 
                max={10}
                value={columnsPerRow}
                onChange={(e) => setColumnsPerRow(Math.min(Math.max(2, Number(e.target.value)), 10))}
              />
            </div>
            <div className="download-actions">
              <button onClick={() => handleDownload('png')} className="download-btn-png">
                📷 PNG
              </button>
              <button onClick={() => handleDownload('jpeg')} className="download-btn-jpg">
                🖼️ JPEG
              </button>
              <button onClick={() => setShowDownloadDialog(false)} className="download-btn-cancel">
                Annuler
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Grille des molécules */}
      <div 
        className="molecule-grid" 
        ref={gridRef}
        style={{ gridTemplateColumns: `repeat(${columnsPerRow}, 1fr)` }}
      >
        {moleculesData.map((mol) => (
          <div 
            key={mol.source}
            className={`molecule-grid-item ${selectedMolecule?.source === mol.source ? 'selected' : ''}`}
            onClick={() => onMoleculeSelect?.(mol)}
          >
            <div className="mol-svg-container">
              {moleculesSvg[mol.source] ? (
                <div 
                  className="mol-svg"
                  dangerouslySetInnerHTML={{ __html: moleculesSvg[mol.source] }}
                />
              ) : (
                <div className="mol-placeholder">
                  <span>🧪</span>
                </div>
              )}
            </div>
            <div className="mol-info">
              <div className="mol-info-row">
                <span className="mol-rank">Rank: {mol.rank}</span>
                <span className="mol-mw">MW: {mol.mw}</span>
              </div>
              <div className="mol-id">{mol.npaid || mol.source}</div>
              <div className="mol-info-row">
                <span className="mol-score">Score: {mol.score}</span>
              </div>
              <div className="mol-deviation">Deviation: {mol.deviation} ppm</div>
            </div>
          </div>
        ))}
      </div>
    </div>
  );
};

export default MoleculeGrid;
