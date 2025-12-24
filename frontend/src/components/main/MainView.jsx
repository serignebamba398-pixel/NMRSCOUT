import React, { useState, useEffect, useCallback, useMemo } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import Viewer from '../../Viewer';
import { normalizeStaticUrl } from '../../utils';
import './MainView.css';

/**
 * SIRIUS-style main view
 * Displays the list of annotated molecules on the left and the viewer on the right
 */
export const MainView = ({ onMoleculeSelect }) => {
  const { viewerData, customColors, setCustomColors } = useAppState();
  const [molecules, setMolecules] = useState([]);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [loading, setLoading] = useState(false);
  const [filterQuery, setFilterQuery] = useState('');

  const normalizedCsvs = useMemo(() => {
    const csvs = viewerData?.csvs;
    if (!csvs || typeof csvs !== 'object') return {};
    return Object.fromEntries(
      Object.entries(csvs).map(([k, v]) => [k, normalizeStaticUrl(v)])
    );
  }, [viewerData?.csvs]);

  const handleColorChange = useCallback(
    (method, color) => setCustomColors((prev) => ({ ...prev, [method]: color })),
    [setCustomColors]
  );

  const hasCsvs = useMemo(() => {
    return !!viewerData?.csvs && Object.keys(viewerData.csvs).length > 0;
  }, [viewerData]);

  const detectDelimiter = useCallback((headerLine) => {
    const line = String(headerLine ?? '');
    const counts = {
      ';': (line.match(/;/g) || []).length,
      ',': (line.match(/,/g) || []).length,
      '\t': (line.match(/\t/g) || []).length,
    };
    const best = Object.entries(counts).sort((a, b) => b[1] - a[1])[0];
    return best && best[1] > 0 ? best[0] : ';';
  }, []);

  const splitDelimitedLine = useCallback((line, delimiter) => {
    const s = String(line ?? '');
    const fields = [];
    let current = '';
    let inQuotes = false;

    for (let i = 0; i < s.length; i++) {
      const ch = s[i];

      if (ch === '"') {
        // Handle escaped double quote "" inside quoted field
        if (inQuotes && s[i + 1] === '"') {
          current += '"';
          i++;
          continue;
        }
        inQuotes = !inQuotes;
        continue;
      }

      if (!inQuotes && ch === delimiter) {
        fields.push(current);
        current = '';
        continue;
      }

      current += ch;
    }

    fields.push(current);
    return fields;
  }, []);

  const stripMethodPrefix = useCallback((name) => {
    if (!name) return name;
    return String(name).replace(/^(?:1H|13C|COSY|HSQC)[\s_\-:]+/i, '');
  }, []);

  // Filter molecules based on search
  const filteredMolecules = useMemo(() => {
    if (!filterQuery.trim()) return molecules;
    const q = filterQuery.toLowerCase();
    return molecules.filter(m => m.source.toLowerCase().includes(q));
  }, [molecules, filterQuery]);

  // Load molecules from annotation CSVs
  useEffect(() => {
    const loadMolecules = async () => {
      console.log('[MainView] === loadMolecules called ===');
      console.log('[MainView] viewerData:', viewerData);
      console.log('[MainView] csvs:', viewerData?.csvs);
      console.log('[MainView] csvs keys:', viewerData?.csvs ? Object.keys(viewerData.csvs) : 'none');
      
      if (!viewerData?.csvs || Object.keys(viewerData.csvs).length === 0) {
        console.log('[MainView] No CSVs available - exiting');
        setMolecules([]);
        return;
      }

      setLoading(true);
      try {
        const allMolecules = new Map();

        // Iterate through each CSV type (1H, 13C, COSY, HSQC)
        for (const [type, csvUrl] of Object.entries(viewerData.csvs)) {
          try {
            const url = normalizeStaticUrl(csvUrl);
            console.log(`[MainView] Fetching ${type} from:`, url);
            const response = await fetch(url);
            console.log(`[MainView] ${type} response status:`, response.status);
            if (!response.ok) {
              console.warn(`[MainView] Failed to fetch ${type}: ${response.status}`);
              continue;
            }
            
            const text = await response.text();
            console.log(`[MainView] ${type} CSV text length:`, text.length);
            
            // Remove UTF-8 BOM if present
            const cleanText = text.replace(/^\uFEFF/, '');
            
            const lines = cleanText.split(/\r?\n/).filter(l => l.trim());
            console.log(`[MainView] ${type} lines count:`, lines.length);
            if (lines.length < 2) {
              console.log(`[MainView] ${type} not enough lines, skipping`);
              continue;
            }

            const delimiter = detectDelimiter(lines[0]);

            // Nettoyer les headers (retirer guillemets)
            const headers = splitDelimitedLine(lines[0], delimiter).map(h => h.replace(/^["']|["']$/g, '').trim());
            console.log(`[MainView] ${type} headers:`, headers);
            
            // Chercher Source ou Filename (insensible à la casse)
            let sourceIdx = headers.findIndex(h => h.toLowerCase() === 'source');
            if (sourceIdx < 0) {
              sourceIdx = headers.findIndex(h => h.toLowerCase() === 'filename');
            }
            console.log(`[MainView] ${type} sourceIdx:`, sourceIdx);
            
            const lowerHeaders = headers.map(h => h.toLowerCase());
            // Chercher coverage_ransy (plusieurs variantes)
            const ransyCandidates = [
              'coverage_ransy',
              'coverage_ransy_mean',
              'coverage_ransy_avg',
              'coverage_ransy_percent',
            ];
            let ransyIdx = ransyCandidates
              .map(k => lowerHeaders.indexOf(k))
              .find(i => i >= 0);
            if (ransyIdx === undefined) ransyIdx = -1;
            if (ransyIdx < 0) {
              // fallback: any column that contains both "ransy" and "coverage"
              ransyIdx = lowerHeaders.findIndex(h => h.includes('ransy') && h.includes('coverage'));
            }

            // Chercher coverage_sub (plusieurs variantes)
            const subCandidates = [
              'coverage_sub',
              'coverage_sub_ratio',
              'coverage_sub_mean',
              'coverage_sub_avg',
            ];
            let subIdx = subCandidates
              .map(k => lowerHeaders.indexOf(k))
              .find(i => i >= 0);
            if (subIdx === undefined) subIdx = -1;
            if (subIdx < 0) {
              subIdx = lowerHeaders.findIndex(h => h.includes('sub') && h.includes('coverage'));
            }
            console.log(`[MainView] ${type} ransyIdx:`, ransyIdx, 'subIdx:', subIdx);

            for (let i = 1; i < lines.length; i++) {
              const values = splitDelimitedLine(lines[i], delimiter);
              // Nettoyer les quotes des valeurs
              const cleanValues = values.map(v => v?.replace(/^["']|["']$/g, '').trim());
              const source = sourceIdx >= 0 ? cleanValues[sourceIdx] : null;
              
              if (i === 1) {
                console.log(`[MainView] ${type} first row values:`, cleanValues);
                console.log(`[MainView] ${type} first row source:`, source);
              }
              
              if (!source) continue;

              // Nettoyer le nom: retirer _predicted, .sdf, etc.
              const cleanSource = stripMethodPrefix(source
                .replace(/_predicted$/i, '')
                .replace(/\.sdf$/i, '')
                .replace(/^NMR_/i, ''));
              
              if (!allMolecules.has(cleanSource)) {
                allMolecules.set(cleanSource, {
                  source: cleanSource,
                  originalSource: source,
                  types: {},
                  ransy: {},
                  sub: {},
                });
              }

              const mol = allMolecules.get(cleanSource);
              mol.types[type] = true;
              
              if (ransyIdx >= 0) {
                const raw = cleanValues[ransyIdx];
                const ransy = parseFloat(String(raw ?? '').replace(',', '.'));
                if (!isNaN(ransy)) mol.ransy[type] = ransy;
              }
              
              if (subIdx >= 0) {
                const raw = cleanValues[subIdx];
                const sub = parseFloat(String(raw ?? '').replace(',', '.'));
                if (!isNaN(sub)) mol.sub[type] = sub;
              }
            }
            
            console.log(`[MainView] ${type} molecules found:`, allMolecules.size);
          } catch (err) {
            console.warn(`Error loading CSV ${type}:`, err);
          }
        }

        // Convertir en tableau et calculer les moyennes
        const moleculeList = Array.from(allMolecules.values()).map(mol => ({
          ...mol,
          avgRansy: Object.values(mol.ransy).length > 0 
            ? Object.values(mol.ransy).reduce((a, b) => a + b, 0) / Object.values(mol.ransy).length 
            : 0,
          avgSub: Object.values(mol.sub).length > 0 
            ? Object.values(mol.sub).reduce((a, b) => a + b, 0) / Object.values(mol.sub).length 
            : 0,
        }));

        // Sort by descending RANSY average (stable tie-break)
        moleculeList.sort((a, b) => (b.avgRansy - a.avgRansy) || a.source.localeCompare(b.source));
        
        setMolecules(moleculeList);
        
        // Automatically select first molecule by default
        if (moleculeList.length > 0) {
          setSelectedMolecule(moleculeList[0]);
        }
      } catch (err) {
        console.error('Error loading molecules:', err);
      } finally {
        setLoading(false);
      }
    };

    loadMolecules();
  }, [viewerData, stripMethodPrefix, detectDelimiter, splitDelimitedLine]);

  const handleMoleculeClick = useCallback((mol) => {
    setSelectedMolecule(mol);
    onMoleculeSelect?.(mol);
  }, [onMoleculeSelect]);

  return (
    <div className="main-view">
      {/* Left panel - Molecules list */}
      <div className="molecule-list-panel">
        <div className="panel-header">
          <span className="panel-title">Annotated Molecules</span>
          <span className="molecule-count">{filteredMolecules.length}</span>
        </div>
        
        {/* SIRIUS-style search field */}
        <div className="filter-container">
          <input
            type="text"
            placeholder="Filter molecules..."
            value={filterQuery}
            onChange={(e) => setFilterQuery(e.target.value)}
            className="filter-input"
          />
          {filterQuery && (
            <button 
              className="filter-clear"
              type="button"
              onClick={() => setFilterQuery('')}
              title="Clear filter"
            >
              ×
            </button>
          )}
        </div>
        
        <div className="molecule-list">
          {loading ? (
            <div className="loading-state">Loading...</div>
          ) : filteredMolecules.length === 0 ? (
            <div className="empty-list">
              {molecules.length === 0 ? 'No molecules' : 'No results'}
            </div>
          ) : (
            filteredMolecules.map((mol, idx) => (
              <div
                key={mol.source}
                className={`molecule-item ${selectedMolecule?.source === mol.source ? 'selected' : ''}`}
                onClick={() => handleMoleculeClick(mol)}
              >
                <div className="molecule-info">
                  <div className="molecule-name">{mol.source}</div>
                </div>
              </div>
            ))
          )}
        </div>
      </div>

      {/* Panel droit - Viewer */}
      <div className="viewer-panel">
        {hasCsvs && selectedMolecule ? (
          <Viewer
            csvs={normalizedCsvs}
            taskId={viewerData.taskId}
            defaultShow3D={true}
            annotationResults={viewerData.annotationResults}
            pdataDirs={viewerData.pdataDirs || {}}
            bgOptions={viewerData.bgOptions || {}}
            externalSource={selectedMolecule.originalSource || selectedMolecule.source}
            hideSidebar={true}
            customColors={customColors}
            onColorChange={handleColorChange}
          />
        ) : (
          <Viewer
            csvs={{}}
            taskId={null}
            defaultShow3D={true}
            annotationResults={{}}
            pdataDirs={{}}
            bgOptions={{}}
            externalSource={null}
            hideSidebar={true}
            customColors={customColors}
            onColorChange={handleColorChange}
          />
        )}
      </div>
    </div>
  );
};

export default MainView;
