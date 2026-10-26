import { isMethodAvailable } from '../../constants/releaseFeatures';
import React, { useState, useEffect, useCallback, useMemo } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import Viewer from '../../Viewer';
import { normalizeStaticUrl } from '../../utils';
import './MainView.css';


export const MainView = ({ onMoleculeSelect }) => {
  const { viewerData, customColors, setCustomColors } = useAppState();
  const [moleculesRaw, setMoleculesRaw] = useState([]);
  const [selectedMolecule, setSelectedMolecule] = useState(null);
  const [loading, setLoading] = useState(false);
  const [filterQuery, setFilterQuery] = useState('');
  const [coverageMin, setCoverageMin] = useState('');
  const [coverageMax, setCoverageMax] = useState('');
  const [signalPerMethod, setSignalPerMethod] = useState({});
  const [signalModalOpen, setSignalModalOpen] = useState(false);
  const [filterModalOpen, setFilterModalOpen] = useState(false);
  const [selectedMethodsForCoverage, setSelectedMethodsForCoverage] = useState(["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable));

  const normalizedCsvs = useMemo(() => {
    const csvs = viewerData?.csvs;
    if (!csvs || typeof csvs !== 'object') return {};
    return Object.fromEntries(
      Object.entries(csvs).filter(([type]) => isMethodAvailable(type)).map(([k, v]) => [k, normalizeStaticUrl(v)])
    );
  }, [viewerData?.csvs]);

  const handleColorChange = useCallback(
    (method, color) => setCustomColors((prev) => ({ ...prev, [method]: color })),
    [setCustomColors]
  );

  const hasCsvs = useMemo(() => {
    return !!viewerData?.csvs && Object.keys(viewerData.csvs).some(isMethodAvailable);
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
    return String(name).replace(/^(?:1H|13C|COSY|TOCSY|HSQC)[\s_\-:]+/i, '');
  }, []);

  // Recompute avgRansy for the selected methods.
  const molecules = useMemo(() => {
    const methodsToUse = selectedMethodsForCoverage.length > 0
      ? selectedMethodsForCoverage
      : ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable);

    const list = moleculesRaw.map(mol => {
      let sumRansy = 0, countRansy = 0;
      let sumSub = 0, countSub = 0;

      for (const type of methodsToUse) {
        if (mol.ransy[type] != null) { sumRansy += mol.ransy[type]; countRansy++; }
        if (mol.sub[type] != null)   { sumSub   += mol.sub[type];   countSub++;   }
      }

      return {
        ...mol,
        avgRansy: countRansy > 0 ? sumRansy / countRansy : 0,
        avgSub:   countSub   > 0 ? sumSub   / countSub   : 0,
      };
    });

    list.sort((a, b) => (b.avgRansy - a.avgRansy) || a.source.localeCompare(b.source));
    return list;
  }, [moleculesRaw, selectedMethodsForCoverage]);

  // Filter molecules based on search, coverage range and signal count
  const filteredMolecules = useMemo(() => {
    let filtered = molecules;

    // Text filter
    if (filterQuery.trim()) {
      const q = filterQuery.toLowerCase();
      filtered = filtered.filter(m => m.source.toLowerCase().includes(q));
    }

    // Coverage filter (uses avgRansy as the coverage value)
    const minCov = parseFloat(coverageMin);
    const maxCov = parseFloat(coverageMax);
    if (!isNaN(minCov) || !isNaN(maxCov)) {
      filtered = filtered.filter(m => {
        const coverage = m.avgRansy ?? 0;
        if (!isNaN(minCov) && coverage < minCov) return false;
        if (!isNaN(maxCov) && coverage > maxCov) return false;
        return true;
      });
    }

    // Per-method signal count filter
    const activeSignalFilters = Object.entries(signalPerMethod).filter(([, v]) => v !== '' && !isNaN(parseInt(v, 10)));
    if (activeSignalFilters.length > 0) {
      filtered = filtered.filter(m => {
        for (const [method, minVal] of activeSignalFilters) {
          const min = parseInt(minVal, 10);
          const count = m.signals?.[method] ?? 0;
          if (count < min) return false;
        }
        return true;
      });
    }

    return filtered;
  }, [molecules, filterQuery, coverageMin, coverageMax, signalPerMethod]);

  // Load molecules from annotation CSVs
  useEffect(() => {
    const loadMolecules = async () => {
      console.log('[MainView] === loadMolecules called ===');
      console.log('[MainView] viewerData:', viewerData);
      console.log('[MainView] csvs:', viewerData?.csvs);
      console.log('[MainView] csvs keys:', viewerData?.csvs ? Object.keys(viewerData.csvs) : 'none');
      
      if (!viewerData?.csvs || Object.keys(viewerData.csvs).length === 0) {
        console.log('[MainView] No CSVs available - exiting');
        setMoleculesRaw([]);
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

            // Remove quotes from column headers.
            const headers = splitDelimitedLine(lines[0], delimiter).map(h => h.replace(/^["']|["']$/g, '').trim());
            console.log(`[MainView] ${type} headers:`, headers);
            
            // Match Source or Filename case-insensitively.
            let sourceIdx = headers.findIndex(h => h.toLowerCase() === 'source');
            if (sourceIdx < 0) {
              sourceIdx = headers.findIndex(h => h.toLowerCase() === 'filename');
            }
            console.log(`[MainView] ${type} sourceIdx:`, sourceIdx);
            
            const lowerHeaders = headers.map(h => h.toLowerCase());
            // Use coverage_ransy for all methods (1H, 13C, COSY, TOCSY, HSQC).
            // For 1H/13C, intersection_coverage is STOCSY-based and is often NA
            // when STOCSY is disabled or finds no intersections — coverage_ransy
            // is always populated (default 0) in all annotation modes.
            const ransyCandidates = [
              'coverage_ransy',
              'coverage_ransy_mean',
              'coverage_ransy_avg',
              'coverage_ransy_percent',
              'intersection_coverage', // fallback if coverage_ransy absent
            ];
            let ransyIdx = ransyCandidates
              .map(k => lowerHeaders.indexOf(k))
              .find(i => i >= 0) ?? -1;
            if (ransyIdx < 0) {
              // last resort: any column containing both "ransy" and "coverage"
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

            // ppm_molecule column: count number of signals
            const ppmMolIdx = lowerHeaders.indexOf('ppm_molecule');
            const is2D = type === 'COSY' || type === 'TOCSY' || type === 'HSQC';

            for (let i = 1; i < lines.length; i++) {
              const values = splitDelimitedLine(lines[i], delimiter);
              // Remove quotes from values.
              const cleanValues = values.map(v => v?.replace(/^["']|["']$/g, '').trim());
              const source = sourceIdx >= 0 ? cleanValues[sourceIdx] : null;
              
              if (i === 1) {
                console.log(`[MainView] ${type} first row values:`, cleanValues);
                console.log(`[MainView] ${type} first row source:`, source);
              }
              
              if (!source) continue;

              // Remove prediction suffixes and file extensions from the display name.
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
                  signals: {},
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

              // Count predicted signals from ppm_molecule
              if (ppmMolIdx >= 0) {
                const pmol = cleanValues[ppmMolIdx] ?? '';
                if (pmol && pmol !== 'NA') {
                  const count = is2D
                    ? pmol.split('||').filter(Boolean).length
                    : pmol.split('|').filter(Boolean).length;
                  if (count > 0) mol.signals[type] = count;
                }
              }
            }
            
            console.log(`[MainView] ${type} molecules found:`, allMolecules.size);
          } catch (err) {
            console.warn(`Error loading CSV ${type}:`, err);
          }
        }

        // Store raw RANSY and Sub metrics by spectrum type before applying weights.
        const rawList = Array.from(allMolecules.values());
        setMoleculesRaw(rawList);

        // Select the first molecule after the initial sort.
        if (rawList.length > 0) {
          // Initially sort by the mean score across available methods.
          const sorted = [...rawList].sort((a, b) => {
            const avgA = Object.values(a.ransy).reduce((s, v) => s + v, 0) / (Object.keys(a.ransy).length || 1);
            const avgB = Object.values(b.ransy).reduce((s, v) => s + v, 0) / (Object.keys(b.ransy).length || 1);
            return (avgB - avgA) || a.source.localeCompare(b.source);
          });
          setSelectedMolecule(sorted[0]);
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
        
        {/* Coverage range filter */}
        <div className="coverage-filter-container">
          <label className="coverage-filter-label">Coverage %</label>
          <div className="coverage-inputs">
            <input
              type="number"
              placeholder="Min"
              value={coverageMin}
              onChange={(e) => setCoverageMin(e.target.value)}
              className="coverage-input"
              min="0"
              max="100"
              step="1"
            />
            <span className="coverage-separator">-</span>
            <input
              type="number"
              placeholder="Max"
              value={coverageMax}
              onChange={(e) => setCoverageMax(e.target.value)}
              className="coverage-input"
              min="0"
              max="100"
              step="1"
            />
            {(coverageMin || coverageMax) && (
              <button 
                className="coverage-clear"
                type="button"
                onClick={() => { setCoverageMin(''); setCoverageMax(''); }}
                title="Clear coverage filter"
              >
                ×
              </button>
            )}
          </div>
        </div>
        
        {/* Signal count filter */}
        <div className="coverage-filter-container">
          <label className="coverage-filter-label">Nb signals</label>
          <div className="coverage-inputs">
            <button
              type="button"
              onClick={() => setSignalModalOpen(true)}
              style={{
                fontSize: 11, padding: "2px 10px", borderRadius: 4,
                border: `1px solid ${Object.values(signalPerMethod).some(v => v !== '') ? '#1d4ed8' : '#cbd5e1'}`,
                background: Object.values(signalPerMethod).some(v => v !== '') ? '#1d4ed8' : '#f8fafc',
                color: Object.values(signalPerMethod).some(v => v !== '') ? '#fff' : '#475569',
                cursor: 'pointer', fontWeight: 600,
              }}
            >
              {Object.values(signalPerMethod).filter(v => v !== '').length > 0
                ? `${Object.values(signalPerMethod).filter(v => v !== '').length} filter(s)`
                : 'Configure'}
            </button>
            {Object.values(signalPerMethod).some(v => v !== '') && (
              <button
                className="coverage-clear"
                type="button"
                onClick={() => setSignalPerMethod({})}
                title="Clear signal filter"
              >
                ×
              </button>
            )}
          </div>
        </div>

        {/* Method filter for molecule ranking. */}
        <div style={{ display: "flex", alignItems: "center", padding: "4px 8px 2px", gap: 6 }}>
          <span style={{ fontSize: 10, color: "#64748b", flex: 1 }}>
            {filteredMolecules.length} / {molecules.length}
          </span>
          <button
            onClick={() => setFilterModalOpen(true)}
            style={{
              fontSize: 10, padding: "2px 8px", borderRadius: 4,
              border: `1px solid #4a90d9`,
              background: selectedMethodsForCoverage.length < 4 ? "#4a90d9" : "#f8fafc",
              color: selectedMethodsForCoverage.length < 4 ? "#fff" : "#4a90d9",
              cursor: "pointer", fontWeight: 600,
            }}
            title="Choisir les méthodes pour le tri par coverage"
          >
            Sort by
          </button>
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
            coverageMethods={selectedMethodsForCoverage}
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

      {/* Method selection for coverage ranking. */}
      {filterModalOpen && (
        <div
          style={{
            position: "fixed", inset: 0, zIndex: 11000,
            background: "rgba(0,0,0,0.45)",
            display: "flex", alignItems: "center", justifyContent: "center",
          }}
          onClick={() => setFilterModalOpen(false)}
        >
          <div
            style={{
              background: "#fff", borderRadius: 12, padding: 24,
              width: 320, boxShadow: "0 8px 40px rgba(0,0,0,0.25)",
            }}
            onClick={e => e.stopPropagation()}
          >
            <div style={{ fontWeight: 700, fontSize: 14, marginBottom: 4, color: "#1e293b" }}>
              Sort by combined coverage
            </div>
            <div style={{ fontSize: 11, color: "#64748b", marginBottom: 16 }}>
              Select the methods to include in the calculation
            </div>
            {["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).map(method => {
              const hasData = moleculesRaw.some(m => m.ransy[method] != null);
              const checked = selectedMethodsForCoverage.includes(method);
              return (
                <label key={method} style={{
                  display: "flex", alignItems: "center", gap: 8,
                  padding: "7px 0", fontSize: 13,
                  color: hasData ? "#1e293b" : "#94a3b8",
                  cursor: hasData ? "pointer" : "not-allowed",
                  borderBottom: "1px solid #f1f5f9",
                }}>
                  <input
                    type="checkbox"
                    checked={checked}
                    disabled={!hasData}
                    onChange={() => {
                      setSelectedMethodsForCoverage(prev =>
                        prev.includes(method)
                          ? prev.filter(m => m !== method)
                          : [...prev, method]
                      );
                    }}
                    style={{ width: 14, height: 14 }}
                  />
                  <span style={{ fontWeight: 600 }}>{method}</span>
                  {!hasData && (
                    <span style={{ fontSize: 10, color: "#94a3b8" }}>(not available)</span>
                  )}
                </label>
              );
            })}
            <div style={{ marginTop: 16, display: "flex", gap: 8, justifyContent: "flex-end" }}>
              <button
                onClick={() => setSelectedMethodsForCoverage(["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable))}
                style={{
                  padding: "5px 14px", borderRadius: 6, fontSize: 12,
                  border: "1px solid #cbd5e1", background: "#f8fafc",
                  cursor: "pointer", color: "#475569",
                }}
              >
                Reset
              </button>
              <button
                onClick={() => setFilterModalOpen(false)}
                style={{
                  padding: "5px 14px", borderRadius: 6, fontSize: 12,
                  border: "1px solid #0ea5e9", background: "#0ea5e9",
                  cursor: "pointer", color: "#fff", fontWeight: 600,
                }}
              >
                Apply
              </button>
            </div>
          </div>
        </div>
      )}

      {/* Per-method signal count filter modal */}
      {signalModalOpen && (
        <div
          style={{
            position: "fixed", inset: 0, zIndex: 11000,
            background: "rgba(0,0,0,0.45)",
            display: "flex", alignItems: "center", justifyContent: "center",
          }}
          onClick={() => setSignalModalOpen(false)}
        >
          <div
            style={{
              background: "#fff", borderRadius: 12, padding: 24,
              width: 320, boxShadow: "0 8px 40px rgba(0,0,0,0.25)",
            }}
            onClick={e => e.stopPropagation()}
          >
            <div style={{ fontWeight: 700, fontSize: 14, marginBottom: 4, color: "#1e293b" }}>
              Minimum signals per method
            </div>
            <div style={{ fontSize: 11, color: "#64748b", marginBottom: 16 }}>
              A molecule must satisfy all active conditions to be displayed.
            </div>
            {["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).map(method => {
              const hasData = moleculesRaw.some(m => (m.signals?.[method] ?? 0) > 0);
              return (
                <div key={method} style={{
                  display: "flex", alignItems: "center", gap: 10,
                  padding: "8px 0", borderBottom: "1px solid #f1f5f9",
                }}>
                  <span style={{
                    fontWeight: 600, fontSize: 13, width: 60,
                    color: hasData ? "#1e293b" : "#94a3b8",
                  }}>{method}</span>
                  {!hasData && (
                    <span style={{ fontSize: 10, color: "#94a3b8", flex: 1 }}>(not available)</span>
                  )}
                  {hasData && (
                    <>
                      <span style={{ fontSize: 11, color: "#64748b" }}>min ≥</span>
                      <input
                        type="number"
                        min="0"
                        step="1"
                        placeholder="—"
                        value={signalPerMethod[method] ?? ''}
                        onChange={e => setSignalPerMethod(prev => ({
                          ...prev,
                          [method]: e.target.value,
                        }))}
                        style={{
                          width: 60, padding: "3px 6px", borderRadius: 4,
                          border: "1px solid #cbd5e1", fontSize: 12,
                          textAlign: "center",
                        }}
                      />
                    </>
                  )}
                </div>
              );
            })}
            <div style={{ marginTop: 16, display: "flex", gap: 8, justifyContent: "flex-end" }}>
              <button
                onClick={() => setSignalPerMethod({})}
                style={{
                  padding: "5px 14px", borderRadius: 6, fontSize: 12,
                  border: "1px solid #cbd5e1", background: "#f8fafc",
                  cursor: "pointer", color: "#475569",
                }}
              >
                Reset
              </button>
              <button
                onClick={() => setSignalModalOpen(false)}
                style={{
                  padding: "5px 14px", borderRadius: 6, fontSize: 12,
                  border: "1px solid #1d4ed8", background: "#1d4ed8",
                  cursor: "pointer", color: "#fff", fontWeight: 600,
                }}
              >
                Apply
              </button>
            </div>
          </div>
        </div>
      )}
    </div>
  );
};

export default MainView;
