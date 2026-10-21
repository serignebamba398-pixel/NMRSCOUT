import { isMethodAvailable } from '../../constants/releaseFeatures';
import { TocsyPredictionFields } from './TocsyPredictionFields';
import React, { useState, useEffect, useCallback } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import { LockOverlay, ProgressBar, FileUpload } from '../common';
import { formatFileSize, formatDuration } from '../../utils';
import './TabbedPredictionInterface.css';

/**
 * Tabbed prediction interface MS-FINDER style
 */
export const TabbedPredictionInterface = ({ inModal = false, onClose }) => {
  const {
    predictionState,
    setPredictionState,
    uiLock,
    setUiLock,
    addLog,
    // Use persistent prediction functions from context
    prediction,
  } = useAppState();

  const {
    running,
    progress,
    progressMessage,
    error,
    startTime,
    endTime,
    totalMolecules,
    completedMolecules,
    csvUploads,
    loadingCsvUploads,
    refreshCsvUploads,
    uploadCsv,
    loadMoleculeCountFromServerCsv,
    launchSmiles,
    launchCSV,
    handleCancel,
  } = prediction;

  const {
    mode,
    smiles,
    molId,
    csvLocalFile,
    csvUploadedName,
    options,
    selectedPredMethods,
    artifacts,
  } = predictionState;

  // Active tab: 'source' | '1H' | '13C' | 'COSY' | 'TOCSY' | 'HSQC'
  const [activeTab, setActiveTab] = useState('source');

  const dimmed = uiLock.lockedBy === "ANNO" || uiLock.lockedBy === "OPTIM";

  // Extracted parameters
  const ex = options?.extractors || {};
  const h1 = ex.h1 || {};
  const c13 = ex.c13 || {};
  const hsqc = ex.hsqc || {};
  const keepExchangeableUI = !!h1.keep_exchangeable;
  const carbonModeUI = c13.carbon_mode || "all";
  const hsqcFoldUI = !!hsqc.fold;
  const hsqcLowUI = Number.isFinite(+hsqc.fold_low) ? +hsqc.fold_low : 4.8618;
  const hsqcHighUI = Number.isFinite(+hsqc.fold_high) ? +hsqc.fold_high : 69.8395;

  // Lock UI during execution
  useEffect(() => {
    if (running) {
      setUiLock({ lockedBy: "PRED", message: "Prediction in progress..." });
    } else if (uiLock.lockedBy === "PRED") {
      setUiLock({ lockedBy: null, message: "" });
    }
  }, [running, setUiLock, uiLock.lockedBy]);

  // Load CSVs on mount if CSV mode
  useEffect(() => {
    if (mode === "CSV" && csvUploads.length === 0) {
      refreshCsvUploads();
    }
  }, [mode, csvUploads.length, refreshCsvUploads]);

  const updateOption = useCallback((key, value) => {
    setPredictionState((prev) => ({
      ...prev,
      options: { ...prev.options, [key]: value },
    }));
  }, [setPredictionState]);

  const updateCosyOption = useCallback((key, value) => {
    setPredictionState((prev) => ({
      ...prev,
      options: {
        ...prev.options,
        cosy: { ...(prev.options.cosy || {}), [key]: value },
      },
    }));
  }, [setPredictionState]);

  const updateTocsyOption = useCallback((key, value) => {
    setPredictionState((prev) => ({
      ...prev,
      options: {
        ...prev.options,
        tocsy: { ...(prev.options.tocsy || {}), [key]: value },
      },
    }));
  }, [setPredictionState]);

  const updateExtractorOption = useCallback((extractor, key, value) => {
    setPredictionState((prev) => ({
      ...prev,
      options: {
        ...prev.options,
        extractors: {
          ...(prev.options?.extractors || {}),
          [extractor]: {
            ...(prev.options?.extractors?.[extractor] || {}),
            [key]: value,
          },
        },
      },
    }));
  }, [setPredictionState]);

  const toggleMethod = useCallback((m) => {
    if (running) return;
    setPredictionState((prev) => {
      const set = new Set(prev.selectedPredMethods);
      if (set.has(m)) set.delete(m);
      else set.add(m);
      return { ...prev, selectedPredMethods: Array.from(set) };
    });
  }, [running, setPredictionState]);

  // Handler to launch prediction AND close modal
  const handleLaunchAndClose = useCallback(() => {
    const launcher = mode === "SMILES" ? launchSmiles : launchCSV;
    launcher();
    // Close modal after launch
    if (inModal && onClose) {
      onClose();
    }
  }, [mode, launchSmiles, launchCSV, inModal, onClose]);

  // Number of selected methods
  const methodCount = selectedPredMethods.length;

  // Source & Methods tab
  // Experimental parameters (for prediction-time JSON generation)
  const freqMhz = options?.freq_mhz ?? 700;
  const solvent = options?.solvent ?? 'CDCl3';
  const temperatureK = options?.temperature_k ?? 298;
  const pH = options?.ph ?? 7.4;

  const renderSourceTab = () => (
    <div className="tab-content">
      {/* Mode selector */}
      <div className="settings-section">
        <h4 className="section-title">Input Mode</h4>
        <div className="mode-selector">
          {["SMILES", "CSV"].map((m) => (
            <button
              key={m}
              onClick={() => !running && setPredictionState((prev) => ({ ...prev, mode: m }))}
              className={`mode-btn ${mode === m ? "active" : ""}`}
              disabled={running}
            >
              {m}
            </button>
          ))}
        </div>
      </div>

      {/* SMILES Mode */}
      {mode === "SMILES" && (
        <div className="settings-section">
          <h4 className="section-title">SMILES Molecule</h4>
          <div className="form-grid">
            <div className="form-field">
              <label>SMILES</label>
              <input
                type="text"
                value={smiles}
                onChange={(e) => setPredictionState((prev) => ({ ...prev, smiles: e.target.value }))}
                disabled={running}
                placeholder="E.g.: CCO"
              />
            </div>
            <div className="form-field">
              <label>Molecule ID</label>
              <input
                type="text"
                value={molId}
                onChange={(e) => setPredictionState((prev) => ({ ...prev, molId: e.target.value }))}
                disabled={running}
                placeholder="ethanol"
              />
            </div>
          </div>
        </div>
      )}

      {/* CSV Mode */}
      {mode === "CSV" && (
        <div className="settings-section">
          <h4 className="section-title">CSV File</h4>
          
          <div className="upload-area">
            <FileUpload
              label="CSV file (columns: id, smiles)"
              file={csvLocalFile}
              onChange={(e) => uploadCsv(e.target.files?.[0])}
              accept=".csv"
              disabled={running}
            />
          </div>

          <div className="csv-selector">
            <label>Or select an existing CSV</label>
            <div className="selector-row">
              <select
                disabled={running || loadingCsvUploads}
                value={csvUploadedName || ""}
                onChange={async (e) => {
                  const name = e.target.value || "";
                  setPredictionState((prev) => ({
                    ...prev,
                    csvUploadedName: name,
                    csvLocalFile: null,
                  }));
                  if (name) {
                    addLog(`CSV selected: ${name}`);
                    await loadMoleculeCountFromServerCsv(name);
                  }
                }}
              >
                <option value="">-- None --</option>
                {csvUploads.map((f) => (
                  <option key={f.filename} value={f.filename}>
                    {f.filename} ({formatFileSize(f.size)})
                  </option>
                ))}
              </select>
              <button
                className="btn-icon"
                onClick={refreshCsvUploads}
                disabled={running || loadingCsvUploads}
                title="Refresh"
              >
                <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                  <path d="M23 4v6h-6M1 20v-6h6M3.51 9a9 9 0 0114.85-3.36L23 10M1 14l4.64 4.36A9 9 0 0020.49 15"/>
                </svg>
              </button>
            </div>
          </div>

          {csvUploadedName && (
            <div className="csv-status success">
              CSV selected: {csvUploadedName}
              {totalMolecules ? ` (${totalMolecules} molecules)` : ''}
            </div>
          )}

          <div className="form-grid" style={{ marginTop: 12 }}>
            <div className="form-field">
              <label>ID Column</label>
              <input
                type="text"
                value={options.idColumn || 'id'}
                onChange={(e) => updateOption('idColumn', e.target.value)}
                disabled={running}
              />
            </div>
            <div className="form-field">
              <label>SMILES Column</label>
              <input
                type="text"
                value={options.smilesColumn || 'smiles'}
                onChange={(e) => updateOption('smilesColumn', e.target.value)}
                disabled={running}
              />
            </div>
          </div>

          <p className="hint" style={{ marginTop: 10 }}>
            Experimental Conditions below apply to this CSV prediction and are stored in the generated spectra metadata.
          </p>

          {/* Incremental prediction - resume from existing task */}
          <div className="settings-section" style={{ marginTop: 16, padding: '12px', backgroundColor: '#f0f7ff', borderRadius: '6px', border: '1px solid #cce0ff' }}>
            <h4 className="section-title" style={{ fontSize: '13px', marginBottom: 8, color: '#1a73e8' }}>
              Incremental Prediction (optional)
            </h4>
            <p className="hint" style={{ marginBottom: 8, fontSize: '11px' }}>
              Resume from a previous task: molecules already predicted with the same parameters (freq, solvent, T, pH) will be reused.
            </p>
            <div className="form-field">
              <label>Previous Task Folder</label>
              <input
                type="text"
                value={options.resumeFromDir || ''}
                onChange={(e) => updateOption('resumeFromDir', e.target.value)}
                disabled={running}
                placeholder="e.g., C:\...\outputs\22c5f5d9-ee0c-4a5f-..."
                style={{ fontFamily: 'monospace', fontSize: '11px' }}
              />
              <p className="hint" style={{ marginTop: 4 }}>
                Paste the full path to a previous output folder containing checkpoint.json and SDF files.
              </p>
            </div>
            {options.resumeFromDir && (
              <div className="csv-status success" style={{ marginTop: 8 }}>
                Incremental mode: will reuse molecules from previous task
              </div>
            )}
          </div>
        </div>
      )}

      {/* Prediction methods */}
      <div className="settings-section">
        <h4 className="section-title">Prediction Methods</h4>
        <p className="section-hint">Check the methods then configure each tab if needed.</p>
        <div className="methods-list">
          { ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).map((m) => {
            const checked = selectedPredMethods.includes(m);
            const hasParams = m === "13C" || m === "COSY" || m === "TOCSY" || m === "HSQC";
            
            return (
              <div key={m} className={`method-row ${checked ? 'selected' : ''}`}>
                <label className="method-checkbox">
                  <input
                    type="checkbox"
                    checked={checked}
                    onChange={() => toggleMethod(m)}
                    disabled={running}
                  />
                  <span className="method-name">{m}</span>
                </label>
                <div className="method-info">
                  {hasParams && (
                    <span className="has-params">Parameters available</span>
                  )}
                </div>
                {hasParams && (
                  <button
                    className="btn-configure"
                    onClick={() => setActiveTab(m)}
                    disabled={running}
                    title={`Configure ${m}`}
                  >
                    Configure
                  </button>
                )}
              </div>
            );
          })}
        </div>

        {/* nb_cores global */}
        <div className="form-field" style={{ marginTop: 12 }}>
          <label>Number of cores (nb_cores)</label>
          <input
            type="number"
            value={options.nb_cores || 4}
            onChange={(e) => updateOption('nb_cores', parseInt(e.target.value) || 4)}
            disabled={running}
            min={1}
            max={32}
            style={{ width: 80 }}
          />
        </div>
      </div>

      <div className="settings-section">
        <h4 className="section-title" style={{ fontSize: '13px', marginTop: 16 }}>Experimental Conditions (SMILES and CSV)</h4>
        <p className="hint" style={{ marginBottom: 12 }}>
          These parameters are used during prediction to generate the multiplicity JSON with accurate shifts and J-couplings.
        </p>

        <div className="form-grid">
          <div className="form-field">
            <label>Spectrometer Frequency (MHz)</label>
            <input
              type="number"
              value={freqMhz}
              onChange={(e) => updateOption('freq_mhz', parseFloat(e.target.value) || 700)}
              disabled={running}
              min={100}
              max={1200}
              step={50}
              style={{ width: 100 }}
            />
          </div>

          <div className="form-field">
            <label>Solvent</label>
            <select
              value={solvent}
              onChange={(e) => updateOption('solvent', e.target.value)}
              disabled={running}
            >
              <option value="CDCl3">CDCl3 (Chloroform-d)</option>
              <option value="DMSO-d6">DMSO-d6 (Dimethylsulfoxide-d6)</option>
              <option value="D2O">D2O (Deuterium oxide)</option>
              <option value="CD3OD">CD3OD (Methanol-d4)</option>
              <option value="acetone-d6">Acetone-d6</option>
              <option value="C6D6">C6D6 (Benzene-d6)</option>
              <option value="pyridine-d5">Pyridine-d5</option>
              <option value="CD3CN">CD3CN (Acetonitrile-d3)</option>
              <option value="THF-d8">THF-d8</option>
            </select>
            <p className="hint">
              Exchangeable protons (OH, NH, SH) will be included/excluded based on the solvent.
              In protic solvents (D2O, CD3OD), they exchange rapidly and are typically not visible.
            </p>
          </div>
        </div>

        <div className="form-grid">
          <div className="form-field">
            <label>Temperature (K)</label>
            <input
              type="number"
              value={temperatureK}
              onChange={(e) => updateOption('temperature_k', parseFloat(e.target.value) || 298)}
              disabled={running}
              min={200}
              max={400}
              step={1}
              style={{ width: 80 }}
            />
          </div>

          <div className="form-field">
            <label>pH</label>
            <input
              type="number"
              value={pH}
              onChange={(e) => updateOption('ph', parseFloat(e.target.value) || 7.4)}
              disabled={running}
              min={0}
              max={14}
              step={0.1}
              style={{ width: 80 }}
            />
          </div>
        </div>
      </div>

      {/* Summary */}
      <div className="methods-summary">
        {methodCount === 0 ? (
          <span className="summary-warning">Select at least one method</span>
        ) : (
          <span className="summary-ready">{methodCount} method(s) selected</span>
        )}
      </div>
    </div>
  );

  // 1H Tab
  const render1HTab = () => {

    return (
      <div className="tab-content">
        <div className="settings-section">
          <h4 className="section-title">1H NMR Parameters</h4>

          <div className="method-status-box">
            <span className={selectedPredMethods.includes("1H") ? "enabled" : "disabled"}>
              {selectedPredMethods.includes("1H") ? "Method enabled" : "Method disabled"}
            </span>
          </div>

          <p className="hint">Experimental conditions are configured in Source &amp; Methods.</p>
        </div>
      </div>
    );
  };

  // 13C Tab
  const render13CTab = () => (
    <div className="tab-content">
      <div className="settings-section">
        <h4 className="section-title">13C NMR Parameters</h4>
        
        <div className="method-status-box">
          <span className={selectedPredMethods.includes("13C") ? "enabled" : "disabled"}>
            {selectedPredMethods.includes("13C") ? "Method enabled" : "Method disabled"}
          </span>
        </div>

        <div className="form-field">
          <label>Carbons to predict</label>
          <select
            value={carbonModeUI}
            onChange={(e) => updateExtractorOption('c13', 'carbon_mode', 
              e.target.value.toLowerCase().startsWith("all") ? "all" : "nonquaternary"
            )}
            disabled={running}
          >
            <option value="all">All carbons (C, CH, CH2, CH3)</option>
            <option value="nonquaternary">Non-quaternary carbons (CH/CH2/CH3)</option>
          </select>
          <p className="hint">
            The "non-quaternary" mode excludes carbons without attached hydrogen, 
            which is often more suitable for HSQC.
          </p>
        </div>
      </div>
    </div>
  );

  // COSY Tab
  const renderCOSYTab = () => (
    <div className="tab-content">
      <div className="settings-section">
        <h4 className="section-title">COSY Parameters</h4>
        
        <div className="method-status-box">
          <span className={selectedPredMethods.includes("COSY") ? "enabled" : "disabled"}>
            {selectedPredMethods.includes("COSY") ? "Method enabled" : "Method disabled"}
          </span>
        </div>

        <div className="params-grid">
          <div className="param-row">
            <label>min_bonds:</label>
            <input
              type="number"
              value={options.cosy?.min_bonds || 2}
              onChange={(e) => updateCosyOption('min_bonds', parseInt(e.target.value) || 2)}
              disabled={running}
              min={1}
              max={10}
            />
          </div>
          <div className="param-row">
            <label>max_bonds:</label>
            <input
              type="number"
              value={options.cosy?.max_bonds || 3}
              onChange={(e) => updateCosyOption('max_bonds', parseInt(e.target.value) || 3)}
              disabled={running}
              min={1}
              max={10}
            />
          </div>
        </div>
        <p className="hint">
          The number of bonds determines which H-H correlations are predicted.
          Typically 2-3 for standard COSY.{isMethodAvailable('TOCSY') && ' TOCSY uses its own spin-system settings.'}
        </p>
      </div>
    </div>
  );

  const renderTOCSYTab = () => (
    <div className="tab-content">
      <div className="settings-section">
        <h4 className="section-title">TOCSY Parameters</h4>
        <div className="method-status-box">
          <span className={selectedPredMethods.includes("TOCSY") ? "enabled" : "disabled"}>
            {selectedPredMethods.includes("TOCSY") ? "Method enabled" : "Method disabled"}
          </span>
        </div>
        <TocsyPredictionFields options={options.tocsy} onChange={updateTocsyOption} disabled={running} />
      </div>
    </div>
  );

  // HSQC Tab
  const renderHSQCTab = () => (
    <div className="tab-content">
      <div className="settings-section">
        <h4 className="section-title">HSQC Parameters</h4>
        
        <div className="method-status-box">
          <span className={selectedPredMethods.includes("HSQC") ? "enabled" : "disabled"}>
            {selectedPredMethods.includes("HSQC") ? "Method enabled" : "Method disabled"}
          </span>
        </div>

        <div className="form-field">
          <label className="checkbox-label">
            <input
              type="checkbox"
              checked={hsqcFoldUI}
              onChange={(e) => updateExtractorOption('hsqc', 'fold', e.target.checked)}
              disabled={running}
            />
            <span>Folded spectrum</span>
          </label>
          <p className="hint">
            If checked, the F1 dimension (13C) will be folded between the specified bounds.
          </p>
        </div>

        {hsqcFoldUI && (
          <div className="params-grid">
            <div className="param-row">
              <label>F1 (13C) lower bound (ppm):</label>
              <input
                type="number"
                step="0.1"
                value={hsqcLowUI}
                onChange={(e) => updateExtractorOption('hsqc', 'fold_low', parseFloat(e.target.value))}
                disabled={running}
              />
            </div>
            <div className="param-row">
              <label>F1 (13C) upper bound (ppm):</label>
              <input
                type="number"
                step="0.1"
                value={hsqcHighUI}
                onChange={(e) => updateExtractorOption('hsqc', 'fold_high', parseFloat(e.target.value))}
                disabled={running}
              />
            </div>
          </div>
        )}
      </div>
    </div>
  );

  const tabs = [
    { id: 'source', label: 'Source & Methods' },
    { id: '1H', label: '¹H' },
    { id: '13C', label: '¹³C' },
    { id: 'COSY', label: 'COSY' },
    { id: 'TOCSY', label: 'TOCSY' },
    { id: 'HSQC', label: 'HSQC' },
  ].filter(tab => isMethodAvailable(tab.id));

  const elapsedMs = startTime ? (endTime || Date.now()) - startTime : 0;
  const canLaunch = methodCount > 0 && (mode === "SMILES" ? smiles.trim() : csvUploadedName);

  return (
    <div className="tabbed-prediction" style={{ position: "relative" }}>
      {dimmed && (
        <LockOverlay message="Annotation or optimization in progress..." />
      )}

      {/* MS-FINDER style tab bar */}
      <div className="tab-bar">
        {tabs.map((tab) => (
          <button
            key={tab.id}
            className={`tab-button ${activeTab === tab.id ? 'active' : ''}`}
            onClick={() => setActiveTab(tab.id)}
            disabled={running}
          >
            {tab.label}
          </button>
        ))}
      </div>

      {/* Active tab content */}
      <div className="tab-panel">
        {activeTab === 'source' && renderSourceTab()}
        {activeTab === '1H' && render1HTab()}
        {activeTab === '13C' && render13CTab()}
        {activeTab === 'COSY' && renderCOSYTab()}
        {isMethodAvailable('TOCSY') && activeTab === 'TOCSY' && renderTOCSYTab()}
        {activeTab === 'HSQC' && renderHSQCTab()}
      </div>

      {/* Progress */}
      {running && (
        <div className="progress-section">
          <ProgressBar progress={progress} message={progressMessage} />
          {totalMolecules > 0 && (
            <div className="progress-details">
              {completedMolecules} / {totalMolecules} molecules
            </div>
          )}
        </div>
      )}

      {/* Error */}
      {error && (
        <div className="error-box">
          Error: {error}
        </div>
      )}

      {/* Results */}
      {!running && artifacts?.length > 0 && (
        <div className="results-box">
          <span className="success-icon">OK</span>
          Prediction completed! {artifacts.length} file(s) generated.
        </div>
      )}

      {/* Action bar at bottom */}
      <div className="action-bar">
        {running ? (
          <>
            <div className="progress-inline">
              <span className="running-text">Prediction in progress...</span>
              {elapsedMs > 0 && <span className="elapsed">{formatDuration(elapsedMs)}</span>}
            </div>
            <button onClick={handleCancel} className="btn btn-cancel">
              Cancel
            </button>
          </>
        ) : (
          <>
            <button
              onClick={handleLaunchAndClose}
              disabled={!canLaunch}
              className="btn btn-finish"
            >
              {!canLaunch 
                ? "Configure first" 
                : `Launch prediction (${methodCount})`
              }
            </button>
            {inModal && onClose && (
              <button onClick={onClose} className="btn btn-secondary">
                Close
              </button>
            )}
          </>
        )}
      </div>
    </div>
  );
};

export default TabbedPredictionInterface;
