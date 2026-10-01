import { isMethodAvailable } from '../../constants/releaseFeatures';
import React, { useState, useEffect, useCallback } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import { LockOverlay, ProgressBar, FileUpload } from '../common';
import { formatFileSize, normalizeStaticUrl } from '../../utils';
import './TabbedAnnotationInterface.css';

/**
 * Tabbed annotation interface MS-FINDER style
 */
export const TabbedAnnotationInterface = ({ inModal = false, onClose }) => {
  const {
    annotationState,
    setAnnotationState,
    uiLock,
    setUiLock,
    addLog,
    // Use persistent annotation functions from context
    annotation,
  } = useAppState();

  const {
    matrixUploads,
    loadingMatrixUploads,
    cancelRequested,
    refreshMatrixUploads,
    loadPredictionTasks,
    uploadMatrix,
    hasCsvForMethod,
    canRunAnnotation,
    runSelectedIntegrated,
    handleCancelAll,
  } = annotation;

  const {
    selectedTask,
    allRunning,
    globalProgress,
    globalProgressMessage,
    readyToView,
    selectedMethods,
    uploadedMatrixNames,
    matrixFiles,
    uploadingMatrix,
    predictionTasks,
    parameters,
    bgPdataDirs,
    bgOptions,
    bgStatus,
  } = annotationState;

  // Active tab: 'source' | '1H' | '13C' | 'COSY' | 'TOCSY' | 'HSQC'
  const [activeTab, setActiveTab] = useState('source');

  const dimmed = uiLock.lockedBy === "PRED" || uiLock.lockedBy === "OPTIM";

  // UI Lock
  useEffect(() => {
    if (allRunning) {
      setUiLock({ lockedBy: "ANNO", message: "Annotations in progress..." });
    } else if (uiLock.lockedBy === "ANNO") {
      setUiLock({ lockedBy: null, message: "" });
    }
  }, [allRunning, setUiLock, uiLock.lockedBy]);

  // Load tasks on mount
  useEffect(() => {
    if (predictionTasks.length === 0) {
      loadPredictionTasks();
    }
    if (matrixUploads.length === 0) {
      refreshMatrixUploads();
    }
  }, [loadPredictionTasks, refreshMatrixUploads, predictionTasks.length, matrixUploads.length]);

  // Calculate number of ready methods
  const readySelected = selectedMethods.filter((t) => canRunAnnotation(t));
  const n = readySelected.length;

  const updateParameter = useCallback((type, key, value) => {
    const num = parseFloat(value);
    setAnnotationState((prev) => ({
      ...prev,
      parameters: {
        ...prev.parameters,
        [type]: {
          ...prev.parameters[type],
          [key]: !isNaN(num) ? num : value,
        },
      },
    }));
  }, [setAnnotationState]);

  // Source & Methods tab
  const renderSourceTab = () => {
    return (
      <div className="tab-content">
        {/* Task selection */}
        <div className="settings-section">
          <h4 className="section-title">Prediction Task</h4>
          <div className="form-row compact">
            <select
              className="form-select"
              value={selectedTask || ""}
              onChange={(e) =>
                setAnnotationState((prev) => ({
                  ...prev,
                  selectedTask: e.target.value,
                  readyToView: false,
                }))
              }
              disabled={allRunning}
            >
              <option value="">Select a task...</option>
              {predictionTasks.map((task) => (
                <option key={task.taskId} value={task.taskId}>
                  {task.taskId.slice(0, 8)}... - {new Date(task.createdAt).toLocaleString()}
                </option>
              ))}
            </select>
            <button
              className="btn-icon"
              onClick={loadPredictionTasks}
              disabled={allRunning}
              title="Refresh"
            >
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M23 4v6h-6M1 20v-6h6M3.51 9a9 9 0 0114.85-3.36L23 10M1 14l4.64 4.36A9 9 0 0020.49 15"/>
              </svg>
            </button>
          </div>
        </div>

        {/* Available CSVs */}
        {selectedTask && (
          <div className="settings-section">
            <h4 className="section-title">Available CSVs</h4>
            <div className="csv-grid">
              {["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).map((t) => {
                const valid = hasCsvForMethod(t);
                return (
                  <div key={t} className={`csv-item ${valid ? 'valid' : 'invalid'}`}>
                    <span className="csv-type">{t}</span>
                    <span className="csv-status">{valid ? '✓' : '✗'}</span>
                  </div>
                );
              })}
            </div>
          </div>
        )}

        {/* Methods to launch */}
        <div className="settings-section">
          <h4 className="section-title">Methods to Launch</h4>
          <p className="section-hint">Check the methods then configure each tab with the corresponding spectral matrix.</p>
          <div className="methods-list">
            {["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).map((t) => {
              const hasCsv = hasCsvForMethod(t);
              const hasMatrix = !!uploadedMatrixNames[t];
              const checked = selectedMethods.includes(t);
              
              return (
                <div key={t} className={`method-row ${!hasCsv ? 'disabled' : ''} ${checked ? 'selected' : ''}`}>
                  <label className="method-checkbox">
                    <input
                      type="checkbox"
                      checked={checked}
                      onChange={() => {
                        if (!hasCsv || allRunning) return;
                        setAnnotationState((prev) => {
                          const current = new Set(prev.selectedMethods);
                          if (current.has(t)) {
                            current.delete(t);
                          } else {
                            current.add(t);
                          }
                          return {
                            ...prev,
                            selectedMethods: Array.from(current),
                            activeType: t,
                            readyToView: false,
                          };
                        });
                      }}
                      disabled={!hasCsv || allRunning}
                    />
                    <span className="method-name">{t}</span>
                  </label>
                  <div className="method-statuses">
                    <span className={`status-badge ${hasCsv ? 'ok' : 'missing'}`}>
                      CSV {hasCsv ? 'OK' : 'X'}
                    </span>
                    <span className={`status-badge ${hasMatrix ? 'ok' : 'missing'}`}>
                      Matrix {hasMatrix ? 'OK' : 'X'}
                    </span>
                  </div>
                  <button
                    className="btn-configure"
                    onClick={() => setActiveTab(t)}
                    disabled={!hasCsv || allRunning}
                    title={`Configure ${t}`}
                  >
                    Configure
                  </button>
                </div>
              );
            })}
          </div>
          
          {/* Summary */}
          <div className="methods-summary">
            {n === 0 ? (
              <span className="summary-warning">Select methods and add missing matrices</span>
            ) : (
              <span className="summary-ok">{n} method{n > 1 ? 's' : ''} ready to launch</span>
            )}
          </div>
        </div>
      </div>
    );
  };

  // Tab for a specific type (1H, 13C, COSY, HSQC)
  const renderTypeTab = (type) => {
    const typeParams = parameters[type] || {};
    const opts = bgOptions[type] || {};
    const pdataDir = bgPdataDirs[type] || '';
    const status = bgStatus[type];

    const is1D = type === "1H" || type === "13C";

    return (
      <div className="tab-content">
        {/* Spectral matrix */}
        <div className="settings-section">
          <h4 className="section-title">Spectral Matrix {type}</h4>
          <div className="form-row compact">
            <select
              className="form-select"
              disabled={allRunning || loadingMatrixUploads}
              value={uploadedMatrixNames[type] || ""}
              onChange={(e) => {
                const name = e.target.value || "";
                setAnnotationState((prev) => ({
                  ...prev,
                  uploadedMatrixNames: {
                    ...prev.uploadedMatrixNames,
                    [type]: name,
                  },
                  readyToView: false,
                }));
                if (name) addLog(`[${type}] matrix selected: ${name}`);
              }}
            >
              <option value="">-- None --</option>
              {matrixUploads.map((u) => (
                <option key={u.filename} value={u.filename}>
                  {u.filename} ({formatFileSize(u.size)})
                </option>
              ))}
            </select>
            <button
              className="btn-icon"
              onClick={refreshMatrixUploads}
              disabled={allRunning}
              title="Refresh"
            >
              <svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
                <path d="M23 4v6h-6M1 20v-6h6M3.51 9a9 9 0 0114.85-3.36L23 10M1 14l4.64 4.36A9 9 0 0020.49 15"/>
              </svg>
            </button>
          </div>
          
          <div className="upload-zone">
            <FileUpload
              label={`Upload a CSV file for ${type}`}
              file={matrixFiles?.[type] || null}
              onChange={(e) => {
                const file = e.target.files?.[0];
                if (file) {
                  uploadMatrix(type, file);
                }
              }}
              accept=".csv"
              disabled={allRunning}
              uploading={!!uploadingMatrix?.[type]}
            />
          </div>
        </div>

        {/* Tolerance parameters */}
        <div className="settings-section">
          <h4 className="section-title">Tolerance Parameters</h4>
          <div className="params-grid">
            {is1D ? (
              <>
                <div className="param-row">
                  <label>Correspondence tolerance (ppm):</label>
                  <input
                    type="number"
                    step="0.01"
                    value={typeParams.tol_correspondence ?? (type === "1H" ? 0.1 : 2)}
                    onChange={(e) => updateParameter(type, 'tol_correspondence', e.target.value)}
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>Match tolerance (ppm):</label>
                  <input
                    type="number"
                    step="0.01"
                    value={typeParams.tol_match ?? (type === "1H" ? 0.1 : 2)}
                    onChange={(e) => updateParameter(type, 'tol_match', e.target.value)}
                    disabled={allRunning}
                  />
                </div>
              </>
            ) : (
              <>
                <div className="param-row">
                  <label>F1 corresp. tolerance:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={typeParams.tol_correspondance_F1 ?? (type === "COSY" ? 0.2 : 4)}
                    onChange={(e) => updateParameter(type, 'tol_correspondance_F1', e.target.value)}
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>F1 match tolerance:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={typeParams.tol_match_F1 ?? (type === "COSY" ? 0.2 : 4)}
                    onChange={(e) => updateParameter(type, 'tol_match_F1', e.target.value)}
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>F2 corresp. tolerance:</label>
                  <input
                    type="number"
                    step="0.01"
                    value={typeParams.tol_correspondance_F2 ?? 0.2}
                    onChange={(e) => updateParameter(type, 'tol_correspondance_F2', e.target.value)}
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>F2 match tolerance:</label>
                  <input
                    type="number"
                    step="0.01"
                    value={typeParams.tol_match_F2 ?? 0.2}
                    onChange={(e) => updateParameter(type, 'tol_match_F2', e.target.value)}
                    disabled={allRunning}
                  />
                </div>
              </>
            )}
          </div>
        </div>

        {/* Advanced parameters */}
        <div className="settings-section">
          <h4 className="section-title">Advanced Parameters</h4>
          <div className="params-grid">
            <div className="param-row">
              <label>Correlation threshold:</label>
              <input
                type="number"
                step="0.1"
                min="0"
                max="1"
                value={typeParams.correlation_threshold ?? (is1D ? 0.9 : 0.8)}
                onChange={(e) => updateParameter(type, 'correlation_threshold', e.target.value)}
                disabled={allRunning}
              />
            </div>
            <div className="param-row">
              <label>Correlation method:</label>
              <select
                value={typeParams.correlation_method ?? "pearson"}
                onChange={(e) => updateParameter(type, 'correlation_method', e.target.value)}
                disabled={allRunning}
              >
                <option value="pearson">Pearson</option>
                <option value="spearman">Spearman</option>
              </select>
            </div>
            <div className="param-row">
              <label>Ratio min:</label>
              <input
                type="number"
                step="0.01"
                min="0"
                value={typeParams.ratio_min ?? (type === "1H" || type === "COSY" ? 0.01 : 0.1)}
                onChange={(e) => updateParameter(type, 'ratio_min', e.target.value)}
                disabled={allRunning}
              />
            </div>
            <div className="param-row">
              <label>RANSY CV threshold:</label>
              <input
                type="number"
                step="0.01"
                min="0"
                value={typeParams.ransy_cv_threshold ?? (type === "1H" ? 0.05 : type === "13C" ? 0.1 : 0.15)}
                onChange={(e) => updateParameter(type, 'ransy_cv_threshold', e.target.value)}
                disabled={allRunning}
              />
            </div>
            <div className="param-row">
              <label>Max missing RANSY:</label>
              <input
                type="number"
                step="1"
                min="0"
                value={typeParams.max_missing_ransy ?? 50}
                onChange={(e) => updateParameter(type, 'max_missing_ransy', e.target.value)}
                disabled={allRunning}
              />
            </div>
            <div className="param-row">
              <label>Coverage threshold group:</label>
              <input
                type="number"
                step="0.1"
                min="0"
                value={typeParams.coverage_threshold_group ?? 0}
                onChange={(e) => updateParameter(type, 'coverage_threshold_group', e.target.value)}
                disabled={allRunning}
              />
            </div>
            <div className="param-row">
              <label>Coverage threshold subchain:</label>
              <input
                type="number"
                step="0.1"
                min="0"
                value={typeParams.coverage_threshold_subchain ?? 0}
                onChange={(e) => updateParameter(type, 'coverage_threshold_subchain', e.target.value)}
                disabled={allRunning}
              />
            </div>
            <div className="param-row">
              <label>Coverage threshold RANSY:</label>
              <input
                type="number"
                step="0.1"
                min="0"
                value={typeParams.coverage_threshold_ransy ?? 0}
                onChange={(e) => updateParameter(type, 'coverage_threshold_ransy', e.target.value)}
                disabled={allRunning}
              />
            </div>
            <div className="param-row">
              <label>Number of cores (nb_cores):</label>
              <input
                type="number"
                step="1"
                min="1"
                max="16"
                value={typeParams.nb_cores ?? 4}
                onChange={(e) => updateParameter(type, 'nb_cores', parseInt(e.target.value) || 4)}
                disabled={allRunning}
              />
            </div>
          </div>
        </div>

        {/* Spectral Normalization */}
        <div className="settings-section">
          <h4 className="section-title">Spectral Normalization</h4>
          <p className="section-hint" style={{ fontSize: '11px', color: '#666', marginBottom: '8px' }}>
            Corrects dilution effects between samples. PQN is recommended for metabolomics.
          </p>
          <div className="params-grid">
            <div className="param-row">
              <label>Normalization method:</label>
              <select
                value={
                  typeParams.use_pqn ? "pqn" :
                  typeParams.use_csn ? "csn" : "none"
                }
                onChange={(e) => {
                  const val = e.target.value;
                  updateParameter(type, 'use_pqn', val === 'pqn');
                  updateParameter(type, 'use_csn', val === 'csn');
                }}
                disabled={allRunning}
              >
                <option value="none">None (raw intensities)</option>
                <option value="csn">CSN (Constant Sum Normalization)</option>
                <option value="pqn">PQN (Probabilistic Quotient Normalization)</option>
              </select>
            </div>
          </div>
        </div>

        {/* Bruker backgrounds (optional) */}
        <div className="settings-section">
          <h4 className="section-title">Bruker Backgrounds (optional)</h4>
          <div className="form-row compact">
            <label>pdata folder:</label>
            <input
              type="text"
              className="form-input"
              placeholder="C:/data/exp/pdata/1"
              value={pdataDir}
              onChange={(e) =>
                setAnnotationState((prev) => ({
                  ...prev,
                  bgPdataDirs: {
                    ...prev.bgPdataDirs,
                    [type]: e.target.value,
                  },
                  // Clear the cached background status for this type when path changes
                  bgStatus: {
                    ...prev.bgStatus,
                    [type]: null,
                  },
                }))
              }
              disabled={allRunning}
            />
          </div>

          {/* Spectral window */}
          <div className="params-grid" style={{ marginTop: 8 }}>
            {is1D ? (
              <>
                <div className="param-row">
                  <label>ppm min:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={opts.ppm_min ?? 0}
                    onChange={(e) =>
                      setAnnotationState((prev) => ({
                        ...prev,
                        bgOptions: {
                          ...prev.bgOptions,
                          [type]: { ...opts, ppm_min: parseFloat(e.target.value) },
                        },
                      }))
                    }
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>ppm max:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={opts.ppm_max ?? 10}
                    onChange={(e) =>
                      setAnnotationState((prev) => ({
                        ...prev,
                        bgOptions: {
                          ...prev.bgOptions,
                          [type]: { ...opts, ppm_max: parseFloat(e.target.value) },
                        },
                      }))
                    }
                    disabled={allRunning}
                  />
                </div>
              </>
            ) : (
              <>
                <div className="param-row">
                  <label>F2 min:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={opts.f2_ppm_min ?? 0}
                    onChange={(e) =>
                      setAnnotationState((prev) => ({
                        ...prev,
                        bgOptions: {
                          ...prev.bgOptions,
                          [type]: { ...opts, f2_ppm_min: parseFloat(e.target.value) },
                        },
                      }))
                    }
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>F2 max:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={opts.f2_ppm_max ?? 10}
                    onChange={(e) =>
                      setAnnotationState((prev) => ({
                        ...prev,
                        bgOptions: {
                          ...prev.bgOptions,
                          [type]: { ...opts, f2_ppm_max: parseFloat(e.target.value) },
                        },
                      }))
                    }
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>F1 min:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={opts.f1_ppm_min ?? 0}
                    onChange={(e) =>
                      setAnnotationState((prev) => ({
                        ...prev,
                        bgOptions: {
                          ...prev.bgOptions,
                          [type]: { ...opts, f1_ppm_min: parseFloat(e.target.value) },
                        },
                      }))
                    }
                    disabled={allRunning}
                  />
                </div>
                <div className="param-row">
                  <label>F1 max:</label>
                  <input
                    type="number"
                    step="0.1"
                    value={opts.f1_ppm_max ?? 200}
                    onChange={(e) =>
                      setAnnotationState((prev) => ({
                        ...prev,
                        bgOptions: {
                          ...prev.bgOptions,
                          [type]: { ...opts, f1_ppm_max: parseFloat(e.target.value) },
                        },
                      }))
                    }
                    disabled={allRunning}
                  />
                </div>
              </>
            )}
          </div>

          {/* Background status */}
          {status && (
            <div className={`bg-status ${status.pngUrl ? 'success' : 'error'}`}>
              {status.pngUrl ? (
                <span>Background exported: <a href={normalizeStaticUrl(status.pngUrl)} target="_blank" rel="noreferrer">View PNG</a></span>
              ) : status.error ? (
                <span>Error: {status.error}</span>
              ) : null}
            </div>
          )}
        </div>
      </div>
    );
  };

  const tabs = [
    { id: 'source', label: 'Source & Methods' },
    { id: '1H', label: '¹H' },
    { id: '13C', label: '¹³C' },
    { id: 'COSY', label: 'COSY' },
    { id: 'TOCSY', label: 'TOCSY' },
    { id: 'HSQC', label: 'HSQC' },
  ].filter(tab => isMethodAvailable(tab.id));

  return (
    <div className="tabbed-annotation" style={{ position: "relative" }}>
      {dimmed && (
        <LockOverlay message="Prediction or optimization in progress..." />
      )}

      {/* MS-FINDER style tab bar */}
      <div className="tab-bar">
        {tabs.map((tab) => (
          <button
            key={tab.id}
            className={`tab-button ${activeTab === tab.id ? 'active' : ''}`}
            onClick={() => setActiveTab(tab.id)}
            disabled={allRunning}
          >
            {tab.label}
          </button>
        ))}
      </div>

      {/* Active tab content */}
      <div className="tab-panel">
        {activeTab === 'source' && renderSourceTab()}
        {activeTab === '1H' && renderTypeTab('1H')}
        {activeTab === '13C' && renderTypeTab('13C')}
        {activeTab === 'COSY' && renderTypeTab('COSY')}
        {isMethodAvailable('TOCSY') && activeTab === 'TOCSY' && renderTypeTab('TOCSY')}
        {activeTab === 'HSQC' && renderTypeTab('HSQC')}
      </div>

      {/* Action bar at bottom */}
      <div className="action-bar">
        {allRunning ? (
          <>
            <div className="progress-inline">
              <ProgressBar progress={globalProgress} message={globalProgressMessage} />
            </div>
            <button
              onClick={handleCancelAll}
              className="btn btn-cancel"
              disabled={cancelRequested}
            >
              {cancelRequested ? "Canceling..." : "Cancel"}
            </button>
          </>
        ) : (
          <>
            <button
              onClick={() => {
                runSelectedIntegrated();
                if (inModal && onClose) onClose();
              }}
              disabled={n === 0}
              className="btn btn-finish"
            >
              {n === 0 ? "Select methods" : `Launch (${n})`}
            </button>
            {inModal && onClose && (
              <button onClick={onClose} className="btn btn-cancel">
                Close
              </button>
            )}
          </>
        )}
      </div>

      {readyToView && !allRunning && (
        <div className="success-message">
          Annotation completed! Close this window to see results.
        </div>
      )}
    </div>
  );
};

export default TabbedAnnotationInterface;
