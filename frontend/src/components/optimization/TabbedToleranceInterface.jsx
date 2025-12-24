import React, { useState, useEffect, useCallback } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import { LockOverlay, ProgressBar } from '../common';
import { ToleranceRangeInputs } from './ToleranceRangeInputs';
import { ToleranceResultsCharts } from './ToleranceResultsCharts';
import { ToleranceSummaryTable } from './ToleranceSummaryTable';
import { 
  getPredictionTasks, 
  listUploads, 
  uploadFile,
  isValidCsvForType,
  API_BASE_URL
} from '../../api';
import { ProgressTracker, parseCSVSmart, normalizeStaticUrl, formatFileSize } from '../../utils';
import { createReconnectingSSE } from '../../createReconnectingSSE';
import './TabbedToleranceInterface.css';

/**
 * Tolerance optimization interface with MS-FINDER style tabs
 */
export const TabbedToleranceInterface = ({ onClose }) => {
  const {
    annotationState,
    setAnnotationState,
    addLog,
    setShowLogs,
    uiLock,
    setUiLock,
  } = useAppState();

  const { toleranceRanges, toleranceResults, parameters } = annotationState;

  // Onglet actif: 'source' ou '1H', '13C', 'COSY', 'HSQC', 'results'
  const [activeTab, setActiveTab] = useState('source');
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState(0);
  const [currentTest, setCurrentTest] = useState("");
  const [matrixUploads, setMatrixUploads] = useState([]);
  const [loadingMatrixUploads, setLoadingMatrixUploads] = useState(false);
  // Type sélectionné pour l'optimisation
  const [selectedType, setSelectedType] = useState('1H');

  const dimmed = uiLock.lockedBy === "PRED" || uiLock.lockedBy === "ANNO";

  // Clés de tolérance à masquer des paramètres
  const toleranceKeysByType = {
    "1H": ["tol_correspondence", "tol_correspondance", "tol_match"],
    "13C": ["tol_correspondence", "tol_correspondance", "tol_match"],
    COSY: [
      "tol_correspondence_f1", "tol_correspondance_f1", "tol_match_f1",
      "tol_correspondence_f2", "tol_correspondance_f2", "tol_match_f2",
    ],
    HSQC: [
      "tol_correspondence_f1", "tol_correspondance_f1", "tol_match_f1",
      "tol_correspondence_f2", "tol_correspondance_f2", "tol_match_f2",
    ],
  };

  useEffect(() => {
    if (running) {
      setUiLock({ lockedBy: "OPTIM", message: "Tolerance optimization in progress..." });
    } else if (uiLock.lockedBy === "OPTIM") {
      setUiLock({ lockedBy: null, message: "" });
    }
  }, [running, setUiLock, uiLock.lockedBy]);

  const refreshMatrixUploads = useCallback(async () => {
    setLoadingMatrixUploads(true);
    try {
      const list = await listUploads("matrix");
      setMatrixUploads(list);
    } catch (e) {
      addLog(`Error listUploads(matrix): ${e.message}`, "error");
    } finally {
      setLoadingMatrixUploads(false);
    }
  }, [addLog]);

  const loadPredictionTasks = useCallback(async () => {
    try {
      const tasks = await getPredictionTasks();
      setAnnotationState((prev) => ({ ...prev, predictionTasks: tasks }));
      if (tasks.length > 0 && !annotationState.selectedTask) {
        setAnnotationState((prev) => ({ ...prev, selectedTask: tasks[0].taskId }));
      }
    } catch (e) {
      addLog(`Error loading prediction tasks: ${e.message}`, "error");
    }
  }, [addLog, annotationState.selectedTask, setAnnotationState]);

  useEffect(() => {
    if (annotationState.predictionTasks.length === 0) {
      loadPredictionTasks();
    }
    if (matrixUploads.length === 0) {
      refreshMatrixUploads();
    }
  }, [loadPredictionTasks, refreshMatrixUploads, annotationState.predictionTasks.length, matrixUploads.length]);

  const uploadMatrix = async (type, file) => {
    if (!file) return;
    setAnnotationState((prev) => ({
      ...prev,
      uploadingMatrix: { ...prev.uploadingMatrix, [type]: true },
    }));
    try {
      const { filename } = await uploadFile(file);
      setAnnotationState((prev) => ({
        ...prev,
        uploadedMatrixNames: { ...prev.uploadedMatrixNames, [type]: filename },
        matrixFiles: { ...prev.matrixFiles, [type]: file },
        uploadingMatrix: { ...prev.uploadingMatrix, [type]: false },
      }));
      addLog(`[${type}] Matrix uploaded: ${filename}`);
    } catch (e) {
      addLog(`[${type}] Matrix upload failed: ${e.message}`, "error");
      setAnnotationState((prev) => ({
        ...prev,
        uploadingMatrix: { ...prev.uploadingMatrix, [type]: false },
      }));
    }
  };

  const updateParameter = useCallback(
    (type, key, value) => {
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
    },
    [setAnnotationState]
  );

  const selectedTaskData = annotationState.predictionTasks.find(
    (t) => t.taskId === annotationState.selectedTask
  );

  const hasCsvForMethod = useCallback(
    (type) => {
      const csv = selectedTaskData?.csvFiles?.[type];
      return !!csv && isValidCsvForType(csv, type) && !!annotationState.selectedTask;
    },
    [annotationState.selectedTask, selectedTaskData]
  );

  const canRunTest = () => {
    return hasCsvForMethod(selectedType) && !!annotationState.uploadedMatrixNames[selectedType];
  };

  // Analyser résultat CSV annotation
  const analyzeAnnotationResult = async (csvUrl, toleranceParams, totalExperimentalSignals) => {
    try {
      const normalized = normalizeStaticUrl(csvUrl);
      const response = await fetch(normalized);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);

      const csvText = await response.text();
      const { data } = parseCSVSmart(csvText);

      if (!data || data.length === 0) {
        return {
          toleranceParams,
          coverageDistribution: Array(5).fill(0),
          totalSignalsMatched: 0,
          totalExperimentalSignals: totalExperimentalSignals ?? null,
          signalsMatchedPct: 0,
          totalMolecules: 0,
          avgCoverage: 0,
        };
      }

      const coverageDistribution = Array(5).fill(0);
      let totalCoverageSum = 0;
      const matchedExpSignals = new Set();

      data.forEach((row) => {
        let coverage = parseFloat(row.coverage_ransy ?? row.coverage_sub ?? 0);
        if (!Number.isFinite(coverage)) coverage = 0;
        coverage = Math.max(0, Math.min(100, coverage));

        const bin = Math.min(4, Math.floor(coverage / 20));
        coverageDistribution[bin]++;
        totalCoverageSum += coverage;

        if (row.RANSY_signals && String(row.RANSY_signals).trim() !== "") {
          const expSignals = String(row.RANSY_signals)
            .split("|")
            .map((s) => s.trim())
            .filter((s) => s !== "");

          expSignals.forEach((sig) => {
            matchedExpSignals.add(sig);
          });
        }
      });

      const totalMolecules = data.length;
      const avgCoverage = totalCoverageSum / totalMolecules;
      const uniqueMatchedCount = matchedExpSignals.size;
      const signalsMatchedPct =
        totalExperimentalSignals && totalExperimentalSignals > 0
          ? (uniqueMatchedCount / totalExperimentalSignals) * 100
          : 0;

      return {
        toleranceParams,
        coverageDistribution,
        totalSignalsMatched: uniqueMatchedCount,
        totalExperimentalSignals: totalExperimentalSignals ?? null,
        signalsMatchedPct,
        totalMolecules,
        avgCoverage,
        matchedSignals: Array.from(matchedExpSignals).map(s => parseFloat(s)).sort((a, b) => a - b),
      };
    } catch (e) {
      console.error("Error analyzing CSV:", e);
      return {
        toleranceParams,
        coverageDistribution: Array(5).fill(0),
        totalSignalsMatched: 0,
        totalExperimentalSignals: totalExperimentalSignals ?? null,
        signalsMatchedPct: 0,
        totalMolecules: 0,
        avgCoverage: 0,
        matchedSignals: [],
        error: e.message,
      };
    }
  };

  const runSingleAnnotation = (predTaskId, matrixFile, annotationType, params) => {
    return new Promise((resolve, reject) => {
      createReconnectingSSE({
        endpoint: "/run-annotation-from-prediction",
        body: {
          predictionTaskId: predTaskId,
          matrixFile: matrixFile,
          annotationType: annotationType,
          parameters: params,
        },
        onProgress: () => {},
        onLog: () => {},
        onError: (err) => reject(err),
        onComplete: (res) => resolve(res),
      });
    });
  };

  const runToleranceTests = async () => {
    const type = selectedType;
    const ranges = toleranceRanges[type];
    const baseParams = parameters[type];

    const computeTotalExperimentalSignalsForMatrix = async (filename) => {
      if (!filename || typeof filename !== "string") return null;
      const cleanName = filename.trim();
      if (!cleanName) return null;

      try {
        const url = `${API_BASE_URL}/uploads/${encodeURIComponent(cleanName)}`;
        const resp = await fetch(url);
        if (!resp.ok) {
          throw new Error(`HTTP ${resp.status} (${url})`);
        }

        const rawText = await resp.text();
        const lines = rawText
          .replace(/^\uFEFF/, "")
          .split(/\r?\n/)
          .map((l) => l.trim())
          .filter((l) => l.length > 0);

        if (lines.length === 0) return 0;
        const totalSignals = lines.length - 1;
        addLog(`Matrix "${cleanName}": ${lines.length} lines, ${totalSignals} signals detected`, "info");
        return totalSignals > 0 ? totalSignals : 0;
      } catch (e) {
        addLog(`Error reading experimental matrix (${filename}): ${e.message}`, "error");
        return null;
      }
    };

    const matrixName = annotationState.uploadedMatrixNames?.[type];
    const totalExperimentalSignals = matrixName
      ? await computeTotalExperimentalSignalsForMatrix(matrixName)
      : null;

    if (totalExperimentalSignals == null) {
      addLog(`Unable to compute total experimental signals for matrix "${matrixName}".`, "error");
    } else if (totalExperimentalSignals === 0) {
      addLog(`Warning: Matrix "${matrixName}" is empty or invalid (0 signals detected).`, "warning");
    } else {
      addLog(`Matrix "${matrixName}" → ${totalExperimentalSignals} experimental signals detected.`, "info");
    }

    // Générer grille de tolérances
    const testCombinations = [];

    if (type === "1H" || type === "13C") {
      for (
        let tol = ranges.tol_correspondence.min;
        tol <= ranges.tol_correspondence.max + 1e-9;
        tol += ranges.tol_correspondence.step
      ) {
        const roundedTol = parseFloat(tol.toFixed(3));
        testCombinations.push({
          tol_correspondence: roundedTol,
          tol_match: roundedTol,
        });
      }
    } else if (type === "COSY" || type === "HSQC") {
      for (
        let tolF1 = ranges.tol_correspondence_f1.min;
        tolF1 <= ranges.tol_correspondence_f1.max + 1e-9;
        tolF1 += ranges.tol_correspondence_f1.step
      ) {
        for (
          let tolF2 = ranges.tol_correspondence_f2.min;
          tolF2 <= ranges.tol_correspondence_f2.max + 1e-9;
          tolF2 += ranges.tol_correspondence_f2.step
        ) {
          const roundedF1 = parseFloat(tolF1.toFixed(3));
          const roundedF2 = parseFloat(tolF2.toFixed(3));

          testCombinations.push({
            tol_correspondence_f1: roundedF1,
            tol_correspondence_f2: roundedF2,
            tol_match_f1: roundedF1,
            tol_match_f2: roundedF2,
          });
        }
      }
    }

    addLog(`Launching ${testCombinations.length} tests for ${type}...`);
    setRunning(true);
    setShowLogs(true);

    const tracker = new ProgressTracker();
    const phases = testCombinations.map((_, i) => ({
      name: `Test ${i + 1}`,
      weight: 100 / testCombinations.length,
    }));
    tracker.definePhases(phases);

    const results = [];

    for (let i = 0; i < testCombinations.length; i++) {
      const tolParams = testCombinations[i];

      const fullParams = {
        ...baseParams,
        ...tolParams,
      };

      let label;
      if (type === "1H" || type === "13C") {
        label = `tol=${tolParams.tol_correspondence} ppm`;
      } else {
        label = `F1=${tolParams.tol_correspondence_f1} / F2=${tolParams.tol_correspondence_f2}`;
      }

      setCurrentTest(`Test ${i + 1}/${testCombinations.length} : ${label}`);
      addLog(`[${i + 1}/${testCombinations.length}] ${label}...`);

      const updateTestProgress = (subPhase) => {
        const phaseProgress = (subPhase / 3) * 100;
        const globalPercent = tracker.updatePhase(i, phaseProgress);
        setProgress(globalPercent);
      };

      updateTestProgress(0);

      try {
        updateTestProgress(1);

        const result = await runSingleAnnotation(
          annotationState.selectedTask,
          annotationState.uploadedMatrixNames[type],
          type,
          fullParams
        );

        updateTestProgress(2);

        if (result.success && result.downloadUrl) {
          const metrics = await analyzeAnnotationResult(
            result.downloadUrl,
            tolParams,
            totalExperimentalSignals
          );

          results.push(metrics);

          addLog(
            `OK ${label} : ${metrics.totalSignalsMatched} signals matched` +
              (metrics.totalExperimentalSignals
                ? ` (${metrics.signalsMatchedPct.toFixed(2)}% of experimental matrix)`
                : "") +
              `, avg coverage = ${metrics.avgCoverage.toFixed(2)}%`
          );
        } else {
          addLog(`Failed ${label} : ${result.error}`, "error");
          results.push({
            toleranceParams: tolParams,
            error: result.error,
            totalSignalsMatched: 0,
            totalExperimentalSignals: totalExperimentalSignals ?? null,
            signalsMatchedPct: 0,
            totalMolecules: 0,
            avgCoverage: 0,
            matchedSignals: [],
          });
        }

        updateTestProgress(3);
      } catch (e) {
        addLog(`Error ${label} : ${e.message}`, "error");
        results.push({
          toleranceParams: tolParams,
          error: e.message,
          totalSignalsMatched: 0,
          totalExperimentalSignals: totalExperimentalSignals ?? null,
          signalsMatchedPct: 0,
          totalMolecules: 0,
          avgCoverage: 0,
          matchedSignals: [],
        });
      }
    }

    const resultSet = {
      id: Date.now(),
      type,
      matrix: matrixName,
      timestamp: new Date().toISOString(),
      label: `${type} @ ${new Date().toLocaleTimeString()} - ${matrixName}`,
      results,
      totalExperimentalSignals,
    };

    setAnnotationState((prev) => ({
      ...prev,
      toleranceResults: [...(prev.toleranceResults || []), resultSet],
    }));

    setRunning(false);
    setProgress(0);
    setCurrentTest("");
    addLog(`Tests completed for ${type}. ${results.length} results recorded.`);
    
    // Switch to results tab
    setActiveTab('results');
  };

  // Render Source & Method tab
  const renderSourceTab = () => (
    <div className="tab-content">
      {/* Prediction task */}
      <div className="param-section">
        <div className="section-title">Source Prediction Task</div>
        <div className="param-row">
          <label>Prediction task</label>
          <div className="input-with-button">
            <select
              value={annotationState.selectedTask || ""}
              onChange={(e) =>
                setAnnotationState((prev) => ({
                  ...prev,
                  selectedTask: e.target.value,
                }))
              }
              disabled={running}
            >
              <option value="">Select a task...</option>
              {annotationState.predictionTasks.map((task) => (
                <option key={task.taskId} value={task.taskId}>
                  {task.taskId.slice(0, 8)}... - {new Date(task.createdAt).toLocaleString()}
                  {task.canRunAnnotations ? " ✓" : ""}
                </option>
              ))}
            </select>
            <button
              className="btn-icon"
              onClick={loadPredictionTasks}
              disabled={running}
              title="Refresh"
            >
              ↻
            </button>
          </div>
        </div>

        {annotationState.selectedTask && (
          <div className="csv-status-grid">
            {["1H", "13C", "COSY", "HSQC"].map((t) => {
              const csv = selectedTaskData?.csvFiles?.[t];
              const valid = csv && isValidCsvForType(csv, t);
              return (
                <div key={t} className={`csv-status-item ${valid ? 'valid' : 'invalid'}`}>
                  <span className="status-icon">{valid ? '✓' : '✗'}</span>
                  <span className="status-label">{t}</span>
                </div>
              );
            })}
          </div>
        )}
      </div>

      {/* Method to optimize */}
      <div className="param-section">
        <div className="section-title">Method to Optimize</div>
        <div className="method-grid">
          {["1H", "13C", "COSY", "HSQC"].map((t) => {
            const available = hasCsvForMethod(t);
            return (
              <label key={t} className={`method-option ${!available ? 'disabled' : ''}`}>
                <input
                  type="radio"
                  name="optimizeType"
                  value={t}
                  checked={selectedType === t}
                  onChange={() => setSelectedType(t)}
                  disabled={running || !available}
                />
                <span className="method-name">{t}</span>
                <span className={`method-status ${available ? 'available' : 'missing'}`}>
                  {available ? 'Available' : 'Missing'}
                </span>
              </label>
            );
          })}
        </div>
      </div>
    </div>
  );

  // Render configuration tab for a type
  const renderTypeTab = (type) => {
    const typeParams = parameters[type] || {};
    const typeHiddenKeys = toleranceKeysByType[type] || [];

    return (
      <div className="tab-content">
        {/* Spectral matrix */}
        <div className="param-section">
          <div className="section-title">Spectral Matrix {type}</div>
          <div className="param-row">
            <label>Select a matrix</label>
            <div className="input-with-button">
              <select
                disabled={running || loadingMatrixUploads}
                value={annotationState.uploadedMatrixNames[type] || ""}
                onChange={(e) => {
                  const name = e.target.value || "";
                  setAnnotationState((prev) => ({
                    ...prev,
                    uploadedMatrixNames: {
                      ...prev.uploadedMatrixNames,
                      [type]: name,
                    },
                  }));
                  if (name) addLog(`[${type}] matrix selected: ${name}`);
                }}
              >
                <option value="">-- None --</option>
                {matrixUploads.map((f) => (
                  <option key={f.filename} value={f.filename}>
                    {f.filename} {Number.isFinite(f.size) ? `(${formatFileSize(f.size)})` : ""}
                  </option>
                ))}
              </select>
              <button
                className="btn-icon"
                onClick={refreshMatrixUploads}
                disabled={running || loadingMatrixUploads}
                title="Refresh"
              >
                ↻
              </button>
            </div>
          </div>

          <div className="upload-zone">
            <span>Or drag a CSV here</span>
            <input
              type="file"
              accept=".csv"
              onChange={(e) => uploadMatrix(type, e.target.files?.[0])}
              disabled={running}
            />
          </div>
        </div>

        {/* Tolerance ranges */}
        <div className="param-section">
          <div className="section-title">Tolerance Ranges to Test</div>
          <ToleranceRangeInputs
            type={type}
            ranges={toleranceRanges[type]}
            setAnnotationState={setAnnotationState}
            running={running}
          />
        </div>

        {/* Other parameters */}
        <div className="param-section">
          <div className="section-title">Other Parameters</div>
          <div className="params-grid">
            {Object.entries(typeParams)
              .filter(([key]) => !typeHiddenKeys.includes(key))
              .map(([key, value]) => (
                <div key={key} className="param-item">
                  <label>{key.replace(/_/g, " ")}</label>
                  <input
                    type="number"
                    step="any"
                    value={value}
                    onChange={(e) => updateParameter(type, key, e.target.value)}
                    disabled={running}
                  />
                </div>
              ))}
          </div>
        </div>
      </div>
    );
  };

  // Rendu de l'onglet Résultats
  const renderResultsTab = () => (
    <div className="tab-content">
      {toleranceResults && toleranceResults.length > 0 ? (
        <div className="results-container">
          {toleranceResults.map((resultSet) => (
            <div key={resultSet.id} className="result-set">
              <ToleranceSummaryTable resultSet={resultSet} />
              <div className="charts-container">
                <ToleranceResultsCharts
                  resultSet={resultSet}
                  onDelete={(id) => {
                    const result = toleranceResults.find((r) => r.id === id);
                    const label = result?.label || `Result ${id}`;
                    if (window.confirm(`Delete "${label}"?`)) {
                      setAnnotationState((prev) => ({
                        ...prev,
                        toleranceResults: prev.toleranceResults.filter((r) => r.id !== id),
                      }));
                      addLog(`Result deleted: ${label}`, "info");
                    }
                  }}
                />
              </div>
            </div>
          ))}
        </div>
      ) : (
        <div className="empty-results">
          <div className="empty-icon">
            <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
              <path d="M3 3v18h18" />
              <path d="M18 17V9" />
              <path d="M13 17V5" />
              <path d="M8 17v-3" />
            </svg>
          </div>
          <p>No optimization results</p>
          <p className="hint">Run tolerance tests to see results here</p>
        </div>
      )}
    </div>
  );

  return (
    <div className="tabbed-tolerance-interface">
      {dimmed && (
        <LockOverlay message="Prediction or annotation in progress..." />
      )}

      {/* Tab bar */}
      <div className="tabs-bar">
        <button
          className={`tab-btn ${activeTab === 'source' ? 'active' : ''}`}
          onClick={() => setActiveTab('source')}
          disabled={running}
        >
          Source & Method
        </button>
        <button
          className={`tab-btn ${activeTab === '1H' ? 'active' : ''}`}
          onClick={() => setActiveTab('1H')}
          disabled={running}
        >
          ¹H
        </button>
        <button
          className={`tab-btn ${activeTab === '13C' ? 'active' : ''}`}
          onClick={() => setActiveTab('13C')}
          disabled={running}
        >
          ¹³C
        </button>
        <button
          className={`tab-btn ${activeTab === 'COSY' ? 'active' : ''}`}
          onClick={() => setActiveTab('COSY')}
          disabled={running}
        >
          COSY
        </button>
        <button
          className={`tab-btn ${activeTab === 'HSQC' ? 'active' : ''}`}
          onClick={() => setActiveTab('HSQC')}
          disabled={running}
        >
          HSQC
        </button>
        <button
          className={`tab-btn ${activeTab === 'results' ? 'active' : ''}`}
          onClick={() => setActiveTab('results')}
          disabled={running}
        >
          Results {toleranceResults?.length > 0 && `(${toleranceResults.length})`}
        </button>
      </div>

      {/* Active tab content */}
      <div className="tab-panel">
        {activeTab === 'source' && renderSourceTab()}
        {activeTab === '1H' && renderTypeTab('1H')}
        {activeTab === '13C' && renderTypeTab('13C')}
        {activeTab === 'COSY' && renderTypeTab('COSY')}
        {activeTab === 'HSQC' && renderTypeTab('HSQC')}
        {activeTab === 'results' && renderResultsTab()}
      </div>

      {/* Action bar */}
      <div className="action-bar">
        {running && (
          <div className="progress-info">
            <ProgressBar progress={progress} message={currentTest} />
          </div>
        )}
        <div className="action-buttons">
          <button
            className="btn-primary"
            onClick={runToleranceTests}
            disabled={running || !canRunTest()}
          >
            {running ? `Tests in progress... (${Math.round(progress)}%)` : `Run tests ${selectedType}`}
          </button>
          <button
            className="btn-secondary"
            onClick={onClose}
            disabled={running}
          >
            Close
          </button>
        </div>
      </div>
    </div>
  );
};

export default TabbedToleranceInterface;
