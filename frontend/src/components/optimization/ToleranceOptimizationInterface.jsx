import { isMethodAvailable } from '../../constants/releaseFeatures';
import React, { useState, useEffect, useCallback } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import { LockOverlay, ProgressBar, ParameterInput } from '../common';
import { ToleranceRangeInputs } from './ToleranceRangeInputs';
import { ToleranceResultsCharts } from './ToleranceResultsCharts';
import { ToleranceSummaryTable } from './ToleranceSummaryTable';
import { 
  getPredictionTasks, 
  listUploads, 
  isValidCsvForType
} from '../../api';
import { formatFileSize, shortCsvLabel } from '../../utils';

/**
 * Main interface for tolerance optimization
 * @param {boolean} inModal - If true, adapts the style for a modal
 * @param {function} onClose - Callback to close the modal
 */
export const ToleranceOptimizationInterface = ({ inModal = false, onClose }) => {
  const {
    annotationState,
    setAnnotationState,
    addLog,
    uiLock,
    currentProject,
    toleranceOptimization,
  } = useAppState();

  const { toleranceRanges, toleranceResults, activeType, parameters } = annotationState;

  const {
    running = false,
    progress = 0,
    currentTest = "",
    runToleranceTests,
    canRunTest,
  } = toleranceOptimization || {};
  const [matrixUploads, setMatrixUploads] = useState([]);
  const [loadingMatrixUploads, setLoadingMatrixUploads] = useState(false);

  const dimmed = uiLock.lockedBy === "PRED" || uiLock.lockedBy === "ANNO";

  // Tolerance fields hidden from the parameter controls.
  const toleranceKeysByType = {
    "1H": ["tol_correspondence", "tol_correspondance", "tol_match"],
    "13C": ["tol_correspondence", "tol_correspondance", "tol_match"],
    COSY: [
      "tol_correspondence_f1",
      "tol_correspondance_f1",
      "tol_match_f1",
      "tol_correspondence_f2",
      "tol_correspondance_f2",
      "tol_match_f2",
    ],
    HSQC: [
      "tol_correspondence_f1",
      "tol_correspondance_f1",
      "tol_match_f1",
      "tol_correspondence_f2",
      "tol_correspondance_f2",
      "tol_match_f2",
    ],
  };

  const hiddenKeys = toleranceKeysByType[activeType] || [];
  const activeParams = parameters[activeType] || {};

  // Lock UI is handled by the global hook (toleranceOptimization)

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
  }, [
    loadPredictionTasks,
    refreshMatrixUploads,
    annotationState.predictionTasks.length,
    matrixUploads.length,
  ]);

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

  const renderCsvStatus = () => {
    if (!annotationState.selectedTask) return null;
    const csvMap = selectedTaskData?.csvFiles || {};
    const order = ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable);
    return (
      <div className="csv-status">
        {order.map((t) => {
          const f = csvMap[t];
          if (!f)
            return (
              <div key={t} className="status-item missing">
                {t}: Missing
              </div>
            );
          if (!isValidCsvForType(f, t))
            return (
              <div key={t} className="status-item missing">
                {t}: Invalid (should contain {"prediction_" + t + "_"}
              </div>
            );
          return (
            <div key={t} className="status-item success">
              {t}: {shortCsvLabel(f)}
            </div>
          );
        })}
      </div>
    );
  };

  // Check that a project is open (after hooks to respect React rules)
  if (!currentProject) {
    return (
      <div className="tolerance-interface">
        <div className="warning-box">
          <span className="warning-icon" aria-hidden="true">!</span>
          <div>
            <h3 className="warning-title">No project open</h3>
            <p className="warning-text">
              You must create or open a project from the home page before using tolerance optimization.
              Files will be stored in your project folder.
            </p>
            <button 
              onClick={() => window.location.href = '/'} 
              className="btn btn-primary"
            >
              Return to home
            </button>
          </div>
        </div>
      </div>
    );
  }

  return (
    <div className="annotation-interface" style={{ position: "relative" }}>
      {dimmed && (
        <LockOverlay message="Prediction or annotation in progress... Optimization is temporarily unavailable." />
      )}

      <div className="header">
        <h1>PPM Tolerance Optimization</h1>
        <p>Test multiple tolerance values to identify optimal parameters</p>
      </div>

      <div className="info-box">
        <strong>Note:</strong> This interface runs annotations with different
        tolerance values to generate optimization curves. The{" "}
        <code>tol_match</code> parameters will automatically follow the{" "}
        <code>tol_correspondence</code> values.
      </div>

      {/* Source prediction task */}
      <div className="card">
        <h3>Source Prediction Task</h3>
        <div className="form-row">
          <label className="form-label">Prediction task</label>
          <select
            className="form-input"
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
                {task.taskId} - {new Date(task.createdAt).toLocaleString()}
                {task.canRunAnnotations ? " ✓" : " !"}
              </option>
            ))}
          </select>
          <button
            className="btn btn-secondary"
            onClick={loadPredictionTasks}
            disabled={running}
          >
            Refresh
          </button>
        </div>

        {annotationState.selectedTask && (
          <div className="task-details">
            <h4>CSV available in {annotationState.selectedTask}:</h4>
            {renderCsvStatus()}
          </div>
        )}
      </div>

      {/* Method selection */}
      <div className="type-selector">
        <h3>Method to Optimize</h3>
        <div className="type-tabs">
          {["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).map((t) => {
            const available = hasCsvForMethod(t);
            return (
              <button
                key={t}
                onClick={() =>
                  running
                    ? null
                    : setAnnotationState((prev) => ({
                        ...prev,
                        activeType: t,
                      }))
                }
                className={`type-tab ${activeType === t ? "active" : ""} ${
                  !available ? "disabled" : ""
                }`}
                disabled={running || !available}
              >
                <div className="tab-title">{t}</div>
                <div className="tab-description">
                  {available ? "✓ Available" : "✗ Invalid / missing"}
                </div>
              </button>
            );
          })}
        </div>
      </div>

      {/* Matrice spectrale */}
      <div className="card">
        <h3>Matrice spectrale</h3>

        <div className="form-row">
          <label className="form-label">Sélectionner une matrice existante (uploads)</label>
          <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
            <select
              className="form-input"
              disabled={running || loadingMatrixUploads}
              value={annotationState.uploadedMatrixNames[activeType] || ""}
              onChange={(e) => {
                const name = e.target.value || "";
                setAnnotationState((prev) => ({
                  ...prev,
                  uploadedMatrixNames: {
                    ...prev.uploadedMatrixNames,
                    [activeType]: name,
                  },
                  matrixFiles: {
                    ...prev.matrixFiles,
                    [activeType]: null,
                  },
                }));
                if (name) addLog(`[${activeType}] matrice sélectionnée: ${name}`);
              }}
            >
              <option value="">— Choisir dans uploads —</option>
              {matrixUploads.map((f) => (
                <option key={f.filename} value={f.filename}>
                  {f.filename}{" "}
                  {Number.isFinite(f.size) ? `(${formatFileSize(f.size)})` : ""}{" "}
                  {f.mtime ? `- ${new Date(f.mtime).toLocaleString()}` : ""}
                </option>
              ))}
            </select>
            <button
              className="btn btn-secondary"
              onClick={refreshMatrixUploads}
              disabled={running || loadingMatrixUploads}
            >
              {loadingMatrixUploads ? "Loading..." : "Refresh"}
            </button>
          </div>
        </div>
      </div>

      {/* Annotation parameters (except tolerances) */}
      <div className="card">
        <h3>Parameters for {activeType}</h3>
        <div className="parameters-grid">
          {Object.entries(activeParams)
            .filter(([key]) => !hiddenKeys.includes(key))
            .map(([key, value]) => (
              <ParameterInput
                key={key}
                label={key.replace(/_/g, " ")}
                value={value}
                onChange={(e) => updateParameter(activeType, key, e.target.value)}
                disabled={running}
              />
            ))}
        </div>
      </div>

      {/* Tolerance windows */}
      <ToleranceRangeInputs
        type={activeType}
        ranges={toleranceRanges[activeType]}
        setAnnotationState={setAnnotationState}
        running={running}
      />

      {/* Launch tests button */}
      <div className="action-section">
        <button
          onClick={() => {
            if (typeof runToleranceTests === "function") runToleranceTests();
            if (inModal && onClose) onClose();
          }}
          disabled={running || !(typeof canRunTest === "function" ? canRunTest() : false)}
          className="btn btn-primary btn-large"
        >
          {running ? `Tests in progress... (${progress}%)` : `Run tests for ${activeType}`}
        </button>
        {running && (
          <div className="progress-section" style={{ marginTop: 12 }}>
            <ProgressBar progress={progress} message={currentTest} />
          </div>
        )}
      </div>

      {/* Test results */}
      {toleranceResults && toleranceResults.length > 0 && (
        <div className="all-tolerance-results">
          {toleranceResults.filter(result => isMethodAvailable(result.type)).map((resultSet) => (
            <div key={resultSet.id} className="tolerance-result-section">
              {/* Summary table */}
              <ToleranceSummaryTable resultSet={resultSet} />
              
              {/* Results charts */}
              <div style={{ marginTop: 30 }}>
                <ToleranceResultsCharts
                  resultSet={resultSet}
                  onDelete={(id) => {
                    const result = annotationState.toleranceResults.find((r) => r.id === id);
                    const label = result?.label || `Result ${id}`;
                    const confirmed = window.confirm(
                      `Are you sure you want to delete:\n"${label}"?\n\nThe table and charts will be removed.`
                    );
                    
                    if (confirmed) {
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
      )}
    </div>
  );
};