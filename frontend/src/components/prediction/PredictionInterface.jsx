import React, { useEffect } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import { LockOverlay } from '../common';
import { PredictionForm } from './PredictionForm';
import { PredictionProgress } from './PredictionProgress';
import { PredictionResults } from './PredictionResults';
import { formatDuration } from '../../utils';

/**
 * Main interface for NMR predictions (SMILES/CSV)
 * @param {boolean} inModal - If true, adapts style for modal
 * @param {function} onClose - Callback to close the modal
 */
export const PredictionInterface = ({ inModal = false, onClose }) => {
  const {
    predictionState,
    setPredictionState,
    uiLock,
    setUiLock,
    currentProject,
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

  const dimmed = uiLock.lockedBy === "ANNO" || uiLock.lockedBy === "OPTIM";

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
    if (predictionState.mode === "CSV" && csvUploads.length === 0) {
      refreshCsvUploads();
    }
  }, [predictionState.mode, csvUploads.length, refreshCsvUploads]);

  // Check that a project is open (after hooks to respect React rules)
  if (!currentProject) {
    return (
      <div className="prediction-interface">
        <div className="warning-box">
          <span className="warning-icon">Warning</span>
          <div>
            <h3 className="warning-title">No project open</h3>
            <p className="warning-text">
              You must create or open a project from the home page before using predictions.
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

  const elapsedMs = startTime ? (endTime || Date.now()) - startTime : 0;

  return (
    <div className="prediction-interface" style={{ position: "relative" }}>
      {dimmed && (
        <LockOverlay message="Annotation or optimization in progress... Predictions section is temporarily unavailable." />
      )}

      <div className="header">
        <h1>NMR Predictions (SMILES/CSV)</h1>
        <p>
          Generates SDF + CSV (1H/13C/COSY/HSQC) usable for annotations
          and visualization
        </p>
      </div>

      <div className="info-box">
        <strong>Integrated workflow:</strong> after a successful prediction,
        you can visualize the molecule and peaks, then launch the{" "}
        <em>Integrated Annotations</em> tab.
      </div>

      <PredictionForm
        predictionState={predictionState}
        setPredictionState={setPredictionState}
        running={running}
        csvUploads={csvUploads}
        loadingCsvUploads={loadingCsvUploads}
        refreshCsvUploads={refreshCsvUploads}
        uploadCsv={uploadCsv}
        loadMoleculeCountFromServerCsv={loadMoleculeCountFromServerCsv}
        totalMolecules={totalMolecules}
        addLog={addLog}
      />

      <div className="action-section">
        {predictionState.mode === "SMILES" ? (
          <>
            <button
              className="btn btn-primary btn-large"
              onClick={() => {
                launchSmiles();
                if (inModal && onClose) onClose();
              }}
              disabled={running}
            >
              {running ? "Prediction in progress..." : "Launch prediction (SMILES)"}
            </button>
            {running && (
              <button
                className="btn btn-danger btn-large"
                onClick={handleCancel}
                style={{ marginLeft: 8 }}
              >
                Cancel
              </button>
            )}
          </>
        ) : (
          <>
            <button
              className="btn btn-primary btn-large"
              onClick={() => {
                launchCSV();
                if (inModal && onClose) onClose();
              }}
              disabled={running || !predictionState.csvUploadedName}
            >
              {running ? "Prediction in progress..." : "Launch prediction (CSV)"}
            </button>
            {running && (
              <button
                className="btn btn-danger btn-large"
                onClick={handleCancel}
                style={{ marginLeft: 8 }}
              >
                Cancel
              </button>
            )}
          </>
        )}
        {elapsedMs > 0 && (
          <div className="duration-info">
            Duration: {formatDuration(elapsedMs)}
          </div>
        )}
      </div>

      {running && (
        <PredictionProgress
          progress={progress}
          progressMessage={progressMessage}
          completedMolecules={completedMolecules}
          totalMolecules={totalMolecules}
        />
      )}

      {error && (
        <div className="error-section">
          <div className="error-message">Error: {error}</div>
        </div>
      )}

      {!running && predictionState.artifacts?.length > 0 && (
        <PredictionResults
          artifacts={predictionState.artifacts}
          predictionTaskId={predictionState.predictionTaskId}
          showCombinedViewer={predictionState.showCombinedViewer}
          setPredictionState={setPredictionState}
          addLog={addLog}
        />
      )}
    </div>
  );
};