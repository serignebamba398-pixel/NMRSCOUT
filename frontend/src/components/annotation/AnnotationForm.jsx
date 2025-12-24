import React from 'react';
import { FileUpload, ParameterInput } from '../common';
import { formatFileSize } from '../../utils';

/**
 * Annotation configuration form
 */
export const AnnotationForm = ({
  annotationState,
  setAnnotationState,
  allRunning,
  hasCsvForMethod,
  matrixUploads,
  loadingMatrixUploads,
  refreshMatrixUploads,
  uploadMatrix,
  addLog,
}) => {
  const { activeType, uploadedMatrixNames, matrixFiles, uploadingMatrix, parameters } =
    annotationState;

  const updateParameter = (type, key, value) => {
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
  };

  // Tabs to select active type
  return (
    <>
      <div className="type-selector">
        <h3>Annotation Type (parameters to edit)</h3>
        <div className="type-tabs">
          {["1H", "13C", "COSY", "HSQC"].map((t) => {
            const available = hasCsvForMethod(t);
            return (
              <button
                key={t}
                onClick={() =>
                  allRunning
                    ? null
                    : setAnnotationState((prev) => ({
                        ...prev,
                        activeType: t,
                        readyToView: false,
                      }))
                }
                className={`type-tab ${activeType === t ? "active" : ""} ${
                  !available ? "disabled" : ""
                }`}
                disabled={allRunning || !available}
              >
                <div className="tab-title">{t}</div>
                <div className="tab-description">
                  {available ? "Available" : "Invalid / missing"}
                </div>
              </button>
            );
          })}
        </div>
      </div>

      {/* Matrix upload */}
      <div className="card">
        <h3>Spectral Matrix</h3>

        <div className="form-row">
          <label className="form-label">
            Select an existing matrix (uploads)
          </label>
          <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
            <select
              className="form-input"
              disabled={allRunning || loadingMatrixUploads}
              value={uploadedMatrixNames[activeType] || ""}
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
                  readyToView: false,
                }));
                if (name) addLog(`[${activeType}] matrix selected: ${name}`);
              }}
            >
              <option value="">-- Select from uploads --</option>
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
              disabled={allRunning || loadingMatrixUploads}
            >
              {loadingMatrixUploads ? "Loading..." : "Refresh"}
            </button>
          </div>
          <div className="hint" style={{ marginTop: 6 }}>
            You can also upload a new file below.
          </div>
        </div>

        <FileUpload
          label={`Matrix file for ${activeType}`}
          file={matrixFiles[activeType] || null}
          onChange={(e) => uploadMatrix(activeType, e.target.files?.[0])}
          accept=".csv"
          disabled={allRunning}
          uploading={!!uploadingMatrix[activeType]}
        />
        {uploadedMatrixNames[activeType] && (
          <div className="file-status-summary">
            <div className="status-item success">
              Matrice {activeType}: {uploadedMatrixNames[activeType]}
            </div>
          </div>
        )}
      </div>

      {/* Annotation parameters */}
      <div className="parameters-section">
        <h3>Parameters {activeType}</h3>
        <div className="parameters-grid">
          {Object.entries(parameters[activeType] || {}).map(([k, v]) => (
            <ParameterInput
              key={k}
              label={k.replace(/_/g, " ")}
              value={v}
              onChange={(e) => updateParameter(activeType, k, e.target.value)}
              disabled={allRunning}
            />
          ))}
        </div>
      </div>
    </>
  );
};