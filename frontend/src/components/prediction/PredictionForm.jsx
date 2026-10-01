import { isMethodAvailable } from '../../constants/releaseFeatures';
import { TocsyPredictionFields } from './TocsyPredictionFields';
import React from 'react';
import { FileUpload, ParameterInput } from '../common';
import { formatFileSize } from '../../utils';

/**
 * Prediction configuration form
 */
export const PredictionForm = ({
  predictionState,
  setPredictionState,
  running,
  csvUploads,
  loadingCsvUploads,
  refreshCsvUploads,
  uploadCsv,
  loadMoleculeCountFromServerCsv,
  totalMolecules,
  addLog,
}) => {
  const { mode, smiles, molId, csvLocalFile, csvUploadedName, options, selectedPredMethods } =
    predictionState;

  const ex = options?.extractors || {};
  const c13 = ex.c13 || {};
  const hsqc = ex.hsqc || {};
  const carbonModeUI = c13.carbon_mode || "nonquaternary";
  const hsqcFoldUI = !!hsqc.fold;
  const hsqcLowUI = Number.isFinite(+hsqc.fold_low) ? +hsqc.fold_low : 4.8618;
  const hsqcHighUI = Number.isFinite(+hsqc.fold_high) ? +hsqc.fold_high : 69.8395;

  return (
    <>
      {/* Mode selector */}
      <div className="type-selector">
        <h3>Mode</h3>
        <div className="type-tabs">
          {["SMILES", "CSV"].map((m) => (
            <button
              key={m}
              onClick={() =>
                running
                  ? null
                  : setPredictionState((prev) => ({ ...prev, mode: m }))
              }
              className={`type-tab ${mode === m ? "active" : ""}`}
              disabled={running}
            >
              <div className="tab-title">{m}</div>
            </button>
          ))}
        </div>
      </div>

      {/* SMILES Mode */}
      {mode === "SMILES" && (
        <div className="card">
          <div className="form-row">
            <label className="form-label">SMILES</label>
            <input
              className="form-input"
              type="text"
              value={smiles}
              onChange={(e) =>
                setPredictionState((prev) => ({
                  ...prev,
                  smiles: e.target.value,
                }))
              }
              disabled={running}
              placeholder="E.g.: CCO"
            />
          </div>
          <div className="form-row">
            <label className="form-label">Molecule ID</label>
            <input
              className="form-input"
              type="text"
              value={molId}
              onChange={(e) =>
                setPredictionState((prev) => ({
                  ...prev,
                  molId: e.target.value,
                }))
              }
              disabled={running}
              placeholder="ethanol"
            />
          </div>
        </div>
      )}

      {/* CSV Mode */}
      {mode === "CSV" && (
        <div className="card">
          <div className="file-uploads">
            <FileUpload
              label="CSV file (columns: id, smiles)"
              file={csvLocalFile}
              onChange={(e) => uploadCsv(e.target.files?.[0])}
              accept=".csv"
              disabled={running}
            />
          </div>

          <div className="form-row" style={{ marginTop: 12 }}>
            <label className="form-label">
              Or select an existing CSV (uploads)
            </label>
            <div style={{ display: "flex", gap: 8, alignItems: "center" }}>
              <select
                className="form-input"
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
                    addLog(`CSV selected from uploads: ${name}`);
                    await loadMoleculeCountFromServerCsv(name);
                  }
                }}
              >
                <option value="">-- Select from uploads --</option>
                {csvUploads.map((f) => (
                  <option key={f.filename} value={f.filename}>
                    {f.filename}{" "}
                    {Number.isFinite(f.size) ? `(${formatFileSize(f.size)})` : ""}{" "}
                    {f.mtime ? `- ${new Date(f.mtime).toLocaleString()}` : ""}
                  </option>
                ))}
              </select>
              <button
                className="btn btn-secondary"
                onClick={refreshCsvUploads}
                disabled={running || loadingCsvUploads}
              >
                {loadingCsvUploads ? "Loading..." : "Refresh"}
              </button>
            </div>
            <div className="hint" style={{ marginTop: 6 }}>
              Tip: no need to re-upload if the CSV already exists.
            </div>
          </div>

          {csvUploadedName && (
            <div className="file-status-summary">
              <div className="status-item success">
                CSV selected: {csvUploadedName}{" "}
                {totalMolecules ? `(${totalMolecules} molecules)` : null}
              </div>
            </div>
          )}

          <div className="parameters-grid" style={{ marginTop: 12 }}>
            <ParameterInput
              label="id column"
              value={options.idColumn}
              type="text"
              onChange={(e) =>
                setPredictionState((prev) => ({
                  ...prev,
                  options: { ...prev.options, idColumn: e.target.value },
                }))
              }
              disabled={running}
            />
            <ParameterInput
              label="smiles column"
              value={options.smilesColumn}
              type="text"
              onChange={(e) =>
                setPredictionState((prev) => ({
                  ...prev,
                  options: { ...prev.options, smilesColumn: e.target.value },
                }))
              }
              disabled={running}
            />
          </div>
        </div>
      )}

      {/* Prediction parameters */}
      <div className="card">
        <h3>Prediction methods to run</h3>
        <div className="parameters-grid">
          {["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable).map((m) => {
            const checked = selectedPredMethods.includes(m);
            return (
              <label
                key={m}
                className="parameter-input"
                style={{ flexDirection: "row", alignItems: "center", gap: 8 }}
              >
                <input
                  type="checkbox"
                  checked={checked}
                  onChange={() =>
                    setPredictionState((prev) => {
                      const set = new Set(prev.selectedPredMethods);
                      if (set.has(m)) set.delete(m);
                      else set.add(m);
                      return { ...prev, selectedPredMethods: Array.from(set) };
                    })
                  }
                  disabled={running}
                  style={{ width: 18, height: 18 }}
                />
                <span className="parameter-label" style={{ marginBottom: 0 }}>
                  {m}
                </span>
              </label>
            );
          })}
        </div>

        <div className="divider" style={{ margin: "12px 0" }} />

        <div className="parameters-grid">
          <ParameterInput
            label="nb_cores (global)"
            value={options.nb_cores}
            onChange={(e) =>
              setPredictionState((prev) => ({
                ...prev,
                options: { ...prev.options, nb_cores: e.target.value },
              }))
            }
            disabled={running}
          />
        </div>

        {/* COSY Parameters */}
        {selectedPredMethods.includes("COSY") && (
          <>
            <h4 style={{ marginTop: 12 }}>COSY Parameters</h4>
            <div className="parameters-grid">
              <ParameterInput
                label="min_bonds"
                value={options.cosy?.min_bonds || 2}
                onChange={(e) =>
                  setPredictionState((prev) => ({
                    ...prev,
                    options: {
                      ...prev.options,
                      cosy: {
                        ...(prev.options.cosy || {}),
                        min_bonds: e.target.value,
                      },
                    },
                  }))
                }
                disabled={running}
              />
              <ParameterInput
                label="max_bonds"
                value={options.cosy?.max_bonds || 3}
                onChange={(e) =>
                  setPredictionState((prev) => ({
                    ...prev,
                    options: {
                      ...prev.options,
                      cosy: {
                        ...(prev.options.cosy || {}),
                        max_bonds: e.target.value,
                      },
                    },
                  }))
                }
                disabled={running}
              />
            </div>
          </>
        )}

        {/* TOCSY Parameters */}
        {isMethodAvailable("TOCSY") && selectedPredMethods.includes("TOCSY") && (
          <>
            <h4 style={{ marginTop: 12 }}>TOCSY Parameters</h4>
            <TocsyPredictionFields
              options={options.tocsy}
              disabled={running}
              onChange={(key, value) => setPredictionState(prev => ({
                ...prev,
                options: { ...prev.options, tocsy: { ...(prev.options.tocsy || {}), [key]: value } },
              }))}
            />
          </>
        )}

        {/* 13C Parameters */}
        {selectedPredMethods.includes("13C") && (
          <>
            <h4 style={{ marginTop: 12 }}>13C Parameters</h4>
            <div className="parameters-grid">
              <div className="parameter-input">
                <label className="parameter-label">carbons to predict</label>
                <select
                  className="parameter-field"
                  value={carbonModeUI}
                  onChange={(e) =>
                    setPredictionState((prev) => ({
                      ...prev,
                      options: {
                        ...prev.options,
                        extractors: {
                          ...(prev.options?.extractors || {}),
                          c13: {
                            ...(prev.options?.extractors?.c13 || {}),
                            carbon_mode: String(e.target.value)
                              .toLowerCase()
                              .startsWith("all")
                              ? "all"
                              : "nonquaternary",
                          },
                        },
                      },
                    }))
                  }
                  disabled={running}
                >
                  <option value="all">All carbons (C, CH, CH2, CH3)</option>
                  <option value="nonquaternary">
                    Non-quaternary carbons (CH/CH2/CH3)
                  </option>
                </select>
              </div>
            </div>
          </>
        )}

        {/* HSQC Parameters */}
        {selectedPredMethods.includes("HSQC") && (
          <>
            <h4 style={{ marginTop: 12 }}>HSQC Parameters</h4>
            <div className="parameter-input" style={{ gridColumn: "1 / -1" }}>
              <label className="parameter-label">HSQC - folded spectrum?</label>
              <input
                type="checkbox"
                checked={hsqcFoldUI}
                onChange={(e) =>
                  setPredictionState((prev) => ({
                    ...prev,
                    options: {
                      ...prev.options,
                      extractors: {
                        ...(prev.options?.extractors || {}),
                        hsqc: {
                          ...(prev.options?.extractors?.hsqc || {}),
                          fold: !!e.target.checked,
                        },
                      },
                    },
                  }))
                }
                disabled={running}
              />
              <div className="hint" style={{ marginTop: 4 }}>
                If checked, the F1 dimension (13C) will be folded between the bounds below.
              </div>
            </div>

            <div className="parameters-grid">
              <ParameterInput
                label="F1 (13C) lower bound (ppm)"
                value={hsqcLowUI}
                onChange={(e) =>
                  setPredictionState((prev) => ({
                    ...prev,
                    options: {
                      ...prev.options,
                      extractors: {
                        ...(prev.options?.extractors || {}),
                        hsqc: {
                          ...(prev.options?.extractors?.hsqc || {}),
                          fold_low: +e.target.value,
                        },
                      },
                    },
                  }))
                }
                disabled={running || !hsqcFoldUI}
              />
              <ParameterInput
                label="F1 (13C) upper bound (ppm)"
                value={hsqcHighUI}
                onChange={(e) =>
                  setPredictionState((prev) => ({
                    ...prev,
                    options: {
                      ...prev.options,
                      extractors: {
                        ...(prev.options?.extractors || {}),
                        hsqc: {
                          ...(prev.options?.extractors?.hsqc || {}),
                          fold_high: +e.target.value,
                        },
                      },
                    },
                  }))
                }
                disabled={running || !hsqcFoldUI}
              />
            </div>
          </>
        )}
      </div>
    </>
  );
};