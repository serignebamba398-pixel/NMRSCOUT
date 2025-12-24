import React, { useMemo } from 'react';
import Viewer from '../../Viewer';
import { normalizeStaticUrl, downloadBlob } from '../../utils';
import { API_BASE_URL } from '../../api';

/**
 * Prediction results display
 */
export const PredictionResults = ({
  artifacts,
  predictionTaskId,
  showCombinedViewer,
  setPredictionState,
  addLog,
}) => {
  const sdfArtifacts = artifacts.filter((a) => a.url?.endsWith(".sdf"));

  const csvArtifactToNmrType = (t) =>
    t === "CSV_1H"
      ? "1H"
      : t === "CSV_13C"
      ? "13C"
      : t === "CSV_COSY"
      ? "COSY"
      : t === "CSV_HSQC"
      ? "HSQC"
      : null;

  const combinedCSVs = useMemo(() => {
    const csvs = {};
    artifacts
      .filter((a) => a.type.startsWith("CSV_"))
      .forEach((a) => {
        const nmrType = csvArtifactToNmrType(a.type);
        if (nmrType) csvs[nmrType] = normalizeStaticUrl(a.url);
      });
    return csvs;
  }, [artifacts]);

  const downloadArtifact = async (a) => {
    try {
      const normalized = normalizeStaticUrl(a.url);
      const isStatic =
        typeof normalized === "string" && normalized.startsWith("/outputs/");
      const href = isStatic
        ? normalized
        : `${API_BASE_URL}${normalized.startsWith("/") ? "" : "/"}${normalized}`;
      const resp = await fetch(href);
      if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
      const blob = await resp.blob();
      const filename = `${a.type}_${a.file || "artifact"}`;
      downloadBlob(blob, filename);
      addLog(`Download: ${filename}`);
    } catch (e) {
      addLog(`Download error: ${e.message}`, "error");
    }
  };

  return (
    <>
      <div className="results-section">
        <h4>Generated Artifacts</h4>
        <div className="artifact-list">
          {artifacts.map((a, i) => (
            <div key={`${a.type}-${i}`} className="artifact-item">
              <div className="artifact-info">
                <div className="artifact-type">{a.type}</div>
                <div className="artifact-file">{a.file}</div>
              </div>
              <div className="artifact-actions">
                <a
                  className="btn btn-small"
                  href={normalizeStaticUrl(a.url)}
                  target="_blank"
                  rel="noreferrer"
                  title="Open"
                >
                  Open
                </a>
                <button className="btn btn-small" onClick={() => downloadArtifact(a)}>
                  Download
                </button>
              </div>
            </div>
          ))}
        </div>

        {sdfArtifacts.length > 0 && (
          <div className="hint">
            SDF files generated: {sdfArtifacts.length}. Download them for inspection (RDKit,
            Avogadro...).
          </div>
        )}
      </div>

      {Object.keys(combinedCSVs).length > 0 && (
        <div className="results-section">
          <div
            style={{
              display: "flex",
              alignItems: "center",
              justifyContent: "space-between",
              marginBottom: 12,
            }}
          >
            <h4>Combined visualization of predicted spectra</h4>
            <button
              className="btn btn-secondary"
              onClick={() =>
                setPredictionState((prev) => ({
                  ...prev,
                  showCombinedViewer: !prev.showCombinedViewer,
                }))
              }
            >
              {showCombinedViewer
                ? "Hide visualization"
                : "Show visualization"}
            </button>
          </div>

          {showCombinedViewer && (
            <Viewer
              csvs={combinedCSVs}
              taskId={predictionTaskId}
              defaultShow3D
              pdataDirs={{}}
            />
          )}
        </div>
      )}
    </>
  );
};