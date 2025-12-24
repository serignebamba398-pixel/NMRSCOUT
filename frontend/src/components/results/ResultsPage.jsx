import React, { useCallback, useMemo } from 'react';
import { useNavigate } from 'react-router-dom';
import { useAppState } from '../../contexts/AppStateContext';
import Viewer from '../../Viewer';
import { LockOverlay } from '../common';
import { normalizeStaticUrl, downloadBlob } from '../../utils';
import { API_BASE_URL } from '../../api';

/**
 * Annotation results display page
 */
export const ResultsPage = () => {
  const { viewerData, uiLock, customColors, setCustomColors, currentProject } = useAppState();
  const navigate = useNavigate();

  const safeViewer = viewerData ?? {
    csvs: {},
    taskId: null,
    title: "—",
    annotationResults: null,
    downloads: {},
    pdataDirs: {},
  };

  const hasCSV = safeViewer.csvs && Object.keys(safeViewer.csvs).length > 0;
  const dimmed = uiLock.lockedBy !== null;

  const normalizedCsvs = useMemo(() => {
    const csvs = safeViewer.csvs;
    if (!csvs || typeof csvs !== 'object') return {};
    return Object.fromEntries(
      Object.entries(csvs).map(([k, v]) => [k, normalizeStaticUrl(v)])
    );
  }, [safeViewer.csvs]);

  const handleDownload = async (label, url) => {
    try {
      const normalized = normalizeStaticUrl(url);
      const isStatic =
        typeof normalized === "string" && normalized.startsWith("/outputs/");
      const href = isStatic
        ? normalized
        : `${API_BASE_URL}${normalized.startsWith("/") ? "" : "/"}${normalized}`;
      const resp = await fetch(href);
      if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
      const blob = await resp.blob();
      const filename = `${label}_${safeViewer.taskId || "task"}_${Date.now()}.csv`;
      downloadBlob(blob, filename);
    } catch (e) {
      console.error("Download failed:", e.message);
    }
  };

  const handleColorChange = useCallback(
    (method, color) => {
      setCustomColors((prev) => ({ ...prev, [method]: color }));
    },
    [setCustomColors]
  );

  // Check that a project is open (after hooks to respect React rules)
  if (!currentProject) {
    return (
      <div className="results-page">
        <div className="warning-box">
          <span className="warning-icon">⚠️</span>
          <div>
            <h3 className="warning-title">No project open</h3>
            <p className="warning-text">
              You must create or open a project from the home page before viewing results.
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
    <div className="app">
      <header className="top-tabs">
        <button className="top-tab" onClick={() => navigate(-1)}>
          ← Back
        </button>
        <div className="top-tab active">Results - {safeViewer.title}</div>
      </header>

      <main className="main-content" style={{ padding: "1rem", position: "relative" }}>
        {dimmed && (
          <LockOverlay
            message={uiLock.message || "Interface locked during execution..."}
          />
        )}

        {!hasCSV ? (
          <div className="card" style={{ padding: "1rem" }}>
            <h3>No results loaded</h3>
            <p>
              Run annotations, then click <em>View results</em>.
            </p>
            <div style={{ marginTop: 12 }}>
              <button
                className="btn btn-primary"
                onClick={() => navigate("/annotations")}
              >
                Go to annotations
              </button>
            </div>
          </div>
        ) : (
          <>
            <div className="card" style={{ marginBottom: 16 }}>
              <h3>Available downloads</h3>
              <div style={{ display: "flex", gap: 8, flexWrap: "wrap" }}>
                {Object.entries(safeViewer.downloads || {}).map(([type, url]) => (
                  <button
                    key={type}
                    className="btn btn-success"
                    onClick={() => handleDownload(type, url)}
                  >
                    Download {type}
                  </button>
                ))}
              </div>
            </div>

            <div
              style={{
                background: "var(--surface-elevated)",
                borderRadius: "var(--border-radius-large)",
                padding: "1rem",
                minHeight: "70vh",
                position: "relative",
              }}
            >
              <Viewer
                csvs={normalizedCsvs}
                taskId={safeViewer.taskId}
                defaultShow3D={true}
                annotationResults={safeViewer.annotationResults}
                pdataDirs={safeViewer.pdataDirs || {}}
                bgOptions={safeViewer.bgOptions || {}}
                customColors={customColors}
                onColorChange={handleColorChange}
              />
            </div>
          </>
        )}
      </main>
    </div>
  );
};