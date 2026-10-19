import React from 'react';
import { useLogs } from '../../contexts/AppStateContext';

/**
 * Controls for showing, exporting, and clearing logs.
 */
export const LogControls = () => {
  const { logs, showLogs, setShowLogs, clearLogs, addLog } = useLogs();

  const isElectron = typeof window !== 'undefined' && window.electronAPI?.isElectron === true;

  const handleExportLogs = async () => {
    if (!isElectron || !window.electronAPI?.exportLogs) {
      addLog("Export logs unavailable (browser mode)", "error");
      return;
    }
    addLog("Exporting logs...", "info");
    try {
      const res = await window.electronAPI.exportLogs(logs);
      if (res?.success && res?.filePath) {
        addLog(`Logs exported: ${res.filePath}`, "success");
      } else if (res?.canceled) {
        addLog("Log export canceled.", "info");
      } else {
        addLog(`Export logs failed (${res?.error || 'unknown reason'})`, "error");
      }
    } catch (e) {
      addLog(`Export logs error (${e?.message || e})`, "error");
    }
  };

  return (
    <div className="log-controls">
      <button
        onClick={() => setShowLogs(!showLogs)}
        className="btn btn-secondary"
        aria-expanded={showLogs}
        aria-controls="log-console"
      >
        {showLogs ? "Hide logs" : "View logs"} ({logs.length})
      </button>
      <button
        onClick={handleExportLogs}
        className="btn btn-secondary"
        disabled={!isElectron}
        title={isElectron ? "Export logs (zip)" : "Available only in the installed app"}
      >
        Export logs
      </button>
      {logs.length > 0 && (
        <button 
          className="btn btn-secondary" 
          onClick={clearLogs}
          aria-label="Clear all logs"
        >
          Clear logs
        </button>
      )}
    </div>
  );
};