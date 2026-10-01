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
      addLog("Export logs indisponible (mode navigateur)", "error");
      return;
    }
    addLog("Export des logs en cours…", "info");
    try {
      const res = await window.electronAPI.exportLogs(logs);
      if (res?.success && res?.filePath) {
        addLog(`Logs exportés: ${res.filePath}`, "success");
      } else if (res?.canceled) {
        addLog("Export des logs annulé.", "info");
      } else {
        addLog(`Export logs: échec (${res?.error || 'raison inconnue'})`, "error");
      }
    } catch (e) {
      addLog(`Export logs: erreur (${e?.message || e})`, "error");
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
        {showLogs ? "Masquer logs" : "Voir logs"} ({logs.length})
      </button>
      <button
        onClick={handleExportLogs}
        className="btn btn-secondary"
        disabled={!isElectron}
        title={isElectron ? "Exporter les logs (zip)" : "Disponible uniquement dans l'app installée"}
      >
        Export logs
      </button>
      {logs.length > 0 && (
        <button 
          className="btn btn-secondary" 
          onClick={clearLogs}
          aria-label="Vider tous les logs"
        >
          Vider logs
        </button>
      )}
    </div>
  );
};