import React, { useRef, useEffect } from 'react';
import { useLogs } from '../../contexts/AppStateContext';

/**
 * Log console with auto-scroll and type filtering
 */
export const LogConsole = () => {
  const { logs, showLogs, clearLogs } = useLogs();
  const contentRef = useRef(null);

  // Auto-scroll to bottom when new log arrives
  useEffect(() => {
    if (contentRef.current) {
      contentRef.current.scrollTop = contentRef.current.scrollHeight;
    }
  }, [logs]);

  if (!showLogs) return null;

  return (
    <div className="log-console">
      <div className="log-console-header">
        <h3>Log Console</h3>
        <button 
          onClick={clearLogs} 
          className="btn btn-small btn-danger"
          aria-label="Clear logs"
        >
          Clear
        </button>
      </div>
      
      <div ref={contentRef} className="log-console-content">
        {logs.length === 0 ? (
          <div className="log-empty">No logs yet...</div>
        ) : (
          logs.map((log, i) => (
            <div key={i} className={`log-entry log-${log.type}`}>
              <span className="log-timestamp">[{log.timestamp}]</span>{" "}
              <span className="log-message">{log.message}</span>
            </div>
          ))
        )}
      </div>
    </div>
  );
};