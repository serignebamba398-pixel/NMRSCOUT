import React from 'react';

/**
 * Progress bar with percentage and optional message
 * @param {number} progress - Progress from 0 to 100
 * @param {string} message - Optional descriptive message
 */
export const ProgressBar = ({ progress = 0, message }) => {
  const pct = Number.isFinite(progress)
    ? Math.max(0, Math.min(100, Math.floor(progress)))
    : 0;

  return (
    <div className="progress-container" style={{ position: "relative" }}>
      <div
        className="progress-bar"
        style={{ width: `${pct}%`, height: "100%" }}
        role="progressbar"
        aria-valuenow={pct}
        aria-valuemin={0}
        aria-valuemax={100}
      />
      <div className="progress-label">{pct}%</div>
      {message && <div className="progress-message">{message}</div>}
    </div>
  );
};