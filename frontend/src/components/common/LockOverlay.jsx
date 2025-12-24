import React from 'react';

/**
 * Overlay that blocks the interface with a spinner and message
 * @param {string} message - Message to display
 */
export const LockOverlay = ({ message }) => {
  return (
    <div className="dimmer-overlay" aria-live="polite" role="alert">
      <div className="dimmer-content">
        <span className="spinner" aria-hidden="true" />
        <span>{message || "Operation in progress..."}</span>
      </div>
    </div>
  );
};