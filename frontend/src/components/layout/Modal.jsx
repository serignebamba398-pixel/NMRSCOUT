import React, { useState, useCallback } from 'react';
import './Modal.css';

/**
 * Shared modal for prediction, annotation, and tolerance settings.
 */
export const Modal = ({
  isOpen,
  onClose,
  title,
  children,
  size = 'large', // 'small', 'medium', 'large', 'fullscreen'
  showCloseButton = true,
  allowFullscreen = true, // Whether to show the fullscreen button.
}) => {
  const [isFullscreen, setIsFullscreen] = useState(false);

  const toggleFullscreen = useCallback(() => {
    setIsFullscreen(prev => !prev);
  }, []);

  if (!isOpen) return null;

  const handleBackdropClick = (e) => {
    if (e.target === e.currentTarget) {
      onClose();
    }
  };

  // Choose the size class for the current display mode.
  const effectiveSize = isFullscreen ? 'fullscreen' : size;

  return (
    <div className="modal-backdrop" onClick={handleBackdropClick}>
      <div className={`modal-container modal-${effectiveSize}`}>
        <div className="modal-header">
          <h2 className="modal-title">{title}</h2>
          <div className="modal-header-buttons">
            {allowFullscreen && (
              <button 
                className="modal-fullscreen-btn" 
                onClick={toggleFullscreen}
                title={isFullscreen ? "Réduire" : "Plein écran"}
              >
                {isFullscreen ? (
                  // Exit fullscreen icon.
                  <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor">
                    <path d="M5.5 0a.5.5 0 0 1 .5.5v4A1.5 1.5 0 0 1 4.5 6h-4a.5.5 0 0 1 0-1h4a.5.5 0 0 0 .5-.5v-4a.5.5 0 0 1 .5-.5zm5 0a.5.5 0 0 1 .5.5v4a.5.5 0 0 0 .5.5h4a.5.5 0 0 1 0 1h-4A1.5 1.5 0 0 1 10 4.5v-4a.5.5 0 0 1 .5-.5zM0 10.5a.5.5 0 0 1 .5-.5h4A1.5 1.5 0 0 1 6 11.5v4a.5.5 0 0 1-1 0v-4a.5.5 0 0 0-.5-.5h-4a.5.5 0 0 1-.5-.5zm10 0a.5.5 0 0 1 .5-.5h4a.5.5 0 0 1 0 1h-4a.5.5 0 0 0-.5.5v4a.5.5 0 0 1-1 0v-4z"/>
                  </svg>
                ) : (
                  // Enter fullscreen icon.
                  <svg width="16" height="16" viewBox="0 0 16 16" fill="currentColor">
                    <path d="M1.5 1a.5.5 0 0 0-.5.5v4a.5.5 0 0 1-1 0v-4A1.5 1.5 0 0 1 1.5 0h4a.5.5 0 0 1 0 1h-4zM10 .5a.5.5 0 0 1 .5-.5h4A1.5 1.5 0 0 1 16 1.5v4a.5.5 0 0 1-1 0v-4a.5.5 0 0 0-.5-.5h-4a.5.5 0 0 1-.5-.5zM.5 10a.5.5 0 0 1 .5.5v4a.5.5 0 0 0 .5.5h4a.5.5 0 0 1 0 1h-4A1.5 1.5 0 0 1 0 14.5v-4a.5.5 0 0 1 .5-.5zm15 0a.5.5 0 0 1 .5.5v4a1.5 1.5 0 0 1-1.5 1.5h-4a.5.5 0 0 1 0-1h4a.5.5 0 0 0 .5-.5v-4a.5.5 0 0 1 .5-.5z"/>
                  </svg>
                )}
              </button>
            )}
            {showCloseButton && (
              <button 
                className="modal-close-btn" 
                onClick={onClose}
                title="Fermer"
              >
                ✕
              </button>
            )}
          </div>
        </div>
        <div className="modal-content">
          {children}
        </div>
      </div>
    </div>
  );
};

export default Modal;
