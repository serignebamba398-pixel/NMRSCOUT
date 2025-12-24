import React from 'react';
import './Modal.css';

/**
 * Composant Modal générique style SIRIUS
 * Utilisé pour Prediction, Annotation, Tolerance
 */
export const Modal = ({
  isOpen,
  onClose,
  title,
  children,
  size = 'large', // 'small', 'medium', 'large', 'fullscreen'
  showCloseButton = true,
}) => {
  if (!isOpen) return null;

  const handleBackdropClick = (e) => {
    if (e.target === e.currentTarget) {
      onClose();
    }
  };

  return (
    <div className="modal-backdrop" onClick={handleBackdropClick}>
      <div className={`modal-container modal-${size}`}>
        <div className="modal-header">
          <h2 className="modal-title">{title}</h2>
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
        <div className="modal-content">
          {children}
        </div>
      </div>
    </div>
  );
};

export default Modal;
