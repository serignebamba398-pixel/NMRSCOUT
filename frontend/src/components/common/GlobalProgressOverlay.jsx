import { isMethodAvailable } from '../../constants/releaseFeatures';
import React, { useState, useRef, useCallback, useEffect } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import './GlobalProgressOverlay.css';

/**
 * Progress overlay - floating draggable and collapsible window
 * Displays when a long operation is running (prediction, annotation, optimization)
 * No blur backdrop - the page remains interactive
 */
export const GlobalProgressOverlay = () => {
  const { 
    prediction, 
    annotationState,
    annotation,
    toleranceOptimization,
  } = useAppState();

  // State for minimize/maximize
  const [isMinimized, setIsMinimized] = useState(false);
  const [isDismissed, setIsDismissed] = useState(false);
  
  // State for drag
  const [position, setPosition] = useState({ x: null, y: null });
  const [isDragging, setIsDragging] = useState(false);
  const dragOffset = useRef({ x: 0, y: 0 });
  const cardRef = useRef(null);

  // Reset position and minimized state when a new operation starts
  const predictionRunning = prediction?.running || false;
  const annotationRunning = annotationState?.allRunning || false;
  const optimizationRunning = toleranceOptimization?.running || false;
  const anyRunning = predictionRunning || annotationRunning || optimizationRunning;

  // If the user cancels a prediction from the overlay, the operation stops immediately
  // and the overlay would disappear before the user sees the final status.
  // Keep the overlay visible briefly to show the cancellation message.
  const [, forceUpdate] = useState(0);
  const predictionCancelledRecently = (() => {
    if (predictionRunning) return false;
    if (prediction?.error !== 'Cancelled by user') return false;
    const end = Number(prediction?.endTime);
    if (!Number.isFinite(end)) return false;
    return Date.now() - end < 4000;
  })();

  // Force re-render after 4 seconds to hide the "cancelled" message
  useEffect(() => {
    if (predictionCancelledRecently) {
      const timer = setTimeout(() => forceUpdate(n => n + 1), 4100);
      return () => clearTimeout(timer);
    }
  }, [predictionCancelledRecently, prediction?.endTime]);

  useEffect(() => {
    if (anyRunning) {
      // New operation: reset
      setPosition({ x: null, y: null });
      setIsMinimized(false);
      setIsDismissed(false);
    }
  }, [anyRunning]);

  useEffect(() => {
    const showProgress = () => {
      setIsDismissed(false);
      setIsMinimized(false);
    };
    window.addEventListener('nmrscout:show-progress', showProgress);
    return () => window.removeEventListener('nmrscout:show-progress', showProgress);
  }, []);

  // Drag handling
  const handleMouseDown = useCallback((e) => {
    if (e.target.closest('button') || e.target.closest('.minimize-btn')) return;
    if (!cardRef.current) return;
    
    const rect = cardRef.current.getBoundingClientRect();
    dragOffset.current = {
      x: e.clientX - rect.left,
      y: e.clientY - rect.top,
    };
    setIsDragging(true);
    e.preventDefault();
  }, []);

  const handleMouseMove = useCallback((e) => {
    if (!isDragging) return;
    
    const newX = e.clientX - dragOffset.current.x;
    const newY = e.clientY - dragOffset.current.y;
    
    // Limit to screen bounds
    const maxX = window.innerWidth - (cardRef.current?.offsetWidth || 400);
    const maxY = window.innerHeight - (cardRef.current?.offsetHeight || 200);
    
    setPosition({
      x: Math.max(0, Math.min(newX, maxX)),
      y: Math.max(0, Math.min(newY, maxY)),
    });
  }, [isDragging]);

  const handleMouseUp = useCallback(() => {
    setIsDragging(false);
  }, []);

  useEffect(() => {
    if (isDragging) {
      document.addEventListener('mousemove', handleMouseMove);
      document.addEventListener('mouseup', handleMouseUp);
      return () => {
        document.removeEventListener('mousemove', handleMouseMove);
        document.removeEventListener('mouseup', handleMouseUp);
      };
    }
  }, [isDragging, handleMouseMove, handleMouseUp]);

  // Determine which operation is running (or recently cancelled)
  let isVisible = false;
  let title = '';
  let type = null;

  if (predictionRunning) {
    isVisible = true;
    type = 'prediction';
    title = 'Prediction in progress...';
  } else if (predictionCancelledRecently) {
    isVisible = true;
    type = 'prediction';
    title = 'Prediction cancelled';
  } else if (annotationRunning) {
    isVisible = true;
    type = 'annotation';
    title = 'Annotation in progress...';
  } else if (optimizationRunning) {
    isVisible = true;
    type = 'optimization';
    title = 'Optimization in progress...'
  }

  if (!isVisible || isDismissed) return null;

  // Dynamic style for position
  const cardStyle = {
    ...(position.x !== null && position.y !== null ? {
      position: 'fixed',
      left: position.x,
      top: position.y,
      right: 'auto',
      bottom: 'auto',
      transform: 'none',
    } : {}),
    cursor: isDragging ? 'grabbing' : 'grab',
  };

  // Colors by method
  const methodColors = {
    '1H': '#1d4ed8',    // Blue
    '13C': '#166534',   // Green
    'COSY': '#b45309',  // Orange
    'TOCSY': '#7c3aed', // Purple
    'HSQC': '#a21caf',  // Pink
  };

  // For annotation, display progress bars by method
  if (type === 'annotation') {
    const progressByMethod = annotationState?.progressByMethod || {};
    const selectedMethods = annotationState?.selectedMethods || [];
    const message = annotationState?.globalProgressMessage || 'Initializing...';
    
    // DEBUG: Log to understand the issue
    console.log('[GlobalProgressOverlay] selectedMethods:', selectedMethods, 'progressByMethod keys:', Object.keys(progressByMethod));
    
    // Always show all selected methods, even if progress is 0 or undefined
    // Fallback to checking progressByMethod keys if selectedMethods is empty
    const activeMethods = selectedMethods.length > 0 
      ? selectedMethods.filter(m => ['1H', '13C', 'COSY', 'TOCSY', 'HSQC'].filter(isMethodAvailable).includes(m))
      : ['1H', '13C', 'COSY', 'TOCSY', 'HSQC'].filter(isMethodAvailable).filter(m => Object.prototype.hasOwnProperty.call(progressByMethod, m));

    // Calculate overall progress
    const totalProgress = activeMethods.length > 0 
      ? activeMethods.reduce((sum, m) => sum + (progressByMethod[m] || 0), 0) / activeMethods.length
      : 0;

    // Minimized view
    if (isMinimized) {
      return (
        <div className="global-progress-floating">
          <div 
            ref={cardRef}
            className="global-progress-minimized"
            style={cardStyle}
            onMouseDown={handleMouseDown}
          >
            <div className="minimized-content">
              <span className="minimized-title">Annotation</span>
              <span className="minimized-percent">{Math.round(totalProgress)}%</span>
              <div className="minimized-bar-container">
                <div className="minimized-bar annotation" style={{ width: `${totalProgress}%` }} />
              </div>
            </div>
            <button className="minimize-btn expand" onClick={() => setIsMinimized(false)} title="Expand">
              &#x25A1;
            </button>
            <button
              className="btn-cancel-small"
              onClick={annotation?.handleCancelAll}
              disabled={!annotationRunning || annotation?.cancelRequested}
              title="Cancel"
            >
              &#x2715;
            </button>
          </div>
        </div>
      );
    }

    return (
      <div className="global-progress-floating">
        <div 
          ref={cardRef}
          className="global-progress-card"
          style={cardStyle}
          onMouseDown={handleMouseDown}
        >
          <div className="global-progress-header">
            <h3 className="global-progress-title">{title}</h3>
            <div className="header-buttons">
              <button className="minimize-btn" onClick={() => setIsMinimized(true)} title="Minimize">
                &#x2212;
              </button>
              <div className="global-progress-spinner">
                <svg viewBox="0 0 50 50" className="spinner-svg">
                  <circle cx="25" cy="25" r="20" fill="none" strokeWidth="4" stroke="#e0e0e0" />
                  <circle
                    cx="25" cy="25" r="20" fill="none" strokeWidth="4"
                    stroke="#28a745" strokeLinecap="round"
                    strokeDasharray="125.6" strokeDashoffset="31.4"
                    className="spinner-circle"
                  />
                </svg>
              </div>
            </div>
          </div>
          
          {/* Progress bars by method */}
          <div className="method-progress-container">
            {activeMethods.map(method => {
              const progress = progressByMethod[method] || 0;
              const color = methodColors[method];
              const isComplete = progress >= 100;
              
              return (
                <div key={method} className="method-progress-row">
                  <div className="method-progress-label" style={{ color }}>
                    {method}
                  </div>
                  <div className="method-progress-bar-container">
                    <div 
                      className={`method-progress-bar ${isComplete ? 'complete' : ''}`}
                      style={{ 
                        width: `${Math.min(100, Math.max(0, progress))}%`,
                        backgroundColor: color,
                      }}
                    />
                  </div>
                  <div className="method-progress-percent">
                    {Math.round(progress)}%
                  </div>
                </div>
              );
            })}
          </div>
          
          <div className="global-progress-message-text">
            {message}
          </div>

          <div style={{ marginTop: 12, display: "flex", justifyContent: "center", gap: 8 }}>
            <button
              className="btn btn-danger"
              onClick={annotation?.handleCancelAll}
              disabled={!annotationRunning || annotation?.cancelRequested}
            >
              {annotation?.cancelRequested ? "Canceling..." : "Cancel"}
            </button>
            <button
              className="btn btn-secondary"
              onClick={() => setIsDismissed(true)}
            >
              Close
            </button>
          </div>
          
          <p className="global-progress-hint">
            You can close this window, the calculation will continue in the background.
          </p>
        </div>
      </div>
    );
  }

  if (type === 'optimization') {
    const optProgress = toleranceOptimization?.progress || 0;
    const optMessage = toleranceOptimization?.currentTest || 'Initializing...';
    const optCancelRequested = toleranceOptimization?.cancelRequested || false;

    // Minimized view
    if (isMinimized) {
      return (
        <div className="global-progress-floating">
          <div 
            ref={cardRef}
            className="global-progress-minimized"
            style={cardStyle}
            onMouseDown={handleMouseDown}
          >
            <div className="minimized-content">
              <span className="minimized-title">Optimization</span>
              <span className="minimized-percent">{Math.round(optProgress)}%</span>
              <div className="minimized-bar-container">
                <div className="minimized-bar optimization" style={{ width: `${optProgress}%` }} />
              </div>
            </div>
            <button className="minimize-btn expand" onClick={() => setIsMinimized(false)} title="Expand">
              &#x25A1;
            </button>
            <button
              className="btn-cancel-small"
              onClick={toleranceOptimization?.handleCancel}
              disabled={!optimizationRunning || optCancelRequested}
              title="Cancel"
            >
              &#x2715;
            </button>
          </div>
        </div>
      );
    }

    return (
      <div className="global-progress-floating">
        <div 
          ref={cardRef}
          className="global-progress-card"
          style={cardStyle}
          onMouseDown={handleMouseDown}
        >
          <div className="global-progress-header">
            <h3 className="global-progress-title">{title}</h3>
            <div className="header-buttons">
              <button className="minimize-btn" onClick={() => setIsMinimized(true)} title="Minimize">
                &#x2212;
              </button>
              <div className="global-progress-spinner">
                <svg viewBox="0 0 50 50" className="spinner-svg">
                  <circle cx="25" cy="25" r="20" fill="none" strokeWidth="4" stroke="#e0e0e0" />
                  <circle
                    cx="25" cy="25" r="20" fill="none" strokeWidth="4"
                    stroke="#6b7280" strokeLinecap="round"
                    strokeDasharray="125.6" strokeDashoffset="31.4"
                    className="spinner-circle"
                  />
                </svg>
              </div>
            </div>
          </div>

          <div className="global-progress-bar-container">
            <div
              className="global-progress-bar optimization"
              style={{ width: `${Math.min(100, Math.max(0, optProgress))}%` }}
            />
          </div>

          <div className="global-progress-info">
            <span className="global-progress-percent">{Math.round(optProgress)}%</span>
            <span className="global-progress-message">{optMessage}</span>
          </div>

          <div style={{ marginTop: 12, display: "flex", justifyContent: "center" }}>
            <button
              className="btn btn-danger"
              onClick={toleranceOptimization?.handleCancel}
              disabled={!optimizationRunning || optCancelRequested}
            >
              {optCancelRequested ? "Canceling..." : "Cancel"}
            </button>
          </div>

          <p className="global-progress-hint">
            You can close this window, the calculation will continue in the background.
          </p>
        </div>
      </div>
    );
  }

  // For prediction, display a single bar
  const progress = prediction?.progress || 0;
  const message = predictionCancelledRecently
    ? (prediction?.error || 'Cancelled by user')
    : (prediction?.progressMessage || 'Initializing...');
  const completedMolecules = prediction?.completedMolecules || 0;
  const totalMolecules = prediction?.totalMolecules || 0;
  
  // Prefer the backend/hook message; only fall back to a counter if missing.
  const displayMessage = message || (
    totalMolecules > 0 ? `Molecule ${completedMolecules}/${totalMolecules} (SDF)` : 'Initializing...'
  );

  // Minimized view
  if (isMinimized) {
    return (
      <div className="global-progress-floating">
        <div 
          ref={cardRef}
          className="global-progress-minimized"
          style={cardStyle}
          onMouseDown={handleMouseDown}
        >
          <div className="minimized-content">
            <span className="minimized-title">Prediction</span>
            <span className="minimized-percent">{Math.round(progress)}%</span>
            <div className="minimized-bar-container">
              <div className="minimized-bar prediction" style={{ width: `${progress}%` }} />
            </div>
          </div>
          <button className="minimize-btn expand" onClick={() => setIsMinimized(false)} title="Expand">
            &#x25A1;
          </button>
          <button
            className="btn-cancel-small"
            onClick={prediction?.handleCancel}
            disabled={!predictionRunning}
            title="Cancel"
          >
            &#x2715;
          </button>
        </div>
      </div>
    );
  }

  return (
    <div className="global-progress-floating">
      <div 
        ref={cardRef}
        className="global-progress-card"
        style={cardStyle}
        onMouseDown={handleMouseDown}
      >
        <div className="global-progress-header">
          <h3 className="global-progress-title">{title}</h3>
          <div className="header-buttons">
            <button className="minimize-btn" onClick={() => setIsMinimized(true)} title="Minimize">
              &#x2212;
            </button>
            <div className="global-progress-spinner">
              <svg viewBox="0 0 50 50" className="spinner-svg">
                <circle cx="25" cy="25" r="20" fill="none" strokeWidth="4" stroke="#e0e0e0" />
                <circle
                  cx="25" cy="25" r="20" fill="none" strokeWidth="4"
                  stroke="#4a90d9" strokeLinecap="round"
                  strokeDasharray="125.6" strokeDashoffset="31.4"
                  className="spinner-circle"
                />
              </svg>
            </div>
          </div>
        </div>
        
        <div className="global-progress-bar-container">
          <div 
            className="global-progress-bar prediction"
            style={{ width: `${Math.min(100, Math.max(0, progress))}%` }}
          />
        </div>
        
        <div className="global-progress-info">
          <span className="global-progress-percent">{Math.round(progress)}%</span>
          <span className="global-progress-message">{displayMessage}</span>
        </div>

        <div style={{ marginTop: 12, display: "flex", justifyContent: "center", gap: 8 }}>
          <button
            className="btn btn-danger"
            onClick={prediction?.handleCancel}
            disabled={!predictionRunning}
          >
            Cancel
          </button>
          <button
            className="btn btn-secondary"
            onClick={() => setIsDismissed(true)}
          >
            Close
          </button>
        </div>
        
        <p className="global-progress-hint">
          You can close this window, the calculation will continue in the background.
        </p>
      </div>
    </div>
  );
};

// For optimization, display a single bar


export default GlobalProgressOverlay;
