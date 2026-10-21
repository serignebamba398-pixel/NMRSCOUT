import React from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import './Toolbar.css';

// Toolbar SVG icons.
const Icons = {
  newFile: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <path d="M14 2H6a2 2 0 0 0-2 2v16a2 2 0 0 0 2 2h12a2 2 0 0 0 2-2V8z"/>
      <polyline points="14 2 14 8 20 8"/>
      <line x1="12" y1="18" x2="12" y2="12"/>
      <line x1="9" y1="15" x2="15" y2="15"/>
    </svg>
  ),
  openFolder: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <path d="M22 19a2 2 0 0 1-2 2H4a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h5l2 3h9a2 2 0 0 1 2 2z"/>
    </svg>
  ),
  save: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <path d="M19 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h11l5 5v11a2 2 0 0 1-2 2z"/>
      <polyline points="17 21 17 13 7 13 7 21"/>
      <polyline points="7 3 7 8 15 8"/>
    </svg>
  ),
  saveAs: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <path d="M19 21H5a2 2 0 0 1-2-2V5a2 2 0 0 1 2-2h11l5 5v11a2 2 0 0 1-2 2z"/>
      <polyline points="17 21 17 13 7 13 7 21"/>
      <polyline points="7 3 7 8 15 8"/>
      <circle cx="19" cy="19" r="4" fill="white" stroke="currentColor"/>
      <path d="M19 17v4M17 19h4" strokeWidth="1.5"/>
    </svg>
  ),
  prediction: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <circle cx="12" cy="12" r="10"/>
      <path d="M12 6v6l4 2"/>
      <path d="M8 14l-2 4M16 14l2 4"/>
    </svg>
  ),
  annotation: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <path d="M11 4H4a2 2 0 0 0-2 2v14a2 2 0 0 0 2 2h14a2 2 0 0 0 2-2v-7"/>
      <path d="M18.5 2.5a2.121 2.121 0 0 1 3 3L12 15l-4 1 1-4 9.5-9.5z"/>
    </svg>
  ),
  settings: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2">
      <circle cx="12" cy="12" r="3"/>
      <path d="M19.4 15a1.65 1.65 0 0 0 .33 1.82l.06.06a2 2 0 0 1 0 2.83 2 2 0 0 1-2.83 0l-.06-.06a1.65 1.65 0 0 0-1.82-.33 1.65 1.65 0 0 0-1 1.51V21a2 2 0 0 1-2 2 2 2 0 0 1-2-2v-.09A1.65 1.65 0 0 0 9 19.4a1.65 1.65 0 0 0-1.82.33l-.06.06a2 2 0 0 1-2.83 0 2 2 0 0 1 0-2.83l.06-.06a1.65 1.65 0 0 0 .33-1.82 1.65 1.65 0 0 0-1.51-1H3a2 2 0 0 1-2-2 2 2 0 0 1 2-2h.09A1.65 1.65 0 0 0 4.6 9a1.65 1.65 0 0 0-.33-1.82l-.06-.06a2 2 0 0 1 0-2.83 2 2 0 0 1 2.83 0l.06.06a1.65 1.65 0 0 0 1.82.33H9a1.65 1.65 0 0 0 1-1.51V3a2 2 0 0 1 2-2 2 2 0 0 1 2 2v.09a1.65 1.65 0 0 0 1 1.51 1.65 1.65 0 0 0 1.82-.33l.06-.06a2 2 0 0 1 2.83 0 2 2 0 0 1 0 2.83l-.06.06a1.65 1.65 0 0 0-.33 1.82V9a1.65 1.65 0 0 0 1.51 1H21a2 2 0 0 1 2 2 2 2 0 0 1-2 2h-.09a1.65 1.65 0 0 0-1.51 1z"/>
    </svg>
  ),
  play: (
    <svg viewBox="0 0 24 24" fill="currentColor" stroke="none">
      <polygon points="5 3 19 12 5 21 5 3"/>
    </svg>
  ),
  spinner: (
    <svg viewBox="0 0 24 24" fill="none" stroke="currentColor" strokeWidth="2.5" className="spinner-icon">
      <circle cx="12" cy="12" r="10" strokeOpacity="0.25"/>
      <path d="M12 2a10 10 0 0 1 10 10" strokeLinecap="round"/>
    </svg>
  ),
};

/**
 * SIRIUS-style toolbar
 * Displays project icons and action buttons
 */
export const Toolbar = ({
  onNewProject,
  onOpenProject,
  onSaveProject,
  onSaveAs,
  onPrediction,
  onAnnotation,
  onComputePrediction,
  onComputeAnnotation,
  currentProject,
  loading = false,
}) => {
  // Read task execution state from context.
  const { prediction, annotationState } = useAppState();
  
  const predictionRunning = prediction?.running || false;
  const annotationRunning = annotationState?.allRunning || false;
  const anyRunning = predictionRunning || annotationRunning;

  return (
    <div className="sirius-toolbar">
      {/* Section Projet */}
      <div className="toolbar-section">
        <button 
          className="toolbar-btn" 
          data-action="new"
          onClick={onNewProject}
          disabled={loading || anyRunning}
          title="New project"
        >
          <span className="toolbar-icon">{Icons.newFile}</span>
          <span className="toolbar-label">New</span>
        </button>
        
        <button 
          className="toolbar-btn" 
          data-action="open"
          onClick={onOpenProject}
          disabled={loading || anyRunning}
          title="Open project"
        >
          <span className="toolbar-icon">{Icons.openFolder}</span>
          <span className="toolbar-label">Open</span>
        </button>
        
        <button 
          className="toolbar-btn" 
          data-action="save"
          onClick={onSaveProject}
          disabled={loading || !currentProject || anyRunning}
          title="Sauvegarder le projet"
        >
          <span className="toolbar-icon">{Icons.save}</span>
          <span className="toolbar-label">Save</span>
        </button>
        
        <button 
          className="toolbar-btn" 
          data-action="saveas"
          onClick={onSaveAs}
          disabled={loading || !currentProject || anyRunning}
          title="Sauvegarder sous..."
        >
          <span className="toolbar-icon">{Icons.saveAs}</span>
          <span className="toolbar-label">Save As</span>
        </button>
      </div>

      <div className="toolbar-separator" />

      {/* Settings dialogs. */}
      <div className="toolbar-section">
        <button 
          className="toolbar-btn" 
          data-action="prediction"
          onClick={onPrediction}
          disabled={loading || anyRunning}
          title="Paramètres de prédiction"
        >
          <span className="toolbar-icon">{Icons.prediction}</span>
          <span className="toolbar-label">Prediction</span>
        </button>
        
        <button 
          className="toolbar-btn" 
          data-action="annotation"
          onClick={onAnnotation}
          disabled={loading || anyRunning}
          title="Paramètres d'annotation"
        >
          <span className="toolbar-icon">{Icons.annotation}</span>
          <span className="toolbar-label">Annotation</span>
        </button>
      </div>

      <div className="toolbar-separator" />

      {/* Task launch controls. */}
      <div className="toolbar-section toolbar-compute-section">
        <button 
          className={`toolbar-btn compute-btn ${predictionRunning ? 'running' : ''}`}
          data-action="compute-prediction"
          onClick={() => predictionRunning
            ? window.dispatchEvent(new Event('nmrscout:show-progress'))
            : onComputePrediction?.()}
          disabled={loading || annotationRunning}
          title={predictionRunning ? "Show prediction progress" : "Start prediction"}
        >
          <span className={`toolbar-icon compute-icon ${predictionRunning ? 'spinning' : ''}`}>
            {predictionRunning ? Icons.spinner : Icons.play}
          </span>
          <span className="toolbar-label">
            {predictionRunning ? 'Running...' : 'Compute Pred.'}
          </span>
        </button>
        
        <button 
          className={`toolbar-btn compute-btn ${annotationRunning ? 'running' : ''}`}
          data-action="compute-annotation"
          onClick={() => annotationRunning
            ? window.dispatchEvent(new Event('nmrscout:show-progress'))
            : onComputeAnnotation?.()}
          disabled={loading || predictionRunning}
          title={annotationRunning ? "Show annotation progress" : "Start annotation"}
        >
          <span className={`toolbar-icon compute-icon ${annotationRunning ? 'spinning' : ''}`}>
            {annotationRunning ? Icons.spinner : Icons.play}
          </span>
          <span className="toolbar-label">
            {annotationRunning ? 'Running...' : 'Compute Anno.'}
          </span>
        </button>
        
      </div>

      <div className="toolbar-spacer" />

      {/* Info projet */}
      <div className="toolbar-project-info">
        {currentProject ? (
          <span className="project-name" title={currentProject.name || currentProject.workDir}>
            {currentProject.name || 'Projet sans nom'}
          </span>
        ) : (
          <span className="no-project">Aucun projet ouvert</span>
        )}
      </div>
    </div>
  );
};

export default Toolbar;
