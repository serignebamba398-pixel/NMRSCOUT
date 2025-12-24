import React, { useState, useCallback } from 'react';
import { AppStateProvider, useAppState } from './contexts/AppStateContext';
import { Toolbar, Modal, LogControls } from './components/layout';
import { LogConsole, GlobalProgressOverlay } from './components/common';
import { MainView } from './components/main';
import { TabbedPredictionInterface } from './components/prediction';
import { TabbedToleranceInterface } from './components/optimization';
import { TabbedAnnotationInterface } from './components/annotation';
import './App.css';

/**
 * Main application content - SIRIUS-style interface
 */
const AppContent = () => {
  const { 
    currentProject, 
    setCurrentProject,
    predictionState,
    annotationState,
    viewerData,
    customColors,
    prediction,
    annotation,
    addLog,
    setShowLogs,
  } = useAppState();
  const [loading, setLoading] = useState(false);
  
  // Modal states
  const [showPrediction, setShowPrediction] = useState(false);
  const [showAnnotation, setShowAnnotation] = useState(false);
  const [showTolerance, setShowTolerance] = useState(false);

  // Handlers projet - Mode Electron
  const handleNewProject = useCallback(async () => {
    setLoading(true);
    try {
      const project = await window.electronAPI.createProject();
      if (project) {
        setCurrentProject(project);
      }
    } catch (err) {
      console.error('Project creation error:', err);
    } finally {
      setLoading(false);
    }
  }, [setCurrentProject]);

  const handleOpenProject = useCallback(async () => {
    setLoading(true);
    try {
      const project = await window.electronAPI.openProject();
      if (project) {
        setCurrentProject(project);
      }
    } catch (err) {
      console.error('Project opening error:', err);
    } finally {
      setLoading(false);
    }
  }, [setCurrentProject]);

  const handleSaveProject = useCallback(async () => {
    if (!currentProject) return;
    setLoading(true);
    try {
      const state = {
        predictionState,
        annotationState,
        viewerData,
        customColors,
      };
      const r = await window.electronAPI.saveProjectState(state);
      if (!r?.success) {
        console.warn('Project save: failed', r?.error);
      } else {
        console.log('Project saved');
      }
    } catch (err) {
      console.error('Project save error:', err);
    } finally {
      setLoading(false);
    }
  }, [
    currentProject,
    predictionState,
    annotationState,
    viewerData,
    customColors,
  ]);

  const handleSaveAs = useCallback(async () => {
    if (!currentProject) return;
    setLoading(true);
    try {
      const state = {
        predictionState,
        annotationState,
        viewerData,
        customColors,
      };
      const result = await window.electronAPI.saveProjectAs(state);
      if (result) {
        console.log('Project saved as:', result);
        setCurrentProject(result);
      }
    } catch (err) {
      console.error('Save as error:', err);
    } finally {
      setLoading(false);
    }
  }, [
    currentProject,
    predictionState,
    annotationState,
    viewerData,
    customColors,
    setCurrentProject,
  ]);

  // ========== HANDLERS COMPUTE ==========
  const handleComputePrediction = useCallback(() => {
    if (!prediction?.runPrediction) {
      addLog("Prediction not available. Configure parameters first.", "error");
      return;
    }
    setShowLogs(true);
    prediction.runPrediction();
  }, [prediction, addLog, setShowLogs]);

  const handleComputeAnnotation = useCallback(() => {
    if (!annotation?.runSelectedIntegrated) {
      addLog("Annotation not available. Configure parameters first.", "error");
      return;
    }
    setShowLogs(true);
    annotation.runSelectedIntegrated();
  }, [annotation, addLog, setShowLogs]);

  const handleComputeOptimization = useCallback(() => {
    // TODO: Implement when optimization hook is ready
    addLog("Optimization - To be implemented", "info");
  }, [addLog]);

  return (
    <div className="app sirius-layout">
      {/* Overlay de progression globale */}
      <GlobalProgressOverlay />

      {/* Toolbar en haut */}
      <Toolbar
        onNewProject={handleNewProject}
        onOpenProject={handleOpenProject}
        onSaveProject={handleSaveProject}
        onSaveAs={handleSaveAs}
        onPrediction={() => setShowPrediction(true)}
        onAnnotation={() => setShowAnnotation(true)}
        onTolerance={() => setShowTolerance(true)}
        onComputePrediction={handleComputePrediction}
        onComputeAnnotation={handleComputeAnnotation}
        onComputeOptimization={handleComputeOptimization}
        currentProject={currentProject}
        loading={loading}
      />

      {/* Main view with molecules list + viewer */}
      <MainView />

      {/* Prediction Modal - Parameters only */}
      <Modal
        isOpen={showPrediction}
        onClose={() => setShowPrediction(false)}
        title="Prediction Parameters"
        size="large"
      >
        <TabbedPredictionInterface inModal onClose={() => setShowPrediction(false)} />
      </Modal>

      {/* Annotation Modal - Parameters only */}
      <Modal
        isOpen={showAnnotation}
        onClose={() => setShowAnnotation(false)}
        title="Annotation Parameters"
        size="large"
      >
        <TabbedAnnotationInterface inModal onClose={() => setShowAnnotation(false)} />
      </Modal>

      {/* Tolerance Modal - Parameters only */}
      <Modal
        isOpen={showTolerance}
        onClose={() => setShowTolerance(false)}
        title="Tolerance Optimization Parameters"
        size="large"
      >
        <TabbedToleranceInterface onClose={() => setShowTolerance(false)} />
      </Modal>

      {/* Project creation dialog (web/Docker mode) */}
      {/* Logs */}
      <LogControls />
      <LogConsole />
    </div>
  );
};

/**
 * Composant racine de l'application
 */
const App = () => (
  <AppStateProvider>
    <AppContent />
  </AppStateProvider>
);

export default App;