import React, { useState, useCallback } from 'react';
import { AppStateProvider, useAppState } from './contexts/AppStateContext';
import { Toolbar, Modal, LogControls } from './components/layout';
import { LogConsole, GlobalProgressOverlay, PendingTasksNotification } from './components/common';
import { MainView } from './components/main';
import { TabbedPredictionInterface } from './components/prediction';
import { TabbedToleranceInterface } from './components/optimization';
import { TabbedAnnotationInterface } from './components/annotation';
import { API_BASE_URL } from './api';
import './App.css';

const IS_ELECTRON = window.electronAPI?.isElectron === true;

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
  // Web mode: project list dialog
  const [showProjectList, setShowProjectList] = useState(false);
  const [projectList, setProjectList] = useState([]);

  // ========== PROJECT HANDLERS ==========
  const handleNewProject = useCallback(async () => {
    if (IS_ELECTRON) {
      setLoading(true);
      try {
        const project = await window.electronAPI.createProject();
        if (project) setCurrentProject(project);
      } catch (err) {
        console.error('Project creation error:', err);
      } finally {
        setLoading(false);
      }
      return;
    }
    // Web/Docker mode: prompt for name, create via API
    const name = prompt('Project name:');
    if (!name) return;
    setLoading(true);
    try {
      const res = await fetch(`${API_BASE_URL}/projects`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name }),
      });
      const project = await res.json();
      if (project?.id) {
        await fetch(`${API_BASE_URL}/projects/${project.id}/set-current`, { method: 'POST' });
        setCurrentProject(project);
        addLog(`Project "${name}" created`, 'info');
      }
    } catch (err) {
      console.error('Project creation error:', err);
      addLog(`Project creation failed: ${err.message}`, 'error');
    } finally {
      setLoading(false);
    }
  }, [setCurrentProject, addLog]);

  const handleOpenProject = useCallback(async () => {
    if (IS_ELECTRON) {
      setLoading(true);
      try {
        const project = await window.electronAPI.openProject();
        if (project) setCurrentProject(project);
      } catch (err) {
        console.error('Project opening error:', err);
      } finally {
        setLoading(false);
      }
      return;
    }
    // Web/Docker mode: fetch project list and show picker
    try {
      const res = await fetch(`${API_BASE_URL}/projects`);
      const projects = await res.json();
      if (!projects || projects.length === 0) {
        addLog('No projects found. Create a new one first.', 'info');
        return;
      }
      setProjectList(projects);
      setShowProjectList(true);
    } catch (err) {
      console.error('Project list error:', err);
      addLog(`Failed to list projects: ${err.message}`, 'error');
    }
  }, [setCurrentProject, addLog]);

  const handleSelectProject = useCallback(async (project) => {
    setShowProjectList(false);
    setLoading(true);
    try {
      await fetch(`${API_BASE_URL}/projects/${project.id}/set-current`, { method: 'POST' });
      setCurrentProject(project);
      addLog(`Project "${project.name}" opened`, 'info');
    } catch (err) {
      addLog(`Failed to open project: ${err.message}`, 'error');
    } finally {
      setLoading(false);
    }
  }, [setCurrentProject, addLog]);

  const handleDeleteProject = useCallback(async (projectId) => {
    try {
      await fetch(`${API_BASE_URL}/projects/${projectId}`, { method: 'DELETE' });
      setProjectList(prev => prev.filter(p => p.id !== projectId));
      if (currentProject?.id === projectId) setCurrentProject(null);
      addLog('Project deleted', 'info');
    } catch (err) {
      addLog(`Failed to delete project: ${err.message}`, 'error');
    }
  }, [currentProject, setCurrentProject, addLog]);

  const handleSaveProject = useCallback(async () => {
    if (!currentProject) return;
    setLoading(true);
    try {
      if (IS_ELECTRON) {
        const state = { predictionState, annotationState, viewerData, customColors };
        const r = await window.electronAPI.saveProjectState(state);
        if (!r?.success) console.warn('Project save: failed', r?.error);
        else console.log('Project saved');
      } else {
        // Web/Docker mode: save state via API
        const stateKeys = { predictionState, annotationState, viewerData, customColors };
        for (const [key, value] of Object.entries(stateKeys)) {
          await fetch(`${API_BASE_URL}/projects/${currentProject.id}/state`, {
            method: 'PUT',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ key, value }),
          });
        }
        addLog('Project saved', 'info');
      }
    } catch (err) {
      console.error('Project save error:', err);
      addLog(`Save failed: ${err.message}`, 'error');
    } finally {
      setLoading(false);
    }
  }, [currentProject, predictionState, annotationState, viewerData, customColors, addLog]);

  const handleSaveAs = useCallback(async () => {
    if (!currentProject) return;
    if (IS_ELECTRON) {
      setLoading(true);
      try {
        const state = { predictionState, annotationState, viewerData, customColors };
        const result = await window.electronAPI.saveProjectAs(state);
        if (!result) { addLog('Save As: No response from system', 'error'); return; }
        if (result.cancelled) return;
        if (!result.success) { addLog(`Save As failed: ${result.error || 'Unknown error'}`, 'error'); setShowLogs(true); return; }
        if (result.warning) addLog(`Save As completed with warning: ${result.warning}`, 'warn');
        addLog(`Project saved to: ${result.path}`, 'info');
        setCurrentProject(result);
      } catch (err) {
        console.error('Save as error:', err);
        addLog(`Save As error: ${err.message}`, 'error');
        setShowLogs(true);
      } finally {
        setLoading(false);
      }
      return;
    }
    // Web/Docker mode: create a copy with new name
    const newName = prompt('Save as (new name):', currentProject.name + ' (copy)');
    if (!newName) return;
    setLoading(true);
    try {
      const res = await fetch(`${API_BASE_URL}/projects`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name: newName }),
      });
      const newProject = await res.json();
      if (newProject?.id) {
        // Save current state to the new project
        const stateKeys = { predictionState, annotationState, viewerData, customColors };
        for (const [key, value] of Object.entries(stateKeys)) {
          await fetch(`${API_BASE_URL}/projects/${newProject.id}/state`, {
            method: 'PUT',
            headers: { 'Content-Type': 'application/json' },
            body: JSON.stringify({ key, value }),
          });
        }
        await fetch(`${API_BASE_URL}/projects/${newProject.id}/set-current`, { method: 'POST' });
        setCurrentProject(newProject);
        addLog(`Project saved as "${newName}"`, 'info');
      }
    } catch (err) {
      addLog(`Save As failed: ${err.message}`, 'error');
    } finally {
      setLoading(false);
    }
  }, [currentProject, predictionState, annotationState, viewerData, customColors, setCurrentProject, addLog, setShowLogs]);

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
        title="Tolerance Optimization"
        size="large"
      >
        <TabbedToleranceInterface onClose={() => setShowTolerance(false)} />
      </Modal>

      {/* Project list dialog (web/Docker mode) */}
      <Modal
        isOpen={showProjectList}
        onClose={() => setShowProjectList(false)}
        title="Open Project"
        size="small"
      >
        <div style={{ maxHeight: '400px', overflowY: 'auto' }}>
          {projectList.length === 0 ? (
            <p style={{ textAlign: 'center', color: '#888' }}>No projects found</p>
          ) : (
            projectList.map(p => (
              <div key={p.id} style={{
                display: 'flex', alignItems: 'center', justifyContent: 'space-between',
                padding: '10px 12px', borderBottom: '1px solid #eee', cursor: 'pointer',
              }}>
                <div style={{ flex: 1 }} onClick={() => handleSelectProject(p)}>
                  <div style={{ fontWeight: 600 }}>{p.name}</div>
                  <div style={{ fontSize: '0.8em', color: '#888' }}>
                    {new Date(p.updatedAt || p.createdAt).toLocaleDateString()}
                  </div>
                </div>
                <button
                  onClick={(e) => { e.stopPropagation(); if (window.confirm(`Delete "${p.name}"?`)) handleDeleteProject(p.id); }}
                  style={{ background: 'none', border: 'none', color: '#e44', cursor: 'pointer', fontSize: '1.1em', padding: '4px 8px' }}
                  title="Delete project"
                >
                  ✕
                </button>
              </div>
            ))
          )}
        </div>
      </Modal>

      {/* Pending Tasks Notification - shows tasks that survived disconnection */}
      <PendingTasksNotification />
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