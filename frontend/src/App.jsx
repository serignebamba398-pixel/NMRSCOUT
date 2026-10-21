import React, { useState, useCallback, useEffect, useRef } from 'react';
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
  // Synchronous guard for rapid clicks before React exposes running=true.
  const taskLaunchLockRef = useRef(null);

  // Cancellation can close the progress streams before their promises settle.
  // Release the click guard as soon as both task states report idle.
  useEffect(() => {
    if (!prediction?.running && !annotationState?.allRunning) {
      taskLaunchLockRef.current = null;
    }
  }, [prediction?.running, annotationState?.allRunning]);

  // Modal states
  const [showPrediction, setShowPrediction] = useState(false);
  const [showAnnotation, setShowAnnotation] = useState(false);
  const [showTolerance, setShowTolerance] = useState(false);
  // Web mode: project list dialog
  const [showProjectList, setShowProjectList] = useState(false);
  const [projectList, setProjectList] = useState([]);
  // Web mode: new project name dialog
  const [showNewProjectDialog, setShowNewProjectDialog] = useState(false);
  const [newProjectName, setNewProjectName] = useState('');
  // Web mode: save-as dialog
  const [showSaveAsDialog, setShowSaveAsDialog] = useState(false);
  const [saveAsName, setSaveAsName] = useState('');

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
    // Web/Docker mode: show name dialog
    setNewProjectName('');
    setShowNewProjectDialog(true);
  }, [setCurrentProject]);

  const handleCreateProjectConfirm = useCallback(async () => {
    const name = newProjectName.trim();
    if (!name) return;
    setShowNewProjectDialog(false);
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
  }, [newProjectName, setCurrentProject, addLog]);

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
    // Web/Docker mode: show save-as dialog
    setSaveAsName(currentProject.name + ' (copy)');
    setShowSaveAsDialog(true);
  }, [currentProject, predictionState, annotationState, viewerData, customColors, setCurrentProject, addLog, setShowLogs]);

  const handleSaveAsConfirm = useCallback(async () => {
    const newName = saveAsName.trim();
    if (!newName || !currentProject) return;
    setShowSaveAsDialog(false);
    setLoading(true);
    try {
      const res = await fetch(`${API_BASE_URL}/projects`, {
        method: 'POST',
        headers: { 'Content-Type': 'application/json' },
        body: JSON.stringify({ name: newName }),
      });
      const newProject = await res.json();
      if (newProject?.id) {
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
  }, [saveAsName, currentProject, predictionState, annotationState, viewerData, customColors, setCurrentProject, addLog]);

  // ========== HANDLERS COMPUTE ==========
  const handleComputePrediction = useCallback(() => {
    if (prediction?.running) {
      window.dispatchEvent(new Event('nmrscout:show-progress'));
      return;
    }

    if (annotationState?.allRunning || taskLaunchLockRef.current) return;

    if (!prediction?.runPrediction) {
      addLog("Prediction not available. Configure parameters first.", "error");
      return;
    }
    taskLaunchLockRef.current = 'prediction';
    setShowLogs(true);
    Promise.resolve(prediction.runPrediction()).finally(() => {
      if (taskLaunchLockRef.current === 'prediction') taskLaunchLockRef.current = null;
    });
  }, [prediction, annotationState?.allRunning, addLog, setShowLogs]);

  const handleComputeAnnotation = useCallback(() => {
    if (annotationState?.allRunning) {
      window.dispatchEvent(new Event('nmrscout:show-progress'));
      return;
    }

    if (prediction?.running || taskLaunchLockRef.current) return;

    if (!annotation?.runSelectedIntegrated) {
      addLog("Annotation not available. Configure parameters first.", "error");
      return;
    }
    taskLaunchLockRef.current = 'annotation';
    setShowLogs(true);
    Promise.resolve(annotation.runSelectedIntegrated()).finally(() => {
      if (taskLaunchLockRef.current === 'annotation') taskLaunchLockRef.current = null;
    });
  }, [annotation, annotationState?.allRunning, prediction?.running, addLog, setShowLogs]);

  return (
    <div className="app sirius-layout">
      {/* Global progress overlay */}
      <GlobalProgressOverlay />

      {/* Top toolbar */}
      <Toolbar
        onNewProject={handleNewProject}
        onOpenProject={handleOpenProject}
        onSaveProject={handleSaveProject}
        onSaveAs={handleSaveAs}
        onPrediction={() => setShowPrediction(true)}
        onAnnotation={() => setShowAnnotation(true)}
        onComputePrediction={handleComputePrediction}
        onComputeAnnotation={handleComputeAnnotation}
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
        size="fullscreen"
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

      {/* New project dialog (web/Docker mode) */}
      <Modal
        isOpen={showNewProjectDialog}
        onClose={() => setShowNewProjectDialog(false)}
        title="New Project"
        size="small"
      >
        <div style={{ padding: '16px' }}>
          <label style={{ display: 'block', marginBottom: 8, fontWeight: 600 }}>Project name</label>
          <input
            type="text"
            value={newProjectName}
            onChange={(e) => setNewProjectName(e.target.value)}
            onKeyDown={(e) => { if (e.key === 'Enter') handleCreateProjectConfirm(); }}
            placeholder="My NMR Project"
            autoFocus
            style={{ width: '100%', padding: '8px 12px', fontSize: '14px', border: '1px solid #ccc', borderRadius: 4, boxSizing: 'border-box' }}
          />
          <div style={{ display: 'flex', justifyContent: 'flex-end', gap: 8, marginTop: 16 }}>
            <button onClick={() => setShowNewProjectDialog(false)} style={{ padding: '8px 16px', border: '1px solid #ccc', borderRadius: 4, background: '#fff', cursor: 'pointer' }}>Cancel</button>
            <button onClick={handleCreateProjectConfirm} disabled={!newProjectName.trim()} style={{ padding: '8px 16px', border: 'none', borderRadius: 4, background: '#2563eb', color: '#fff', cursor: 'pointer', opacity: newProjectName.trim() ? 1 : 0.5 }}>Create</button>
          </div>
        </div>
      </Modal>

      {/* Save As dialog (web/Docker mode) */}
      <Modal
        isOpen={showSaveAsDialog}
        onClose={() => setShowSaveAsDialog(false)}
        title="Save As"
        size="small"
      >
        <div style={{ padding: '16px' }}>
          <label style={{ display: 'block', marginBottom: 8, fontWeight: 600 }}>New project name</label>
          <input
            type="text"
            value={saveAsName}
            onChange={(e) => setSaveAsName(e.target.value)}
            onKeyDown={(e) => { if (e.key === 'Enter') handleSaveAsConfirm(); }}
            autoFocus
            style={{ width: '100%', padding: '8px 12px', fontSize: '14px', border: '1px solid #ccc', borderRadius: 4, boxSizing: 'border-box' }}
          />
          <div style={{ display: 'flex', justifyContent: 'flex-end', gap: 8, marginTop: 16 }}>
            <button onClick={() => setShowSaveAsDialog(false)} style={{ padding: '8px 16px', border: '1px solid #ccc', borderRadius: 4, background: '#fff', cursor: 'pointer' }}>Cancel</button>
            <button onClick={handleSaveAsConfirm} disabled={!saveAsName.trim()} style={{ padding: '8px 16px', border: 'none', borderRadius: 4, background: '#2563eb', color: '#fff', cursor: 'pointer', opacity: saveAsName.trim() ? 1 : 0.5 }}>Save</button>
          </div>
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
 * Application root component
 */
const App = () => (
  <AppStateProvider>
    <AppContent />
  </AppStateProvider>
);

export default App;
