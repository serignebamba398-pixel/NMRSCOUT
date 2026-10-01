import { useState, useEffect } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import './ProjectManager.css';

/**
 * Project Manager - Allows creating/opening/closing projects
 * Electron application only
 */
export const ProjectManager = () => {
  const { currentProject, setCurrentProject } = useAppState();
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [backendReady, setBackendReady] = useState(false);

  // Load current project at startup
  useEffect(() => {
    const loadCurrentProject = async () => {
      try {
        const project = await window.electronAPI.getCurrentProject();
        console.log('[ProjectManager] Project loaded from Electron:', project);
        if (project) {
          setCurrentProject(project);
        }
      } catch (err) {
        console.error('Error loading project:', err);
      }
    };

    loadCurrentProject();

    // Listen for backend events
    const handleBackendReady = (data) => {
      console.log('Backend ready:', data);
      setBackendReady(true);
      setError(null);
    };

    const handleBackendError = (data) => {
      console.error('Backend error:', data);
      setError(data.error || 'Backend connection error');
      setBackendReady(false);
    };

    window.electronAPI.onBackendReady(handleBackendReady);
    window.electronAPI.onBackendError(handleBackendError);
  }, [setCurrentProject]);

  const handleCreateProject = async () => {
    setLoading(true);
    setError(null);
    try {
      const project = await window.electronAPI.createProject();
      console.log('[ProjectManager] Project created:', project);
      if (project) {
        setCurrentProject(project);
      }
    } catch (err) {
      setError(err.message || 'Error creating project');
    } finally {
      setLoading(false);
    }
  };

  const handleOpenProject = async () => {
    setLoading(true);
    setError(null);
    try {
      const project = await window.electronAPI.openProject();
      console.log('[ProjectManager] Project opened:', project);
      if (project) {
        setCurrentProject(project);
      }
    } catch (err) {
      setError(err.message || 'Error opening project');
    } finally {
      setLoading(false);
    }
  };

  const handleCloseProject = async () => {
    try {
      await window.electronAPI.closeProject();
      console.log('[ProjectManager] Project closed');
      setCurrentProject(null);
    } catch (err) {
      setError(err.message || 'Error closing project');
    }
  };

  return (
    <div className="project-manager">
      <div className="project-toolbar">
        <h2 className="project-toolbar-title">Project Management</h2>
        <div className="project-toolbar-actions">
          <button
            onClick={handleCreateProject}
            disabled={loading}
            className="btn btn-primary"
          >
            New project
          </button>
          <button
            onClick={handleOpenProject}
            disabled={loading}
            className="btn btn-secondary"
          >
            Open project
          </button>
          {currentProject && (
            <button
              onClick={handleCloseProject}
              disabled={loading}
              className="btn btn-danger"
            >
              ✕ Close project
            </button>
          )}
        </div>
      </div>

      {error && (
        <div className="error-box">
          <strong>Error:</strong> {error}
        </div>
      )}

      {loading && (
        <div className="info-box">
          Loading...
        </div>
      )}

      {currentProject && (
        <div className="project-info">
          <div className="project-info-header">
            <h3 className="project-info-name">{currentProject.name}</h3>
            <span className={`backend-status ${backendReady ? 'ready' : 'loading'}`}>
              {backendReady ? '✓ Backend ready' : 'Backend starting...'}
            </span>
          </div>
          <p className="project-info-path">{currentProject.path}</p>
          {currentProject.createdAt && (
            <p className="project-info-date">
              Created on {new Date(currentProject.createdAt).toLocaleString()}
            </p>
          )}
        </div>
      )}
    </div>
  );
};
