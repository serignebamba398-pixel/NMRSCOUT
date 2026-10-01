import { isMethodAvailable } from '../../constants/releaseFeatures';
import React from 'react';
import { useNavigate } from 'react-router-dom';
import { ProjectManager } from '../project';
import { useAppState } from '../../contexts/AppStateContext';
import './HomePage.css';

/**
 * Home page - Displays the project manager and shortcuts
 */
export const HomePage = () => {
  const navigate = useNavigate();
  const { currentProject } = useAppState();

  // Log project state for debugging.
  console.log('[HomePage] currentProject:', currentProject);

  const shortcuts = [
    {
      title: 'NMR Prediction',
      description: `Predict NMR spectra (${['1H', '13C', 'COSY', 'TOCSY', 'HSQC'].filter(isMethodAvailable).join(', ')}) from SMILES or CSV`,
      icon: '¹H',
      path: '/prediction',
      color: 'blue',
      disabled: !currentProject,
    },
    {
      title: 'Integrated Annotation',
      description: 'Annotate NMR spectra with predicted data',
      icon: 'ID',
      path: '/annotation',
      color: 'green',
      disabled: !currentProject,
    },
    {
      title: 'Tolerance Optimization',
      description: 'Optimize tolerance parameters for annotation',
      icon: '±',
      path: '/optimization',
      color: 'purple',
      disabled: !currentProject,
    },
    {
      title: 'Results',
      description: 'View results and annotated spectra',
      icon: '≡',
      path: '/results',
      color: 'orange',
      disabled: !currentProject,
    },
  ];

  return (
    <div className="home-page">
      {/* Header */}
      <div className="home-header">
        <h1 className="home-title">NMR Annotation Desktop</h1>
        <p className="home-subtitle">
          NMR spectrum annotation and prediction platform
        </p>
      </div>

      {/* Project manager */}
      <ProjectManager />

      {/* Message if no project */}
      {!currentProject && (
        <div className="warning-box">
          <span className="warning-icon" aria-hidden="true">!</span>
          <div>
            <h3 className="warning-title">No project open</h3>
            <p className="warning-text">
              You must create or open a project before using prediction and annotation features.
              Files (uploads, outputs, data) will be stored in your project folder.
            </p>
          </div>
        </div>
      )}

      {/* Feature shortcuts */}
      <div className="shortcuts-grid">
        {shortcuts.map((shortcut) => (
          <button
            key={shortcut.path}
            onClick={() => !shortcut.disabled && navigate(shortcut.path)}
            disabled={shortcut.disabled}
            className={`shortcut-card ${shortcut.disabled ? 'disabled' : ''} ${shortcut.color}`}
          >
            <span className="shortcut-icon">{shortcut.icon}</span>
            <div className="shortcut-content">
              <h3 className="shortcut-title">{shortcut.title}</h3>
              <p className="shortcut-description">{shortcut.description}</p>
              {shortcut.disabled && (
                <p className="shortcut-warning">
                  Open a project to access this feature
                </p>
              )}
            </div>
          </button>
        ))}
      </div>

      {/* Current project information */}
      {currentProject && (
        <div className="project-info-box">
          <h3 className="project-info-title">
            Active project: {currentProject.name}
          </h3>
          <div className="project-info-content">
            <p className="project-path">{currentProject.path}</p>
            <p className="project-info-label">
              Files will be saved in:
            </p>
            <ul className="project-folders">
              <li><code>{currentProject.path}/uploads</code> - Your uploaded files</li>
              <li><code>{currentProject.path}/outputs</code> - Generated results</li>
              <li><code>{currentProject.path}/data</code> - NMR data</li>
            </ul>
          </div>
        </div>
      )}
    </div>
  );
};
