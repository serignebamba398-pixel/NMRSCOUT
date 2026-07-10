import React, { useState, useEffect, useCallback } from 'react';
import { useAppState } from '../../contexts/AppStateContext';
import { API_BASE_URL } from '../../api';

/**
 * PendingTasksNotification
 * 
 * Affiche une notification si des tâches sont en attente après
 * une déconnexion (fermeture de l'ordinateur, refresh du navigateur, etc.)
 */
export const PendingTasksNotification = () => {
  const [pendingTasks, setPendingTasks] = useState([]);
  const [showNotification, setShowNotification] = useState(false);
  const [loading, setLoading] = useState(true);
  const { addLog } = useAppState();

  // Vérifier les tâches en attente au chargement
  const checkPendingTasks = useCallback(async () => {
    try {
      const response = await fetch(`${API_BASE_URL}/tasks/pending`);
      if (response.ok) {
        const data = await response.json();
        if (data.tasks && data.tasks.length > 0) {
          setPendingTasks(data.tasks);
          setShowNotification(true);
          addLog(`${data.tasks.length} tâche(s) en attente détectée(s)`, 'info');
        }
      }
    } catch (error) {
      console.warn('Could not check pending tasks:', error.message);
    } finally {
      setLoading(false);
    }
  }, [addLog]);

  useEffect(() => {
    // Vérifier au chargement initial
    checkPendingTasks();
    
    // Vérifier périodiquement (toutes les 30 secondes)
    const interval = setInterval(() => {
      checkPendingTasks();
    }, 30000);

    return () => clearInterval(interval);
  }, [checkPendingTasks]);

  // Ignorer une tâche (ne pas la reprendre)
  const handleDismissTask = async (taskId) => {
    try {
      await fetch(`${API_BASE_URL}/tasks/${taskId}/dismiss`, {
        method: 'POST',
      });
      setPendingTasks(prev => prev.filter(t => t.taskId !== taskId));
      addLog(`Tâche ${taskId.slice(0, 8)}... ignorée`, 'info');
      
      if (pendingTasks.length <= 1) {
        setShowNotification(false);
      }
    } catch (error) {
      console.error('Failed to dismiss task:', error);
    }
  };

  // Ignorer toutes les tâches
  const handleDismissAll = async () => {
    for (const task of pendingTasks) {
      try {
        await fetch(`${API_BASE_URL}/tasks/${task.taskId}/dismiss`, {
          method: 'POST',
        });
      } catch (error) {
        console.warn(`Failed to dismiss task ${task.taskId}:`, error);
      }
    }
    setPendingTasks([]);
    setShowNotification(false);
    addLog('Toutes les tâches en attente ont été ignorées', 'info');
  };

  // Consulter les résultats d'une tâche
  const handleViewTask = async (task) => {
    try {
      const response = await fetch(`${API_BASE_URL}/tasks/${task.taskId}/checkpoint`);
      if (response.ok) {
        const data = await response.json();
        if (data.checkpoint?.status === 'completed' && data.checkpoint?.outputs) {
          addLog(`Tâche ${task.taskId.slice(0, 8)} terminée - fichiers disponibles`, 'success');
          // TODO: Charger les fichiers dans le viewer
        } else {
          addLog(`Tâche ${task.taskId.slice(0, 8)}: progression ${data.checkpoint?.progress?.percent || 0}%`, 'info');
        }
      }
    } catch (error) {
      console.error('Failed to get task checkpoint:', error);
    }
  };

  if (loading || !showNotification || pendingTasks.length === 0) {
    return null;
  }

  return (
    <div className="pending-tasks-notification">
      <div className="notification-header">
        <span className="notification-icon">⏳</span>
        <span className="notification-title">
          {pendingTasks.length === 1 
            ? '1 tâche en attente' 
            : `${pendingTasks.length} tâches en attente`}
        </span>
        <button 
          className="notification-close"
          onClick={() => setShowNotification(false)}
          title="Masquer (ne supprime pas les tâches)"
        >
          ×
        </button>
      </div>
      
      <div className="notification-body">
        <p className="notification-description">
          Ces tâches étaient en cours lors de la dernière session et peuvent avoir terminé.
        </p>
        
        <ul className="pending-tasks-list">
          {pendingTasks.map(task => (
            <li key={task.taskId} className="pending-task-item">
              <div className="task-info">
                <span className="task-type">
                  {task.type === 'prediction' ? '🧪' : task.type === 'annotation' ? '📋' : '⚙️'}
                  {task.type}
                </span>
                <span className="task-progress">
                  {task.progress ? `${task.progress.percent}%` : '?'}
                </span>
                <span className="task-time">
                  {task.interruptedAt 
                    ? new Date(task.interruptedAt).toLocaleString() 
                    : ''}
                </span>
              </div>
              <div className="task-actions">
                <button 
                  className="btn-small btn-primary"
                  onClick={() => handleViewTask(task)}
                  title="Voir le statut/résultats"
                >
                  👁️ Voir
                </button>
                <button 
                  className="btn-small btn-secondary"
                  onClick={() => handleDismissTask(task.taskId)}
                  title="Ignorer cette tâche"
                >
                  ✕ Ignorer
                </button>
              </div>
            </li>
          ))}
        </ul>
      </div>
      
      <div className="notification-footer">
        <button 
          className="btn-text"
          onClick={handleDismissAll}
        >
          Ignorer toutes les tâches
        </button>
      </div>

      <style jsx="true">{`
        .pending-tasks-notification {
          position: fixed;
          bottom: 80px;
          right: 20px;
          width: 360px;
          background: #1e2530;
          border: 1px solid #3a4558;
          border-radius: 8px;
          box-shadow: 0 4px 20px rgba(0, 0, 0, 0.3);
          z-index: 1000;
          overflow: hidden;
          animation: slideIn 0.3s ease;
        }

        @keyframes slideIn {
          from {
            opacity: 0;
            transform: translateY(20px);
          }
          to {
            opacity: 1;
            transform: translateY(0);
          }
        }

        .notification-header {
          display: flex;
          align-items: center;
          padding: 12px 16px;
          background: #2a3444;
          border-bottom: 1px solid #3a4558;
        }

        .notification-icon {
          font-size: 18px;
          margin-right: 10px;
        }

        .notification-title {
          flex: 1;
          font-weight: 600;
          color: #ffc107;
        }

        .notification-close {
          background: none;
          border: none;
          color: #8899aa;
          font-size: 20px;
          cursor: pointer;
          padding: 0 4px;
        }

        .notification-close:hover {
          color: #fff;
        }

        .notification-body {
          padding: 12px 16px;
        }

        .notification-description {
          font-size: 12px;
          color: #8899aa;
          margin: 0 0 12px 0;
        }

        .pending-tasks-list {
          list-style: none;
          padding: 0;
          margin: 0;
          max-height: 200px;
          overflow-y: auto;
        }

        .pending-task-item {
          display: flex;
          justify-content: space-between;
          align-items: center;
          padding: 8px;
          background: #2a3444;
          border-radius: 4px;
          margin-bottom: 8px;
        }

        .pending-task-item:last-child {
          margin-bottom: 0;
        }

        .task-info {
          display: flex;
          flex-direction: column;
          gap: 2px;
        }

        .task-type {
          font-size: 13px;
          font-weight: 500;
          text-transform: capitalize;
        }

        .task-progress {
          font-size: 11px;
          color: #4da6ff;
        }

        .task-time {
          font-size: 10px;
          color: #6b7a8c;
        }

        .task-actions {
          display: flex;
          gap: 6px;
        }

        .btn-small {
          padding: 4px 8px;
          font-size: 11px;
          border: none;
          border-radius: 4px;
          cursor: pointer;
          display: flex;
          align-items: center;
          gap: 4px;
        }

        .btn-primary {
          background: #4da6ff;
          color: #fff;
        }

        .btn-primary:hover {
          background: #3d8ad9;
        }

        .btn-secondary {
          background: #3a4558;
          color: #ccc;
        }

        .btn-secondary:hover {
          background: #4a5568;
        }

        .notification-footer {
          padding: 8px 16px;
          border-top: 1px solid #3a4558;
          text-align: center;
        }

        .btn-text {
          background: none;
          border: none;
          color: #8899aa;
          font-size: 12px;
          cursor: pointer;
        }

        .btn-text:hover {
          color: #ff6b6b;
          text-decoration: underline;
        }
      `}</style>
    </div>
  );
};

export default PendingTasksNotification;
