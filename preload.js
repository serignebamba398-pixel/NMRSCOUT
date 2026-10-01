const { contextBridge, ipcRenderer } = require('electron');

// Debug: preload loaded
console.log('[Preload] preload.js loaded');

contextBridge.exposeInMainWorld('electronAPI', {
  // Project management.
  createProject: async () => { console.log('[Preload] createProject invoked'); return await ipcRenderer.invoke('create-project'); },
  openProject: async () => { console.log('[Preload] openProject invoked'); return await ipcRenderer.invoke('open-project'); },
  saveProjectState: async (state) => { console.log('[Preload] saveProjectState invoked'); return await ipcRenderer.invoke('save-project-state', state); },
  saveProjectAs: async (state) => { console.log('[Preload] saveProjectAs invoked'); return await ipcRenderer.invoke('save-project-as', state); },
  getCurrentProject: async () => { console.log('[Preload] getCurrentProject invoked'); return await ipcRenderer.invoke('get-current-project'); },
  getLastProject: async () => { console.log('[Preload] getLastProject invoked'); return await ipcRenderer.invoke('get-last-project'); },
  closeProject: async () => { console.log('[Preload] closeProject invoked'); return await ipcRenderer.invoke('close-project'); },
  
  // Directory selection.
  selectDirectory: async () => { console.log('[Preload] selectDirectory invoked'); return await ipcRenderer.invoke('select-directory'); },
  
  // System paths.
  getUserDataPath: async () => { console.log('[Preload] getUserDataPath invoked'); return await ipcRenderer.invoke('get-user-data-path'); },
  getDocumentsPath: async () => { console.log('[Preload] getDocumentsPath invoked'); return await ipcRenderer.invoke('get-documents-path'); },

  // Export logs (zip)
  exportLogs: async (rendererLogs = []) => {
    console.log('[Preload] exportLogs invoked');
    return await ipcRenderer.invoke('export-logs', { rendererLogs });
  },
  
  // Events.
  onBackendReady: (callback) => ipcRenderer.on('backend-ready', (event, data) => callback(data)),
  onBackendError: (callback) => ipcRenderer.on('backend-error', (event, data) => callback(data)),
  
  // Expose the Electron runtime flag.
  isElectron: true,
});

