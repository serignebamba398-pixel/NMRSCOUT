import { isMethodAvailable, restoreAvailableSelections } from '../constants/releaseFeatures';
import React, { createContext, useContext, useState, useEffect, useCallback, useRef, useMemo, useSyncExternalStore } from 'react';
import { getDefaultParameters, API_BASE_URL } from '../api';
import { useAnnotation, usePrediction, useToleranceOptimization } from '../hooks';
import {
  DEFAULT_COLORS,
  sessionStorage as ss,
  DEFAULT_PREDICTION_OPTIONS,
  DEFAULT_TOLERANCE_RANGES,
  DEFAULT_BG_OPTIONS
} from '../constants';
import { normalizeStaticUrl } from '../utils';
import { logStore } from '../stores/logStore';

const IS_ELECTRON = window.electronAPI?.isElectron === true;

const AppStateContext = createContext();

// Initial states for reset when project changes
const INITIAL_PREDICTION_STATE = {
  mode: "SMILES",
  smiles: "CCO",
  molId: "ethanol",
  csvLocalFile: null,
  csvUploadedName: null,
  uploading: false,
  options: DEFAULT_PREDICTION_OPTIONS,
  selectedPredMethods: ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable),
  artifacts: [],
  predictionTaskId: null,
  showCombinedViewer: false,
};

// Initial state for global progress (central bar)
const INITIAL_GLOBAL_PROGRESS = {
  visible: false,
  type: null, // 'prediction' | 'annotation' | 'optimization'
  progress: 0,
  message: '',
  title: '',
};

const getInitialAnnotationState = () => ({
  predictionTasks: [],
  selectedTask: null,
  activeType: "1H",
  matrixFiles: {},
  uploadedMatrixNames: {},
  uploadingMatrix: {},
  parameters: {
    "1H": getDefaultParameters("1H"),
    "13C": getDefaultParameters("13C"),
    COSY: getDefaultParameters("COSY"),
    TOCSY: getDefaultParameters("TOCSY"),
    HSQC: getDefaultParameters("HSQC"),
  },
  lastResultByType: {},
  viewerTaskId: null,
  allRunning: false,
  showCombined: false,
  globalProgress: 0,
  globalProgressMessage: "",
  readyToView: false,
  selectedMethods: ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable),
  bgPdataDirs: { "1H": "", "13C": "", "COSY": "", "TOCSY": "", "HSQC": "" },
  bgOptions: DEFAULT_BG_OPTIONS,
  bgStatus: {},
  toleranceRanges: DEFAULT_TOLERANCE_RANGES,
  toleranceResults: [],
  // Tolerance optimization state (preserved when modal closes)
  toleranceOptimProgress: 0,
  toleranceOptimRunning: false,
  toleranceOptimCurrentTest: "",
  toleranceOptimSelectedType: "1H",
});

/**
 * Global application context provider
 * Manages all shared states across different pages
 */
export const AppStateProvider = ({ children }) => {
  // ========== PREDICTIONS STATE ==========
  const [predictionState, setPredictionState] = useState(INITIAL_PREDICTION_STATE);

  // ========== ANNOTATIONS STATE ==========
  const [annotationState, setAnnotationState] = useState(getInitialAnnotationState());

  // ========== VIEWER STATE ==========
  const [viewerData, setViewerData] = useState(null);

  // ========== PROGRESSION GLOBALE (barre centrale) ==========
  const [globalProgress, setGlobalProgress] = useState(INITIAL_GLOBAL_PROGRESS);

  // ========== LOGS ==========
  // Logs are stored outside of React state so frequent log appends do not
  // re-render the whole application tree.
  const [showLogs, setShowLogs] = useState(false);

  // ========== PERSONNALISATION ==========
  const [customColors, setCustomColors] = useState(DEFAULT_COLORS);

  // ========== VERROUILLAGE UI ==========
  const [uiLock, setUiLock] = useState({ lockedBy: null, message: "" });

  // ========== PROJET ACTUEL ==========
  const [currentProject, setCurrentProject] = useState(null);
  const previousProjectRef = useRef(null);

  const normalizeViewerData = useCallback((raw) => {
    if (!raw || typeof raw !== 'object') return raw;
    const fixed = {
      ...raw,
      csvs: Object.fromEntries(
        Object.entries(raw.csvs || {}).map(([k, v]) => [k, normalizeStaticUrl(v)])
      ),
      downloads: Object.fromEntries(
        Object.entries(raw.downloads || {}).map(([k, v]) => [k, normalizeStaticUrl(v)])
      ),
    };

    if (fixed.annotationResults && typeof fixed.annotationResults === 'object') {
      const next = { ...fixed.annotationResults };
      for (const [k, v] of Object.entries(next)) {
        if (!v || typeof v !== 'object') continue;
        next[k] = {
          ...v,
          downloadUrl: normalizeStaticUrl(v.downloadUrl),
          pngUrl: normalizeStaticUrl(v.pngUrl),
          jsonUrl: normalizeStaticUrl(v.jsonUrl),
        };
      }
      fixed.annotationResults = next;
    }

    return fixed;
  }, []);

  const getProjectKey = useCallback((p) => {
    // Electron projects use `path` as stable identity; backend projects may have `id`.
    return p?.id || p?.path || null;
  }, []);

  const hydrateFromProjectState = useCallback((projectState) => {
    if (!projectState || typeof projectState !== 'object') return;

    if (projectState.predictionState && typeof projectState.predictionState === 'object') {
      setPredictionState((prev) => {
        const base = INITIAL_PREDICTION_STATE;
        const saved = restoreAvailableSelections(projectState.predictionState);
        return {
          ...base,
          ...saved,
          // avoid non-serializable fields
          csvLocalFile: null,
          uploading: false,
          options: {
            ...base.options,
            ...(saved.options || {}),
          },
        };
      });
    }

    if (projectState.annotationState && typeof projectState.annotationState === 'object') {
      setAnnotationState(() => {
        const base = getInitialAnnotationState();
        const saved = restoreAvailableSelections(projectState.annotationState);

        // Deep merge toleranceRanges to preserve structure
        const mergedToleranceRanges = { ...base.toleranceRanges };
        if (saved.toleranceRanges && typeof saved.toleranceRanges === 'object') {
          for (const type of Object.keys(base.toleranceRanges)) {
            if (saved.toleranceRanges[type] && typeof saved.toleranceRanges[type] === 'object') {
              mergedToleranceRanges[type] = {
                ...base.toleranceRanges[type],
                ...saved.toleranceRanges[type],
              };
              // Deep merge each dimension (min/max/step)
              for (const dim of Object.keys(base.toleranceRanges[type])) {
                if (saved.toleranceRanges[type][dim] && typeof saved.toleranceRanges[type][dim] === 'object') {
                  mergedToleranceRanges[type][dim] = {
                    ...base.toleranceRanges[type][dim],
                    ...saved.toleranceRanges[type][dim],
                  };
                }
              }
            }
          }
        }

        // Deep merge parameters by type to preserve new default params
        const mergedParameters = { ...base.parameters };
        if (saved.parameters && typeof saved.parameters === 'object') {
          for (const type of Object.keys(base.parameters)) {
            if (saved.parameters[type] && typeof saved.parameters[type] === 'object') {
              mergedParameters[type] = {
                ...base.parameters[type],
                ...saved.parameters[type],
              };
            }
          }
        }

        return {
          ...base,
          ...saved,
          // avoid non-serializable fields
          matrixFiles: {},
          uploadingMatrix: {},
          allRunning: false,
          globalProgress: 0,
          globalProgressMessage: "",
          parameters: mergedParameters,
          bgPdataDirs: {
            ...base.bgPdataDirs,
            ...(saved.bgPdataDirs || {}),
          },
          bgOptions: {
            ...base.bgOptions,
            ...(saved.bgOptions || {}),
          },
          toleranceRanges: mergedToleranceRanges,
        };
      });
    }

    if (projectState.viewerData && typeof projectState.viewerData === 'object') {
      setViewerData(normalizeViewerData(projectState.viewerData));
    }

    if (projectState.customColors && typeof projectState.customColors === 'object') {
      setCustomColors(projectState.customColors);
    }
  }, [normalizeViewerData]);

  // ========== RESET ON PROJECT CHANGE ==========
  useEffect(() => {
    const newProjectKey = getProjectKey(currentProject);
    const prevProjectKey = previousProjectRef.current;

    // First mount: hydrate if a project is already loaded
    if (prevProjectKey === null) {
      previousProjectRef.current = newProjectKey;
      if (currentProject?.state) {
        hydrateFromProjectState(currentProject.state);
      }
      return;
    }

    if (newProjectKey !== prevProjectKey) {
      console.log('[AppState] Project change detected:', prevProjectKey, '->', newProjectKey);
      
      // Reset all project-related states
      setPredictionState(INITIAL_PREDICTION_STATE);
      setAnnotationState(getInitialAnnotationState());
      setViewerData(null);
      logStore.clear();
      setUiLock({ lockedBy: null, message: "" });
      
      // Clean sessionStorage
      ss.remove("viewerData");
      
      console.log('[AppState] States reset for new project');

      // Then hydrate from project state (if available)
      if (currentProject?.state) {
        hydrateFromProjectState(currentProject.state);
      }
    }

    previousProjectRef.current = newProjectKey;
  }, [currentProject, getProjectKey, hydrateFromProjectState]);

  // Log state
  const addLog = useCallback((message, type = "info") => {
    const timestamp = new Date().toLocaleTimeString();
    logStore.addLog({ message, type, timestamp });
  }, []);

  const clearLogs = useCallback(() => logStore.clear(), []);

  // ========== PERSISTENT ANNOTATION HOOK ==========
  // This hook lives at context level so annotations
  // continue even if modal is closed
  const annotation = useAnnotation({
    annotationState,
    setAnnotationState,
    setViewerData,
    addLog,
    setShowLogs,
    setUiLock,
  });

  // ========== PERSISTENT PREDICTION HOOK ==========
  // This hook lives at context level so predictions
  // continue even if modal is closed
  const prediction = usePrediction({
    predictionState,
    setPredictionState,
    addLog,
    setShowLogs,
  });

  // ========== PERSISTENT TOLERANCE OPTIMIZATION HOOK ==========
  // This hook lives at context level so optimization
  // continues even if modal is closed
  const toleranceOptimization = useToleranceOptimization({
    annotationState,
    setAnnotationState,
    addLog,
    setShowLogs,
    uiLock,
    setUiLock,
  });

  // ========== PERSISTENCE (sessionStorage) ==========
  useEffect(() => {
    // Restore viewer data from session storage in development or browser mode.
    const savedViewer = ss.get("viewerData", null);
    if (savedViewer && typeof savedViewer === "object") {
      const fixed = {
        ...savedViewer,
        csvs: Object.fromEntries(
          Object.entries(savedViewer.csvs || {}).map(([k, v]) => [
            k,
            normalizeStaticUrl(v),
          ])
        ),
        downloads: Object.fromEntries(
          Object.entries(savedViewer.downloads || {}).map(([k, v]) => [
            k,
            normalizeStaticUrl(v),
          ])
        ),
      };
      setViewerData(fixed);
    }

    // Restore custom colors from session storage.
    const savedColors = ss.get("customColors", null);
    if (savedColors && typeof savedColors === "object") {
      setCustomColors(savedColors);
    }
  }, []);

  // Persist viewerData changes, removing the saved value when null.
  useEffect(() => {
    if (viewerData) {
      ss.set("viewerData", viewerData);
    } else {
      ss.remove("viewerData");
    }
  }, [viewerData]);

  // Persist custom color changes.
  useEffect(() => {
    ss.set("customColors", customColors);
  }, [customColors]);

  // ========== WEB/DOCKER: load last project on startup ==========
  useEffect(() => {
    if (IS_ELECTRON) return;
    const loadCurrentProject = async () => {
      try {
        const res = await fetch(`${API_BASE_URL}/projects/current`);
        const project = await res.json();
        if (project?.id) {
          setCurrentProject(project);
        }
      } catch {
        // No current project or backend not ready yet
      }
    };
    loadCurrentProject();
  }, []); // eslint-disable-line react-hooks/exhaustive-deps

  // ========== PERSISTENCE PROJET (Electron) ==========
  // Do not overwrite `project.nmrproj` unless the user explicitly
  // clicks Save / Save As. Annotation/prediction updates must not auto-persist
  // to disk because it can overwrite a previous project state without consent.
  //
  // We still persist lightweight state to sessionStorage above (viewer/customColors)
  // for the current session, but disk writes are manual-only.

  // Normalize prediction settings.
  useEffect(() => {
    setPredictionState((prev) => {
      const o = prev?.options || {};
      const ex = o.extractors || {};
      const legacyC13 = o.c13 || o.C13 || {};
      const legacyHSQC = o.hsqc || o.HSQC || {};
      const newC13 = ex.c13 || {};
      const newHSQC = ex.hsqc || {};

      const carbon_mode =
        newC13.carbon_mode ??
        legacyC13.carbon_mode ??
        legacyC13.carbonMode ??
        "all";

      const fold = newHSQC.fold ?? legacyHSQC.fold ?? false;
      const fold_low = Number.isFinite(+newHSQC.fold_low)
        ? +newHSQC.fold_low
        : Number.isFinite(+legacyHSQC.fold_low)
        ? +legacyHSQC.fold_low
        : 4.8618;
      const fold_high = Number.isFinite(+newHSQC.fold_high)
        ? +newHSQC.fold_high
        : Number.isFinite(+legacyHSQC.fold_high)
        ? +legacyHSQC.fold_high
        : 69.8395;

      return {
        ...prev,
        options: {
          ...o,
          extractors: {
            ...ex,
            c13: { carbon_mode },
            hsqc: { fold, fold_low, fold_high },
          },
        },
      };
    });
  }, []);

  // Keep axis limits synchronized even when the Bruker path is unchanged.
  useEffect(() => {
    const { bgPdataDirs, bgOptions } = annotationState;
    setViewerData(prev => {
      if (!prev) return prev;
      const pdataDirs = { ...prev.pdataDirs, ...bgPdataDirs };
      const options = { ...prev.bgOptions, ...bgOptions };
      if (JSON.stringify(pdataDirs) === JSON.stringify(prev.pdataDirs || {}) &&
          JSON.stringify(options) === JSON.stringify(prev.bgOptions || {})) return prev;
      return { ...prev, pdataDirs, bgOptions: options };
    });
  }, [annotationState.bgPdataDirs, annotationState.bgOptions]);

  // ========== CONTEXT VALUES ==========
  const value = useMemo(() => ({
    predictionState,
    setPredictionState,
    annotationState,
    setAnnotationState,
    viewerData,
    setViewerData,
    // Global progress (central bar)
    globalProgress,
    setGlobalProgress,
    addLog,
    clearLogs,
    showLogs,
    setShowLogs,
    uiLock,
    setUiLock,
    customColors,
    setCustomColors,
    currentProject,
    setCurrentProject,
    // Persistent annotation functions (don't stop when modal closes)
    annotation,
    // Persistent prediction functions
    prediction,
    // Persistent tolerance optimization functions
    toleranceOptimization,
  }), [
    predictionState,
    setPredictionState,
    annotationState,
    setAnnotationState,
    viewerData,
    setViewerData,
    globalProgress,
    setGlobalProgress,
    addLog,
    clearLogs,
    showLogs,
    setShowLogs,
    uiLock,
    setUiLock,
    customColors,
    setCustomColors,
    currentProject,
    setCurrentProject,
    annotation,
    prediction,
    toleranceOptimization,
  ]);

  return (
    <AppStateContext.Provider value={value}>
      {children}
    </AppStateContext.Provider>
  );
};

/**
 * Hook to access global context
 * @throws {Error} If used outside AppStateProvider
 */
export const useAppState = () => {
  const context = useContext(AppStateContext);
  if (!context) {
    throw new Error("useAppState must be used within AppStateProvider");
  }
  return context;
};

/**
 * Hook to access logs context
 * @throws {Error} If used outside AppStateProvider
 */
export const useLogs = () => {
  const { showLogs, setShowLogs, addLog, clearLogs } = useAppState();
  const logs = useSyncExternalStore(logStore.subscribe, logStore.getSnapshot);
  return { logs, showLogs, setShowLogs, addLog, clearLogs };
};