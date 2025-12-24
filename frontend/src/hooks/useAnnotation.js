import { useState, useRef, useCallback } from 'react';
import { 
  getPredictionTasks, 
  uploadFile, 
  listUploads,
  isValidCsvForType,
  cancelTask,
  normalizeServerPath,
  exportAllSpectralBackgroundsWithRetry
} from '../api';
import { createReconnectingSSE } from '../createReconnectingSSE';
import { normalizeStaticUrl } from '../utils';

/**
 * Hook to manage integrated annotation logic
 */
export const useAnnotation = ({
  annotationState,
  setAnnotationState,
  setViewerData,
  addLog,
  setShowLogs,
}) => {
  const [matrixUploads, setMatrixUploads] = useState([]);
  const [loadingMatrixUploads, setLoadingMatrixUploads] = useState(false);
  const [cancelRequested, setCancelRequested] = useState(false);
  
  const runningSSERef = useRef(new Map());

  /**
   * Refreshes the list of uploaded matrices
   */
  const refreshMatrixUploads = useCallback(async () => {
    setLoadingMatrixUploads(true);
    try {
      const list = await listUploads("matrix");
      setMatrixUploads(list);
    } catch (e) {
      addLog(`Error listUploads(matrix): ${e.message}`, "error");
    } finally {
      setLoadingMatrixUploads(false);
    }
  }, [addLog]);

  /**
   * Loads available prediction tasks
   */
  const loadPredictionTasks = useCallback(async () => {
    try {
      const tasks = await getPredictionTasks();
      setAnnotationState((prev) => ({ ...prev, predictionTasks: tasks }));
      
      if (tasks.length > 0 && !annotationState.selectedTask) {
        setAnnotationState((prev) => ({
          ...prev,
          selectedTask: tasks[0].taskId,
        }));
      }
    } catch (e) {
      addLog(`Error loading prediction tasks: ${e.message}`, "error");
    }
  }, [addLog, annotationState.selectedTask, setAnnotationState]);

  /**
   * Uploads a spectral matrix
   */
  const uploadMatrix = async (type, file) => {
    if (!file) return;
    
    setAnnotationState((prev) => ({
      ...prev,
      uploadingMatrix: { ...prev.uploadingMatrix, [type]: true },
    }));

    try {
      const { filename } = await uploadFile(file);
      setAnnotationState((prev) => ({
        ...prev,
        uploadedMatrixNames: { ...prev.uploadedMatrixNames, [type]: filename },
        matrixFiles: { ...prev.matrixFiles, [type]: file },
        uploadingMatrix: { ...prev.uploadingMatrix, [type]: false },
      }));
      addLog(`[${type}] Matrix uploaded: ${filename}`);
    } catch (e) {
      addLog(`[${type}] Matrix upload failed: ${e.message}`, "error");
      setAnnotationState((prev) => ({
        ...prev,
        uploadingMatrix: { ...prev.uploadingMatrix, [type]: false },
      }));
    }
  };

  /**
   * Checks if a CSV is available for a method
   */
  const hasCsvForMethod = useCallback(
    (type) => {
      const selectedTaskData = annotationState.predictionTasks.find(
        (t) => t.taskId === annotationState.selectedTask
      );
      const csv = selectedTaskData?.csvFiles?.[type];
      return !!csv && isValidCsvForType(csv, type) && !!annotationState.selectedTask;
    },
    [annotationState.selectedTask, annotationState.predictionTasks]
  );

  /**
   * Checks if an annotation can be launched
   */
  const canRunAnnotation = useCallback(
    (type) => {
      return hasCsvForMethod(type) && !!annotationState.uploadedMatrixNames[type];
    },
    [hasCsvForMethod, annotationState.uploadedMatrixNames]
  );

  /**
   * Launches all selected annotations
   */
  const runSelectedIntegrated = useCallback(async () => {
    const types = annotationState.selectedMethods.filter((t) =>
      canRunAnnotation(t)
    );

    if (types.length === 0) {
      addLog("No method selected or ready.", "error");
      return;
    }

    setAnnotationState((prev) => ({
      ...prev,
      showCombined: false,
      allRunning: true,
      globalProgress: 0,
      globalProgressMessage: "",
      readyToView: false,
    }));
    setShowLogs(true);

    // ===== NEW EXACT PROGRESSION LOGIC =====
    // Total = nb_molecules × nb_methods
    // Each completed annotation = +1
    // Progression = (counter / total) × 100
    const totalSourcesByMethod = {};  // nb of molecules per method (received from R)
    const completedByMethod = {};     // nb of annotated molecules per method
    const progressByMethod = {};      // progression 0-100 per method
    types.forEach(t => { 
      totalSourcesByMethod[t] = 0;    // Will be updated when R sends total_sources
      completedByMethod[t] = 0;       // Counter of completed annotations
      progressByMethod[t] = 0;        // Progression 0-100
    });
    let totalAnnotations = 0;         // Global total (will be calculated once all totals received)

    const updateGlobalProgress = (message, methodType = null, methodProgress = null) => {
      // Update method progression if provided
      if (methodType && methodProgress !== null) {
        progressByMethod[methodType] = Math.min(100, Math.max(0, methodProgress));
      }
      
      // Calculate total if we received all totals per method
      const receivedTotals = Object.values(totalSourcesByMethod).filter(v => v > 0);
      if (receivedTotals.length > 0) {
        totalAnnotations = Object.values(totalSourcesByMethod).reduce((a, b) => a + b, 0);
      }
      
      // Calculate total number of completed annotations
      const completedTotal = Object.values(completedByMethod).reduce((a, b) => a + b, 0);
      
      // Progression = (completed / total) × 100
      let globalPercent = 0;
      if (totalAnnotations > 0) {
        globalPercent = Math.round((completedTotal / totalAnnotations) * 100);
      }
      
      setAnnotationState((prev) => ({
        ...prev,
        globalProgress: Math.min(100, Math.max(0, globalPercent)),
        globalProgressMessage: message,
        // Add progression per method for display
        progressByMethod: { ...progressByMethod },
        selectedMethods: types, // Keep list of active methods
      }));
    };

    addLog(`Starting annotations (${types.join(", ")}) – task ${annotationState.selectedTask}…`);
    updateGlobalProgress('Initializing…');

    // ===== PHASE 0 : EXPORT BACKGROUNDS (no progression, just log) =====
    addLog("🎨 Exporting spectral backgrounds (grouped with automatic retry)...");

    const bgResults = {};

    try {
      const pdataPaths = {};
      const bgOptionsForExport = {};

      for (const type of types) {
        const dir = (annotationState.bgPdataDirs?.[type] || "").trim();
        if (dir) {
          try {
            const normalized = await normalizeServerPath(dir);
            pdataPaths[type] = normalized || dir;

            const opts = annotationState.bgOptions?.[type] || {};
            const cleanedOpts = {};

            const asNumber = (v) => {
              const n = Number(v);
              return Number.isFinite(n) ? n : undefined;
            };
            const pushNum = (k, v) => {
              const n = asNumber(v);
              if (Number.isFinite(n)) cleanedOpts[k] = n;
            };

            if (type === "1H" || type === "13C") {
              pushNum("ppm_min", opts.ppm_min);
              pushNum("ppm_max", opts.ppm_max);
            } else {
              pushNum("f2_ppm_min", opts.f2_ppm_min);
              pushNum("f2_ppm_max", opts.f2_ppm_max);
              pushNum("f1_ppm_min", opts.f1_ppm_min);
              pushNum("f1_ppm_max", opts.f1_ppm_max);
            }

            bgOptionsForExport[type] = cleanedOpts;
          } catch (e) {
            addLog(`[${type}] Invalid pdata path: ${e.message}`, "error");
          }
        } else {
          addLog(`[${type}] No pdata folder provided — background not exported (OK).`);
        }
      }

      if (Object.keys(pdataPaths).length > 0) {
        addLog(`[INFO] Grouped export of ${Object.keys(pdataPaths).length} spectrum(s)...`);

        const options = {
          oneH: bgOptionsForExport["1H"] || {},
          c13: bgOptionsForExport["13C"] || {},
          cosy: bgOptionsForExport["COSY"] || {},
          hsqc: bgOptionsForExport["HSQC"] || {},
        };

        const results = await exportAllSpectralBackgroundsWithRetry(
          pdataPaths,
          options,
          {
            persistKey: annotationState.selectedTask,
            overwrite: true,
          }
        );

        for (const [method, output] of Object.entries(results.outputs || {})) {
          if (output?.ok) {
            const normOut = {
              ...output,
              pngUrl: normalizeStaticUrl(output.pngUrl),
              jsonUrl: normalizeStaticUrl(output.jsonUrl),
            };
            bgResults[method] = normOut;

            setAnnotationState((prev) => ({
              ...prev,
              bgStatus: {
                ...prev.bgStatus,
                [method]: normOut,
              },
            }));

            addLog(`[${method}] Background export OK → PNG: ${normOut.pngUrl}`);
          } else {
            const errMsg = output?.error || "Export failed";
            addLog(`[${method}] Background export failed: ${errMsg}`, "error");

            setAnnotationState((prev) => ({
              ...prev,
              bgStatus: {
                ...prev.bgStatus,
                [method]: { ok: false, error: errMsg },
              },
            }));
          }
        }

        addLog(`[OK] Grouped export complete (${Object.keys(bgResults).length} backgrounds generated)`);
      } else {
        addLog("[INFO] No backgrounds to export (no pdata paths configured)");
      }
    } catch (e) {
      addLog(`[ERROR] Error during grouped background export: ${e.message}`, "error");
    }

    // ===== ANNOTATIONS =====
    addLog("[START] Launching annotations...");;

    const completed = new Set();
    const results = [];

    const promises = types.map((type) => {
      return new Promise((resolve) => {
        const lightParams = {
          ...annotationState.parameters[type],
          nb_cores: 2,
          // batch_size managed automatically on R side (adaptive based on nb molecules)
        };

        const sse = createReconnectingSSE({
          endpoint: "/run-annotation-from-prediction",
          body: {
            predictionTaskId: annotationState.selectedTask,
            matrixFile: annotationState.uploadedMatrixNames[type],
            annotationType: type,
            parameters: lightParams,
          },
          onProgress: (data) => {
            // Handle total_sources message sent by R at the beginning
            if (data?.type === 'total_sources' && Number.isFinite(Number(data?.count))) {
              totalSourcesByMethod[type] = Number(data.count);
              addLog(`[${type}] ${data.count} molecules to annotate`);
              updateGlobalProgress(`${type}: Starting...`, type, 0);
              return;
            }

            // Calculate number of completed molecules for this method
            // R sends current (source index) and total (total sources)
            const cur = Number(data?.current);
            const tot = Number(data?.total);
            
            if (Number.isFinite(cur) && Number.isFinite(tot) && tot > 0) {
              // Update total if not yet received via total_sources
              if (totalSourcesByMethod[type] === 0) {
                totalSourcesByMethod[type] = tot;
              }
              // Calculate number of completed annotations (cur represents % from 0-100)
              // So we convert to number of sources
              const completedCount = Math.round((cur / 100) * totalSourcesByMethod[type]);
              completedByMethod[type] = completedCount;
              progressByMethod[type] = cur; // cur is already in % (0-100)
            }
            
            // Recalculate and display global progression
            const methodProgress = Number.isFinite(cur) ? cur : progressByMethod[type];
            updateGlobalProgress(`${type}: ${data?.message || 'In progress…'}`, type, methodProgress);
          },
          onLog: (msg) => addLog(`[${type}] ${msg}`),
          onError: (err) => {
            addLog(`[${type}] Error: ${err}`, "error");
            completed.add(type);
            // Mark all annotations of this method as completed
            completedByMethod[type] = totalSourcesByMethod[type] || 0;
            progressByMethod[type] = 100; // Mark as completed (even on error)
            updateGlobalProgress(`${type}: Error`, type, 100);

            const bgOut = bgResults[type] || null;

            results.push({
              success: false,
              error: String(err),
              type,
              ...(bgOut ? {
                pngUrl: bgOut.pngUrl,
                jsonUrl: bgOut.jsonUrl,
                background: bgOut,
              } : {}),
            });
            runningSSERef.current.delete(type);
            resolve();
          },
          onComplete: (res) => {
            if (res?.success) {
              setAnnotationState((prev) => ({
                ...prev,
                lastResultByType: { ...prev.lastResultByType, [type]: res },
                viewerTaskId: annotationState.selectedTask,
              }));
            }
            addLog(
              res?.success ? `[${type}] OK` : `[${type}] Failed: ${res?.error}`,
              res?.success ? "info" : "error"
            );
            completed.add(type);
            // Mark all annotations of this method as completed
            completedByMethod[type] = totalSourcesByMethod[type] || 0;
            progressByMethod[type] = 100; // Mark as completed
            updateGlobalProgress(`${type}: Completed`, type, 100);

            const bgOut = bgResults[type] || null;

            results.push({
              ...(res || {}),
              type,
              ...(bgOut ? {
                pngUrl: bgOut.pngUrl,
                jsonUrl: bgOut.jsonUrl,
                background: bgOut,
              } : {}),
            });
            runningSSERef.current.delete(type);
            resolve();
          },
        });

        runningSSERef.current.set(type, sse);
      });
    });

    await Promise.all(promises);

    // Prepare viewer data
    const downloads = {};
    const csvs = {};
    const annotByType = {};

    results.forEach((r) => {
      if (r && r.success && r.downloadUrl) {
        const u = normalizeStaticUrl(r.downloadUrl);
        downloads[r.type] = u;
        csvs[r.type] = u;
      }
      annotByType[r.type] = r;
    });

    setViewerData({
      title: `Task ${annotationState.selectedTask}`,
      taskId: annotationState.selectedTask,
      csvs,
      annotationResults: annotByType,
      downloads,
      pdataDirs: annotationState.bgPdataDirs,
      bgOptions: annotationState.bgOptions,
    });

    setAnnotationState((prev) => ({
      ...prev,
      readyToView: true,
      allRunning: false,
      globalProgress: 100,
      globalProgressMessage: 'Completed',
    }));

    addLog('[OK] Workflow completed');
  }, [
    annotationState,
    setAnnotationState,
    setViewerData,
    canRunAnnotation,
    addLog,
    setShowLogs,
  ]);

  /**
   * Cancels all running annotations
   */
  const handleCancelAll = useCallback(async () => {
    setCancelRequested(true);
    addLog("Cancelling all running annotations...", "error");

    const entries = Array.from(runningSSERef.current.entries());

    if (entries.length === 0) {
      addLog("No running annotation to cancel", "error");
      setAnnotationState((prev) => ({
        ...prev,
        allRunning: false,
        globalProgressMessage: "No running task",
        globalProgress: 0,
      }));
      setCancelRequested(false);
      return;
    }

    addLog(`Cancelling ${entries.length} annotation(s)...`, "error");

    for (const [type, sse] of entries) {
      try {
        addLog(`Cancelling ${type}...`, "error");
        const taskId = sse.getTaskId?.();
        if (taskId) {
          await cancelTask(taskId);
          addLog(`${type}: task stopped (taskId: ${taskId})`, "error");
        }

        if (typeof sse.cancel === "function") {
          sse.cancel();
        }
      } catch (e) {
        addLog(`Error cancelling ${type}: ${e.message}`, "error");
      }
    }

    runningSSERef.current.clear();

    setAnnotationState((prev) => ({
      ...prev,
      allRunning: false,
      globalProgressMessage: `Cancelled (${entries.length} task(s) stopped)`,
      globalProgress: 0,
    }));

    setCancelRequested(false);
  }, [setAnnotationState, addLog]);

  // NOTE: No cleanup on unmount!
  // The hook now lives in the global context (AppStateContext)
  // so annotations continue even if the modal is closed.
  // To cancel, the user must click the "Cancel" button.

  return {
    matrixUploads,
    loadingMatrixUploads,
    cancelRequested,
    refreshMatrixUploads,
    loadPredictionTasks,
    uploadMatrix,
    hasCsvForMethod,
    canRunAnnotation,
    runSelectedIntegrated,
    handleCancelAll,
  };
};