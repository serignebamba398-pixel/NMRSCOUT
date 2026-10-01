import { isMethodAvailable } from '../constants/releaseFeatures';
import { useCallback, useEffect, useRef, useState } from "react";
import {
  API_BASE_URL,
  getPredictionTasks,
  isValidCsvForType,
  listUploads,
  uploadFile,
  saveToleranceResults,
  loadToleranceResults,
} from "../api";
import { createReconnectingSSE } from "../createReconnectingSSE";
import { ProgressTracker, normalizeStaticUrl, parseCSVSmart, shortCsvLabel } from "../utils";

/**
 * Global (persistent) hook for tolerance optimization.
 * Lives in AppStateContext so that the calculation continues even if the modal is closed.
 */
export const useToleranceOptimization = ({
  annotationState,
  setAnnotationState,
  addLog,
  setShowLogs,
  uiLock,
  setUiLock,
}) => {
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState(0);
  const [currentTest, setCurrentTest] = useState("");
  const [cancelRequested, setCancelRequested] = useState(false);

  const cancelFlagRef = useRef(false);
  const currentSseRef = useRef(null);
  const currentRejectRef = useRef(null);

  // Lock UI globally while optimization runs (even if modal is closed)
  useEffect(() => {
    if (running) {
      setUiLock({ lockedBy: "OPTIM", message: "Tolerance optimization in progress…" });
    } else if (uiLock.lockedBy === "OPTIM") {
      setUiLock({ lockedBy: null, message: "" });
    }
  }, [running, setUiLock, uiLock.lockedBy]);

  const selectedTaskData = annotationState.predictionTasks.find(
    (t) => t.taskId === annotationState.selectedTask
  );

  const hasCsvForMethod = useCallback(
    (type) => {
      const csv = selectedTaskData?.csvFiles?.[type];
      return !!csv && isValidCsvForType(csv, type) && !!annotationState.selectedTask;
    },
    [annotationState.selectedTask, selectedTaskData]
  );

  const canRunTest = useCallback(() => {
    const type = annotationState.activeType;
    if (!isMethodAvailable(type)) return false;
    return hasCsvForMethod(type) && !!annotationState.uploadedMatrixNames[type];
  }, [annotationState.activeType, annotationState.uploadedMatrixNames, hasCsvForMethod]);

  const refreshMatrixUploads = useCallback(async () => {
    // kept for parity with UI, but not required for background run
    try {
      await listUploads("matrix");
    } catch {
      // no-op
    }
  }, []);

  const loadPredictionTasks = useCallback(async () => {
    try {
      const tasks = await getPredictionTasks();
      setAnnotationState((prev) => ({ ...prev, predictionTasks: tasks }));
      if (tasks.length > 0 && !annotationState.selectedTask) {
        setAnnotationState((prev) => ({ ...prev, selectedTask: tasks[0].taskId }));
      }
    } catch (e) {
      addLog(`Error loading prediction tasks: ${e.message}`, "error");
    }
  }, [addLog, annotationState.selectedTask, setAnnotationState]);

  const loadSavedToleranceResults = useCallback(async (taskId) => {
    if (!taskId) return;
    try {
      const { results } = await loadToleranceResults(taskId);
      if (results && results.length > 0) {
        setAnnotationState((prev) => ({
          ...prev,
          toleranceResults: results,
        }));
        addLog(`[OK] Loaded ${results.length} saved tolerance optimization result(s) from task ${taskId}`, "info");
      }
    } catch (e) {
      // Don't log errors if file doesn't exist (normal case)
      if (!e.message.includes('not found') && !e.message.includes('404')) {
        addLog(`[WARN] Failed to load tolerance results: ${e.message}`, "warn");
      }
    }
  }, [addLog, setAnnotationState]);

  const uploadMatrix = useCallback(
    async (type, file) => {
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
    },
    [addLog, setAnnotationState]
  );

  const analyzeAnnotationResult = useCallback(async (csvUrl, toleranceParams, totalExperimentalSignals) => {
    try {
      const normalized = normalizeStaticUrl(csvUrl);
      const response = await fetch(normalized);
      if (!response.ok) throw new Error(`HTTP ${response.status}`);

      const csvText = await response.text();
      const { data } = parseCSVSmart(csvText);

      if (!data || data.length === 0) {
        return {
          toleranceParams,
          coverageDistribution: Array(5).fill(0),
          totalSignalsMatched: 0,
          totalExperimentalSignals: totalExperimentalSignals ?? null,
          signalsMatchedPct: 0,
          totalMolecules: 0,
          avgCoverage: 0,
          matchedSignals: [],
        };
      }

      const coverageDistribution = Array(5).fill(0);
      let totalCoverageSum = 0;
      const matchedExpSignals = new Set();

      data.forEach((row) => {
        let coverage = parseFloat(row.coverage_ransy ?? row.coverage_sub ?? 0);
        if (!Number.isFinite(coverage)) coverage = 0;
        coverage = Math.max(0, Math.min(100, coverage));

        const bin = Math.min(4, Math.floor(coverage / 20));
        coverageDistribution[bin]++;
        totalCoverageSum += coverage;

        if (row.RANSY_signals && String(row.RANSY_signals).trim() !== "") {
          const rawSignals = String(row.RANSY_signals).trim();
          // Detect 2D format (COSY/HSQC uses "||" to separate signal pairs)
          const is2D = rawSignals.includes("||");
          const expSignals = is2D
            ? rawSignals.split("||").map((s) => s.trim()).filter((s) => s !== "")
            : rawSignals.split("|").map((s) => s.trim()).filter((s) => s !== "");

          expSignals.forEach((sig) => matchedExpSignals.add(sig));
        }
      });

      const totalMolecules = data.length;
      const avgCoverage = totalCoverageSum / totalMolecules;
      const uniqueMatchedCount = matchedExpSignals.size;
      const signalsMatchedPct =
        totalExperimentalSignals && totalExperimentalSignals > 0
          ? (uniqueMatchedCount / totalExperimentalSignals) * 100
          : 0;

      return {
        toleranceParams,
        coverageDistribution,
        totalSignalsMatched: uniqueMatchedCount,
        totalExperimentalSignals: totalExperimentalSignals ?? null,
        signalsMatchedPct,
        totalMolecules,
        avgCoverage,
        matchedSignals: Array.from(matchedExpSignals)
          .map((s) => parseFloat(s))
          .filter((n) => Number.isFinite(n))
          .sort((a, b) => a - b),
      };
    } catch (e) {
      return {
        toleranceParams,
        coverageDistribution: Array(5).fill(0),
        totalSignalsMatched: 0,
        totalExperimentalSignals: totalExperimentalSignals ?? null,
        signalsMatchedPct: 0,
        totalMolecules: 0,
        avgCoverage: 0,
        matchedSignals: [],
        error: e.message,
      };
    }
  }, []);

  const runSingleAnnotation = useCallback((predTaskId, matrixFile, annotationType, parameters) => {
    return new Promise((resolve, reject) => {
      currentRejectRef.current = reject;

      const sse = createReconnectingSSE({
        endpoint: "/run-annotation-from-prediction",
        body: {
          predictionTaskId: predTaskId,
          matrixFile,
          annotationType,
          parameters,
        },
        onProgress: () => {},
        onLog: () => {},
        onTaskId: (taskId) => {
          // Store taskId for cancellation
          if (sse) {
            sse._taskId = taskId;
          }
        },
        onError: (err) => {
          if (currentRejectRef.current === reject) currentRejectRef.current = null;
          reject(new Error(String(err)));
        },
        onComplete: (res) => {
          if (currentRejectRef.current === reject) currentRejectRef.current = null;
          resolve(res);
        },
      });

      currentSseRef.current = sse;
    });
  }, []);

  const handleCancel = useCallback(async () => {
    if (!running) return;

    cancelFlagRef.current = true;
    setCancelRequested(true);
    addLog("Cancelling optimization in progress...", "error");

    try {
      // Try to cancel server-side first, then close SSE.
      const taskId = currentSseRef.current?.getTaskId?.() || currentSseRef.current?._taskId;
      if (taskId) {
        addLog(`Killing server task: ${taskId}`, "warn");
        if (currentSseRef.current?.cancelTask) {
          await currentSseRef.current.cancelTask();
          addLog(`Task ${taskId} cancelled on server`, "info");
        }
      } else {
        addLog("No taskId found - cannot kill server process", "warn");
      }

      if (currentSseRef.current?.cancel) {
        currentSseRef.current.cancel();
        addLog("SSE connection closed", "info");
      }
    } catch (e) {
      addLog(`Optimization cancellation error: ${e.message}`, "error");
    } finally {
      // unblock the await loop if necessary
      if (typeof currentRejectRef.current === "function") {
        try {
          currentRejectRef.current(new Error("Cancelled"));
        } catch {
          // ignore
        }
        currentRejectRef.current = null;
      }
    }
  }, [addLog, running]);

  const runToleranceTests = useCallback(async () => {
    if (running) return;

    const type = annotationState.activeType;
    if (!isMethodAvailable(type)) return false;
    const ranges = annotationState.toleranceRanges[type];
    const baseParams = annotationState.parameters[type];

    if (!canRunTest()) {
      addLog(`Cannot start: CSV/matrix missing for ${type}`, "error");
      return;
    }

    cancelFlagRef.current = false;
    setCancelRequested(false);
    setRunning(true);
    setProgress(0);
    setCurrentTest("Initializing…");
    setShowLogs(true);

    const computeTotalExperimentalSignalsForMatrix = async (filename) => {
      if (!filename || typeof filename !== "string") return null;
      const cleanName = filename.trim();
      if (!cleanName) return null;

      try {
        const url = `${API_BASE_URL}/uploads/${encodeURIComponent(cleanName)}`;
        const resp = await fetch(url);
        if (!resp.ok) throw new Error(`HTTP ${resp.status} (${url})`);

        const rawText = await resp.text();
        const lines = rawText
          .replace(/^\uFEFF/, "")
          .split(/\r?\n/)
          .map((l) => l.trim())
          .filter((l) => l.length > 0);

        if (lines.length === 0) return 0;

        const totalSignals = lines.length - 1;
        addLog(
          `Matrix "${cleanName}": ${lines.length} lines, ${totalSignals} signals detected`,
          "info"
        );
        return totalSignals > 0 ? totalSignals : 0;
      } catch (e) {
        addLog(`Error reading experimental matrix (${filename}): ${e.message}`, "error");
        return null;
      }
    };

    const matrixName = annotationState.uploadedMatrixNames?.[type];
    const totalExperimentalSignals = matrixName
      ? await computeTotalExperimentalSignalsForMatrix(matrixName)
      : null;

    if (totalExperimentalSignals == null) {
      addLog(
        `Cannot calculate total experimental signals for matrix "${matrixName}".`,
        "error"
      );
    } else if (totalExperimentalSignals === 0) {
      addLog(`[WARN] Matrix "${matrixName}" empty or invalid (0 signals detected).`, "warning");
    } else {
      addLog(`[OK] Matrix "${matrixName}" → ${totalExperimentalSignals} experimental signals detected.`, "info");
    }

    const testCombinations = [];

    if (type === "1H" || type === "13C") {
      for (
        let tol = ranges.tol_correspondence.min;
        tol <= ranges.tol_correspondence.max + 1e-9;
        tol += ranges.tol_correspondence.step
      ) {
        const roundedTol = parseFloat(tol.toFixed(3));
        testCombinations.push({
          tol_correspondence: roundedTol,
          tol_correspondance: roundedTol,  // spelling used by R scripts
          tol_match: roundedTol,
        });
      }
    } else if (type === "COSY" || type === "TOCSY" || type === "HSQC") {
      for (
        let tolF1 = ranges.tol_correspondence_f1.min;
        tolF1 <= ranges.tol_correspondence_f1.max + 1e-9;
        tolF1 += ranges.tol_correspondence_f1.step
      ) {
        for (
          let tolF2 = ranges.tol_correspondence_f2.min;
          tolF2 <= ranges.tol_correspondence_f2.max + 1e-9;
          tolF2 += ranges.tol_correspondence_f2.step
        ) {
          const roundedF1 = parseFloat(tolF1.toFixed(3));
          const roundedF2 = parseFloat(tolF2.toFixed(3));

          testCombinations.push({
            tol_correspondence_f1: roundedF1,
            tol_correspondence_f2: roundedF2,
            tol_correspondance_F1: roundedF1,  // spelling + case used by R scripts
            tol_correspondance_F2: roundedF2,
            tol_match_f1: roundedF1,
            tol_match_f2: roundedF2,
            tol_match_F1: roundedF1,
            tol_match_F2: roundedF2,
          });
        }
      }
    }

    addLog(`Launching ${testCombinations.length} tests for ${type}...`);

    const tracker = new ProgressTracker();
    const phases = testCombinations.map((_, i) => ({
      name: `Test ${i + 1}`,
      weight: 100 / testCombinations.length,
    }));
    tracker.definePhases(phases);

    const results = [];

    try {
      for (let i = 0; i < testCombinations.length; i++) {
        if (cancelFlagRef.current) break;

        const tolParams = testCombinations[i];
        const fullParams = { ...baseParams, ...tolParams };

        let label;
        if (type === "1H" || type === "13C") {
          label = `tol=${tolParams.tol_correspondence} ppm (match=${tolParams.tol_match})`;
        } else {
          label = `F1=${tolParams.tol_correspondence_f1} (match=${tolParams.tol_match_f1}) / F2=${tolParams.tol_correspondence_f2} (match=${tolParams.tol_match_f2})`;
        }

        setCurrentTest(`Test ${i + 1}/${testCombinations.length} : ${label}`);
        addLog(`[${i + 1}/${testCombinations.length}] ${label}...`);

        const updateTestProgress = (subPhase) => {
          const phaseProgress = (subPhase / 3) * 100;
          const globalPercent = tracker.updatePhase(i, phaseProgress);
          setProgress(globalPercent);
        };

        updateTestProgress(0);

        try {
          updateTestProgress(1);

          const result = await runSingleAnnotation(
            annotationState.selectedTask,
            annotationState.uploadedMatrixNames[type],
            type,
            fullParams
          );

          if (cancelFlagRef.current) break;

          updateTestProgress(2);

          if (result.success && result.downloadUrl) {
            const metrics = await analyzeAnnotationResult(
              result.downloadUrl,
              tolParams,
              totalExperimentalSignals
            );
            results.push(metrics);

            addLog(
              `[OK] ${label} : ${metrics.totalSignalsMatched} unique signals matched` +
                (metrics.totalExperimentalSignals
                  ? ` (${metrics.signalsMatchedPct.toFixed(2)}% of experimental matrix)`
                  : "") +
                `, avg coverage = ${metrics.avgCoverage.toFixed(2)}%`
            );
          } else {
            addLog(`✗ ${label} : Failed - ${result.error}`, "error");
            results.push({
              toleranceParams: tolParams,
              error: result.error,
              totalSignalsMatched: 0,
              totalExperimentalSignals: totalExperimentalSignals ?? null,
              signalsMatchedPct: 0,
              totalMolecules: 0,
              avgCoverage: 0,
              matchedSignals: [],
            });
          }

          updateTestProgress(3);
        } catch (e) {
          if (cancelFlagRef.current) break;
          addLog(`✗ ${label} : Error - ${e.message}`, "error");
          results.push({
            toleranceParams: tolParams,
            error: e.message,
            totalSignalsMatched: 0,
            totalExperimentalSignals: totalExperimentalSignals ?? null,
            signalsMatchedPct: 0,
            totalMolecules: 0,
            avgCoverage: 0,
            matchedSignals: [],
          });
        }
      }

      if (cancelFlagRef.current) {
        setCurrentTest("Cancelled");
        setProgress(0);
        addLog("Optimization cancelled.", "error");
        return;
      }

      const newResult = {
        id: `${type}_${Date.now()}`,
        type,
        results,
        timestamp: new Date().toISOString(),
        label: `${type} - ${new Date().toLocaleString()}`,
      };

      const updatedResults = [...annotationState.toleranceResults, newResult];

      setAnnotationState((prev) => ({
        ...prev,
        toleranceResults: updatedResults,
      }));

      // Save results to disk
      if (annotationState.selectedTask) {
        try {
          await saveToleranceResults(annotationState.selectedTask, updatedResults);
          addLog(`[OK] Tolerance results saved to outputs/${annotationState.selectedTask}/tolerance_results.json`, "info");
        } catch (e) {
          addLog(`[WARN] Failed to save tolerance results: ${e.message}`, "warn");
        }
      }

      setProgress(100);
      setCurrentTest("Completed");
      addLog(`[OK] Tests completed: ${results.length} configurations tested`);
    } finally {
      setRunning(false);
      currentSseRef.current = null;
      currentRejectRef.current = null;
      cancelFlagRef.current = false;
      setCancelRequested(false);
    }
  }, [
    addLog,
    analyzeAnnotationResult,
    annotationState,
    canRunTest,
    runSingleAnnotation,
    setAnnotationState,
    setShowLogs,
    running,
  ]);

  const getStatusLabel = useCallback(() => {
    if (!annotationState.selectedTask) return "";
    const csvMap = selectedTaskData?.csvFiles || {};
    const order = ["1H", "13C", "COSY", "TOCSY", "HSQC"].filter(isMethodAvailable);
    const ok = order.filter((t) => {
      const f = csvMap[t];
      return f && isValidCsvForType(f, t);
    });
    if (ok.length === 0) return "";
    return ok.map((t) => `${t}: ${shortCsvLabel(csvMap[t])}`).join(" • ");
  }, [annotationState.selectedTask, selectedTaskData]);

  // Load saved tolerance results when selectedTask changes
  useEffect(() => {
    if (annotationState.selectedTask) {
      loadSavedToleranceResults(annotationState.selectedTask);
    }
  }, [annotationState.selectedTask, loadSavedToleranceResults]);

  return {
    running,
    progress,
    currentTest,
    cancelRequested,
    runToleranceTests,
    handleCancel,
    // helpers used by UI
    canRunTest,
    hasCsvForMethod,
    loadPredictionTasks,
    refreshMatrixUploads,
    uploadMatrix,
    getStatusLabel,
    loadSavedToleranceResults,
  };
};
