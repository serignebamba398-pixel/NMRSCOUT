import { useState, useRef, useCallback } from 'react';
import { 
  runPredictionFromSmiles, 
  runPredictionFromCSV,
  uploadFile,
  listUploads,
  API_BASE_URL
} from '../api';

/**
 * Hook to manage prediction logic (SMILES and CSV)
 */
export const usePrediction = ({ 
  predictionState, 
  setPredictionState, 
  addLog, 
  setShowLogs 
}) => {
  // Local states
  const [running, setRunning] = useState(false);
  const [progress, setProgress] = useState(0);
  const [progressMessage, setProgressMessage] = useState("");
  const [error, setError] = useState(null);
  const [startTime, setStartTime] = useState(null);
  const [endTime, setEndTime] = useState(null);
  const [csvUploads, setCsvUploads] = useState([]);
  const [loadingCsvUploads, setLoadingCsvUploads] = useState(false);

  // Refs to avoid closure problem
  const totalMoleculesRef = useRef(1);
  const completedMoleculesRef = useRef(0);
  const currentRunRef = useRef(null);

  // States synchronized with refs
  const [totalMolecules, setTotalMolecules] = useState(1);
  const [completedMolecules, setCompletedMolecules] = useState(0);

  /**
   * Refreshes the list of uploaded CSVs
   */
  const refreshCsvUploads = useCallback(async () => {
    setLoadingCsvUploads(true);
    try {
      const list = await listUploads("csv");
      setCsvUploads(list);
    } catch (e) {
      addLog(`Error listUploads(csv): ${e.message}`, "error");
    } finally {
      setLoadingCsvUploads(false);
    }
  }, [addLog]);

  /**
   * Estimates the number of molecules in a raw CSV
   */
  const estimateMoleculeCountFromText = (text) => {
    if (!text) return 0;
    const lines = text
      .split(/\r?\n/)
      .map((l) => l.trim())
      .filter((l) => l !== "");
    if (lines.length === 0) return 0;
    return Math.max(0, lines.length - 1); // -1 pour le header
  };

  /**
   * Loads molecule count from a CSV on the server
   */
  const loadMoleculeCountFromServerCsv = async (filename) => {
    try {
      const url = `${API_BASE_URL}/uploads/${encodeURIComponent(filename)}`;
      const resp = await fetch(url);
      if (!resp.ok) throw new Error(`HTTP ${resp.status}`);
      const text = await resp.text();
      const count = estimateMoleculeCountFromText(text);
      if (!count) {
        addLog(`CSV ${filename}: no molecules detected`, "error");
      } else {
        addLog(`CSV ${filename}: ${count} molecules detected`, "info");
      }
      return count || 1;
    } catch (e) {
      addLog(`Error counting molecules in ${filename}: ${e.message}`, "error");
      return 1;
    }
  };

  /**
   * Uploads a CSV file
   */
  const uploadCsv = async (file) => {
    if (!file) return;
    setPredictionState((prev) => ({ ...prev, uploading: true }));
    try {
      const { filename } = await uploadFile(file);
      const text = await file.text();
      const nMol = estimateMoleculeCountFromText(text) || 1;

      totalMoleculesRef.current = nMol;
      completedMoleculesRef.current = 0;
      setTotalMolecules(nMol);
      setCompletedMolecules(0);

      setPredictionState((prev) => ({
        ...prev,
        csvUploadedName: filename,
        csvLocalFile: file,
      }));
      addLog(`CSV uploaded: ${filename} (${nMol} molecules detected)`);
    } catch (e) {
      addLog(`CSV upload failed: ${e.message}`, "error");
    } finally {
      setPredictionState((prev) => ({ ...prev, uploading: false }));
    }
  };

  /**
   * Builds the list of artifacts from outputs
   */
  const buildArtifactsFromOutputs = (outputs) => {
    const arr = [];
    const types = [
      { key: 'sdf', type: 'SDF' },
      { key: 'csv', type: 'CSV' },
      { key: 'json', type: 'JSON' },
      { key: 'csv_1H', type: 'CSV_1H' },
      { key: 'csv_13C', type: 'CSV_13C' },
      { key: 'csv_COSY', type: 'CSV_COSY' },
      { key: 'csv_HSQC', type: 'CSV_HSQC' },
    ];

    types.forEach(({ key, type }) => {
      if (outputs?.[key]) {
        arr.push({
          type,
          file: outputs[key].split("/").pop(),
          url: outputs[key],
        });
      }
    });

    return arr;
  };

  /**
   * Launches a SMILES prediction
   */
  const launchSmiles = useCallback(async () => {
    const { smiles, molId, options, selectedPredMethods } = predictionState;

    if (!smiles.trim()) {
      addLog("Empty SMILES", "error");
      return;
    }

    if (!selectedPredMethods.length) {
      addLog("Select at least one method.", "error");
      return;
    }

    // Initialisation
    totalMoleculesRef.current = 1;
    completedMoleculesRef.current = 0;
    setTotalMolecules(1);
    setCompletedMolecules(0);
    setProgress(0);
    setProgressMessage("Starting…");

    setPredictionState((prev) => ({
      ...prev,
      artifacts: [],
      predictionTaskId: null,
      showCombinedViewer: false,
    }));

    setShowLogs(true);
    setRunning(true);
    setError(null);
    setStartTime(Date.now());
    setEndTime(null);

    // Normalisation min/max bonds
    const minB = Number(options.cosy?.min_bonds) || 2;
    const maxB = Number(options.cosy?.max_bonds) || 3;
    const normalizedMin = Math.min(minB, maxB);
    const normalizedMax = Math.max(minB, maxB);

    const carbonModeUI = options.extractors?.c13?.carbon_mode || "nonquaternary";
    const hsqcFoldUI = options.extractors?.hsqc?.fold ?? true;
    const hsqcLowUI = options.extractors?.hsqc?.fold_low ?? 4.8618;
    const hsqcHighUI = options.extractors?.hsqc?.fold_high ?? 69.8395;

    // Experimental parameters for multiplicity JSON generation
    const freqMhz = Number(options.freq_mhz) || 700;
    const solvent = options.solvent || 'CDCl3';
    const temperatureK = Number(options.temperature_k) || 298;
    
    // ✅ VALIDATION pH ULTRA-ROBUSTE - Prévention bug pH 26 (SMILES mode)
    let pH = Number(options.ph) || 7.4;
    
    // Log pour debugging - voir toutes les valeurs pH reçues
    console.log(`[DEBUG usePrediction SMILES] pH reçu: ${options.ph} (type: ${typeof options.ph}) → Converti: ${pH}`);
    
    if (pH === 26) {
      console.warn('[usePrediction SMILES] DANGER: pH=26 détecté - probablement confusion avec année 2026 ou température! Fallback vers 7.4');
      addLog('⚠️ pH=26 détecté et corrigé vers 7.4 (confusion probable avec année/température)', 'warning');
      pH = 7.4;
    } else if (pH === 2026) {
      console.warn('[usePrediction SMILES] DANGER: pH=2026 détecté - probablement confusion avec année! Fallback vers 7.4');
      addLog('⚠️ pH=2026 détecté et corrigé vers 7.4 (confusion probable avec année)', 'warning');
      pH = 7.4;
    } else if (pH === 15) {
      console.warn('[usePrediction SMILES] DANGER: pH=15 détecté - probablement confusion avec température! Fallback vers 7.4');
      addLog('⚠️ pH=15 détecté et corrigé vers 7.4 (confusion probable avec température)', 'warning');
      pH = 7.4;
    } else if (pH === 16) {
      console.warn('[usePrediction SMILES] DANGER: pH=16 détecté - probablement confusion avec température! Fallback vers 7.4');
      addLog('⚠️ pH=16 détecté et corrigé vers 7.4 (confusion probable avec température)', 'warning');
      pH = 7.4;
    } else if (pH >= 15 && pH <= 30) {
      console.warn(`[usePrediction SMILES] DANGER: pH=${pH} dans plage suspecte [15-30] - probablement confusion! Fallback vers 7.4`);
      addLog(`⚠️ pH=${pH} détecté et corrigé vers 7.4 (valeur suspecte)`, 'warning');
      pH = 7.4;
    } else if (pH < 1.0 || pH > 14.0) {
      console.warn(`[usePrediction SMILES] pH=${pH} hors plage valide [1.0-14.0] - Fallback vers 7.4`);
      addLog(`⚠️ pH=${pH} hors plage et corrigé vers 7.4`, 'warning');
      pH = 7.4;
    }

    try {
      const runPromise = runPredictionFromSmiles(
        smiles,
        {
          molId: molId || "mol",
          methods: selectedPredMethods,
          global: {
            nb_cores: Number(options.nb_cores) || 2,
            min_bonds: normalizedMin,
            max_bonds: normalizedMax,
            // Experimental parameters
            freq_mhz: freqMhz,
            solvent: solvent,
            temperature_k: temperatureK,
            ph: pH,
            extractors: {
              c13: { carbon_mode: carbonModeUI },
              hsqc: {
                fold: hsqcFoldUI,
                fold_low: hsqcLowUI,
                fold_high: hsqcHighUI,
              },
            },
          },
        },
        // onProgress
        (data) => {
          // Total number of molecules (from totalMolecules, not total which is 100)
          const tFromData = Number(data?.totalMolecules ?? data?.total);
          if (Number.isFinite(tFromData) && tFromData > 0) {
            totalMoleculesRef.current = tFromData;
            setTotalMolecules(tFromData);
          }

          // Number of processed molecules (from completedMolecules, not current which is %)
          const cFromData = Number(data?.completedMolecules ?? data?.current);
          if (Number.isFinite(cFromData) && cFromData >= 0) {
            completedMoleculesRef.current = cFromData;
            setCompletedMolecules(cFromData);
          }

          // Progression percentage (0-100)
          const pFromData = Number(data?.percent);
          if (Number.isFinite(pFromData) && pFromData >= 0 && pFromData <= 100) {
            // Always use backend progress - it properly manages SDF (0-95%) and CSV extractors (95-100%)
            setProgress(pFromData);
          }

          // Phase-specific messages for CSV extractors
          if (data?.phase === 'extractors') {
            const extractorMsg = data?.message || 'Generating CSV files (1H, 13C, COSY, HSQC)...';
            setProgressMessage(extractorMsg);
          } else if (data?.message) {
            setProgressMessage(data.message);
          }
        },
        // onLog
        (msg) => {
          addLog(`[PRED] ${msg}`);

          // Detect "[OK] Molecule X/Y processed"
          const molMatch = msg.match(/\[OK\]\s+molécule\s+(\d+)\/(\d+)\s+traitée/i)
                        || msg.match(/\[OK\]\s+molecule\s+(\d+)\/(\d+)\s+processed/i);
          if (molMatch) {
            const current = parseInt(molMatch[1], 10);
            const total = parseInt(molMatch[2], 10);

            totalMoleculesRef.current = total;
            completedMoleculesRef.current = current;
            setTotalMolecules(total);
            setCompletedMolecules(current);

            // Only update progress from logs if backend hasn't sent progress events
            // Backend sends: SDF (0-98%), CSV extractors (98-100%)
            const frac = total > 0 ? current / total : 0;
            const newPercent = Math.min(98, frac * 98); // Cap at 98% for SDF phase

            setProgress((prev) => Math.max(prev, newPercent));
            setProgressMessage(`Molecule ${current}/${total} (SDF created)`);
            return;
          }

          // Detect CSV phase messages (for logging only, backend handles progress 95-100%)
          const isCsvPhase =
            /extraction.*csv/i.test(msg) ||
            /post.*traitement/i.test(msg) ||
            /generating csv/i.test(msg) ||
            /predict_1h|predict_13c|predict_cosy|predict_hsqc/i.test(msg) ||
            /csv.*en cours/i.test(msg);

          if (isCsvPhase) {
            // Just update message, let backend handle progress percentage
            setProgressMessage("Generating CSV files (1H, 13C, COSY, HSQC)...");
          }
        },
        // onError
        (err) => {
          addLog(`[PRED] Error: ${err}`, "error");
          setError(err);
          setRunning(false);
          setEndTime(Date.now());
        },
        // onComplete
        (res) => {
          setRunning(false);
          setEndTime(Date.now());

          if (res?.success) {
            setProgress(100);
            setProgressMessage("[OK] Completed (SDF + CSV)");

            const arts = buildArtifactsFromOutputs(res.outputs);
            setPredictionState((prev) => ({
              ...prev,
              predictionTaskId: res.taskId || null,
              artifacts: arts,
              showCombinedViewer: true,
            }));
          } else {
            setError(res?.error || "Failed");
          }
        }
      );

      currentRunRef.current = runPromise;
      await runPromise;
    } catch (e) {
      addLog(`Prediction error: ${e.message}`, "error");
      setError(e.message);
      setRunning(false);
      setEndTime(Date.now());
    }
  }, [predictionState, setPredictionState, addLog, setShowLogs]);

  /**
   * Launches a CSV batch prediction
   */
  const launchCSV = useCallback(async () => {
    const { csvUploadedName, options, selectedPredMethods } = predictionState;

    if (!csvUploadedName) {
      addLog("No CSV selected", "error");
      return;
    }

    if (!selectedPredMethods.length) {
      addLog("Select at least one method.", "error");
      return;
    }

    let nMol = totalMoleculesRef.current;
    if (!Number.isFinite(nMol) || nMol <= 0) {
      nMol = 1;
    }

    totalMoleculesRef.current = nMol;
    completedMoleculesRef.current = 0;
    setTotalMolecules(nMol);
    setCompletedMolecules(0);
    setProgress(0);
    setProgressMessage(`Starting… (~${nMol} molecules)`);

    setPredictionState((prev) => ({
      ...prev,
      artifacts: [],
      predictionTaskId: null,
      showCombinedViewer: false,
    }));

    setShowLogs(true);
    setRunning(true);
    setError(null);
    setStartTime(Date.now());
    setEndTime(null);

    const minB = Number(options.cosy?.min_bonds) || 2;
    const maxB = Number(options.cosy?.max_bonds) || 3;
    const normalizedMin = Math.min(minB, maxB);
    const normalizedMax = Math.max(minB, maxB);

    const carbonModeUI = options.extractors?.c13?.carbon_mode || "nonquaternary";
    const hsqcFoldUI = options.extractors?.hsqc?.fold ?? true;
    const hsqcLowUI = options.extractors?.hsqc?.fold_low ?? 4.8618;
    const hsqcHighUI = options.extractors?.hsqc?.fold_high ?? 69.8395;

    // Experimental parameters for multiplicity JSON generation
    const freqMhz = Number(options.freq_mhz) || 700;
    const solvent = options.solvent || 'CDCl3';
    const temperatureK = Number(options.temperature_k) || 298;
    
    // ✅ VALIDATION pH ULTRA-ROBUSTE - Prévention bug pH 26 (CSV mode)
    let pH = Number(options.ph) || 7.4;
    
    // Log pour debugging - voir toutes les valeurs pH reçues
    console.log(`[DEBUG usePrediction CSV] pH reçu: ${options.ph} (type: ${typeof options.ph}) → Converti: ${pH}`);
    
    if (pH === 26) {
      console.warn('[usePrediction CSV] DANGER: pH=26 détecté - probablement confusion avec année 2026 ou température! Fallback vers 7.4');
      addLog('⚠️ pH=26 détecté et corrigé vers 7.4 (confusion probable avec année/température)', 'warning');
      pH = 7.4;
    } else if (pH === 2026) {
      console.warn('[usePrediction CSV] DANGER: pH=2026 détecté - probablement confusion avec année! Fallback vers 7.4');
      addLog('⚠️ pH=2026 détecté et corrigé vers 7.4 (confusion probable avec année)', 'warning');
      pH = 7.4;
    } else if (pH === 15) {
      console.warn('[usePrediction CSV] DANGER: pH=15 détecté - probablement confusion avec température! Fallback vers 7.4');
      addLog('⚠️ pH=15 détecté et corrigé vers 7.4 (confusion probable avec température)', 'warning');
      pH = 7.4;
    } else if (pH === 16) {
      console.warn('[usePrediction CSV] DANGER: pH=16 détecté - probablement confusion avec température! Fallback vers 7.4');
      addLog('⚠️ pH=16 détecté et corrigé vers 7.4 (confusion probable avec température)', 'warning');
      pH = 7.4;
    } else if (pH >= 15 && pH <= 30) {
      console.warn(`[usePrediction CSV] DANGER: pH=${pH} dans plage suspecte [15-30] - probablement confusion! Fallback vers 7.4`);
      addLog(`⚠️ pH=${pH} détecté et corrigé vers 7.4 (valeur suspecte)`, 'warning');
      pH = 7.4;
    } else if (pH < 1.0 || pH > 14.0) {
      console.warn(`[usePrediction CSV] pH=${pH} hors plage valide [1.0-14.0] - Fallback vers 7.4`);
      addLog(`⚠️ pH=${pH} hors plage et corrigé vers 7.4`, 'warning');
      pH = 7.4;
    }

    try {
      const runPromise = runPredictionFromCSV(
        csvUploadedName,
        {
          idColumn: options.idColumn || "id",
          smilesColumn: options.smilesColumn || "smiles",
          methods: selectedPredMethods,
          // Reprise incrémentale: réutiliser les molécules déjà prédites
          resumeFromDir: options.resumeFromDir || null,
          global: {
            nb_cores: Number(options.nb_cores) || 2,
            min_bonds: normalizedMin,
            max_bonds: normalizedMax,
            // Experimental parameters
            freq_mhz: freqMhz,
            solvent: solvent,
            temperature_k: temperatureK,
            ph: pH,
            extractors: {
              c13: { carbon_mode: carbonModeUI },
              hsqc: {
                fold: hsqcFoldUI,
                fold_low: hsqcLowUI,
                fold_high: hsqcHighUI,
              },
            },
          },
        },
        // onProgress
        (data) => {
          // Total number of molecules (from totalMolecules, not total which is 100)
          const tFromData = Number(data?.totalMolecules ?? data?.total);
          if (Number.isFinite(tFromData) && tFromData > 0) {
            totalMoleculesRef.current = tFromData;
            setTotalMolecules(tFromData);
          }

          // Number of processed molecules (from completedMolecules, not current which is %)
          const cFromData = Number(data?.completedMolecules ?? data?.current);
          if (Number.isFinite(cFromData) && cFromData >= 0) {
            completedMoleculesRef.current = cFromData;
            setCompletedMolecules(cFromData);
          }

          // Progression percentage (0-100)
          const pFromData = Number(data?.percent);
          if (Number.isFinite(pFromData) && pFromData >= 0 && pFromData <= 100) {
            setProgress(pFromData);
          }

          if (data?.message) {
            setProgressMessage(data.message);
          }
        },
        // onLog (same logic as SMILES)
        (msg) => {
          addLog(`[PRED] ${msg}`);

          const molMatch = msg.match(/\[OK\]\s+molécule\s+(\d+)\/(\d+)\s+traitée/i)
                        || msg.match(/\[OK\]\s+molecule\s+(\d+)\/(\d+)\s+processed/i);
          if (molMatch) {
            const current = parseInt(molMatch[1], 10);
            const total = parseInt(molMatch[2], 10);

            totalMoleculesRef.current = total;
            completedMoleculesRef.current = current;
            setTotalMolecules(total);
            setCompletedMolecules(current);

            const frac = total > 0 ? current / total : 0;
            const newPercent = Math.min(98, frac * 98);

            setProgress((prev) => Math.max(prev, newPercent));
            setProgressMessage(`Molecule ${current}/${total} (SDF created)`);
            return;
          }

          // Detect CSV/extractor phase - just update message, backend sends percent via onProgress
          const isExtractorPhase = /extraction.*csv|post.*traitement|csv.*en cours|predict_.*_real|1H.*csv|13C.*csv|COSY.*csv|HSQC.*csv/i.test(msg);
          if (isExtractorPhase) {
            setProgressMessage("Generating CSV files...");
            // Don't modify percent here - backend sends correct value via onProgress with phase='extractors'
          }
        },
        // onError
        (err) => {
          addLog(`[PRED] Error: ${err}`, "error");
          setError(err);
          setRunning(false);
          setEndTime(Date.now());
        },
        // onComplete
        (res) => {
          setRunning(false);
          setEndTime(Date.now());

          if (res?.success) {
            setProgress(100);
            setProgressMessage("[OK] Completed (SDF + CSV)");

            const arts = buildArtifactsFromOutputs(res.outputs);
            setPredictionState((prev) => ({
              ...prev,
              predictionTaskId: res.taskId || null,
              artifacts: arts,
              showCombinedViewer: true,
            }));
          } else {
            setError(res?.error || "Failed");
          }
        }
      );

      currentRunRef.current = runPromise;
      await runPromise;
    } catch (e) {
      addLog(`Batch prediction error: ${e.message}`, "error");
      setError(e.message);
      setRunning(false);
      setEndTime(Date.now());
    }
  }, [predictionState, setPredictionState, addLog, setShowLogs]);

  /**
   * Cancels the current prediction
   */
  const handleCancel = useCallback(async () => {
    const taskId = currentRunRef.current?.getTaskId?.();
    if (!taskId) {
      addLog("Cannot cancel: taskId not yet received from server", "error");
      return;
    }

    addLog(`Cancelling task ${taskId}...`, "error");
    try {
      if (typeof currentRunRef.current.cancelTask === "function") {
        await currentRunRef.current.cancelTask();
        addLog(`Server process stopped (taskId: ${taskId})`, "error");
      }
      if (typeof currentRunRef.current.cancel === "function") {
        currentRunRef.current.cancel();
        addLog("SSE connection closed", "error");
      }
      setError("Cancelled by user");
      setRunning(false);
      setEndTime(Date.now());
    } catch (e) {
      addLog(`Cancellation error: ${e.message}`, "error");
    }
  }, [addLog]);

  return {
    running,
    progress,
    progressMessage,
    error,
    startTime,
    endTime,
    totalMolecules,
    completedMolecules,
    csvUploads,
    loadingCsvUploads,
    refreshCsvUploads,
    uploadCsv,
    loadMoleculeCountFromServerCsv,
    launchSmiles,
    launchCSV,
    handleCancel,
    // Unified function to launch prediction (automatically chooses SMILES or CSV)
    runPrediction: predictionState.mode === "CSV" ? launchCSV : launchSmiles,
  };
};