// src/components/viewer/SpectrumDetailModal.js
// Detailed simulated 1H spectrum with line shapes.

import React, { useRef, useEffect, useState, useCallback } from "react";
import { simulateSpectrum, getMultiplicityData, extractBrukerSpectrum, getTaskParams } from "../../api";

const MODAL_STYLES = {
  overlay: {
    position: "fixed",
    top: 0,
    left: 0,
    right: 0,
    bottom: 0,
    backgroundColor: "rgba(0, 0, 0, 0.85)",
    display: "flex",
    alignItems: "center",
    justifyContent: "center",
    zIndex: 10000,
  },
  modal: {
    backgroundColor: "#fff",
    borderRadius: 4,
    padding: 20,
    width: "99vw",
    height: "98vh",
    maxWidth: "none",
    maxHeight: "none",
    overflow: "hidden",
    display: "flex",
    flexDirection: "column",
    boxShadow: "0 20px 60px rgba(0,0,0,0.5)",
  },
  header: {
    display: "flex",
    justifyContent: "space-between",
    alignItems: "center",
    marginBottom: 16,
    paddingBottom: 12,
    borderBottom: "1px solid #e5e7eb",
  },
  title: {
    fontSize: 18,
    fontWeight: 700,
    color: "#1e293b",
  },
  closeBtn: {
    background: "none",
    border: "none",
    fontSize: 28,
    cursor: "pointer",
    color: "#64748b",
    padding: "4px 12px",
    borderRadius: 4,
    lineHeight: 1,
  },
  controls: {
    display: "flex",
    gap: 16,
    marginBottom: 12,
    flexWrap: "wrap",
    alignItems: "center",
  },
  controlGroup: {
    display: "flex",
    alignItems: "center",
    gap: 6,
  },
  label: {
    fontSize: 12,
    fontWeight: 500,
    color: "#64748b",
  },
  select: {
    padding: "4px 8px",
    borderRadius: 4,
    border: "1px solid #d1d5db",
    fontSize: 12,
  },
  input: {
    padding: "4px 8px",
    borderRadius: 4,
    border: "1px solid #d1d5db",
    fontSize: 12,
    width: 70,
  },
  canvasContainer: {
    flex: 1,
    minHeight: 0,
    position: "relative",
    backgroundColor: "#f8fafc",
    borderRadius: 8,
    overflow: "hidden",
  },
  loading: {
    position: "absolute",
    top: "50%",
    left: "50%",
    transform: "translate(-50%, -50%)",
    display: "flex",
    flexDirection: "column",
    alignItems: "center",
    gap: 12,
  },
  spinner: {
    width: 40,
    height: 40,
    border: "3px solid #e5e7eb",
    borderTop: "3px solid #3b82f6",
    borderRadius: "50%",
    animation: "spin 1s linear infinite",
  },
  error: {
    color: "#dc2626",
    textAlign: "center",
    padding: 20,
  },
  peaksTable: {
    marginTop: 8,
    maxHeight: 100,
    overflow: "auto",
    fontSize: 11,
    border: "1px solid #e5e7eb",
    borderRadius: 6,
  },
  table: {
    width: "100%",
    borderCollapse: "collapse",
  },
  th: {
    padding: "6px 8px",
    textAlign: "left",
    backgroundColor: "#f1f5f9",
    fontWeight: 600,
    position: "sticky",
    top: 0,
  },
  td: {
    padding: "4px 8px",
    borderTop: "1px solid #e5e7eb",
  },
};

// CSS keyframes for spinner
const spinnerKeyframes = `
  @keyframes spin {
    0% { transform: rotate(0deg); }
    100% { transform: rotate(360deg); }
  }
`;

function SpectrumDetailModal({
  isOpen,
  onClose,
  molId,
  taskId,
  theoryShifts = [],
  expShifts = [],
  pdataDir1H = null,
  matchTolerancePpm = 0.05,  // Default annotation tolerance.
  multiplicityJson = null,   // Preloaded multiplicity data.
  preloadedSpectrumData = null,  // Preloaded simulated spectrum.
  isPreloading = false       // Whether preloading is in progress.
}) {
  const canvasSimRef = useRef(null);  // Simulated spectrum canvas (top).
  const canvasExpRef = useRef(null);  // Experimental spectrum canvas (bottom).
  const [loading, setLoading] = useState(false);
  const [error, setError] = useState(null);
  const [spectrumData, setSpectrumData] = useState(null);
  const [fetchedMultiplicityData, setFetchedMultiplicityData] = useState(null);
  const [loadingMultiplicity, setLoadingMultiplicity] = useState(false);
  const [multiplicityError, setMultiplicityError] = useState(null);

  // Experimental Bruker data.
  const [expSpectrumData, setExpSpectrumData] = useState(null);
  const [loadingExpSpectrum, setLoadingExpSpectrum] = useState(false);
  const [expSpectrumError, setExpSpectrumError] = useState(null);

  // Simulation options
  const [freqMhz, setFreqMhzInternal] = useState(700);
  const [lineShape, setLineShapeInternal] = useState("lorentzian");

  // Choose the default line width from the spectrometer frequency.
  const getDefaultLineWidth = useCallback((freq) => {
    // Keep these values aligned with the backend FREQ_LINEWIDTH_TABLE.
    if (freq <= 150) return 1.80;
    if (freq <= 250) return 1.50;
    if (freq <= 350) return 1.30;
    if (freq <= 450) return 1.20;
    if (freq <= 600) return 1.00;
    if (freq <= 800) return 0.90;
    return 0.80;
  }, []);

  const [lineWidthHz, setLineWidthHzInternal] = useState(() => getDefaultLineWidth(700));
  const [isLineWidthAuto, setIsLineWidthAuto] = useState(true); // Track whether the value was set automatically or manually.

  // Dynamics settings: solvent, temperature, and pH.
  // The dynamics module handles solvent-dependent visibility of OH and NH protons.
  const [solvent, setSolventInternal] = useState('CDCl3');  // Default solvent.
  const [temperature, setTemperatureInternal] = useState(298);  // Temperature in kelvin.
  const [ph, setPhInternal] = useState(7.4);

  // Track changes that require a new simulation.
  const [userChangedParams, setUserChangedParams] = useState(false);

  // Update settings without automatically rerunning the simulation.
  // Simulation starts when the user clicks Refresh.
  const setFreqMhz = useCallback((v) => {
    setFreqMhzInternal(v);
    // Recalculate the line width in automatic mode.
    if (isLineWidthAuto) {
      setLineWidthHzInternal(getDefaultLineWidth(v));
    }
    
  }, [isLineWidthAuto, getDefaultLineWidth]);

  const setLineShape = useCallback((v) => { setLineShapeInternal(v); }, []);

  const setLineWidthHz = useCallback((v) => {
    setLineWidthHzInternal(v);
    setIsLineWidthAuto(false); // Switch to manual mode.
    
  }, []);

  const resetLineWidthToAuto = useCallback(() => {
    setLineWidthHzInternal(getDefaultLineWidth(freqMhz));
    setIsLineWidthAuto(true);
    
  }, [freqMhz, getDefaultLineWidth]);

  const setSolvent = useCallback((v) => { setSolventInternal(v); }, []);
  const setTemperature = useCallback((v) => { setTemperatureInternal(v); }, []);
  const setPh = useCallback((v) => { setPhInternal(v); }, []);

  // Zoom state
  const [zoomRange, setZoomRange] = useState(null);
  const [verticalScaleSim, setVerticalScaleSim] = useState(1.0); // Simulated spectrum intensity scale.
  const [verticalScaleExp, setVerticalScaleExp] = useState(1.0); // Experimental spectrum intensity scale.
  const [ppmOffsetExp, setPpmOffsetExp] = useState(0.0); // Experimental spectrum offset in ppm.
  const [isDragging, setIsDragging] = useState(false);
  const [dragStart, setDragStart] = useState(null);

  // Reset userChangedParams when the modal closes.
  useEffect(() => {
    if (!isOpen) {
      setUserChangedParams(false);
    }
  }, [isOpen]);

  // Reset data when the selected molecule changes.
  useEffect(() => {
    console.log('[SpectrumDetailModal] molId changed to:', molId, '- resetting data');
    setSpectrumData(null);
    setFetchedMultiplicityData(null);
    setMultiplicityError(null);
    setZoomRange(null);
    setError(null);
  }, [molId]);

  // Annotation tolerance loaded from the task settings.
  const [taskTolerancePpm, setTaskTolerancePpm] = useState(matchTolerancePpm);

  // Load annotation settings for this task.
  useEffect(() => {
    if (isOpen && taskId) {
      getTaskParams(taskId)
        .then(params => {
          // Read the tolerance from deltaPpm or ppmTol.
          const tolerance = params.deltaPpm || params.ppmTol || params.match_tolerance_ppm || matchTolerancePpm;
          setTaskTolerancePpm(tolerance);
          console.log('[SpectrumDetailModal] Tolérance d\'annotation:', tolerance, 'ppm');
        })
        .catch(err => {
          console.warn('[SpectrumDetailModal] Failed to get task params:', err);
          setTaskTolerancePpm(matchTolerancePpm); // Use the default if no task-specific tolerance is available.
        });
    }
  }, [isOpen, taskId, matchTolerancePpm]);

  // Fetch multiplicity data from task directory when modal opens
  useEffect(() => {
    if (isOpen && taskId && molId && !multiplicityJson) {
      console.log('[SpectrumDetailModal] Fetching multiplicity data for:', { taskId, molId });
      setLoadingMultiplicity(true);
      setMultiplicityError(null);

      getMultiplicityData(taskId, molId)
        .then(data => {
          console.log('[SpectrumDetailModal] Fetched multiplicity data:', data);
          setFetchedMultiplicityData(data.multiplicity_data);
          setLoadingMultiplicity(false);
          setMultiplicityError(null);
        })
        .catch(err => {
          console.warn('[SpectrumDetailModal] Failed to fetch multiplicity data:', err);
          const errorMsg = err.response?.data?.error || err.message || 'Unknown error';
          const errorDetails = err.response?.data ? JSON.stringify(err.response.data, null, 2) : '';
          setMultiplicityError(`${errorMsg}${errorDetails ? '\n\nDetails:\n' + errorDetails : ''}`);
          setFetchedMultiplicityData(null);
          setLoadingMultiplicity(false);
        });
    }
  }, [isOpen, taskId, molId, multiplicityJson]);

  // Initialize parameters from multiplicity JSON
  useEffect(() => {
    const multData = multiplicityJson || fetchedMultiplicityData;
    if (multData && !userChangedParams) {
      // Load experimental parameters from multiplicity JSON
      if (multData.freq_MHz) setFreqMhzInternal(multData.freq_MHz);
      if (multData.solvent) setSolventInternal(multData.solvent);
      if (multData.temperature_K) setTemperatureInternal(multData.temperature_K);
      // pH not in JSON, keep default
      console.log('[SpectrumDetailModal] Loaded parameters from multiplicity JSON:', {
        freq: multData.freq_MHz,
        solvent: multData.solvent,
        temp: multData.temperature_K
      });
    }
  }, [multiplicityJson, fetchedMultiplicityData, userChangedParams]);

  // Fetch experimental Bruker spectrum data when modal opens
  useEffect(() => {
    if (isOpen && pdataDir1H) {
      console.log('[SpectrumDetailModal] Extracting Bruker spectrum from:', pdataDir1H);
      setLoadingExpSpectrum(true);
      setExpSpectrumError(null);

      extractBrukerSpectrum(pdataDir1H)
        .then(data => {
          // Use reduce to find min/max to avoid stack overflow with large arrays
          const ppmMin = data.ppm.reduce((min, val) => Math.min(min, val), Infinity);
          const ppmMax = data.ppm.reduce((max, val) => Math.max(max, val), -Infinity);

          console.log('[SpectrumDetailModal] Extracted experimental spectrum:', {
            n_points: data.n_points,
            sf: data.sf,
            ppm_range: [ppmMin, ppmMax]
          });
          setExpSpectrumData(data);
          setLoadingExpSpectrum(false);
        })
        .catch(err => {
          console.warn('[SpectrumDetailModal] Failed to extract Bruker spectrum:', err);
          setExpSpectrumError(err.message || 'Failed to load experimental spectrum');
          setExpSpectrumData(null);
          setLoadingExpSpectrum(false);
        });
    }
  }, [isOpen, pdataDir1H]);

  // Fetch spectrum data
  const fetchSpectrum = useCallback(async () => {
    // Prefer multiplicity data passed as a prop over fetched data.
    const multData = multiplicityJson || fetchedMultiplicityData;

    // Wait for multiplicity data before simulating.
    if (loadingMultiplicity) {
      console.log('[SpectrumDetailModal] Still loading multiplicity data, waiting...');
      return;
    }

    // Simulation requires multiplicity JSON.
    if (!multData) {
      const errorMessage = multiplicityError
        ? `Failed to load multiplicity data: ${multiplicityError}`
        : "No multiplicity data available. Please load a prediction task.";
      setError(errorMessage);
      return;
    }

    setLoading(true);
    setError(null);

    try {
      const data = await simulateSpectrum({
        // Pass multiplicity JSON to the simulator.
        multiplicity_json: multData,
        mol_id: molId,
        freq_mhz: freqMhz,
        line_shape: lineShape,
        // ppm_range is auto-calculated by backend from molecule shifts
        digital_resolution_hz: 0.1,
        line_width_hz: lineWidthHz, // Use the supplied line width or let the backend choose it.
        // The dynamics module handles solvent-dependent OH and NH visibility.
        enable_dynamics: true,  // Keep enabled so the corrections are applied.
        solvent: solvent,
        temperature_k: temperature,
        ph: ph,
      });

      setSpectrumData(data);
      setZoomRange(null);  // Reset zoom on new data
    } catch (err) {
      console.error("[SpectrumDetailModal] Error:", err);
      // Show detailed error for debugging in production
      let errorMsg = err.message || "Failed to simulate spectrum";
      if (err.response?.data) {
        console.error("[SpectrumDetailModal] Server error details:", err.response.data);
        errorMsg = `${errorMsg}\n\nDetails: ${JSON.stringify(err.response.data, null, 2)}`;
      }
      setError(errorMsg);
    } finally {
      setLoading(false);
    }
  }, [molId, freqMhz, lineShape, multiplicityJson, fetchedMultiplicityData, lineWidthHz, solvent, temperature, ph, loadingMultiplicity, multiplicityError]);

  // Trigger simulation ONLY on modal open (not on param changes)
  useEffect(() => {
    // Wait for data loading to finish before starting the simulation.
    if (loadingMultiplicity || isPreloading) {
      console.log('[SpectrumDetailModal] Waiting for multiplicity data to load...');
      return;
    }

    //  Run automatically only when the molecule is first loaded.
    // Changing settings does not automatically run the simulation.
    // Click Refresh to rerun the simulation.
    if (isOpen && (multiplicityJson || fetchedMultiplicityData)) {
      if (!spectrumData && preloadedSpectrumData) {
        // Use preloaded data only for the initial render.
        console.log('[SpectrumDetailModal] Using preloaded data (first load)');
        setSpectrumData(preloadedSpectrumData);
        setLoading(false);
        setError(null);
      } else if (!spectrumData) {
        // Run the initial simulation when no preloaded spectrum is available.
        console.log('[SpectrumDetailModal] Initial spectrum simulation');
        fetchSpectrum();
      }
      
      // Apply parameter changes when the user clicks Refresh.
    }
  }, [isOpen, loadingMultiplicity, isPreloading, multiplicityJson, fetchedMultiplicityData, preloadedSpectrumData, spectrumData, fetchSpectrum]);
  // Parameter changes are intentionally excluded from this effect; Refresh applies them.
  // Parameter changes should not trigger simulation automatically.

  // Draw spectrum on canvas
  useEffect(() => {
    if (!spectrumData || !canvasSimRef.current) return;
    
    const canvas = canvasSimRef.current;
    const ctx = canvas.getContext("2d");
    const DPR = window.devicePixelRatio || 1;
    
    const container = canvas.parentElement;
    const W = container.clientWidth || 800;
    const H = container.clientHeight || 400;
    
    canvas.width = W * DPR;
    canvas.height = H * DPR;
    canvas.style.width = `${W}px`;
    canvas.style.height = `${H}px`;
    ctx.scale(DPR, DPR);
    
    const PAD = { top: 20, right: 30, bottom: 40, left: 50 };
    const plotW = W - PAD.left - PAD.right;
    const plotH = H - PAD.top - PAD.bottom;
    
    // Clear
    ctx.fillStyle = "#fff";
    ctx.fillRect(0, 0, W, H);
    
    // Determine PPM range (use zoom if set, otherwise use data range)
    const defaultRange = spectrumData.ppm_range
      ? { min: spectrumData.ppm_range[0], max: spectrumData.ppm_range[1] }
      : { min: Math.min(...spectrumData.ppm), max: Math.max(...spectrumData.ppm) };
    const range = zoomRange || defaultRange;
    const ppmToX = (ppm) => PAD.left + ((range.max - ppm) / (range.max - range.min)) * plotW;
    
    // Find data indices for current range
    const ppmArr = spectrumData.ppm;
    const intArr = spectrumData.intensity;
    
    // Filter and get max intensity in range
    let maxInt = 0;
    const visibleData = [];
    for (let i = 0; i < ppmArr.length; i++) {
      const ppm = ppmArr[i];
      if (ppm >= range.min && ppm <= range.max) {
        visibleData.push({ ppm, int: intArr[i] });
        if (intArr[i] > maxInt) maxInt = intArr[i];
      }
    }
    if (maxInt === 0) maxInt = 1;
    
    // Apply the simulated spectrum intensity scale.
    const intToY = (int) => PAD.top + plotH - (int / maxInt) * plotH * 0.9 * verticalScaleSim;
    
    // Draw grid with adaptive step based on zoom level
    ctx.strokeStyle = "#e5e7eb";
    ctx.lineWidth = 1;
    
    // Vertical grid lines (ppm) - adaptive step
    const ppmRange = range.max - range.min;
    let ppmStep;
    if (ppmRange < 0.005) ppmStep = 0.0005;     // 0.0005 ppm ticks correspond to 0.5 Hz at 1000 MHz.
    else if (ppmRange < 0.01) ppmStep = 0.001;  // 0.001 ppm tick spacing.
    else if (ppmRange < 0.05) ppmStep = 0.005;  // 0.005 ppm tick spacing.
    else if (ppmRange < 0.1) ppmStep = 0.01;    // 0.01 ppm tick spacing.
    else if (ppmRange < 0.5) ppmStep = 0.05;    // 0.05 ppm tick spacing.
    else if (ppmRange < 1) ppmStep = 0.1;       // Zoom moyen: 0.1 ppm
    else if (ppmRange < 3) ppmStep = 0.2;       // 0.2 ppm tick spacing.
    else if (ppmRange < 6) ppmStep = 0.5;       // Vue large: 0.5 ppm
    else ppmStep = 1;                            // 1 ppm tick spacing.
    
    // Align ticks to multiples of the step size.
    const startPpm = Math.floor(range.min / ppmStep) * ppmStep;
    for (let ppm = startPpm; ppm <= range.max; ppm += ppmStep) {
      if (ppm < range.min) continue;
      const x = ppmToX(ppm);
      ctx.beginPath();
      ctx.moveTo(x, PAD.top);
      ctx.lineTo(x, PAD.top + plotH);
      ctx.stroke();
    }
    
    // Draw spectrum line
    ctx.strokeStyle = "#1e40af";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    let started = false;
    for (const pt of visibleData) {
      const x = ppmToX(pt.ppm);
      const y = intToY(pt.int);
      if (!started) {
        ctx.moveTo(x, y);
        started = true;
      } else {
        ctx.lineTo(x, y);
      }
    }
    ctx.stroke();
    
    // Fill under curve (light blue)
    ctx.fillStyle = "rgba(59, 130, 246, 0.1)";
    ctx.beginPath();
    ctx.moveTo(ppmToX(range.max), intToY(0));
    for (const pt of visibleData) {
      ctx.lineTo(ppmToX(pt.ppm), intToY(pt.int));
    }
    ctx.lineTo(ppmToX(range.min), intToY(0));
    ctx.closePath();
    ctx.fill();
    
    // Draw theory markers (blue sticks)
    if (theoryShifts && theoryShifts.length > 0) {
      ctx.strokeStyle = "#3b82f6";
      ctx.lineWidth = 2;
      for (const ppm of theoryShifts) {
        if (ppm >= range.min && ppm <= range.max) {
          const x = ppmToX(ppm);
          ctx.beginPath();
          ctx.moveTo(x, PAD.top + plotH);
          ctx.lineTo(x, PAD.top + plotH - 30);
          ctx.stroke();
        }
      }
    }
    
    // Draw experimental markers (red sticks)
    if (expShifts && expShifts.length > 0) {
      ctx.strokeStyle = "#dc2626";
      ctx.lineWidth = 2;
      for (const ppm of expShifts) {
        if (ppm >= range.min && ppm <= range.max) {
          const x = ppmToX(ppm);
          ctx.beginPath();
          ctx.moveTo(x, PAD.top + plotH);
          ctx.lineTo(x, PAD.top + plotH - 20);
          ctx.stroke();
        }
      }
    }
    
    // Draw axes
    ctx.strokeStyle = "#374151";
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(PAD.left, PAD.top);
    ctx.lineTo(PAD.left, PAD.top + plotH);
    ctx.lineTo(PAD.left + plotW, PAD.top + plotH);
    ctx.stroke();
    
    // Draw ppm labels - adapt precision to zoom level
    ctx.fillStyle = "#374151";
    ctx.font = "11px system-ui, sans-serif";
    ctx.textAlign = "center";
    // Match label precision to the zoom level.
    const precision = ppmRange < 0.01 ? 4 : ppmRange < 0.1 ? 3 : ppmRange < 1 ? 2 : 1;
    const startPpmLabel = Math.floor(range.min / ppmStep) * ppmStep;
    for (let ppm = startPpmLabel; ppm <= range.max; ppm += ppmStep) {
      if (ppm < range.min) continue;
      const x = ppmToX(ppm);
      ctx.fillText(ppm.toFixed(precision), x, PAD.top + plotH + 15);
    }
    
    // X-axis label
    ctx.fillStyle = "#64748b";
    ctx.font = "12px system-ui, sans-serif";
    ctx.fillText("Chemical Shift (ppm)", PAD.left + plotW / 2, H - 8);

    // Draw zoom selection box if dragging
    if (isDragging && dragStart) {
      const currentPos = { x: dragStart.currentX, y: dragStart.currentY };
      ctx.strokeStyle = "rgba(59, 130, 246, 0.8)";
      ctx.fillStyle = "rgba(59, 130, 246, 0.1)";
      ctx.lineWidth = 2;
      ctx.setLineDash([5, 5]);

      const rectX = Math.min(dragStart.x, currentPos.x);
      const rectY = PAD.top;
      const rectW = Math.abs(currentPos.x - dragStart.x);
      const rectH = plotH;

      ctx.fillRect(rectX, rectY, rectW, rectH);
      ctx.strokeRect(rectX, rectY, rectW, rectH);
      ctx.setLineDash([]);
    }

    // Draw ONE global rectangle around all annotated signals with match tolerance
    if (expShifts && expShifts.length > 0) {
      const minSignal = Math.min(...expShifts);
      const maxSignal = Math.max(...expShifts);
      
      const ppmMin = minSignal - taskTolerancePpm;
      const ppmMax = maxSignal + taskTolerancePpm;

      const xMin = ppmToX(ppmMax); // Reversed for NMR
      const xMax = ppmToX(ppmMin);
      const rectX = xMin;
      const rectW = xMax - xMin;
      const rectY = PAD.top;
      const rectH = plotH;

      ctx.strokeStyle = "#f59e0b";
      ctx.fillStyle = "rgba(245, 158, 11, 0.08)";
      ctx.lineWidth = 2;
      ctx.setLineDash([5, 3]);
      ctx.fillRect(rectX, rectY, rectW, rectH);
      ctx.strokeRect(rectX, rectY, rectW, rectH);
      ctx.setLineDash([]);
    }
    
  }, [spectrumData, theoryShifts, expShifts, zoomRange, molId, verticalScaleSim, isDragging, dragStart, taskTolerancePpm]);

  // Draw experimental spectrum on canvas with zoom rectangles
  useEffect(() => {
    if (!expSpectrumData || !canvasExpRef.current || !spectrumData) return;

    const canvas = canvasExpRef.current;
    const ctx = canvas.getContext("2d");
    const DPR = window.devicePixelRatio || 1;

    const container = canvas.parentElement;
    const W = container.clientWidth || 800;
    const H = container.clientHeight || 400;

    canvas.width = W * DPR;
    canvas.height = H * DPR;
    canvas.style.width = `${W}px`;
    canvas.style.height = `${H}px`;
    ctx.scale(DPR, DPR);

    const PAD = { top: 20, right: 30, bottom: 40, left: 50 };
    const plotW = W - PAD.left - PAD.right;
    const plotH = H - PAD.top - PAD.bottom;

    // Clear
    ctx.fillStyle = "#fff";
    ctx.fillRect(0, 0, W, H);

    // Use same PPM range as simulated spectrum for alignment
    const defaultRange = spectrumData.ppm_range
      ? { min: spectrumData.ppm_range[0], max: spectrumData.ppm_range[1] }
      : { min: Math.min(...spectrumData.ppm), max: Math.max(...spectrumData.ppm) };
    const range = zoomRange || defaultRange;
    const ppmToX = (ppm) => PAD.left + ((range.max - ppm) / (range.max - range.min)) * plotW;

    // Filter experimental data in range (with offset applied)
    const ppmArr = expSpectrumData.ppm;
    const intArr = expSpectrumData.intensity;

    let maxInt = 0;
    const visibleData = [];
    for (let i = 0; i < ppmArr.length; i++) {
      // Apply ppm offset to shift experimental spectrum
      const ppm = ppmArr[i] + ppmOffsetExp;
      if (ppm >= range.min && ppm <= range.max) {
        visibleData.push({ ppm, int: intArr[i] });
        if (Math.abs(intArr[i]) > maxInt) maxInt = Math.abs(intArr[i]);
      }
    }
    if (maxInt === 0) maxInt = 1;

    // Apply vertical scale for experimental spectrum (independent)
    const intToY = (int) => PAD.top + plotH - (int / maxInt) * plotH * 0.9 * verticalScaleExp;

    // Draw grid with adaptive step
    ctx.strokeStyle = "#e5e7eb";
    ctx.lineWidth = 1;

    const ppmRange = range.max - range.min;
    let ppmStep;
    if (ppmRange < 0.005) ppmStep = 0.0005;     // 0.0005 ppm tick spacing.
    else if (ppmRange < 0.01) ppmStep = 0.001;  // 0.001 ppm tick spacing.
    else if (ppmRange < 0.05) ppmStep = 0.005;  // 0.005 ppm tick spacing.
    else if (ppmRange < 0.1) ppmStep = 0.01;    // 0.01 ppm tick spacing.
    else if (ppmRange < 0.5) ppmStep = 0.05;    // 0.05 ppm tick spacing.
    else if (ppmRange < 1) ppmStep = 0.1;       // Zoom moyen: 0.1 ppm
    else if (ppmRange < 3) ppmStep = 0.2;       // 0.2 ppm tick spacing.
    else if (ppmRange < 6) ppmStep = 0.5;       // Vue large: 0.5 ppm
    else ppmStep = 1;                            // 1 ppm tick spacing.
    
    const startPpm = Math.floor(range.min / ppmStep) * ppmStep;
    for (let ppm = startPpm; ppm <= range.max; ppm += ppmStep) {
      if (ppm < range.min) continue;
      const x = ppmToX(ppm);
      ctx.beginPath();
      ctx.moveTo(x, PAD.top);
      ctx.lineTo(x, PAD.top + plotH);
      ctx.stroke();
    }

    // Draw experimental spectrum (red)
    ctx.strokeStyle = "#dc2626";
    ctx.lineWidth = 1.5;
    ctx.beginPath();
    let started = false;
    for (const pt of visibleData) {
      const x = ppmToX(pt.ppm);
      const y = intToY(pt.int);
      if (!started) {
        ctx.moveTo(x, y);
        started = true;
      } else {
        ctx.lineTo(x, y);
      }
    }
    ctx.stroke();

    // Draw ONE global rectangle around all annotated signals with match tolerance  
    if (expShifts && expShifts.length > 0) {
      const minSignal = Math.min(...expShifts);
      const maxSignal = Math.max(...expShifts);
      
      const ppmMin = minSignal - taskTolerancePpm;
      const ppmMax = maxSignal + taskTolerancePpm;

      const xMin = ppmToX(ppmMax); // Reversed for NMR
      const xMax = ppmToX(ppmMin);
      const rectX = xMin;
      const rectW = xMax - xMin;
      const rectY = PAD.top;
      const rectH = plotH;

      ctx.strokeStyle = "#f59e0b";
      ctx.fillStyle = "rgba(245, 158, 11, 0.08)";
      ctx.lineWidth = 2;
      ctx.setLineDash([5, 3]);
      ctx.fillRect(rectX, rectY, rectW, rectH);
      ctx.strokeRect(rectX, rectY, rectW, rectH);
      ctx.setLineDash([]);
    }

    // X-axis labels - show ORIGINAL ppm values (without offset) for experimental spectrum
    // This way user sees the real experimental ppm values
    ctx.fillStyle = "#64748b";
    ctx.font = "11px system-ui, sans-serif";
    ctx.textAlign = "center";
    // Match label precision to the zoom level.
    const precision = ppmRange < 0.01 ? 4 : ppmRange < 0.1 ? 3 : ppmRange < 1 ? 2 : 1;
    const startPpmLabel = Math.floor(range.min / ppmStep) * ppmStep;
    for (let ppm = startPpmLabel; ppm <= range.max; ppm += ppmStep) {
      if (ppm < range.min) continue;
      const x = ppmToX(ppm);
      // Display the original experimental ppm (subtract the offset to show what ppm the data really is)
      const originalPpm = ppm - ppmOffsetExp;
      ctx.fillText(originalPpm.toFixed(precision), x, H - 20);
    }

    ctx.fillText("Chemical Shift (ppm)", PAD.left + plotW / 2, H - 8);

  }, [expSpectrumData, spectrumData, expShifts, zoomRange, taskTolerancePpm, verticalScaleExp, ppmOffsetExp]);

  // Handle mouse wheel for vertical scale (intensity)
  const handleWheel = useCallback((e) => {
    e.preventDefault();
    if (!spectrumData) return;

    // Ctrl/Cmd + wheel = horizontal zoom, otherwise = vertical scale
    if (e.ctrlKey || e.metaKey) {
      // Horizontal zoom (PPM axis)
      // Use the canvas that received the event (works for both simulated and experimental)
      const canvas = e.currentTarget;
      const rect = canvas.getBoundingClientRect();
      const mouseX = e.clientX - rect.left;

      const PAD = { left: 50, right: 30 };
      const plotW = rect.width - PAD.left - PAD.right;

      const defaultRange = spectrumData.ppm_range
        ? { min: spectrumData.ppm_range[0], max: spectrumData.ppm_range[1] }
        : { min: Math.min(...spectrumData.ppm), max: Math.max(...spectrumData.ppm) };
      const range = zoomRange || defaultRange;
      const centerPpm = range.max - ((mouseX - PAD.left) / plotW) * (range.max - range.min);

      const zoomFactor = e.deltaY < 0 ? 0.8 : 1.25;
      const newWidth = (range.max - range.min) * zoomFactor;

      // Allow zooming down to 0.001 ppm, equivalent to 1 Hz at 1000 MHz.
      if (newWidth < 0.001 || newWidth > 15) return;

      let newMin = centerPpm - newWidth / 2;
      let newMax = centerPpm + newWidth / 2;

      if (newMin < -1) { newMin = -1; newMax = newMin + newWidth; }
      if (newMax > 15) { newMax = 15; newMin = 15 - newWidth; }

      setZoomRange({ min: newMin, max: newMax });
    } else {
      // Vertical scale (intensity axis) - independent for each spectrum
      const scaleDelta = e.deltaY < 0 ? 1.1 : 0.9;
      const canvas = e.currentTarget;
      const isExpCanvas = canvas === canvasExpRef.current;
      if (isExpCanvas) {
        const newScale = Math.max(0.1, Math.min(10, verticalScaleExp * scaleDelta));
        setVerticalScaleExp(newScale);
      } else {
        const newScale = Math.max(0.1, Math.min(10, verticalScaleSim * scaleDelta));
        setVerticalScaleSim(newScale);
      }
    }
  }, [spectrumData, zoomRange, verticalScaleSim, verticalScaleExp]);

  // Mouse handlers for drag-to-zoom selection and horizontal panning
  const handleMouseDown = useCallback((e) => {
    if (!spectrumData) return;
    // Use the canvas that received the event (works for both simulated and experimental)
    const canvas = e.currentTarget;
    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;
    const isExpCanvas = canvas === canvasExpRef.current;

    const PAD = { top: 20, left: 50, bottom: 40, right: 30 };

    // Only start drag in plot area
    if (x >= PAD.left && x <= rect.width - PAD.right &&
        y >= PAD.top && y <= rect.height - PAD.bottom) {
      // Shift+drag on experimental canvas = horizontal panning
      const isPanning = e.shiftKey && isExpCanvas;
      setIsDragging(true);
      setDragStart({ x, y, currentX: x, currentY: y, canvas, isPanning, startOffset: ppmOffsetExp });
    }
  }, [spectrumData, ppmOffsetExp]);

  const handleMouseMove = useCallback((e) => {
    if (!isDragging || !dragStart) return;
    const canvas = dragStart.canvas;
    if (!canvas) return;
    const rect = canvas.getBoundingClientRect();
    const x = e.clientX - rect.left;
    const y = e.clientY - rect.top;

    // If panning the experimental spectrum
    if (dragStart.isPanning) {
      const PAD = { left: 50, right: 30 };
      const plotW = rect.width - PAD.left - PAD.right;
      
      const defaultRange = spectrumData.ppm_range
        ? { min: spectrumData.ppm_range[0], max: spectrumData.ppm_range[1] }
        : { min: Math.min(...spectrumData.ppm), max: Math.max(...spectrumData.ppm) };
      const range = zoomRange || defaultRange;
      
      // Calculate ppm delta from pixel movement
      const ppmPerPixel = (range.max - range.min) / plotW;
      const deltaX = x - dragStart.x;
      // Moving right (positive deltaX) = shift spectrum to lower ppm (negative offset)
      const deltaPpm = -deltaX * ppmPerPixel;
      
      setPpmOffsetExp(dragStart.startOffset + deltaPpm);
    }

    setDragStart({ ...dragStart, currentX: x, currentY: y });
  }, [isDragging, dragStart, spectrumData, zoomRange]);

  const handleMouseUp = useCallback(() => {
    if (!isDragging || !dragStart) return;

    // If it was a panning operation, just end it (offset already applied in handleMouseMove)
    if (dragStart.isPanning) {
      setIsDragging(false);
      setDragStart(null);
      return;
    }

    // Use the canvas that was stored in dragStart (works for both simulated and experimental)
    const canvas = dragStart.canvas;
    if (!canvas) {
      setIsDragging(false);
      setDragStart(null);
      return;
    }
    const rect = canvas.getBoundingClientRect();
    const PAD = { left: 50, right: 30 };
    const plotW = rect.width - PAD.left - PAD.right;

    const defaultRange = spectrumData.ppm_range
      ? { min: spectrumData.ppm_range[0], max: spectrumData.ppm_range[1] }
      : { min: Math.min(...spectrumData.ppm), max: Math.max(...spectrumData.ppm) };
    const range = zoomRange || defaultRange;

    // Convert pixel positions to PPM
    const xToPpm = (x) => range.max - ((x - PAD.left) / plotW) * (range.max - range.min);

    const ppm1 = xToPpm(dragStart.x);
    const ppm2 = xToPpm(dragStart.currentX);

    const newMin = Math.min(ppm1, ppm2);
    const newMax = Math.max(ppm1, ppm2);

    // Only zoom if selection is wide enough (0.001 ppm = ~0.7 Hz at 700 MHz)
    if (Math.abs(newMax - newMin) > 0.001) {
      setZoomRange({ min: newMin, max: newMax });
    }

    setIsDragging(false);
    setDragStart(null);
  }, [isDragging, dragStart, zoomRange, spectrumData]);

  // Reset zoom
  const resetZoom = () => {
    setZoomRange(null);
    setVerticalScaleSim(1.0);
    setVerticalScaleExp(1.0);
    setPpmOffsetExp(0.0);
  };

  if (!isOpen) return null;

  return (
    <>
      <style>{spinnerKeyframes}</style>
      <div style={MODAL_STYLES.overlay} onClick={onClose}>
        <div style={MODAL_STYLES.modal} onClick={(e) => e.stopPropagation()}>
          <div style={MODAL_STYLES.header}>
            <span style={MODAL_STYLES.title}>
              ¹H NMR Spectrum Detail - {molId}
            </span>
            <button style={MODAL_STYLES.closeBtn} onClick={onClose}>×</button>
          </div>
          
          <div style={MODAL_STYLES.controls}>
            <div style={MODAL_STYLES.controlGroup}>
              <span style={MODAL_STYLES.label}>Line Shape:</span>
              <select
                style={MODAL_STYLES.select}
                value={lineShape}
                onChange={(e) => setLineShape(e.target.value)}
              >
                <option value="lorentzian">Lorentzian</option>
                <option value="gaussian">Gaussian</option>
                <option value="voigt">Voigt</option>
              </select>
            </div>

            <div style={MODAL_STYLES.controlGroup}>
              <span style={MODAL_STYLES.label}>Line Width:</span>
              <input
                type="range"
                style={{ width: 120, marginRight: 8 }}
                value={lineWidthHz}
                onChange={(e) => setLineWidthHz(Number(e.target.value))}
                min={0.1}
                max={5}
                step={0.1}
              />
              <input
                type="number"
                style={{ ...MODAL_STYLES.input, width: 50 }}
                value={lineWidthHz}
                onChange={(e) => {
                  // Update value immediately for display
                  setLineWidthHzInternal(Number(e.target.value));
                }}
                onBlur={(e) => {
                  // Only trigger simulation when user finishes typing (blur)
                  const newValue = Number(e.target.value);
                  if (newValue >= 0.1 && newValue <= 10) {
                    setLineWidthHz(newValue);
                  } else {
                    // Reset to valid value if out of range
                    setLineWidthHzInternal(lineWidthHz);
                  }
                }}
                min={0.1}
                max={10}
                step={0.1}
              />
              <span style={{ fontSize: 11, color: "#6b7280" }}>Hz</span>
              <button
                onClick={resetLineWidthToAuto}
                style={{
                  marginLeft: 8,
                  padding: "3px 8px",
                  fontSize: 11,
                  borderRadius: 4,
                  border: "1px solid #d1d5db",
                  background: isLineWidthAuto ? "#e0f2fe" : "#fff",
                  color: isLineWidthAuto ? "#0369a1" : "#64748b",
                  cursor: "pointer",
                  fontWeight: isLineWidthAuto ? 600 : 400,
                }}
                title={`Auto: ${getDefaultLineWidth(freqMhz).toFixed(2)} Hz at ${freqMhz} MHz`}
              >
                Auto
              </button>
            </div>

            <div style={{ display: 'flex', gap: 12, fontSize: 11, color: '#6b7280', alignItems: 'center', marginLeft: 12, paddingLeft: 12, borderLeft: '1px solid #e5e7eb' }}>
              <span><strong>Conditions:</strong> {freqMhz} MHz, {solvent}, {temperature}K</span>
            </div>

            <button
              onClick={fetchSpectrum}
              disabled={loading}
              style={{
                padding: "6px 12px",
                backgroundColor: "#3b82f6",
                color: "#fff",
                border: "none",
                borderRadius: 4,
                cursor: loading ? "not-allowed" : "pointer",
                fontSize: 12,
                fontWeight: 500,
              }}
            >
              {loading ? "Simulating..." : "Refresh"}
            </button>
            
            <div style={MODAL_STYLES.controlGroup}>
              <span style={MODAL_STYLES.label}>Vertical Scale:</span>
              <input
                type="range"
                style={{ width: 100, marginRight: 8 }}
                value={verticalScaleSim}
                onChange={(e) => { setVerticalScaleSim(Number(e.target.value)); setVerticalScaleExp(Number(e.target.value)); }}
                min={0.1}
                max={10}
                step={0.1}
              />
              <input
                type="number"
                style={{ ...MODAL_STYLES.input, width: 50 }}
                value={verticalScaleSim.toFixed(1)}
                onChange={(e) => { setVerticalScaleSim(Number(e.target.value)); setVerticalScaleExp(Number(e.target.value)); }}
                min={0.1}
                max={10}
                step={0.1}
              />
              <span style={{ fontSize: 11, color: "#6b7280" }}>×</span>
            </div>

            <div style={{ display: "flex", gap: 8, alignItems: "center", flexWrap: "wrap" }}>
              <button
                onClick={resetZoom}
                disabled={!zoomRange && verticalScaleSim === 1.0 && verticalScaleExp === 1.0 && ppmOffsetExp === 0.0}
                style={{
                  padding: "6px 12px",
                  backgroundColor: zoomRange || verticalScaleSim !== 1.0 || verticalScaleExp !== 1.0 || ppmOffsetExp !== 0.0 ? "#f1f5f9" : "#f8fafc",
                  color: "#374151",
                  border: "1px solid #d1d5db",
                  borderRadius: 4,
                  cursor: zoomRange || verticalScaleSim !== 1.0 || verticalScaleExp !== 1.0 || ppmOffsetExp !== 0.0 ? "pointer" : "not-allowed",
                  fontSize: 12,
                  opacity: zoomRange || verticalScaleSim !== 1.0 || verticalScaleExp !== 1.0 || ppmOffsetExp !== 0.0 ? 1 : 0.5,
                }}
              >
                Reset View
              </button>
              <span style={{ fontSize: 10, color: "#64748b", fontStyle: "italic" }}>
                Drag to zoom | Scroll to scale | Ctrl+Scroll for H-zoom | Shift+Drag on Exp. to align
              </span>
              <button
                onClick={onClose}
                style={{
                  padding: "6px 12px",
                  backgroundColor: "#3b82f6",
                  color: "#fff",
                  border: "none",
                  borderRadius: 4,
                  cursor: "pointer",
                  fontSize: 12,
                  fontWeight: 500,
                }}
              >
                View Detail
              </button>
            </div>
          </div>

          <div style={MODAL_STYLES.canvasContainer}>
            {loadingMultiplicity && (
              <div style={MODAL_STYLES.loading}>
                <div style={MODAL_STYLES.spinner}></div>
                <span>Loading multiplicity data...</span>
              </div>
            )}
            {!loadingMultiplicity && loading && (
              <div style={MODAL_STYLES.loading}>
                <div style={MODAL_STYLES.spinner}></div>
                <span>Simulating spectrum...</span>
              </div>
            )}
            
            {error && (
              <div style={MODAL_STYLES.error}>
                <p>{error}</p>
                <button
                  onClick={fetchSpectrum}
                  style={{
                    marginTop: 12,
                    padding: "8px 16px",
                    backgroundColor: "#3b82f6",
                    color: "#fff",
                    border: "none",
                    borderRadius: 4,
                    cursor: "pointer",
                  }}
                >
                  Retry
                </button>
              </div>
            )}
            
            {!loading && !error && (
              <div style={{ display: "flex", flexDirection: "column", height: "100%", gap: 8 }}>
                {/* Simulated spectrum (top). */}
                <div style={{ flex: 1, position: "relative", border: "1px solid #e5e7eb", borderRadius: 8, minHeight: 0 }}>
                  <div style={{
                    position: "absolute",
                    top: 4,
                    left: 8,
                    fontSize: 11,
                    fontWeight: 600,
                    color: "#3b82f6",
                    zIndex: 10
                  }}>
                    Simulated Spectrum
                  </div>
                  <canvas
                    ref={canvasSimRef}
                    onWheel={handleWheel}
                    onMouseDown={handleMouseDown}
                    onMouseMove={handleMouseMove}
                    onMouseUp={handleMouseUp}
                    onMouseLeave={handleMouseUp}
                    style={{
                      display: "block",
                      width: "100%",
                      height: "100%",
                      cursor: isDragging ? "crosshair" : "default"
                    }}
                  />
                </div>

                {/* Experimental spectrum (bottom). */}
                <div style={{ flex: 1, position: "relative", border: "1px solid #e5e7eb", borderRadius: 8, minHeight: 0 }}>
                  <div style={{
                    position: "absolute",
                    top: 4,
                    left: 8,
                    fontSize: 11,
                    fontWeight: 600,
                    color: "#ef4444",
                    zIndex: 10
                  }}>
                    Experimental Spectrum
                    {ppmOffsetExp !== 0 && (
                      <span style={{ marginLeft: 8, fontSize: 10, color: "#f97316", fontWeight: 500 }}>
                        (shifted {ppmOffsetExp > 0 ? '+' : ''}{ppmOffsetExp.toFixed(2)} ppm)
                      </span>
                    )}
                  </div>
                  {loadingExpSpectrum && (
                    <div style={{
                      position: "absolute",
                      top: "50%",
                      left: "50%",
                      transform: "translate(-50%, -50%)",
                      textAlign: "center",
                      color: "#64748b"
                    }}>
                      Loading experimental data...
                    </div>
                  )}
                  {expSpectrumError && (
                    <div style={{
                      position: "absolute",
                      top: "50%",
                      left: "50%",
                      transform: "translate(-50%, -50%)",
                      textAlign: "center",
                      color: "#ef4444",
                      fontSize: 12
                    }}>
                      {expSpectrumError}
                    </div>
                  )}
                  {!loadingExpSpectrum && !expSpectrumError && expSpectrumData && (
                    <canvas
                      ref={canvasExpRef}
                      onWheel={handleWheel}
                      onMouseDown={handleMouseDown}
                      onMouseMove={handleMouseMove}
                      onMouseUp={handleMouseUp}
                      onMouseLeave={handleMouseUp}
                      style={{
                        display: "block",
                        width: "100%",
                        height: "100%",
                        cursor: isDragging ? "crosshair" : "default"
                      }}
                    />
                  )}
                  {!loadingExpSpectrum && !expSpectrumError && !expSpectrumData && (
                    <div style={{
                      position: "absolute",
                      top: "50%",
                      left: "50%",
                      transform: "translate(-50%, -50%)",
                      textAlign: "center",
                      color: "#94a3b8",
                      fontSize: 12
                    }}>
                      No experimental data available
                    </div>
                  )}
                </div>
              </div>
            )}
          </div>
          
          {spectrumData?.peaks && spectrumData.peaks.length > 0 && (
            <div style={MODAL_STYLES.peaksTable}>
              <table style={MODAL_STYLES.table}>
                <thead>
                  <tr>
                    <th style={MODAL_STYLES.th}>#</th>
                    <th style={MODAL_STYLES.th}>δ (ppm)</th>
                    <th style={MODAL_STYLES.th}>Multiplicity</th>
                    <th style={MODAL_STYLES.th}>J (Hz)</th>
                  </tr>
                </thead>
                <tbody>
                  {spectrumData.peaks.map((peak, i) => (
                    <tr key={i}>
                      <td style={MODAL_STYLES.td}>{i + 1}</td>
                      <td style={MODAL_STYLES.td}>{peak.ppm.toFixed(2)}</td>
                      <td style={MODAL_STYLES.td}>{peak.multiplicity || peak.label}</td>
                      <td style={MODAL_STYLES.td}>
                        {peak.J_hz?.length ? peak.J_hz.map(j => j.toFixed(1)).join(", ") : "-"}
                      </td>
                    </tr>
                  ))}
                </tbody>
              </table>
            </div>
          )}
        </div>
      </div>
    </>
  );
}

export default SpectrumDetailModal;
