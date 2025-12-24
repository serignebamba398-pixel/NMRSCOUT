import React, { useRef, useEffect, useState } from 'react';

/**
 * Display of tolerance optimization result charts
 */
export const ToleranceResultsCharts = ({ resultSet, onDelete }) => {
  const { id, type, results, label } = resultSet;
  
  const canvasRef1 = useRef(null);
  const canvasRef2 = useRef(null);
  const canvasRef3 = useRef(null);
  const wrapperRef = useRef(null);

  const [canvasSize, setCanvasSize] = useState({
    width: 700,
    height: 500,
  });

  // Responsive canvas sizing
  useEffect(() => {
    const updateSize = () => {
      if (!wrapperRef.current) return;
      const rect = wrapperRef.current.getBoundingClientRect();

      const isSingleColumn = window.innerWidth < 900;
      const columns = isSingleColumn ? 1 : 2;

      const columnWidth = rect.width / columns;
      const innerWidth = Math.max(320, columnWidth - 32);
      const innerHeight = Math.max(260, innerWidth * 0.7);

      setCanvasSize({
        width: Math.round(innerWidth),
        height: Math.round(innerHeight),
      });
    };

    updateSize();
    window.addEventListener("resize", updateSize);
    return () => window.removeEventListener("resize", updateSize);
  }, []);

  // Chart 1: Coverage distribution (stacked area)
  const drawCoverageChart = (canvas, data) => {
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const width = canvas.width;
    const height = canvas.height;
    const PAD = 80;
    const legendWidth = 120;
    const graphW = width - 2 * PAD - legendWidth;
    const graphH = height - 2 * PAD;

    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);

    // Title
    ctx.fillStyle = "#000000";
    ctx.font = "bold 14px sans-serif";
    ctx.fillText(
      "Molecules by coverage class (% of molecules) vs Tolerance",
      PAD,
      26
    );

    // Axes
    ctx.strokeStyle = "#333333";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(PAD, PAD);
    ctx.lineTo(PAD, PAD + graphH);
    ctx.lineTo(PAD + graphW, PAD + graphH);
    ctx.stroke();

    const resultsFiltered = data.results.filter((r) => !r.error);
    if (resultsFiltered.length === 0) return;

    const nbPoints = resultsFiltered.length;
    const colors = [
      "rgba(31, 119, 180, 0.7)",
      "rgba(255, 127, 14, 0.7)",
      "rgba(44, 160, 44, 0.7)",
      "rgba(214, 39, 40, 0.7)",
      "rgba(148, 103, 189, 0.7)",
    ];
    const binLabels = ["0–20%", "20–40%", "40–60%", "60–80%", "80–100%"];
    const maxY = 100;

    // Calcul cumuls
    const stacked = [];
    resultsFiltered.forEach((r) => {
      const totalMol = r.totalMolecules || 0;
      const dist = r.coverageDistribution || [];
      const binsPct = [];

      for (let b = 0; b < 5; b++) {
        const countBin = dist[b] || 0;
        const pct = totalMol > 0 ? (countBin / totalMol) * 100 : 0;
        binsPct.push(pct);
      }

      const cumBase = [];
      const cumTop = [];
      let acc = 0;
      for (let b = 0; b < 5; b++) {
        cumBase[b] = acc;
        acc += binsPct[b];
        cumTop[b] = acc;
      }

      stacked.push({ binsPct, cumBase, cumTop });
    });

    // Dessin aires empilées
    for (let bin = 0; bin < 5; bin++) {
      ctx.beginPath();

      for (let i = 0; i < nbPoints; i++) {
        const x =
          nbPoints === 1 ? PAD + graphW / 2 : PAD + (i / (nbPoints - 1)) * graphW;
        const topVal = stacked[i].cumTop[bin];
        const yTop = PAD + graphH - (topVal / maxY) * graphH;

        if (i === 0) ctx.moveTo(x, yTop);
        else ctx.lineTo(x, yTop);
      }

      for (let i = nbPoints - 1; i >= 0; i--) {
        const x =
          nbPoints === 1 ? PAD + graphW / 2 : PAD + (i / (nbPoints - 1)) * graphW;
        const baseVal = stacked[i].cumBase[bin];
        const yBase = PAD + graphH - (baseVal / maxY) * graphH;
        ctx.lineTo(x, yBase);
      }

      ctx.closePath();
      ctx.fillStyle = colors[bin];
      ctx.fill();
      ctx.strokeStyle = "rgba(0,0,0,0.2)";
      ctx.stroke();
    }

    // Légende
    binLabels.forEach((label, bin) => {
      const legendX = PAD + graphW + 20;
      const legendY = PAD + bin * 25;
      ctx.fillStyle = colors[bin];
      ctx.fillRect(legendX, legendY, 15, 10);
      ctx.fillStyle = "#000000";
      ctx.font = "11px sans-serif";
      ctx.fillText(label, legendX + 20, legendY + 9);
    });

    // Axis labels
    ctx.fillStyle = "#000000";
    ctx.font = "12px sans-serif";
    ctx.fillText("Tolerance (ppm)", PAD + graphW / 2 - 50, height - 10);

    ctx.save();
    ctx.translate(15, PAD + graphH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("% of molecules", 0, 0);
    ctx.restore();

    // Ticks
    ctx.font = "10px sans-serif";
    resultsFiltered.forEach((r, i) => {
      const x = nbPoints === 1 ? PAD + graphW / 2 : PAD + (i / (nbPoints - 1)) * graphW;

      let label;
      if (r.toleranceParams.tol_correspondence !== undefined) {
        label = r.toleranceParams.tol_correspondence.toFixed(2);
      } else {
        label = `${r.toleranceParams.tol_correspondence_f1.toFixed(
          2
        )},${r.toleranceParams.tol_correspondence_f2.toFixed(2)}`;
      }
      ctx.save();
      ctx.translate(x, PAD + graphH + 15);
      ctx.rotate(-Math.PI / 4);
      ctx.fillText(label, 0, 0);
      ctx.restore();
    });

    const ySteps = [0, 20, 40, 60, 80, 100];
    ySteps.forEach((val) => {
      const y = PAD + graphH - (val / maxY) * graphH;
      ctx.fillText(val.toString(), PAD - 28, y + 5);
    });
  };

  // Chart 2: Matched signals
  const drawMatchedSignalsChart = (canvas, data) => {
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const width = canvas.width;
    const height = canvas.height;
    const PAD = 80;
    const graphW = width - 2 * PAD;
    const graphH = height - 2 * PAD;

    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);

    const resultsFiltered = data.results.filter((r) => !r.error);
    if (resultsFiltered.length === 0) return;

    const usePct = resultsFiltered.some(
      (r) =>
        typeof r.signalsMatchedPct === "number" &&
        Number.isFinite(r.signalsMatchedPct)
    );

    const maxSignals = usePct
      ? 100
      : Math.max(...resultsFiltered.map((r) => r.totalSignalsMatched || 0)) || 1;

    ctx.fillStyle = "#000000";
    ctx.font = "bold 14px sans-serif";
    ctx.fillText(
      usePct
        ? "Experimental signals matched (% of total signals) vs Tolerance"
        : "Number of experimental signals matched vs Tolerance",
      PAD,
      26
    );

    ctx.strokeStyle = "#333333";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(PAD, PAD);
    ctx.lineTo(PAD, PAD + graphH);
    ctx.lineTo(PAD + graphW, PAD + graphH);
    ctx.stroke();

    ctx.strokeStyle = "#2ca02c";
    ctx.fillStyle = "#2ca02c";
    ctx.lineWidth = 3;
    ctx.beginPath();

    resultsFiltered.forEach((r, i) => {
      const x =
        resultsFiltered.length === 1
          ? PAD + graphW / 2
          : PAD + (i / (resultsFiltered.length - 1)) * graphW;

      const value = usePct ? r.signalsMatchedPct || 0 : r.totalSignalsMatched || 0;
      const y = PAD + graphH - (value / maxSignals) * graphH;

      if (i === 0) ctx.moveTo(x, y);
      else ctx.lineTo(x, y);

      ctx.save();
      ctx.beginPath();
      ctx.arc(x, y, 4, 0, Math.PI * 2);
      ctx.fill();
      ctx.restore();
    });

    ctx.stroke();

    ctx.fillStyle = "#000000";
    ctx.font = "12px sans-serif";
    ctx.fillText("Tolerance (ppm)", PAD + graphW / 2 - 50, height - 10);

    ctx.save();
    ctx.translate(15, PAD + graphH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText(
      usePct ? "Signals matched (% of exp. signals)" : "Signals matched",
      0,
      0
    );
    ctx.restore();

    ctx.font = "10px sans-serif";
    resultsFiltered.forEach((r, i) => {
      const x =
        resultsFiltered.length === 1
          ? PAD + graphW / 2
          : PAD + (i / (resultsFiltered.length - 1)) * graphW;

      let label;
      if (r.toleranceParams.tol_correspondence !== undefined) {
        label = r.toleranceParams.tol_correspondence.toFixed(2);
      } else {
        label = `${r.toleranceParams.tol_correspondence_f1.toFixed(
          2
        )},${r.toleranceParams.tol_correspondence_f2.toFixed(2)}`;
      }
      ctx.save();
      ctx.translate(x, PAD + graphH + 15);
      ctx.rotate(-Math.PI / 4);
      ctx.fillText(label, 0, 0);
      ctx.restore();
    });

    ctx.font = "10px sans-serif";
    if (usePct) {
      const steps = [0, 20, 40, 60, 80, 100];
      steps.forEach((val) => {
        const y = PAD + graphH - (val / 100) * graphH;
        ctx.fillText(`${val}`, PAD - 35, y + 5);
      });
    } else {
      const ySteps = 5;
      for (let i = 0; i <= ySteps; i++) {
        const y = PAD + graphH - (i / ySteps) * graphH;
        const val = Math.round((i / ySteps) * maxSignals);
        ctx.fillText(val, PAD - 40, y + 5);
      }
    }
  };

  // Chart 3: Signal overlap
  const drawSignalOverlapChart = (canvas, data) => {
    if (!canvas) return;
    const ctx = canvas.getContext("2d");
    const width = canvas.width;
    const height = canvas.height;
    const PAD = 80;
    const graphW = width - 2 * PAD;
    const graphH = height - 2 * PAD;

    ctx.clearRect(0, 0, width, height);
    ctx.fillStyle = "#ffffff";
    ctx.fillRect(0, 0, width, height);

    ctx.fillStyle = "#000000";
    ctx.font = "bold 16px sans-serif";
    ctx.fillText(
      "Experimental signals matched by tolerance",
      PAD,
      30
    );

    const resultsFiltered = data.results.filter((r) => !r.error);
    if (resultsFiltered.length === 0) return;

    // Calculer recouvrements
    const overlapInfo = resultsFiltered.map((r, idx) => {
      const currSignals = new Set(r.matchedSignals || []);
      const commonWithOthers = new Set();

      resultsFiltered.forEach((other, otherIdx) => {
        if (otherIdx !== idx) {
          const otherSignals = new Set(other.matchedSignals || []);
          currSignals.forEach((sig) => {
            if (otherSignals.has(sig)) {
              commonWithOthers.add(sig);
            }
          });
        }
      });

      return commonWithOthers;
    });

    const allSignals = resultsFiltered.flatMap((r) => r.matchedSignals || []);
    const minPPM = Math.min(...allSignals, 0);
    const maxPPM = Math.max(...allSignals, 10);

    ctx.strokeStyle = "#333333";
    ctx.lineWidth = 2;
    ctx.beginPath();
    ctx.moveTo(PAD, PAD);
    ctx.lineTo(PAD, PAD + graphH);
    ctx.lineTo(PAD + graphW, PAD + graphH);
    ctx.stroke();

    const colSpacing = graphW / Math.max(resultsFiltered.length, 1);
    const rangePPM = maxPPM - minPPM || 10;

    resultsFiltered.forEach((r, idx) => {
      const x = PAD + idx * colSpacing + colSpacing / 2;

      ctx.strokeStyle = "#e5e7eb";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(x, PAD);
      ctx.lineTo(x, PAD + graphH);
      ctx.stroke();

      let tolLabel = "";
      if (r.toleranceParams.tol_correspondence !== undefined) {
        tolLabel = r.toleranceParams.tol_correspondence.toFixed(2);
      } else if (data.type === "COSY") {
        tolLabel = `${r.toleranceParams.tol_correspondence_f1.toFixed(
          1
        )}/${r.toleranceParams.tol_correspondence_f2.toFixed(1)}`;
      } else if (data.type === "HSQC") {
        tolLabel = `${r.toleranceParams.tol_correspondence_f1.toFixed(
          1
        )}/${r.toleranceParams.tol_correspondence_f2.toFixed(1)}`;
      }

      ctx.fillStyle = "#000000";
      ctx.font = "12px sans-serif";
      ctx.save();
      ctx.translate(x - 5, PAD + graphH + 40);
      ctx.rotate(-Math.PI / 4);
      ctx.fillText(tolLabel, 0, 0);
      ctx.restore();

      const signals = r.matchedSignals || [];
      const commonWithOthers = overlapInfo[idx];

      signals.forEach((ppm) => {
        const y = PAD + graphH - ((ppm - minPPM) / rangePPM) * graphH;
        const isCommon = commonWithOthers.has(ppm);

        ctx.fillStyle = isCommon ? "#ef4444" : "#3b82f6";
        ctx.beginPath();
        ctx.arc(x, y, isCommon ? 3 : 2, 0, 2 * Math.PI);
        ctx.fill();
      });
    });

    // Ticks PPM
    ctx.fillStyle = "#000000";
    ctx.font = "12px sans-serif";
    const numTicks = 5;
    for (let i = 0; i <= numTicks; i++) {
      const ppm = minPPM + (i / numTicks) * rangePPM;
      const y = PAD + graphH - (i / numTicks) * graphH;

      ctx.strokeStyle = "#666";
      ctx.lineWidth = 1;
      ctx.beginPath();
      ctx.moveTo(PAD - 5, y);
      ctx.lineTo(PAD, y);
      ctx.stroke();

      ctx.fillText(ppm.toFixed(1), PAD - 40, y + 5);
    }

    ctx.fillStyle = "#000000";
    ctx.font = "14px sans-serif";
    ctx.fillText("Tolerances →", PAD + graphW / 2 - 50, height - 12);

    ctx.save();
    ctx.translate(15, PAD + graphH / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillText("PPM (δ) →", 0, 0);
    ctx.restore();

    // Legend
    ctx.fillStyle = "#3b82f6";
    ctx.fillRect(PAD, height - 40, 18, 4);
    ctx.fillStyle = "#000";
    ctx.font = "13px sans-serif";
    ctx.fillText("Signal unique to this tolerance", PAD + 25, height - 34);

    ctx.fillStyle = "#ef4444";
    ctx.fillRect(PAD + 230, height - 40, 18, 4);
    ctx.fillStyle = "#000";
    ctx.fillText("Signal common with other tolerances", PAD + 255, height - 34);
  };

  // Calculate overlap matrix
  const calculateOverlapMatrix = () => {
    const resultsFiltered = results.filter((r) => !r.error);
    if (resultsFiltered.length === 0) return [];

    const matrix = [];

    resultsFiltered.forEach((r1, i) => {
      const signals1 = new Set(r1.matchedSignals || []);
      const row = [];

      resultsFiltered.forEach((r2, j) => {
        const signals2 = new Set(r2.matchedSignals || []);

        if (i === j) {
          row.push(100);
        } else {
          const common = Array.from(signals1).filter((s) => signals2.has(s)).length;
          const percentage = signals1.size > 0 ? (common / signals1.size) * 100 : 0;
          row.push(percentage);
        }
      });

      matrix.push(row);
    });

    return matrix;
  };

  const overlapMatrix = calculateOverlapMatrix();

  useEffect(() => {
    if (!results || results.length === 0) return;
    if (!canvasRef1.current || !canvasRef2.current || !canvasRef3.current) return;

    const adaptedData = { type, results };
    drawCoverageChart(canvasRef1.current, adaptedData);
    drawMatchedSignalsChart(canvasRef2.current, adaptedData);
    drawSignalOverlapChart(canvasRef3.current, adaptedData);
  }, [results, type, canvasSize]);

  return (
    <div
      className="results-section"
      style={{
        marginBottom: 40,
        padding: 20,
        border: "2px solid #e5e7eb",
        borderRadius: 8,
        background: "#f9fafb",
      }}
    >
      <div
        style={{
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
          marginBottom: 20,
        }}
      >
        <h3 style={{ margin: 0 }}>{label}</h3>
        <button
          onClick={() => onDelete(id)}
          className="btn btn-danger btn-small"
          title="Delete this test result"
        >
          🗑️ Delete
        </button>
      </div>

      <div
        ref={wrapperRef}
        style={{
          display: "grid",
          gridTemplateColumns: "1fr 1fr",
          gap: 20,
          marginTop: 20,
        }}
      >
        <div>
          <canvas
            ref={canvasRef1}
            width={canvasSize.width}
            height={canvasSize.height}
            style={{
              border: "1px solid #ddd",
              borderRadius: 8,
              width: "100%",
              height: "auto",
            }}
          />
        </div>
        <div>
          <canvas
            ref={canvasRef2}
            width={canvasSize.width}
            height={canvasSize.height}
            style={{
              border: "1px solid #ddd",
              borderRadius: 8,
              width: "100%",
              height: "auto",
            }}
          />
        </div>
      </div>

      <div style={{ marginTop: 30 }}>
        <div
          style={{
            display: "grid",
            gridTemplateColumns: "55% 45%",
            gap: 20,
            alignItems: "end",
          }}
        >
          <div>
            <canvas
              ref={canvasRef3}
              width={canvasSize.width * 2}
              height={canvasSize.height}
              style={{
                border: "1px solid #ddd",
                borderRadius: 8,
                width: "100%",
                height: "auto",
                display: "block",
              }}
            />
          </div>

          {overlapMatrix.length > 0 && (
            <div
              style={{
                border: "1px solid #ddd",
                borderRadius: 8,
                padding: 15,
                background: "#fff",
                height: "100%",
                display: "flex",
                flexDirection: "column",
              }}
            >
              <h4
                style={{
                  marginBottom: 15,
                  marginTop: 0,
                  fontSize: 16,
                  fontWeight: "bold",
                }}
              >
                Matrice de recouvrement (%)
              </h4>
              <div
                style={{
                  flex: 1,
                  overflowX: "auto",
                  overflowY: "auto",
                  minHeight: 0,
                }}
              >
                <table
                  style={{
                    borderCollapse: "collapse",
                    fontSize: 10,
                    width: "100%",
                    tableLayout: "fixed",
                  }}
                >
                  <thead
                    style={{
                      position: "sticky",
                      top: 0,
                      background: "#fff",
                      zIndex: 1,
                    }}
                  >
                    <tr>
                      <th
                        style={{
                          border: "1px solid #ddd",
                          padding: "4px",
                          background: "#f3f4f6",
                          fontWeight: "bold",
                          fontSize: 9,
                          width: "50px",
                          minWidth: "50px",
                        }}
                      >
                        Tol.
                      </th>
                      {results
                        .filter((r) => !r.error)
                        .map((r, idx) => {
                          let tolLabel = "";
                          if (r.toleranceParams.tol_correspondence !== undefined) {
                            tolLabel = r.toleranceParams.tol_correspondence.toFixed(2);
                          } else if (type === "COSY") {
                            tolLabel = `${r.toleranceParams.tol_correspondence_f1.toFixed(
                              1
                            )}/${r.toleranceParams.tol_correspondence_f2.toFixed(1)}`;
                          } else if (type === "HSQC") {
                            tolLabel = `${r.toleranceParams.tol_correspondence_f1.toFixed(
                              1
                            )}/${r.toleranceParams.tol_correspondence_f2.toFixed(1)}`;
                          }
                          return (
                            <th
                              key={idx}
                              style={{
                                border: "1px solid #ddd",
                                padding: "4px",
                                background: "#f3f4f6",
                                fontWeight: "bold",
                                textAlign: "center",
                                fontSize: 8,
                                whiteSpace: "nowrap",
                              }}
                            >
                              {tolLabel}
                            </th>
                          );
                        })}
                    </tr>
                  </thead>
                  <tbody>
                    {results
                      .filter((r) => !r.error)
                      .map((r1, i) => {
                        let tolLabel1 = "";
                        if (r1.toleranceParams.tol_correspondence !== undefined) {
                          tolLabel1 = r1.toleranceParams.tol_correspondence.toFixed(2);
                        } else if (type === "COSY") {
                          tolLabel1 = `${r1.toleranceParams.tol_correspondence_f1.toFixed(
                            1
                          )}/${r1.toleranceParams.tol_correspondence_f2.toFixed(1)}`;
                        } else if (type === "HSQC") {
                          tolLabel1 = `${r1.toleranceParams.tol_correspondence_f1.toFixed(
                            1
                          )}/${r1.toleranceParams.tol_correspondence_f2.toFixed(1)}`;
                        }
                        return (
                          <tr key={i}>
                            <td
                              style={{
                                border: "1px solid #ddd",
                                padding: "4px",
                                background: "#f9fafb",
                                fontWeight: "bold",
                                fontSize: 8,
                                textAlign: "center",
                                whiteSpace: "nowrap",
                              }}
                            >
                              {tolLabel1}
                            </td>
                            {overlapMatrix[i].map((percentage, j) => {
                              let bgColor = "#ffffff";
                              if (i === j) {
                                bgColor = "#dbeafe";
                              } else if (percentage >= 80) {
                                bgColor = "#fecaca";
                              } else if (percentage >= 50) {
                                bgColor = "#fde68a";
                              } else if (percentage > 0) {
                                bgColor = "#d1fae5";
                              }

                              return (
                                <td
                                  key={j}
                                  style={{
                                    border: "1px solid #ddd",
                                    padding: "4px",
                                    textAlign: "center",
                                    background: bgColor,
                                    fontWeight: i === j ? "bold" : "normal",
                                    fontSize: 8,
                                    whiteSpace: "nowrap",
                                  }}
                                >
                                  {percentage.toFixed(1)}%
                                </td>
                              );
                            })}
                          </tr>
                        );
                      })}
                  </tbody>
                </table>
              </div>
            </div>
          )}
        </div>
      </div>
    </div>
  );
};