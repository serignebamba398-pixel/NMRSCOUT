import React from 'react';
import { ParameterInput } from '../common';

/**
 * Inputs to define tolerance ranges to test
 */
export const ToleranceRangeInputs = ({ type, ranges, setAnnotationState, running }) => {
  const updateRange = (dimension, param, value) => {
    setAnnotationState((prev) => ({
      ...prev,
      toleranceRanges: {
        ...prev.toleranceRanges,
        [type]: {
          ...prev.toleranceRanges[type],
          [dimension]: {
            ...prev.toleranceRanges[type][dimension],
            [param]: parseFloat(value),
          },
        },
      },
    }));
  };

  const calculateTestCount = () => {
    if (type === "1H" || type === "13C") {
      const count =
        Math.ceil(
          (ranges.tol_correspondence.max - ranges.tol_correspondence.min) /
            ranges.tol_correspondence.step
        ) + 1;
      return count;
    } else {
      const countF1 =
        Math.ceil(
          (ranges.tol_correspondence_f1.max - ranges.tol_correspondence_f1.min) /
            ranges.tol_correspondence_f1.step
        ) + 1;
      const countF2 =
        Math.ceil(
          (ranges.tol_correspondence_f2.max - ranges.tol_correspondence_f2.min) /
            ranges.tol_correspondence_f2.step
        ) + 1;
      return countF1 * countF2;
    }
  };

  return (
    <div className="card">
      <h3>Tolerance Windows to Test - {type}</h3>

      {(type === "1H" || type === "13C") && (
        <>
          <h4>tol_correspondence (this value will be tested)</h4>
          <div className="parameters-grid">
            <ParameterInput
              label="Min (ppm)"
              value={ranges.tol_correspondence.min}
              onChange={(e) => updateRange("tol_correspondence", "min", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Max (ppm)"
              value={ranges.tol_correspondence.max}
              onChange={(e) => updateRange("tol_correspondence", "max", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Step (ppm)"
              value={ranges.tol_correspondence.step}
              onChange={(e) => updateRange("tol_correspondence", "step", e.target.value)}
              disabled={running}
            />
          </div>
          <div
            className="hint"
            style={{
              marginTop: 4,
              marginBottom: 12,
              fontSize: 12,
              color: "#475569",
            }}
          >
            <strong>tol_match</strong> will automatically use the same value as{" "}
            <strong>tol_correspondence</strong>
          </div>
        </>
      )}

      {type === "COSY" && (
        <>
          <h4>tol_correspondence_f1 (F1 axis)</h4>
          <div className="parameters-grid">
            <ParameterInput
              label="Min (ppm)"
              value={ranges.tol_correspondence_f1.min}
              onChange={(e) => updateRange("tol_correspondence_f1", "min", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Max (ppm)"
              value={ranges.tol_correspondence_f1.max}
              onChange={(e) => updateRange("tol_correspondence_f1", "max", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Step (ppm)"
              value={ranges.tol_correspondence_f1.step}
              onChange={(e) => updateRange("tol_correspondence_f1", "step", e.target.value)}
              disabled={running}
            />
          </div>
          <div
            className="hint"
            style={{
              marginTop: 4,
              marginBottom: 12,
              fontSize: 12,
              color: "#475569",
            }}
          >
            <strong>tol_match_f1</strong> will automatically use the same value as{" "}
            <strong>tol_correspondence_f1</strong>
          </div>

          <h4>tol_correspondence_f2 (F2 axis)</h4>
          <div className="parameters-grid">
            <ParameterInput
              label="Min (ppm)"
              value={ranges.tol_correspondence_f2.min}
              onChange={(e) => updateRange("tol_correspondence_f2", "min", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Max (ppm)"
              value={ranges.tol_correspondence_f2.max}
              onChange={(e) => updateRange("tol_correspondence_f2", "max", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Step (ppm)"
              value={ranges.tol_correspondence_f2.step}
              onChange={(e) => updateRange("tol_correspondence_f2", "step", e.target.value)}
              disabled={running}
            />
          </div>
          <div
            className="hint"
            style={{
              marginTop: 4,
              marginBottom: 12,
              fontSize: 12,
              color: "#475569",
            }}
          >
            <strong>tol_match_f2</strong> will automatically use the same value as{" "}
            <strong>tol_correspondence_f2</strong>
          </div>
        </>
      )}

      {type === "HSQC" && (
        <>
          <h4>tol_correspondence_f1 (F1 axis - ¹³C)</h4>
          <div className="parameters-grid">
            <ParameterInput
              label="Min (ppm)"
              value={ranges.tol_correspondence_f1.min}
              onChange={(e) => updateRange("tol_correspondence_f1", "min", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Max (ppm)"
              value={ranges.tol_correspondence_f1.max}
              onChange={(e) => updateRange("tol_correspondence_f1", "max", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Step (ppm)"
              value={ranges.tol_correspondence_f1.step}
              onChange={(e) => updateRange("tol_correspondence_f1", "step", e.target.value)}
              disabled={running}
            />
          </div>
          <div
            className="hint"
            style={{
              marginTop: 4,
              marginBottom: 12,
              fontSize: 12,
              color: "#475569",
            }}
          >
            <strong>tol_match_f1</strong> will automatically use the same value as{" "}
            <strong>tol_correspondence_f1</strong>
          </div>

          <h4>tol_correspondence_f2 (F2 axis - ¹H)</h4>
          <div className="parameters-grid">
            <ParameterInput
              label="Min (ppm)"
              value={ranges.tol_correspondence_f2.min}
              onChange={(e) => updateRange("tol_correspondence_f2", "min", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Max (ppm)"
              value={ranges.tol_correspondence_f2.max}
              onChange={(e) => updateRange("tol_correspondence_f2", "max", e.target.value)}
              disabled={running}
            />
            <ParameterInput
              label="Step (ppm)"
              value={ranges.tol_correspondence_f2.step}
              onChange={(e) => updateRange("tol_correspondence_f2", "step", e.target.value)}
              disabled={running}
            />
          </div>
        </>
      )}

      <div className="hint" style={{ marginTop: 12 }}>
        {type === "1H" || type === "13C" ? (
          <>
            Number of tests: <strong>{calculateTestCount()}</strong>
            <br />
            <span style={{ fontSize: 11, color: "#64748b" }}>
              (tol_match will automatically take the {calculateTestCount()} values of
              tol_correspondence)
            </span>
          </>
        ) : (
          <>
            Number of combinations: <strong>{calculateTestCount()}</strong>
            <br />
            <span style={{ fontSize: 11, color: "#64748b" }}>
              (tol_match_f1 and tol_match_f2 will automatically follow correspondence
              values)
            </span>
          </>
        )}
      </div>
    </div>
  );
};