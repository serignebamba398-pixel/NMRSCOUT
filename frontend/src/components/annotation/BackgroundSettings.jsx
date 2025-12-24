import React from 'react';
import { ParameterInput } from '../common';
import { normalizeStaticUrl } from '../../utils';

/**
 * Bruker spectral background configuration (PNG + JSON)
 */
export const BackgroundSettings = ({ annotationState, setAnnotationState, allRunning }) => {
  const { activeType, bgPdataDirs, bgOptions, bgStatus } = annotationState;

  const BgWindowInputs = () => {
    const opts = bgOptions[activeType] || {};

    if (activeType === "1H" || activeType === "13C") {
      return (
        <div className="parameters-grid" style={{ marginTop: 8 }}>
          <ParameterInput
            label="ppm_min"
            value={opts.ppm_min}
            onChange={(e) =>
              setAnnotationState((prev) => ({
                ...prev,
                bgOptions: {
                  ...prev.bgOptions,
                  [activeType]: {
                    ...prev.bgOptions[activeType],
                    ppm_min: parseFloat(e.target.value),
                  },
                },
              }))
            }
            disabled={allRunning}
          />
          <ParameterInput
            label="ppm_max"
            value={opts.ppm_max}
            onChange={(e) =>
              setAnnotationState((prev) => ({
                ...prev,
                bgOptions: {
                  ...prev.bgOptions,
                  [activeType]: {
                    ...prev.bgOptions[activeType],
                    ppm_max: parseFloat(e.target.value),
                  },
                },
              }))
            }
            disabled={allRunning}
          />
        </div>
      );
    }

    return (
      <div className="parameters-grid" style={{ marginTop: 8 }}>
        <ParameterInput
          label="F2 ppm min"
          value={opts.f2_ppm_min}
          onChange={(e) =>
            setAnnotationState((prev) => ({
              ...prev,
              bgOptions: {
                ...prev.bgOptions,
                [activeType]: {
                  ...prev.bgOptions[activeType],
                  f2_ppm_min: parseFloat(e.target.value),
                },
              },
            }))
          }
          disabled={allRunning}
        />
        <ParameterInput
          label="F2 ppm max"
          value={opts.f2_ppm_max}
          onChange={(e) =>
            setAnnotationState((prev) => ({
              ...prev,
              bgOptions: {
                ...prev.bgOptions,
                [activeType]: {
                  ...prev.bgOptions[activeType],
                  f2_ppm_max: parseFloat(e.target.value),
                },
              },
            }))
          }
          disabled={allRunning}
        />
        <ParameterInput
          label="F1 ppm min"
          value={opts.f1_ppm_min}
          onChange={(e) =>
            setAnnotationState((prev) => ({
              ...prev,
              bgOptions: {
                ...prev.bgOptions,
                [activeType]: {
                  ...prev.bgOptions[activeType],
                  f1_ppm_min: parseFloat(e.target.value),
                },
              },
            }))
          }
          disabled={allRunning}
        />
        <ParameterInput
          label="F1 ppm max"
          value={opts.f1_ppm_max}
          onChange={(e) =>
            setAnnotationState((prev) => ({
              ...prev,
              bgOptions: {
                ...prev.bgOptions,
                [activeType]: {
                  ...prev.bgOptions[activeType],
                  f1_ppm_max: parseFloat(e.target.value),
                },
              },
            }))
          }
          disabled={allRunning}
        />
      </div>
    );
  };

  return (
    <div className="card">
      <h3>Bruker Backgrounds (PNG + JSON) - {activeType}</h3>
      <div className="form-row">
        <label className="form-label">
          pdata folder (e.g.: /path/to/.../pdata/1)
        </label>
        <input
          className="form-input"
          type="text"
          value={bgPdataDirs[activeType] || ""}
          onChange={(e) =>
            setAnnotationState((prev) => ({
              ...prev,
              bgPdataDirs: {
                ...prev.bgPdataDirs,
                [activeType]: e.target.value,
              },
            }))
          }
          disabled={allRunning}
          placeholder="Paste your path (e.g.: C:\\...\\pdata\\1 or /Users/.../pdata/1)"
        />
      </div>

      <BgWindowInputs />

      <div className="hint" style={{ marginTop: 8 }}>
        The background (PNG+JSON) will be <strong>generated automatically</strong> for the{" "}
        <em>selected</em> methods on <strong>Launch annotation</strong>.
      </div>

      {bgStatus[activeType] && (
        <div className="file-status-summary" style={{ marginTop: 8 }}>
          {bgStatus[activeType].ok ? (
            <div className="status-item success">
              Last background:{" "}
              <a
                href={normalizeStaticUrl(bgStatus[activeType].pngUrl)}
                target="_blank"
                rel="noreferrer"
              >
                PNG
              </a>{" "}
              {bgStatus[activeType].meta
                ? `- window: ${
                    activeType === "1H" || activeType === "13C"
                      ? `${bgStatus[activeType].meta.ppm_left}..${bgStatus[activeType].meta.ppm_right} ppm`
                      : `F2 ${bgStatus[activeType].meta.f2_left}..${bgStatus[activeType].meta.f2_right} / F1 ${bgStatus[activeType].meta.f1_top}..${bgStatus[activeType].meta.f1_bottom} ppm`
                  }`
                : ""}
            </div>
          ) : (
            <div className="status-item missing">
              Background export failed: {bgStatus[activeType].error}
            </div>
          )}
        </div>
      )}
    </div>
  );
};