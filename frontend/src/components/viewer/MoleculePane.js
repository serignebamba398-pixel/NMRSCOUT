// src/components/viewer/MoleculePane.js
// Component for displaying a molecule in 2D (RDKit) and 3D (3Dmol)

import React, { useRef, useEffect, useState, useCallback, memo } from "react";
import { loadRDKitFromCDN, load3DmolFromCDN, isWebGLAvailable, isProblematicBrowser } from "../../utils/cdnLoaders";
import { extractFirstMolBlockPreserving, resolveMolForSource } from "../../utils/molParser";
import { debugLog, debugWarn, debugError } from "../../utils/debug";
import { API_BASE_URL } from "../../api";

/**
 * Molecule visualization panel (2D + 3D)
 */
function MoleculePane({
  molBlock = null,
  molUrl = null,
  filename = null,
  height2D = 240,
  height3D = 280,
  show3D = true,
  highlightAtoms = [],
  highlightColor = "#ff6b6b",
  onMolLoaded,
  style = {},
}) {
  const container2DRef = useRef(null);
  const container3DRef = useRef(null);
  const viewerRef = useRef(null);

  const [rdkit, setRdkit] = useState(null);
  const [$3Dmol, set$3Dmol] = useState(null);
  const [currentMol, setCurrentMol] = useState(null);
  const [loading2D, setLoading2D] = useState(false);
  const [loading3D, setLoading3D] = useState(false);
  const [error2D, setError2D] = useState(null);
  const [error3D, setError3D] = useState(null);
  const [webglOK, setWebglOK] = useState(true);

  // Load libraries
  useEffect(() => {
    loadRDKitFromCDN()
      .then((kit) => {
        setRdkit(kit);
        debugLog("[OK] RDKit loaded in MoleculePane");
      })
      .catch((e) => {
        debugError("[ERROR] RDKit load failed:", e);
        setError2D("RDKit not available");
      });

    if (show3D && !isProblematicBrowser()) {
      if (!isWebGLAvailable()) {
        setWebglOK(false);
        setError3D("WebGL not available");
        return;
      }
      load3DmolFromCDN()
        .then((mol) => {
          set$3Dmol(mol);
          debugLog("[OK] 3Dmol loaded in MoleculePane");
        })
        .catch((e) => {
          debugError("[ERROR] 3Dmol load failed:", e);
          setError3D("3Dmol not available");
        });
    }
  }, [show3D]);

  // Resolve molBlock
  useEffect(() => {
    let cancelled = false;

    async function resolve() {
      if (molBlock) {
        setCurrentMol(extractFirstMolBlockPreserving(molBlock));
        return;
      }
      if (molUrl || filename) {
        const source = molUrl || filename;
        setLoading2D(true);
        try {
          const mol = await resolveMolForSource(source, API_BASE_URL);
          if (!cancelled) {
            setCurrentMol(mol);
            onMolLoaded && onMolLoaded(mol);
          }
        } catch (e) {
          if (!cancelled) {
            setError2D(e.message);
          }
        } finally {
          if (!cancelled) setLoading2D(false);
        }
      }
    }

    resolve();
    return () => {
      cancelled = true;
    };
  }, [molBlock, molUrl, filename, onMolLoaded]);

  // 2D rendering with RDKit
  useEffect(() => {
    if (!rdkit || !currentMol || !container2DRef.current) return;

    const container = container2DRef.current;
    container.innerHTML = "";
    setError2D(null);

    try {
      const mol = rdkit.get_mol(currentMol);
      if (!mol) {
        setError2D("Invalid molecule");
        return;
      }

      const opts = {
        width: container.clientWidth || 320,
        height: height2D,
        backgroundColour: [1, 1, 1, 1],
        highlightColour: hexToRgbFloat(highlightColor),
        highlightAtoms: highlightAtoms.length ? highlightAtoms : undefined,
        addAtomIndices: false,
        addBondIndices: false,
        addStereoAnnotation: true,
        legend: "",
      };

      const svg = mol.get_svg_with_highlights(JSON.stringify(opts));
      container.innerHTML = svg;

      mol.delete();
      debugLog("[OK] 2D mol rendered");
    } catch (e) {
      debugError("[ERROR] RDKit render error:", e);
      setError2D("2D rendering error");
    }
  }, [rdkit, currentMol, height2D, highlightAtoms, highlightColor]);

  // 3D rendering with 3Dmol
  useEffect(() => {
    if (!$3Dmol || !webglOK || !currentMol || !container3DRef.current) return;

    const container = container3DRef.current;
    setError3D(null);
    setLoading3D(true);

    try {
      // Nettoyage ancien viewer
      if (viewerRef.current) {
        try {
          viewerRef.current.removeAllModels();
          viewerRef.current.removeAllLabels();
          viewerRef.current.removeAllShapes();
        } catch (e) {
          // Ignore cleanup errors
        }
      }

      const viewer = $3Dmol.createViewer(container, {
        backgroundColor: "white",
        antialias: true,
        cartoonQuality: 6,
      });
      viewerRef.current = viewer;

      viewer.addModel(currentMol, "mol");
      viewer.setStyle({}, { stick: { radius: 0.15 }, sphere: { radius: 0.3 } });

      // Highlight atoms si besoin
      if (highlightAtoms.length > 0) {
        for (const idx of highlightAtoms) {
          viewer.setStyle(
            { serial: idx + 1 },
            {
              stick: { radius: 0.2, color: highlightColor },
              sphere: { radius: 0.4, color: highlightColor },
            }
          );
        }
      }

      viewer.zoomTo();
      viewer.render();
      viewer.spin("y", 0.5);

      debugLog("[OK] 3D mol rendered");
    } catch (e) {
      debugError("[ERROR] 3Dmol render error:", e);
      setError3D("3D rendering error");
    } finally {
      setLoading3D(false);
    }

    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.spin(false);
        } catch (e) {
          // Ignore
        }
      }
    };
  }, [$3Dmol, webglOK, currentMol, highlightAtoms, highlightColor]);

  // Cleanup viewer on unmount
  useEffect(() => {
    return () => {
      if (viewerRef.current) {
        try {
          viewerRef.current.spin(false);
          viewerRef.current.removeAllModels();
        } catch (e) {
          // Ignore cleanup errors
        }
      }
    };
  }, []);

  const toggleSpin = useCallback(() => {
    if (!viewerRef.current) return;
    try {
      const v = viewerRef.current;
      // Toggle spin
      v.spin(v.isAnimated() ? false : "y", 0.5);
    } catch (e) {
      debugWarn("Spin toggle failed:", e);
    }
  }, []);

  return (
    <div
      style={{
        display: "flex",
        flexDirection: "column",
        gap: 16,
        ...style,
      }}
    >
      {/* 2D Panel */}
      <div
        style={{
          background: "#fff",
          borderRadius: 10,
          boxShadow: "0 2px 8px rgba(0,0,0,0.08)",
          overflow: "hidden",
          position: "relative",
        }}
      >
        <div
          style={{
            padding: "8px 12px",
            borderBottom: "1px solid #e2e8f0",
            fontWeight: 600,
            fontSize: 13,
            color: "#334155",
          }}
        >
          Vue 2D
        </div>
        <div
          ref={container2DRef}
          style={{
            width: "100%",
            height: height2D,
            display: "flex",
            alignItems: "center",
            justifyContent: "center",
            background: "#fafafa",
          }}
        >
          {loading2D && (
            <span style={{ color: "#64748b", fontSize: 13 }}>
              Loading...
            </span>
          )}
          {error2D && (
            <span style={{ color: "#ef4444", fontSize: 13 }}>{error2D}</span>
          )}
          {!currentMol && !loading2D && !error2D && (
            <span style={{ color: "#94a3b8", fontSize: 13 }}>
              Aucune molécule
            </span>
          )}
        </div>
      </div>

      {/* 3D Panel */}
      {show3D && (
        <div
          style={{
            background: "#fff",
            borderRadius: 10,
            boxShadow: "0 2px 8px rgba(0,0,0,0.08)",
            overflow: "hidden",
            position: "relative",
          }}
        >
          <div
            style={{
              padding: "8px 12px",
              borderBottom: "1px solid #e2e8f0",
              fontWeight: 600,
              fontSize: 13,
              color: "#334155",
              display: "flex",
              justifyContent: "space-between",
              alignItems: "center",
            }}
          >
            <span>Vue 3D</span>
            <button
              onClick={toggleSpin}
              style={{
                padding: "2px 8px",
                borderRadius: 6,
                border: "1px solid #cbd5e1",
                background: "#fff",
                cursor: "pointer",
                fontSize: 11,
              }}
            >
              Rotation
            </button>
          </div>
          <div
            ref={container3DRef}
            style={{
              width: "100%",
              height: height3D,
              position: "relative",
              background: "#fafafa",
            }}
          >
            {loading3D && (
              <div
                style={{
                  position: "absolute",
                  inset: 0,
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "center",
                  color: "#64748b",
                  fontSize: 13,
                }}
              >
                Loading 3D...
              </div>
            )}
            {error3D && (
              <div
                style={{
                  position: "absolute",
                  inset: 0,
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "center",
                  color: "#ef4444",
                  fontSize: 13,
                }}
              >
                {error3D}
              </div>
            )}
            {!webglOK && !error3D && (
              <div
                style={{
                  position: "absolute",
                  inset: 0,
                  display: "flex",
                  alignItems: "center",
                  justifyContent: "center",
                  color: "#94a3b8",
                  fontSize: 13,
                }}
              >
                WebGL non supporté
              </div>
            )}
          </div>
        </div>
      )}
    </div>
  );
}

// Helper: hex color to RGB float array [r, g, b] for RDKit
function hexToRgbFloat(hex) {
  const result = /^#?([a-f\d]{2})([a-f\d]{2})([a-f\d]{2})$/i.exec(hex);
  if (!result) return [1, 0.4, 0.4];
  return [
    parseInt(result[1], 16) / 255,
    parseInt(result[2], 16) / 255,
    parseInt(result[3], 16) / 255,
  ];
}

export default memo(MoleculePane);
