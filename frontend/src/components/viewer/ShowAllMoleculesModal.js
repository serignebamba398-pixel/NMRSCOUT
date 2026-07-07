// src/components/viewer/ShowAllMoleculesModal.js
// Modal to display all molecules of a project in a paginated grid

import { useState, useEffect, useMemo, useCallback, memo } from "react";
import { loadRDKitFromCDN } from "../../utils/cdnLoaders";
import { extractFirstMolBlockPreserving, tryFetchText, stripExplicitHFromMolblock } from "../../utils/molParser";
import { debugLog, debugError } from "../../utils/debug";
import { API_BASE_URL } from "../../api";

const ITEMS_PER_PAGE = 12;

/**
 * Fullscreen modal to display all molecules of a project
 */
function ShowAllMoleculesModal({
  molecules = [],
  currentIndex = 0,
  onSelectMolecule,
  onClose,
}) {
  const [rdkit, setRdkit] = useState(null);
  const [page, setPage] = useState(Math.floor(currentIndex / ITEMS_PER_PAGE));
  const [molBlocks, setMolBlocks] = useState({});
  const [svgCache, setSvgCache] = useState({});
  const [loading, setLoading] = useState(false);

  // Loading RDKit
  useEffect(() => {
    loadRDKitFromCDN()
      .then((kit) => {
        setRdkit(kit);
        debugLog("[OK] RDKit loaded in ShowAllMoleculesModal");
      })
      .catch((e) => {
        debugError("[ERROR] RDKit load failed:", e);
      });
  }, []);

  // Pagination
  const totalPages = Math.ceil(molecules.length / ITEMS_PER_PAGE);
  const pageItems = useMemo(() => {
    const start = page * ITEMS_PER_PAGE;
    return molecules.slice(start, start + ITEMS_PER_PAGE);
  }, [molecules, page]);

  // Loading molBlocks for current page items
  useEffect(() => {
    if (!pageItems.length) return;

    const toFetch = pageItems.filter(
      (m) => !molBlocks[m.filename] && !m.molBlock
    );
    if (toFetch.length === 0) return;

    setLoading(true);
    Promise.all(
      toFetch.map(async (m) => {
        try {
          let mol = null;

          // Try direct URL first
          if (m.sdf_url) {
            const url = m.sdf_url.startsWith("http")
              ? m.sdf_url
              : `${API_BASE_URL}${m.sdf_url}`;
            mol = await tryFetchText(url);
          }

          // Try via filename
          if (!mol && m.filename) {
            const url = `${API_BASE_URL}/outputs/${encodeURIComponent(
              m.filename
            )}`;
            mol = await tryFetchText(url);
          }

          if (mol) {
            return { filename: m.filename, molBlock: extractFirstMolBlockPreserving(mol) };
          }
          return { filename: m.filename, molBlock: null };
        } catch (e) {
          debugError(`Fetch mol failed for ${m.filename}:`, e);
          return { filename: m.filename, molBlock: null };
        }
      })
    ).then((results) => {
      const newBlocks = {};
      for (const r of results) {
        if (r.molBlock) {
          newBlocks[r.filename] = r.molBlock;
        }
      }
      setMolBlocks((prev) => ({ ...prev, ...newBlocks }));
      setLoading(false);
    });
  }, [pageItems, molBlocks]);

  // SVG generation for each molecule
  useEffect(() => {
    if (!rdkit) return;

    for (const m of pageItems) {
      const mol = m.molBlock || molBlocks[m.filename];
      if (!mol || svgCache[m.filename]) continue;

      try {
        const { molblock: cleanMol } = stripExplicitHFromMolblock(mol);
        const tmpMol = rdkit.get_mol(cleanMol);
        if (!tmpMol) continue;
        // Convertir en SMILES puis recréer pour obtenir une 2D propre sans coordonnées 3D
        const smiles = tmpMol.get_smiles();
        tmpMol.delete();
        if (!smiles) continue;
        const rdkMol = rdkit.get_mol(smiles);
        if (!rdkMol) continue;

        const opts = {
          width: 380,
          height: 300,
          backgroundColour: [1, 1, 1, 1],
          addStereoAnnotation: false,
          explicitMethyl: false,
          bondLineWidth: 2.0,
          atomLabelFontSize: 14,
          useMolBlockCoords: false,
          padding: 0.0,
        };
        const svg = rdkMol.get_svg(JSON.stringify(opts));
        rdkMol.delete();

        setSvgCache((prev) => ({ ...prev, [m.filename]: svg }));
      } catch (e) {
        debugError(`SVG gen failed for ${m.filename}:`, e);
      }
    }
  }, [rdkit, pageItems, molBlocks, svgCache]);

  const handleSelect = useCallback(
    (idx) => {
      onSelectMolecule && onSelectMolecule(idx);
      onClose && onClose();
    },
    [onSelectMolecule, onClose]
  );

  const handleKeyDown = useCallback(
    (e) => {
      if (e.key === "Escape") {
        onClose && onClose();
      } else if (e.key === "ArrowLeft" && page > 0) {
        setPage((p) => p - 1);
      } else if (e.key === "ArrowRight" && page < totalPages - 1) {
        setPage((p) => p + 1);
      }
    },
    [onClose, page, totalPages]
  );

  useEffect(() => {
    window.addEventListener("keydown", handleKeyDown);
    return () => window.removeEventListener("keydown", handleKeyDown);
  }, [handleKeyDown]);

  return (
    <div
      style={{
        position: "fixed",
        inset: 0,
        zIndex: 9999,
        background: "rgba(0,0,0,0.6)",
        display: "flex",
        alignItems: "center",
        justifyContent: "center",
        backdropFilter: "blur(4px)",
      }}
      onClick={(e) => {
        if (e.target === e.currentTarget) onClose && onClose();
      }}
    >
      <div
        style={{
          background: "#fff",
          borderRadius: 16,
          maxWidth: "90vw",
          maxHeight: "90vh",
          width: 900,
          display: "flex",
          flexDirection: "column",
          boxShadow: "0 25px 50px -12px rgba(0,0,0,0.25)",
          overflow: "hidden",
        }}
      >
        {/* Header */}
        <div
          style={{
            padding: "16px 24px",
            borderBottom: "1px solid #e2e8f0",
            display: "flex",
            justifyContent: "space-between",
            alignItems: "center",
          }}
        >
          <h2 style={{ margin: 0, fontSize: 18, fontWeight: 600 }}>
            All molecules ({molecules.length})
          </h2>
          <button
            onClick={onClose}
            style={{
              background: "none",
              border: "none",
              fontSize: 24,
              cursor: "pointer",
              color: "#64748b",
              padding: 4,
            }}
          >
            ×
          </button>
        </div>

        {/* Grid */}
        <div
          style={{
            flex: 1,
            overflowY: "auto",
            padding: 20,
          }}
        >
          {loading && pageItems.length > 0 && (
            <div
              style={{
                textAlign: "center",
                color: "#64748b",
                padding: 16,
                fontSize: 14,
              }}
            >
              Loading structures...
            </div>
          )}
          <div
            style={{
              display: "grid",
              gridTemplateColumns: "repeat(auto-fill, minmax(320px, 1fr))",
              gap: 16,
            }}
          >
            {pageItems.map((m, i) => {
              const globalIdx = page * ITEMS_PER_PAGE + i;
              const isSelected = globalIdx === currentIndex;
              const svg = svgCache[m.filename];

              return (
                <div
                  key={m.filename || globalIdx}
                  onClick={() => handleSelect(globalIdx)}
                  style={{
                    background: isSelected ? "#eff6ff" : "#fafafa",
                    border: isSelected
                      ? "2px solid #3b82f6"
                      : "1px solid #e2e8f0",
                    borderRadius: 12,
                    padding: 12,
                    cursor: "pointer",
                    transition: "all 0.15s ease",
                    boxShadow: isSelected
                      ? "0 4px 12px rgba(59,130,246,0.2)"
                      : "0 1px 3px rgba(0,0,0,0.05)",
                  }}
                  onMouseEnter={(e) => {
                    if (!isSelected) {
                      e.currentTarget.style.borderColor = "#93c5fd";
                      e.currentTarget.style.background = "#f8fafc";
                    }
                  }}
                  onMouseLeave={(e) => {
                    if (!isSelected) {
                      e.currentTarget.style.borderColor = "#e2e8f0";
                      e.currentTarget.style.background = "#fafafa";
                    }
                  }}
                >
                  {/* Molecule SVG */}
                  <div
                    style={{
                      width: "100%",
                      height: 260,
                      display: "flex",
                      alignItems: "center",
                      justifyContent: "center",
                      background: "#fff",
                      borderRadius: 8,
                      marginBottom: 8,
                    }}
                  >
                    {svg ? (
                      <div dangerouslySetInnerHTML={{ __html: svg }} />
                    ) : (
                      <span style={{ color: "#94a3b8", fontSize: 12 }}>
                        {loading ? "..." : "No structure"}
                      </span>
                    )}
                  </div>

                  {/* Label */}
                  <div
                    style={{
                      fontSize: 12,
                      fontWeight: 500,
                      color: "#334155",
                      textAlign: "center",
                      overflow: "hidden",
                      textOverflow: "ellipsis",
                      whiteSpace: "nowrap",
                    }}
                    title={m.filename || m.name || `Mol ${globalIdx + 1}`}
                  >
                    {m.name ||
                      (m.filename || "").replace(/\.(sdf|mol)$/i, "") ||
                      `Mol ${globalIdx + 1}`}
                  </div>

                  {/* Index badge */}
                  <div
                    style={{
                      textAlign: "center",
                      marginTop: 4,
                    }}
                  >
                    <span
                      style={{
                        fontSize: 10,
                        background: isSelected ? "#3b82f6" : "#e2e8f0",
                        color: isSelected ? "#fff" : "#64748b",
                        padding: "2px 8px",
                        borderRadius: 10,
                      }}
                    >
                      #{globalIdx + 1}
                    </span>
                  </div>
                </div>
              );
            })}
          </div>

          {molecules.length === 0 && (
            <div
              style={{
                textAlign: "center",
                color: "#94a3b8",
                padding: 40,
              }}
            >
              No molecule in this project
            </div>
          )}
        </div>

        {/* Pagination */}
        {totalPages > 1 && (
          <div
            style={{
              padding: "12px 24px",
              borderTop: "1px solid #e2e8f0",
              display: "flex",
              justifyContent: "center",
              alignItems: "center",
              gap: 12,
            }}
          >
            <button
              onClick={() => setPage((p) => Math.max(0, p - 1))}
              disabled={page === 0}
              style={{
                padding: "6px 16px",
                borderRadius: 8,
                border: "1px solid #cbd5e1",
                background: page === 0 ? "#f1f5f9" : "#fff",
                color: page === 0 ? "#94a3b8" : "#334155",
                cursor: page === 0 ? "not-allowed" : "pointer",
                fontWeight: 500,
              }}
            >
              ← Previous
            </button>
            <span style={{ fontSize: 14, color: "#64748b" }}>
              Page {page + 1} / {totalPages}
            </span>
            <button
              onClick={() => setPage((p) => Math.min(totalPages - 1, p + 1))}
              disabled={page >= totalPages - 1}
              style={{
                padding: "6px 16px",
                borderRadius: 8,
                border: "1px solid #cbd5e1",
                background: page >= totalPages - 1 ? "#f1f5f9" : "#fff",
                color: page >= totalPages - 1 ? "#94a3b8" : "#334155",
                cursor: page >= totalPages - 1 ? "not-allowed" : "pointer",
                fontWeight: 500,
              }}
            >
              Next →
            </button>
          </div>
        )}
      </div>
    </div>
  );
}

export default memo(ShowAllMoleculesModal);
