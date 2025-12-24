// src/components/viewer/Pill.js
// Composant pastille métrique

import React from "react";

/**
 * Petite pastille affichant une métrique (label: value)
 */
const Pill = ({ label, value }) => (
  <div
    style={{
      background: "#f2f7ff",
      border: "1px solid #dbeafe",
      borderRadius: 999,
      padding: "2px 6px",
      fontSize: 10,
      display: "inline-flex",
      gap: 3,
      alignItems: "center",
      marginRight: 4,
      marginBottom: 3,
    }}
  >
    <span style={{ color: "#2563eb", fontWeight: 600 }}>{label}:</span>
    <span style={{ color: "#0f172a", fontWeight: 700 }}>{String(value)}</span>
  </div>
);

export default Pill;
