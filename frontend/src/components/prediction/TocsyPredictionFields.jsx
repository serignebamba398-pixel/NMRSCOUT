import React from 'react';

export const TocsyPredictionFields = ({ options = {}, onChange, disabled }) => (
  <>
    <div className="params-grid">
      <label className="param-row">
        <span>Minimum bonds:</span>
        <input type="number" min="1" step="1" aria-label="Minimum bonds"
          title="Minimum H-H topological distance retained within a detected spin system (diagonal excluded)."
          value={options.min_bonds ?? 2}
          onChange={e => onChange('min_bonds', e.target.value)} disabled={disabled} />
      </label>
      <label className="param-row">
        <span>Maximum bonds:</span>
        <input type="number" min="1" step="1" aria-label="Maximum bonds"
          title="Maximum H-H topological distance retained within a detected spin system."
          value={options.max_bonds ?? 3}
          onChange={e => onChange('max_bonds', e.target.value)} disabled={disabled} />
      </label>
    </div>
    <label className="checkbox-label">
      <input type="checkbox" checked={options.include_diagonal ?? true}
        onChange={e => onChange('include_diagonal', e.target.checked)} disabled={disabled} />
      <span>Include diagonal peaks</span>
    </label>
    <p className="hint">
      Spin systems are detected automatically from the structure. Minimum and maximum bonds filter
      the retained H-H distance within each system, without changing its boundaries.
    </p>
  </>
);
