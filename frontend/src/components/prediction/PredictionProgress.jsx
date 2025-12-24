import React from 'react';
import { ProgressBar } from '../common';

/**
 * Prediction progress display
 */
export const PredictionProgress = ({
  progress,
  progressMessage,
  completedMolecules,
  totalMolecules,
}) => {
  return (
    <div className="progress-section">
      <h4>Progress</h4>
      <ProgressBar progress={progress} message={progressMessage} />
      <div
        style={{
          marginTop: 6,
          fontSize: 12,
          color: "#64748b",
          display: "flex",
          justifyContent: "space-between",
          alignItems: "center",
        }}
      >
        <span>
          {completedMolecules} / {totalMolecules} molecules processed
        </span>
        {totalMolecules > 0 && (
          <span>
            SDF phase: {((completedMolecules / totalMolecules) * 100).toFixed(1)}%
          </span>
        )}
      </div>
    </div>
  );
};