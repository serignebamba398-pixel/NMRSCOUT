import React from 'react';

/**
 * Summary table of tolerance optimization results
 */
export const ToleranceSummaryTable = ({ resultSet }) => {
  const { results, type, label } = resultSet;

  if (!results || results.length === 0) {
    return <p>No results to display</p>;
  }

  // Trouver le meilleur résultat selon les signaux matchés
  const bestResult = results.reduce((best, current) => {
    if (current.error) return best;
    if (!best) return current;
    return (current.totalSignalsMatched || 0) > (best.totalSignalsMatched || 0) ? current : best;
  }, null);

  const formatTolerance = (params) => {
    if (!params) return "N/A";
    
    if (type === "1H" || type === "13C") {
      return `tol: ${params.tol_correspondence?.toFixed(3) || params.tol_match?.toFixed(3) || "N/A"}`;
    } else if (type === "COSY" || type === "HSQC") {
      const f1 = params.tol_correspondence_f1?.toFixed(3) || params.tol_match_f1?.toFixed(3) || "N/A";
      const f2 = params.tol_correspondence_f2?.toFixed(3) || params.tol_match_f2?.toFixed(3) || "N/A";
      return `F1: ${f1}, F2: ${f2}`;
    }
    return "N/A";
  };

  const formatPercent = (value) => {
    return typeof value === 'number' ? `${value.toFixed(2)}%` : 'N/A';
  };

  const formatNumber = (value) => {
    return typeof value === 'number' ? Math.round(value) : 'N/A';
  };

  return (
    <div className="tolerance-summary" style={{ marginTop: 20 }}>
      <div className="summary-header" style={{ marginBottom: 15 }}>
        <h4>{label}</h4>
        {bestResult && !bestResult.error && (
          <div className="best-result-badge" style={{
            backgroundColor: '#e8f5e9',
            border: '1px solid #4caf50',
            borderRadius: 4,
            padding: '8px 12px',
            marginTop: 8,
            fontSize: '0.9em'
          }}>
            <strong>✓ Best result:</strong> {formatTolerance(bestResult.toleranceParams)} 
            → {bestResult.totalSignalsMatched} signals ({formatPercent(bestResult.signalsMatchedPct)})
          </div>
        )}
      </div>

      <div style={{ overflowX: 'auto' }}>
        <table style={{
          width: '100%',
          borderCollapse: 'collapse',
          fontSize: '0.9em',
          border: '1px solid #ddd'
        }}>
          <thead>
            <tr style={{ backgroundColor: '#f5f5f5', borderBottom: '2px solid #ddd' }}>
              <th style={{ padding: 10, textAlign: 'left', borderRight: '1px solid #ddd' }}>#</th>
              <th style={{ padding: 10, textAlign: 'left', borderRight: '1px solid #ddd' }}>Tolerances</th>
              <th style={{ padding: 10, textAlign: 'center', borderRight: '1px solid #ddd' }}>Signals Matched</th>
              <th style={{ padding: 10, textAlign: 'center', borderRight: '1px solid #ddd' }}>% Signals</th>
              <th style={{ padding: 10, textAlign: 'center', borderRight: '1px solid #ddd' }}>Molecules</th>
              <th style={{ padding: 10, textAlign: 'center', borderRight: '1px solid #ddd' }}>Avg Coverage</th>
              <th style={{ padding: 10, textAlign: 'center' }}>Status</th>
            </tr>
          </thead>
          <tbody>
            {results.map((result, index) => {
              const isError = !!result.error;
              const isBest = !isError && bestResult && result === bestResult;
              const rowBg = isBest ? '#e8f5e9' : (index % 2 === 0 ? '#fafafa' : '#fff');
              const borderLeft = isBest ? '4px solid #4caf50' : 'none';

              return (
                <tr
                  key={index}
                  style={{
                    backgroundColor: rowBg,
                    borderBottom: '1px solid #eee',
                    borderLeft
                  }}
                >
                  <td style={{ padding: 10, borderRight: '1px solid #ddd', fontWeight: 'bold' }}>
                    {index + 1}
                  </td>
                  <td style={{ padding: 10, borderRight: '1px solid #ddd', fontFamily: 'monospace', fontSize: '0.85em' }}>
                    {isError ? (
                      <span style={{ color: '#d32f2f' }}>Error</span>
                    ) : (
                      formatTolerance(result.toleranceParams)
                    )}
                  </td>
                  <td style={{ padding: 10, borderRight: '1px solid #ddd', textAlign: 'center', fontWeight: 'bold' }}>
                    {isError ? '—' : formatNumber(result.totalSignalsMatched)}
                  </td>
                  <td style={{ padding: 10, borderRight: '1px solid #ddd', textAlign: 'center' }}>
                    {isError ? '—' : formatPercent(result.signalsMatchedPct)}
                  </td>
                  <td style={{ padding: 10, borderRight: '1px solid #ddd', textAlign: 'center' }}>
                    {isError ? '—' : formatNumber(result.totalMolecules)}
                  </td>
                  <td style={{ padding: 10, borderRight: '1px solid #ddd', textAlign: 'center' }}>
                    {isError ? '—' : formatPercent(result.avgCoverage)}
                  </td>
                  <td style={{ padding: 10, textAlign: 'center' }}>
                    {isError ? (
                      <span style={{ color: '#d32f2f', fontWeight: 'bold' }}>✗ Error</span>
                    ) : (
                      <span style={{ color: '#4caf50', fontWeight: 'bold' }}>✓ OK</span>
                    )}
                  </td>
                </tr>
              );
            })}
          </tbody>
        </table>
      </div>

      {/* Global statistics */}
      <div className="summary-stats" style={{
        marginTop: 20,
        padding: 15,
        backgroundColor: '#f9f9f9',
        borderRadius: 4,
        border: '1px solid #eee'
      }}>
        <h5>Global Statistics</h5>
        <div style={{ display: 'grid', gridTemplateColumns: 'repeat(auto-fit, minmax(200px, 1fr))', gap: 15 }}>
          <div>
            <strong>Total tests:</strong> {results.length}
          </div>
          <div>
            <strong>Successful tests:</strong> {results.filter(r => !r.error).length}
          </div>
          <div>
            <strong>Failed tests:</strong> {results.filter(r => r.error).length}
          </div>
          {bestResult && !bestResult.error && (
            <>
              <div>
                <strong>Max signals matched:</strong> {bestResult.totalSignalsMatched}
              </div>
              <div>
                <strong>Max % signals:</strong> {formatPercent(bestResult.signalsMatchedPct)}
              </div>
              <div>
                <strong>Best coverage:</strong> {formatPercent(bestResult.avgCoverage)}
              </div>
            </>
          )}
        </div>
      </div>
    </div>
  );
};
