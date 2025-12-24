import React from 'react';
import { Link, useLocation } from 'react-router-dom';
import { useAppState } from '../../contexts/AppStateContext';

/**
 * Navigation header with conditional tab locking
 */
export const NavigationHeader = () => {
  const location = useLocation();
  const { uiLock } = useAppState();
  const currentPath = location.pathname;
  const lockedBy = uiLock.lockedBy;

  // Locking rules between sections
  const disablePred = lockedBy === "ANNO" || lockedBy === "OPTIM";
  const disableAnno = lockedBy === "PRED" || lockedBy === "OPTIM";
  const disableOptim = lockedBy === "PRED" || lockedBy === "ANNO";
  const disableRes = lockedBy === "PRED" || lockedBy === "ANNO" || lockedBy === "OPTIM";

  const TabLink = ({ to, disabled, children }) => (
    <Link
      to={disabled ? currentPath : to}
      onClick={disabled ? (e) => e.preventDefault() : undefined}
      className={`top-tab ${currentPath.includes(to) ? "active" : ""} ${
        disabled ? "disabled" : ""
      }`}
      aria-disabled={disabled ? "true" : "false"}
      tabIndex={disabled ? -1 : 0}
    >
      {children}
    </Link>
  );

  return (
    <header className="top-tabs">
      <TabLink to="/" disabled={false}>
        Home
      </TabLink>
      <TabLink to="/prediction" disabled={disablePred}>
        Predictions
      </TabLink>
      <TabLink to="/optimization" disabled={disableOptim}>
        Tolerance Optimization
      </TabLink>
      <TabLink to="/annotation" disabled={disableAnno}>
        Integrated Annotations
      </TabLink>
      <TabLink to="/results" disabled={disableRes}>
        Results
      </TabLink>
    </header>
  );
};