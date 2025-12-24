import React from 'react';
import { API_BASE_URL } from '../../api';

/**
 * Footer with backend information
 */
export const Footer = () => {
  return (
    <footer className="footer">
      Backend: <code>{API_BASE_URL}</code>
    </footer>
  );
};