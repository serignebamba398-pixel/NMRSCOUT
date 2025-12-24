import React, { useRef, useEffect } from 'react';
import { formatFileSize } from '../../utils/formatters';

/**
 * File upload component with status display
 * @param {string} label - Field label
 * @param {File} file - Currently selected file
 * @param {function} onChange - Callback called on selection
 * @param {string} accept - Accepted file types (e.g. ".csv")
 * @param {boolean} disabled - Disable the field
 * @param {boolean} uploading - Show "uploading" state
 */
export const FileUpload = ({
  label,
  file,
  onChange,
  accept = ".csv",
  disabled = false,
  uploading = false,
}) => {
  const inputRef = useRef(null);

  // Reset input if file is removed
  useEffect(() => {
    if (!file && inputRef.current) {
      inputRef.current.value = "";
    }
  }, [file]);

  return (
    <div className="form-group">
      <label className="form-label">{label}</label>
      <input
        ref={inputRef}
        type="file"
        accept={accept}
        onChange={onChange}
        className="form-input"
        disabled={disabled}
        aria-label={label}
      />
      
      {uploading && (
        <div className="file-status uploading">
          <span className="spinner" aria-hidden="true" /> Uploading...
        </div>
      )}
      
      {file && !uploading && (
        <div className="file-status success">
          Uploaded: <strong>{file.name}</strong> ({formatFileSize(file.size)})
        </div>
      )}
    </div>
  );
};