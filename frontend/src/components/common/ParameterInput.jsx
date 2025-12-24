import React from 'react';

/**
 * Parameter input field with label
 * @param {string} label - Parameter label
 * @param {any} value - Current value
 * @param {function} onChange - Change callback
 * @param {string} type - Input type (number, text, etc.)
 * @param {string} step - Step for numbers (e.g. "0.01")
 * @param {boolean} disabled - Disable the field
 * @param {string} placeholder - Optional placeholder
 */
export const ParameterInput = ({
  label,
  value,
  onChange,
  type = "number",
  step = "any",
  disabled = false,
  placeholder = "",
}) => (
  <div className="parameter-input">
    <label className="parameter-label">{label}</label>
    <input
      type={type}
      className="parameter-field"
      step={step}
      value={value}
      onChange={onChange}
      disabled={disabled}
      placeholder={placeholder}
      aria-label={label}
    />
  </div>
);