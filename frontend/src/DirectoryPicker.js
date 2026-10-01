// frontend/src/DirectoryPicker.js
import React, { useEffect, useState } from 'react';
import axios from 'axios';

export default function DirectoryPicker({ value, onSelect }) {
  const [folders, setFolders] = useState([]);

  useEffect(() => {
    // Load output directories when the picker mounts.
    axios.get(`${process.env.REACT_APP_API_URL}/api/list`)
      .then(res => setFolders(res.data.folders))
      .catch(() => setFolders([]));
  }, []);

  return (
    <div className="mb-2">
      <select
        className="border p-2 w-full"
        value={value}
        onChange={e => onSelect(e.target.value)}
      >
        <option value="">-- Choisir un dossier --</option>
        {folders.map(f => (
          <option key={f.path} value={f.path}>
            {f.name}
          </option>
        ))}
      </select>
    </div>
  );
}
