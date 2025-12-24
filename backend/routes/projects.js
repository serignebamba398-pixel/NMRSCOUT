const express = require('express');
const router = express.Router();
const path = require('path');
const fs = require('fs');
const projectService = require('../services/projectService');
const { outputsDir } = require('../config/paths');

// GET /api/projects
router.get('/api/projects', (req, res) => {
  res.json(projectService.getAll());
});

// GET /api/projects/current
router.get('/api/projects/current', (req, res) => {
  const current = projectService.getCurrent();
  if (current) {
    res.json(current);
  } else {
    res.json(null);
  }
});

// GET /api/projects/:id
router.get('/api/projects/:id', (req, res) => {
  const project = projectService.getById(req.params.id);
  if (project) {
    res.json(project);
  } else {
    res.status(404).json({ error: 'Project not found' });
  }
});

// POST /api/projects
router.post('/api/projects', (req, res) => {
  const { name, description, workDir } = req.body;
  if (!name) {
    return res.status(400).json({ error: 'Project name is required' });
  }
  const project = projectService.create(name, description, workDir);
  res.json(project);
});

// PUT /api/projects/:id
router.put('/api/projects/:id', (req, res) => {
  const project = projectService.update(req.params.id, req.body);
  if (project) {
    res.json(project);
  } else {
    res.status(404).json({ error: 'Project not found' });
  }
});

// PUT /api/projects/:id/state
router.put('/api/projects/:id/state', (req, res) => {
  const { key, value } = req.body;
  if (!key) {
    return res.status(400).json({ error: 'State key is required' });
  }
  const project = projectService.updateState(req.params.id, key, value);
  if (project) {
    res.json(project);
  } else {
    res.status(404).json({ error: 'Project not found' });
  }
});

// POST /api/projects/:id/set-current
router.post('/api/projects/:id/set-current', (req, res) => {
  const project = projectService.setCurrent(req.params.id);
  if (project) {
    res.json(project);
  } else {
    res.status(404).json({ error: 'Project not found' });
  }
});

// DELETE /api/projects/:id
router.delete('/api/projects/:id', (req, res) => {
  const success = projectService.delete(req.params.id);
  if (success) {
    res.json({ success: true });
  } else {
    res.status(404).json({ error: 'Project not found' });
  }
});

// POST /api/projects/:id/export
router.post('/api/projects/:id/export', async (req, res) => {
  try {
    const { destPath } = req.body;
    if (!destPath) {
      return res.status(400).json({ error: 'Destination path is required' });
    }
    const result = await projectService.exportProject(req.params.id, destPath);
    res.json(result);
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

// POST /api/projects/import
router.post('/api/projects/import', async (req, res) => {
  try {
    const { srcPath, workDir } = req.body;
    if (!srcPath) {
      return res.status(400).json({ error: 'Source path is required' });
    }
    const project = await projectService.importProject(srcPath, workDir);
    res.json(project);
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

// GET /api/projects/:id/outputs-dir
router.get('/api/projects/:id/outputs-dir', (req, res) => {
  const dir = projectService.getProjectOutputsDir(req.params.id);
  res.json({ outputsDir: dir });
});

// GET /api/outputs-dir - returns current outputs folder
router.get('/api/outputs-dir', (req, res) => {
  res.json({ outputsDir });
});

// GET /api/projects/:id/files - list a project's files
router.get('/api/projects/:id/files', (req, res) => {
  const project = projectService.getById(req.params.id);
  if (!project) {
    return res.status(404).json({ error: 'Project not found' });
  }
  
  const workDir = project.workDir;
  if (!workDir || !fs.existsSync(workDir)) {
    return res.json({ files: [] });
  }

  const listFilesRecursive = (dir, baseDir = dir) => {
    const files = [];
    const items = fs.readdirSync(dir);
    for (const item of items) {
      const fullPath = path.join(dir, item);
      const relativePath = path.relative(baseDir, fullPath);
      const stat = fs.statSync(fullPath);
      if (stat.isDirectory()) {
        files.push(...listFilesRecursive(fullPath, baseDir));
      } else {
        files.push({
          name: item,
          path: relativePath,
          size: stat.size,
          modified: stat.mtime,
        });
      }
    }
    return files;
  };

  try {
    const files = listFilesRecursive(workDir);
    res.json({ files, workDir });
  } catch (e) {
    res.status(500).json({ error: e.message });
  }
});

module.exports = router;
