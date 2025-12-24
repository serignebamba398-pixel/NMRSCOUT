const fs = require('fs');
const path = require('path');
const archiver = require('archiver');
const unzipper = require('unzipper');
const { v4: uuidv4 } = require('uuid');
const { projectsFile, outputsDir, uploadsDir } = require('../config/paths');
const { logWithMemory } = require('../utils/logger');

class ProjectService {
  constructor() {
    this.projects = [];
    this.currentProject = null;
    this.load();
  }

  load() {
    try {
      if (fs.existsSync(projectsFile)) {
        const data = fs.readFileSync(projectsFile, 'utf8');
        this.projects = JSON.parse(data);
      }
    } catch (e) {
      logWithMemory(`Project loading error: ${e.message}`, 'ERROR');
      this.projects = [];
    }
  }

  save() {
    try {
      fs.writeFileSync(projectsFile, JSON.stringify(this.projects, null, 2));
    } catch (e) {
      logWithMemory(`Project save error: ${e.message}`, 'ERROR');
    }
  }

  getAll() {
    return this.projects;
  }

  getById(id) {
    return this.projects.find(p => p.id === id);
  }

  getCurrent() {
    return this.currentProject;
  }

  setCurrent(id) {
    this.currentProject = this.getById(id);
    return this.currentProject;
  }

  create(name, description = '', workDir = null) {
    const id = uuidv4();
    const projectWorkDir = workDir || path.join(outputsDir, `project_${id}`);
    
    // Create the project folder
    if (!fs.existsSync(projectWorkDir)) {
      fs.mkdirSync(projectWorkDir, { recursive: true });
    }

    const newProject = {
      id,
      name,
      description,
      workDir: projectWorkDir,
      createdAt: new Date().toISOString(),
      updatedAt: new Date().toISOString(),
      tasks: [],
      state: {
        predictionState: null,
        annotationState: null,
        viewerData: null,
        customColors: {},
      }
    };
    this.projects.push(newProject);
    this.currentProject = newProject;
    this.save();
    return newProject;
  }

  update(id, updates) {
    const project = this.getById(id);
    if (project) {
      Object.assign(project, updates);
      project.updatedAt = new Date().toISOString();
      this.save();
      return project;
    }
    return null;
  }

  updateState(id, stateKey, stateValue) {
    const project = this.getById(id);
    if (project) {
      if (!project.state) project.state = {};
      project.state[stateKey] = stateValue;
      project.updatedAt = new Date().toISOString();
      this.save();
      return project;
    }
    return null;
  }

  delete(id) {
    const index = this.projects.findIndex(p => p.id === id);
    if (index !== -1) {
      const project = this.projects[index];
      // Delete the project folder if it exists
      if (project.workDir && fs.existsSync(project.workDir)) {
        try {
          fs.rmSync(project.workDir, { recursive: true, force: true });
        } catch (e) {
          logWithMemory(`Error deleting project folder: ${e.message}`, 'WARN');
        }
      }
      this.projects.splice(index, 1);
      if (this.currentProject?.id === id) {
        this.currentProject = null;
      }
      this.save();
      return true;
    }
    return false;
  }

  addTask(projectId, task) {
    const project = this.getById(projectId);
    if (project) {
      project.tasks.push({
        id: task.id,
        type: task.type,
        status: 'pending',
        createdAt: new Date().toISOString(),
        ...task
      });
      project.updatedAt = new Date().toISOString();
      this.save();
      return true;
    }
    return false;
  }

  updateTaskStatus(projectId, taskId, status, result = null) {
    const project = this.getById(projectId);
    if (project) {
      const task = project.tasks.find(t => t.id === taskId);
      if (task) {
        task.status = status;
        if (result) task.result = result;
        task.updatedAt = new Date().toISOString();
        this.save();
        return true;
      }
    }
    return false;
  }

  // Export a project to a .nmrproj file (ZIP)
  async exportProject(id, destPath) {
    const project = this.getById(id);
    if (!project) throw new Error('Project not found');

    return new Promise((resolve, reject) => {
      const output = fs.createWriteStream(destPath);
      const archive = archiver('zip', { zlib: { level: 9 } });

      output.on('close', () => {
        logWithMemory(`Project exported: ${destPath} (${archive.pointer()} bytes)`);
        resolve({ path: destPath, size: archive.pointer() });
      });

      archive.on('error', reject);
      archive.pipe(output);

      // Add project metadata
      archive.append(JSON.stringify(project, null, 2), { name: 'project.json' });

      // Add the project work directory
      if (project.workDir && fs.existsSync(project.workDir)) {
        archive.directory(project.workDir, 'data');
      }

      archive.finalize();
    });
  }

  // Import a project from a .nmrproj file
  async importProject(srcPath, newWorkDir = null) {
    return new Promise(async (resolve, reject) => {
      try {
        const tempDir = path.join(outputsDir, `temp_import_${Date.now()}`);
        fs.mkdirSync(tempDir, { recursive: true });

        // Extract the ZIP
        await fs.createReadStream(srcPath)
          .pipe(unzipper.Extract({ path: tempDir }))
          .promise();

        // Read metadata
        const projectJsonPath = path.join(tempDir, 'project.json');
        if (!fs.existsSync(projectJsonPath)) {
          throw new Error('Missing project.json in imported project');
        }

        const projectData = JSON.parse(fs.readFileSync(projectJsonPath, 'utf8'));
        
        // Generate a new ID
        const newId = uuidv4();
        const projectWorkDir = newWorkDir || path.join(outputsDir, `project_${newId}`);

        // Copy data
        const dataDir = path.join(tempDir, 'data');
        if (fs.existsSync(dataDir)) {
          fs.cpSync(dataDir, projectWorkDir, { recursive: true });
        } else {
          fs.mkdirSync(projectWorkDir, { recursive: true });
        }

        // Create the new project
        const newProject = {
          ...projectData,
          id: newId,
          workDir: projectWorkDir,
          importedAt: new Date().toISOString(),
          originalId: projectData.id,
        };

        this.projects.push(newProject);
        this.currentProject = newProject;
        this.save();

        // Clean up temp folder
        fs.rmSync(tempDir, { recursive: true, force: true });

        resolve(newProject);
      } catch (e) {
        reject(e);
      }
    });
  }

  // Get the outputs folder path for a project
  getProjectOutputsDir(projectId) {
    const project = this.getById(projectId);
    if (project && project.workDir) {
      return project.workDir;
    }
    return outputsDir;
  }
}

module.exports = new ProjectService();
