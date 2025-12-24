/**
 * Class to manage multi-phase progression with interpolation
 */
export class ProgressTracker {
  constructor() {
    this.phases = [];
    this.currentPhase = 0;
    this.currentPhaseProgress = 0;
    this.lastReportedProgress = 0;
    this.interpolationTimer = null;
  }

  /**
   * Defines phases with their relative weights
   * @param {Array<{name: string, weight: number}>} phases
   */
  definePhases(phases) {
    const totalWeight = phases.reduce((sum, p) => sum + p.weight, 0);
    
    this.phases = phases.map((phase, index) => ({
      ...phase,
      startPercent: phases
        .slice(0, index)
        .reduce((sum, p) => sum + (p.weight / totalWeight * 100), 0),
      endPercent: phases
        .slice(0, index + 1)
        .reduce((sum, p) => sum + (p.weight / totalWeight * 100), 0),
    }));
  }

  /**
   * Updates the progression of a phase
   * @param {number} phaseIndex - Phase index (0-based)
   * @param {number} phaseProgress - Progression from 0 to 100 in this phase
   * @returns {number} Global progression (0-100)
   */
  updatePhase(phaseIndex, phaseProgress) {
    if (phaseIndex >= this.phases.length) {
      return this.lastReportedProgress;
    }
    
    this.currentPhase = phaseIndex;
    this.currentPhaseProgress = Math.max(0, Math.min(100, phaseProgress));
    
    const phase = this.phases[phaseIndex];
    const phaseRange = phase.endPercent - phase.startPercent;
    const globalProgress = phase.startPercent + 
      (this.currentPhaseProgress / 100) * phaseRange;
    
    this.lastReportedProgress = Math.floor(globalProgress);
    return this.lastReportedProgress;
  }

  /**
   * Starts a progressive interpolation towards a target percentage
   * @param {number} targetPercent - Target (0-100)
   * @param {number} durationMs - Animation duration
   * @param {function} onUpdate - Callback called on each frame
   */
  startInterpolation(targetPercent, durationMs = 2000, onUpdate) {
    if (this.interpolationTimer) {
      clearInterval(this.interpolationTimer);
    }
    
    const startPercent = this.lastReportedProgress;
    const startTime = Date.now();
    
    this.interpolationTimer = setInterval(() => {
      const elapsed = Date.now() - startTime;
      const progress = Math.min(1, elapsed / durationMs);
      
      // Easing function (ease-out quadratic)
      const eased = 1 - Math.pow(1 - progress, 2);
      const interpolated = startPercent + (targetPercent - startPercent) * eased;
      
      this.lastReportedProgress = Math.floor(interpolated);
      onUpdate(this.lastReportedProgress);
      
      if (progress >= 1) {
        clearInterval(this.interpolationTimer);
        this.interpolationTimer = null;
      }
    }, 50);
  }

  /**
   * Stops the current interpolation
   */
  stopInterpolation() {
    if (this.interpolationTimer) {
      clearInterval(this.interpolationTimer);
      this.interpolationTimer = null;
    }
  }

  /**
   * Resets the tracker
   */
  reset() {
    this.currentPhase = 0;
    this.currentPhaseProgress = 0;
    this.lastReportedProgress = 0;
    this.stopInterpolation();
  }

  /**
   * Gets the current progress
   * @returns {number} Progression globale (0-100)
   */
  getProgress() {
    return this.lastReportedProgress;
  }
}