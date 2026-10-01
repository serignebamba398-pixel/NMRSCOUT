// optimization.js - Complete optimization system for ML NMR predictions
const fs = require('fs');
const path = require('path');
const crypto = require('crypto');
const { Worker } = require('worker_threads');

/**
 * Intelligent cache for ML predictions
 */
class ModelCache {
  constructor(options = {}) {
    this.cache = new Map();
    this.maxCacheSize = options.maxCacheSize || 100;
    this.maxAgeDays = options.maxAgeDays || 7;
    this.cacheDir = options.cacheDir || path.join(__dirname, 'cache');
    this.hitCount = 0;
    this.missCount = 0;
    
    this.ensureCacheDir();
    this.startCleanupTimer();
  }

  ensureCacheDir() {
    if (!fs.existsSync(this.cacheDir)) {
      fs.mkdirSync(this.cacheDir, { recursive: true });
      console.log(`[DIR] Cache folder created: ${this.cacheDir}`);
    }
  }

  startCleanupTimer() {
    // Automatic cleanup every 6 hours
    setInterval(() => {
      this.cleanup();
    }, 6 * 60 * 60 * 1000);
  }

  generateCacheKey(smiles, predictionType, parameters = {}) {
    const normalizedParams = this.normalizeParameters(parameters);
    const data = JSON.stringify({ 
      smiles: smiles.trim().toLowerCase(), 
      predictionType, 
      parameters: normalizedParams 
    });
    return crypto.createHash('sha256').update(data).digest('hex').substring(0, 16);
  }

  normalizeParameters(params) {
    // Normalize parameters to avoid minor variations
    const normalized = {};
    Object.keys(params).sort().forEach(key => {
      const value = params[key];
      if (typeof value === 'number') {
        normalized[key] = Math.round(value * 1000) / 1000; // 3 decimal places
      } else {
        normalized[key] = value;
      }
    });
    return normalized;
  }

  has(cacheKey) {
    // Check in memory
    if (this.cache.has(cacheKey)) {
      return true;
    }
    
    // Check on disk
    const cacheFile = path.join(this.cacheDir, `${cacheKey}.json`);
    return fs.existsSync(cacheFile);
  }

  get(cacheKey) {
    try {
      // Try memory first
      if (this.cache.has(cacheKey)) {
        const cached = this.cache.get(cacheKey);
        cached.lastAccessed = Date.now();
        cached.accessCount = (cached.accessCount || 0) + 1;
        this.hitCount++;
        console.log(`[CACHE] Cache hit (memory): ${cacheKey}`);
        return cached.data;
      }

      // Try disk
      const cacheFile = path.join(this.cacheDir, `${cacheKey}.json`);
      if (fs.existsSync(cacheFile)) {
        const cached = JSON.parse(fs.readFileSync(cacheFile, 'utf8'));
        
        // Check cache age
        const age = Date.now() - cached.cachedAt;
        const maxAge = this.maxAgeDays * 24 * 60 * 60 * 1000;
        
        if (age > maxAge) {
          fs.unlinkSync(cacheFile);
          this.missCount++;
          return null;
        }
        
        cached.lastAccessed = Date.now();
        cached.accessCount = (cached.accessCount || 0) + 1;
        
        // Put back in memory if not too full
        if (this.cache.size < this.maxCacheSize) {
          this.cache.set(cacheKey, {
            data: cached.data,
            cachedAt: cached.cachedAt,
            lastAccessed: cached.lastAccessed,
            accessCount: cached.accessCount
          });
        }
        
        this.hitCount++;
        console.log(`[CACHE] Cache hit (disk): ${cacheKey}`);
        return cached.data;
      }

      this.missCount++;
      return null;
      
    } catch (error) {
      console.error(`[ERROR] Cache read error ${cacheKey}:`, error);
      this.missCount++;
      return null;
    }
  }

  set(cacheKey, data) {
    try {
      const cacheEntry = {
        data: data,
        cachedAt: Date.now(),
        lastAccessed: Date.now(),
        accessCount: 1
      };

      // Save in memory
      if (this.cache.size >= this.maxCacheSize) {
        this.evictOldest();
      }
      this.cache.set(cacheKey, cacheEntry);

      // Save to disk
      const cacheFile = path.join(this.cacheDir, `${cacheKey}.json`);
      fs.writeFileSync(cacheFile, JSON.stringify(cacheEntry, null, 2));
      
      console.log(`[CACHE] Cache saved: ${cacheKey}`);
      
    } catch (error) {
      console.error(`[ERROR] Cache save error ${cacheKey}:`, error);
    }
  }

  evictOldest() {
    let oldestKey = null;
    let oldestTime = Date.now();

    for (const [key, entry] of this.cache.entries()) {
      if (entry.lastAccessed < oldestTime) {
        oldestTime = entry.lastAccessed;
        oldestKey = key;
      }
    }

    if (oldestKey) {
      this.cache.delete(oldestKey);
      console.log(`Cache eviction: ${oldestKey}`);
    }
  }

  cleanup(force = false) {
    const maxAge = this.maxAgeDays * 24 * 60 * 60 * 1000;
    const now = Date.now();
    let deletedMemory = 0;
    let deletedDisk = 0;

    // Clean memory
    for (const [key, entry] of this.cache.entries()) {
      if (force || (now - entry.cachedAt > maxAge)) {
        this.cache.delete(key);
        deletedMemory++;
      }
    }

    // Clean disk
    try {
      const files = fs.readdirSync(this.cacheDir);

      for (const file of files) {
        if (file.endsWith('.json')) {
          const filePath = path.join(this.cacheDir, file);
          const stats = fs.statSync(filePath);
          
          if (force || (now - stats.mtime.getTime() > maxAge)) {
            fs.unlinkSync(filePath);
            deletedDisk++;
          }
        }
      }

      if (deletedMemory > 0 || deletedDisk > 0) {
        console.log(`[CLEAN] Cache cleaned: ${deletedMemory} memory, ${deletedDisk} disk`);
      }
      
    } catch (error) {
      console.error('[ERROR] Cache cleanup error:', error);
    }
  }

  getStats() {
    const memoryEntries = this.cache.size;
    let diskEntries = 0;
    let totalSizeBytes = 0;

    try {
      const files = fs.readdirSync(this.cacheDir);
      for (const file of files) {
        if (file.endsWith('.json')) {
          diskEntries++;
          const filePath = path.join(this.cacheDir, file);
          totalSizeBytes += fs.statSync(filePath).size;
        }
      }
    } catch (error) {
      console.error('[ERROR] Cache stats calculation error:', error);
    }

    const total = this.hitCount + this.missCount;
    const hitRate = total > 0 ? (this.hitCount / total * 100).toFixed(1) : 0;

    return {
      memoryEntries,
      diskEntries,
      totalSizeMB: (totalSizeBytes / 1024 / 1024).toFixed(2),
      hitCount: this.hitCount,
      missCount: this.missCount,
      hitRate: hitRate,
      maxCacheSize: this.maxCacheSize,
      maxAgeDays: this.maxAgeDays
    };
  }

  clear() {
    this.cache.clear();
    this.cleanup(true);
    this.hitCount = 0;
    this.missCount = 0;
    console.log('Cache completely cleared');
  }
}

/**
 * Rate limiting manager
 */
class RateLimiter {
  constructor(options = {}) {
    this.windowMs = options.windowMs || 60000; // 1 minute
    this.maxRequests = options.maxRequests || 30;
    this.requests = new Map();
    this.cleanupInterval = setInterval(() => this.cleanup(), this.windowMs);
  }

  isAllowed(clientId) {
    const now = Date.now();
    const windowStart = now - this.windowMs;

    if (!this.requests.has(clientId)) {
      this.requests.set(clientId, []);
    }

    const clientRequests = this.requests.get(clientId);
    const recentRequests = clientRequests.filter(time => time > windowStart);
    this.requests.set(clientId, recentRequests);

    if (recentRequests.length >= this.maxRequests) {
      return {
        allowed: false,
        remaining: 0,
        resetTime: windowStart + this.windowMs
      };
    }

    recentRequests.push(now);
    return {
      allowed: true,
      remaining: this.maxRequests - recentRequests.length,
      resetTime: windowStart + this.windowMs
    };
  }

  cleanup() {
    const now = Date.now();
    const windowStart = now - this.windowMs;

    for (const [clientId, requests] of this.requests.entries()) {
      const recentRequests = requests.filter(time => time > windowStart);
      if (recentRequests.length === 0) {
        this.requests.delete(clientId);
      } else {
        this.requests.set(clientId, recentRequests);
      }
    }
  }

  getMiddleware() {
    return (req, res, next) => {
      const clientId = req.ip || req.connection.remoteAddress || 'unknown';
      const result = this.isAllowed(clientId);
      
      // Add rate limiting headers
      res.set({
        'X-RateLimit-Limit': this.maxRequests,
        'X-RateLimit-Remaining': result.remaining,
        'X-RateLimit-Reset': new Date(result.resetTime).toISOString()
      });

      if (!result.allowed) {
        return res.status(429).json({
          error: 'Too Many Requests',
          message: `Limit of ${this.maxRequests} requests per ${this.windowMs/1000} seconds reached`,
          retryAfter: Math.ceil((result.resetTime - Date.now()) / 1000)
        });
      }

      next();
    };
  }

  getStats() {
    return {
      windowMs: this.windowMs,
      maxRequests: this.maxRequests,
      activeClients: this.requests.size,
      totalRequests: Array.from(this.requests.values()).reduce((sum, reqs) => sum + reqs.length, 0)
    };
  }

  destroy() {
    if (this.cleanupInterval) {
      clearInterval(this.cleanupInterval);
    }
  }
}

/**
 * Advanced monitoring system
 */
class MonitoringSystem {
  constructor() {
    this.startTime = Date.now();
    this.metrics = {
      requests: {
        total: 0,
        successful: 0,
        failed: 0,
        byEndpoint: {}
      },
      predictions: {
        '1H': { count: 0, totalTime: 0, errors: 0, avgTime: 0 },
        '13C': { count: 0, totalTime: 0, errors: 0, avgTime: 0 },
        'HSQC': { count: 0, totalTime: 0, errors: 0, avgTime: 0 },
        'COSY': { count: 0, totalTime: 0, errors: 0, avgTime: 0 },
        'NMR_Suite': { count: 0, totalTime: 0, errors: 0, avgTime: 0 }
      },
      systemHealth: {
        memoryUsage: 0,
        cpuUsage: 0,
        diskSpace: 0,
        uptime: 0
      },
      performance: {
        recentResponseTimes: [],
        averageResponseTime: 0,
        p95ResponseTime: 0,
        p99ResponseTime: 0
      }
    };
    
    // Collect system metrics every 30 seconds
    this.metricsInterval = setInterval(() => this.collectSystemMetrics(), 30000);
    
    // Clean old metrics every hour
    this.cleanupInterval = setInterval(() => this.cleanupOldMetrics(), 60 * 60 * 1000);
  }

  recordRequest(endpoint, success = true, responseTime = 0) {
    this.metrics.requests.total++;
    
    if (success) {
      this.metrics.requests.successful++;
    } else {
      this.metrics.requests.failed++;
    }

    // Per endpoint
    if (!this.metrics.requests.byEndpoint[endpoint]) {
      this.metrics.requests.byEndpoint[endpoint] = {
        count: 0,
        successful: 0,
        failed: 0,
        avgResponseTime: 0,
        totalResponseTime: 0
      };
    }

    const endpointMetrics = this.metrics.requests.byEndpoint[endpoint];
    endpointMetrics.count++;
    endpointMetrics.totalResponseTime += responseTime;
    endpointMetrics.avgResponseTime = endpointMetrics.totalResponseTime / endpointMetrics.count;

    if (success) {
      endpointMetrics.successful++;
    } else {
      endpointMetrics.failed++;
    }

    // Overall performance
    this.metrics.performance.recentResponseTimes.push(responseTime);
    if (this.metrics.performance.recentResponseTimes.length > 1000) {
      this.metrics.performance.recentResponseTimes.shift();
    }
    
    this.updatePerformanceMetrics();
  }

  recordPrediction(type, duration, success = true) {
    if (!this.metrics.predictions[type]) {
      this.metrics.predictions[type] = { count: 0, totalTime: 0, errors: 0, avgTime: 0 };
    }

    const predictionMetrics = this.metrics.predictions[type];
    predictionMetrics.count++;
    predictionMetrics.totalTime += duration;
    predictionMetrics.avgTime = Math.round(predictionMetrics.totalTime / predictionMetrics.count);
    
    if (!success) {
      predictionMetrics.errors++;
    }
  }

  updatePerformanceMetrics() {
    const times = this.metrics.performance.recentResponseTimes;
    if (times.length === 0) return;

    const sortedTimes = [...times].sort((a, b) => a - b);
    this.metrics.performance.averageResponseTime = Math.round(times.reduce((a, b) => a + b, 0) / times.length);
    
    const p95Index = Math.floor(sortedTimes.length * 0.95);
    const p99Index = Math.floor(sortedTimes.length * 0.99);
    
    this.metrics.performance.p95ResponseTime = sortedTimes[p95Index] || 0;
    this.metrics.performance.p99ResponseTime = sortedTimes[p99Index] || 0;
  }

  collectSystemMetrics() {
    try {
      const usage = process.memoryUsage();
      this.metrics.systemHealth.memoryUsage = Math.round(usage.heapUsed / 1024 / 1024);
      this.metrics.systemHealth.uptime = Math.round(process.uptime());
      
      // Approximate CPU usage (for precise usage, use dedicated libraries)
      const startUsage = process.cpuUsage();
      setTimeout(() => {
        const endUsage = process.cpuUsage(startUsage);
        const totalUsage = endUsage.user + endUsage.system;
        this.metrics.systemHealth.cpuUsage = Math.round(totalUsage / 1000000 * 100) / 100;
      }, 100);
      
    } catch (error) {
      console.error('[ERROR] System metrics collection error:', error);
    }
  }

  cleanupOldMetrics() {
    // Keep only the last 1000 response times
    if (this.metrics.performance.recentResponseTimes.length > 1000) {
      this.metrics.performance.recentResponseTimes = 
        this.metrics.performance.recentResponseTimes.slice(-1000);
    }
  }

  getMetrics() {
    // Compute error rates
    const predictions = {};
    for (const [type, data] of Object.entries(this.metrics.predictions)) {
      predictions[type] = {
        ...data,
        errorRate: data.count > 0 ? ((data.errors / data.count) * 100).toFixed(1) : 0,
        successRate: data.count > 0 ? (((data.count - data.errors) / data.count) * 100).toFixed(1) : 0
      };
    }

    return {
      ...this.metrics,
      predictions,
      uptime: Math.round((Date.now() - this.startTime) / 1000),
      timestamp: new Date().toISOString()
    };
  }

  getPrometheusMetrics() {
    const metrics = this.getMetrics();
    let output = '';
    
    // Request metrics
    output += `# HELP nmr_requests_total Total number of requests\n`;
    output += `# TYPE nmr_requests_total counter\n`;
    output += `nmr_requests_total{status="success"} ${metrics.requests.successful}\n`;
    output += `nmr_requests_total{status="failed"} ${metrics.requests.failed}\n`;
    
    // Prediction metrics
    output += `# HELP nmr_predictions_total Total number of NMR predictions\n`;
    output += `# TYPE nmr_predictions_total counter\n`;
    
    for (const [type, data] of Object.entries(metrics.predictions)) {
      output += `nmr_predictions_total{type="${type}"} ${data.count}\n`;
      output += `nmr_prediction_duration_ms{type="${type}"} ${data.avgTime}\n`;
      output += `nmr_prediction_errors_total{type="${type}"} ${data.errors}\n`;
    }
    
    // System metrics
    output += `# HELP nmr_memory_usage_mb Memory usage in MB\n`;
    output += `# TYPE nmr_memory_usage_mb gauge\n`;
    output += `nmr_memory_usage_mb ${metrics.systemHealth.memoryUsage}\n`;
    
    output += `# HELP nmr_uptime_seconds Uptime in seconds\n`;
    output += `# TYPE nmr_uptime_seconds gauge\n`;
    output += `nmr_uptime_seconds ${metrics.uptime}\n`;
    
    // Performance metrics
    output += `# HELP nmr_response_time_ms Response time metrics\n`;
    output += `# TYPE nmr_response_time_ms gauge\n`;
    output += `nmr_response_time_ms{quantile="0.5"} ${metrics.performance.averageResponseTime}\n`;
    output += `nmr_response_time_ms{quantile="0.95"} ${metrics.performance.p95ResponseTime}\n`;
    output += `nmr_response_time_ms{quantile="0.99"} ${metrics.performance.p99ResponseTime}\n`;
    
    return output;
  }

  getHealthStatus() {
    const metrics = this.getMetrics();
    const memoryMB = metrics.systemHealth.memoryUsage;
    const errorRate = metrics.requests.total > 0 ? 
      (metrics.requests.failed / metrics.requests.total) * 100 : 0;
    
    let status = 'healthy';
    let issues = [];
    
    if (memoryMB > 1024) { // More than 1GB
      status = 'warning';
      issues.push(`High memory usage: ${memoryMB}MB`);
    }
    
    if (errorRate > 10) {
      status = 'warning';
      issues.push(`High error rate: ${errorRate.toFixed(1)}%`);
    }
    
    if (memoryMB > 2048 || errorRate > 25) {
      status = 'unhealthy';
    }
    
    return {
      status,
      issues,
      metrics: {
        memoryUsage: memoryMB,
        errorRate: errorRate.toFixed(1),
        uptime: metrics.uptime,
        totalRequests: metrics.requests.total
      }
    };
  }

  destroy() {
    if (this.metricsInterval) {
      clearInterval(this.metricsInterval);
    }
    if (this.cleanupInterval) {
      clearInterval(this.cleanupInterval);
    }
  }
}

/**
 * Main optimization manager
 */
class OptimizationManager {
  constructor(options = {}) {
    this.cache = new ModelCache(options.cache);
    this.rateLimiter = new RateLimiter(options.rateLimit);
    this.monitoring = new MonitoringSystem();
    
    console.log('[START] Optimization system initialized');
  }

  // Middleware to measure request performance
  getPerformanceMiddleware() {
    return (req, res, next) => {
      const startTime = Date.now();
      const originalSend = res.send;
      
      res.send = function(data) {
        const duration = Date.now() - startTime;
        const endpoint = req.route ? req.route.path : req.path;
        const success = res.statusCode < 400;
        
        this.monitoring.recordRequest(endpoint, success, duration);
        
        originalSend.call(this, data);
      }.bind(this);
      
      next();
    };
  }

  // Middleware for rate limiting
  getRateLimitMiddleware() {
    return this.rateLimiter.getMiddleware();
  }

  // Check cache for a prediction
  getCachedPrediction(smiles, predictionType, parameters = {}) {
    const cacheKey = this.cache.generateCacheKey(smiles, predictionType, parameters);
    return this.cache.get(cacheKey);
  }

  // Save a prediction in cache
  cachePrediction(smiles, predictionType, result, parameters = {}) {
    const cacheKey = this.cache.generateCacheKey(smiles, predictionType, parameters);
    this.cache.set(cacheKey, result);
    return cacheKey;
  }

  // Record a prediction in metrics
  recordPrediction(type, duration, success = true) {
    this.monitoring.recordPrediction(type, duration, success);
  }

  // Get all stats
  getAllStats() {
    return {
      cache: this.cache.getStats(),
      rateLimit: this.rateLimiter.getStats(),
      monitoring: this.monitoring.getMetrics(),
      health: this.monitoring.getHealthStatus()
    };
  }

  // Clean all data
  cleanup() {
    this.cache.cleanup();
    console.log('[CLEAN] Optimization system cleanup completed');
  }

  // Destroy manager
  destroy() {
    this.rateLimiter.destroy();
    this.monitoring.destroy();
    this.cache.cleanup(true);
    console.log('[STOP] Optimization system stopped');
  }
}

// Exports
module.exports = {
  ModelCache,
  RateLimiter,
  MonitoringSystem,
  OptimizationManager
};

// Clean shutdown on process stop
process.on('SIGINT', () => {
  console.log('\n[STOP] Shutting down...');
  // Instances will be destroyed automatically
  process.exit(0);
});

process.on('SIGTERM', () => {
  console.log('\n[STOP] Shutdown requested...');
  process.exit(0);
});