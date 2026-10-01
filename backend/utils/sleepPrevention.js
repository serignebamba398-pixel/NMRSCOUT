/**
 * Sleep Prevention Utility
 * 
 * Keep Windows awake during long-running tasks.
 * Request wakefulness through PowerShell.
 */

const { exec, spawn } = require('child_process');
const { logWithMemory } = require('./logger');

let preventSleepProcess = null;
let activeTaskCount = 0;

/**
 * PowerShell script for sleep prevention.
 * Use SetThreadExecutionState to tell Windows
 * to keep the system awake.
 */
const PREVENT_SLEEP_SCRIPT = `
Add-Type -TypeDefinition @"
using System;
using System.Runtime.InteropServices;
public class SleepPreventer {
    [DllImport("kernel32.dll", CharSet = CharSet.Auto, SetLastError = true)]
    public static extern uint SetThreadExecutionState(uint esFlags);
    
    public const uint ES_CONTINUOUS = 0x80000000;
    public const uint ES_SYSTEM_REQUIRED = 0x00000001;
    public const uint ES_DISPLAY_REQUIRED = 0x00000002;
    
    public static void PreventSleep() {
        SetThreadExecutionState(ES_CONTINUOUS | ES_SYSTEM_REQUIRED);
    }
    
    public static void AllowSleep() {
        SetThreadExecutionState(ES_CONTINUOUS);
    }
}
"@

# Prevent sleep.
[SleepPreventer]::PreventSleep()
Write-Host "Sleep prevention ENABLED"

# Wait until the process is terminated.
while ($true) {
    Start-Sleep -Seconds 60
    # Refresh the signal every minute.
    [SleepPreventer]::PreventSleep()
}
`;

/**
 * Start sleep prevention.
 * Called when a long-running task starts.
 */
function startPreventSleep(taskId) {
    if (process.platform !== 'win32') {
        logWithMemory(`[${taskId}] Sleep prevention not available on this platform`, 'DEBUG');
        return;
    }

    activeTaskCount++;
    logWithMemory(`[${taskId}] Task started - active tasks: ${activeTaskCount}`, 'DEBUG');

    if (preventSleepProcess) {
        // Already active.
        logWithMemory(`[${taskId}] Sleep prevention already active`, 'DEBUG');
        return;
    }

    try {
        preventSleepProcess = spawn('powershell.exe', ['-NoProfile', '-Command', PREVENT_SLEEP_SCRIPT], {
            stdio: ['ignore', 'pipe', 'pipe'],
            detached: false,
            windowsHide: true,
        });

        preventSleepProcess.stdout.on('data', (data) => {
            logWithMemory(`[SleepPreventer] ${data.toString().trim()}`, 'INFO');
        });

        preventSleepProcess.stderr.on('data', (data) => {
            logWithMemory(`[SleepPreventer] Error: ${data.toString().trim()}`, 'WARN');
        });

        preventSleepProcess.on('close', (code) => {
            logWithMemory(`[SleepPreventer] Process exited with code ${code}`, 'DEBUG');
            preventSleepProcess = null;
        });

        preventSleepProcess.on('error', (err) => {
            logWithMemory(`[SleepPreventer] Spawn error: ${err.message}`, 'WARN');
            preventSleepProcess = null;
        });

        logWithMemory(`[${taskId}] Sleep prevention started (PID: ${preventSleepProcess.pid})`, 'INFO');
    } catch (err) {
        logWithMemory(`[${taskId}] Failed to start sleep prevention: ${err.message}`, 'WARN');
    }
}

/**
 * Stop sleep prevention.
 * Called when a long-running task ends.
 */
function stopPreventSleep(taskId) {
    if (process.platform !== 'win32') {
        return;
    }

    activeTaskCount = Math.max(0, activeTaskCount - 1);
    logWithMemory(`[${taskId}] Task finished - active tasks: ${activeTaskCount}`, 'DEBUG');

    if (activeTaskCount > 0) {
        // Other tasks still need sleep prevention.
        logWithMemory(`[${taskId}] Other tasks still running, keeping sleep prevention`, 'DEBUG');
        return;
    }

    if (preventSleepProcess) {
        try {
            preventSleepProcess.kill('SIGTERM');
            logWithMemory(`[${taskId}] Sleep prevention stopped`, 'INFO');
        } catch (err) {
            logWithMemory(`[${taskId}] Error stopping sleep prevention: ${err.message}`, 'WARN');
        }
        preventSleepProcess = null;
    }

    // Release the system wakefulness request.
    exec('powershell.exe -NoProfile -Command "[SleepPreventer]::AllowSleep()" 2>$null', () => {});
}

/**
 * Force cleanup of sleep prevention.
 */
function forceStopPreventSleep() {
    activeTaskCount = 0;
    if (preventSleepProcess) {
        try {
            preventSleepProcess.kill('SIGKILL');
        } catch {}
        preventSleepProcess = null;
    }
}

/**
 * Check whether sleep prevention is active.
 */
function isPreventingSleep() {
    return preventSleepProcess !== null && activeTaskCount > 0;
}

module.exports = {
    startPreventSleep,
    stopPreventSleep,
    forceStopPreventSleep,
    isPreventingSleep,
};
