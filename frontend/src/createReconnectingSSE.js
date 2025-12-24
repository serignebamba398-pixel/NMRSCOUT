// src/createReconnectingSSE.js
import { API_BASE_URL } from "./api";

const RECONNECT_DELAY = 2500; // ms
const IDLE_CHECK_MS   = 8000; // if no data received for 8s, check status

// No hard timeout: 0 => never automatically cut a long task
const HARD_TIMEOUT_MS = 0;

export function createReconnectingSSE({
  endpoint,
  body,
  onProgress,
  onLog,
  onError,
  onComplete,
}) {
  let abortController = new AbortController();
  let taskId = null;
  let hasCompleted = false;

  let reconnectTimeout = null;
  let statusInterval   = null;
  let idleTimer        = null;
  let hardTimeout      = null;

  const stopAllTimers = () => {
    if (reconnectTimeout) { clearTimeout(reconnectTimeout); reconnectTimeout = null; }
    if (statusInterval)   { clearInterval(statusInterval);   statusInterval   = null; }
    if (idleTimer)        { clearTimeout(idleTimer);         idleTimer        = null; }
    if (hardTimeout)      { clearTimeout(hardTimeout);       hardTimeout      = null; }
  };

  const armIdleTimer = () => {
    if (idleTimer) clearTimeout(idleTimer);
    idleTimer = setTimeout(async () => {
      // No events for a while → force a check.
      await checkTaskStatus(true);
    }, IDLE_CHECK_MS);
  };

  const cancel = () => {
    hasCompleted = true;
    abortController.abort();
    stopAllTimers();
  };

  const getTaskId = () => taskId;

  const finalizeSuccess = (payload = {}) => {
    if (hasCompleted) return;
    hasCompleted = true;
    stopAllTimers();
    if (onComplete) onComplete({ success: true, ...payload, taskId });
  };

  const finalizeError = (message = "Unknown error") => {
    if (hasCompleted) return;
    hasCompleted = true;
    stopAllTimers();
    if (onError) onError(message);
  };

  const connect = async () => {
    if (hasCompleted) return;

    try {
      const res = await fetch(`${API_BASE_URL}${endpoint}`, {
        method: "POST",
        headers: {
          "Content-Type": "application/json",
          Accept: "text/event-stream",
        },
        body: JSON.stringify(body),
        signal: abortController.signal,
      });

      if (!res.ok) throw new Error(`HTTP ${res.status}`);

      const reader  = res.body.getReader();
      const decoder = new TextDecoder();
      let buf = "";

      // (Re)arm timers
      armIdleTimer();

      // Only arm hard timeout if enabled (>0)
      if (!hardTimeout && HARD_TIMEOUT_MS && HARD_TIMEOUT_MS > 0) {
        hardTimeout = setTimeout(
          () => finalizeError(`Timeout exceeded (${HARD_TIMEOUT_MS / 60000} min)`),
          HARD_TIMEOUT_MS
        );
      }

      const pump = async () => {
        if (hasCompleted) return;

        try {
          const { done, value } = await reader.read();

          if (done) {
            // Stream closed server-side → retry reconnection (fallback)
            if (!hasCompleted) {
              reconnectTimeout = setTimeout(() => connect(), RECONNECT_DELAY);
            }
            return;
          }

          // Received bytes → reset idle timer
          armIdleTimer();

          buf += decoder.decode(value, { stream: true });
          const lines = buf.split("\n");
          buf = lines.pop() || "";

          for (const raw of lines) {
            const line = raw.trim();
            if (!line.startsWith("data:")) continue;

            try {
              const data = JSON.parse(line.slice(5).trim());

              switch (data.type) {
                case "taskId":
                  taskId = data.taskId;
                  break;

                case "progress":
                  onProgress && onProgress(data);
                  break;

                case "log":
                  onLog && onLog(data.message);
                  break;

                case "error":
                  finalizeError(data.message || "Unknown error");
                  return;

                case "completed":
                  finalizeSuccess(data);
                  return;

                default:
                  // silently ignore
                  break;
              }
            } catch (parseErr) {
              console.error("SSE parse error:", parseErr);
            }
          }

          return pump();
        } catch (e) {
          if (e.name !== "AbortError" && !hasCompleted) {
            // Network issue → retry
            reconnectTimeout = setTimeout(() => connect(), RECONNECT_DELAY);
          }
        }
      };

      await pump();
    } catch (e) {
      if (e.name === "AbortError" || hasCompleted) return;
      reconnectTimeout = setTimeout(() => connect(), RECONNECT_DELAY);
    }
  };

  // Periodic status check (backup outside SSE)
  const checkTaskStatus = async (fromIdle = false) => {
    if (hasCompleted) return;

    try {
      // If we don't have a taskId yet and the call comes from idle timer,
      // we can't query the status API → force a reconnection.
      if (!taskId) {
        if (fromIdle && !hasCompleted) {
          if (onLog) onLog("No taskId received — forced reconnection...");
          abortController.abort();
          abortController = new AbortController();
          reconnectTimeout = setTimeout(() => connect(), 50);
        }
        return;
      }

      // [NOTE] Important: API_BASE_URL already includes "/api" if needed
      const r = await fetch(`${API_BASE_URL}/tasks/${taskId}`);
      if (!r.ok) return;

      const data = await r.json();
      if (!data) return;

      if ((data.status === "completed" && !data.pid) || data.status === "success") {
        finalizeSuccess({
          downloadUrl: data.downloadUrl,
          type: data.type,
        });
      } else if (data.status === "failed" || data.status === "error") {
        finalizeError(data.error || "Task failed");
      }
    } catch (err) {
      // Log it, let the next iteration try its chance
      console.error("Error checking task status:", err);
    }
  };

  // Startup
  connect();
  statusInterval = setInterval(checkTaskStatus, 3000);

  return {
    cancel: () => cancel(),
    cancelTask: async () => {
      if (!taskId) return;
      try {
        const { cancelTask } = await import("./api");
        await cancelTask(taskId);
      } catch (e) {
        console.error("Error cancelTask:", e);
      }
    },
    getTaskId,
  };
}
