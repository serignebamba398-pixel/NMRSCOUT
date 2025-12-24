let logs = [];
const listeners = new Set();

const emit = () => {
  for (const listener of listeners) {
    try {
      listener();
    } catch (_) {
      // ignore listener errors
    }
  }
};

export const logStore = {
  subscribe(listener) {
    listeners.add(listener);
    return () => listeners.delete(listener);
  },

  getSnapshot() {
    return logs;
  },

  addLog(entry) {
    logs = [...logs, entry];
    emit();
  },

  clear() {
    if (logs.length === 0) return;
    logs = [];
    emit();
  },
};
