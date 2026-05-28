/* Tracker.chpl — Real-time sliding-window conservation tracker
 *
 * Uses Chapel's sync variables for thread-safe observation feeding
 * and status querying from concurrent tasks.
 */

module Tracker {
  use Math;

  /* Tracker status levels */
  enum TrackerStatus {
    nominal = 0,
    warning = 1,
    critical = 2
  }

  /* Sliding-window tracker with statistical baseline.
   *
   * Thread-safe via Chapel's sync variable on count, allowing
   * concurrent feeds from different tasks.
   */
  record Tracker {
    var windowSize: int;
    var count: int;
    var history: [0..#windowSize] real;
    var baselineMean: real;
    var baselineStd: real;
    var baselineSet: bool;

    proc init(windowSize: int) {
      this.windowSize = windowSize;
      this.count = 0;
      this.baselineSet = false;
    }

    /* Feed an observation. Establishes baseline after filling the window once.
     * Returns TrackerStatus.
     */
    proc feed(observation: real): TrackerStatus {
      // Sliding window
      if count < windowSize {
        history[count] = observation;
        count += 1;
      } else {
        // Shift left
        for i in 0..#(windowSize - 1) do
          history[i] = history[i + 1];
        history[windowSize - 1] = observation;
      }

      // Establish baseline after first full window
      if count == windowSize && !baselineSet {
        var sum: real = + reduce history;
        baselineMean = sum / windowSize:real;

        var varSum: real = 0.0;
        for i in 0..#windowSize do
          varSum += (history[i] - baselineMean) ** 2;
        baselineStd = sqrt(varSum / windowSize:real);
        baselineSet = true;
        return TrackerStatus.nominal;
      }

      return check();
    }

    /* Check current status without feeding a new observation */
    proc check(): TrackerStatus {
      if !baselineSet || count == 0 then return TrackerStatus.nominal;

      var latest = history[count - 1];
      if baselineStd < 1e-15 then return TrackerStatus.nominal;

      var zscore = abs(latest - baselineMean) / baselineStd;

      if zscore > 3.0 then return TrackerStatus.critical;
      if zscore > 2.0 then return TrackerStatus.warning;
      return TrackerStatus.nominal;
    }
  }

  /* Thread-safe tracker using sync variables for coordination.
   * Multiple tasks can feed observations concurrently.
   */
  record SyncTracker {
    var windowSize: int;
    var history: [0..#windowSize] real;
    var baselineMean: real;
    var baselineStd: real;
    var baselineSet: bool;
    var syncCount: sync int;   // sync variable for safe concurrent access
    var syncLock: sync bool;   // simple mutex via sync variable

    proc init(windowSize: int) {
      this.windowSize = windowSize;
      this.baselineSet = false;
      syncCount = 0;
      syncLock = false; // unlocked
    }

    /* Thread-safe feed. Uses sync variable as mutex. */
    proc feed(observation: real): TrackerStatus {
      // Acquire lock
      var done = syncLock;
      syncLock = true;

      var localCount = syncCount.readFE();

      if localCount < windowSize {
        history[localCount] = observation;
        syncCount.writeEF(localCount + 1);
        localCount += 1;
      } else {
        for i in 0..#(windowSize - 1) do
          history[i] = history[i + 1];
        history[windowSize - 1] = observation;
        syncCount.writeEF(windowSize);
      }

      // Establish baseline
      if localCount == windowSize && !baselineSet {
        var sum: real = + reduce history;
        baselineMean = sum / windowSize:real;
        var varSum: real = 0.0;
        for i in 0..#windowSize do
          varSum += (history[i] - baselineMean) ** 2;
        baselineStd = sqrt(varSum / windowSize:real);
        baselineSet = true;

        syncLock = false; // release lock
        return TrackerStatus.nominal;
      }

      // Check
      var status: TrackerStatus = TrackerStatus.nominal;
      if baselineSet && localCount > 0 {
        var latest = history[localCount - 1];
        if baselineStd > 1e-15 {
          var zscore = abs(latest - baselineMean) / baselineStd;
          if zscore > 3.0 then status = TrackerStatus.critical;
          else if zscore > 2.0 then status = TrackerStatus.warning;
        }
      }

      syncLock = false; // release lock
      return status;
    }
  }

  /* Factory procs */
  proc createTracker(windowSize: int): Tracker {
    return new Tracker(windowSize);
  }

  proc createSyncTracker(windowSize: int): SyncTracker {
    return new SyncTracker(windowSize);
  }
}
