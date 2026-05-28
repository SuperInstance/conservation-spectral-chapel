# Conservation Spectral SDK — Chapel

[![Chapel](https://img.shields.io/badge/language-Chapel-orange.svg)](https://chapel-lang.org/)
[![Version](https://img.shields.io/badge/version-0.1.0-blue.svg)](https://github.com/SuperInstance/conservation-spectral-chapel)

A Chapel implementation of the Conservation Spectral SDK, leveraging Chapel's native parallelism features — `forall`, `coforall`, reductions, domain-based arrays, locale-aware distributions, and sync variables — for high-performance spectral graph analysis.

## Why Chapel?

Chapel (from Cray/HPE) excels at:
- **Data parallelism**: `forall` loops over domains automatically parallelize
- **Task parallelism**: `coforall` spawns concurrent tasks
- **Locality-aware computing**: Block distributions partition data across multi-locale systems
- **Reductions**: Built-in `+ reduce`, `min reduce`, etc. for parallel aggregation
- **Sync variables**: Lightweight coordination between concurrent tasks

This makes Chapel a natural fit for spectral graph algorithms where Laplacian construction, eigenvector computation, and conservation analysis all benefit from massive parallelism.

## Architecture

```
conservation-spectral-chapel/
├── Makefile
├── src/
│   ├── ConservationSpectral.chpl    # Main module (re-exports all)
│   ├── TensionGraph.chpl            # Graph with domain-based arrays
│   ├── Laplacian.chpl               # Parallel Laplacian construction
│   ├── Eigen.chpl                   # Eigendecomposition (power iteration)
│   ├── Conservation.chpl            # Conservation analysis + fingerprinting
│   └── Tracker.chpl                 # Real-time sliding-window tracker
├── tests/
│   └── test_chords.chpl             # Test suite (chord graph)
├── benchmarks/
│   └── bench_scale.chpl             # Scaling benchmarks
└── README.md
```

## Key Chapel Features Used

### 1. Domains and Arrays
All matrices use Chapel's domain-based arrays:
```chapel
var adjMatrix: [0..#n, 0..#n] real;  // Dense n×n matrix
```

### 2. forall Loops (Data Parallelism)
Parallel Laplacian construction:
```chapel
forall (i, j) in {0..#n, 0..#n} with (ref lap) {
  if i == j then lap.matrix[i, j] = degree[i];
  else lap.matrix[i, j] = -adjMatrix[i, j];
}
```

### 3. Reductions
Parallel degree computation:
```chapel
forall i in 0..#n with (ref degree) do
  degree[i] = + reduce adjMatrix[i, ..];
```

### 4. Sync Variables (Thread Coordination)
The `SyncTracker` uses Chapel's sync variables for concurrent observation feeding:
```chapel
var syncCount: sync int;   // Full/empty semantics
var syncLock: sync bool;   // Mutex via sync variable
```

### 5. Locale-Aware Distribution
For large graphs spanning multiple nodes:
```chapel
use BlockDist;
var dom2D = Block.createDomain({0..#n, 0..#n});
var distMatrix: [dom2D] real;
```

## Building

### Prerequisites
- [Chapel compiler](https://chapel-lang.org/download.html) (`chpl` 1.30+)

### Compile
```bash
make                # Build library
make test           # Build and run tests
make bench          # Build and run benchmarks
```

### Manual Compilation
```bash
# Compile tests
chpl -o test_chords src/*.chpl tests/test_chords.chpl
./test_chords

# Compile benchmarks
chpl -o bench_scale src/*.chpl benchmarks/bench_scale.chpl
./bench_scale --minSize=50 --maxSize=500 --step=50

# Multi-locale (requires Chapel built with GASNet)
chpl -nl 4 -o bench_scale_dist src/*.chpl benchmarks/bench_scale.chpl
./bench_scale_dist -nl 4
```

## API Overview

### Graph Construction
```chapel
var graph = createGraph(8);
graph.addVertex(0, 0.5);      // vertex id, attribute
graph.addEdge(0, 1, 1.0);     // from, to, weight
```

### Laplacian
```chapel
var lap = buildLaplacian(graph, normalized=false);
// Or normalized:
var lapNorm = buildLaplacian(graph, normalized=true);
// Or distributed for large graphs:
var (distLap, degree) = buildDistributedLaplacian(graph);
```

### Eigendecomposition
```chapel
var eigen = eigendecompose(lap, k=10);  // top-k or all if k=0
var gap = spectralGap(eigen);
```

### Conservation Analysis
```chapel
// Full pipeline
var report = csAnalyze(graph);
writeln(report.spectralGap);
writeln(report.cheegerConstant);
writeln(report.anomalyCount);

// Individual operations
var ratio = conservationRatio(eigen, graph.vertexAttrs, 1);
var cheeger = cheegerConstant(lap, fiedlerVector);
var fp = computeFingerprint(eigen.eigenvalues, n);
```

### Tracker
```chapel
var tracker = createTracker(windowSize=10);
var status = tracker.feed(42.0);  // TrackerStatus.nominal/warning/critical

// Thread-safe version for concurrent tasks
var syncTracker = createSyncTracker(10);
coforall tid in 0..#numTasks {
  syncTracker.feed(observations[tid]);
}
```

## Running Tests

```bash
make test
```

Expected output:
```
=== Conservation Spectral SDK — Chapel Test Suite ===

Test 1: TensionGraph construction
  Vertices: 8
  Edges: 12
  PASS

Test 2: Laplacian construction
  Max row sum (should be ~0): 0.0
  PASS

Test 3: Eigendecomposition
  λ[0] ≈ 0 ✓
  PASS

...

=== All tests complete ===
```

## Running Benchmarks

```bash
make bench
# Or with custom parameters:
./bench_scale --minSize=100 --maxSize=1000 --step=100
```

## Reference Implementation

This is the Chapel port of the [C reference implementation](https://github.com/SuperInstance/conservation-spectral-c). The API mirrors the C version while leveraging Chapel's native parallelism.

## License

MIT
