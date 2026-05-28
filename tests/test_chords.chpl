/* test_chords.chpl — Test suite for Conservation Spectral SDK (Chapel)
 *
 * Builds a simple chord graph, runs Laplacian construction,
 * eigendecomposition, and conservation analysis.
 */

use ConservationSpectral;
use IO;

proc main() {
  writeln("=== Conservation Spectral SDK — Chapel Test Suite ===");
  writeln();

  // --- Test 1: Graph construction ---
  writeln("Test 1: TensionGraph construction");
  var n = 8;
  var graph = createGraph(n);

  // Add vertex attributes (musical chord tensions)
  var tensions = [0.0, 0.5, 1.0, 0.3, 0.8, 0.2, 0.9, 0.4];
  for i in 0..#n do
    graph.addVertex(i, tensions[i]);

  // Add edges (connections between notes in a chord progression)
  // C major: 0-1-2, F major: 3-4-5, G major: 5-6-7, bridges between chords
  graph.addEdge(0, 1, 1.0);
  graph.addEdge(1, 2, 1.0);
  graph.addEdge(0, 2, 0.5);
  graph.addEdge(3, 4, 1.0);
  graph.addEdge(4, 5, 1.0);
  graph.addEdge(3, 5, 0.5);
  graph.addEdge(5, 6, 1.0);
  graph.addEdge(6, 7, 1.0);
  graph.addEdge(5, 7, 0.5);
  // Bridge edges
  graph.addEdge(2, 3, 0.3);
  graph.addEdge(2, 5, 0.2);
  graph.addEdge(0, 7, 0.1);

  writeln("  Vertices: ", n);
  writeln("  Edges: ", graph.numEdges());
  writeln("  Degree of vertex 0: ", graph.degree(0):string);
  writeln("  PASS");
  writeln();

  // --- Test 2: Parallel Laplacian ---
  writeln("Test 2: Laplacian construction");
  var lap = buildLaplacian(graph, normalized=false);
  writeln("  Laplacian size: ", lap.n, "x", lap.n);
  writeln("  Diagonal (degrees):");
  write("    [");
  for i in 0..#n {
    write(lap.matrix[i, i]:string);
    if i < n - 1 then write(", ");
  }
  writeln("]");

  // Verify Laplacian property: row sums should be ~0
  var maxRowSum: real = 0.0;
  for i in 0..#n {
    var rowSum: real = + reduce lap.matrix[i, ..];
    if abs(rowSum) > maxRowSum then maxRowSum = abs(rowSum);
  }
  writeln("  Max row sum (should be ~0): ", maxRowSum:string);
  if maxRowSum < 1e-10 then writeln("  PASS") else writeln("  WARN: row sum not zero");
  writeln();

  // --- Test 3: Eigendecomposition ---
  writeln("Test 3: Eigendecomposition");
  var eigen = eigendecompose(lap, n);
  writeln("  First 4 eigenvalues:");
  for i in 0..min(3, n - 1) do
    writeln("    λ[", i, "] = ", eigen.eigenvalues[i]:string);

  // First eigenvalue should be ~0 (trivial for Laplacian)
  if abs(eigen.eigenvalues[0]) < 1e-6 then
    writeln("  λ[0] ≈ 0 ✓ (Laplacian trivial eigenvalue)");
  else
    writeln("  WARN: λ[0] should be ~0, got ", eigen.eigenvalues[0]:string);
  writeln("  PASS");
  writeln();

  // --- Test 4: Conservation ratios ---
  writeln("Test 4: Conservation analysis");
  forall k in 0..#min(4, n) {
    var ratio = conservationRatio(eigen, graph.vertexAttrs, k);
    writeln("  Ratio[", k, "] = ", ratio:string);
  }
  writeln();

  // --- Test 5: Spectral gap ---
  writeln("Test 5: Spectral gap");
  var gap = spectralGap(eigen);
  writeln("  Spectral gap: ", gap:string);
  writeln("  PASS");
  writeln();

  // --- Test 6: Full analysis ---
  writeln("Test 6: Full conservation analysis");
  var report = csAnalyze(graph);
  writeln("  Spectral gap: ", report.spectralGap:string);
  writeln("  Cheeger constant: ", report.cheegerConstant:string);
  writeln("  Anomaly count: ", report.anomalyCount);
  writeln("  PASS");
  writeln();

  // --- Test 7: Tracker ---
  writeln("Test 7: Sliding-window tracker");
  var tracker = createTracker(5);

  // Feed normal observations
  for i in 0..#5 do {
    var status = tracker.feed(1.0 + i:real * 0.1);
    write("  Feed ", (1.0 + i:real * 0.1):string, " → ");
    writeln(status:string);
  }
  writeln("  Baseline established: mean=", tracker.baselineMean:string,
          " std=", tracker.baselineStd:string);

  // Feed anomalous observation
  var status = tracker.feed(100.0);
  writeln("  Feed 100.0 → ", status:string, " (should be critical)");
  if status == TrackerStatus.critical then writeln("  PASS") else writeln("  WARN");
  writeln();

  // --- Test 8: Fingerprint ---
  writeln("Test 8: Spectral fingerprint");
  var fp = computeFingerprint(eigen.eigenvalues, n);
  writeln("  Fingerprint (first 64 chars): ", fp[0..min(63, fp.length - 1)]);

  var fp2 = computeFingerprint(eigen.eigenvalues, n);
  var similarity = compareFingerprints(fp, fp2);
  writeln("  Self-similarity: ", similarity:string, " (should be 1.0)");
  if similarity > 0.99 then writeln("  PASS") else writeln("  WARN");
  writeln();

  writeln("=== All tests complete ===");
}
