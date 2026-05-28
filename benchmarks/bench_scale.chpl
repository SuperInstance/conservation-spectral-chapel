/* bench_scale.chpl — Benchmark scaling of parallel Laplacian construction
 *
 * Generates random graphs of increasing size and benchmarks
 * Laplacian construction and eigendecomposition using Chapel's timers.
 */

use ConservationSpectral;
use IO;
use Random;
use Time;

config const minSize = 50;
config const maxSize = 500;
config const step = 50;
config const edgeDensity = 0.3;
config const seed = 42;

proc main() {
  writeln("=== Conservation Spectral SDK — Scaling Benchmark ===");
  writeln();
  writeln("  Graph sizes: ", minSize, " to ", maxSize, " (step ", step, ")");
  writeln("  Edge density: ", edgeDensity);
  writeln();

  writeln("| n     | edges | Laplacian (ms) | Eigendecomp (ms) | Total (ms) |");
  writeln("|-------|-------|----------------|-------------------|------------|");

  var rng = new randomStream(real, seed);

  for size in minSize..maxSize by step {
    // Build random graph
    var graph = createGraph(size);

    // Random vertex attributes
    for i in 0..#size do
      graph.addVertex(i, rng.getNext());

    // Random edges based on density
    for i in 0..#size {
      for j in (i + 1)..#(size - i - 1) {
        if rng.getNext() < edgeDensity then
          graph.addEdge(i, j, rng.getNext());
      }
    }

    // Benchmark Laplacian
    var tLap = new Timer();
    tLap.start();
    var lap = buildLaplacian(graph, normalized=false);
    tLap.stop();

    // Benchmark eigendecomposition (limit k for large graphs)
    var k = if size <= 100 then size else 20;
    var tEigen = new Timer();
    tEigen.start();
    var eigen = eigendecompose(lap, k);
    tEigen.stop();

    var totalMs = tLap.elapsed() * 1000 + tEigen.elapsed() * 1000;
    writeln("| ", size, " | ", graph.numEdges(), " | ",
            (tLap.elapsed() * 1000):int, " | ",
            (tEigen.elapsed() * 1000):int, " | ",
            totalMs:int, " |");
  }

  writeln();
  writeln("=== Benchmark complete ===");
}
