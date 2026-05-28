/* Laplacian.chpl — Parallel Laplacian construction
 *
 * Builds both unnormalized (L = D - W) and symmetric normalized
 * (L = I - D^{-1/2} W D^{-1/2}) Laplacians using forall parallelism
 * and Chapel's reduce operations.
 */

module Laplacian {
  use TensionGraph;

  /* Laplacian matrix stored as a dense 2D Chapel array */
  record LaplacianMatrix {
    var n: int;
    var normalized: bool;
    var matrix: [0..#n, 0..#n] real;

    proc init(n: int, normalized: bool = false) {
      this.n = n;
      this.normalized = normalized;
    }
  }

  /* Build the Laplacian from a TensionGraph.
   *
   * Uses forall loops for parallel degree computation and matrix fill.
   * Degree is computed via Chapel's + reduce over each row.
   */
  proc buildLaplacian(graph: borrowed TensionGraph, normalized: bool = false): LaplacianMatrix {
    var lap = new LaplacianMatrix(n=graph.n, normalized=normalized);

    // Parallel degree computation: one reduction per row
    var degree: [0..#graph.n] real;
    forall i in 0..#graph.n with (ref degree) do
      degree[i] = + reduce graph.adjMatrix[i, ..];

    if !normalized {
      // L = D - W (unnormalized combinatorial Laplacian)
      forall (i, j) in {0..#graph.n, 0..#graph.n} with (ref lap) {
        if i == j then
          lap.matrix[i, j] = degree[i];
        else
          lap.matrix[i, j] = -graph.adjMatrix[i, j];
      }
    } else {
      // Symmetric normalized: L = I - D^{-1/2} W D^{-1/2}
      forall (i, j) in {0..#graph.n, 0..#graph.n} with (ref lap) {
        var diSqrt: real = if degree[i] > 0.0 then 1.0 / sqrt(degree[i]) else 0.0;
        var djSqrt: real = if degree[j] > 0.0 then 1.0 / sqrt(degree[j]) else 0.0;
        var wNorm = diSqrt * graph.adjMatrix[i, j] * djSqrt;
        if i == j then
          lap.matrix[i, j] = 1.0 - wNorm;
        else
          lap.matrix[i, j] = -wNorm;
      }
    }

    return lap;
  }

  /* Distributed Laplacian construction for large graphs.
   * Uses Block distribution to partition the matrix across locales.
   * Requires: use BlockDist;
   */
  proc buildDistributedLaplacian(graph: borrowed TensionGraph, normalized: bool = false) {
    use BlockDist;

    var n = graph.n;
    var dom2D = Block.createDomain({0..#n, 0..#n});
    var lapMatrix: [dom2D] real;

    // Build local degree array (parallel on each locale's portion)
    var degree: [0..#n] real;
    forall i in 0..#n with (ref degree) do
      degree[i] = + reduce graph.adjMatrix[i, ..];

    // Parallel fill — each locale handles its block
    forall (i, j) in dom2D with (ref lapMatrix) {
      if i == j then
        lapMatrix[i, j] = degree[i];
      else
        lapMatrix[i, j] = -graph.adjMatrix[i, j];
    }

    return (lapMatrix, degree);
  }
}
