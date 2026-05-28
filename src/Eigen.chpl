/* Eigen.chpl — Parallel eigendecomposition via power iteration + deflation
 *
 * Computes k eigenvalues/eigenvectors of a Laplacian matrix using
 * shifted power iteration with deflation. Uses coforall for independent
 * eigenvector computations and forall for parallel vector operations.
 */

module Eigen {
  use Laplacian;
  use Math;

  /* Eigen decomposition result */
  record EigenDecomp {
    var n: int;
    var k: int;  // number of eigenpairs computed
    var eigenvalues: [0..#n] real;
    var eigenvectors: [0..#n, 0..#n] real;  // column i = eigenvector i

    proc init(n: int, k: int = 0) {
      this.n = n;
      this.k = if k <= 0 then n else k;
    }
  }

  /* Dense matrix-vector multiply: y = A * x */
  proc matvec(A: [] real, x: [] real, ref y: [] real, n: int) {
    forall i in 0..#n with (ref y) {
      var sum: real = 0.0;
      for j in 0..#n do
        sum += A[i, j] * x[j];
      y[i] = sum;
    }
  }

  /* Compute vector norm */
  proc vecNorm(v: [] real, n: int): real {
    return sqrt(+ reduce (v[0..#n] * v[0..#n]));
  }

  /* Normalize a vector in place, return the norm */
  proc vecNormalize(ref v: [] real, n: int): real {
    var nm = vecNorm(v, n);
    if nm > 1e-30 then
      forall i in 0..#n with (ref v) do v[i] /= nm;
    return nm;
  }

  /* Eigendecompose a Laplacian matrix.
   *
   * Uses shifted power iteration: computes largest eigenvalues of M = shift*I - L,
   * then converts back to eigenvalues of L.
   * Deflation removes found eigenvectors to find subsequent ones.
   *
   * The outer loop over eigenvectors is sequential (deflation dependency),
   * but inner vector operations use forall parallelism.
   */
  proc eigendecompose(lap: borrowed LaplacianMatrix, k: int = 0): EigenDecomp {
    var n = lap.n;
    var kk = if k <= 0 || k > n then n else k;
    var eigen = new EigenDecomp(n=n, k=kk);

    // Find shift = max diagonal (upper bound on max eigenvalue)
    var shift: real = 0.0;
    for i in 0..#n do
      if lap.matrix[i, i] > shift then shift = lap.matrix[i, i];

    // Build M = shift*I - L
    var M: [0..#n, 0..#n] real;
    forall (i, j) in {0..#n, 0..#n} with (ref M) {
      if i == j then
        M[i, j] = shift - lap.matrix[i, j];
      else
        M[i, j] = -lap.matrix[i, j];
    }

    // Residual matrix for deflation
    var R: [0..#n, 0..#n] real = M;

    // Temp vectors
    var v: [0..#n] real;
    var w: [0..#n] real;

    for ev in 0..#kk {
      // Seed vector (varied to avoid convergence to same eigenvector)
      forall i in 0..#n do
        v[i] = 1.0 / (i + 1 + ev):real;

      var lambda: real = 0.0;
      const maxIter = 2000;
      const tol = 1e-12;

      // Power iteration
      for iter in 0..#maxIter {
        // w = R * v
        matvec(R, v, w, n);
        var nm = vecNormalize(w, n);
        if nm < 1e-30 then break;
        forall i in 0..#n do v[i] = w[i];

        // Rayleigh quotient
        matvec(R, v, w, n);
        var rq: real = 0.0;
        for i in 0..#n do rq += v[i] * w[i];

        if abs(rq - lambda) < tol {
          lambda = rq;
          break;
        }
        lambda = rq;
      }

      // Store: eigenvalue of L = shift - lambda_M
      eigen.eigenvalues[ev] = shift - lambda;

      // Store eigenvector (column ev)
      forall i in 0..#n do
        eigen.eigenvectors[i, ev] = v[i];

      // Deflate: R = R - lambda * v * v^T
      forall (i, j) in {0..#n, 0..#n} with (ref R) {
        R[i, j] -= lambda * v[i] * v[j];
      }
    }

    // Sort eigenvalues ascending (and their corresponding eigenvectors)
    for i in 0..#(n - 1) {
      for j in (i + 1)..#(n - i - 1) {
        if eigen.eigenvalues[j] < eigen.eigenvalues[i] {
          // Swap eigenvalues
          var tmp = eigen.eigenvalues[i];
          eigen.eigenvalues[i] = eigen.eigenvalues[j];
          eigen.eigenvalues[j] = tmp;
          // Swap eigenvector columns
          forall r in 0..#n {
            var tv = eigen.eigenvectors[r, i];
            eigen.eigenvectors[r, i] = eigen.eigenvectors[r, j];
            eigen.eigenvectors[r, j] = tv;
          }
        }
      }
    }

    return eigen;
  }

  /* Spectral gap: largest gap between consecutive eigenvalues */
  proc spectralGap(eigen: borrowed EigenDecomp): real {
    if eigen.n < 2 then return 0.0;
    var maxGap: real = 0.0;
    for i in 0..#(eigen.n - 1) {
      var gap = eigen.eigenvalues[i + 1] - eigen.eigenvalues[i];
      if gap > maxGap then maxGap = gap;
    }
    return maxGap;
  }
}
