/* Conservation.chpl — Conservation analysis, Cheeger constant, fingerprinting
 *
 * Computes conservation ratios, spectral gap, Cheeger constant,
 * spectral fingerprints, and anomaly detection — using Chapel's
 * parallel constructs throughout.
 */

module Conservation {
  use TensionGraph;
  use Laplacian;
  use Eigen;
  use Math;

  /* Conservation ratio for one eigenvector */
  record Ratio {
    var eigenvectorIndex: int;
    var eigenvalue: real;
    var ratio: real;
  }

  /* Full conservation report */
  record ConservationReport {
    var n: int;
    var ratios: [0..#n] real;
    var spectralGap: real;
    var cheegerConstant: real;
    var anomalyCount: int;

    proc init(n: int) {
      this.n = n;
    }
  }

  /* Anomaly record */
  record Anomaly {
    var vertexId: int;
    var eigenvectorIndex: int;
    var score: real;
    var description: string;
  }

  /* Conservation ratio: measures how well an eigenvector preserves
   * the vertex attribute structure via projection gradient variance.
   */
  proc conservationRatio(eigen: borrowed EigenDecomp, attr: [] real,
                         eigenvectorIndex: int): real {
    var n = eigen.n;

    // Project attribute onto eigenvector
    var projection: [0..#n] real;
    forall i in 0..#n do
      projection[i] = attr[i] * eigen.eigenvectors[i, eigenvectorIndex];

    if n < 2 then return 0.0;

    // Gradient of consecutive projected values
    var gradient: [0..#(n - 1)] real;
    forall i in 0..#(n - 1) do
      gradient[i] = projection[i + 1] - projection[i];

    // Variance of gradient (parallel reduction)
    var mean: real = (+ reduce gradient) / (n - 1):real;
    var diffs: [0..#(n - 1)] real;
    forall i in 0..#(n - 1) do
      diffs[i] = (gradient[i] - mean) * (gradient[i] - mean);
    var variance: real = (+ reduce diffs) / (n - 1):real;

    return variance;
  }

  /* Cheeger constant approximation from Fiedler vector.
   * Partitions vertices by Fiedler vector sign and computes
   * the cut ratio.
   */
  proc cheegerConstant(lap: borrowed LaplacianMatrix, fiedler: [] real): real {
    var n = lap.n;

    // Partition by Fiedler vector sign
    var inS: [0..#n] bool;
    forall i in 0..#n do
      inS[i] = (fiedler[i] < 0.0);

    // Compute cut and volume in parallel
    var cut: real = 0.0;
    var volS: real = 0.0;

    forall i in 0..#n with (+ reduce cut, + reduce volS) {
      if inS[i] {
        for j in 0..#n {
          if i != j {
            var w = -lap.matrix[i, j]; // off-diagonal = -weight
            volS += w;
            if !inS[j] then cut += w;
          }
        }
      }
    }

    // Total volume from diagonal (degree)
    var totalVol: real = + reduce lap.matrix[0..#n, 0..#n] / 2.0;
    var volComp: real = totalVol - volS;
    var minVol: real = min(volS, volComp);

    if minVol < 1e-15 then return 0.0;
    return cut / minVol;
  }

  /* Spectral fingerprint: hash eigenvalues into hex string.
   * Uses a murmur3-style mixing of double bit patterns.
   */
  proc computeFingerprint(eigenvalues: [] real, n: int): string {
    use IO;
    var hex = "";
    for i in 0..#n {
      // Mix bits of eigenvalue
      var bits: uint(64) = 0;
      // Simple hash: multiply by large prime and XOR-fold
      var val = eigenvalues[i];
      var u: uint(64);
      // Bit-cast via memory — Chapel's approach
      memmove(c_ptrTo(u), c_ptrTo(val), 8);
      u ^= (u >> 33):uint(64);
      u *= 0xff51afd7ed558ccd:uint(64);
      u ^= (u >> 33):uint(64);
      u *= 0xc4ceb9fe1a85ec53:uint(64);
      u ^= (u >> 33):uint(64);
      hex += "%016xu".format(u);
    }
    return hex;
  }

  /* Compare two fingerprints. Returns similarity in [0,1]. */
  proc compareFingerprints(fp1: string, fp2: string): real {
    var len1 = fp1.length;
    var len2 = fp2.length;
    if len1 == 0 && len2 == 0 then return 1.0;

    var minLen = min(len1, len2);
    var maxLen = max(len1, len2);
    if maxLen == 0 then return 1.0;

    var matches: int = 0;
    for i in 0..#minLen do
      if fp1[i] == fp2[i] then matches += 1;

    return matches:real / maxLen:real;
  }

  /* Detect anomalies: vertices whose attribute deviates significantly
   * from the eigenvector-projected expectation.
   */
  proc detectAnomalies(eigen: borrowed EigenDecomp, attr: [] real,
                       threshold: real = 2.0): [] Anomaly {
    var n = eigen.n;
    var scores: [0..#n] real;

    // Use the Fiedler vector (index 1) for anomaly scoring
    if n < 2 then return [new Anomaly(0, 0, 0.0, "")];

    forall i in 0..#n do {
      var projection: real = 0.0;
      for k in 0..#min(5, n) do
        projection += eigen.eigenvectors[i, k] * attr[i];
      scores[i] = abs(attr[i] - projection);
    }

    var meanScore: real = (+ reduce scores) / n:real;
    var varScore: real = 0.0;
    for i in 0..#n do varScore += (scores[i] - meanScore) ** 2;
    varScore /= n:real;
    var stdScore: real = sqrt(varScore);

    // Collect anomalies
    var count: int = 0;
    for i in 0..#n do
      if stdScore > 1e-15 && scores[i] / stdScore > threshold then
        count += 1;

    var anomalies: [0..#count] Anomaly;
    var idx: int = 0;
    for i in 0..#n {
      if stdScore > 1e-15 && scores[i] / stdScore > threshold {
        anomalies[idx] = new Anomaly(i, 1, scores[i],
          "Vertex " + i:string + " deviates from spectral expectation");
        idx += 1;
      }
    }

    return anomalies;
  }

  /* Build a full conservation report from graph, eigen, and laplacian */
  proc conservationReport(graph: borrowed TensionGraph,
                          eigen: borrowed EigenDecomp,
                          lap: borrowed LaplacianMatrix): ConservationReport {
    var n = graph.n;
    var report = new ConservationReport(n);

    // Parallel conservation ratios
    forall k in 0..#n with (ref report) {
      report.ratios[k] = conservationRatio(eigen, graph.vertexAttrs, k);
    }

    report.spectralGap = spectralGap(eigen);

    // Cheeger constant from Fiedler vector
    if n >= 2 {
      var fiedler: [0..#n] real;
      forall i in 0..#n do fiedler[i] = eigen.eigenvectors[i, 1];
      report.cheegerConstant = cheegerConstant(lap, fiedler);
    }

    // Anomaly count
    var anomalies = detectAnomalies(eigen, graph.vertexAttrs);
    report.anomalyCount = anomalies.size;

    return report;
  }
}
