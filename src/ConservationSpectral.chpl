/* ConservationSpectral.chpl — Main module for the Conservation Spectral SDK (Chapel)
 *
 * This module re-exports the public API. Users import this single module.
 *
 * Version: 0.1.0
 */

module ConservationSpectral {
  use TensionGraph;
  use Laplacian;
  use Eigen;
  use Conservation;
  use Tracker;

  public use TensionGraph;
  public use Laplacian;
  public use Eigen;
  public use Conservation;
  public use Tracker;

  /* Version constants */
  param CS_VERSION_MAJOR = 0;
  param CS_VERSION_MINOR = 1;
  param CS_VERSION_PATCH = 0;

  proc csVersion(): string {
    return CS_VERSION_MAJOR:string + "." + CS_VERSION_MINOR:string + "." + CS_VERSION_PATCH:string;
  }

  /* Convenience: full analysis pipeline in one call */
  proc csAnalyze(graph: borrowed TensionGraph, k: int = 0): ConservationReport {
    var n = graph.n;
    var kk = if k <= 0 then n else k;

    var lap = buildLaplacian(graph, normalized=false);
    var eigen = eigendecompose(lap, kk);
    var report = conservationReport(graph, eigen, lap);
    return report;
  }
}
