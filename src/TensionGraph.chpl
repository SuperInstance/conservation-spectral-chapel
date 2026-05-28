/* TensionGraph.chpl — Graph data structure using Chapel domains and arrays
 *
 * Uses a dense adjacency matrix backed by Chapel's domain-based arrays.
 * Supports edge-list construction and parallel adjacency queries.
 */

module TensionGraph {

  /* A single edge in the graph */
  record Edge {
    var from: int;
    var to: int;
    var weight: real;
  }

  /* A graph with n vertices and dense adjacency storage.
   *
   * The adjacency matrix adjMatrix[i,j] holds the edge weight (0 if no edge).
   * vertexAttrs[i] holds the scalar attribute for vertex i.
   */
  record TensionGraph {
    var n: int;
    var directed: bool;
    var vertexAttrs: [0..#n] real;
    var adjMatrix: [0..#n, 0..#n] real;
    var _nEdges: int;

    proc init(n: int, directed: bool = false) {
      this.n = n;
      this.directed = directed;
      // Chapel zero-initializes arrays
    }

    /* Add a vertex attribute */
    proc addVertex(id: int, attr: real) {
      if id < 0 || id >= n then
        halt("TensionGraph.addVertex: index out of bounds");
      vertexAttrs[id] = attr;
    }

    /* Add an edge with weight */
    proc addEdge(from: int, to: int, weight: real) {
      if from < 0 || from >= n || to < 0 || to >= n then
        halt("TensionGraph.addEdge: index out of bounds");
      adjMatrix[from, to] += weight;
      if !directed then
        adjMatrix[to, from] += weight;
      _nEdges += 1;
    }

    /* Number of edges (undirected counted once per addEdge call) */
    proc numEdges(): int {
      return _nEdges;
    }

    /* Degree of vertex i */
    proc degree(i: int): real {
      return + reduce adjMatrix[i, ..];
    }

    /* Parallel degree computation for all vertices */
    proc allDegrees(): [] real {
      var deg: [0..#n] real;
      forall i in 0..#n do
        deg[i] = + reduce adjMatrix[i, ..];
      return deg;
    }
  }

  /* Factory proc */
  proc createGraph(n: int, directed: bool = false): TensionGraph {
    return new TensionGraph(n, directed);
  }
}
