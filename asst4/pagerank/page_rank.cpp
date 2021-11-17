#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"


// pageRank --
//
// g:           graph to process (see common/graph.h)
// solution:    array of per-vertex vertex scores (length of array is num_nodes(g))
// damping:     page-rank algorithm's damping parameter
// convergence: page-rank algorithm's convergence threshold
//
void pageRank(Graph g, double* solution, double damping, double convergence)
{


  // initialize vertex weights to uniform probability. Double
  // precision scores are used to avoid underflow for large graphs

  int numNodes = num_nodes(g);
  double equal_prob = 1.0 / numNodes;
  for (int i = 0; i < numNodes; ++i) {
    solution[i] = equal_prob;
  }
  bool converged = false;
  while (!converged) {
    int global_diff = 0;
    for (int node = 0; node < numNodes; node++){
       int new_score = 0;
       int start_edge = g->outgoing_starts[node];
       int end_edge = (node == g->num_nodes - 1)
                 ? g->num_edges
                 : g->outgoing_starts[node + 1];
       //this is looping through the outgoing edges of node vi
       //sum over all vj  nodes reachable from incoming edges
       //so number of nodes coming from direct incoming nodes of vi?
       int count_edges_leaving_vj = 0;
       int total_score_old = 0;
       for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
         int neighbor_start = g -> outgoing_starts[neighbor];
         int end_edge = (neighbor == g->num_nodes - 1)
                  ? g->num_edges
                  : g->outgoing_starts[neighbor + 1];
         int addition = neighbor_end - neighbor_start;
         count_edges_leaving_vj +=addition;
         total_score_old_vj += solution[neighbor];
       }
       new_score = total_score_old_vj/count_edges_leaving_vj;

       //do other additions
    }

}
  /*
     CS149 students: Implement the page rank algorithm here.  You
     are expected to parallelize the algorithm using openMP.  Your
     solution may need to allocate (and free) temporary arrays.

     Basic page rank pseudocode is provided below to get you started:

     // initialization: see example code above
     score_old[vi] = 1/numNodes;

     while (!converged) {

       // compute score_new[vi] for all nodes vi:
       score_new[vi] = sum over all nodes vj reachable from incoming edges
                          { score_old[vj] / number of edges leaving vj  }
       score_new[vi] = (damping * score_new[vi]) + (1.0-damping) / numNodes;

       score_new[vi] += sum over all nodes v in graph with no outgoing edges
                          { damping * score_old[v] / numNodes }

       // compute how much per-node scores have changed
       // quit once algorithm has converged

       global_diff = sum over all nodes vi { abs(score_new[vi] - score_old[vi]) };
       converged = (global_diff < convergence)
     }

   */
}
