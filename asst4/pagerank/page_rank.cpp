#include "page_rank.h"

#include <stdlib.h>
#include <cmath>
#include <omp.h>
#include <utility>

#include "../common/CycleTimer.h"
#include "../common/graph.h"
#include <iostream>


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
  double* new_score = new double[numNodes];
  double equal_prob = 1.0 / num_nodes(g);
  #pragma omp parallel for
  for (int i = 0; i < numNodes; i++) {
    solution[i] = equal_prob;
  }
  bool converged = false;
  double global_diff = 0;
  double no_out = 0.0;
  //std::cout << "convergence: " << convergence << " numNodes: " << numNodes << std::endl;
  while (!converged) {
    global_diff = 0.0;
    no_out = 0.0;
    #pragma omp parallel for reduction (+:no_out)
    for (int i = 0 ; i < numNodes; i++){
        if (outgoing_size(g,i) == 0 ) {
            no_out += (damping * solution[i] ) / numNodes;
          }
    }
    #pragma omp parallel for reduction (+:global_diff)
    for (int node = 0; node < g->num_nodes; node++){
       new_score[node] = 0.0;
       int start_edge = g->incoming_starts[node];
       int end_edge = (node == g->num_nodes - 1)
                 ? g->num_edges
                : g->incoming_starts[node + 1];
       for (int neighbor = start_edge; neighbor < end_edge; neighbor++) {
         new_score[node] += (solution[g -> incoming_edges[neighbor]] / outgoing_size(g,g-> incoming_edges[neighbor]));
       }
       new_score[node] = (damping * new_score[node]) + ((1.0-damping) / num_nodes(g));
       new_score[node] += no_out;
       global_diff += fabs(new_score[node] - solution[node]);
     //  solution[node] = new_score[node];//might be issue if some instances get here before referencing before
    }
    #pragma omp parallel for
    for (int i = 0; i < numNodes; i++){
      solution[i] = new_score[i];
    }
    converged = (global_diff < convergence);
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
