#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>
#include <iostream>

#include "../common/CycleTimer.h"
#include "../common/graph.h"

#define ROOT_NODE_ID 0
#define NOT_VISITED_MARKER -1

void vertex_set_clear(vertex_set* list) {
    list->count = 0;
}

void vertex_set_init(vertex_set* list, int count) {
    list->max_vertices = count;
    list->vertices = (int*)malloc(sizeof(int) * list->max_vertices);
    vertex_set_clear(list);
}

// Take one step of "top-down" BFS.  For each vertex on the frontier,
// follow all outgoing edges, and add all neighboring vertices to the
// new_frontier.
void top_down_step(
    Graph g,
    vertex_set* frontier,
    vertex_set** new_frontier_array,
    int* distances)
{
    #pragma omp parallel for// schedule(dynamic, 1000)

    for (int i=0; i<frontier->count; i++) {
        int threadIndex =  omp_get_thread_num();

        int node = frontier->vertices[i];

        int start_edge = g->outgoing_starts[node];
        int end_edge = (node == g->num_nodes - 1)
                           ? g->num_edges
                           : g->outgoing_starts[node + 1];

        // attempt to add all neighbors to the new frontier
        for (int neighbor=start_edge; neighbor<end_edge; neighbor++) {
            int outgoing = g->outgoing_edges[neighbor];


            if(distances[outgoing] != NOT_VISITED_MARKER){
              continue;
            }

            if (__sync_bool_compare_and_swap(&distances[outgoing],NOT_VISITED_MARKER, distances[node] + 1)){//__sync_val_compare_and_swap(&distances[outgoing], distances[outgoing], distances[node] + 1) == NOT_VISITED_MARKER) {
                int index = new_frontier_array[threadIndex] -> count++;
        //        std::cout << "distindex and sol: " << outgoing << ", " << distances[outgoing] << std::endl;
                new_frontier_array[threadIndex] -> vertices[index] = outgoing;
            }
        }
    }
}

// Implements top-down BFS.
//
// Result of execution is that, for each node in the graph, the
// distance to the root is stored in sol.distances.
void bfs_top_down(Graph graph, solution* sol) {

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // TODO: reduction for all new frontiers

   // vertex_set** new_frontier_array = (vertex_set**)malloc(sizeof(vertex_set*) * omp_get_max_threads());
    //std::cout << "omp threads: " << omp_get_max_threads() << std::endl;
    vertex_set** new_frontier_array = new vertex_set*[omp_get_max_threads()];

    #pragma omp parallel for
    for (int i = 0; i < omp_get_max_threads(); i++)  {
        new_frontier_array[i] = new vertex_set();
      //  std::cout << "thread index: " << i << std::endl;
        vertex_set_init(new_frontier_array[i], graph->num_nodes);
    }

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {

#ifdef VERBOSE
        double start_time = CycleTimer::currentSeconds();
#endif

        vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier_array, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
    // TODO: how to put all the new frontiers into one frontier
        for (int i = 0; i < omp_get_max_threads(); i++) {
            memcpy(new_frontier -> vertices + new_frontier -> count , new_frontier_array[i] -> vertices, new_frontier_array[i] -> count * sizeof(int));
            new_frontier->count += new_frontier_array[i] -> count;
            vertex_set_clear(new_frontier_array[i]);
        }

        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }
}

void bottom_up_step(
    Graph g,
    int* distances,
    bool* indic_frontier,
    bool* new_frontiers,
    int& count)
{
int numNodes = g -> num_nodes;

//I am having issues with different amounts of threads also ask about the atomic stuff

#pragma omp parallel for schedule(dynamic, 1000) reduction(+:count)
    // KEV TODO: somehow do reduction stuffs on the new frontier???
    for (int node = 0; node < g -> num_nodes; node++){
    // if already in, then break TODO: how to break in parallel for loop? is it just break or continue?
        if (distances[node] == NOT_VISITED_MARKER) {
            // Loop through everything that connects to this node
            // Vertex is typedef'ed to an int. Vertex* points into g.outgoing_edges[]
            int start_edge = g->incoming_starts[node];
            int end_edge = (node == g->num_nodes - 1)
                            ? g->num_edges
                            : g->incoming_starts[node + 1];
            for (int neighbor = start_edge; neighbor < end_edge; neighbor++) {
                 // If in frontier
                 int incoming = g -> incoming_edges[neighbor];
                 if (indic_frontier[incoming]) {
                    // This node is in new frontier
                     new_frontiers[node] = true;
                     //#pragma omp atomic write
                     distances[node] = distances[incoming] + 1;
                     //#pragma omp atomic
                     count ++;
                     break;
                  }
                }
            }
    }
}

void bfs_bottom_up(Graph graph, solution* sol)
{
    // Everything is init as not visited until populated with a distance
    // sol->distances[i] = NOT_VISITED_MARKER is populated with distance

    // initialize all nodes to NOT_VISITED
    //bfs_top_down(graph,sol);
    //retrun;
    int numNodes = graph -> num_nodes;
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }
    bool* indic_frontier = new bool[numNodes];//might have to use new with other arrays which may be iffy 
    memset(indic_frontier, 0, sizeof(bool) * numNodes);
    // setup frontier flags with first as 1
    indic_frontier[ROOT_NODE_ID] = 1;
    sol->distances[ROOT_NODE_ID] = 0;
    int converged = 1;
    bool* new_frontiers = new bool[numNodes];
    memset(new_frontiers, 0, sizeof(bool) * numNodes);
    while(converged != 0) {
        converged = 0;

        memset(new_frontiers, 0, sizeof(bool) * numNodes);
        
        bottom_up_step(graph, sol->distances, indic_frontier, new_frontiers,converged);
        
        //delete(new_frontiers);
        // swap pointers
        bool* tmp = indic_frontier;
        indic_frontier = new_frontiers;
        new_frontiers = tmp;
        //clear new_frontier

    }


    // CS149 students:
    //
    // You will need to implement the "bottom up" BFS here as
    // described in the handout.
    //
    // As a result of your code's execution, sol.distances should be
    // correctly populated for all nodes in the graph.
    //
    // As was done in the top-down case, you may wish to organize your
    // code by creating subroutine bottom_up_step() that is called in
    // each step of the BFS process.
}

void bfs_hybrid(Graph graph, solution* sol)
{
    //bfs_top_down(graph,sol);

    //set up for bottom up

    int numNodes = graph -> num_nodes;
   
    bool* indic_frontier = new bool[numNodes];//might have to use new with other arrays which may be iffy 
    memset(indic_frontier, 0, sizeof(bool) * numNodes);
    // setup frontier flags with first as 1
    indic_frontier[ROOT_NODE_ID] = 1;
    
    int converged = 1;
    bool* new_frontiers = new bool[numNodes];
    memset(new_frontiers, 0, sizeof(bool) * numNodes);

    bool lastWasBottom = false;


    // set up for top down
    int num_unvisited_nodes = numNodes;
    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

    // TODO: reduction for all new frontiers

   // vertex_set** new_frontier_array = (vertex_set**)malloc(sizeof(vertex_set*) * omp_get_max_threads());
    //std::cout << "omp threads: " << omp_get_max_threads() << std::endl;
    vertex_set** new_frontier_array = new vertex_set*[omp_get_max_threads()];
    #pragma omp parallel for
    for (int i = 0; i < omp_get_max_threads(); i++)  {
        new_frontier_array[i] = new vertex_set();
      //  std::cout << "thread index: " << i << std::endl;
        vertex_set_init(new_frontier_array[i], graph->num_nodes);
    }

    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier with the root node
    frontier->vertices[frontier->count++] = ROOT_NODE_ID;
    sol->distances[ROOT_NODE_ID] = 0;

    while (frontier->count != 0) {
        num_unvisited_nodes -= frontier->count;
        vertex_set_clear(new_frontier);
        
        if (frontier->count / (num_unvisited_nodes+1) > 3.25) {
            //std::cout << "going into bottom up" << std::endl;
            if (!lastWasBottom) {
                // last time was top down, must readjust stuffs
                for (int i = 0; i < frontier -> count; i++){
                    indic_frontier[frontier->vertices[i]] = true;
                    
                }
            }
            // do bottom up step
            lastWasBottom = true;
            frontier -> count = 0;
            memset(new_frontiers, 0, sizeof(bool) * numNodes);
            bottom_up_step(graph, sol->distances, indic_frontier, new_frontiers,frontier -> count);
            // swap pointers
            bool* tmp = indic_frontier;
            indic_frontier = new_frontiers;
            new_frontiers = tmp;

        } else {
           // std::cout << "going into top down" << std::endl;
            if (lastWasBottom) {
                // last time was bottom up, need to adjust for top down
                int num_input = 0;
                   for (int i = 0; i < numNodes; i++){
                    if (indic_frontier[i]){
                        frontier->vertices[num_input] = i;
                        num_input++;
                    }
                }
            }
           // std::cout << "finish set up top down step" << std::endl;
            lastWasBottom = false;
        vertex_set_clear(new_frontier);
        top_down_step(graph, frontier, new_frontier_array, sol->distances);

    // TODO: how to put all the new frontiers into one frontier
        for (int i = 0; i < omp_get_max_threads(); i++) {
            memcpy(new_frontier -> vertices + new_frontier -> count , new_frontier_array[i] -> vertices, new_frontier_array[i] -> count * sizeof(int));
            new_frontier->count += new_frontier_array[i] -> count;
            vertex_set_clear(new_frontier_array[i]);
        }
        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
        }
    }
    
    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}
