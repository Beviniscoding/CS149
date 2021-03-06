#include "bfs.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cstddef>
#include <omp.h>

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
    #pragma omp parallel for schedule(dynamic, 1000)

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



            if (__sync_val_compare_and_swap(&distances[outgoing], distances[outgoing], distances[node] + 1) == 1) {
                int index = new_frontier_array[threadIndex]->count++;
                new_frontier_array[i]->vertices[threadIndex] = outgoing;
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

    // TODO: How to init this with new frontiers?
    // omp get thread num for size of vertex set. numthreads ver
    vertex_set** new_frontier_array = new vertex_set*[omp_get_max_threads()];
    //maybe have to max threads* max nodes so that we dont have to do 32 mallocs
    // KEV TODO: reduction for all the new frontiers? we are essentially letting each thread have it's own new frontier then putting it all together (consolidate) into a single new frontier) after all the parallel work is done

    for (int i = 0; i < omp_get_max_threads(); i++) {
        new_frontier_array[i] = new vertex_set();
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

        // vertex_set_clear(new_frontier);

        top_down_step(graph, frontier, new_frontier_array, sol->distances);

#ifdef VERBOSE
    double end_time = CycleTimer::currentSeconds();
    printf("frontier=%-10d %.4f sec\n", frontier->count, end_time - start_time);
#endif
    // TODO: how to put all the new frontiers into one frontier
    // for loop of membpy into the new frontier for the count into vertex list.
    // KEV TODO: somehow use reduction to put the new frontiers together into one actual new frontier.
        for (int i = 0; i < omp_get_max_threads(); i++) {
            memcpy(new_frontier + new_frontier->count, new_frontier_array[i], new_frontier_array[i]->count);
            new_frontier->count += new_frontier_array[i]->count;
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
    vertex_set* frontier,
    vertex_set* new_frontier,
    int* distances)
{
#pragma omp parallel for schedule(dynamic, 1000)
    // KEV TODO: somehow do reduction stuffs on the new frontier???
    for (int i=0; i<g->num_nodes; i++) {
    // if already in, then break TODO: how to break in parallel for loop? is it just break or continue?
        if (distances[i] == NOT_VISITED_MARKER) {
            
            // Loop through everything that connects to this node
            for (int i=0; i<num_nodes(g); i++) {
                // Vertex is typedef'ed to an int. Vertex* points into g.outgoing_edges[]
                const Vertex* start = incoming_begin(g, i);
                const Vertex* end = incoming_end(g, i);
                for (const Vertex* v=start; v!=end; v++) {
                    // If in frontier
                    if (frontier->vertices[*v] == 1) {
                        // This node is in new frontier
                        new_frontier->vertices[new_frontier->count] = i;
                        // Set distance as v's distance + 1 TODO: need CAS here?
                        distances[i] = distances[*v] + 1;
                        // Increment the new frontier counter atomically
                        
                        #pragma omp atomic // capture
                        new_frontier ->count = new_frontier -> count + 1;



                    }
                }
            }
        }
    
                






      
    }






}

void bfs_bottom_up(Graph graph, solution* sol)
{
    // Everything is init as not visited until populated with a distance
    // sol->distances[i] = NOT_VISITED_MARKER is populated with distance


    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes); 
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;
    // initialize all nodes to NOT_VISITED
    for (int i=0; i<graph->num_nodes; i++) {
        sol->distances[i] = NOT_VISITED_MARKER;
    }

    // setup frontier flags with first as 1
    frontier->vertices[ROOT_NODE_ID] = 1; // need to flush out frontier every time and reset
    frontier->count++; // TODO: it's init together in top-down -- how to understand it?
    sol->distances[ROOT_NODE_ID] = 0;

    while(frontier->count != 0) {

        vertex_set_clear(new_frontier); // TODO: does this set things to 0 again?
        bottom_up_step(graph, frontier, new_frontier, sol->distances);
        // TODO: ask if it's better to for loop like this to make frontier hold index's
        for (int i = 0; i < new_frontier->count) {
            // KEV TODO: Make frontier into new frontier. new frontier currently holds index's of which nodes should be in frontier. frontier is a flags list: frontier[index] = 1 means node index is in frontier
            // We do this so that the switch thingy will be easier. We could also implement this such that frontier holds index's rather than flags (this would be best for switch), but Idk if this is better for being fast
        }

         // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;



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

    int num_unvisited_nodes = graph->num_nodes;

    vertex_set list1;
    vertex_set list2;
    vertex_set_init(&list1, graph->num_nodes);
    vertex_set_init(&list2, graph->num_nodes);

    vertex_set* frontier = &list1;
    vertex_set* new_frontier = &list2;

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

        if (frontier->count > num_unvisited_nodes) {
            for (int i = 0; i < frontier->count; i++) {

            }

            bottom_up_step(graph, frontier, new_frontier, sol->distances);

            //num nodes visited vs num nodes in frontier
        } else {
            // TODO: Set up top down reduction stuffs

            //top_down_step(graph, frontier, &new_frontier, sol->distances);
            // TODO: Put together new frontier from array of new frontiers
        }







        // swap pointers
        vertex_set* tmp = frontier;
        frontier = new_frontier;
        new_frontier = tmp;
    }




    // CS149 students:
    //
    // You will need to implement the "hybrid" BFS here as
    // described in the handout.
}
