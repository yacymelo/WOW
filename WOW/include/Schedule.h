#include "GMDP/gmdp.h"
//#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <vector>
#include <utility>
#include <sys/time.h>

#ifdef __ASSERT
#include <assert.h>
#endif

#include "Graph.h"
#include "GraphProgram.h"
#include "SPMV.h"

namespace PIMSecedule {

const int UNTIL_CONVERGENCE = -1;

template<class T, class U, class V>
struct run_graph_program_temp_structure {
PIMSecedule::SpVec<PIMSecedule::DenseSegment<T> >* px;
PIMSecedule::SpVec<PIMSecedule::DenseSegment<U> >* py;
};

template<class T, class U, class V, class E>
struct run_graph_program_temp_structure<T,U,V> graph_program_init(const GraphProgram<T,U,V,E>& gp, const Graph<V, E>& g) {

struct run_graph_program_temp_structure<T,U,V> rgpts;
rgpts.px = new PIMSecedule::SpVec<PIMSecedule::DenseSegment<T> >(g.nvertices, PIMSecedule::get_global_nrank(), PIMSecedule::vector_partition_fn);
T _t;
rgpts.px->setAll(_t);
rgpts.py = new PIMSecedule::SpVec<PIMSecedule::DenseSegment<U> >(g.nvertices, PIMSecedule::get_global_nrank(), PIMSecedule::vector_partition_fn);
U _u;
rgpts.py->setAll(_u);
return rgpts;
}

template<class T, class U, class V>
void graph_program_clear(struct run_graph_program_temp_structure<T,U,V>& rgpts) {
delete rgpts.px;
delete rgpts.py;
}


template <class T,class U, class V, class E>
void send_message(const bool& a, const V& _v, T* b, void* gpv) {
GraphProgram<T,U,V,E>* gp = (GraphProgram<T,U,V,E>*) gpv;
if(a == true) {
gp->send_message(_v, *b);
}
}

template <class T, class U, class V, class E>
void apply_func(const U& y, V* b, void* gpv) {
GraphProgram<T,U,V,E>* gp = (GraphProgram<T,U,V,E>*) gpv;
gp->apply(y, *b);
}

template <class T, typename U, class V, class E>
void run_graph_program(GraphProgram<T,U,V,E>* gp, Graph<V,E>& g, int iterations=1, struct run_graph_program_temp_structure<T,U,V>* rgpts=NULL) { //iterations = -1 ==> until convergence
int it = 0;
int converged = 1;

struct timeval start, end, init_start, init_end, iteration_start, iteration_end;
double time;
int global_myrank = PIMSecedule::get_global_myrank();

gettimeofday(&init_start, 0);


auto act = gp->getActivity();

PIMSecedule::SpVec<PIMSecedule::DenseSegment<T> >* px;
PIMSecedule::SpVec<PIMSecedule::DenseSegment<U> >* py;

if (rgpts == NULL) {
px = new PIMSecedule::SpVec<PIMSecedule::DenseSegment<T> >(g.nvertices, PIMSecedule::get_global_nrank(), PIMSecedule::vector_partition_fn);
T _t;
px->setAll(_t);
py = new PIMSecedule::SpVec<PIMSecedule::DenseSegment<U> >(g.nvertices, PIMSecedule::get_global_nrank(), PIMSecedule::vector_partition_fn);
U _u;
py->setAll(_u);
}

PIMSecedule::SpVec<PIMSecedule::DenseSegment<T> >& x = (rgpts==NULL)?(*px):*(rgpts->px);//*px;
PIMSecedule::SpVec<PIMSecedule::DenseSegment<U> >& y = (rgpts==NULL)?(*py):*(rgpts->py);//*py;

if (act == ALL_VERTICES) {
g.setAllActive();
}
#ifdef __TIMING
printf("Nvertices = %d \n", g.getNumberOfVertices());
#endif

gettimeofday(&init_end, 0);

#ifdef __TIMING
time = (init_end.tv_sec-init_start.tv_sec)*1e3+(init_end.tv_usec-init_start.tv_usec)*1e-3;
printf("PIMSecedule init time = %f ms \n", time);
#endif

while(1) {
gettimeofday(&iteration_start, 0);

PIMSecedule::Clear(&x);
PIMSecedule::Clear(&y);
converged = 1;

gettimeofday(&start, 0);

PIMSecedule::IntersectReduce(g.active, g.vertexproperty, &x, send_message<T,U,V,E>, (void*)gp);

#ifdef __TIMING
printf("x.length = %d \n", x.getNNZ());
#endif
gettimeofday(&end, 0);
#ifdef __TIMING
time = (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3;
printf("Send message time = %.3f ms \n", time);
#endif

gettimeofday(&start, 0);


//do SpMV
if (gp->getOrder() == OUT_EDGES) {

SpMTSpV(g, gp, &x, &y);

} else if (gp->getOrder() == IN_EDGES) {

SpMSpV(g, gp, &x, &y);

} else if (gp->getOrder() == ALL_EDGES) {

SpMTSpV(g, gp, &x, &y);
SpMSpV(g, gp, &x, &y);

} else {
printf("Unrecognized option \n");
exit(1);
}
gettimeofday(&end, 0);
#ifdef __TIMING
time = (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3;
printf("SPMV time = %.3f ms \n", time);
#endif

gettimeofday(&start, 0);
g.setAllInactive();

//update state and activity and check for convergence if needed
int nout = 0;
int total_search = 0;
int local_converged = 1;
converged = 1;

//PIMSecedule::IntersectReduce(g.active, y, &g.vertexproperty, set_y<U,V>);
//auto apply_func  = set_y_apply<U,V>;
//PIMSecedule::Apply(y, &g.vertexproperty, apply_func<T,U,V>, (void*)gp);
for(int segmentId = 0 ; segmentId < y.nsegments ; segmentId++)
{
if(y.nodeIds[segmentId] == global_myrank)
{
auto segment = y.segments[segmentId]->properties;
auto vpValueArray = g.vertexproperty->segments[segmentId]->properties->value;
#pragma omp parallel for reduction(&:local_converged)
for (int i = 0; i < y.segments[segmentId]->num_ints; i++) {
unsigned int value = segment->bit_vector[i];
while (value != 0) {
int last_bit = _bit_scan_forward(value);
int idx = i*32 + last_bit;

V old_prop;
//old_prop = g.vertexproperty.segments[segmentId].properties->value[idx];
old_prop = vpValueArray[idx];

//gp->apply(segment->value[idx], g.vertexproperty.segments[segmentId].properties->value[idx]);
gp->apply(segment->value[idx], vpValueArray[idx]);
if (old_prop != vpValueArray[idx]) {
g.active->segments[segmentId]->properties->value[idx] = true;
PIMSecedule::set_bitvector(idx, g.active->segments[segmentId]->properties->bit_vector);
local_converged = 0;
}

value &= (~(1<<last_bit));
}
}

}
}
MPI_Allreduce(&local_converged, &converged, 1, MPI_INT, MPI_LAND, MPI_COMM_WORLD);

gettimeofday(&end, 0);
#ifdef __TIMING
time = (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3;
printf("Apply time = %.3f ms \n", time);
#endif

gettimeofday(&start, 0);

gp->do_every_iteration(it);

#ifdef __TIMING
gettimeofday(&end, 0);
time = (end.tv_sec-start.tv_sec)*1e3+(end.tv_usec-start.tv_usec)*1e-3;
printf("Do every iteration time = %.3f ms \n", time);
#endif

gettimeofday(&iteration_end, 0);
#ifdef __TIMING
time = (iteration_end.tv_sec-iteration_start.tv_sec)*1e3+(iteration_end.tv_usec-iteration_start.tv_usec)*1e-3;
printf("Iteration %d :: %f msec :: updated %d vertices :: changed %d vertices \n", it, time, y.getNNZ(), g.active->getNNZ());
#endif

if (act == ALL_VERTICES) {
g.setAllActive();
}

it++;
if (it == iterations) {
break;
}
if (iterations <= 0 && converged == 1) {
break;
}
}

struct timeval clear_start, clear_end;
gettimeofday(&clear_start, 0);

if (rgpts == NULL) {
delete px;
delete py;
}

gettimeofday(&clear_end, 0);
#ifdef __TIMING
time = (clear_end.tv_sec-clear_start.tv_sec)*1e3+(clear_end.tv_usec-clear_start.tv_usec)*1e-3;
printf("PIMSecedule clear time = %f msec \n", time);
#endif

printf("Completed %d iterations \n", it);

}

} //namespace PIMSecedule
