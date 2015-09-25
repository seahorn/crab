/*******************************************************************************
 * Iterators and simple types
 ******************************************************************************/

#ifndef DBM_ITER_H__
#define DBM_ITER_H__

#include <crab/domains/dbm/expr.h>

// SEMANTICS:
// An edge (i, w, j) represents the constraint [[i - j <= w]]
typedef struct {
  dbm d;
  int si;
  int di;
  adjlist* slist;
} edge_iter;

edge_iter edge_iterator(dbm d);
bool edges_end(edge_iter& iter);
int src(edge_iter& iter);
int dest(edge_iter& iter);
void next_edge(edge_iter& iter);
val_t iter_val(edge_iter& iter);
val_t edge_val(dbm x, int i, int j);
void next_src(edge_iter&iter);
void next_dest(edge_iter&iter);
bool srcs_end(edge_iter&iter);
bool dests_end(edge_iter&iter);

#endif
