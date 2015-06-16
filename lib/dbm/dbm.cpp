#include <climits>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <algorithm>
#include <vector>

#include <ikos/domains/dbm/dbm.h>
#include <ikos/domains/dbm/expr.h>
#include <ikos/domains/dbm/dbm_iter.h>
#include <ikos/domains/dbm/util/Heap.h>

// #define NDEBUG // disable assert's

#define VAR_MAX ((short) ((1<<(8*sizeof(short)-1))-1))

#define VAL_MAX INT_MAX
// #define VAL_MIN INT_MIN
#define VAL_MIN (-INT_MAX) // We need -VAL_MAX = VAL_MIN.

// These are datastructures used by consistency checking and
// transitive closure. They are global for efficiency.
static int scratch_sz = 0;
static int* _gamma = NULL;
static int* pi_prime = NULL;
static int* fwd_dist = NULL;

// Always initialized to 0.
static char* var_flags = NULL;

// Comparator for use with min-heaps.
typedef struct DistComp {
public:
  DistComp(int* _A)
    : A(_A)
  { }
  bool operator() (int x, int y) const {
    return A[x] < A[y]; 
  }
  int* A;
} DistComp;

static void verify_potentials(dbm abs);

//bool in_graph(dbm x, int i, int j);

adjlist* src_list(dbm abs, int i)
{
  return (adjlist*) (&(abs->srcs[(2 + abs->sz)*i]));
}
adjlist* dest_list(dbm abs, int j)
{
  return (adjlist*) (&(abs->dests[(2 + abs->sz)*j]));
}

static int var_lb(dbm abs, int x);
static int var_ub(dbm abs, int x);

void update_scratch(int sz)
{
  if(scratch_sz < sz)
  {
    var_flags = (char *) realloc(var_flags, sizeof(char)*sz);
    fwd_dist = (int *) realloc(fwd_dist, sizeof(int)*sz);
    _gamma = (int *) realloc(_gamma, sizeof(int)*sz);
    pi_prime = (int *) realloc(pi_prime, sizeof(int)*sz);

    for(; scratch_sz < sz; scratch_sz++)
      var_flags[scratch_sz] = 0;
  }
}

// Check if a given vertex is the source of some edge.
bool src_is_live(dbm abs, int i)
{
  unsigned short inv = src_list(abs, i)->inv;
  return inv < abs->num_srcs && abs->live_srcs[inv] == i;
}

bool dest_is_live(dbm abs, int j)
{
  unsigned short inv = dest_list(abs, j)->inv;
  return inv < abs->num_dests && abs->live_dests[inv] == j;
}

///
// Iterator operations
///

edge_iter edge_iterator(dbm d)
{
  edge_iter iter = { d, 0, 0, NULL};
  if(d->num_srcs > 0)
    iter.slist = src_list(d, d->live_srcs[0]);
  return iter;
}

edge_iter src_iterator(dbm d, int i)
{
  adjlist* ilist = src_list(d, i);
  edge_iter iter = { d, ilist->inv, 0, ilist };
  return iter;
}

edge_iter pred_iterator(dbm d, int j)
{
  adjlist* jlist = dest_list(d, j);
  edge_iter iter = { d, jlist->inv, 0, jlist };
  return iter;
}

// Precondition: iterator allocated with pred_iterator.
int pred(edge_iter& iter)
{
  return iter.slist->elt[iter.di];
}

void next_pred(edge_iter& iter)
{
  iter.di++;
}

bool preds_end(edge_iter& iter)
{
  return iter.di >= iter.slist->sz;
}

int src(edge_iter& iter)
{
  return iter.d->live_srcs[iter.si];
}

bool srcs_end(edge_iter& iter)
{
  return iter.si >= iter.d->num_srcs;
}

void next_src(edge_iter& iter)
{
  iter.si++;
  iter.di = 0;
  if(iter.si < iter.d->num_srcs)
    iter.slist = src_list(iter.d, iter.d->live_srcs[iter.si]);
}

int dest(edge_iter& iter)
{
  return iter.slist->elt[iter.di];
}

bool dests_end(edge_iter& iter)
{
  return iter.di >= iter.slist->sz;
}

void next_dest(edge_iter& iter)
{
  iter.di++;
}

bool edges_end(edge_iter& iter)
{
  return srcs_end(iter);
}
void next_edge(edge_iter& iter)
{
  next_dest(iter);
  if(dests_end(iter))
    next_src(iter); 
}

// FIXME: Need to add allocation checking.
dbm dbm_alloc(int sz)
{

  update_scratch(sz);
  dbm_data* d = new dbm_data;
  
  d->sz = sz;
  d->checked = true;
  d->feasible = true;
  d->closed = true;
  
  d->pi = new int[sz];
  for(int vi = 0; vi < sz; vi++)
    d->pi[vi] = 0;
  
  // Live sources & dests
  d->num_srcs = 0;
  d->live_srcs = new short[sz];
  
  d->num_dests = 0;
  d->live_dests = new short[sz];
  
  // The array of adjacency lists.
  d->srcs = new short[(sz+2)*sz];
  d->dests = new short[(sz+2)*sz];
  
  d->mtx = new einfo[sz*sz];
  
  return d;
}

dbm dbm_top(unsigned int sz)
{
  return dbm_alloc(sz);
}

dbm dbm_bottom()
{
  return NULL;
}

void dbm_dealloc(dbm d)
{
  if(!d)
    return;
  
  delete[] d->pi;

  delete[] d->live_srcs;
  delete[] d->live_dests;

  delete[] d->srcs;
  delete[] d->dests;

  delete[] d->mtx;

  delete d;
}

void dbm_dealloc_ptr(dbm* d)
{
  dbm_dealloc(*d);
}

// Reset counts of sources and dests
void clear_dbm(dbm d)
{
  if(!d)
    return;
  d->num_srcs = 0;
  d->num_dests = 0;
}

// Check the potential function: only for sanity checks
void verify_potentials(dbm abs)
{
#if 1
  return;
#else
  edge_iter iter_check = edge_iterator(abs);
  for(; !edges_end(iter_check); next_edge(iter_check))
  {
    int s = src(iter_check);
    int d = dest(iter_check);
    int w = iter_val(iter_check);
    assert(abs->pi[s] + w - abs->pi[d] >= 0); 

    assert(src_is_live(abs, s));
    assert(dest_is_live(abs, d));
    assert(in_graph(abs, s, d));
  }
#endif 
}

// // Either clear the existing dbm, or allocate
// // a fresh one.
// static dbm maybe_alloc(int sz, dbm* out)
// {
//   if(out)
//   {
//     clear_dbm(*out);
//     return *out;
//   } else {
//     return dbm_alloc(sz);
//   }
// }

dbm dbm_copy(dbm abs)
{
  if(!abs)
    return NULL;
  if(abs->checked && !abs->feasible)
    return NULL;

  int sz = abs->sz;
  dbm ret = dbm_alloc(sz);
  
  // Two options: memcpy, or linear copy.
  // Using memcpy for now, even though it's
  // more or less O(n^2)
  memcpy(ret->pi, abs->pi, sizeof(int)*sz);

  ret->sz = sz;
  ret->checked = abs->checked;
  ret->feasible = abs->feasible;
  ret->closed = abs->closed;

  // // Live sources & dests
  // ret->num_srcs = abs->num_srcs;
  // memcpy(ret->live_srcs, abs->live_srcs, sizeof(short)*sz);

  // ret->num_dests = abs->num_dests;
  // memcpy(ret->live_dests, abs->live_dests, sizeof(short)*sz);

  // // The array of adjacency lists.
  // memcpy(ret->srcs, abs->srcs, sizeof(short)*(sz+2)*sz);
  // memcpy(ret->dests, abs->dests, sizeof(short)*(sz+2)*sz);

  // memcpy(ret->mtx, abs->mtx, sizeof(einfo)*sz*sz);

  edge_iter iter = edge_iterator(abs);
  for(; !edges_end(iter); next_edge(iter))
  {
    dbm_add_edge(ret, src(iter), dest(iter), iter_val(iter));
  }
  return ret;
}

// Domain operations
bool check_feasibility(dbm abs)
{
  if(!abs)
    return false;

  if(abs->checked)
    return abs->feasible;
  
  // Use Bellman-Ford to check feasibility, and restore
  // a valid potential function (if true)
  // Currently not incremental.
  
  // Initialize the potentials.
  for(int vi = 0; vi < abs->sz; vi++)
  {
    abs->pi[vi] = 0;
  }

  // relaxation phase
  // for all vertices
  for(int iter_num = 0; iter_num < abs->sz; iter_num++)
  {
    // for all edges
    edge_iter iter = edge_iterator(abs);
    for(; !srcs_end(iter); next_src(iter))
    {
      int i = src(iter);

      for(; !dests_end(iter); next_dest(iter))
      {
        int j = dest(iter);
        
        int ij_val = iter_val(iter);
        abs->pi[j] = std::min(abs->pi[j], abs->pi[i] + ij_val);
      }
    }
  }

  // check for negative cycle
  edge_iter iter = edge_iterator(abs);
  for(; !edges_end(iter); next_edge(iter))
  {
    int i = src(iter);
    int j = dest(iter);
    int ij_val = edge_val(abs, i, j);
    if(abs->pi[j] > abs->pi[i] + ij_val)
    {
      abs->feasible = false;
      return false;
    }
  }

  verify_potentials(abs);
  abs->checked = true;

  return true;
}

// Update feasibility when a new constraint is added.
// Perform an incremental consistency checking.
bool update_feasibility(dbm x, int i, int j, int c)
{
  // Ensure there's enough scratch space. 
  int sz = x->sz;
  
  for(int vi = 0; vi < sz; vi++)
  {
    _gamma[vi] = 0;
    pi_prime[vi] = x->pi[vi];
  }
  _gamma[j] = x->pi[i] + c - x->pi[j];

  if(_gamma[j] >= 0)
    return true;

  DistComp comp(_gamma);
  Heap<DistComp> heap(comp);

  heap.insert(j);

  while(!heap.empty())
  {
    int es = heap.removeMin();

    pi_prime[es] = x->pi[es] + _gamma[es];

    if(!src_is_live(x, es))
      continue;
    edge_iter iter = src_iterator(x, es);
    for(; !dests_end(iter); next_dest(iter))
    {
      int ed = dest(iter);
      if(pi_prime[ed] == x->pi[ed])
      {
        int gnext_ed = pi_prime[es] + iter_val(iter) - pi_prime[ed];
        if(gnext_ed < _gamma[ed])
        {
          _gamma[ed] = gnext_ed;
          if(heap.inHeap(ed))
          {
            heap.decrease(ed);
          } 
          else {
            heap.insert(ed);
          }
        }
      }
    }
  }
  if(_gamma[i] < 0)
    return false;

  return true;
}


// Assumption: (i -> j) does not yet exist.
void dbm_add_edge(dbm x, int i, int j, int val)
{
  if(!x)
    return;

  assert(!src_is_live(x, i) || !in_graph(x, i, j));
  
  assert(i < x->sz);
  assert(j < x->sz);

  adjlist* ilist = src_list(x, i);
  adjlist* jlist = dest_list(x, j);

  // Check if i or j is not yet live.
  if(((unsigned short) ilist->inv) >= x->num_srcs || 
     x->live_srcs[ilist->inv] != i)
  {
    ilist->sz = 0;
    ilist->inv = x->num_srcs;
    x->live_srcs[x->num_srcs] = i;
    x->num_srcs++;
  }

  if(((unsigned short) jlist->inv) >= x->num_dests || 
     x->live_dests[jlist->inv] != j)
  {
    jlist->sz = 0;
    jlist->inv = x->num_dests;
    x->live_dests[x->num_dests] = j;
    x->num_dests++;
  }

  assert(x->num_srcs <= x->sz);
  assert(x->num_dests <= x->sz);

  // Add i, j to the source/dest lists.
  int eidx = i*x->sz + j;

  x->mtx[eidx].val = val;

  x->mtx[eidx].i_inv = ilist->sz;
  ilist->elt[ilist->sz] = j;
  ilist->sz++;

  x->mtx[eidx].j_inv = jlist->sz;
  jlist->elt[jlist->sz] = i;
  jlist->sz++;

  for(int xi = 0; xi < x->num_srcs; xi++)
  {
    adjlist* xlist = src_list(x, x->live_srcs[xi]);
    assert(xlist->inv == xi);
  }
}

// Precondition:
//   i <- live_srcs
bool in_graph(dbm abs, int i, int j)
{
  assert(abs);
  assert(i < abs->sz);
  assert(j < abs->sz);
  int elt_idx = i*(abs->sz) + j;

  unsigned short i_inv = abs->mtx[elt_idx].i_inv;
  
  adjlist* iadj = src_list(abs, i);

  return (i_inv < iadj->sz) && (iadj->elt[i_inv] == j);
}

int edge_val(dbm x, int i, int j)
{
  assert(x);
  return x->mtx[i*x->sz + j].val;
}

int iter_val(edge_iter& iter)
{
  return edge_val(iter.d, src(iter), dest(iter));
}

int rev_iter_val(edge_iter& iter)
{
  return edge_val(iter.d, pred(iter), iter.d->live_dests[iter.si]);
}

void dijkstra(dbm abs, int svar, 
              std::vector< std::pair<int, int> >& out)
{
  // fwd_dist contains the new non-negative weights
  // fwd_dist[j] = pi[i] + w(i,j) - pi[j] >= 0
  // this is called the "reduced cost"

  // We run dijkstra using fwd_dist rather than original weights. We
  // return out that is a vector of pairs (v,c) where v is the node
  // that has been modified and c is the weight of the shortest path
  // between svar and v but *without* using reducing costs (i.e.,
  // fwd_dist[v] - pi[svar] + pi[v])

  assert(abs);
  assert(src_is_live(abs, svar));
  assert(abs->checked && abs->feasible);
  verify_potentials(abs);
   
  for(int vi = 0; vi < abs->sz; vi++)
  {
    fwd_dist[vi] = INT_MAX;
  }
  fwd_dist[svar] = 0;
  DistComp comp(fwd_dist);
  Heap<DistComp> heap(comp);

  edge_iter iter = src_iterator(abs, svar);
  for(; !dests_end(iter); next_dest(iter))
  {
    int ed = dest(iter);
    fwd_dist[ed] = abs->pi[svar] + iter_val(iter) - abs->pi[ed];
    heap.insert(ed);
  }

  while(!heap.empty())
  {
    int es = heap.removeMin();
    int es_cost = fwd_dist[es] + abs->pi[es];
    int es_val = es_cost - abs->pi[svar];
    if(!in_graph(abs, svar, es) || edge_val(abs, svar, es) > es_val)
      out.push_back( std::pair<int, int>(es, es_val) );

    if(!src_is_live(abs, es))
      continue;

    iter = src_iterator(abs, es);
    for(; !dests_end(iter); next_dest(iter))
    {
      int ed = dest(iter);
      int v = es_cost + iter_val(iter) - abs->pi[ed];
      if(v < fwd_dist[ed])
      {
        fwd_dist[ed] = v;
        if(heap.inHeap(ed))
        {
          heap.decrease(ed);
        } 
        else {
          heap.insert(ed);
        }
      }
    }
  }
}

// Single-dest shortest path algorithm.
// Yes, I should probably templatify it so I can just swap out the iterators.
// But it's not worth the effort for now.
void update_closure(dbm x, int ii, int jj, int c)
{
  assert(x);
  assert(x->checked && x->feasible);
  // scratch holds the distances so far.
  assert(scratch_sz >= x->sz);
    
  // There may be a cheaper way to do this.
  if(dest_is_live(x, ii))
  {
    edge_iter piter = pred_iterator(x, ii);
    for(; !preds_end(piter); next_pred(piter))
    {
      int se = pred(piter);
      int sval = rev_iter_val(piter);

      assert(src_is_live(x, se));
      if(se != jj)
      {
        if(in_graph(x, se, jj))
        {
          x->mtx[se*x->sz + jj].val = std::min(edge_val(x, se, jj), sval + c);
        } else {
          dbm_add_edge(x, se, jj, sval + c);
        }
      
        if(src_is_live(x, jj))
        {
          edge_iter siter = src_iterator(x, jj);
          for(; !dests_end(siter); next_dest(siter))
          {
            int de = dest(siter);
            
            // We know src_is_live(x, src).
            if(se != de)
            {
              int val = rev_iter_val(piter) + c + iter_val(siter);
              if(in_graph(x, se, de))
              {
                x->mtx[se*x->sz + de].val = std::min(edge_val(x, se, de), val);
              } else {
                dbm_add_edge(x, se, de, val);
              }
            }
          }
        }
      }
    }
  }

  if(src_is_live(x, jj))
  {
    edge_iter siter = src_iterator(x, jj);
    for(; !dests_end(siter); next_dest(siter))
    {
      int de = dest(siter);
      int val = iter_val(siter) + c;
      if(de != ii)
      {
        if(in_graph(x, ii, de))
        {
          x->mtx[ii*x->sz + de].val = std::min(edge_val(x, ii, de), val);
        } else {
          dbm_add_edge(x, ii, de, val);
        }
      }
    }
  }
  // Closure is now updated.
  x->closed = true; 
}

// The transitive closure follows the Johnson's algorithm which finds
// efficiently for sparse graphs all shortest paths between all pairs
// allowing negative edges but *no* negative cycles. It uses
// Bellman-ford (check_feasibility) to check that there is no negative
// cycles and update the potential functions. Then, we transform the
// graph using those potential functions to remove negative edges so
// we can use Dijkstra's algorithm.
void dbm_canonical(dbm abs)
{
  if(!abs)
    return;
  
  if(!abs->checked)
    check_feasibility(abs);
  if(!abs->feasible)
    return;
  
  verify_potentials(abs);
  // Iterate through the live sources, and use Dijkstra's
  // algorithm on each.
  // At this point, abs->pi should be up to date.
  edge_iter iter = edge_iterator(abs);
  std::vector< std::vector< std::pair<int, int> > > changed;
  for(; !srcs_end(iter); next_src(iter))
  {
    changed.push_back( std::vector< std::pair<int, int> >() );
    dijkstra(abs, src(iter), changed.back());
  }
  for(int si = 0; si < abs->num_srcs; si++)
  {
    int s = abs->live_srcs[si];
    std::vector< std::pair<int, int> >& s_changed(changed[si]);
    for(unsigned int ei = 0; ei < s_changed.size(); ei++)
    {
      int d = s_changed[ei].first;
      int v = s_changed[ei].second;
      if(!in_graph(abs, s, d))
      {
        dbm_add_edge(abs, s, d, v);
      } else {
        abs->mtx[s*abs->sz + d].val = v;
      }
    }
  }
  abs->closed = true;
  
  return;
}

// out = x /\ y (destructive if out != null)
dbm dbm_meet(dbm x, dbm y)
{
  if(!x || !y)
    return NULL;
  assert(x->sz == y->sz);

  // Collect the min of all edges
  assert(x->sz == y->sz);
  if(!x->checked || !x->closed)
    dbm_canonical(x);
  if(!y->checked || !y->closed)
    dbm_canonical(y);

  if(!x->feasible || !y->feasible)
    return NULL;

  assert(x->checked && x->feasible);
  assert(y->checked && y->feasible);

  dbm ret = dbm_alloc(x->sz);

  bool xchange = false;
  bool ychange = false;

  edge_iter iter = edge_iterator(x);
  for(; !edges_end(iter); next_edge(iter))
  {
    dbm_add_edge(ret, src(iter), dest(iter), iter_val(iter));
    if(!src_is_live(y, src(iter)) || !in_graph(y, src(iter), dest(iter)))
      ychange = true;
  }

  edge_iter yiter = edge_iterator(y);
  for(; !edges_end(yiter); next_edge(yiter))
  {
    int i = src(yiter);
    int j = dest(yiter);

    if(src_is_live(ret, i) && in_graph(ret, i, j))
    {
      //ret->mtx[i*ret->sz + j].val = std::min(edge_val(ret, i, j), iter_val(yiter));
      int vx = edge_val(ret, i, j);
      int vy = iter_val(yiter);
      if(vy < vx)
      {
        xchange = true;
        ret->mtx[i*ret->sz + j].val = vy;
      } else if(vx < vy) {
        ychange = true;
      }
    } else {
      xchange = true;
      dbm_add_edge(ret, i, j, iter_val(yiter));
    }
  }
  
  if((x->closed && !xchange))
  {
    memcpy(ret->pi, x->pi, sizeof(int)*x->sz);
    ret->checked = true;
    ret->feasible = true;
    ret->closed = true;
  } else if((y->closed && !ychange)) {
    memcpy(ret->pi, y->pi, sizeof(int)*y->sz);
    ret->checked = true;
    ret->feasible = true;
    ret->closed = true;
  } else {
    ret->checked = false;
    ret->closed = false;
  }
  /*
  ret->checked = false;
  ret->closed = false;
  */

  return ret;
}

void dbm_print(dbm x)
{
  return dbm_print_to(std::cout, x);
}

dbm dbm_join(dbm x, dbm y)
{
  if(!x)
    return dbm_copy(y);
  if(!y)
    return dbm_copy(x);

  assert(x->sz == y->sz);
  if(!x->checked || !x->closed)
    dbm_canonical(x);
  if(!y->checked || !y->closed)
    dbm_canonical(y);

  if(!x->feasible)
    return dbm_copy(y);
  if(!y->feasible)
    return dbm_copy(x);

  dbm ret = dbm_alloc(x->sz);

  memcpy(ret->pi, x->pi, sizeof(int)*x->sz);
  
  edge_iter xiter = edge_iterator(x);
  for(; !srcs_end(xiter); next_src(xiter))
  {
    int i = src(xiter);
    if(!src_is_live(y, i))
      continue;

    while(!dests_end(xiter) && !in_graph(y, i, dest(xiter)))
      next_dest(xiter);

    for(; !dests_end(xiter); next_dest(xiter))
    {
      int j = dest(xiter);

      if(!in_graph(y, i, j))
        continue;

      int ij_val = std::max(iter_val(xiter), edge_val(y, i, j));

      dbm_add_edge(ret, i, j, ij_val);
    }
  }

  // The join of closed elements is also closed.
  ret->checked = true;
  ret->feasible = true;
  ret->closed = true;

  return ret;
}

// For DBM widening, you can close the left
// argument, but not the right.
// x is the previous iteration, y is the current.
dbm dbm_widen(dbm x, dbm y)
{
  /*
  if(!x || !y || dbm_is_top(y) || !x->feasible)
    return dbm_copy(y);
    */
  if(!y)
    return NULL;
  if(!y->checked || !y->closed)
    dbm_canonical(y);

  if(!y->feasible)
    return NULL;

  /*
  // No longer guarantees termination.
  if(!x->checked || !x->closed)
    dbm_canonical(x);
  */
  if(!x || (x->checked && !x->feasible))
    return dbm_copy(y);

  dbm ret = dbm_alloc(y->sz);
  assert(x->sz == y->sz);

  // Widen is an upper-bound function, so
  // the potential function should still
  // be fine.
  memcpy(ret->pi, x->pi, sizeof(int)*x->sz);
  
  // Walk through all the edges of y, and throw
  // away anything that isn't stable in x.
  edge_iter yiter = edge_iterator(y);
  for(; !srcs_end(yiter); next_src(yiter))
  {
    int i = src(yiter);
    // The source should be live in x, assuming
    // we've been moving up.
    if(!src_is_live(x, i))
      continue;

    // Skip past anything 
    while(!dests_end(yiter) && !in_graph(x, i, dest(yiter)))
      next_dest(yiter);

    for(; !dests_end(yiter); next_dest(yiter))
    {
      int j = dest(yiter);

      int yval = edge_val(y, i, j);

//    if(!in_graph(x, i, j))
//      continue;
      if(!in_graph(x, i, j)){
        dbm_add_edge(ret, i, j, yval);
      }
      if(in_graph(x, i, j)){ // FIX 
        int xval = edge_val(x, i, j);
        if(yval <= xval){
          dbm_add_edge(ret, i, j, xval);
        }
      }
    }
  }

  // The join of closed elements is also closed.
  ret->checked = true;
  ret->feasible = true;
  ret->closed = false;
  return ret;
}

// Existentially project a variable.
dbm dbm_forget(int v, dbm x)
{
  if(!x)
    return x;
  dbm_canonical(x);
  if(!x->feasible)
    return NULL;   
  int sz = x->sz;
  dbm ret = dbm_alloc(x->sz);
  memcpy(ret->pi, x->pi, sizeof(int)*sz);

  ret->sz = sz;
  ret->checked = true;
  ret->feasible = true;
  ret->closed = true;


  edge_iter iter = edge_iterator(x);
  for(; !srcs_end(iter); next_src(iter))
  {
    int s = src(iter);
    if(s == v)
      continue;

    for(; !dests_end(iter); next_dest(iter))
    {
      int d = dest(iter);
      if(d == v)
        continue;
      assert(!src_is_live(ret, s) || !in_graph(ret, s, d));
      dbm_add_edge(ret, s, d, iter_val(iter));
      /*
      if(in_graph(ret, s, d))
      {
        ret->mtx[s*ret->sz + d].val = std::min(edge_val(ret, s, d), iter_val(iter));
      } else {
        dbm_add_edge(ret, s, d, iter_val(iter));
      }
      */
    }
  }

  return ret;
}

dbm dbm_forget_array(int* vs, int vs_len, dbm x)
{
  if(!x)
    return x;
  dbm_canonical(x);
  if(!x->feasible)
    return NULL;   

  // Mark the var flags.
  for(int vi = 0; vi < vs_len; vi++)
    var_flags[vs[vi]] = 1;

  int sz = x->sz;
  dbm ret = dbm_alloc(x->sz);
  memcpy(ret->pi, x->pi, sizeof(int)*sz);

  ret->sz = sz;
  ret->checked = true;
  ret->feasible = true;
  ret->closed = true;

  edge_iter iter = edge_iterator(x);
  for(; !srcs_end(iter); next_src(iter))
  {
    int s = src(iter);
    if(var_flags[s])
      continue;

    for(; !dests_end(iter); next_dest(iter))
    {
      int d = dest(iter);
      if(var_flags[d])
        continue;
      assert(!src_is_live(ret, s) || !in_graph(ret, s, d));
      dbm_add_edge(ret, s, d, iter_val(iter));
    }
  }

  // Clear the var flags.
  for(int vi = 0; vi < vs_len; vi++)
    var_flags[vs[vi]] = 0;

  return ret;
}

dbm dbm_extract(int* vs, int vs_len, dbm x)
{
  if(!x)
    return x;
  dbm_canonical(x);
  if(!x->feasible)
    return NULL;   

  // Mark the var flags.
  for(int vi = 0; vi < vs_len; vi++)
    var_flags[vs[vi]] = 1;

  int sz = x->sz;
  dbm ret = dbm_alloc(x->sz);
  memcpy(ret->pi, x->pi, sizeof(int)*sz);

  ret->sz = sz;
  ret->checked = true;
  ret->feasible = true;
  ret->closed = true;

  edge_iter iter = edge_iterator(x);
  for(; !srcs_end(iter); next_src(iter))
  {
    int s = src(iter);

    for(; !dests_end(iter); next_dest(iter))
    {
      int d = dest(iter);
      if(!var_flags[s] && !var_flags[d])
        continue;

      assert(!src_is_live(ret, s) || !in_graph(ret, s, d));
      dbm_add_edge(ret, s, d, iter_val(iter));
    }
  }

  // Clear the var flags.
  for(int vi = 0; vi < vs_len; vi++)
    var_flags[vs[vi]] = 0;

  return ret;
}

// Copy x, but apply variable substitutions in subs.
// A variable may occur as each a source and a destination
// at most once.
// If a variable occurs as a destination but not a source,
// existing constraints are forgotten.
dbm dbm_rename(rmap* subs, int slen, dbm x)
{
  if(!x)
    return NULL;

  dbm_canonical(x); 
  if(!x || !x->feasible)
    return NULL;

  int sz = x->sz;
  dbm ret = dbm_alloc(x->sz);

  // Set up the mapping.
  // Abusing _gamma for temporary storage.
  bool changed = false;
  for(int vi = 0; vi < sz; vi++)
  {
    _gamma[vi] = vi;
  }
  for(int si = 0; si < slen; si++)
  {
    if(subs[si].r_to == subs[si].r_from)
      continue;

    changed = true;
    var_flags[subs[si].r_to] = 1; // Possibly overwritten
    _gamma[subs[si].r_from] = subs[si].r_to;
  }
  if(!changed)
  {
    for(int vi = 0; vi < sz; vi++)
      var_flags[vi] = 0;
    return dbm_copy(x);
  }

  for(int vi = 0; vi < sz; vi++)
  {
    if(_gamma[vi] != vi)
      var_flags[vi] = 0;
  }

  memcpy(ret->pi, x->pi, sizeof(int)*sz);
  // Permute the value of reassigned variables.
  for(int si = 0; si < slen; si++)
  {
    ret->pi[subs[si].r_to] = x->pi[subs[si].r_from];
  }

  // Copy edges with var_flags[{source, dest}] = 0.
  edge_iter iter = edge_iterator(x);
  for(; !srcs_end(iter); next_src(iter))
  {
    int i = src(iter);
    if(var_flags[i])
      continue;

    for(; !dests_end(iter); next_dest(iter))
    {
      int j = dest(iter);
      if(var_flags[j])
        continue;

      dbm_add_edge(ret, _gamma[i], _gamma[j], iter_val(iter));
    }
  }
  for(int vi = 0; vi < sz; vi++)
    var_flags[vi] = 0;

  ret->checked = true;
  ret->feasible = true;
  ret->closed = true;

  return ret;
}

// As dbm_rename, but:
// (1) Variables which do not occur as destinations are eliminated
// (2) A variable may occur as a source multiple times. In this case,
// phi [x -> y, x -> z] == phi [x->y] TT (y = z).
dbm dbm_rename_strict(rmap* subs, int slen, dbm x)
{
  // assert(0 && "dbm: dbm_rename_strict not implemented yet"); 

  if(!x)
    return NULL;

  dbm_canonical(x); 
  if(!x || !x->feasible)
    return NULL;

  int sz = x->sz;
  dbm ret = dbm_alloc(x->sz);

  memcpy(ret->pi, x->pi, sizeof(int)*sz);
  // Permute the value of reassigned variables.
  for(int si = 0; si < slen; si++)
  {
    ret->pi[subs[si].r_to] = x->pi[subs[si].r_from];
  }
  return ret;

}

// The graph is bottom iff it's infeasible.
int dbm_is_bottom(dbm x)
{
  if(!x)
    return true;

  if(x->checked)
    return !x->feasible;
  return !check_feasibility(x);
}

// The graph is top iff there are no edges.
int dbm_is_top(dbm x)
{
  if(!x)
    return false;

  return x->num_srcs == 0;
}

int dbm_is_leq(dbm x, dbm y)
{
  // FIXME: Need to normalize before doing the cutoff.
  dbm_canonical(x);
  dbm_canonical(y);

  if(dbm_is_bottom(x) || dbm_is_top(y))
    return true;
  if(dbm_is_top(x) || dbm_is_bottom(y))
    return false;
  
  //  dbm_print_to(std::cout, x);
  //  dbm_print_to(std::cout, y);

  edge_iter yiter = edge_iterator(y);
  for(; !edges_end(yiter); next_edge(yiter))
  {
    int s = src(yiter);
    int d = dest(yiter);

    if(!src_is_live(x, s) || !in_graph(x, s, d) || edge_val(x, s, d) > iter_val(yiter))
      return false;
  }
  return true;
}

// We can also do this without
// computing the closure.
int dbm_implies(dbm x, dexpr con)
{
  dbm_canonical(x);
  if(!x || !x->feasible)
    return true;

  // Convert the ucon to an edge.
  // Factor this out.
  int i;
  int j;
  int w;

  switch(con.kind)
  {
    case D_LB:
      i = con.args[0];
      j = x->sz-1;
      w = con.konst;
      break;
    case D_UB:
      i = x->sz-1;
      j = con.args[0];
      w = -con.konst;
      break;
    default: // D_DIFF
      i = con.args[0];
      j = con.args[1];
      w = con.konst;
      break;
  }
  if(i == j)
    return w >= 0;

  return src_is_live(x, i) && in_graph(x, i, j) && edge_val(x, i, j) <= w;
}

void dbm_print_to(std::ostream &o, dbm x)
{

  dbm_canonical(x);

  if(!x || !x->feasible)
  {
    o << "_|_" << std::endl;
    return;
  }

#if 0
  o <<"pi: [";
  if(x->sz > 0)
    o << x->pi[0];
  for(int vi = 1; vi < x->sz; vi++)
  {
    o << ", " << x->pi[vi];
  }
  o << "]" << std::endl;
#endif 

  edge_iter iter = edge_iterator(x);
  for(; !srcs_end(iter); next_src(iter))
  {
    int ii = src(iter);
    o << ii << " ->";    
    for(; !dests_end(iter); next_dest(iter))
    {
      int jj = dest(iter);
      o << " {" << jj << " : " << edge_val(x, ii, jj) << "}";
    }
    o << std::endl;
  }
}


// Preconditions:
// - x is in canonical form
// - x is feasible
static int var_lb(dbm x, int i)
{
  int v0 = x->sz - 1;
  if(src_is_live(x, v0) && in_graph(x, v0, i))
  {
    return -1*edge_val(x, v0, i); 
  } else {
    return VAL_MIN;
  }
}

static int var_ub(dbm x, int i)
{
  int v0 = x->sz - 1;
  if(src_is_live(x, i) && in_graph(x, i, v0))
  {
    return edge_val(x, i, v0);
  } else {
    return VAL_MAX;
  }
}

void exp_collect_vars(exp_t e, std::vector<int>& vs)
{
  switch(e->kind)
  {
    case E_VAR:
      if(!var_flags[(size_t) e->args[0]])
      {
        var_flags[(size_t) e->args[0]] = 1;
        vs.push_back((size_t) e->args[0]);
      }
      break;
    case E_CONST:
      break;
    default:
      exp_collect_vars((exp_t) e->args[0], vs);
      exp_collect_vars((exp_t) e->args[1], vs);
      break;
  }
}

typedef struct {
  int lb;
  int ub;
} interval;

int max(int x, int y, int z, int t) {
  return std::max(x, std::max(y, std::max(z,t)));
}

int min(int x, int y, int z, int t) {
  return std::min(x, std::min(y, std::min(z,t)));
}

// cx + k
typedef struct {
  interval coeff;
  interval konst;
} linterm;

bool is_zero(interval i)
{
  return i.lb == 0 && i.ub == 0;
}

bool is_top(interval i)
{
  return i.lb == VAL_MIN && i.ub == VAL_MAX;
}

bool is_bot(interval i)
{
  return i.lb == VAL_MAX && i.ub == VAL_MIN;
}

bool is_plus_infinite(interval i)
{
  return (i.ub == VAL_MAX);
}

bool is_minus_infinite(interval i)
{
  return (i.lb == VAL_MIN);
}

bool is_bounded(interval i)
{
  return (!is_top(i) && !is_minus_infinite(i) && !is_plus_infinite(i));
}

bool is_singleton(interval i)
{
  return i.lb == i.ub;
}

interval mk_interval(int lb, int ub)
{
  interval r = { lb, ub };
  return r;
}

interval mk_top()
{
  interval r = { VAL_MIN, VAL_MAX };
  return r;
}

interval mk_bot()
{
  interval r = {VAL_MAX,  VAL_MIN};
  return r;
}

void print_interval_to(std::ostream &o, interval i){
  o << "[" << i.lb << "," << i.ub << "]";
}

void print_linterm_to(std::ostream &o, linterm term){
  print_interval_to(o, term.coeff);
  o << "x +";
  print_interval_to(o, term.konst);
}

// pre: i is not zero
interval trim_zero(interval i)
{
  if (i.lb == 0)
    return mk_interval(i.lb + 1, i.ub);
  else if (i.ub == 0)
    return mk_interval(i.lb, i.ub - 1);
  else 
    return i;
}

interval i_const(int v)
{
  interval r = { v, v };
  return r;
}

interval i_neg(interval x)
{
  interval r = { -x.ub, -x.lb };
  return r;
}

interval i_add(interval x, interval y)
{
  interval r = {
    (x.lb == VAL_MIN || y.lb == VAL_MIN) ? VAL_MIN : x.lb + y.lb,
    (x.ub == VAL_MAX || y.ub == VAL_MAX) ? VAL_MAX : x.ub + y.ub
  };
  return r;
}

interval i_mul(interval x, interval y)
{

  int ll = x.lb * y.lb;
  int lu = x.lb * y.ub;
  int ul = x.ub * y.lb;
  int uu = x.ub * y.ub;
  
  int lb = (x.lb == VAL_MIN || y.lb == VAL_MIN) ? VAL_MIN : min(ll,lu,ul,uu);
  int ub = (x.ub == VAL_MAX || y.ub == VAL_MAX) ? VAL_MAX : max(ll,lu,ul,uu);
  interval r = { lb, ub };
  return r; 
}

interval i_div(interval x, interval y)
{
  // assert(!is_zero(y) && "dbm: interval cannot be [0,0]"); 
  if (is_zero (y))
    return mk_bot ();

  interval y1 = trim_zero(y);

  if ( (!is_singleton(y1)) & (y1.lb <= 0) && (0 <= y1.ub)){
    // we give up if 0 \in y
    return mk_top();
  }
  else{
    int ll = x.lb / y1.lb;
    int lu = x.lb / y1.ub;
    int ul = x.ub / y1.lb;
    int uu = x.ub / y1.ub;
    
    int lb = (x.lb == VAL_MIN || y1.lb == VAL_MIN) ? VAL_MIN : min(ll,lu,ul,uu);
    int ub = (x.ub == VAL_MAX || y1.ub == VAL_MAX) ? VAL_MAX : max(ll,lu,ul,uu);
    interval r = { lb, ub };
    return r; 
  }
}

bool eval_exp(dbm d, int v, exp_t e, linterm& term)
{
  linterm left;
  linterm right;
  switch(e->kind)
  {
    case E_VAR:
      if((size_t) e->args[0] == v)
      {
        term.coeff = i_const(1);
        term.konst = i_const(0);
      } 
      else 
      {
        term.coeff = i_const(0);
        term.konst = mk_interval(var_lb(d, (size_t) e->args[0]),var_ub(d, (size_t) e->args[0]));
      }
      return true;
      break;
    case E_CONST:
      term.coeff = i_const(0);
      term.konst = i_const((size_t) e->args[0]);
      return true;
      break;
    default:
      break;
  }
  if(!eval_exp(d, v, (exp_t) e->args[0], left))
    return false;
  if(!eval_exp(d, v, (exp_t) e->args[1], right))
    return false;
  
  switch(e->kind)
  {
    case E_ADD:
      // (ax + b) + (cx + d) = (a+c)x + (b+d)
      term.coeff = i_add(left.coeff, right.coeff);
      term.konst = i_add(left.konst, right.konst);
      break;
    case E_SUB:
      // (ax + b) - (cx + d) = (a-c)x + (b-d)
      term.coeff = i_add(left.coeff, i_neg(right.coeff));
      term.konst = i_add(left.konst, i_neg(right.konst));
      break;
    case E_MUL:
      if(is_zero(left.coeff))
      {
        // b * (cx + d) = (b*c)x + (b*d)
        term.coeff = i_mul(left.konst, right.coeff); 
        term.konst = i_mul(left.konst, right.konst);
      } 
      else if(is_zero(right.coeff)) 
      {
        // (ax + b) * d = (a*d)x + (b*d)
        term.coeff = i_mul(left.coeff, right.konst);
        term.konst = i_mul(left.konst, right.konst); 
      } 
      else 
      {
        // the resulting term is not linear
        return false;
      }
      //if (is_top(term.coeff) || is_top(term.konst))
      if (!is_bounded(term.coeff) || !is_bounded(term.konst))
      {
        return false;
      }
      break;
    case E_DIV:
      if(is_zero(right.coeff)) 
      {
        // (ax + b) / d = (a/d)*x + (b/d)
        term.coeff = i_div(left.coeff, right.konst);
        term.konst = i_div(left.konst, right.konst); 
      }
      else
      {
        return false;
      }
      if (!is_bounded(term.coeff) || !is_bounded(term.konst))
      {
        return false;
      }
      break;
    default:
      assert(0 && "dbm: unreachable");
  }

  // FIXME
  // if ( (term.coeff.lb == VAL_MIN || term.coeff.ub == VAL_MAX) ||
  //      (term.konst.lb == VAL_MIN || term.konst.ub == VAL_MAX) )
  //   return false;
  // else 
  return true;
}

// typedef struct {
//   int val;
//   int is_bot;
// } pi_t;

// pi_t mk_pi(int val){
//   pi_t r = { val, 0 /*false*/};
//   return r;
// }

// pi_t mk_bot_pi(){
//   pi_t r = { 0 /*any val here*/, 1 /*true*/};
//   return r;
// }

// Compute a potential value for a freshly introduced variable.
int eval_pi_rec(dbm x, exp_t e, int v0, int* is_bot)
{

  *is_bot = 0;
  switch(e->kind)
  {
    case E_CONST:
      return ((size_t) e->args[0]);
    case E_VAR:
      return v0 - x->pi[(size_t) e->args[0]];
    case E_ADD:
      {
        int is_bot1, is_bot2;
        int pi = eval_pi_rec(x, (exp_t) e->args[0], v0, &is_bot1) 
            + eval_pi_rec(x, (exp_t) e->args[1], v0, &is_bot2);
        *is_bot = is_bot1 + is_bot2;
        return pi;
      }
    case E_SUB:
      {
        int is_bot1, is_bot2;
        int pi = eval_pi_rec(x, (exp_t) e->args[0], v0, &is_bot1) 
            - eval_pi_rec(x, (exp_t) e->args[1], v0, &is_bot2);
        *is_bot = is_bot1 + is_bot2;
        return pi;
      }
    case E_MUL:
      {
        int is_bot1, is_bot2;
        int pi = eval_pi_rec(x, (exp_t) e->args[0], v0, &is_bot1) 
            * eval_pi_rec(x, (exp_t) e->args[1], v0, &is_bot2);
        *is_bot = is_bot1 + is_bot2;
        return pi;
      }
    case E_DIV:
      {
        int d = eval_pi_rec(x, (exp_t) e->args[1], v0, is_bot);
        //assert(d != 0 && "dbm: division by 0 in pi");
        if (*is_bot || d == 0)
        {
          *is_bot = 1;
          return VAL_MIN;
        }
        else
        {
          return eval_pi_rec(x, (exp_t) e->args[0], v0, is_bot) / d;
        }
      }
    default:
      assert(0 && "dbm: unreachable"); 
      return VAL_MIN;
  }
}

int eval_pi(dbm x, exp_t e, int* is_bot)
{
  int v0 = x->pi[x->sz - 1];
  return v0 - eval_pi_rec(x, e, v0, is_bot);
}

void post_linterm(dbm x, int v, int y, linterm t)
{
  if(y < 0)
  {
    // Bounds on v.
    int v0 = x->sz - 1;

    if (!is_zero(t.coeff)){
      std::cerr << "DBM ERROR: linear term ";
      print_linterm_to(std::cerr , t);
      std::cerr << " should have coefficient zero" << std::endl;
      assert(0);
    }

    if(t.konst.lb != VAL_MIN)
      dbm_add_edge(x, v0, v, -t.konst.lb);
    if(t.konst.ub != VAL_MAX)
      dbm_add_edge(x, v, v0, t.konst.ub);
  }
  else {
    // Should be able to compute some bounds for other coeffs.
    assert(t.coeff.lb == 1 && t.coeff.ub == 1);
    if(t.konst.lb != VAL_MIN)
      dbm_add_edge(x, y, v, -t.konst.lb);
    if(t.konst.ub != VAL_MAX)
      dbm_add_edge(x, v, y, t.konst.ub);
  }
}

dbm dbm_assign(int v, exp_t e, dbm x)
{
  if (!x || (x->checked && !x->feasible))
    return NULL;

  // std::cout << "dbm_assign BEGIN: \n";
  // dbm_print_to(std::cout, x);

  dbm ret = dbm_forget(v, x);
  if (!ret)
    return NULL;

  // std::cout << "After forgetting " << v  << std::endl;
  // dbm_print_to(std::cout, ret);
      
  std::vector<int> ys;
  exp_collect_vars(e, ys);

  linterm term;
  if (!eval_exp(x, -1, e, term))
    return ret;

  if (is_bot (term.coeff) || is_bot (term.konst))
  {
    return dbm_bottom ();
  }
  else
  {
    post_linterm(ret, v, -1, term);
  }

  for(unsigned int yi = 0; yi < ys.size(); yi++)
  {
    int y = ys[yi];
    var_flags[y] = 0;
    if(eval_exp(x, y, e, term))
    {
      // v = [c_l, c_u].y + [k_l, k_u]
      if (is_bot (term.coeff) || is_bot (term.konst))
      {
        return dbm_bottom ();
      }
      else if(term.coeff.lb == 1 && term.coeff.ub == 1)
      {
        post_linterm(ret, v, y, term);
      }
    }
  }

  // std::cout << "After applying assignment " << std::endl;
  // dbm_print_to(std::cout, ret);


  if(x->checked)
  {
    ret->checked = true;
    ret->feasible = true;
    int is_bot;
    int pi = eval_pi(x, e, &is_bot);
    if (is_bot)
    {
      return dbm_bottom ();
    }
    
    ret->pi[v] = pi;
    verify_potentials(ret);

#if 0
    if(x->closed)
    {
      // We can update the closure by following edges adjacent to v.
      if(src_is_live(ret, v))
      {
        edge_iter viter = src_iterator(ret, v);
        for(; !dests_end(viter); next_dest(viter))
        {
          int k = dest(viter);
          int kval = iter_val(viter);
          if(!src_is_live(x, k))
            continue;
          edge_iter diter = src_iterator(x, k);
          for(; !dests_end(diter); next_dest(diter))
          {
            int d = dest(diter);
            int dval = iter_val(diter);
            if(v == d)
              continue;
            if(in_graph(ret, v, d))
            {
              ret->mtx[v*x->sz + d].val = std::min(edge_val(ret, v, d), kval + dval);
            }
            else 
            {
              dbm_add_edge(ret, v, d, kval + dval);
            }
          }
        }

        viter = pred_iterator(ret, v);
        for(; !preds_end(viter); next_pred(viter))
        {
          int k = pred(viter);
          int kval = rev_iter_val(viter);
          if(!dest_is_live(x, k))
              continue;
          edge_iter siter = pred_iterator(x, k);
          for(; !preds_end(siter); next_pred(siter))
          {
            int s = pred(siter);
            int sval = rev_iter_val(siter);
            if(s == v)
              continue;
            if(src_is_live(ret, s) && in_graph(ret, s, v))
            {
              ret->mtx[s*x->sz + v].val = std::min(edge_val(ret, s, v), kval + sval);
            } 
            else 
            {
              dbm_add_edge(ret, s, v, kval + sval);
            }
          }
        }
      }
      // std::cout << "After transitive closure: \n" << std::endl;
      // dbm_print_to(std::cout, ret);

      verify_potentials(ret);
      ret->closed = true;
    }
#else
  dbm_canonical(ret);
#endif 
  }
  return ret;
}

dbm dbm_store(int a, uterm i, uterm v, dbm x)
{
  return dbm_copy(x);
}

dbm dbm_load(int v, int a, uterm i, dbm x)
{
  return dbm_forget(v, x);
}

dbm dbm_apply_edge(dbm x, int i, int j, int w)
{
  if (!x || (x->checked && !x->feasible))
  {
    return NULL;
  }

  dbm ret = NULL;

  // Edge is added, but isn't redundant.
  if(x->checked)
  {
    verify_potentials(x);
    if(!update_feasibility(x, i, j, w)){
      //assert(!src_is_live(x, i) || !in_graph(x, i, j) || edge_val(x, i, j) > w);
      ///// Updated constraint isn't feasible.
      ///// fprintf(stderr, "PING: B (%d -> %d: %d.\n", i, j, w);
      return NULL;
    } 
    else {
      // Copy x, and update pi.
      ret = dbm_copy(x);
      for(int vi = 0; vi < ret->sz; vi++)
        ret->pi[vi] = pi_prime[vi];
      ret->checked = true;
      ret->feasible = true;
      verify_potentials(ret);
    }
  } 
  else {
    ret = dbm_copy(x);
  }

  assert(ret);

  // Check if the edge is (i, w, j) is redundant.
  if(src_is_live(ret, i) && in_graph(ret, i, j))
  {
    // Redundant; no need to do anything.
    if(edge_val(ret, i, j) <= w)
      return ret;
    ret->mtx[i*ret->sz + j].val = w;
  } else {
    dbm_add_edge(ret, i, j, w);
  }
  if(x->checked)
    verify_potentials(ret);

  // If x is closed, we can update the closure.
  if(x->checked && x->closed)
  {
    update_closure(ret, i, j, w);
    verify_potentials(ret);
  }

  assert(ret);
  return ret;
}

dbm dbm_cond(ucon con, dbm x)
{
  if(!x || (x->checked && !x->feasible))
    return NULL;

  if(con.kind == U_DIS){
    if(con.args[0].kind == U_CONST && con.args[1].kind == U_CONST)
    {
      if (con.args[0].val != con.args[1].val)
        return dbm_copy(x);
      else 
        return NULL;
    } 
    else {
      // TODO: x != k when x >= k1 or x <= k2 where k,k1,k2 are constants
      return dbm_copy(x);
    }
  }

  int i;
  int j;
  int w;
  
  // assert(con.kind == U_LEQ || con.kind == U_LT);

  // C_1 op C_2.
  if(con.args[0].kind == U_CONST && con.args[1].kind == U_CONST)
  {
    if((con.kind == U_LEQ && con.args[0].val <= con.args[1].val) ||
       (con.kind == U_LT && con.args[0].val < con.args[1].val)){
      return dbm_copy(x);
    } 
    else {
      return NULL;
    }
  } 
  else if(con.args[0].kind == U_CONST) {
    // C_0 op V
    i = x->sz - 1;
    j = con.args[1].val;
    w = -1*con.args[0].val;
  } 
  else if(con.args[1].kind == U_CONST) {
    // V op C_0
    i = con.args[0].val;
    j = x->sz - 1;
    w = con.args[1].val;
  } 
  else {
    // V_1 op V_2
    i = con.args[0].val;
    j = con.args[1].val;
    w = 0;
    if(i == j){
      if(con.kind == U_LEQ)
        return dbm_copy(x);
      else
        return NULL;
    }
  }
  if(con.kind == U_LT)
    w -= 1;

  /*
  fprintf(stderr, "Cond: [%d %s %d] :: (%d, %d, %d)\n",
      con.args[0].val, con.kind == U_LEQ ? "<=" : "<", con.args[1].val, 
      i, j, w);
  */
  if(con.kind == U_EQ){
    // To avoid leaking x1
    dbm tmp = dbm_apply_edge(x, i, j, w);
    dbm res = dbm_apply_edge(tmp, j, i, -w);
    dbm_dealloc(tmp);
    return res;
    //return dbm_apply_edge(dbm_apply_edge(x, i, j, w), j, i, -w);
  } 
  else {
    return dbm_apply_edge(x, i, j, w);
  }
}

dbm dbm_apply_dexpr(dexpr con, dbm x)
{
  dbm_canonical(x);
  if(!x || !x->feasible)
    return NULL;

  // Convert the ucon to an edge.
  // Factor this out.
  int i;
  int j;
  int w;

  switch(con.kind)
  {
    case D_LB:
      i = con.args[0];
      j = x->sz-1;
      w = con.konst;
      break;
    case D_UB:
      i = x->sz-1;
      j = con.args[0];
      w = -con.konst;
      break;
    default: // D_DIFF
      i = con.args[0];
      j = con.args[1];
      w = con.konst;
      break;
  }

  if(i == j)
  {
    if(w < 0)
      return NULL;
    else
      return dbm_copy(x);
  } 
  else {
    return dbm_apply_edge(x, i, j, w);
  }
}

// copy the variable v to a new one new_v.
dbm dbm_expand (int v, int new_v, dbm x)
{
  if (!x || (x->checked && !x->feasible))
    return NULL;

  rmap subs[1];
  subs[0].r_from = v;
  subs[0].r_to = new_v;

  dbm tmp = dbm_rename (subs, 1, x);
  dbm res = dbm_meet (x,tmp);
  dbm_dealloc(tmp);
  return res;
}
