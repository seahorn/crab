#ifndef CRAB_GRAPH_OPS_HPP
#define CRAB_GRAPH_OPS_HPP
#include <crab/domains/dbm/util/Heap.h>
//============================
// A set of utility algorithms for manipulating graphs.

using namespace std;

namespace crab {
  // Graph views - for when we want to traverse some mutation
  // of the graph without actually constructing it.
  // ============
  
  // Processing a graph under a (possibly incomplete)
  // permutation of vertices.
  // We assume perm[x] is unique; otherwise, we'd have
  // to introduce a edges for induced equivalence classes.
  template<class G>
  class GraphPerm {
  public: 
    typedef typename G::vert_id vert_id;
    typedef typename G::Wt Wt;
    typedef typename G::adj_list g_adj_list;

    GraphPerm(vector<vert_id>& _perm, G& _g) 
      : perm(_perm), g(_g)
    {
      vector<vert_id> inv(g.size(), -1);
      for(unsigned int vi = 0; vi < perm.size(); vi++)
      {
        assert(inv[perm[vi]] == -1);
        inv[perm[vi]] = vi;
      }
    }

    // Check whether an edge is live
    bool elem(vert_id x, vert_id y) const {
      if(perm[x] < 0 || perm[y] < 0)
        return false;
      return g.elem(perm[x], perm[y]);
    }

    // Precondition: elem(x, y) is true.
    Wt& edge_val(vert_id x, vert_id y) const {
      assert(perm[x] >= 0 && perm[y] >= 0);
      return g.edge_val(perm[x], perm[y]);
    }

    // Precondition: elem(x, y) is true.
    Wt operator()(vert_id x, vert_id y) const {
      assert(perm[x] >= 0 && perm[y] >= 0);
      return g(perm[x], perm[y]);
    }

    // Number of allocated vertices
    int size(void) const {
      return perm.size();
    }

    class adj_iterator {
    public:
      adj_iterator(vector<vert_id>& _inv, vert_id* _v)
        : inv(_inv), v(_v)
      { }

      vert_id operator*(void)
      {
        while(inv[*v] < 0)
          ++v;
        return *v;
      }

      adj_iterator& operator++(void) {
        ++v;
        return this;
      }

      bool operator!=(adj_iterator& other)
      {
        return v != other.v;
      }
    protected:
      vector<vert_id>& inv;
      vert_id* v;
    };

    class adj_list {
    public:
      adj_list(vector<vert_id>& _inv, g_adj_list& _adj)
        : inv(_inv), adj(_adj)
      {
        vert_id* adj_end(adj.end());
        while(adj_end != adj.begin() && inv[adj_end[-1]] < 0)
          --adj_end;
      }
      adj_iterator begin(void) const { return adj_iterator(adj.begin()); } 
      adj_iterator end(void) const { return adj_iterator(adj_end); }

      bool mem(unsigned int v) const {
        if(perm[v] < 0)
          return false;
        return adj.mem(perm[v]); 
      }

    protected:
      vector<vert_id>& inv;
      g_adj_list adj;
      vert_id* adj_end;
    };

    adj_list succs(vert_id v)
    {
      return adj_list(inv, g.succs(perm[v]));
    }
    adj_list preds(vert_id v)
    {
      return adj_list(inv, g.preds(perm[v]));
    } 
    G& g;
    vector<vert_id> perm;
    vector<vert_id> inv;
  };

  // Viewing a graph with all edges reversed.
  // Useful if we want to run single-dest shortest paths,
  // for updating bounds and incremental closure.
  template<class G>
  class GraphRev {
  public: 
    typedef typename G::vert_id vert_id;
    typedef typename G::Wt Wt;
    typedef typename G::adj_list g_adj_list;

    GraphRev(vector<vert_id>& _perm, G& _g) 
      : g(_g)
    { }

    // Check whether an edge is live
    bool elem(vert_id x, vert_id y) const {
      return g.elem(y, x);
    }

    // Precondition: elem(x, y) is true.
    Wt& edge_val(vert_id x, vert_id y) const {
      return g.edge_val(y, x);
    }

    // Precondition: elem(x, y) is true.
    Wt operator()(vert_id x, vert_id y) const {
      return g(y, x);
    }

    // Number of allocated vertices
    int size(void) const {
      return g.size();
    }

    typedef typename G::adj_list adj_list;

    adj_list succs(vert_id v)
    {
      return g.preds(v);
    }
    adj_list preds(vert_id v)
    {
      return g.succs(v);
    } 
    G& g;
  };

  // Comparator for use with min-heaps.
  template<class V>
  class DistComp {
  public:
    DistComp(V& _A)
      : A(_A)
    { }
    bool operator() (int x, int y) const {
      return A[x] < A[y]; 
    }
    V& A;
  };

  // GKG - What's the best way to split this out?
  template<class Number>
  class GraphOps {
  public:
    typedef Number Wt;
    // The following code assumes vert_id is an integer.
    typedef SparseWtGraph<Wt> graph_t;
    typedef typename graph_t::vert_id vert_id;

    typedef vector< pair< pair<vert_id, vert_id>, Wt > > edge_vector;

    typedef DistComp< vector<Wt> > WtComp;
    typedef Heap<WtComp> WtHeap;

    // ===========================================
    // Scratch space needed by the graph algorithms.
    // Should really switch to some kind of arena allocator, rather
    // than having all these static structures.
    // ===========================================
    // Used for marking edge colour during meet
    enum MarkT { E_NONE = 0, E_LEFT = 1, E_RIGHT = 2, E_BOTH = 3 };
    static char* edge_marks;

    // Used for Bellman-Ford queueing
    static vert_id* dual_queue;
    static int* vert_marks;
    static unsigned int scratch_sz;

    // For locality, should combine dists & dist_ts.
    // Wt must have an empty constructor, but does _not_
    // need a top or infty element.
    // dist_ts tells us which distances are current,
    // and ts_idx prevents wraparound problems, in the unlikely
    // circumstance that we have more than 2^sizeof(uint) iterations.
    static vector<Wt> dists;
    static vector<unsigned int> dist_ts;
    static unsigned int ts; 
    static unsigned int ts_idx;

    static void grow_scratch(unsigned int sz) {
      if(sz <= scratch_sz)
        return;

      if(scratch_sz == 0)
        scratch_sz = 10; // Introduce enums for init_sz and growth_factor
      while(scratch_sz < sz)
        scratch_sz *= 1.5;

      edge_marks = (char*) realloc(edge_marks, sizeof(char)*scratch_sz*scratch_sz);
      dual_queue = (vert_id*) realloc(dual_queue, sizeof(vert_id)*2*scratch_sz);
      vert_marks = (int*) realloc(vert_marks, sizeof(int)*scratch_sz);

      // Initialize new elements as necessary.
      while(dists.size() < scratch_sz)
      {
        dists.push_back(Wt());
        dist_ts.push_back(ts-1);
      }
    }

    // Syntactic join.
    template<class G1, class G2>
    static graph_t join(G1& l, G2& r)
    {
      // For the join, potentials are preserved 
      assert(l.size() == r.size());
      int sz = l.size();

      graph_t g(sz); 
      
      for(vert_id s : l.verts())
      {
        for(vert_id d : l.succs(s))
        {
          if(r.elem(s, d))
            g.add_edge(s, max(l.edge_val(s, d), r.edge_val(s, d)), d);
        }
      }
      return g;
    }

    // Syntactic meet
    template<class G1, class G2>
    static graph_t meet(G1& l, G2& r)
    {
      assert(l.size() == r.size());
      int sz = l.size();

      graph_t g(l); 
      
      for(vert_id s : r.verts())
      {
        for(vert_id d : r.succs(s))
        {
          if(!l.elem(s, d))
            g.add_edge(s, r.edge_val(s, d), d);
          else
            g.edge_val(s, d) = min(g.edge_val(s, d), r.edge_val(s, d));
        }
      }
      return g;
    }

    // Compute the strongly connected components
    // Duped pretty much verbatim from Wikipedia
    // Abuses 'dual_queue' to store indices.
    template<class G>
    static void strong_connect(G& x, vector<vert_id>& stack, int& index, vert_id v, vector< vector<vert_id> >& sccs)
    {
      vert_marks[v] = (index<<1)|1;
      assert(vert_marks[v]&1);
      dual_queue[v] = index;
      index++;

      stack.push_back(v);

      // Consider successors of v
      for(vert_id w : x.succs(v))
      {
        if(!vert_marks[w])
        {
          strong_connect(x, stack, index, w, sccs);
          dual_queue[v] = min(dual_queue[v], dual_queue[w]);
        } else if(vert_marks[w]&1) {
          // W is on the stack
          dual_queue[v] = min(dual_queue[v], (vert_id) (vert_marks[w]>>1));
        }
      }

      // If v is a root node, pop the stack and generate an SCC
      if(dual_queue[v] == (vert_marks[v]>>1))
      {
        sccs.push_back(vector<vert_id>()); 
        vector<vert_id>& scc(sccs.back());
        int w;
        do 
        {
          w = stack.back();
          stack.pop_back();
          vert_marks[w] &= (~1);
          scc.push_back(w);
        } while (v != w);
      }
    }

    template<class G>
    static void compute_sccs(G& x, vector< vector<vert_id> >& out_scc)
    {
      int sz = x.size();
      grow_scratch(sz);

      for(vert_id v : x.verts())
        vert_marks[v] = 0;
      int index = 1;
      vector<vert_id> stack;
      for(vert_id v : x.verts())
      {
        if(!vert_marks[v])
          strong_connect(x, stack, index, v, out_scc);
      }
      /*
      printf("[");
      for(int ii = 0; ii < out_scc.size(); ii++)
      {
        printf("[");
        for(int jj = 0; jj < out_scc[ii].size(); jj++)
          printf(" %d", out_scc[ii][jj]);
        printf("]");
      }
      printf("]\n");
      */

      for(vert_id v : x.verts())
        vert_marks[v] = 0;
    }

    // Run Bellman-Ford to compute a valid model of a set of difference constraints.
    // Returns false if there is some negative cycle.
    template<class G, class P>
    static bool select_potentials(G& g, P& potentials)
    {
      int sz = g.size();
      assert(potentials.size() >= sz);
      grow_scratch(sz);

      vector< vector<vert_id> > sccs;
      compute_sccs(g, sccs);

      // Zero existing potentials.
      // Not strictly necessary, but means we're less
      // likely to run into over/underflow.
      //
      // Though hurts our chances of early cutoff.
      for(vert_id v : g.verts())
        potentials[v] = 0;

      // Run Bellman-ford on each SCC.
      //for(vector<vert_id>& scc : sccs)
      // Current implementation returns sccs in reverse topological order.
      for(auto it = sccs.rbegin(); it != sccs.rend(); ++it)
      {
        vector<vert_id>& scc(*it);

        vert_id* qhead = dual_queue;
        vert_id* qtail = qhead;

        vert_id* next_head = dual_queue+sz;     
        vert_id* next_tail = next_head;

        for(vert_id v : scc)
        {
          *qtail = v;
          vert_marks[v] = 3;
          qtail++;
        }

        for(int iter = 0; iter < scc.size(); iter++)
        {
          for(; qtail != qhead; )
          {
            vert_id s = *(--qtail); 
            vert_marks[s] = 2; 
            
            Wt s_pot = potentials[s];

            for(vert_id d : g.succs(s))
            {
              Wt sd_pot = s_pot + g.edge_val(s, d);
              if(sd_pot < potentials[d])
              {
                potentials[d] = sd_pot;
                if(vert_marks[d] == 2)
                {
                  *next_tail = d;
                  vert_marks[d] = 3;
                  next_tail++;
                }
              }
            }
          }
          // Prepare for the next iteration
          swap(qhead, next_head);
          qtail = next_tail;
          next_tail = next_head;
          if(qhead == qtail)
            break;
        }
        // Check if the SCC is feasible.
        for(; qtail != qhead; qtail--)
        { 
          vert_id s = *qtail;
          Wt s_pot = potentials[s];
          for(vert_id d : g.succs(s))
          {
            if(s_pot + g.edge_val(s, d) < potentials[d])
            {
              // Cleanup vertex marks
              for(vert_id v : g.verts())
                vert_marks[v] = 0;
              return false;
            }
          }
        }
      }
      return true;
    }

    template<class G, class B, class P>
    static bool update_bounds(G& g, B& bounds, P& pots)
    {
      // We compute upper bounds by running Dijkstra from (the imaginary) v0.
      // Lower bounds are obtained by reversing the graph (running a single-dest
      // shortest path problem).
       
      return true;
    }

    template<class G1, class G2, class P>
    static void close_after_meet(graph_t& g, P& pots, G1& l, G2& r, edge_vector& delta)
    {
      // We assume the syntactic meet has already been computed,
      // and potentials have been initialized.
      // We just want to restore closure.
      assert(l.size() == r.size());
      unsigned int sz = l.size();
      delta.clear();
      
      vector< vector<vert_id> > colour_succs(2*sz);
      // Partition edges into r-only/rb/b-only.
      for(vert_id s : g.verts())
      {
        unsigned int g_count = 0;
        unsigned int r_count = 0;
        for(vert_id d : g.succs(s))
        {
          char mark = 0;
          if(l.elem(s, d) && l.edge_val(s, d) == g.edge_val(s, d))
            mark |= E_LEFT;
          if(r.elem(s, d) && r.edge_val(s, d) == g.edge_val(s, d))
            mark |= E_RIGHT;
          // Add them to the appropriate coloured successor list 
          // Could do it inline, but this'll do.
          assert(mark != 0);
          switch(mark)
          {
            case E_LEFT:
              colour_succs[2*s].push_back(d);
              break;
            case E_RIGHT:
              colour_succs[2*s+1].push_back(d);
              break;
            default:
              break;
          }
          edge_marks[sz*s + d] = mark;
        }
      }

      // We can run the chromatic Dijkstra variant
      // on each source.
      for(vert_id v = 0; v < sz; v++)
      {
        delta.push_back(); 
        chrome_dijkstra(g, pots, v, delta.back());
      }

      // Should actually return this set.
      for(vert_id v = 0; v < sz; v++)
      {
        for(auto p : delta[v])
        {
          if(g.elem(v, p.first)) 
          {
            g.edge_val(v, p.first) = p.second;
          } else {
            g.add_edge(v, p.second, p.first);
          }
        }
      }
    }

    static void apply_delta(graph_t& g, edge_vector& delta)
    {
      for(vert_id v = 0; v < delta.size(); v++)
      {
        for(auto p : delta[v])
        {
          if(g.elem(v, p.first))
          {
            g.edge_val(v, p.first) = p.second;
          } else {
            g.add_dge(v, p.second, p.first);
          } 
        }
      }
    }

    // P is some vector-alike holding a valid system of potentials.
    // Don't need to clear/initialize 
    template<class P>
    static void chrome_dijkstra(graph_t& g, P& p, vector< vector<vert_id> >& colour_succs, vert_id src, vector< pair<vert_id, Wt> >& out)
    {
      unsigned int sz = g.size();
      if(sz == 0)
        return;
      grow_scratch(sz);

      // Reset all vertices to infty.
      dist_ts[ts_idx] = ts++;
      ts_idx = (ts_idx+1) % dists.size();

      dists[src] = Wt(0);
      dist_ts[src] = ts;

      WtComp comp(dists);
      WtHeap heap(comp);

      for(vert_id dest : succs(src))
      {
        dists[dest] = p[src] + g.edge_val(src, dest) - p[dest];
        dist_ts[dest] = ts;

        vert_marks[dest] = edge_marks[sz*src + dest];
        heap.insert(dest);
      }

      while(!heap.empty())
      {
        int es = heap.removeMin();
        int es_cost = dists[es] + p[es]; // If it's on the queue, distance is not infinite.
        int es_val = es_cost - p(src);
        if(!g.elem(src, es) || g.edge_val(src, es) > es_val)
          out.push_back( make_pair(es, es_val) );

//        if(!src_is_live(abs, es))
//          continue;

        if(vert_marks[es] == (E_LEFT|E_RIGHT))
          continue;

        // Pick the appropriate set of successors
        vector<vert_id>& es_succs = (vert_marks[es] == E_LEFT) ?
          colour_succs[2*es] : colour_succs[2*es+1];
        for(vert_id ed : es_succs)
        {
          int v = es_cost + g.edge_val(es, ed) - p[ed];
          if(dist_ts[ed] != ts || v < dists[ed])
          {
            dists[ed] = v;
            dist_ts[ed] = ts;
            vert_marks[ed] = edge_marks[sz*es+ed];

            if(heap.inHeap(ed))
            {
              heap.decrease(ed);
            } else {
              heap.insert(ed);
            }
          } else if(v == dists[ed]) {
            vert_marks[ed] |= edge_marks[sz*es+ed];
          }
        }
      }
    }

    template<class G, class P>
    static void close_after_widen(graph_t& g, G& orig, vector< vector< pair<vert_id, Wt> > >& delta)
    {
      unsigned int sz = g.size();
      assert(orig.size() == sz);
      
      for(vert_id v : g.verts())
      {
        vert_marks[v] = 0;
        // Assumption: stable iff |G(v)| = |H(v)|.
        if(g.succs(v).size() == orig.succs(v).size())
          vert_marks[v] = 1;
      }
    }
  };

  // Static data allocation
  template<class Wt>
  char* GraphOps<Wt>::edge_marks = NULL;

  // Used for Bellman-Ford queueing
  template<class Wt>
  typename GraphOps<Wt>::vert_id* GraphOps<Wt>::dual_queue = NULL;

  template<class Wt>
  int* GraphOps<Wt>::vert_marks = NULL;

  template<class Wt>
  unsigned int GraphOps<Wt>::scratch_sz = 0;

  template<class Wt>
  vector<Wt> GraphOps<Wt>::dists;
  template<class Wt>
  vector<unsigned int> GraphOps<Wt>::dist_ts;
  template<class Wt>
  unsigned int GraphOps<Wt>::ts = 0;
  template<class Wt>
  unsigned int GraphOps<Wt>::ts_idx = 0;

} // namespace crab

#endif
