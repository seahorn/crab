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
  template<class V>
  class range_iterator {
  public:
    class value_ref {
    public:
      value_ref(const V& _v)
        : v(_v)
      { }
      const V& operator*(void) const { return v; }
      value_ref& operator++(void) { v++; return *this; }
      bool operator!=(const value_ref& o) const {
        return v < o.v;
      }
    protected: 
      V v;
    };
    range_iterator(const V& _after)
      : after(_after)
    { }

    value_ref begin(void) const { return value_ref((V) 0); }
    value_ref end(void) const { return value_ref(after); }
  protected:     
    V after;
  };
 
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
      : g(_g), perm(_perm)
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

    typedef range_iterator<vert_id> vert_list;
    vert_list verts(void) const { return vert_list(perm.size()); }

    // GKG: Should probably modify this to handle cases where
    // the vertex iterator isn't just a vert_id*.
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
        return *this;
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
      adj_list(vector<vert_id>& _inv, const g_adj_list& _adj)
        : inv(_inv), adj(_adj)
      {
        vert_id* adj_end(adj.end());
        while(adj_end != adj.begin() && inv[adj_end[-1]] < 0)
          --adj_end;
      }
      adj_iterator begin(void) const { return adj_iterator(inv, adj.begin()); } 
      adj_iterator end(void) const { return adj_iterator(inv, adj_end); }

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

    GraphRev(G& _g) 
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

    typename G::vert_list verts(void) const { return g.verts(); }

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

    typedef pair<pair< vert_id, vert_id >, Wt> edge_ref;

    //===========================================
    // Enums used to mark vertices/edges during algorithms
    //===========================================
    // Edge colour during chromatic Dijkstra
    enum CMarkT { E_NONE = 0, E_LEFT = 1, E_RIGHT = 2, E_BOTH = 3 };
    // Whether a vertex is 'stable' during widening
    enum SMarkT { V_UNSTABLE = 0, V_STABLE = 1 };
    // Whether a vertex is in the current SCC/queue for Bellman-Ford.
    enum QMarkT { BF_NONE = 0, BF_SCC = 1, BF_QUEUED = 2 };
    // ===========================================
    // Scratch space needed by the graph algorithms.
    // Should really switch to some kind of arena allocator, rather
    // than having all these static structures.
    // ===========================================
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

      graph_t g;
      g.growTo(sz);
      
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

      graph_t g(graph_t::copy(l)); 
      
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

    template<class G1, class G2>
    static graph_t widen(G1& l, G2& r)
    {
      assert(l.size() == r.size());
      size_t sz = l.size();
      graph_t g;
      g.growTo(sz);
      /*
       * GKG: Check correctness
       */
      for(vert_id s : r.verts())
      {
        for(vert_id d : r.succs(s))
        {
          if(l.elem(s, d) && l.edge_val(s, d) <= r.edge_val(s, d))
            g.add_edge(s, l.edge_val(s, d), d);  
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
      // (gmp_z won't overflow, but there's a performance penalty
      // if magnitudes become large)
      //
      // Though this hurts our chances of early cutoff.
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
          vert_marks[v] = BF_SCC|BF_QUEUED;
          qtail++;
        }

        for(int iter = 0; iter < scc.size(); iter++)
        {
          for(; qtail != qhead; )
          {
            vert_id s = *(--qtail); 
            // If it _was_ on the queue, it must be in the SCC
            vert_marks[s] = BF_SCC; 
            
            Wt s_pot = potentials[s];

            for(vert_id d : g.succs(s))
            {
              Wt sd_pot = s_pot + g.edge_val(s, d);
              if(sd_pot < potentials[d])
              {
                potentials[d] = sd_pot;
                if(vert_marks[d] == BF_SCC)
                {
                  *next_tail = d;
                  vert_marks[d] = (BF_SCC|BF_QUEUED);
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
                vert_marks[v] = BF_NONE;
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
//        unsigned int g_count = 0;
//        unsigned int r_count = 0;
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
      vector< pair<vert_id, Wt> > adjs;
//      for(vert_id v = 0; v < sz; v++)
      for(vert_id v : g.verts())
      {
        adjs.clear();
        chrome_dijkstra(g, pots, colour_succs, v, adjs);

        for(pair<vert_id, Wt>& p : adjs)
          delta.push_back( make_pair(make_pair(v, p.first), p.second) );
      }

    }

    static void apply_delta(graph_t& g, edge_vector& delta)
    {
      for(pair< pair<vert_id, vert_id>, Wt>& e : delta)
      {
        g.set_edge(e.first.first, e.second, e.first.second);
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

      for(vert_id dest : g.succs(src))
      {
        dists[dest] = p[src] + g.edge_val(src, dest) - p[dest];
        dist_ts[dest] = ts;

        vert_marks[dest] = edge_marks[sz*src + dest];
        heap.insert(dest);
      }

      while(!heap.empty())
      {
        int es = heap.removeMin();
        Wt es_cost = dists[es] + p[es]; // If it's on the queue, distance is not infinite.
        Wt es_val = es_cost - p[src];
        if(!g.elem(src, es) || g.edge_val(src, es) > es_val)
          out.push_back( make_pair(es, es_val) );

        if(vert_marks[es] == (E_LEFT|E_RIGHT))
          continue;

        // Pick the appropriate set of successors
        vector<vert_id>& es_succs = (vert_marks[es] == E_LEFT) ?
          colour_succs[2*es+1] : colour_succs[2*es];
        for(vert_id ed : es_succs)
        {
          Wt v = es_cost + g.edge_val(es, ed) - p[ed];
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

    // Run Dijkstra's algorithm, but similar to the chromatic algorithm, avoid expanding
    // anything that _was_ stable.
    // GKG: Factor out common elements of this & the previous algorithm.
    template<class G, class P, class S>
    static void dijkstra_recover(G& g, P& p, const S& is_stable, vert_id src, vector< pair<vert_id, Wt> >& out)
    {
      unsigned int sz = g.size();
      if(sz == 0)
        return;
      if(is_stable[src])
        return;

      grow_scratch(sz);

      // Reset all vertices to infty.
      dist_ts[ts_idx] = ts++;
      ts_idx = (ts_idx+1) % dists.size();

      dists[src] = Wt(0);
      dist_ts[src] = ts;

      WtComp comp(dists);
      WtHeap heap(comp);

      for(vert_id dest : g.succs(src))
      {
        dists[dest] = p[src] + g.edge_val(src, dest) - p[dest];
        dist_ts[dest] = ts;

        vert_marks[dest] = V_UNSTABLE;
        heap.insert(dest);
      }

      while(!heap.empty())
      {
        int es = heap.removeMin();
        Wt es_cost = dists[es] + p[es]; // If it's on the queue, distance is not infinite.
        Wt es_val = es_cost - p[src];
        if(!g.elem(src, es) || g.edge_val(src, es) > es_val)
          out.push_back( make_pair(es, es_val) );

        if(vert_marks[es] == V_STABLE)
          continue;

        char es_mark = is_stable[es] ? V_STABLE : V_UNSTABLE;

        // Pick the appropriate set of successors
        for(vert_id ed : g.succs(es))
        {
          Wt v = es_cost + g.edge_val(es, ed) - p[ed];
          if(dist_ts[ed] != ts || v < dists[ed])
          {
            dists[ed] = v;
            dist_ts[ed] = ts;
            vert_marks[ed] = es_mark;

            if(heap.inHeap(ed))
            {
              heap.decrease(ed);
            } else {
              heap.insert(ed);
            }
          } else if(v == dists[ed]) {
            vert_marks[ed] |= es_mark;
          }
        }
      }
    }

    class forall_except {
    public:
      forall_except(vert_id _v)
        : v(_v)
      { }      

      bool operator[](vert_id w) const { return w != v; }
    protected:
      vert_id v;
    };

    /*
    template<class V>
    class vec_neg {
    public:
      vec_neg(const V& _vec)
        : vec(_vec)
      { }

      auto operator[](size_t idx) const { return -(vec[idx]); }
    protected:
      const V& vec;
    };
    template<class V>
    vec_neg<V> negate_vec(const V& vec) { return vec_neg<V>(vec); }
    */

    template<class G, class P>
    static void close_after_widen(graph_t& g, P& p, G& orig, edge_vector& delta)
    {
      unsigned int sz = g.size();
      assert(orig.size() == sz);
      
      for(vert_id v : g.verts())
      {
        // We're abusing edge_marks to store _vertex_ flags.
        // Should really just switch this to allocating regions of a fixed-size buffer.
        edge_marks[v] = V_UNSTABLE;
        // Assumption: stable iff |G(v)| = |H(v)|.
        if(g.succs(v).size() == orig.succs(v).size())
          edge_marks[v] = V_STABLE;
      }
      
      vector< pair<vert_id, Wt> > aux;
      for(vert_id v : g.verts())
      {
        if(!edge_marks[v])
        {
          aux.clear();
          dijkstra_recover(g, p, edge_marks, aux); 
          for(auto p : aux)
            delta.push_back( make_pair( make_pair(v, p.first), p.second ) );   
        }
      }
    }

    // Used for sorting successors of some vertex by increasing slack.
    // operator() may only be called on vertices for which
    // dists is initialized.
    template<class P>
    class AdjCmp {
    public: 
      AdjCmp(const P& _p)
        : p(_p)
      { }

      bool operator()(vert_id d1, vert_id d2) const {
       return (dists[d1] - p[d1]) < (dists[d2] - p[d2]);
      }
    protected:
      const P& p;
    };

    template<class P>
    static AdjCmp<P> make_adjcmp(const P& p)
    {
      return AdjCmp<P>(p);
    }

    // Compute the transitive closure of edges reachable from v, assuming
    // (1) the subgraph G \ {v} is closed, and (2) P is a valid model of G.
    template<class G, class P>
    static void close_after_assign_fwd(G& g, P& p, vert_id v, vector< pair<vert_id, Wt> >& aux)
    {
      // Initialize the queue and distances.
      vert_marks[v] = BF_QUEUED;
      dists[v] = Wt(0);
      vert_id* adj_head = dual_queue;
      vert_id* adj_tail = adj_head;
      for(vert_id d : g.succs(v))
      {
        vert_marks[d] = BF_QUEUED;
        dists[d] = g.edge_val(v, d);
        *adj_tail = d;
        adj_tail++;
      }

      // Sort the immediate edges by increasing slack.
      std::sort(adj_head, adj_tail, make_adjcmp(p));

      vert_id* reach_tail = adj_tail;
      for(; adj_head < adj_tail; adj_head++)
      {
        vert_id d = *adj_head;  
        
        Wt d_wt = dists[d];
        for(vert_id e : g.succs(d))
        {
          Wt e_wt = d_wt + g.edge_val(d, e);
          if(!vert_marks[e])
          {
            dists[e] = e_wt;
            vert_marks[e] = BF_QUEUED;
            *reach_tail = e;
            reach_tail++;
          } else {
            dists[e] = min(e_wt, dists[e]);
          }
        }
      }

      // Now collect the adjacencies, and clear vertex flags
      for(adj_head = dual_queue; adj_head < reach_tail; adj_head++)
      {
        aux.push_back(make_pair(*adj_head, dists[*adj_head]));
        vert_marks[*adj_head] = 0;
      }
    }
    
    template<class P>
    static void close_after_assign(graph_t& g, P& p, vert_id v, edge_vector& delta)
    {
      unsigned int sz = g.size();
      grow_scratch(sz);

      vector< pair<vert_id, Wt> > aux;

      close_after_assign_fwd(g, p, v, aux);
      for(auto p : aux)
        delta.push_back( make_pair( make_pair(v, p.first), p.second ) );   

      GraphRev<graph_t> g_rev(g);
      close_after_assign_fwd(g_rev, p, v, aux);
      for(auto p : aux)
        delta.push_back( make_pair( make_pair(p.first, v), p.second ) );
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
