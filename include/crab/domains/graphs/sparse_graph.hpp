#pragma once

#include <crab/domains/graphs/graph_iterators.hpp>
#include <crab/domains/graphs/util/reference_wrapper.hpp>
#include <crab/support/os.hpp>
#include <memory>
#include <vector>

/*
 * An implementation of a sparse weighted graph.
 */
namespace crab {

template <class Weight> class SparseWtGraph {
public:
  using Wt = Weight;
  using graph_t = SparseWtGraph<Wt>;
  using vert_id = graph_iterators::vert_id;
  using wt_ref_t = crab::reference_wrapper<const Wt>;

private:
  class adj_iterator {
  public:
    adj_iterator(void) : ptr(nullptr) {}
    adj_iterator(uint16_t *_p) : ptr(_p) {}
    // XXX: to make sure that we always return the same address
    // for the "empty" iterator, otherwise we can trigger
    // undefined behavior.
    static adj_iterator empty_iterator() {
      static std::unique_ptr<adj_iterator> it = nullptr;
      if (!it)
        it = std::unique_ptr<adj_iterator>(new adj_iterator());
      return *it;
    }
    vert_id operator*(void) const { return (vert_id)*ptr; }
    adj_iterator &operator++(void) {
      ptr++;
      return *this;
    }
    bool operator!=(const adj_iterator &o) const { return ptr < o.ptr; }

  protected:
    uint16_t *ptr;
  };

  class adj_list {
  public:
    using iterator = adj_iterator;

    adj_list(uint16_t *_ptr, unsigned int max_sz)
        : ptr(_ptr), sparseptr(_ptr + 1 + max_sz) {}
    adj_iterator begin(void) const { return adj_iterator(ptr + 1); }
    adj_iterator end(void) const { return adj_iterator(ptr + 1 + size()); }
    vert_id operator[](unsigned int idx) const { return (vert_id)ptr[1 + idx]; }
    uint16_t *sparse(void) const { return sparseptr; }
    uint16_t *dense(void) const { return (uint16_t *)(ptr + 1); }
    unsigned int size(void) const { return *ptr; }

    bool mem(unsigned int v) const {
      unsigned int idx = sparse()[v];
      return idx < size() && dense()[idx] == v;
    }
    void add(unsigned int v) {
      //       assert(!mem(v));
      //        assert(dense()+v < sparse());
      unsigned int idx = *ptr;
      dense()[idx] = v;
      sparse()[v] = idx;
      (*ptr)++;
    }
    void remove(unsigned int v) {
      //        assert(mem(v));
      (*ptr)--;
      unsigned int idx = sparse()[v];
      unsigned int w = dense()[*ptr];
      dense()[idx] = w;
      sparse()[w] = idx;
    }
    void clear() { *ptr = 0; }

  protected:
    uint16_t *ptr;
    uint16_t *sparseptr;
  };

  using edge_ref_t = graph_iterators::edge_ref_t<Wt>;
  using const_edge_ref_t = graph_iterators::const_edge_ref_t<Wt>;

public:
  using vert_iterator = graph_iterators::vert_iterator;
  using vert_range = graph_iterators::vert_range;
  using const_succ_range = adj_list;
  using const_pred_range = adj_list;
  using succ_range = const_succ_range;
  using pred_range = const_pred_range;
  using e_succ_range =
      graph_iterators::fwd_edge_range<graph_t, adj_iterator, edge_ref_t>;
  using e_pred_range =
      graph_iterators::rev_edge_range<graph_t, adj_iterator, edge_ref_t>;
  using const_e_succ_range =
      graph_iterators::const_fwd_edge_range<graph_t, adj_iterator,
                                            const_edge_ref_t>;
  using const_e_pred_range =
      graph_iterators::const_rev_edge_range<graph_t, adj_iterator,
                                            const_edge_ref_t>;

  SparseWtGraph(unsigned int _maxsz = 10, float _growth_rate = 1.4)
      : max_sz(_maxsz), sz(0), growth_rate(_growth_rate), edge_count(0),
        fwd_adjs(
            (uint16_t *)malloc(sizeof(uint16_t) * max_sz * (2 * max_sz + 1))),
        rev_adjs(
            (uint16_t *)malloc(sizeof(uint16_t) * max_sz * (2 * max_sz + 1))),
        mtx((Wt *)malloc(sizeof(Wt) * max_sz * max_sz)) {
    /*
      for(vert_id v = 0 ; v < sz; v++) {
        succs(v).clear(); preds(v).clear();
      }
      check_adjs();
    */
  }

  template <class Wo>
  SparseWtGraph(const SparseWtGraph<Wo> &o)
      : max_sz(o.max_sz), sz(o.sz), growth_rate(o.growth_rate), edge_count(0),
        fwd_adjs(
            (uint16_t *)malloc(sizeof(uint16_t) * max_sz * (2 * max_sz + 1))),
        rev_adjs(
            (uint16_t *)malloc(sizeof(uint16_t) * max_sz * (2 * max_sz + 1))),
        mtx((Wt *)malloc(sizeof(Wt) * max_sz * max_sz)), is_free(o.is_free),
        free_id(o.free_id) {
    assert(sz <= max_sz);
    for (vert_id v = 0; v < sz; v++) {
      succs(v).clear();
      preds(v).clear();
    }

    for (vert_id v : o.verts())
      for (vert_id d : o.succs(v))
        add_edge(v, o.edge_val(v, d), d);

    //      check_adjs();
  }

  SparseWtGraph(const SparseWtGraph<Wt> &o)
      : max_sz(o.max_sz), sz(o.sz), growth_rate(o.growth_rate), edge_count(0),
        fwd_adjs(
            (uint16_t *)malloc(sizeof(uint16_t) * max_sz * (2 * max_sz + 1))),
        rev_adjs(
            (uint16_t *)malloc(sizeof(uint16_t) * max_sz * (2 * max_sz + 1))),
        mtx((Wt *)malloc(sizeof(Wt) * max_sz * max_sz)), is_free(o.is_free),
        free_id(o.free_id) {
    assert(sz <= max_sz);
    for (vert_id v = 0; v < sz; v++) {
      succs(v).clear();
      preds(v).clear();
    }

    for (vert_id v : o.verts())
      for (vert_id d : o.succs(v))
        add_edge(v, o.edge_val(v, d), d);

    //      check_adjs();
  }

  SparseWtGraph(SparseWtGraph<Wt> &&o)
      : max_sz(o.max_sz), sz(o.sz), growth_rate(o.growth_rate),
        edge_count(o.edge_count), fwd_adjs(o.fwd_adjs), rev_adjs(o.rev_adjs),
        mtx(o.mtx), is_free(std::move(o.is_free)),
        free_id(std::move(o.free_id)) {
    o.max_sz = 0;
    o.sz = 0;
    o.fwd_adjs = NULL;
    o.rev_adjs = NULL;
    o.mtx = NULL;
  }

  SparseWtGraph &operator=(const SparseWtGraph<Wt> &o) {
    if ((&o) == this)
      return *this;

    // Make sure the number of vertices matches,
    // and active adjacency lists are initialized.
    /*
    clear_edges();
    growCap(o.sz);
    for(; sz < o.sz; sz++)
    {
      // Initialize extra adj_lists
      succs(sz).clear();
      preds(sz).clear();
    }
    sz = o.sz;
    */
    clear();
    growTo(o.sz);

    // Now copy the edges
    for (vert_id s : o.verts()) {
      for (vert_id d : o.succs(s)) {
        add_edge(s, o.edge_val(s, d), d);
      }
    }

    // Copy the free vertices
    free_id = o.free_id;
    is_free = o.is_free;

    //      check_adjs();

    return *this;
  }

  SparseWtGraph &operator=(SparseWtGraph<Wt> &&o) {
    if (max_sz > 0) {
      free(mtx);
      free(fwd_adjs);
      free(rev_adjs);
    }

    max_sz = o.max_sz;
    sz = o.sz;
    growth_rate = o.growth_rate;
    edge_count = o.edge_count;
    fwd_adjs = o.fwd_adjs;
    rev_adjs = o.rev_adjs;
    mtx = o.mtx;
    is_free = std::move(o.is_free);
    free_id = std::move(o.free_id);

    o.max_sz = 0;
    o.sz = 0;
    o.fwd_adjs = NULL;
    o.rev_adjs = NULL;
    o.mtx = NULL;

    return *this;
  }

  // void check_adjs(void) {
  //   for(vert_id v : verts()) {
  //     assert(succs(v).size() <= sz);
  //     for(vert_id s : succs(v)) {
  //       assert(s < sz);
  //       assert(preds(s).mem(v));
  //     }

  //     assert(preds(v).size() <= sz);
  //     for(vert_id p : preds(v)) {
  //       assert(p < sz);
  //       assert(succs(p).mem(v));
  //     }
  //   }
  // }

  ~SparseWtGraph() {
    if (max_sz > 0) {
      for (vert_id v = 0; v < sz; v++) {
        for (vert_id d : succs(v))
          (&(mtx[max_sz * v + d]))->~Wt();
      }
      free(mtx);
      free(fwd_adjs);
      free(rev_adjs);
    }
  }

  template <class G> static graph_t copy(G &g) {
    graph_t ret;
    ret.growTo(g.size());

    for (vert_id s : g.verts()) {
      for (vert_id d : g.succs(s)) {
        ret.add_edge(s, g.edge_val(s, d), d);
      }
    }
    // Unforunately, views don't have
    // convenient 'is_free' arrays.
    //      ret.free_id = g.free_id;
    //      ret.is_free = g.is_free;
    //      ret.check_adjs();

    return ret;
  }

  bool is_empty(void) const { return edge_count == 0; }

  vert_id new_vertex(void) {
    vert_id v;
    if (free_id.size() > 0) {
      v = free_id.back();
      assert(v < sz);
      free_id.pop_back();
      is_free[v] = false;
    } else {
      if (max_sz <= sz) {
        // Make sure max_sz is strictly increasing
        growCap(1 + max_sz * growth_rate);
      }
      assert(sz < max_sz);
      v = sz++;
      is_free.push_back(false);
    }
    preds(v).clear();
    succs(v).clear();

    return v;
  }

  void forget(vert_id v) {
    assert(v < sz);
    // FIXME: currently assumes v is not free
    if (is_free[v])
      return;

    free_id.push_back(v);
    is_free[v] = true;

    // Remove (s -> v) from preds.
    edge_count -= succs(v).size();
    for (vert_id s : succs(v)) {
      (&(mtx[max_sz * v + s]))->~Wt();
      preds(s).remove(v);
    }
    succs(v).clear();

    // Remove (v -> p) from succs
    edge_count -= preds(v).size();
    for (vert_id p : preds(v)) {
      (&(mtx[max_sz * p + v]))->~Wt();
      succs(p).remove(v);
    }
    preds(v).clear();
  }

  // Check whether an edge is live
  bool elem(vert_id x, vert_id y) const { return succs(x).mem(y); }

  bool lookup(vert_id x, vert_id y, wt_ref_t &w) const {
    if (!succs(x).mem(y))
      return false;
    w = wt_ref_t(mtx[max_sz * x + y]);
    return true;
  }

  Wt &edge_val(vert_id x, vert_id y) { return mtx[max_sz * x + y]; }

  const Wt &edge_val(vert_id x, vert_id y) const { return mtx[max_sz * x + y]; }

  // Precondition: elem(x, y) is true.
  Wt operator()(vert_id x, vert_id y) const { return mtx[max_sz * x + y]; }

  void clear_edges(void) {
    for (vert_id v : verts()) {
      // Make sure the matrix entries are destructed
      for (vert_id d : succs(v))
        (&(mtx[max_sz * v + d]))->~Wt();

      preds(v).clear();
      succs(v).clear();
    }
    edge_count = 0;
  }

  void clear(void) {
    clear_edges();
    is_free.clear();
    free_id.clear();
    sz = 0;
  }

  // Number of allocated vertices
  int size(void) const { return sz; }

  // Number of edges
  size_t num_edges(void) const { return edge_count; }

  // Assumption: (x, y) not in mtx
  void add_edge(vert_id x, Wt wt, vert_id y) {
    // assert(x < size() && y < size());
    // assert(x != y);
    assert(!elem(x, y));
    succs(x).add(y);
    preds(y).add(x);
    // Allocate a new Wt at the given offset.
    new (&(mtx[x * max_sz + y])) Wt(wt);
    edge_count++;
  }

  void set_edge(vert_id s, Wt w, vert_id d) {
    //  assert(s < size() && d < size());
    if (!elem(s, d))
      add_edge(s, w, d);
    else
      edge_val(s, d) = w;
  }

  void set_edge_if_less_than(vert_id s, Wt w, vert_id d) {
    if (!elem(s, d)) {
      add_edge(s, w, d);
    } else {
      Wt &v = edge_val(s, d);
      if (w < v) {
        v = w;
      }
    }
  }

  template <class Op> void update_edge(vert_id s, Wt w, vert_id d, Op &op) {
    if (elem(s, d)) {
      edge_val(s, d) = op.apply(edge_val(s, d), w);
      return;
    }

    if (!op.default_is_absorbing())
      add_edge(s, w, d);
  }

  // FIXME: Verts currently iterates over free vertices,
  // as well as existing ones
  vert_range verts(void) const { return vert_range(sz, is_free); }
  succ_range succs(vert_id v) {
    return succ_range(fwd_adjs + v * (2 * max_sz + 1), max_sz);
  }
  pred_range preds(vert_id v) {
    return pred_range(rev_adjs + v * (2 * max_sz + 1), max_sz);
  }
  const_succ_range succs(vert_id v) const {
    return const_succ_range(fwd_adjs + v * (2 * max_sz + 1), max_sz);
  }
  const_pred_range preds(vert_id v) const {
    return const_pred_range(rev_adjs + v * (2 * max_sz + 1), max_sz);
  }
  e_succ_range e_succs(vert_id v) { return e_succ_range(*this, v); }
  e_pred_range e_preds(vert_id v) { return e_pred_range(*this, v); }
  const_e_succ_range e_succs(vert_id v) const {
    return const_e_succ_range(*this, v);
  }
  const_e_pred_range e_preds(vert_id v) const {
    return const_e_pred_range(*this, v);
  }

  // growTo shouldn't be used after forget
  void growTo(unsigned int new_sz) {
    growCap(new_sz);
    for (; sz < new_sz; sz++) {
      succs(sz).clear();
      preds(sz).clear();
      is_free.push_back(false);
    }
  }

  void write(crab_os &o) const {
    o << "[|";
    bool first = true;
    for (vert_id v = 0; v < sz; v++) {
      auto it = succs(v).begin();
      auto end = succs(v).end();

      if (it != end) {
        if (first)
          first = false;
        else
          o << ", ";

        o << "[v" << v << " -> ";
        o << "(" << edge_val(v, *it) << ":" << *it << ")";
        for (++it; it != end; ++it) {
          o << ", (" << edge_val(v, *it) << ":" << *it << ")";
        }
        o << "]";
      }
    }
    o << "|]";
  }

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const SparseWtGraph<Weight> &g) {
    g.write(o);
    return o;
  }

protected:
  // Allocate new memory, and duplicate
  // the content.
  // Add an element
  void _adj_add(uint16_t *adj, unsigned int max, unsigned int val) {
    adj[1 + *adj] = val;
    adj[1 + max + val] = *adj;
    (*adj)++;
  }

  void growCap(unsigned int new_max) {
    //      check_adjs();
    if (new_max <= max_sz)
      return;

    unsigned int new_mtxsz = new_max * new_max;
    unsigned int new_adjsz = (2 * new_max) + 1;

    Wt *new_mtx = (Wt *)malloc(sizeof(Wt) * new_mtxsz);
    uint16_t *new_fwd =
        (uint16_t *)malloc(sizeof(uint16_t) * new_max * new_adjsz);
    assert(new_fwd);
    uint16_t *new_rev =
        (uint16_t *)malloc(sizeof(uint16_t) * new_max * new_adjsz);
    assert(new_rev);

    for (vert_id v = 0; v < sz; v++) {
      if (is_free[v])
        continue;
      assert(v < new_max);

      uint16_t *new_fwd_ptr = new_fwd + v * new_adjsz;
      *new_fwd_ptr = 0;
      for (vert_id d : succs(v)) {
        assert(d < new_max);
        _adj_add(new_fwd_ptr, new_max, d);
        //          new_mtx[v*new_max + d] = mtx[v*max_sz + d];
        new (&(new_mtx[v * new_max + d])) Wt(mtx[v * max_sz + d]);
        (&(mtx[max_sz * v + d]))->~Wt();
      }

      uint16_t *new_rev_ptr = new_rev + v * new_adjsz;
      *new_rev_ptr = 0;
      for (vert_id s : preds(v)) {
        assert(s < new_max);
        _adj_add(new_rev_ptr, new_max, s);
      }
    }
    free(mtx);
    free(fwd_adjs);
    free(rev_adjs);

    mtx = new_mtx;
    fwd_adjs = new_fwd;
    rev_adjs = new_rev;

    max_sz = new_max;

    //      check_adjs();
  }

  unsigned int max_sz;
  unsigned int sz;
  float growth_rate;
  unsigned int edge_count;

  // Each element of fwd/rev adjs:
  // [ sz/1 | dense/max_sz | inv/max_sz ]
  // Total size: sizeof(uint) * max_sz * (1 + 2*max_sz)
  uint16_t *fwd_adjs;
  uint16_t *rev_adjs;
  Wt *mtx;

  std::vector<bool> is_free;
  std::vector<int> free_id;
};

} // namespace crab
