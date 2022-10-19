#pragma once

#include <crab/domains/graphs/graph_iterators.hpp>
#include <crab/domains/graphs/util/reference_wrapper.hpp>
#include <crab/support/os.hpp>

#include <unordered_map>
#include <unordered_set>
#include <vector>

/*
 * Hash table backed implementation of a sparse, weighted graph.
 */
namespace crab {

template <class Weight> class HtGraph {
  enum { table_init_sz = 17 };

public:
  using Wt = Weight;
  using graph_t = HtGraph<Wt>;
  using vert_id = graph_iterators::vert_id;
  using wt_ref_t = crab::reference_wrapper<const Wt>;

private:
  using pred_t = std::unordered_set<vert_id>;
  using succ_t = std::unordered_map<vert_id, Wt>;

  class pred_iterator {
  public:
    using ItP = typename pred_t::const_iterator;
    using type = pred_iterator;

    pred_iterator(ItP _it) : it(_it) {}
    pred_iterator(void) : it() {}

    bool operator!=(const type &o) const { return it != o.it; }
    type &operator++(void) {
      ++it;
      return *this;
    }
    vert_id operator*(void) const { return (*it); }

  protected:
    ItP it;
  };

  class succ_iterator {
  public:
    using ItS = typename succ_t::const_iterator;
    using type = succ_iterator;

    succ_iterator(ItS _it) : it(_it) {}
    succ_iterator(void) : it() {}

    // XXX: to make sure that we always return the same address
    // for the "empty" iterator, otherwise we can trigger
    // undefined behavior.
    static type empty_iterator() {
      static std::unique_ptr<type> it = nullptr;
      if (!it)
        it = std::unique_ptr<type>(new type());
      return *it;
    }
    bool operator!=(const type &o) const { return it != o.it; }
    type &operator++(void) {
      ++it;
      return *this;
    }
    vert_id operator*(void) const { return (*it).first; }

  protected:
    ItS it;
  };

  using edge_ref_t = graph_iterators::edge_ref_t<Wt>;
  using const_edge_ref_t = graph_iterators::const_edge_ref_t<Wt>;

public:
  using vert_iterator = graph_iterators::vert_iterator;
  using vert_range = graph_iterators::vert_range;
  using pred_range = graph_iterators::pred_range<pred_t, pred_iterator>;
  using const_pred_range =
      graph_iterators::const_pred_range<pred_t, pred_iterator>;
  using succ_range = graph_iterators::succ_range<succ_t, succ_iterator>;
  using const_succ_range =
      graph_iterators::const_succ_range<succ_t, succ_iterator>;
  using e_succ_range =
      graph_iterators::fwd_edge_range<graph_t, succ_iterator, edge_ref_t>;
  using e_pred_range =
      graph_iterators::rev_edge_range<graph_t, pred_iterator, edge_ref_t>;
  using const_e_succ_range =
      graph_iterators::const_fwd_edge_range<graph_t, succ_iterator,
                                            const_edge_ref_t>;
  using const_e_pred_range =
      graph_iterators::const_rev_edge_range<graph_t, pred_iterator,
                                            const_edge_ref_t>;

  HtGraph() : edge_count(0), _succs(), _preds(), is_free(), free_id() {}

  template <class Wo> HtGraph(const HtGraph<Wo> &o) : edge_count(0) {
    for (vert_id v : o.verts())
      for (vert_id d : o.succs(v))
        add_edge(v, o.edge_val(v, d), d);
  }

  HtGraph(const HtGraph<Wt> &o)
      : edge_count(o.edge_count), _succs(o._succs), _preds(o._preds),
        is_free(o.is_free), free_id(o.free_id) {}

  HtGraph(HtGraph<Wt> &&o)
      : edge_count(o.edge_count), _succs(std::move(o._succs)),
        _preds(std::move(o._preds)), is_free(std::move(o.is_free)),
        free_id(std::move(o.free_id)) {
    o.edge_count = 0;
  }

  HtGraph &operator=(const HtGraph<Wt> &o) {
    if ((&o) == this)
      return *this;

    edge_count = o.edge_count;
    _succs = o._succs;
    _preds = o._preds;
    is_free = o.is_free;
    free_id = o.free_id;

    return *this;
  }

  HtGraph &operator=(HtGraph<Wt> &&o) {
    edge_count = o.edge_count;
    _succs = std::move(o._succs);
    _preds = std::move(o._preds);
    is_free = std::move(o.is_free);
    free_id = std::move(o.free_id);

    return *this;
  }

  // void check_adjs(void) {
  //   for(vert_id v : verts()) {
  //     assert(succs(v).size() <= _succs.size());
  //     for(vert_id s : succs(v)) {
  //       assert(s < _succs.size());
  //       assert(preds(s).mem(v));
  //     }
  //     assert(preds(v).size() <= _succs.size());
  //     for(vert_id p : preds(v)) {
  //       assert(p < _succs.size());
  //       assert(succs(p).mem(v));
  //     }
  //   }
  // }

  ~HtGraph() {}

  // GKG: Can do this more efficiently
  template <class G> static graph_t copy(G &g) {
    graph_t ret;
    ret.growTo(g.size());

    for (vert_id s : g.verts()) {
      for (vert_id d : g.succs(s)) {
        ret.add_edge(s, g.edge_val(s, d), d);
      }
    }

    return ret;
  }

  bool is_empty(void) const { return edge_count == 0; }

  vert_id new_vertex(void) {
    vert_id v;
    if (free_id.size() > 0) {
      v = free_id.back();
      assert(v < _succs.size());
      free_id.pop_back();
      is_free[v] = false;
    } else {
      v = is_free.size();
      _succs.push_back(succ_t(table_init_sz));
      _preds.push_back(pred_t(table_init_sz));
      is_free.push_back(false);
    }

    return v;
  }

  void forget(vert_id v) {
    assert(v < _succs.size());
    if (is_free[v])
      return;

    free_id.push_back(v);
    is_free[v] = true;

    // Remove (s -> v) from preds.
    edge_count -= succs(v).size();
    for (vert_id d : succs(v))
      preds(d).remove(v);
    _succs[v].clear();

    // Remove (v -> p) from succs
    edge_count -= preds(v).size();
    for (vert_id p : preds(v))
      succs(p).remove(v);
    _preds[v].clear();
  }

  // Check whether an edge is live
  bool elem(vert_id x, vert_id y) const { return succs(x).mem(y); }

  bool lookup(vert_id x, vert_id y, wt_ref_t &w) const {
    if (!succs(x).mem(y))
      return false;
    w = succs(x).value(y);
    return true;
  }

  Wt &edge_val(vert_id x, vert_id y) { return succs(x).value(y); }

  const Wt &edge_val(vert_id x, vert_id y) const { return succs(x).value(y); }

  // Precondition: elem(x, y) is true.
  Wt operator()(vert_id x, vert_id y) {
    // return mtx[max_sz*x + y];
    return succs(x).value(y);
  }

  void clear_edges(void) {
    edge_count = 0;
    for (vert_id v : verts()) {
      _succs[v].clear();
      _preds[v].clear();
    }
  }

  void clear(void) {
    edge_count = 0;
    is_free.clear();
    free_id.clear();
    _succs.clear();
    _preds.clear();
  }

  // Number of allocated vertices
  int size(void) const { return is_free.size(); }

  // Number of edges
  size_t num_edges(void) const { return edge_count; }

  // Assumption: (x, y) not in mtx

  void add_edge(vert_id x, Wt wt, vert_id y) {
    assert(!elem(x, y));
    succs(x).add(y, wt);
    preds(y).add(x);
    edge_count++;
  }

  void set_edge(vert_id s, Wt w, vert_id d) {
    // assert(s < size() && d < size());
    if (!elem(s, d))
      add_edge(s, w, d);
    else {
      edge_val(s, d) = w;
    }
  }

  void set_edge_if_less_than(vert_id s, Wt w, vert_id d) {
    // assert(s < size() && d < size());
    if (!elem(s, d))
      add_edge(s, w, d);
    else {
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
  vert_range verts(void) const { return vert_range(is_free.size(), is_free); }
  succ_range succs(vert_id v) { return succ_range(_succs[v]); }
  pred_range preds(vert_id v) { return pred_range(_preds[v]); }
  const_succ_range succs(vert_id v) const {
    return const_succ_range(_succs[v]);
  }
  const_pred_range preds(vert_id v) const {
    return const_pred_range(_preds[v]);
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
    size_t sz = is_free.size();
    for (; sz < new_sz; sz++) {
      is_free.push_back(false);
      _preds.push_back(pred_t(table_init_sz));
      _succs.push_back(succ_t(table_init_sz));
    }
    assert(free_id.size() == 0);
  }

  void write(crab_os &o) const {
    o << "[|";
    bool first = true;
    for (vert_id v = 0; v < _succs.size(); v++) {
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

  friend crab::crab_os &operator<<(crab::crab_os &o, const HtGraph<Weight> &g) {
    g.write(o);
    return o;
  }

private:
  unsigned int edge_count;
  std::vector<succ_t> _succs;
  std::vector<pred_t> _preds;
  std::vector<bool> is_free;
  std::vector<int> free_id;
};

} // namespace crab
