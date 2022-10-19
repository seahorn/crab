#pragma once

#include <crab/domains/graphs/graph_iterators.hpp>
#include <crab/domains/patricia_trees.hpp>
#include <crab/support/os.hpp>
#include <crab/types/indexable.hpp>

#include <vector>

/*
 * Patricia-tree backed sparse weighted graph.
 * Trades some time penalty for much lower memory consumption.
 *
 * This is an implementation of a "persistent" graph, meaning that
 * after a graph update the new and old version share all the common
 * nodes and edges (i.e., structural sharing) between them. This is
 * what it makes this implementation very memory efficient. The price
 * to pay is that the lookups are done in logarithmic time.
 */
namespace crab {

template <class Weight> class PtGraph {
public:
  using Wt = Weight;
  using graph_t = PtGraph<Wt>;
  using vert_id = graph_iterators::vert_id;

  class wt_ref_t {
    const graph_t *g;
    vert_id s;
    vert_id d;
    // potentially expensive with patricia trees because we cannot
    // store a reference
    Wt w;

  public:
    wt_ref_t() : g(nullptr) {}

    wt_ref_t(const graph_t &_g, vert_id _s, vert_id _d)
        : g(&_g), s(_s), d(_d), w(g->succs(s).value(d)) {}

    const Wt &get() const {
      assert(g);
      return w;
    }
  };

private:
  class vert_idx : public indexable {
  public:
    vert_idx(vert_id _v) : v(_v) {}
    virtual ikos::index_t index(void) const override {
      return (ikos::index_t)v;
    }
    virtual void write(crab_os &o) const override { o << v; }

    vert_id v;
  };

  using pred_t = ikos::patricia_tree_set<vert_idx>;
  using succ_t = ikos::patricia_tree<vert_idx, Wt>;

  class pred_iterator {
  public:
    using ItP = typename pred_t::iterator;
    using iter_t = pred_iterator;

    pred_iterator(const ItP &_it) : it(_it) {}
    pred_iterator(void) : it() {}
    bool operator!=(const iter_t &o) { return it != o.it; }
    iter_t &operator++(void) {
      ++it;
      return *this;
    }
    vert_id operator*(void) const { return (*it).v; }

  private:
    ItP it;
  };

  class succ_iterator {
  public:
    using ItS = typename succ_t::iterator;
    using iter_t = succ_iterator;
    succ_iterator(const ItS &_it) : it(_it) {}
    succ_iterator(void) : it() {}
    // XXX: to make sure that we always return the same address
    // for the "empty" iterator, otherwise we can trigger
    // undefined behavior.
    static iter_t empty_iterator() {
      static std::unique_ptr<iter_t> it = nullptr;
      if (!it)
        it = std::unique_ptr<iter_t>(new iter_t());
      return *it;
    }
    bool operator!=(const iter_t &o) { return it != o.it; }
    iter_t &operator++(void) {
      ++it;
      return *this;
    }
    vert_id operator*(void) const { return (*it).first.v; }

  private:
    ItS it;
  };

  using edge_val_t = graph_iterators::edge_val_t<Wt>;

public:
  class pred_range {
  public:
    using iterator = pred_iterator;

    pred_range(pred_t &_p) : p(_p) {}
    iterator begin(void) const { return iterator(p.begin()); }
    iterator end(void) const { return iterator(p.end()); }
    size_t size(void) const { return p.size(); }
    bool mem(unsigned int v) const { return p[v]; }
    void add(unsigned int v) { p += v; }
    void remove(unsigned int v) { p -= v; }
    void clear() { p.clear(); }

  private:
    pred_t &p;
  };

  class const_pred_range {
  public:
    using iterator = pred_iterator;

    const_pred_range(const pred_t &_p) : p(_p) {}
    iterator begin(void) const { return iterator(p.begin()); }
    iterator end(void) const { return iterator(p.end()); }
    size_t size(void) const { return p.size(); }
    bool mem(unsigned int v) const { return p[v]; }

  private:
    const pred_t &p;
  };

  class succ_range {
  public:
    using iterator = succ_iterator;

    succ_range(succ_t &_p) : p(_p) {}
    iterator begin(void) const { return iterator(p.begin()); }
    iterator end(void) const { return iterator(p.end()); }
    size_t size(void) const { return p.size(); }
    bool mem(unsigned int v) const {
      if (p.lookup(v))
        return true;
      else
        return false;
    }
    void add(unsigned int v, const Wt &w) { p.insert(v, w); }
    Wt value(unsigned int v) const { return *(p.lookup(v)); }
    void remove(unsigned int v) { p.remove(v); }
    void clear() { p.clear(); }

  private:
    succ_t &p;
  };

  class const_succ_range {
  public:
    using iterator = succ_iterator;

    const_succ_range(const succ_t &_p) : p(_p) {}
    iterator begin(void) const { return iterator(p.begin()); }
    iterator end(void) const { return iterator(p.end()); }
    size_t size(void) const { return p.size(); }
    bool mem(unsigned int v) const {
      if (p.lookup(v))
        return true;
      else
        return false;
    }
    Wt value(unsigned int v) const { return *(p.lookup(v)); }

  private:
    const succ_t &p;
  };

  using vert_iterator = graph_iterators::vert_iterator;
  using vert_range = graph_iterators::vert_range;
  using const_e_succ_range =
      graph_iterators::const_fwd_edge_range<graph_t, succ_iterator, edge_val_t>;
  using const_e_pred_range =
      graph_iterators::const_rev_edge_range<graph_t, pred_iterator, edge_val_t>;

  PtGraph() : edge_count(0), _succs(), _preds(), is_free(), free_id() {}

  // PtGraph is a persistent graph so the copy is cheap.
  PtGraph(const PtGraph<Wt> &o)
      : edge_count(o.edge_count), _succs(o._succs), _preds(o._preds),
        is_free(o.is_free), free_id(o.free_id) {}

  PtGraph(PtGraph<Wt> &&o)
      : edge_count(o.edge_count), _succs(std::move(o._succs)),
        _preds(std::move(o._preds)), is_free(std::move(o.is_free)),
        free_id(std::move(o.free_id)) {
    o.edge_count = 0;
  }

  // PtGraph is a persistent graph so the assignment is cheap  
  PtGraph &operator=(const PtGraph<Wt> &o) {
    if (this != &o) {
      edge_count = o.edge_count;
      _succs = o._succs;
      _preds = o._preds;
      is_free = o.is_free;
      free_id = o.free_id;
    }
    return *this;
  }

  PtGraph &operator=(PtGraph<Wt> &&o) {
    if (this != &o) {
      edge_count = o.edge_count;
      _succs = std::move(o._succs);
      _preds = std::move(o._preds);
      is_free = std::move(o.is_free);
      free_id = std::move(o.free_id);
    }

    return *this;
  }

  ~PtGraph() {}

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
      _succs.push_back(succ_t());
      _preds.push_back(pred_t());
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
    if (!succs(x).mem(y)) {
      return false;
    }
    w = wt_ref_t(*this, x, y);
    return true;
  }

  // Important: no longer a ref
  Wt edge_val(vert_id x, vert_id y) const { return succs(x).value(y); }

  // Precondition: elem(x, y) is true.
  Wt operator()(vert_id x, vert_id y) { return succs(x).value(y); }

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
    else
      _succs[s].insert(vert_idx(d), w);
  }

  void set_edge_if_less_than(vert_id s, Wt w, vert_id d) {
    // assert(s < size() && d < size());
    if (!elem(s, d)) {
      add_edge(s, w, d);
    } else {
      // Expensive with patricia trees because it requires
      // lookup+insert
      // 
      // TODO: extend patricia tree API to insert a new pair key-value
      // if some user-given binary predicate holds on the old and new
      // value associated to the key
      if (w < edge_val(s, d)) {
        _succs[s].insert(vert_idx(d), w);
      }
    }
  }

  template <class Op> void update_edge(vert_id s, Wt w, vert_id d, Op &op) {
    if (elem(s, d)) {
      // _succs[s].insert(vert_idx(d), w, op);
      _succs[s].insert(vert_idx(d), op.apply(edge_val(s, d), w));
      return;
    }

    if (!op.default_is_absorbing())
      add_edge(s, w, d);
  }

  // FIXME: Verts currently iterates over free vertices,
  // as well as existing ones
  vert_range verts(void) const { return vert_range(is_free.size(), is_free); }

  succ_range succs(vert_id v) { return succ_range(_succs[v]); }

  const_succ_range succs(vert_id v) const {
    return const_succ_range(_succs[v]);
  }

  pred_range preds(vert_id v) { return pred_range(_preds[v]); }

  const_pred_range preds(vert_id v) const {
    return const_pred_range(_preds[v]);
  }

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
      _preds.push_back(pred_t());
      _succs.push_back(succ_t());
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

  friend crab::crab_os &operator<<(crab::crab_os &o, const PtGraph<Weight> &g) {
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
