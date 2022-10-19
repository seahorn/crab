#pragma once

#include <crab/domains/graphs/adapt_smap.hpp>
#include <crab/domains/graphs/graph_iterators.hpp>
#include <crab/domains/graphs/util/Vec.h>
#include <crab/domains/graphs/util/reference_wrapper.hpp>
#include <crab/support/os.hpp>

// Adaptive sparse-set based weighted graph implementation

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {

template <class Weight> class AdaptGraph {
  using smap_t = AdaptSMap<size_t>;

public:
  using vert_id = graph_iterators::vert_id;
  using Wt = Weight;
  using wt_ref_t = crab::reference_wrapper<const Wt>;

private:
  vec<smap_t> _preds;
  vec<smap_t> _succs;
  vec<Wt> _ws;
  int edge_count;
  std::vector<bool> is_free;
  vec<vert_id> free_id;
  vec<size_t> free_widx;

  class edge_iter {
  public:
    using edge_wrapper = typename graph_iterators::edge_ref_t<Wt>;
    edge_iter(const smap_t::elt_iter_t &_it, vec<Wt> &_ws)
        : it(_it), ws(&_ws) {}
    edge_iter(const edge_iter &o) : it(o.it), ws(o.ws) {}
    edge_iter(void) : ws(nullptr) {}

    // XXX: to make sure that we always return the same address
    // for the "empty" iterator, otherwise we can trigger
    // undefined behavior.
    static edge_iter empty_iterator() {
      static std::unique_ptr<edge_iter> it = nullptr;
      if (!it)
        it = std::unique_ptr<edge_iter>(new edge_iter());
      return *it;
    }

    edge_wrapper operator*(void) const {
      return edge_wrapper((*it).key, (*ws)[(*it).val]);
    }
    edge_iter operator++(void) {
      ++it;
      return *this;
    }
    bool operator!=(const edge_iter &o) const { return it != o.it; }

    smap_t::elt_iter_t it;
    vec<Wt> *ws;
  };

  class const_edge_iter {
  public:
    using edge_wrapper = typename graph_iterators::const_edge_ref_t<Wt>;
    const_edge_iter(const smap_t::elt_iter_t &_it, const vec<Wt> &_ws)
        : it(_it), ws(&_ws) {}
    const_edge_iter(const edge_iter &o) : it(o.it), ws(o.ws) {}
    const_edge_iter(void) : ws(nullptr) {}

    static const_edge_iter empty_iterator() {
      static std::unique_ptr<const_edge_iter> it = nullptr;
      if (!it)
        it = std::unique_ptr<const_edge_iter>(new const_edge_iter());
      return *it;
    }

    edge_wrapper operator*(void) const {
      return edge_wrapper((*it).key, (*ws)[(*it).val]);
    }

    const_edge_iter operator++(void) {
      ++it;
      return *this;
    }
    bool operator!=(const const_edge_iter &o) const { return it != o.it; }

    smap_t::elt_iter_t it;
    const vec<Wt> *ws;
  };

  class edge_range_t {
  public:
    using elt_range_t = typename smap_t::elt_range_t;
    using iterator = edge_iter;
    edge_range_t(const edge_range_t &o) : r(o.r), ws(o.ws) {}
    edge_range_t(const elt_range_t &_r, vec<Wt> &_ws) : r(_r), ws(_ws) {}
    iterator begin(void) { return iterator(r.begin(), ws); }
    iterator end(void) { return iterator(r.end(), ws); }
    size_t size(void) const { return r.size(); }

    elt_range_t r;
    vec<Wt> &ws;
  };

  class const_edge_range_t {
  public:
    using elt_range_t = typename smap_t::elt_range_t;
    using iterator = const_edge_iter;
    const_edge_range_t(const const_edge_range_t &o) : r(o.r), ws(o.ws) {}
    const_edge_range_t(const elt_range_t &_r, const vec<Wt> &_ws)
        : r(_r), ws(_ws) {}
    iterator begin(void) const { return iterator(r.begin(), ws); }
    iterator end(void) const { return iterator(r.end(), ws); }
    size_t size(void) const { return r.size(); }

    elt_range_t r;
    const vec<Wt> &ws;
  };

public:
  using vert_iterator = graph_iterators::vert_iterator;
  using vert_range = graph_iterators::vert_range;
  using const_pred_range = typename smap_t::key_range_t;
  using const_succ_range = typename smap_t::key_range_t;
  using pred_range = const_pred_range;
  using succ_range = const_succ_range;
  using e_pred_range = edge_range_t;
  using e_succ_range = edge_range_t;
  using const_e_pred_range = const_edge_range_t;
  using const_e_succ_range = const_edge_range_t;

  AdaptGraph(void) : edge_count(0) {}
  AdaptGraph(AdaptGraph<Wt> &&o) = default;
  AdaptGraph(const AdaptGraph<Wt> &o) = default;
  AdaptGraph<Wt> &operator=(const AdaptGraph<Wt> &o) = default;
  AdaptGraph<Wt> &operator=(AdaptGraph<Wt> &&o) = default;

  template <class G> static AdaptGraph<Wt> copy(const G &o) {
    AdaptGraph<Wt> g;
    g.growTo(o.size());

    for (vert_id s : o.verts()) {
      for (auto e : const_cast<G &>(o).e_succs(s)) {
        g.add_edge(s, e.val, e.vert);
      }
    }
    return g;
  }

  vert_range verts(void) const { return vert_range(is_free.size(), is_free); }

  succ_range succs(vert_id v) { return _succs[v].keys(); }

  pred_range preds(vert_id v) { return _preds[v].keys(); }

  const_succ_range succs(vert_id v) const { return _succs[v].keys(); }

  const_pred_range preds(vert_id v) const { return _preds[v].keys(); }

  e_succ_range e_succs(vert_id v) {
    return e_succ_range(_succs[v].elts(), _ws);
  }

  e_pred_range e_preds(vert_id v) {
    return e_pred_range(_preds[v].elts(), _ws);
  }

  const_e_succ_range e_succs(vert_id v) const {
    return const_e_succ_range(_succs[v].elts(), _ws);
  }

  const_e_pred_range e_preds(vert_id v) const {
    return const_e_pred_range(_preds[v].elts(), _ws);
  }

  bool is_empty(void) const { return edge_count == 0; }

  size_t size(void) const { return _succs.size(); }

  size_t num_edges(void) const { return edge_count; }

  vert_id new_vertex(void) {
    vert_id v;
    if (free_id.size() > 0) {
      v = free_id.last();
      assert(v < _succs.size());
      free_id.pop();
      is_free[v] = false;
    } else {
      v = _succs.size();
      is_free.push_back(false);
      _succs.push();
      _preds.push();
    }

    return v;
  }

  void growTo(vert_id v) {
    while (size() < v)
      new_vertex();
  }

  void forget(vert_id v) {
    if (is_free[v])
      return;

    for (smap_t::elt_t &e : _succs[v].elts()) {
      free_widx.push(e.val);
      _preds[e.key].remove(v);
    }
    edge_count -= _succs[v].size();
    _succs[v].clear();

    for (smap_t::key_t k : _preds[v].keys())
      _succs[k].remove(v);
    edge_count -= _preds[v].size();
    _preds[v].clear();

    is_free[v] = true;
    free_id.push(v);
  }

  void clear_edges(void) {
    _ws.clear();
    for (vert_id v : verts()) {
      _succs[v].clear();
      _preds[v].clear();
    }
    edge_count = 0;
  }
  void clear(void) {
    _ws.clear();
    _succs.clear();
    _preds.clear();
    is_free.clear();
    free_id.clear();
    free_widx.clear();

    edge_count = 0;
  }

  bool elem(vert_id s, vert_id d) const { return _succs[s].elem(d); }

  Wt &edge_val(vert_id s, vert_id d) {
    size_t idx;
    _succs[s].lookup(d, &idx);
    return _ws[idx];
  }

  const Wt &edge_val(vert_id s, vert_id d) const {
    size_t idx;
    _succs[s].lookup(d, &idx);
    return _ws[idx];
  }

  bool lookup(vert_id s, vert_id d, wt_ref_t &w) const {
    size_t idx;
    if (_succs[s].lookup(d, &idx)) {
      w = wt_ref_t(_ws[idx]);
      return true;
    }
    return false;
  }

  void add_edge(vert_id s, Wt w, vert_id d) {
    assert(!elem(s, d));
    size_t idx;
    if (free_widx.size() > 0) {
      idx = free_widx.last();
      free_widx.pop();
      _ws[idx] = w;
    } else {
      idx = _ws.size();
      _ws.push(w);
    }

    _succs[s].add(d, idx);
    _preds[d].add(s, idx);
    edge_count++;
  }

  template <class Op> void update_edge(vert_id s, Wt w, vert_id d, Op &op) {
    size_t idx;
    if (_succs[s].lookup(d, &idx)) {
      _ws[idx] = op.apply(_ws[idx], w);
    } else {
      if (!op.default_is_absorbing())
        add_edge(s, w, d);
    }
  }

  void set_edge(vert_id s, Wt w, vert_id d) {
    size_t idx;
    if (_succs[s].lookup(d, &idx)) {
      _ws[idx] = w;
    } else {
      add_edge(s, w, d);
    }
  }

  void set_edge_if_less_than(vert_id s, Wt w, vert_id d) {
    size_t idx;
    if (_succs[s].lookup(d, &idx)) {
      if (w < _ws[idx]) {
        _ws[idx] = w;
      }
    } else {
      add_edge(s, w, d);
    }
  }

  void write(crab_os &o) const {
    o << "[|";
    bool first = true;
    for (vert_id v : verts()) {
      auto it = e_succs(v).begin();
      auto end = e_succs(v).end();

      if (it != end) {
        if (first)
          first = false;
        else
          o << ", ";

        o << "[v" << v << " -> ";
        o << "(" << (*it).val << ":" << (*it).vert << ")";
        for (++it; it != end; ++it) {
          o << ", (" << (*it).val << ":" << (*it).vert << ")";
        }
        o << "]";
      }
    }
    o << "|]";
  }

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const AdaptGraph<Weight> &g) {
    g.write(o);
    return o;
  }
};
} // namespace crab
#pragma GCC diagnostic pop
