#pragma once

#include <boost/optional.hpp>
#include <vector>

/**
 * Graph views - to traverse some mutation of the graph without
 * actually constructing it.
 **/

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {

template <class V> class num_range {
public:
  class value_ref {
  public:
    value_ref(const V &_v) : v(_v) {}
    const V &operator*(void) const { return v; }
    value_ref &operator++(void) {
      v++;
      return *this;
    }
    value_ref &operator--(void) {
      v--;
      return *this;
    }
    bool operator!=(const value_ref &o) const { return v < o.v; }

  protected:
    V v;
  };
  num_range(const V &_after) : after(_after) {}

  value_ref begin(void) const { return value_ref((V)0); }
  value_ref end(void) const { return value_ref(after); }

protected:
  V after;
};

// Processing a graph under a (possibly incomplete)
// permutation of vertices.
// We assume perm[x] is unique; otherwise, we'd have
// to introduce a edges for induced equivalence classes.
template <class G> class GraphPerm {
public:
  using vert_id = typename G::vert_id;
  using Wt = typename G::Wt;
  using wt_ref_t = typename G::wt_ref_t;
  // These defined below
  // using const_pred_range = ...
  // using const_succ_range = ...
  // using const_e_pred_range = ...
  // using const_e_succ_range = ...

  GraphPerm(std::vector<vert_id> &_perm, const G &_g)
      : g(_g), perm(_perm), inv(_g.size(), -1) {
    for (unsigned int vi = 0; vi < perm.size(); vi++) {
      if (perm[vi] == -1)
        continue;
      assert(inv[perm[vi]] == -1);
      inv[perm[vi]] = vi;
    }
  }

  // Check whether an edge is live
  bool elem(vert_id x, vert_id y) const {
    if (perm[x] > g.size() || perm[y] > g.size())
      return false;
    return g.elem(perm[x], perm[y]);
  }

  bool lookup(vert_id x, vert_id y, wt_ref_t &w) const {
    if (perm[x] > g.size() || perm[y] > g.size())
      return false;
    return g.lookup(perm[x], perm[y], w);
  }

  // Precondition: elem(x, y) is true.
  Wt edge_val(vert_id x, vert_id y) const {
    // assert(perm[x] < g.size() && perm[y] < g.size());
    return g.edge_val(perm[x], perm[y]);
  }

  // Precondition: elem(x, y) is true.
  Wt operator()(vert_id x, vert_id y) const {
    // assert(perm[x] < g.size() && perm[y] < g.size());
    return g(perm[x], perm[y]);
  }

  // Number of allocated vertices
  int size(void) const { return perm.size(); }

  using vert_range = num_range<vert_id>;
  using vert_iterator = typename num_range<vert_id>::value_ref;
  vert_range verts(void) const { return vert_range(perm.size()); }

  // GKG: Should probably modify this to handle cases where
  // the vertex iterator isn't just a vert_id*.
  template <class ItG> class adj_iterator {
  public:
    adj_iterator(const std::vector<vert_id> &_inv, const ItG &_v)
        : inv(_inv), v(_v) {}

    vert_id operator*(void) const { return inv[*v]; }

    adj_iterator &operator++(void) {
      ++v;
      return *this;
    }

    bool operator!=(const adj_iterator &other) {
      while (v != other.v && inv[*v] == (-1))
        ++v;
      return v != other.v;
    }

  protected:
    const std::vector<vert_id> &inv;
    ItG v;
  };

  template <class ItG> class e_adj_iterator {
  public:
    using edge_wrapper = typename ItG::edge_wrapper;

    e_adj_iterator(const std::vector<vert_id> &_inv, const ItG &_v)
        : inv(_inv), v(_v) {}

    edge_wrapper operator*(void) {
      return edge_wrapper(inv[(*v).vert], (*v).val);
    }

    e_adj_iterator &operator++(void) {
      ++v;
      return *this;
    }

    bool operator!=(const e_adj_iterator &other) {
      while (v != other.v && inv[(*v).vert] == (-1))
        ++v;
      return v != other.v;
    }

  protected:
    const std::vector<vert_id> &inv;
    ItG v;
  };

  template <class RG, class It> class adj_list {
  public:
    using ItG = typename RG::iterator;
    using adj_list_t = adj_list<RG, It>;
    using iterator = It;

    adj_list(const std::vector<vert_id> &_perm,
             const std::vector<vert_id> &_inv, const RG &_adj)
        : perm(_perm), inv(_inv), adj(_adj) {}

    adj_list(const std::vector<vert_id> &_perm,
             const std::vector<vert_id> &_inv)
        : perm(_perm), inv(_inv), adj() {}

    iterator begin(void) const {
      if (adj)
        return iterator(inv, (*adj).begin());
      else
        return iterator(inv, ItG::empty_iterator());
    }
    iterator end(void) const {
      if (adj)
        return iterator(inv, (*adj).end());
      else
        return iterator(inv, ItG::empty_iterator());
    }

    bool mem(unsigned int v) const {
      if (!adj || perm[v] == (-1))
        return false;
      return (*adj).mem(perm[v]);
    }

  protected:
    const std::vector<vert_id> &perm;
    const std::vector<vert_id> &inv;
    boost::optional<RG> adj;
  };

  // public typedefs
  using const_pred_range =
      adj_list<typename G::const_pred_range,
               adj_iterator<typename G::const_pred_range::iterator>>;
  using const_succ_range =
      adj_list<typename G::const_succ_range,
               adj_iterator<typename G::const_succ_range::iterator>>;
  using const_e_pred_range =
      adj_list<typename G::const_e_pred_range,
               e_adj_iterator<typename G::const_e_pred_range::iterator>>;
  using const_e_succ_range =
      adj_list<typename G::const_e_succ_range,
               e_adj_iterator<typename G::const_e_succ_range::iterator>>;

  const_succ_range succs(vert_id v) const {
    if (perm[v] == (-1))
      return const_succ_range(perm, inv);
    else
      return const_succ_range(perm, inv, g.succs(perm[v]));
  }
  const_pred_range preds(vert_id v) const {
    if (perm[v] == (-1))
      return const_pred_range(perm, inv);
    else
      return const_pred_range(perm, inv, g.preds(perm[v]));
  }

  const_e_succ_range e_succs(vert_id v) const {
    if (perm[v] == (-1))
      return const_e_succ_range(perm, inv);
    else
      return const_e_succ_range(perm, inv, g.e_succs(perm[v]));
  }

  const_e_pred_range e_preds(vert_id v) const {
    if (perm[v] == (-1))
      return const_e_pred_range(perm, inv);
    else
      return const_e_pred_range(perm, inv, g.e_preds(perm[v]));
  }

  void write(crab_os &o) const {
    o << "[|";
    bool first = true;
    for (vert_id v : verts()) {
      auto it = e_succs(v).begin();
      auto end = e_succs(v).end();
      if (it != end) {
        if (first) {
          first = false;
        } else {
          o << ", ";
        }
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

  const G &g;
  std::vector<vert_id> perm;
  std::vector<vert_id> inv;
};

// View of a graph, omitting a given vertex
template <class G> class SubGraph {
public:
  using vert_id = typename G::vert_id;
  using Wt = typename G::Wt;
  using wt_ref_t = typename G::wt_ref_t;
  // These defined below
  // using const_pred_range = ...
  // using const_succ_range = ...
  // using const_e_pred_range = ...
  // using const_e_succ_range = ...

  SubGraph(const G &_g, vert_id _v_ex) : g(_g), v_ex(_v_ex) {}

  bool elem(vert_id x, vert_id y) const {
    return (x != v_ex && y != v_ex && g.elem(x, y));
  }

  bool lookup(vert_id x, vert_id y, wt_ref_t &w) const {
    return (x != v_ex && y != v_ex && g.lookup(x, y, w));
  }

  Wt edge_val(vert_id x, vert_id y) const { return g.edge_val(x, y); }

  // Precondition: elem(x, y) is true.
  Wt operator()(vert_id x, vert_id y) const { return g(x, y); }

  // Number of allocated vertices
  int size(void) const { return g.size(); }

  class vert_iterator {
  public:
    vert_iterator(const typename G::vert_iterator &_iG, vert_id _v_ex)
        : v_ex(_v_ex), iG(_iG) {}

    // Skipping of v_ex is done entirely by !=.
    // So we _MUST_ test it != verts.end() before dereferencing.
    vert_id operator*(void) { return *iG; }
    vert_iterator operator++(void) {
      ++iG;
      return *this;
    }
    bool operator!=(const vert_iterator &o) {
      if (iG != o.iG && (*iG) == v_ex)
        ++iG;
      return iG != o.iG;
    }

    vert_id v_ex;
    typename G::vert_iterator iG;
  };
  class vert_range {
  public:
    vert_range(const typename G::vert_range &_rG, vert_id _v_ex)
        : rG(_rG), v_ex(_v_ex) {}

    vert_iterator begin(void) const { return vert_iterator(rG.begin(), v_ex); }
    vert_iterator end(void) const { return vert_iterator(rG.end(), v_ex); }

    typename G::vert_range rG;
    vert_id v_ex;
  };
  vert_range verts(void) const { return vert_range(g.verts(), v_ex); }

  template <class It> class adj_iterator {
  public:
    adj_iterator(const It &_iG, vert_id _v_ex) : iG(_iG), v_ex(_v_ex) {}
    vert_id operator*(void) const { return *iG; }
    adj_iterator &operator++(void) {
      ++iG;
      return *this;
    }
    bool operator!=(const adj_iterator &o) {
      if (iG != o.iG && (*iG) == v_ex)
        ++iG;
      return iG != o.iG;
    }

    It iG;
    vert_id v_ex;
  };

  template <class It> class e_adj_iterator {
  public:
    using edge_wrapper = typename It::edge_wrapper;

    e_adj_iterator(const It &_iG, vert_id _v_ex) : iG(_iG), v_ex(_v_ex) {}
    edge_wrapper operator*(void) const { return *iG; }
    e_adj_iterator &operator++(void) {
      ++iG;
      return *this;
    }
    bool operator!=(const e_adj_iterator &o) {
      if (iG != o.iG && (*iG).vert == v_ex)
        ++iG;
      return iG != o.iG;
    }

    It iG;
    vert_id v_ex;
  };

  template <class R, class It> class adj_list {
  public:
    using g_iter = typename R::iterator;
    using iterator = It;

    adj_list(const R &_rG, vert_id _v_ex) : rG(_rG), v_ex(_v_ex) {}
    iterator begin() const { return iterator(rG.begin(), v_ex); }
    iterator end() const { return iterator(rG.end(), v_ex); }

  protected:
    R rG;
    vert_id v_ex;
  };

  // public typedefs
  using const_pred_range =
      adj_list<typename G::const_pred_range,
               adj_iterator<typename G::const_pred_range::iterator>>;
  using const_succ_range =
      adj_list<typename G::const_succ_range,
               adj_iterator<typename G::const_succ_range::iterator>>;
  using const_e_pred_range =
      adj_list<typename G::const_e_pred_range,
               e_adj_iterator<typename G::const_e_pred_range::iterator>>;
  using const_e_succ_range =
      adj_list<typename G::const_e_succ_range,
               e_adj_iterator<typename G::const_e_succ_range::iterator>>;

  const_succ_range succs(vert_id v) const {
    // assert(v != v_ex);
    return const_succ_range(g.succs(v), v_ex);
  }
  const_pred_range preds(vert_id v) const {
    // assert(v != v_ex);
    return const_pred_range(g.preds(v), v_ex);
  }
  const_e_succ_range e_succs(vert_id v) const {
    return const_e_succ_range(g.e_succs(v), v_ex);
  }
  const_e_pred_range e_preds(vert_id v) const {
    return const_e_pred_range(g.e_preds(v), v_ex);
  }

  const G &g;
  vert_id v_ex;
};

// Viewing a graph with all edges reversed.
// Useful if we want to run single-dest shortest paths,
// for updating bounds and incremental closure.
template <class G> class GraphRev {
public:
  using vert_id = typename G::vert_id;
  using Wt = typename G::Wt;
  using wt_ref_t = typename G::wt_ref_t;
  using const_pred_range = typename G::const_succ_range;
  using const_succ_range = typename G::const_pred_range;
  using const_e_pred_range = typename G::const_e_succ_range;
  using const_e_succ_range = typename G::const_e_pred_range;

  GraphRev(const G &_g) : g(_g) {}

  // Check whether an edge is live
  bool elem(vert_id x, vert_id y) const { return g.elem(y, x); }

  bool lookup(vert_id x, vert_id y, wt_ref_t &w) const {
    return g.lookup(y, x, w);
  }

  // Precondition: elem(x, y) is true.
  Wt edge_val(vert_id x, vert_id y) const { return g.edge_val(y, x); }

  // Precondition: elem(x, y) is true.
  Wt operator()(vert_id x, vert_id y) const { return g(y, x); }

  // Number of allocated vertices
  int size(void) const { return g.size(); }

  typename G::vert_range verts(void) const { return g.verts(); }

  const_succ_range succs(vert_id v) const { return g.preds(v); }
  const_succ_range preds(vert_id v) const { return g.succs(v); }

  const_e_succ_range e_succs(vert_id v) const { return g.e_preds(v); }
  const_e_pred_range e_preds(vert_id v) const { return g.e_succs(v); }

  const G &g;
};

} // namespace crab
#pragma GCC diagnostic pop
