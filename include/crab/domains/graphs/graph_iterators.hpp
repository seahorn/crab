#pragma once

#include <memory>
#include <vector>

/** Graph iterators **/

namespace crab {
namespace graph_iterators {

// This can be a template parameter of ver_iterator
using vert_id = unsigned int;

/** Iterators over graph nodes **/

class vert_iterator {
public:
  vert_iterator(vert_id _v, const std::vector<bool> &_is_free)
      : v(_v), is_free(_is_free) {}
  vert_id operator*(void) const { return v; }
  vert_iterator &operator++(void) {
    ++v;
    return *this;
  }
  vert_iterator &operator--(void) {
    --v;
    return *this;
  }
  bool operator!=(const vert_iterator &o) {
    while (v < o.v && is_free[v])
      ++v;
    return v < o.v;
  }

private:
  vert_id v;
  const std::vector<bool> &is_free;
};

class vert_range {
public:
  vert_range(vert_id _sz, const std::vector<bool> &_is_free)
      : sz(_sz), is_free(_is_free) {}
  vert_iterator begin(void) const { return vert_iterator(0, is_free); }
  vert_iterator end(void) const { return vert_iterator(sz, is_free); }
  unsigned int size(void) const { return (unsigned int)sz; }

private:
  vert_id sz;
  const std::vector<bool> &is_free;
};

// FIXME(more reusable): it assumes that Pred is a std container but
// patricia trees/sets have different API
template <class Pred, class PredIt> class pred_range {
public:
  using iterator = PredIt;

  pred_range(Pred &_p) : p(_p) {}
  iterator begin(void) const { return iterator(p.begin()); }
  iterator end(void) const { return iterator(p.end()); }
  size_t size(void) const { return p.size(); }

  bool mem(unsigned int v) const { return p.find(v) != p.end(); }
  void add(unsigned int v) { p.insert(v); }
  void remove(unsigned int v) { p.erase(v); }
  void clear() { p.clear(); }

private:
  Pred &p;
};

template <class Pred, class PredIt> class const_pred_range {
public:
  using iterator = PredIt;
  const_pred_range(const Pred &_p) : p(_p) {}
  iterator begin(void) const { return iterator(p.begin()); }
  iterator end(void) const { return iterator(p.end()); }
  size_t size(void) const { return p.size(); }
  bool mem(unsigned int v) const { return p.find(v) != p.end(); }

private:
  const Pred &p;
};

// FIXME(more reusable): it assumes that Succ is a std container but
// patricia trees/sets have different API
template <class Succ, class SuccIt> class succ_range {
  using succ_elt_t = typename Succ::value_type;
  using Wt = typename Succ::mapped_type;

public:
  using iterator = SuccIt;

  succ_range(Succ &_p) : p(_p) {}
  iterator begin(void) const { return iterator(p.begin()); }
  iterator end(void) const { return iterator(p.end()); }
  size_t size(void) const { return p.size(); }

  bool mem(unsigned int v) const { return p.find(v) != p.end(); }
  void add(unsigned int v, const Wt &w) { p.insert(succ_elt_t(v, w)); }
  const Wt &value(unsigned int v) const { return (*(p.find(v))).second; }
  Wt &value(unsigned int v) { return (*(p.find(v))).second; }
  void remove(unsigned int v) { p.erase(v); }
  void clear() { p.clear(); }

private:
  Succ &p;
};

template <class Succ, class SuccIt> class const_succ_range {
  using succ_elt_t = typename Succ::value_type;
  using Wt = typename Succ::mapped_type;

public:
  using iterator = SuccIt;
  const_succ_range(const Succ &_p) : p(_p) {}
  iterator begin(void) const { return iterator(p.begin()); }
  iterator end(void) const { return iterator(p.end()); }
  size_t size(void) const { return p.size(); }
  bool mem(unsigned int v) const { return p.find(v) != p.end(); }
  void add(unsigned int v, const Wt &w) { p.insert(succ_elt_t(v, w)); }
  const Wt &value(unsigned int v) const { return (*(p.find(v))).second; }

private:
  const Succ &p;
};

/** Iterators over graph edges **/

template <class Wt> class edge_ref_t {
public:
  edge_ref_t(vert_id _v, Wt &_w) : vert(_v), val(_w) {}
  vert_id vert;
  Wt &val;
};

template <class Wt> class const_edge_ref_t {
public:
  const_edge_ref_t(vert_id _v, const Wt &_w) : vert(_v), val(_w) {}
  vert_id vert;
  const Wt &val;
};

template <class Wt> class edge_val_t {
public:
  edge_val_t(vert_id _v, Wt _w) : vert(_v), val(_w) {}
  vert_id vert;
  Wt val; // no longer a ref
};

template <class G, class It, class E> class fwd_edge_iterator {
public:
  using Wt = typename G::Wt;
  using edge_wrapper = E;
  using type = fwd_edge_iterator<G, It, E>;

  fwd_edge_iterator(void) : g(nullptr) {}
  fwd_edge_iterator(G &_g, vert_id _s, It _it) : g(&_g), s(_s), it(_it) {}

  // XXX: to make sure that we always return the same address
  // for the "empty" iterator, otherwise we can trigger
  // undefined behavior.
  static type empty_iterator() {
    static std::unique_ptr<type> it = nullptr;
    if (!it)
      it = std::unique_ptr<type>(new type());
    return *it;
  }
  edge_wrapper operator*(void) const {
    return edge_wrapper((*it), g->edge_val(s, (*it)));
  }
  type &operator++(void) {
    ++it;
    return *this;
  }
  bool operator!=(const type &o) { return it != o.it; }

  G *g;
  vert_id s;
  It it;
};

template <class G, class It, class E> class const_fwd_edge_iterator {
public:
  using Wt = typename G::Wt;
  using edge_wrapper = E;
  using type = const_fwd_edge_iterator<G, It, E>;

  const_fwd_edge_iterator(void) : g(nullptr) {}
  const_fwd_edge_iterator(const G &_g, vert_id _s, It _it)
      : g(&_g), s(_s), it(_it) {}
  static type empty_iterator() {
    static std::unique_ptr<type> it = nullptr;
    if (!it)
      it = std::unique_ptr<type>(new type());
    return *it;
  }
  edge_wrapper operator*(void) const {
    return edge_wrapper((*it), g->edge_val(s, (*it)));
  }
  type &operator++(void) {
    ++it;
    return *this;
  }
  bool operator!=(const type &o) { return it != o.it; }

  const G *g;
  vert_id s;
  It it;
};

template <class G, class It, class E> class rev_edge_iterator {
public:
  using Wt = typename G::Wt;
  using edge_wrapper = E;
  using type = rev_edge_iterator<G, It, E>;

  rev_edge_iterator(void) : g(nullptr) {}
  rev_edge_iterator(G &_g, vert_id _d, It _it) : g(&_g), d(_d), it(_it) {}

  static type empty_iterator() {
    static std::unique_ptr<type> it = nullptr;
    if (!it)
      it = std::unique_ptr<type>(new type());
    return *it;
  }

  edge_wrapper operator*(void) const {
    return edge_wrapper((*it), g->edge_val((*it), d));
  }
  type &operator++(void) {
    ++it;
    return *this;
  }
  bool operator!=(const type &o) { return it != o.it; }

  G *g;
  vert_id d;
  It it;
};

template <class G, class It, class E> class const_rev_edge_iterator {
public:
  using Wt = typename G::Wt;
  using edge_wrapper = E;
  using type = const_rev_edge_iterator<G, It, E>;

  const_rev_edge_iterator(void) : g(nullptr) {}
  const_rev_edge_iterator(const G &_g, vert_id _d, It _it)
      : g(&_g), d(_d), it(_it) {}

  static type empty_iterator() {
    static std::unique_ptr<type> it = nullptr;
    if (!it)
      it = std::unique_ptr<type>(new type());
    return *it;
  }

  edge_wrapper operator*(void) const {
    return edge_wrapper((*it), g->edge_val((*it), d));
  }
  type &operator++(void) {
    ++it;
    return *this;
  }
  bool operator!=(const type &o) { return it != o.it; }

  const G *g;
  vert_id d;
  It it;
};

template <class G, class It, class E> class fwd_edge_range {
public:
  using iterator = fwd_edge_iterator<G, It, E>;
  fwd_edge_range(G &_g, vert_id _s) : g(_g), s(_s) {}

  iterator begin(void) const { return iterator(g, s, g.succs(s).begin()); }
  iterator end(void) const { return iterator(g, s, g.succs(s).end()); }
  G &g;
  vert_id s;
};

template <class G, class It, class E> class const_fwd_edge_range {
public:
  using iterator = const_fwd_edge_iterator<G, It, E>;

  const_fwd_edge_range(const G &_g, vert_id _s) : g(_g), s(_s) {}

  iterator begin(void) const { return iterator(g, s, g.succs(s).begin()); }
  iterator end(void) const { return iterator(g, s, g.succs(s).end()); }
  const G &g;
  vert_id s;
};

template <class G, class It, class E> class rev_edge_range {
public:
  using iterator = rev_edge_iterator<G, It, E>;
  rev_edge_range(G &_g, vert_id _d) : g(_g), d(_d) {}

  iterator begin(void) const { return iterator(g, d, g.preds(d).begin()); }
  iterator end(void) const { return iterator(g, d, g.preds(d).end()); }
  G &g;
  vert_id d;
};

template <class G, class It, class E> class const_rev_edge_range {
public:
  using iterator = const_rev_edge_iterator<G, It, E>;
  const_rev_edge_range(const G &_g, vert_id _d) : g(_g), d(_d) {}

  iterator begin(void) const { return iterator(g, d, g.preds(d).begin()); }
  iterator end(void) const { return iterator(g, d, g.preds(d).end()); }
  const G &g;
  vert_id d;
};

} // end namespace graph_iterators
} // end namespace crab
