#pragma once

#include <crab/domains/graphs/graph_ops.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {

namespace domains {

// This class generalizes GrOps from graph_ops.hpp to a more general
// semiring (GrOps only deals with integers).
template <class Gr, bool IsDistWt> class SemiringGrOps {
  //
  // IsDistWt must be true iff
  //   Wt_meet.apply(x, Wt_join.apply(y,z)) =
  //   Wt_join.apply(Wt_meet.apply(x,y),Wt_meet.apply(x,z))
  //
  using graph_t = Gr;
  using Wt = typename graph_t::Wt;
  using GrRev = GraphRev<graph_t>;

  using vert_id = typename graph_t::vert_id;
  using mut_val_ref_t = typename graph_t::mut_val_ref_t;
  using edge_vector = std::vector<std::pair<std::pair<vert_id, vert_id>, Wt>>;

  static Wt make_wt_top() {
    Wt w;
    return w.make_top();
  }

  static Wt make_wt_bottom() {
    Wt w;
    return w.make_bottom();
  }

public:
  static bool less_than(const Wt &w1, const Wt &w2) { return (!(w2 <= w1)); }

  struct Wt_join {
    Wt apply(const Wt &x, const Wt &y) { return x | y; }
    bool default_is_absorbing() { return true; }
  };

  struct Wt_meet {
    Wt apply(const Wt &x, const Wt &y) { return x & y; }
    bool default_is_absorbing() { return false; }
  };

  template <class G> static void set_edge(G &g, vert_id i, Wt w, vert_id j) {
    if (w.is_top())
      return;
    g.set_edge(i, w, j);
  }

  template <class G> static void add_edge(G &g, vert_id i, Wt w, vert_id j) {
    if (w.is_top())
      return;
    g.add_edge(i, w, j);
  }

  template <class G> static void apply_delta(G &g, edge_vector &delta) {
    for (auto &e : delta)
      g.set_edge(e.first.first, e.second, e.first.second);
  }

  // Syntactic join.
  template <class G1, class G2> static graph_t join(const G1 &l, const G2 &r) {
    assert(l.size() == r.size());
    int sz = l.size();

    graph_t g;
    g.growTo(sz);

    mut_val_ref_t wr;
    for (vert_id s : l.verts()) {
      for (auto e : l.e_succs(s)) {
        vert_id d = e.vert;
        if (r.lookup(s, d, &wr))
          add_edge(g, s, e.val | wr.get(), d);
      }
    }
    return g;
  }

  // Syntactic meet/narrowing
  // apply meet if is_meet is true, otherwise apply narrowing
  template <class G1, class G2>
  static graph_t meet_or_narrowing(const G1 &l, const G2 &r, const bool is_meet,
                                   std::vector<vert_id> &changed) {
    assert(l.size() == r.size());
    graph_t g(graph_t::copy(l));
    mut_val_ref_t wg;
    for (vert_id s : r.verts()) {
      for (auto e : r.e_succs(s)) {
        if (!g.lookup(s, e.vert, &wg)) {
          // XXX: this is correct even if narrowing since we just
          // refine from a top weight
          add_edge(g, s, e.val, e.vert);
        } else {
          wg = (is_meet ? e.val & wg.get() : e.val && wg.get());
        }
      }
      // Check if this vertex is stable
      for (auto e : l.e_succs(s)) {
        vert_id d = e.vert;
        if (g.lookup(s, d, &wg)) {
          if (less_than(wg.get(), e.val)) {
            changed.push_back(s);
            break;
          }
        } else {
          changed.push_back(s);
          break;
        }
      }
      for (auto e : r.e_succs(s)) {
        vert_id d = e.vert;
        if (g.lookup(s, d, &wg)) {
          if (less_than(wg.get(), e.val)) {
            changed.push_back(s);
            break;
          }
        } else {
          changed.push_back(s);
          break;
        }
      }
    }
    return g;
  }

  // Syntactic widening
  template <class G1, class G2>
  static graph_t widen(const G1 &l, const G2 &r,
                       std::vector<vert_id> &unstable) {
    assert(l.size() == r.size());
    size_t sz = l.size();
    graph_t g;
    g.growTo(sz);
    mut_val_ref_t wl;
    for (vert_id s : r.verts()) {
      for (auto e : r.e_succs(s)) {
        vert_id d = e.vert;
        if (l.lookup(s, d, &wl))
          add_edge(g, s, wl.get() || e.val, d);
      }
      // Check if this vertex is stable
      mut_val_ref_t widen_wl;
      for (auto e : l.e_succs(s)) {
        vert_id d = e.vert;
        if (g.lookup(s, d, &widen_wl)) {
          if (less_than(e.val, widen_wl.get())) {
            unstable.push_back(s);
            break;
          }
        } else {
          unstable.push_back(s);
          break;
        }
      }
    }
    return g;
  }

  template <class G>
  static bool maybe_extend(G &g, vert_id i, vert_id j, Wt w_ik, Wt w_kj) {
    mut_val_ref_t w0;
    if (g.lookup(i, j, &w0)) {
      // w := w0 & (w_ik | w_kj)
      Wt w = (w0.get() & w_ik) | (w0.get() & w_kj);
      if (less_than(w, w0.get())) {
        w0 = w;
        return true;
      }
    } else {
      auto w_ij = w_ik | w_kj;
      if (!w_ij.is_top()) {
        set_edge(g, i, w_ij, j);
        return true;
      }
    }
    return false;
  }

  // Produce a closed graph even if weight is not distributive.
  // XXX: we only use it for debugging
  template <class G> static void closure(G &g, unsigned max_cycles = 10) {
    if (IsDistWt) {
      floyd_warshall(g);
    } else {
      bool change = false;
      unsigned cycle = 0;
      do {
        ++cycle;
        change = floyd_warshall(g);
      } while (change && cycle <= max_cycles);
    }
  }

  // XXX: we only use it for debugging
  template <class G> static bool floyd_warshall(G &g) {
    bool change = false;
    mut_val_ref_t w_ij, w_ik, w_kj;
    for (vert_id k : g.verts()) {
      for (vert_id i : g.verts()) {
        for (vert_id j : g.verts()) {
          bool has_ik = g.lookup(i, k, &w_ik);
          bool has_kj = g.lookup(k, j, &w_kj);
          change = maybe_extend(g, i, j, (has_ik ? w_ik.get() : make_wt_top()),
                                (has_kj ? w_kj.get() : make_wt_top()));
        }
      }
    }
    return change;
  }

  // restore closure after adding edge (x,y)
  // XXX: if weight is not distributive then closure is not
  // guaranteed.
  template <class G> static bool close_after_edge(G &g, vert_id x, vert_id y) {

    mut_val_ref_t w;
    if (!g.lookup(x, y, &w))
      return false;

    // find vertices that change
    std::vector<vert_id> Q1, Q2;
    Q1.push_back(x);
    Q2.push_back(y);
    for (vert_id i : g.verts()) {
      mut_val_ref_t w_ix, w_iy, w_yi, w_xi;
      if (g.lookup(i, x, &w_ix)) {
        if (g.lookup(i, y, &w_iy)) {
          if (less_than(w_ix.get() | w.get(), w_iy.get()))
            Q1.push_back(i);
        } else
          Q1.push_back(i);
      }
      if (g.lookup(y, i, &w_yi)) {
        if (g.lookup(x, i, &w_xi)) {
          if (less_than(w_yi.get() | w.get(), w_xi.get()))
            Q2.push_back(i);
        } else
          Q2.push_back(i);
      }
    }

    bool change = false;
    mut_val_ref_t w_ij, w_ix, w_yj;
    for (auto i : Q1) {
      for (auto j : Q2) {
        if (y == j && g.lookup(i, x, &w_ix)) {
          change = maybe_extend(g, i, j, w_ix.get(), w.get());
        } else if (i == x && g.lookup(y, j, &w_yj)) {
          change = maybe_extend(g, i, j, w.get(), w_yj.get());
        } else if (g.lookup(i, x, &w_ix) && g.lookup(y, j, &w_yj)) {
          change = maybe_extend(g, i, j, w_ix.get() | w.get(), w_yj.get());
        }
      }
    }
    return change;
  }

  // Comparator for use with min-heaps.
  template <class V> class DistComp {
  public:
    DistComp(V &_A) : A(_A) {}
    bool operator()(int x, int y) const { return !(A[y] <= A[x]); }
    V &A;
  };
  using WtComp = DistComp<std::vector<Wt>>;
  using WtHeap = Heap<WtComp>;

  /// FIXME: Avoid expanding stable vertices
  template <typename G>
  static void dijkstra(const G &g, vert_id s,
                       std::vector<std::pair<vert_id, Wt>> &out) {

    std::vector<Wt> dists;
    while (dists.size() < g.size())
      dists.push_back(make_wt_top());
    dists[s] = make_wt_bottom();

    WtComp comp(dists);
    WtHeap heap(comp);
    heap.insert(s);

    while (!heap.empty()) {
      vert_id q = heap.removeMin();
      for (auto e : g.e_succs(q)) {
        vert_id d = e.vert;
        Wt w = dists[q] | e.val;
        if (less_than(w, dists[d])) {
          dists[d] = w;
          if (heap.inHeap(d))
            heap.decrease(d);
          else
            heap.insert(d);
        }
      }
    }

    for (vert_id i = 0, e = dists.size(); i != e; ++i) {
      if (!dists[i].is_top())
        out.push_back(std::make_pair(i, dists[i]));
    }
  }

  template <class G>
  static void close_over_vertex(G &g, vert_id v, edge_vector &delta) {
    std::vector<std::pair<vert_id, Wt>> aux;
    dijkstra(g, v, aux);
    for (auto p : aux)
      delta.push_back(std::make_pair(std::make_pair(v, p.first), p.second));
    aux.clear();
    GraphRev<Gr> g_rev(g);
    dijkstra(g_rev, v, aux);
    for (auto p : aux)
      delta.push_back(std::make_pair(std::make_pair(p.first, v), p.second));
  }

  // restore closure after widening
  // XXX: if weight is not distributive then closure is not
  // guaranteed.
  template <class G, class V>
  static void close_after_widen(G &g, const V &is_stable) {
    edge_vector delta;
    for (vert_id v : g.verts()) {
      if (is_stable[v])
        continue;
      close_over_vertex(g, v, delta);
    }
    apply_delta(g, delta);
  }

  // restore closure after meet/narrowing
  // XXX: if weight is not distributive then closure is not
  // guaranteed.
  template <class G, class V>
  static void close_after_meet_or_narrowing(G &g, const V &is_stable) {
    edge_vector delta;
    for (vert_id v : g.verts()) {
      if (is_stable[v])
        continue;
      close_over_vertex(g, v, delta);
    }
    apply_delta(g, delta);
  }
};
} // end namespace domains
} // end namespace crab
#pragma GCC diagnostic pop
