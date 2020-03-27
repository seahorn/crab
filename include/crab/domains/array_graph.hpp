/*******************************************************************************
 * An array content domain.
 *
 * This domain is a simplified implementation of the paper:
 * "A Partial-Order Approach to Array Content Analysis" by
 *  Gange, Navas, Schachte, Sondergaard, and Stuckey
 *  available here http://arxiv.org/pdf/1408.1754v1.pdf.
 *
 * It keeps a graph where vertices are array indexes and edges are
 * labelled with weights.  An edge (i,j) with weight w denotes that
 * the property w holds for the array segment [i,j). A weight is an
 * arbitrary lattice that can relate multiple array variables as well
 * as array with scalar variables.
 ******************************************************************************/

/**
 * TODOs:
 *
 * The implementation works for toy programs but we need to fix the
 * following issues for being able to analyze real programs:
 *
 * - landmarks must be kept as local state as part of each abstract
 *   state.
 *
 * - reduction between scalar and weight domains must be done
 *   incrementally. For that, we need some assumptions about the
 *   underlying scalar domain. For instance, if we assume zones then
 *   after each operation we know which are the indexes affected by
 *   the operation. We can use that information for doing reduction
 *   only on those indexes. This would remove the need of having
 *   methods such as array_graph_domain_helper_traits::is_unsat and
 *   array_graph_domain_helper_traits::active_variables which are
 *   anyway domain dependent.
 **/

#pragma once

#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/array_graph/array_graph_ops.hpp>
#include <crab/domains/array_graph/array_segmentation.hpp>
#include <crab/domains/graphs/adapt_sgraph.hpp>
#include <crab/domains/graphs/sparse_graph.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/linear_constraints.hpp>

// XXX: if expression domain is a template parameter no need to include
#include <crab/domains/term_equiv.hpp>
// XXX: for customized propagations between weight and scalar domains
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/nullity.hpp>

#include <boost/container/flat_map.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/optional.hpp>
#include <boost/range/iterator_range.hpp>
#include <functional>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <unordered_set>

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {
namespace domains {
/*
   A weighted directed graph where the weight is an abstract
   domain. The graph should be always kept in a consistent form,
   i.e., for all i,j,k:: weight(i,j) <= join(weight(i,k), weight(k,j))
*/
template <typename Vertex, typename Weight, bool IsDistWeight>
class array_graph_ : public writeable {

public:
  // XXX: make this a template parameter later
  // typedef AdaptGraph<Weight> graph_t;
  typedef SparseWtGraph<Weight> graph_t;
  typedef typename graph_t::vert_id _vert_id;
  typedef typename graph_t::edge_ref_t edge_ref_t;
  typedef typename graph_t::Wt Wt;
  typedef typename graph_t::mut_val_ref_t mut_val_ref_t;
  // XXX: needs to have this typedef so we can use GraphRev
  typedef Vertex vert_id;

private:
  typedef GraphPerm<graph_t> GrPerm;
  typedef ArrayGrOps<graph_t, IsDistWeight> GrOps;
  typedef typename GrOps::Wt_join Wt_join;
  typedef typename GrOps::Wt_meet Wt_meet;

  typedef boost::container::flat_map<Vertex, _vert_id> vert_map_t;
  typedef typename vert_map_t::value_type vmap_elt_t;
  typedef std::vector<boost::optional<Vertex>> rev_map_t;
  typedef std::unordered_set<_vert_id> vert_set_t;
  typedef array_graph_<Vertex, Weight, IsDistWeight> array_graph_t;

  vert_map_t _vert_map;
  rev_map_t _rev_map;
  graph_t _g;
  vert_set_t _unstable;
  bool _is_bottom;

  struct vert_set_wrap_t {
    vert_set_wrap_t(const vert_set_t &_vs) : vs(_vs) {}

    bool operator[](_vert_id v) const { return vs.find(v) != vs.end(); }
    const vert_set_t &vs;
  };

  _vert_id get_vert(Vertex v) {
    auto it = _vert_map.find(v);
    if (it != _vert_map.end())
      return (*it).second;

    _vert_id vert(_g.new_vertex());
    assert(vert <= _rev_map.size());

    if (vert < _rev_map.size()) {
      assert(!_rev_map[vert]);
      _rev_map[vert] = v;
    } else {
      _rev_map.push_back(v);
    }
    _vert_map.insert(vmap_elt_t(v, vert));

    return vert;
  }

public:
  template <class ItS> class iterator {
  public:
    typedef iterator<ItS> iter_t;
    iterator(const ItS &_it, const rev_map_t &_rev_map)
        : it(_it), rev_map(_rev_map) {}
    bool operator!=(const iter_t &o) { return it != o.it; }
    iter_t &operator++(void) {
      ++it;
      return *this;
    }
    Vertex operator*(void)const {
      if (!rev_map[*it])
        CRAB_ERROR("Reverse map failed");
      return *(rev_map[*it]);
    }

  protected:
    ItS it;
    const rev_map_t &rev_map;
  };

  struct edge_t {
    edge_t(Vertex _v, Wt &_w) : vert(_v), val(_w) {}
    Vertex vert;
    Wt &val;
  };

  template <class ItS> class edge_iterator {
  public:
    typedef edge_iterator<ItS> iter_t;
    edge_iterator(const ItS &_it, const rev_map_t &_rev_map)
        : it(_it), rev_map(_rev_map) {}
    bool operator!=(const iter_t &o) { return it != o.it; }
    iter_t &operator++(void) {
      ++it;
      return *this;
    }
    edge_t operator*(void)const {
      edge_ref_t e = *it;
      if (!rev_map[e.vert])
        CRAB_ERROR("Reverse map failed");
      Vertex v = *(rev_map[e.vert]);
      return edge_t(v, e.val);
    }

  protected:
    ItS it;
    const rev_map_t &rev_map;
  };

  template <class Range, class ItS> class iterator_range {
  public:
    typedef ItS iterator;
    iterator_range(const Range &r, const rev_map_t &rev_map)
        : _r(r), _rev_map(rev_map) {}
    iterator begin(void) const { return iterator(_r.begin(), _rev_map); }
    iterator end(void) const { return iterator(_r.end(), _rev_map); }

  protected:
    Range _r;
    const rev_map_t &_rev_map;
  };

  typedef iterator<typename graph_t::succ_iterator> succ_transform_iterator;
  typedef iterator<typename graph_t::pred_iterator> pred_transform_iterator;
  typedef iterator<typename graph_t::vert_iterator> vert_transform_iterator;
  typedef edge_iterator<typename graph_t::fwd_edge_iterator>
      fwd_edge_transform_iterator;
  typedef edge_iterator<typename graph_t::rev_edge_iterator>
      rev_edge_transform_iterator;

  typedef iterator_range<typename graph_t::succ_range, succ_transform_iterator>
      succ_range;
  typedef iterator_range<typename graph_t::pred_range, pred_transform_iterator>
      pred_range;
  typedef iterator_range<typename graph_t::vert_range, vert_transform_iterator>
      vert_range;
  typedef iterator_range<typename graph_t::e_succ_range,
                         fwd_edge_transform_iterator>
      e_succ_range;
  typedef iterator_range<typename graph_t::e_pred_range,
                         rev_edge_transform_iterator>
      e_pred_range;

  vert_range verts() {
    typename graph_t::vert_range p = _g.verts();
    return vert_range(p, _rev_map);
  }

  succ_range succs(Vertex v) {
    typename graph_t::succ_range p = _g.succs(get_vert(v));
    return succ_range(p, _rev_map);
  }

  pred_range preds(Vertex v) {
    typename graph_t::pred_range p = _g.preds(get_vert(v));
    return pred_range(p, _rev_map);
  }

  e_succ_range e_succs(Vertex v) {
    typename graph_t::e_succ_range p = _g.e_succs(get_vert(v));
    return e_succ_range(p, _rev_map);
  }

  e_pred_range e_preds(Vertex v) {
    typename graph_t::e_pred_range p = _g.e_preds(get_vert(v));
    return e_pred_range(p, _rev_map);
  }

public:
  array_graph_(bool is_bottom = false) : _is_bottom(is_bottom) {}

  array_graph_(const array_graph_t &o)
      : _vert_map(o._vert_map), _rev_map(o._rev_map), _g(o._g),
        _unstable(o._unstable), _is_bottom(false) {
    if (o._is_bottom)
      set_to_bottom();
  }

  array_graph_(array_graph_t &&o)
      : _vert_map(std::move(o._vert_map)), _rev_map(std::move(o._rev_map)),
        _g(std::move(o._g)), _unstable(std::move(o._unstable)),
        _is_bottom(o._is_bottom) {}

  array_graph_(vert_map_t &vert_map, rev_map_t &rev_map, graph_t &g,
               vert_set_t unstable)
      : _vert_map(vert_map), _rev_map(rev_map), _g(g), _unstable(unstable),
        _is_bottom(false) {}

  array_graph_(vert_map_t &&vert_map, rev_map_t &&rev_map, graph_t &&g,
               vert_set_t &&unstable)
      : _vert_map(std::move(vert_map)), _rev_map(std::move(rev_map)),
        _g(std::move(g)), _unstable(std::move(unstable)), _is_bottom(false) {}

  array_graph_t &operator=(const array_graph_t &o) {
    if (this != &o) {
      if (o._is_bottom)
        set_to_bottom();
      else {
        _is_bottom = false;
        _vert_map = o._vert_map;
        _rev_map = o._rev_map;
        _g = o._g;
        _unstable = o._unstable;
      }
    }
    return *this;
  }

  array_graph_t &operator=(array_graph_t &&o) {
    if (o._is_bottom) {
      set_to_bottom();
    } else {
      _is_bottom = false;
      _vert_map = std::move(o._vert_map);
      _rev_map = std::move(o._rev_map);
      _unstable = std::move(o._unstable);
      _g = std::move(o._g);
    }
    return *this;
  }

public:
  void set_to_bottom() {
    _vert_map.clear();
    _rev_map.clear();
    _g.clear();
    _unstable.clear();
    _is_bottom = true;
  }

  static array_graph_t top() { return array_graph_t(false); }

  static array_graph_t bottom() { return array_graph_t(true); }

  bool is_bottom() const { return _is_bottom; }

  bool is_top() {
    if (_is_bottom)
      return false;
    return _g.is_empty();
  }

  bool lookup_edge(Vertex s, Vertex d, mut_val_ref_t *w) {
    if (is_bottom())
      return false;
    auto se = get_vert(s);
    auto de = get_vert(d);
    return _g.lookup(se, de, w);
  }

  // // update edge but do not close graph
  // void update_edge_unclosed(Vertex s, Weight w, Vertex d) {
  //   if (w.is_top()) return;
  //   normalize();
  //   if (is_bottom()) return;
  //   auto se = get_vert(s);
  //   auto de = get_vert(d);
  //   Wt_meet op;
  //   _g.update_edge(se, w, de, op);
  // }

  // close the graph after edge (s,d) has been updated
  void close_edge(Vertex s, Vertex d) {
    normalize();
    if (is_bottom())
      return;
    auto se = get_vert(s);
    auto de = get_vert(d);
    GrOps::close_after_edge(_g, se, de);
  }

  void update_edge(Vertex s, Weight w, Vertex d) {
    if (w.is_top())
      return;

    normalize();

    if (is_bottom())
      return;

    auto se = get_vert(s);
    auto de = get_vert(d);
    Wt_meet op;
    _g.update_edge(se, w, de, op);
    GrOps::close_after_edge(_g, se, de);
  }

  // void full_close() { // for debugging
  //   if (is_bottom()) return;
  //   GrOps::floyd_warshall(_g);
  // }

  void expand(Vertex s, Vertex d) {
    if (is_bottom())
      return;

    auto it = _vert_map.find(d);
    if (it != _vert_map.end()) {
      CRAB_ERROR("array_graph expand failed because vertex ", d,
                 " already exists");
    }

    auto se = get_vert(s);
    auto de = get_vert(d);

    for (auto edge : _g.e_preds(se))
      _g.add_edge(edge.vert, edge.val, de);

    for (auto edge : _g.e_succs(se))
      _g.add_edge(de, edge.val, edge.vert);
  }

  void normalize() {
#if 0
    GrOps::closure(_g); // only for debugging purposes
#else
    // Always maintained in closed form except for widening
    if (_unstable.size() == 0)
      return;
    GrOps::close_after_widen(_g, vert_set_wrap_t(_unstable));
    _unstable.clear();
#endif
  }

  void operator|=(array_graph_t &o) { *this = *this | o; }

  bool operator<=(array_graph_t &o) {
    if (is_bottom())
      return true;
    else if (o.is_bottom())
      return false;
    else if (o.is_top())
      return true;
    else if (is_top())
      return false;
    else {
      normalize();

      if (_vert_map.size() < o._vert_map.size())
        return false;

      // Set up a mapping from o to this.
      std::vector<unsigned int> vert_renaming(o._g.size(), -1);
      for (auto p : o._vert_map) {
        auto it = _vert_map.find(p.first);
        // We can't have this <= o if we're missing some
        // vertex.
        if (it == _vert_map.end())
          return false;
        vert_renaming[p.second] = (*it).second;
      }

      assert(_g.size() > 0);
      mut_val_ref_t wx;

      for (_vert_id ox : o._g.verts()) {
        assert(vert_renaming[ox] != -1);
        _vert_id x = vert_renaming[ox];
        for (auto edge : o._g.e_succs(ox)) {
          _vert_id oy = edge.vert;
          assert(vert_renaming[ox] != -1);
          _vert_id y = vert_renaming[oy];
          auto ow = (Weight)edge.val;
          if (!_g.lookup(x, y, &wx) || (!((Weight)wx <= ow)))
            return false;
        }
      }
      return true;
    }
  }

  array_graph_t operator|(array_graph_t &o) {

    if (is_bottom() || o.is_top())
      return o;
    else if (is_top() || o.is_bottom())
      return *this;
    else {
      CRAB_LOG("array-sgraph", crab::outs() << "Before join:\n"
                                            << "Graph 1\n"
                                            << *this << "\n"
                                            << "Graph 2\n"
                                            << o << "\n");

      normalize();
      o.normalize();

      // Figure out the common renaming.
      std::vector<_vert_id> perm_x;
      std::vector<_vert_id> perm_y;

      vert_map_t out_vmap;
      rev_map_t out_revmap;

      for (auto p : _vert_map) {
        auto it = o._vert_map.find(p.first);
        // Vertex exists in both
        if (it != o._vert_map.end()) {
          out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
          out_revmap.push_back(p.first);

          perm_x.push_back(p.second);
          perm_y.push_back((*it).second);
        }
      }

      // Build the permuted view of x and y.
      assert(_g.size() > 0);
      GrPerm gx(perm_x, _g);
      assert(o._g.size() > 0);
      GrPerm gy(perm_y, o._g);

      // We now have the relevant set of relations. Because g_rx
      // and g_ry are closed, the result is also closed.
      graph_t join_g(GrOps::join(gx, gy));

      // Now garbage collect any unused vertices
      for (_vert_id v : join_g.verts()) {
        if (join_g.succs(v).size() == 0 && join_g.preds(v).size() == 0) {
          join_g.forget(v);
          if (out_revmap[v]) {
            out_vmap.erase(*(out_revmap[v]));
            out_revmap[v] = boost::none;
          }
        }
      }

      array_graph_t res(std::move(out_vmap), std::move(out_revmap),
                        std::move(join_g), vert_set_t());
      CRAB_LOG("array-sgraph", crab::outs() << "Result join:\n"
                                            << res << "\n";);
      return res;
    }
  }

  template <typename Thresholds>
  array_graph_t widening_thresholds(array_graph_t &o, const Thresholds &) {
    return (*this || o);
  }

  array_graph_t operator||(array_graph_t &o) {
    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else {
      CRAB_LOG("array-sgraph", crab::outs() << "Before widening:\n"
                                            << "Graph 1\n"
                                            << *this << "\n"
                                            << "Graph 2\n"
                                            << o << "\n";);
      o.normalize();

      // Figure out the common renaming
      std::vector<_vert_id> perm_x;
      std::vector<_vert_id> perm_y;
      vert_map_t out_vmap;
      rev_map_t out_revmap;
      vert_set_t widen_unstable(_unstable);

      for (auto p : _vert_map) {
        auto it = o._vert_map.find(p.first);
        // Vertex exists in both
        if (it != o._vert_map.end()) {
          out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
          out_revmap.push_back(p.first);

          perm_x.push_back(p.second);
          perm_y.push_back((*it).second);
        }
      }

      // Build the permuted view of x and y.
      // assert(_g.size() > 0);
      GrPerm gx(perm_x, _g);
      // assert(o._g.size() > 0);
      GrPerm gy(perm_y, o._g);

      // Now perform the widening
      std::vector<_vert_id> destabilized;
      graph_t widen_g(GrOps::widen(gx, gy, destabilized));
      for (_vert_id v : destabilized) {
        widen_unstable.insert(v);
      }

      array_graph_t res(std::move(out_vmap), std::move(out_revmap),
                        std::move(widen_g), std::move(widen_unstable));
      CRAB_LOG("array-sgraph", crab::outs() << "Result widening:\n"
                                            << res << "\n";);
      return res;
    }
  }

  array_graph_t meet_or_narrowing(array_graph_t &o, bool is_meet,
                                  const std::string op) {

    if (is_bottom() || o.is_bottom())
      return bottom();
    else if (is_top())
      return o;
    else if (o.is_top())
      return *this;
    else {
      CRAB_LOG("array-sgraph", crab::outs() << "Before " << op << ":\n"
                                            << "Graph 1\n"
                                            << *this << "\n"
                                            << "Graph 2\n"
                                            << o << "\n");

      normalize();
      o.normalize();

      // Figure out the common renaming.
      std::vector<_vert_id> perm_x;
      std::vector<_vert_id> perm_y;

      vert_map_t out_vmap;
      rev_map_t out_revmap;

      for (auto p : _vert_map) {
        _vert_id vv = perm_x.size();
        out_vmap.insert(vmap_elt_t(p.first, vv));
        out_revmap.push_back(p.first);

        perm_x.push_back(p.second);
        perm_y.push_back(-1);
      }

      // Add missing mappings from the right operand.
      for (auto p : o._vert_map) {
        auto it = out_vmap.find(p.first);
        if (it == out_vmap.end()) {
          _vert_id vv = perm_y.size();
          out_revmap.push_back(p.first);

          perm_y.push_back(p.second);
          perm_x.push_back(-1);
          out_vmap.insert(vmap_elt_t(p.first, vv));
        } else {
          perm_y[(*it).second] = p.second;
        }
      }

      // Build the permuted view of x and y.
      GrPerm gx(perm_x, _g);
      GrPerm gy(perm_y, o._g);

      // Compute the syntactic meet/narrowing of the permuted graphs.
      std::vector<_vert_id> changes;
      graph_t out_g(GrOps::meet_or_narrowing(gx, gy, is_meet, changes));
      vert_set_t unstable;
      for (_vert_id v : changes)
        unstable.insert(v);

      GrOps::close_after_meet_or_narrowing(_g, vert_set_wrap_t(unstable));

      array_graph_t res(std::move(out_vmap), std::move(out_revmap),
                        std::move(out_g), vert_set_t());
      CRAB_LOG("array-sgraph", crab::outs() << "Result " << op << ":\n"
                                            << res << "\n";);
      return res;
    }
  }

  array_graph_t operator&(array_graph_t &o) {
    return meet_or_narrowing(o, true, "meet");
  }

  array_graph_t operator&&(array_graph_t &o) {
    return meet_or_narrowing(o, false, "narrowing");
  }

  // array_graph_t operator&(array_graph_t& o) {

  //   if (is_bottom() || o.is_bottom())
  //     return bottom();
  //   else if (is_top())
  //     return o;
  //   else if (o.is_top())
  //     return *this;
  //   else {
  //     CRAB_LOG("array-sgraph",
  //               crab::outs() << "Before meet:\n"<<"Graph 1\n"<<*this<<"\n"
  //                            <<"Graph 2\n"<<o << "\n");

  //     normalize();
  //     o.normalize();

  //     // Figure out the common renaming.
  //     std::vector<_vert_id> perm_x;
  //     std::vector<_vert_id> perm_y;

  //     vert_map_t out_vmap;
  //     rev_map_t out_revmap;

  //     for(auto p : _vert_map)
  //     {
  //       _vert_id vv = perm_x.size();
  //       out_vmap.insert(vmap_elt_t(p.first, vv));
  //       out_revmap.push_back(p.first);

  //       perm_x.push_back(p.second);
  //       perm_y.push_back(-1);
  //     }

  //     // Add missing mappings from the right operand.
  //     for(auto p : o._vert_map)
  //     {
  //       auto it = out_vmap.find(p.first);
  //       if(it == out_vmap.end())
  //       {
  //         _vert_id vv = perm_y.size();
  //         out_revmap.push_back(p.first);

  //         perm_y.push_back(p.second);
  //         perm_x.push_back(-1);
  //         out_vmap.insert(vmap_elt_t(p.first, vv));
  //       } else {
  //         perm_y[(*it).second] = p.second;
  //       }
  //     }

  //     // Build the permuted view of x and y.
  //     //assert(_g.size() > 0);
  //     GrPerm gx(perm_x, _g);
  //     //assert(o._g.size() > 0);
  //     GrPerm gy(perm_y, o._g);

  //     // Compute the syntactic meet of the permuted graphs.
  //     std::vector<_vert_id> changes;
  //     graph_t meet_g(GrOps::meet_or_narrowing(gx, gy, true /*meet*/,
  //     changes)); vert_set_t unstable; for(_vert_id v : changes)
  //       unstable.insert(v);

  //     GrOps::close_after_meet_or_narrowing(_g, vert_set_wrap_t(unstable));

  //     array_graph_t res(std::move(out_vmap), std::move(out_revmap),
  //                              std::move(meet_g), vert_set_t());
  //     CRAB_LOG("array-sgraph", crab::outs() << "Result meet:\n"<< res
  //     <<"\n";); return res;
  //   }
  // }

  // array_graph_t operator&&(array_graph_t& o) {
  //   if (is_bottom() || o.is_bottom())
  //     return bottom();
  //   else if (is_top())
  //     return o;
  //   else{
  //     CRAB_LOG("array-sgraph",
  //               crab::outs() << "Before narrowing:\n"<<"Graph
  //               1\n"<<*this<<"\n"
  //                            <<"Graph 2\n"<<o<<"\n";);

  //     // Narrowing as a no-op should be sound.
  //     normalize();
  //     array_graph_t res(*this);

  //     CRAB_LOG("array-sgraph",
  //               crab::outs() << "Result narrowing:\n" << res<<"\n";);
  //     return res;
  //   }
  // }

  void operator-=(Vertex v) {
    if (is_bottom())
      return;
    auto it = _vert_map.find(v);
    if (it != _vert_map.end()) {
      normalize();
      _g.forget(it->second);
      _rev_map[it->second] = boost::none;
      _vert_map.erase(v);
    }
  }

  void remove_from_weights(typename Weight::variable_t v) {
    mut_val_ref_t w_pq;
    for (auto p : _g.verts())
      for (auto e : _g.e_succs(p)) {
        auto q = e.vert;
        if (_g.lookup(p, q, &w_pq)) {
          Weight w = (Weight)w_pq;
          w -= v;
          Wt_meet op;
          _g.update_edge(p, w, q, op);
          GrOps::close_after_edge(_g, p, q);
        }
      }
  }

  void write(crab_os &o) { write(o, true); }

  void write(crab_os &o, bool print_bottom_edges) {

    normalize();

    if (is_bottom()) {
      o << "_|_";
      return;
    } else if (is_top()) {
      o << "{}";
      return;
    } else {
      bool first = true;
      o << "{";

      // sorting only for helping debugging
      std::vector<_vert_id> verts;
      for (auto v : _g.verts())
        verts.push_back(v);
      std::sort(verts.begin(), verts.end());

      for (_vert_id s : verts) {
        if (!_rev_map[s])
          continue;

        auto vs = *_rev_map[s];

        // sorting only for helping debugging
        std::vector<_vert_id> succs;
        for (auto d : _g.succs(s))
          succs.push_back(d);
        std::sort(succs.begin(), succs.end());

        for (_vert_id d : succs) {
          if (!_rev_map[d])
            continue;

          auto w = _g.edge_val(s, d);
          if (!print_bottom_edges && w.is_bottom())
            continue; // do not print bottom edges

          auto vd = *_rev_map[d];

          if (first)
            first = false;
          else
            o << ", ";
          o << "[" << vs << "," << vd << ")=>" << w;
        }
      }
      o << "}";
    }
  }
};

namespace array_graph_impl {

// JN: I do not know how to propagate arbitrary invariants
// between weight and scalar domains in a domain-independent
// manner. Here, we define propagations between specific
// domains.

template <typename Dom>
interval<typename Dom::number_t>
eval_interval(Dom dom, typename Dom::linear_expression_t e) {
  interval<typename Dom::number_t> r = e.constant();
  for (auto p : e)
    r += p.first * dom[p.second];
  return r;
}

template <typename Dom1, typename Dom2>
void propagate_between_weight_and_scalar(
    Dom1 src, typename Dom1::linear_expression_t src_e, variable_type ty,
    Dom2 &dst, typename Dom2::variable_t dst_var) {

  if (ty == ARR_INT_TYPE || ty == ARR_REAL_TYPE) {
    // --- XXX: simplification wrt Gange et.al.:
    //     Only non-relational numerical invariants are
    //     propagated from the graph domain to the scalar domain.
    dst.set(dst_var, eval_interval(src, src_e));
  } else {
    CRAB_WARN("Unsupported array type ", __LINE__, ":",
              "missing propagation between weight and scalar domains");
  }
}

template <typename BaseDom>
void propagate_between_weight_and_scalar(
    numerical_nullity_domain<BaseDom> src,
    typename BaseDom::linear_expression_t src_e, variable_type ty,
    numerical_nullity_domain<BaseDom> &dst,
    typename BaseDom::variable_t dst_var) {

  if (ty == ARR_INT_TYPE || ARR_REAL_TYPE) {
    // --- XXX: simplification wrt Gange et.al.:
    //     Only non-relational numerical invariants are
    //     propagated from the graph domain to the scalar domain.
    dst.set(dst_var, eval_interval(src, src_e));
  } else if (ty == ARR_PTR_TYPE) {
    if (auto src_var = src_e.get_variable()) {
      auto &null_dom = dst.second();
      null_dom.set_nullity(dst_var, src.get_nullity(*src_var));
    }
  } else {
    CRAB_WARN("Unsupported array type ", __LINE__, ":",
              "missing propagation between weight and scalar domains");
  }
}

template <typename BaseDom>
void propagate_between_weight_and_scalar(
    nullity_domain<typename BaseDom::number_t, typename BaseDom::varname_t> src,
    typename BaseDom::linear_expression_t src_e, variable_type ty,
    numerical_nullity_domain<BaseDom> &dst,
    typename BaseDom::variable_t dst_var) {

  if (ty == ARR_INT_TYPE || ty == ARR_REAL_TYPE) {
    // do nothing
  } else if (ty == ARR_PTR_TYPE) {
    if (auto src_var = src_e.get_variable()) {
      auto &null_dom = dst.second();
      null_dom.set_nullity(dst_var, src.get_nullity(*src_var));
    }
  } else {
    CRAB_WARN("Unsupported array type ", __LINE__, ":",
              "missing propagation between weight and scalar domains");
  }
}

template <typename BaseDom>
void propagate_between_weight_and_scalar(
    numerical_nullity_domain<BaseDom> src,
    typename BaseDom::linear_expression_t src_e, variable_type ty,
    nullity_domain<typename BaseDom::number_t, typename BaseDom::varname_t>
        &dst,
    typename BaseDom::variable_t dst_var) {

  if (ty == ARR_INT_TYPE || ty == ARR_REAL_TYPE) {
    // do nothing
  } else if (ty == ARR_PTR_TYPE) {
    if (auto src_var = src_e.get_variable()) {
      dst.set_nullity(dst_var, src.second().get_nullity(*src_var));
    }
  } else {
    CRAB_WARN("Unsupported array type ", __LINE__, ":",
              "missing propagation between weight and scalar domains");
  }
}

} // namespace array_graph_impl

#if 1
template <class Vertex, class Weight, bool IsDistWeight>
using array_graph = array_graph_<Vertex, Weight, IsDistWeight>;
#else
// Wrapper which uses shared references with copy-on-write.
template <class Vertex, class Weight, bool IsDistWeight>
class array_graph : public writeable {
public:
  typedef array_graph_<Vertex, Weight, IsDistWeight> array_graph_impl_t;
  typedef std::shared_ptr<array_graph_impl_t> array_graph_ref_t;
  typedef array_graph<Vertex, Weight, IsDistWeight> array_graph_t;

  typedef typename array_graph_impl_t::Wt Wt;
  typedef typename array_graph_impl_t::mut_val_ref_t mut_val_ref_t;
  typedef typename array_graph_impl_t::graph_t graph_t;
  typedef typename array_graph_impl_t::vert_id vert_id;
  typedef typename array_graph_impl_t::succ_range succ_range;
  typedef typename array_graph_impl_t::pred_range pred_range;
  typedef typename array_graph_impl_t::vert_range vert_range;
  typedef typename array_graph_impl_t::e_succ_range e_succ_range;
  typedef typename array_graph_impl_t::e_pred_range e_pred_range;

  array_graph(array_graph_ref_t _ref) : norm_ref(_ref) {}

  array_graph(array_graph_ref_t _base, array_graph_ref_t _norm)
      : base_ref(_base), norm_ref(_norm) {}

  array_graph_t create(array_graph_impl_t &&t) {
    return std::make_shared<array_graph_impl_t>(std::move(t));
  }

  array_graph_t create_base(array_graph_impl_t &&t) {
    array_graph_ref_t base = std::make_shared<array_graph_impl_t>(t);
    array_graph_ref_t norm = std::make_shared<array_graph_impl_t>(std::move(t));
    return array_graph_t(base, norm);
  }

  void lock(void) { // Allocate a fresh copy.
    if (!norm_ref.unique())
      norm_ref = std::make_shared<array_graph_impl_t>(*norm_ref);
    base_ref.reset();
  }

public:
  static array_graph_t top() { return array_graph(false); }

  static array_graph_t bottom() { return array_graph(true); }

  array_graph(bool is_bottom = false)
      : norm_ref(std::make_shared<array_graph_impl_t>(is_bottom)) {}

  array_graph(const array_graph_t &o)
      : base_ref(o.base_ref), norm_ref(o.norm_ref) {}

  array_graph_t &operator=(const array_graph_t &o) {
    base_ref = o.base_ref;
    norm_ref = o.norm_ref;
    return *this;
  }

  array_graph_impl_t &base(void) {
    if (base_ref)
      return *base_ref;
    else
      return *norm_ref;
  }

  array_graph_impl_t &norm(void) { return *norm_ref; }
  const array_graph_impl_t &norm(void) const { return *norm_ref; }

  bool is_bottom() { return norm().is_bottom(); }

  bool is_top() { return norm().is_top(); }

  bool operator<=(array_graph_t &o) { return norm() <= o.norm(); }

  void operator|=(array_graph_t o) {
    lock();
    norm() |= o.norm();
  }

  array_graph_t operator|(array_graph_t o) { return create(norm() | o.norm()); }

  array_graph_t operator||(array_graph_t o) {
    return create_base(base() || o.norm());
  }

  array_graph_t operator&(array_graph_t o) { return create(norm() & o.norm()); }

  array_graph_t operator&&(array_graph_t o) {
    return create(norm() && o.norm());
  }

  template <typename Thresholds>
  array_graph_t widening_thresholds(array_graph_t o, const Thresholds &ts) {
    return create_base(
        base().template widening_thresholds<Thresholds>(o.norm(), ts));
  }

  void normalize() {
    lock();
    norm().normalize();
  }

  vert_range verts() { return norm().verts(); }

  succ_range succs(Vertex v) { return norm().succs(v); }

  pred_range preds(Vertex v) { return norm().preds(v); }

  e_succ_range e_succs(Vertex v) { return norm().e_succs(v); }

  e_pred_range e_preds(Vertex v) { return norm().e_preds(v); }

  void set_to_bottom() {
    lock();
    norm().set_to_bottom();
  }

  bool lookup_edge(Vertex s, Vertex d, mut_val_ref_t *w) {
    lock();
    return norm().lookup_edge(s, d, w);
  }

  void expand(Vertex s, Vertex d) {
    lock();
    norm().expand(s, d);
  }

  void update_edge(Vertex s, Weight w, Vertex d) {
    lock();
    norm().update_edge(s, w, d);
  }

  // void full_close() { lock(); norm().full_close(); }
  // void update_edge_unclosed(Vertex s, Weight w, Vertex d) {
  // lock(); norm().update_edge_unclosed(s,w,d); }
  void close_edge(Vertex s, Vertex d) {
    lock();
    norm().close_edge(s, d);
  }

  void remove_from_weights(typename Weight::variable_t v) {
    lock();
    norm().remove_from_weights(v);
  }

  void operator-=(Vertex v) {
    lock();
    norm() -= v;
  }

  void write(crab_os &o) { norm().write(o); }
  void write(crab_os &o, bool print_bottom_edges) {
    norm().write(o, print_bottom_edges);
  }

protected:
  array_graph_ref_t base_ref;
  array_graph_ref_t norm_ref;
};
#endif

// Landmark: another C++ datatype to wrap variables and numbers as
// graph vertices.

enum landmark_kind_t { LMC, LMV, LMVP };

template <class Variable, class Number> class landmark {
protected:
  landmark_kind_t _kind;
  landmark(landmark_kind_t kind) : _kind(kind) {}

public:
  virtual ~landmark() {}

  landmark_kind_t kind() const { return _kind; }

  virtual bool operator==(const landmark<Variable, Number> &o) const = 0;

  virtual bool operator<(const landmark<Variable, Number> &o) const = 0;

  virtual void write(crab_os &o) const = 0;

  virtual std::size_t hash() const = 0;
};

template <class Variable, class Number>
class landmark_cst : public landmark<Variable, Number> {

  Number _n;
  typedef landmark<Variable, Number> landmark_t;
  typedef landmark_cst<Variable, Number> landmark_cst_t;

public:
  landmark_cst(Number n) : landmark_t(landmark_kind_t::LMC), _n(n) {}

  bool operator==(const landmark_t &o) const {
    if (this->_kind != o.kind())
      return false;

    assert(o.kind() == landmark_kind_t::LMC);
    auto o_ptr = static_cast<const landmark_cst_t *>(&o);
    return (_n == o_ptr->_n);
  }

  bool operator<(const landmark_t &o) const {
    if (this->_kind != o.kind())
      return true;

    assert(o.kind() == landmark_kind_t::LMC);
    return (_n < static_cast<const landmark_cst_t *>(&o)->_n);
  }

  void write(crab_os &o) const { o << _n; }

  std::size_t hash() const { return std::hash<Number>{}(_n); }

  Number get_cst() const { return _n; }
};

template <class Variable, class Number>
class landmark_var : public landmark<Variable, Number> {

  Variable _v;
  typedef landmark<Variable, Number> landmark_t;
  typedef landmark_var<Variable, Number> landmark_var_t;

public:
  landmark_var(Variable v) : landmark_t(landmark_kind_t::LMV), _v(v) {}

  bool operator==(const landmark_t &o) const {
    if (this->_kind != o.kind())
      return false;

    assert(o.kind() == landmark_kind_t::LMV);
    auto o_ptr = static_cast<const landmark_var_t *>(&o);
    return (_v == o_ptr->_v);
  }

  bool operator<(const landmark_t &o) const {
    if (this->_kind == o.kind()) {
      assert(o.kind() == landmark_kind_t::LMV);
      auto o_ptr = static_cast<const landmark_var_t *>(&o);
      return (_v < o_ptr->_v);
    } else if (o.kind() == LMC) {
      return false;
    } else if (o.kind() == LMVP) {
      return true;
    } else
      CRAB_ERROR("unreachable!");
  }

  void write(crab_os &o) const { o << _v; }

  std::size_t hash() const { return _v.hash(); }

  Variable get_var() const { return _v; }
};

template <class Variable, class Number>
class landmark_varprime : public landmark<Variable, Number> {
  std::string _lm;
  Variable _v;

  typedef landmark<Variable, Number> landmark_t;
  typedef landmark_varprime<Variable, Number> landmark_var_prime_t;

public:
  landmark_varprime(std::string lm, Variable v)
      : landmark_t(landmark_kind_t::LMVP), _lm(lm), _v(v) {}

  bool operator==(const landmark_t &o) const {
    if (this->_kind != o.kind())
      return false;

    assert(o.kind() == landmark_kind_t::LMVP);
    auto o_ptr = static_cast<const landmark_var_prime_t *>(&o);
    return (_v == o_ptr->_v);
  }

  bool operator<(const landmark_t &o) const {
    if (this->_kind != o.kind())
      return false;

    assert(o.kind() == landmark_kind_t::LMVP);
    auto o_ptr = static_cast<const landmark_var_prime_t *>(&o);
    return (_v < o_ptr->_v);
  }

  void write(crab_os &o) const { o << _lm << "'"; }

  std::size_t hash() const { return _v.hash(); }

  Variable get_var() const { return _v; }
};

// Wrapper for landmark
template <class Variable, class Number> class landmark_ref {
  typedef landmark<Variable, Number> landmark_t;
  typedef landmark_ref<Variable, Number> landmark_ref_t;

public:
  // yeah I know this is bad ...
  std::shared_ptr<landmark_t> _ref;

  landmark_ref(Variable v, std::string name = "") : _ref(nullptr) {
    if (name == "") {
      _ref = std::static_pointer_cast<landmark_t>(
          std::make_shared<landmark_var<Variable, Number>>(
              landmark_var<Variable, Number>(v)));
    } else {
      _ref = std::static_pointer_cast<landmark_t>(
          std::make_shared<landmark_varprime<Variable, Number>>(
              landmark_varprime<Variable, Number>(name, v)));
    }
  }
  landmark_ref(Number n)
      : _ref(std::static_pointer_cast<landmark_t>(
            std::make_shared<landmark_cst<Variable, Number>>(
                landmark_cst<Variable, Number>(n)))) {}

  landmark_kind_t kind() const { return _ref->kind(); }

  bool operator==(const landmark_ref &o) const { return (*_ref == *(o._ref)); }

  bool operator<(const landmark_ref &o) const { return (*_ref < *(o._ref)); }

  void write(crab_os &o) const { _ref->write(o); }

  std::size_t hash() const { return _ref->hash(); }
};

// super unsafe!
template <class Variable, class Number>
inline Variable get_var(const landmark_ref<Variable, Number> &lm) {
  assert(lm.kind() == LMVP);
  return std::static_pointer_cast<const landmark_varprime<Variable, Number>>(
             lm._ref)
      ->get_var();
}

// super unsafe!
template <class Variable, class Number>
inline Number get_cst(const landmark_ref<Variable, Number> &lm) {
  assert(lm.kind() == LMC);
  return std::static_pointer_cast<const landmark_cst<Variable, Number>>(lm._ref)
      ->get_cst();
}

template <class Variable, class Number>
inline crab_os &operator<<(crab_os &o,
                           const landmark_ref<Variable, Number> &lm) {
  lm.write(o);
  return o;
}
} // end namespace domains
} // end namespace crab

namespace std {
template <class Variable, class Number>
struct hash<crab::domains::landmark_ref<Variable, Number>> {
  using landmark_ref_t = crab::domains::landmark_ref<Variable, Number>;
  size_t operator()(const landmark_ref_t &lm) { return lm.hash(); }
};
} // namespace std

namespace crab {
namespace domains {

template <typename Variable, typename Number> struct landmark_ref_hasher {
  size_t operator()(const landmark_ref<Variable, Number> &lm) const {
    return lm.hash();
  }
};

template <typename Variable, typename Number> struct landmark_ref_equal {
  bool operator()(const landmark_ref<Variable, Number> &lm1,
                  const landmark_ref<Variable, Number> &lm2) const {
    return lm1 == lm2;
  }
};

template <typename Variable, class Number, class MappedVal>
using landmark_ref_unordered_map =
    std::unordered_map<landmark_ref<Variable, Number>, MappedVal,
                       landmark_ref_hasher<Variable, Number>,
                       landmark_ref_equal<Variable, Number>>;

/*
  Reduced product of a numerical domain with a weighted array
  graph.

  FIXME: the set of array landmarks are chosen statically before
  starting the analysis (do_initialization). However, at
  anytime only those alive (using scalar domain) are
  considered.

  The main issue is that the landmarks are kept as global
  state. This is really error-prone. For instance,
  landmarks are reset each time a new CFG is analyzed. For
  a summary-based inter-procedural analysis might be ok but
  not, e.g., for inlining.
*/
template <typename NumDom, typename Content, bool IsDistContent = false>
class array_graph_domain final
    : public abstract_domain<
          array_graph_domain<NumDom, Content, IsDistContent>> {
public:
  typedef typename NumDom::number_t number_t;
  typedef typename NumDom::varname_t varname_t;

private:
  static_assert(std::is_same<number_t, typename Content::number_t>::value,
                "Scalar and Content domains must have the same number type");
  static_assert(std::is_same<varname_t, typename Content::varname_t>::value,
                "Scalar and Content domains must have the same varname type");

  typedef array_graph_domain<NumDom, Content, IsDistContent>
      array_graph_domain_t;
  typedef abstract_domain<array_graph_domain_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;

  typedef typename NumDom::variable_t variable_t;
  typedef typename NumDom::variable_vector_t variable_vector_t;
  typedef interval<number_t> interval_t;

  typedef landmark_cst<variable_t, number_t> landmark_cst_t;
  typedef landmark_var<variable_t, number_t> landmark_var_t;
  typedef landmark_varprime<variable_t, number_t> landmark_var_prime_t;
  typedef landmark_ref<variable_t, number_t> landmark_ref_t;
  typedef array_graph<landmark_ref_t, Content, IsDistContent> array_graph_t;

  //// XXX: make this a template parameter later
  typedef crab::cfg::var_factory_impl::str_var_alloc_col::varname_t
      str_varname_t;
  typedef interval_domain<z_number, str_varname_t> str_interval_dom_t;
  typedef term::TDomInfo<z_number, varname_t, str_interval_dom_t> idom_info;
  typedef term_domain<idom_info> expression_domain_t;

private:
  typedef typename array_graph_t::mut_val_ref_t mut_val_ref_t;

  // Quick wrapper to perform efficient unsat queries on the
  // scalar domain.
  struct solver_wrapper {
    // XXX: do not pass by reference
    NumDom _inv;
    solver_wrapper(NumDom inv) : _inv(inv) {}
    bool is_unsat(linear_constraint_t cst) {
      // XXX: it might modify _inv so that's why we make a copy in
      // the constructor.
      return array_graph_domain_helper_traits<NumDom>::is_unsat(_inv, cst);
    }
  };

  NumDom _scalar;
  expression_domain_t
      _expressions; // map each program variable to a symbolic expression
  array_graph_t _g;

  // A landmark is either a variable or number that may appear as
  // an array index. In addition, for each landmark l we keep
  // track of a prime landmark l' whose meaning is l'=l+1.

  /// === Static data
  using lm_map_t =
      landmark_ref_unordered_map<variable_t, number_t, landmark_ref_t>;
  static lm_map_t s_var_landmarks;
  static lm_map_t s_cst_landmarks;

  // --- landmark iterators
  struct get_first : public std::unary_function<typename lm_map_t::value_type,
                                                landmark_ref_t> {
    get_first() {}
    landmark_ref_t operator()(const typename lm_map_t::value_type &p) const {
      return p.first;
    }
  };
  struct get_second : public std::unary_function<typename lm_map_t::value_type,
                                                 landmark_ref_t> {
    get_second() {}
    landmark_ref_t operator()(const typename lm_map_t::value_type &p) const {
      return p.second;
    }
  };
  typedef boost::transform_iterator<get_first, typename lm_map_t::iterator>
      lm_iterator;
  typedef boost::transform_iterator<get_second, typename lm_map_t::iterator>
      lm_prime_iterator;
  typedef boost::iterator_range<lm_iterator> lm_range;
  typedef boost::iterator_range<lm_prime_iterator> lm_prime_range;

  lm_prime_iterator var_lm_prime_begin() {
    return boost::make_transform_iterator(s_var_landmarks.begin(),
                                          get_second());
  }
  lm_prime_iterator var_lm_prime_end() {
    return boost::make_transform_iterator(s_var_landmarks.end(), get_second());
  }
  lm_prime_range var_lm_primes() {
    return boost::make_iterator_range(var_lm_prime_begin(), var_lm_prime_end());
  }

  lm_prime_iterator cst_lm_prime_begin() {
    return boost::make_transform_iterator(s_cst_landmarks.begin(),
                                          get_second());
  }
  lm_prime_iterator cst_lm_prime_end() {
    return boost::make_transform_iterator(s_cst_landmarks.end(), get_second());
  }
  lm_prime_range cst_lm_primes() {
    return boost::make_iterator_range(cst_lm_prime_begin(), cst_lm_prime_end());
  }

public:
  template <class CFG> static void do_initialization(CFG cfg) {

    typedef crab::analyzer::array_segmentation<CFG> array_segment_analysis_t;
    typedef typename array_segment_analysis_t::array_segment_domain_t
        array_segment_domain_t;
    typedef crab::analyzer::array_constant_segment_visitor<
        CFG, array_segment_domain_t>
        array_cst_segment_visitor_t;

    std::set<landmark_ref_t> lms;

    // add variables
    array_segment_analysis_t analysis(cfg);
    analysis.exec();
    auto var_indexes = analysis.get_variables(cfg.entry());

    if (var_indexes.begin() == var_indexes.end()) {
      CRAB_WARN("No variables found in the cfg. No array graph landmarks will "
                "be added\n");
      return;
    }
    lms.insert(var_indexes.begin(), var_indexes.end());

    // get variable factory
    auto &vfac = (*var_indexes.begin()).name().get_var_factory();

    // add constants
    // make sure 0 is always considered as an array index
    lms.insert(landmark_ref_t(number_t(0)));
    typename array_cst_segment_visitor_t::constant_set_t constants;
    for (auto &bb : boost::make_iterator_range(cfg.begin(), cfg.end())) {
      auto var_indexes = analysis.get_variables(bb.label());
      // XXX: use some heuristics to choose "relevant" constants
      array_cst_segment_visitor_t vis(var_indexes);
      for (auto &s : boost::make_iterator_range(bb.begin(), bb.end()))
        s.accept(&vis);
      auto cst_indexes = vis.get_constants();
      lms.insert(cst_indexes.begin(), cst_indexes.end());
    }

    set_landmarks(lms, vfac);
  }

  template <class Range, class VarFactory>
  static void set_landmarks(const Range &lms, VarFactory &vfac) {

    s_var_landmarks.clear();
    s_cst_landmarks.clear();

    unsigned num_vl = 0;
    unsigned num_cl = 0;

    for (auto lm : lms) {
      switch (lm.kind()) {
      case LMV: {
        auto v =
            std::static_pointer_cast<const landmark_var_t>(lm._ref)->get_var();
        variable_t v_prime(vfac.get(v.index()));
        landmark_ref_t lm_prime(v_prime, v.name().str());
        s_var_landmarks.insert(std::make_pair(lm, lm_prime));
        num_vl++;
        break;
      }
      case LMC: {
        auto n =
            std::static_pointer_cast<const landmark_cst_t>(lm._ref)->get_cst();
        variable_t v_prime(vfac.get());
        landmark_ref_t lm_prime(v_prime, n.get_str());
        s_cst_landmarks.insert(std::make_pair(lm, lm_prime));
        num_cl++;
        break;
      }
      default:
        CRAB_ERROR("A landmark can only be either variable or constant");
      }
    }
    CRAB_LOG(
        "array-sgraph-domain-landmark",
        crab::outs() << "Added " << num_vl << " variable landmarks "
                     << "and " << num_cl << " constant landmarks={";
        bool first = true; for (auto &l
                                : s_var_landmarks) {
          if (!first)
            crab::outs() << ",";
          first = false;
          crab::outs() << l.first;
        } for (auto &l
               : s_cst_landmarks) {
          if (!first)
            crab::outs() << ",";
          first = false;
          crab::outs() << l.first;
        } crab::outs() << "}\n";);
  }

public: // public only for tests
  void add_landmark(variable_t v) {
    landmark_ref_t lm_v(v);
    landmark_ref_t lm_v_prime(
        variable_t(v.name().get_var_factory().get(v.name().index())),
        v.name().str());
    // add pair  x -> x'
    s_var_landmarks.insert(std::make_pair(lm_v, lm_v_prime));
    // x' = x + 1
    _scalar += make_prime_relation(lm_v_prime, lm_v);

    // reduce between _scalar and the array graph
    if (!reduce(_scalar, _g)) {
      // FIXME: incremental version
      // TODO: we can assume that the scalar domain is zones so
      // that we can return the affected edges after each
      // operation and apply reduction only on those edges. That
      // would suffice for now. If the scalar domain is not zones
      // then we don't reduce incrementally.
      set_to_bottom();
    }

    CRAB_LOG("array-sgraph-domain-landmark",
             crab::outs() << "Added landmark " << v << "\n";);
  }

  void remove_landmark(variable_t v) {
    array_forget(v);
    forget_prime_var(v);
    s_var_landmarks.erase(landmark_ref_t(v));

    CRAB_LOG("array-sgraph-domain-landmark",
             crab::outs() << "Removed landmark " << v << "\n";);
  }

private:
  // By active we mean current variables that are kept track by
  // the scalar domain.
  void get_active_landmarks(NumDom &scalar,
                            std::vector<landmark_ref_t> &landmarks) const {
    landmarks.reserve(s_cst_landmarks.size());
    for (auto p : s_cst_landmarks) {
      landmarks.push_back(p.first);
      landmarks.push_back(p.second);
    }
    std::vector<variable_t> active_vars;
    array_graph_domain_helper_traits<NumDom>::active_variables(scalar,
                                                               active_vars);
    for (auto v : active_vars) {
      auto it = s_var_landmarks.find(landmark_ref_t(v));
      if (it != s_var_landmarks.end()) {
        landmarks.push_back(landmark_ref_t(v));
        landmarks.push_back(it->second);
      }
    }
  }

  linear_expression_t make_expr(landmark_ref_t x) {
    switch (x.kind()) {
    case LMC:
      return std::static_pointer_cast<landmark_cst_t>(x._ref)->get_cst();
    case LMV:
      return variable_t(
          std::static_pointer_cast<landmark_var_t>(x._ref)->get_var());
    case LMVP:
      return variable_t(
          std::static_pointer_cast<landmark_var_prime_t>(x._ref)->get_var());
    default:
      CRAB_ERROR("unreachable!");
    }
  }

  // make constraint x < y
  linear_constraint_t make_lt_cst(landmark_ref_t x, landmark_ref_t y) {
    return linear_constraint_t(make_expr(x) <= make_expr(y) - 1);
  }

  // make constraint x <= y
  linear_constraint_t make_leq_cst(landmark_ref_t x, landmark_ref_t y) {
    return linear_constraint_t(make_expr(x) <= make_expr(y));
  }

  // make constraint x == y
  linear_constraint_t make_eq_cst(landmark_ref_t x, landmark_ref_t y) {
    return linear_constraint_t(make_expr(x) == make_expr(y));
  }

  // make constraint x' == x+1
  linear_constraint_t make_prime_relation(landmark_ref_t x_prime,
                                          landmark_ref_t x) {
    return linear_constraint_t(make_expr(x_prime) == make_expr(x) + 1);
  }

  // return true if v is a landmark in the graph
  bool is_landmark(variable_t v) const {
    landmark_ref_t lm_v(v);
    auto it = s_var_landmarks.find(lm_v);
    return (it != s_var_landmarks.end());
  }

  // return true if n is a landmark in the graph
  bool is_landmark(z_number n) const {
    landmark_ref_t lm_n(n);
    auto it = s_cst_landmarks.find(lm_n);
    return (it != s_cst_landmarks.end());
  }

  // return the prime landmark of v
  landmark_ref_t get_landmark_prime(variable_t v) const {
    landmark_ref_t lm_v(v);
    auto it = s_var_landmarks.find(lm_v);
    assert(it != s_var_landmarks.end());
    return it->second;
  }

  // Return the weight from the edge (i, i') otherwise top
  Content array_edge(variable_t i) {
    if (is_bottom())
      return Content::bottom();
    if (is_top() || !is_landmark(i))
      return Content::top();

    mut_val_ref_t wi;
    if (_g.lookup_edge(landmark_ref_t(i), get_landmark_prime(i), &wi))
      return (Content)wi;
    else
      return Content::top();
  }

  // Remove v from the edge (i,i')
  void array_edge_forget(variable_t i, variable_t v) {
    if (is_bottom())
      return;

    if (!is_landmark(i))
      return;

    mut_val_ref_t wi;
    landmark_ref_t lm_i(i);
    landmark_ref_t lm_i_prime = get_landmark_prime(i);
    if (_g.lookup_edge(lm_i, lm_i_prime, &wi)) {
      Content w = (Content)wi;
      w -= v;
      // XXX: update_edge closes the array graph
      _g.update_edge(lm_i, w, lm_i_prime);
    }
  }

  // Remove v from all vertices and edges
  void array_forget(variable_t v) {
    if (is_bottom())
      return;
    if (!is_landmark(v))
      return;

    _g -= landmark_ref_t(v);
    _g.remove_from_weights(v);
  }

  // Update the weight from the edge (i, i')
  void array_edge_update(variable_t i, Content w) {
    if (is_bottom())
      return;

    //--- strong update
    if (!is_landmark(i))
      return;

    landmark_ref_t lm_i(i);
    landmark_ref_t lm_i_prime = get_landmark_prime(i);

    _g.update_edge(lm_i, w, lm_i_prime);
    mut_val_ref_t wi;
    if (!_g.lookup_edge(lm_i, lm_i_prime, &wi))
      return;

    //--- weak update
    // An edge (p,q) must be weakened if p <= i <= q and p < q
    solver_wrapper solve(_scalar);
    mut_val_ref_t w_pq;
    for (auto p : _g.verts()) {
      for (auto e : _g.e_succs(p)) {
        auto q = e.vert;
        if ((p == lm_i) && (q == lm_i_prime))
          continue;
        if (_g.lookup_edge(p, q, &w_pq) && ((Content)w_pq).is_bottom())
          continue;
        // we know already that p < q in the array graph

        // check p <= i
        if (solve.is_unsat(make_leq_cst(p, lm_i)))
          continue;
        // check i' <= q
        if (solve.is_unsat(make_leq_cst(lm_i_prime, q)))
          continue;

        w_pq = (Content)w_pq | (Content)wi;
      }
    }
  }

  // x := x op k
  template <typename VarOrNum>
  void apply_one_variable(operation_t op, variable_t x, VarOrNum k) {
    if (is_bottom())
      return;

    if (!is_landmark(x)) {
      // If x is not a landmark we just apply the operation on the
      // scalar domain and return.
      apply_only_scalar(op, x, x, k);
      return;
    }

    landmark_ref_t lm_x(x);
    landmark_ref_t lm_x_prime = get_landmark_prime(x);

    /// --- Add x_old and x_old' to store old values of x and x'

    variable_t x_old(x.name().get_var_factory().get());
    variable_t x_old_prime(x.name().get_var_factory().get());
    landmark_ref_t lm_x_old(x_old);
    landmark_ref_t lm_x_old_prime(x_old_prime, x_old.name().str());
    s_var_landmarks.insert(std::make_pair(lm_x_old, lm_x_old_prime));
    // x_old = x
    _scalar.assign(x_old, x);
    // relation between x_old and x'
    _scalar += make_prime_relation(lm_x_old_prime, lm_x_old);
    //_scalar += make_eq_cst(lm_x_old_prime, lm_x_prime);

    /*** Incremental graph reduction ***/
    //// x_old  has all the x predecessors and successors
    _g.expand(lm_x, lm_x_old);
    //// x_old' has all the x' predecessors and successors
    _g.expand(lm_x_prime, lm_x_old_prime);
    //// edges between x and x_old
    _g.update_edge(lm_x, Content::bottom(), lm_x_old);
    _g.update_edge(lm_x_old, Content::bottom(), lm_x);
    //// edges between x' and x_old'
    _g.update_edge(lm_x_prime, Content::bottom(), lm_x_old_prime);
    _g.update_edge(lm_x_old_prime, Content::bottom(), lm_x_prime);
    //// edges between x_old and x_old'
    mut_val_ref_t w;
    if (_g.lookup_edge(lm_x, lm_x_prime, &w))
      _g.update_edge(lm_x_old, (Content)w, lm_x_old_prime);
    _g.update_edge(lm_x_old_prime, Content::bottom(), lm_x_old);

    /// --- Remove x and x'
    _g -= lm_x;
    _g -= lm_x_prime;

    /// --- Perform operation in the scalar domain
    _scalar.apply(op, x, x, k);

    // restore relation between x and x'
    _scalar.apply(OP_ADDITION, get_var(lm_x_prime), x, 1);
    //_scalar -= get_var(lm_x_prime);
    //_scalar += make_prime_relation(lm_x_prime, lm_x);

    if (!reduce(_scalar, _g)) { // FIXME: incremental version
      set_to_bottom();
      return;
    }

    /// --- Remove x_old and x_old'
    _g -= lm_x_old;
    _g -= lm_x_old_prime;
    _scalar -= x_old;
    _scalar -= x_old_prime;
    s_var_landmarks.erase(lm_x_old);
  }

  // remove v' from scalar and array graph
  void forget_prime_var(variable_t v) {
    if (!is_landmark(v))
      return;

    landmark_ref_t lm_v_prime = get_landmark_prime(v);
    _scalar -= get_var(lm_v_prime);
    // XXX: v' cannot appear in the array weights so we do not
    //      need to call array_forget.
    _g -= lm_v_prime;
  }

  // perform the operation in the scalar domain assuming that
  // nothing can be done in the graph domain.
  template <class Op, class K>
  void apply_only_scalar(Op op, variable_t x, variable_t y, K k) {
    _scalar.apply(op, x, y, k);

    // Abstract x in the array graph
    if (is_landmark(x)) {
      array_forget(x);     // remove x from the array graph
      forget_prime_var(x); // remove x' from scalar and array graph
      /// XXX: I think no need to reduce here
    }
  }

  // return a pair with the normalized offset and a bool that is
  // true if a new landmark was added in the array graph
  std::pair<variable_t, bool> normalize_offset(variable_t o, z_number n) {
    CRAB_LOG("array-sgraph-domain-norm",
             crab::outs() << "BEFORE NORMALIZE OFFSET: expressions="
                          << _expressions << "\n");

    // --- create a fresh variable no such that no := o;
    variable_t no(o.name().get_var_factory().get());
    _expressions.assign(no, o);

    // -- apply no := no / n; in the expressions domain
    _expressions.apply(operation_t::OP_SDIV, no, no, n);

    // -- simplify the expression domain
    bool simp_done = _expressions.simplify(no);

    CRAB_LOG("array-sgraph-domain-norm",
             crab::outs() << "AFTER NORMALIZE OFFSET: expressions="
                          << _expressions << "\n");

    if (!simp_done) {
      CRAB_LOG(
          "array-sgraph-domain-norm",
          crab::outs()
              << "NO NORMALIZATION done using the expression abstraction\n");

      // cleanup of the expression abstraction
      _expressions -= no;

      bool added_lm = false;
      if (!is_landmark(o)) {
        add_landmark(o);
        added_lm = true;
      }

      return std::make_pair(o, added_lm);
    }

    CRAB_LOG("array-sgraph-domain-norm",
             crab::outs()
                 << "NORMALIZATION DONE! using the expression abstraction\n");

    // -- propagate equalities from _expressions to _scalar
    linear_constraint_system_t e_csts;
    reduced_domain_traits<expression_domain_t>::extract(
        _expressions, no, e_csts, /*only_equalities=*/false);
    _scalar += e_csts;

    // -- add landmark for the new array index
    add_landmark(no);

    // cleanup of the expression abstraction
    _expressions -= no;

    return std::make_pair(no, true);
  }

  // T1 and T2 are either variable_t or z_number
  template <typename T1, typename T2>
  void array_init(variable_t arr, T1 src, T2 dst, linear_expression_t val) {
    if (!is_landmark(src)) {
      crab::outs() << "WARNING no landmark found for " << src << "\n";
      return;
    }

    if (!is_landmark(dst)) {
      crab::outs() << "WARNING no landmark found for " << dst << "\n";
      return;
    }

    landmark_ref_t lm_src(src);
    landmark_ref_t lm_dst(dst);

    Content w = Content::top();
    array_graph_impl::propagate_between_weight_and_scalar(
        _scalar, val, arr.get_type(), w, arr);
    _g.update_edge(lm_src, w, lm_dst);
  }

  interval_t to_interval(linear_expression_t expr) {
    interval_t r(expr.constant());
    for (typename linear_expression_t::iterator it = expr.begin();
         it != expr.end(); ++it) {
      interval_t c(it->first);
      r += c * _scalar[it->second];
    }
    return r;
  }

public:
  static void clear_global_state() {
    s_var_landmarks.clear();
    s_cst_landmarks.clear();
  }

  void set_to_top() {
    array_graph_domain_t abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    _scalar = NumDom::bottom();
    _expressions = expression_domain_t::bottom();
    _g.set_to_bottom();
  }

  // void set_to_bottom() {
  //   array_graph_domain_t abs(true);
  // 	std::swap(*this, abs);
  // }

  array_graph_domain(bool is_bottom = false)
      : _scalar(NumDom::top()), _expressions(expression_domain_t::top()),
        _g(array_graph_t::top()) {
    if (is_bottom)
      set_to_bottom();
  }

  array_graph_domain(const NumDom &s, const expression_domain_t &e,
                     const array_graph_t &g)
      : _scalar(s), _expressions(e), _g(g) {
    if (_scalar.is_bottom() || _expressions.is_bottom() || _g.is_bottom())
      set_to_bottom();
  }

  array_graph_domain(NumDom &&s, expression_domain_t &&e, array_graph_t &&g)
      : _scalar(std::move(s)), _expressions(std::move(e)), _g(std::move(g)) {
    if (_scalar.is_bottom() || _expressions.is_bottom() || _g.is_bottom())
      set_to_bottom();
  }

  array_graph_domain(const array_graph_domain_t &o)
      : _scalar(o._scalar), _expressions(o._expressions), _g(o._g) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
  }

  array_graph_domain(array_graph_domain_t &&o)
      : _scalar(std::move(o._scalar)), _expressions(std::move(o._expressions)),
        _g(std::move(o._g)) {}

  array_graph_domain_t &operator=(const array_graph_domain_t &o) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");
    if (this != &o) {
      _scalar = o._scalar;
      _expressions = o._expressions;
      _g = o._g;
    }
    return *this;
  }

  array_graph_domain_t &operator=(array_graph_domain_t &&o) {
    _scalar = std::move(o._scalar);
    _expressions = std::move(o._expressions);
    _g = std::move(o._g);
    return *this;
  }

  bool is_top() {
    return _scalar.is_top() && _expressions.is_top() && _g.is_top();
  }

  bool is_bottom() {
    return _scalar.is_bottom() || _expressions.is_bottom() || _g.is_bottom();
  }

  bool operator<=(array_graph_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");

    CRAB_LOG("array-sgraph-domain", crab::outs() << "Leq " << *this << " and\n"
                                                 << o << "=\n";);
    bool res = (_scalar <= o._scalar) && (_expressions <= o._expressions) &&
               (_g <= o._g);
    CRAB_LOG("array-sgraph-domain", crab::outs() << res << "\n";);
    return res;
  }

  void operator|=(array_graph_domain_t o) { *this = (*this | o); }

  array_graph_domain_t operator|(array_graph_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");

    CRAB_LOG("array-sgraph-domain",
             crab::outs() << "Join " << *this << " and " << o << "=\n");
    array_graph_domain_t join(_scalar | o._scalar,
                              _expressions | o._expressions, _g | o._g);
    CRAB_LOG("array-sgraph-domain", crab::outs() << join << "\n";);
    return join;
  }

  array_graph_domain_t
  widening_thresholds(array_graph_domain_t o,
                      const iterators::thresholds<number_t> &ts) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");

    CRAB_LOG("array-sgraph-domain", crab::outs()
                                        << "Widening (w/ thresholds) " << *this
                                        << " and " << o << "=\n";);
    auto widen_scalar(_scalar.widening_thresholds(o._scalar, ts));
    auto widen_expr(_expressions.widening_thresholds(o._expressions, ts));
    auto widen_g(_g.widening_thresholds(o._g, ts));
    if (!reduce(widen_scalar, widen_g)) {
      CRAB_LOG("array-sgraph-domain", crab::outs() << "_|_\n";);
      return array_graph_domain_t::bottom();
    } else {
      array_graph_domain_t widen(widen_scalar, widen_expr, widen_g);
      CRAB_LOG("array-sgraph-domain", crab::outs() << widen << "\n";);
      return widen;
    }
  }

  array_graph_domain_t operator||(array_graph_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");

    CRAB_LOG("array-sgraph-domain",
             crab::outs() << "Widening " << *this << " and " << o << "=\n");
    auto widen_scalar(_scalar || o._scalar);
    auto widen_expr(_expressions || o._expressions);
    auto widen_g(_g || o._g);
    if (!reduce(widen_scalar, widen_g)) {
      CRAB_LOG("array-sgraph-domain", crab::outs() << "_|_\n";);
      return array_graph_domain_t::bottom();
    } else {
      array_graph_domain_t widen(widen_scalar, widen_expr, widen_g);
      CRAB_LOG("array-sgraph-domain", crab::outs() << widen << "\n";);
      return widen;
    }
  }

  array_graph_domain_t operator&(array_graph_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");

    CRAB_LOG("array-sgraph-domain",
             crab::outs() << "Meet " << *this << " and " << o << "=\n");
    auto meet_scalar(_scalar & o._scalar);
    auto meet_expr(_expressions & o._expressions);
    auto meet_g(_g & o._g);
    if (!reduce(meet_scalar, meet_g)) {
      CRAB_LOG("array-sgraph-domain", crab::outs() << "_|_\n";);
      return array_graph_domain_t::bottom();
    } else {
      array_graph_domain_t meet(meet_scalar, meet_expr, meet_g);
      CRAB_LOG("array-sgraph-domain", crab::outs() << meet << "\n";);
      return meet;
    }
  }

  array_graph_domain_t operator&&(array_graph_domain_t o) {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

    CRAB_LOG("array-sgraph-domain",
             crab::outs() << "Narrowing " << *this << " and " << o << "=\n");
    auto narrow_scalar(_scalar && o._scalar);
    auto narrow_expr(_expressions && o._expressions);
    auto narrow_g(_g && o._g);
    if (!reduce(narrow_scalar, narrow_g)) {
      CRAB_LOG("array-sgraph-domain", crab::outs() << "_|_\n";);
      return array_graph_domain_t::bottom();
    } else {
      array_graph_domain_t narrow(narrow_scalar, narrow_expr, narrow_g);
      CRAB_LOG("array-sgraph-domain", crab::outs() << narrow << "\n";);
      return narrow;
    }
  }

  // The reduction consists of detecting dead segments so it is
  // done only in one direction (scalar -> array graph). Note that
  // whenever an edge becomes bottom closure is also happening.
  // Return false if bottom is detected during the reduction.
  bool reduce(NumDom &scalar, array_graph_t &g) {
    crab::CrabStats::count(getDomainName() + ".count.reduce");
    crab::ScopedCrabStats __st__(getDomainName() + ".reduce");

    scalar.normalize();
    g.normalize();

    if (scalar.is_bottom() || g.is_bottom())
      return false;

    if (!scalar.is_top()) {
      std::vector<landmark_ref_t> active_landmarks;
      get_active_landmarks(scalar, active_landmarks);
      solver_wrapper solve(scalar);
      for (auto lm_s : active_landmarks)
        for (auto lm_d : active_landmarks) {
          // XXX: we do not exploit the following facts:
          //   - i < i' is always sat
          //   - i' < i is always unsat
          //   - if i < j  unsat then i' < j unsat.
          //   - if i < j' unsat then i' < j unsat.
          if ((lm_s == lm_d) || solve.is_unsat(make_lt_cst(lm_s, lm_d))) {
            g.update_edge(lm_s, Content::bottom(), lm_d);
          }
        }
    }
    return (!g.is_bottom());
  }

  void operator-=(variable_t v) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");

    if (is_bottom())
      return;

    // remove v from scalar and array graph
    _scalar -= v;
    // remove v from expressions
    _expressions -= v;

    if (is_landmark(v)) {
      array_forget(v);
      // remove v' from scalar and array graph
      forget_prime_var(v);
    }
  }

  void project(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.project");
    crab::ScopedCrabStats __st__(getDomainName() + ".project");

    if (is_bottom() || is_top()) {
      return;
    }
    if (variables.empty()) {
      set_to_top();
      return;
    }

    std::set<variable_t> keep_vars(variables.begin(), variables.end());
    std::vector<variable_t> active_vars;
    array_graph_domain_helper_traits<NumDom>::active_variables(_scalar,
                                                               active_vars);
    for (auto v : active_vars) {
      if (!keep_vars.count(v)) {
        array_forget(v);
        forget_prime_var(v);
      }
    }

    _scalar.project(variables);
    _expressions.project(variables);
  }

  void forget(const variable_vector_t &variables) {
    if (is_bottom() || is_top()) {
      return;
    }
    for (variable_t v : variables) {
      this->operator-=(v);
    }
  }

  void expand(variable_t var, variable_t new_var) {
    CRAB_WARN("array_graph_domain expand not implemented");
  }

  void normalize() {
    CRAB_WARN("array_graph_domain normalize not implemented");
  }

  void minimize() {
    _scalar.minimize();
    _expressions.minimize();
  }

  /* begin intrinsics operations */    
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  array_graph_domain_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    assert(from.size() == to.size());

    if (is_top() || is_bottom())
      return;

    CRAB_WARN(getDomainName(), "::rename not implemented");
  }

  void operator+=(linear_constraint_system_t csts) {
    crab::CrabStats::count(getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

    if (is_bottom())
      return;

    _scalar += csts;
    _expressions += csts;

    if (!reduce(_scalar, _g)) { // FIXME: incremental version
      set_to_bottom();
      return;
    }
    CRAB_LOG("array-sgraph-domain",
             crab::outs() << "Assume(" << csts << ") --- " << *this << "\n";);
  }

  void assign(variable_t x, linear_expression_t e) { assign(x, e, true); }

  // Perform the operation in the scalar (optionally expression)
  // domain and reduce.
  //
  // NOTE: if the assignment is something like i = i + k then we
  // will lose precision in the array graph. This kind of
  // assignments should be managed by the apply methods instead.
  void assign(variable_t x, linear_expression_t e, bool update_expressions) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    if (is_bottom())
      return;

    if (auto y = e.get_variable()) {
      // skip x:=x
      if ((*y) == x)
        return;
    }

    _scalar.assign(x, e);
    if (update_expressions)
      _expressions.assign(x, e);

    if (is_landmark(x)) {
      array_forget(x);
      // remove x' from scalar and array graph
      forget_prime_var(x);
      // restore the relationship between x and x'
      _scalar.apply(OP_ADDITION, get_var(get_landmark_prime(x)), x, 1);
      // XXX: is it needed ??
      //_g.close_edge(landmark_ref_t(x), get_landmark_prime(x));
    }

    if (!reduce(_scalar, _g)) { // FIXME: incremental version
      set_to_bottom();
      return;
    }

    CRAB_LOG("array-sgraph-domain", crab::outs()
                                        << "Assign " << x << " := " << e
                                        << " ==> " << *this << "\n";);
  }

  void apply(operation_t op, variable_t x, variable_t y, number_t z) {
    if (x == y) {
      crab::CrabStats::count(getDomainName() + ".count.apply");
      crab::ScopedCrabStats __st__(getDomainName() + ".apply");
      _expressions.apply(op, x, y, z);
      apply_one_variable<number_t>(op, x, z);
      CRAB_LOG("array-sgraph-domain",
               crab::outs() << "Apply " << x << " := " << y << " " << op << " "
                            << z << " ==> " << *this << "\n";);
    } else {
      switch (op) {
      case OP_ADDITION:
        assign(x, y + z);
        break;
      case OP_SUBTRACTION:
        assign(x, y - z);
        break;
      case OP_MULTIPLICATION:
        assign(x, y * z);
        break;
      default:
        CRAB_WARN(op, "not implemented in array-sgraph-domain\n");
      }
    }
  }

  void apply(operation_t op, variable_t x, variable_t y, variable_t z) {
    if (x == y) {
      crab::CrabStats::count(getDomainName() + ".count.apply");
      crab::ScopedCrabStats __st__(getDomainName() + ".apply");
      _expressions.apply(op, x, y, z);
      apply_one_variable<variable_t>(op, x, z);
      CRAB_LOG("array-sgraph-domain",
               crab::outs() << "Apply " << x << " := " << y << " " << op << " "
                            << z << " ==> " << *this << "\n";);
    } else {
      switch (op) {
      case OP_ADDITION:
        assign(x, y + z);
        break;
      case OP_SUBTRACTION:
        assign(x, y - z);
        break;
      default:
        CRAB_WARN(op, "not implemented in array-sgraph-domain");
      }
    }
  }
  // backward arithmetic operations

  void backward_assign(variable_t x, linear_expression_t e,
                       array_graph_domain_t invariant) {
    operator-=(x);
    CRAB_WARN("backward assign not implemented");
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      array_graph_domain_t invariant) {
    operator-=(x);
    CRAB_WARN("backward apply not implemented");
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      array_graph_domain_t invariant) {

    operator-=(x);
    CRAB_WARN("backward apply not implemented");
  }

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {
    operator-=(lhs);
    CRAB_WARN("assign_bool_cst not implemented in array-sgraph-domain");
  }

  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {
    operator-=(lhs);
    CRAB_WARN("assign_bool_var not implemented in array-sgraph-domain");
  }

  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {
    operator-=(x);
    CRAB_WARN("apply_binary_bool not implemented in array-sgraph-domain");
  }

  void assume_bool(variable_t v, bool is_negated) {
    CRAB_WARN("assume_bool not implemented in array-sgraph-domain");
  }

  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                array_graph_domain_t invariant) {
    operator-=(lhs);
    CRAB_WARN(
        "backward_assign_bool_cst not implemented in array-sgraph-domain");
  }

  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                array_graph_domain_t invariant) {
    operator-=(lhs);
    CRAB_WARN(
        "backward_assign_bool_var not implemented in array-sgraph-domain");
  }

  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z,
                                  array_graph_domain_t invariant) {
    operator-=(x);
    CRAB_WARN(
        "backward_apply_binary_bool not implemented in array-sgraph-domain");
  }

  // cast operations

  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    _expressions.apply(op, dst, src);
    // assume unlimited precision so widths are ignored.
    assign(dst, src, false);
  }

  // bitwise operations

  void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    _expressions.apply(op, x, y, z);
    // XXX: we give up soundly in the graph domain
    apply_only_scalar(op, x, y, z);
  }

  void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    _expressions.apply(op, x, y, k);
    // XXX: we give up soundly in the graph domain
    apply_only_scalar(op, x, y, k);
  }

  interval_t operator[](variable_t v) { return _scalar[v]; }

  // pointer_operators_api
  virtual void pointer_load(variable_t lhs, variable_t rhs) override {
    _scalar.pointer_load(lhs, rhs);
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs) override {
    _scalar.pointer_store(lhs, rhs);
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs,
                              linear_expression_t offset) override {
    _scalar.pointer_assign(lhs, rhs, offset);
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    _scalar.pointer_mk_obj(lhs, address);
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    _scalar.pointer_function(lhs, func);
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    _scalar.pointer_mk_null(lhs);
  }

  virtual void pointer_assume(pointer_constraint_t cst) override {
    _scalar.pointer_assume(cst);
  }

  virtual void pointer_assert(pointer_constraint_t cst) override {
    _scalar.pointer_assert(cst);
  }

  // array_operators_api

  virtual void array_init(variable_t a, linear_expression_t /*elem_size*/,
                          linear_expression_t lb_idx,
                          linear_expression_t ub_idx,
                          linear_expression_t val) override {

    auto lb_var_opt = lb_idx.get_variable();
    auto ub_var_opt = ub_idx.get_variable();

    if (lb_idx.is_constant() && ub_idx.is_constant())
      array_init(a, lb_idx.constant(), ub_idx.constant(), val);
    else if (lb_idx.is_constant() && ub_var_opt)
      array_init(a, lb_idx.constant(), *ub_var_opt, val);
    else if (lb_var_opt && ub_idx.is_constant())
      array_init(a, *lb_var_opt, ub_idx.constant(), val);
    else if (lb_var_opt && ub_var_opt)
      array_init(a, *lb_var_opt, *ub_var_opt, val);
    else
      CRAB_ERROR("unreachable");
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t elem_size,
                          linear_expression_t i) override {
    crab::CrabStats::count(getDomainName() + ".count.load");
    crab::ScopedCrabStats __st__(getDomainName() + ".load");

    auto vi = i.get_variable();
    if (!vi) {
      CRAB_WARN("TODO: array load index must be a variable");
      return;
    }

    // -- normalization ensures that closure and reduction have
    // -- been applied.
    interval_t i_elem_size = to_interval(elem_size);
    boost::optional<number_t> n_elem_size = i_elem_size.singleton();
    if (!n_elem_size) {
      CRAB_WARN("array_graph ignored array load because element size is not "
                "constant");
      return;
    }
    assert(static_cast<int64_t>(*n_elem_size) > 0 &&
           static_cast<int64_t>(*n_elem_size) <=
               std::numeric_limits<uint64_t>::max());
    int64_t num_bytes = static_cast<int64_t>(*n_elem_size);
    auto p = normalize_offset(*vi, num_bytes);
    variable_t norm_idx = p.first;

    // #if 0
    // if (is_landmark(norm_idx)) {
    //   landmark_ref_t lm_norm_idx(norm_idx);
    //   landmark_ref_t lm_norm_idx_prime = get_landmark_prime(norm_idx);

    //   _g.close_edge(lm_norm_idx, lm_norm_idx_prime);
    //   crab::outs() << "#### 1 " << _g << "\n";

    //   Content w;
    //   w += linear_constraint_t(linear_expression_t(lhs) ==
    //   linear_expression_t(a)); crab::outs() << "#### 2 " << w << "\n";
    //   //_g.update_edge_unclosed(lm_norm_idx, w, lm_norm_idx_prime);
    //   _g.update_edge(lm_norm_idx, w, lm_norm_idx_prime);
    // }
    // #endif

    Content w = array_edge(norm_idx);

    if (a.get_type() == ARR_INT_TYPE || a.get_type() == ARR_REAL_TYPE) {
      // Only non-relational numerical invariants are
      // propagated from the graph domain to the expressions domain.
      _expressions.set(lhs, w[a]);
    }

    array_graph_impl::propagate_between_weight_and_scalar(w, a, a.get_type(),
                                                          _scalar, lhs);

    // if normalize_offset created a landmark we remove it here to
    // keep smaller array graph
    if (p.second)
      remove_landmark(norm_idx);

    /// XXX: due to the above simplification we need to reduce
    /// only if the content of an array cell can be an index.
    if (is_landmark(lhs))
      if (!reduce(_scalar, _g)) // FIXME: incremental version
        set_to_bottom();

    CRAB_LOG("array-sgraph-domain", crab::outs() << "Array read " << lhs
                                                 << " := " << a << "[" << i
                                                 << "] ==> " << *this << "\n";);
  }

  virtual void array_store(variable_t a, linear_expression_t elem_size,
                           linear_expression_t i, linear_expression_t val,
                           bool /*is_strong_update*/) override {
    crab::CrabStats::count(getDomainName() + ".count.store");
    crab::ScopedCrabStats __st__(getDomainName() + ".store");

    auto vi = i.get_variable();
    if (!vi) {
      CRAB_WARN("TODO: array store index must be a variable");
      return;
    }

    Content w = Content::top();
    array_graph_impl::propagate_between_weight_and_scalar(_scalar, val,
                                                          a.get_type(), w, a);

    interval_t i_elem_size = to_interval(elem_size);
    boost::optional<number_t> n_elem_size = i_elem_size.singleton();
    if (!n_elem_size) {
      CRAB_WARN("array_graph ignored array store because element size is not "
                "constant");
      return;
    }
    assert(static_cast<int64_t>(*n_elem_size) > 0 &&
           static_cast<int64_t>(*n_elem_size) <=
               std::numeric_limits<uint64_t>::max());
    int64_t num_bytes = static_cast<int64_t>(*n_elem_size);
    auto p = normalize_offset(*vi, num_bytes);
    variable_t norm_idx = p.first;

    array_edge_forget(norm_idx, a);
    array_edge_update(norm_idx, w);

    // XXX: since we do not propagate from the array weights to
    // the scalar domain I think we don't need to reduce here.

    CRAB_LOG("array-sgraph-domain", crab::outs() << "Array write " << a << "["
                                                 << i << "] := " << val
                                                 << " ==> " << *this << "\n";);
  }

  virtual void array_store(variable_t a_new, variable_t a_old,
                           linear_expression_t elem_size, linear_expression_t i,
                           linear_expression_t val,
                           bool /*is_strong_update*/) override {
    CRAB_WARN("array_store in array_graph not implemented");
  }

  virtual void array_store_range(variable_t a, linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    CRAB_WARN("array_store_range in array_graph not implemented");
  }

  virtual void array_store_range(variable_t a_new, variable_t a_old,
                                 linear_expression_t elem_size,
                                 linear_expression_t i, linear_expression_t j,
                                 linear_expression_t val) override {
    CRAB_WARN("array_store_range in array_graph not implemented");
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    CRAB_WARN("array_assign in array_graph not implemented");
  }

  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           array_graph_domain_t invariant) {
    CRAB_WARN("backward_array_init in array_graph domain not implemented");
  }
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           array_graph_domain_t invariant) {
    CRAB_WARN("backward_array_load in array_graph domain not implemented");
    this->operator-=(lhs);
  }
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            array_graph_domain_t invariant) {
    CRAB_WARN("backward_array_store in array_graph domain not implemented");
  }
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update,
                            array_graph_domain_t invariant) {
    CRAB_WARN("backward_array_store in array_graph domain not implemented");
  }
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  array_graph_domain_t invariant) {
    CRAB_WARN(
        "backward_array_store_range in array_graph domain not implemented");
  }
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v,
                                  array_graph_domain_t invariant) {
    CRAB_WARN(
        "backward_array_store_range in array_graph domain not implemented");
  }
  void backward_array_assign(variable_t lhs, variable_t rhs,
                             array_graph_domain_t invariant) {
    CRAB_WARN("backward_array_assign in array_graph domain not implemented");
  }

  void write(crab_os &o) {
#if 1
    NumDom copy_scalar(_scalar);
    array_graph_t copy_g(_g);
    // Remove all primed variables for pretty printing
    for (auto lm : var_lm_primes()) {
      copy_scalar -= get_var(lm);
      copy_g -= lm;
    }
    for (auto lm : cst_lm_primes()) {
      copy_scalar -= get_var(lm);
      copy_g -= lm;
    }
    o << "(" << copy_scalar << ",";
    copy_g.write(o, false); // we do not print bottom edges
    o << ")";
    // o << "##" << _expressions;
#else
    o << "("
      << "S=" << _scalar << ","
      << "E=" << _expressions << ","
      << "G=" << _g << ")";
#endif
  }

  // XXX: the array domain is disjunctive so it is not really
  // useful to express it through a conjunction of linear
  // constraints
  linear_constraint_system_t to_linear_constraint_system() {
    CRAB_ERROR("array-sgraph does not implement to_linear_constraint_system");
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    CRAB_ERROR("TODO: array-sgraph does not implement "
               "to_disjunctive_linear_constraint_system");
  }

  static std::string getDomainName() {
    std::string name("ArrayGraph(" + NumDom::getDomainName() + "," +
                     Content::getDomainName() + ")");
    return name;
  }
};

template <typename Dom, typename Content, bool IsDistContent>
struct abstract_domain_traits<array_graph_domain<Dom, Content, IsDistContent>> {
  // assume Dom::variable_t = Content::variable_t
  typedef typename Dom::number_t number_t;
  typedef typename Dom::varname_t varname_t;
};

template <typename Dom, typename Content, bool IsDistContent>
class array_graph_domain_traits<
    array_graph_domain<Dom, Content, IsDistContent>> {
public:
  template <class CFG> static void do_initialization(CFG cfg) {
    array_graph_domain<Dom, Content, IsDistContent>::do_initialization(cfg);
  }
};

template <typename Dom, typename Content, bool IsDistContent>
class special_domain_traits<array_graph_domain<Dom, Content, IsDistContent>> {
public:
  static void clear_global_state(void) {
    array_graph_domain<Dom, Content, IsDistContent>::clear_global_state();
  }
};

// Static data allocation
template <class Dom, class Content, bool IsDistContent>
landmark_ref_unordered_map<
    typename Dom::variable_t, typename Dom::number_t,
    landmark_ref<typename Dom::variable_t, typename Dom::number_t>>
    array_graph_domain<Dom, Content, IsDistContent>::s_var_landmarks;

template <class Dom, class Content, bool IsDistContent>
landmark_ref_unordered_map<
    typename Dom::variable_t, typename Dom::number_t,
    landmark_ref<typename Dom::variable_t, typename Dom::number_t>>
    array_graph_domain<Dom, Content, IsDistContent>::s_cst_landmarks;

} // end namespace domains
} // end namespace crab
#pragma GCC diagnostic pop
