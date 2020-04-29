/*******************************************************************************
 *
 * Sparse DBM implementation, with the same underlying architecture
 * as SplitDBM
 *
 * Graeme Gange (gkgange@unimelb.edu.au)
 * Jorge A. Navas (jorge.navas@sri.com)
 ******************************************************************************/

#pragma once

#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/common/types.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/graphs/graph_config.hpp>
#include <crab/domains/graphs/graph_ops.hpp>
#include <crab/domains/interval.hpp>
#include <crab/domains/linear_constraints.hpp>

#include <boost/container/flat_map.hpp>
#include <boost/optional.hpp>
#include <unordered_set>

//#define CHECK_POTENTIAL
//#define SDBM_NO_NORMALIZE

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {

namespace domains {

template <class Number, class VariableName,
          class Params = DBM_impl::DefaultParams<Number>>
class SparseDBM_ final
    : public abstract_domain<SparseDBM_<Number, VariableName, Params>> {
  typedef SparseDBM_<Number, VariableName, Params> DBM_t;
  typedef abstract_domain<DBM_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  typedef typename linear_constraint_t::kind_t constraint_kind_t;
  typedef ikos::interval<number_t> interval_t;

private:
  typedef ikos::bound<number_t> bound_t;
  typedef typename Params::Wt Wt;
  typedef typename Params::graph_t graph_t;
  typedef DBM_impl::NtoW<number_t, Wt> ntow;
  typedef typename graph_t::vert_id vert_id;
  typedef boost::container::flat_map<variable_t, vert_id> vert_map_t;
  typedef typename vert_map_t::value_type vmap_elt_t;
  typedef std::vector<boost::optional<variable_t>> rev_map_t;
  typedef GraphOps<graph_t> GrOps;
  typedef GraphPerm<graph_t> GrPerm;
  typedef typename GrOps::edge_vector edge_vector;
  // < <x, y>, k> == x - y <= k.
  typedef std::pair<std::pair<variable_t, variable_t>, Wt> diffcst_t;
  typedef std::unordered_set<vert_id> vert_set_t;

protected:
  //================
  // Domain data
  //================
  vert_map_t vert_map; // Mapping from variables to vertices
  rev_map_t rev_map;
  graph_t g;                 // The underlying relation graph
  std::vector<Wt> potential; // Stored potential for the vertex
  vert_set_t unstable;
  bool _is_bottom;

  /*
  void forget(std::vector<int> idxs) {
    dbm ret = NULL;
    ret = dbm_forget_array(&idxs[0], idxs.size(), _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
  }
  */

  class Wt_max {
  public:
    Wt_max() {}
    Wt apply(const Wt &x, const Wt &y) { return max(x, y); }
    bool default_is_absorbing() { return true; }
  };

  class Wt_min {
  public:
    Wt_min() {}
    Wt apply(const Wt &x, const Wt &y) { return std::min(x, y); }
    bool default_is_absorbing() { return false; }
  };

  vert_id get_vert(variable_t v) {
    auto it = vert_map.find(v);
    if (it != vert_map.end())
      return (*it).second;

    vert_id vert(g.new_vertex());
    // Initialize
    assert(vert <= rev_map.size());
    if (vert < rev_map.size()) {
      assert(!rev_map[vert]);
      potential[vert] = Wt(0);
      rev_map[vert] = v;
    } else {
      potential.push_back(Wt(0));
      rev_map.push_back(v);
    }
    vert_map.insert(vmap_elt_t(v, vert));

    assert(vert != 0);

    return vert;
  }

  vert_id get_vert(graph_t &g, vert_map_t &vmap, rev_map_t &rmap,
                   std::vector<Wt> &pot, variable_t v) {
    auto it = vmap.find(v);
    if (it != vmap.end())
      return (*it).second;

    vert_id vert(g.new_vertex());
    // vmap.insert(vmap_elt_t(v, vert));
    // Initialize
    assert(vert <= rmap.size());
    if (vert < rmap.size()) {
      assert(!rmap[vert]);
      pot[vert] = Wt(0);
      rmap[vert] = v;
    } else {
      pot.push_back(Wt(0));
      rmap.push_back(v);
    }
    vmap.insert(vmap_elt_t(v, vert));

    return vert;
  }

  template <class G, class P> inline bool check_potential(G &g, P &p) {
#ifdef CHECK_POTENTIAL
    for (vert_id v : g.verts()) {
      for (vert_id d : g.succs(v)) {
        if (p[v] + g.edge_val(v, d) - p[d] < Wt(0)) {
          assert(0 && "Invalid potential.");
          return false;
        }
      }
    }
#endif
    return true;
  }

  class vert_set_wrap_t {
  public:
    vert_set_wrap_t(const vert_set_t &_vs) : vs(_vs) {}

    bool operator[](vert_id v) const { return vs.find(v) != vs.end(); }
    const vert_set_t &vs;
  };

  // Evaluate the potential value of a variable.
  Wt pot_value(variable_t v) {
    auto it = vert_map.find(v);
    if (it != vert_map.end())
      return potential[(*it).second];

    return ((Wt)0);
  }

  Wt pot_value(variable_t v, std::vector<Wt> &potential) {
    auto it = vert_map.find(v);
    if (it != vert_map.end())
      return potential[(*it).second];

    return ((Wt)0);
  }

  // Evaluate an expression under the chosen potentials
  Wt eval_expression(linear_expression_t e, bool overflow) {
    Wt v(ntow::convert(e.constant(), overflow));
    if (overflow) {
      return Wt(0);
    }

    for (auto p : e) {
      Wt coef = ntow::convert(p.first, overflow);
      if (overflow) {
        return Wt(0);
      }
      v += (pot_value(p.second) - potential[0]) * coef;
    }
    return v;
  }

  interval_t eval_interval(linear_expression_t e) {
    interval_t r = e.constant();
    for (auto p : e) {
      r += p.first * operator[](p.second);
    }
    return r;
  }

  interval_t compute_residual(linear_expression_t e, variable_t pivot) {
    interval_t residual(-e.constant());
    for (typename linear_expression_t::iterator it = e.begin(); it != e.end();
         ++it) {
      variable_t v = it->second;
      if (v.index() != pivot.index()) {
        residual = residual - (interval_t(it->first) * this->operator[](v));
      }
    }
    return residual;
  }

  interval_t get_interval(variable_t x) { return get_interval(vert_map, g, x); }

  interval_t get_interval(vert_map_t &m, graph_t &r, variable_t x) {
    auto it = m.find(x);
    if (it == m.end()) {
      return interval_t::top();
    }
    vert_id v = (*it).second;
    interval_t x_out = interval_t(
        r.elem(v, 0) ? -number_t(r.edge_val(v, 0)) : bound_t::minus_infinity(),
        r.elem(0, v) ? number_t(r.edge_val(0, v)) : bound_t::plus_infinity());
    return x_out;
    /*
    boost::optional< interval_t > v = r.lookup(x);
    if(v)
      return *v;
    else
      return interval_t::top();
      */
  }

  // Turn an assignment into a set of difference constraints.
  void diffcsts_of_assign(variable_t x, linear_expression_t exp,
                          std::vector<std::pair<variable_t, Wt>> &lb,
                          std::vector<std::pair<variable_t, Wt>> &ub) {
    {
      // Process upper bounds.
      boost::optional<variable_t> unbounded_ubvar;
      bool overflow;

      Wt exp_ub(ntow::convert(exp.constant(), overflow));
      if (overflow) {
        return;
      }

      std::vector<std::pair<variable_t, Wt>> ub_terms;
      for (auto p : exp) {
        Wt coeff(ntow::convert(p.first, overflow));
        if (overflow) {
          continue;
        }
        if (coeff < Wt(0)) {
          // Can't do anything with negative coefficients.
          bound_t y_lb = operator[](p.second).lb();
          if (y_lb.is_infinite())
            goto assign_ub_finish;
          exp_ub += ntow::convert(*(y_lb.number()), overflow) * coeff;
          if (overflow) {
            continue;
          }
        } else {
          variable_t y(p.second);
          bound_t y_ub = operator[](y).ub();
          if (y_ub.is_infinite()) {
            if (unbounded_ubvar || coeff != Wt(1))
              goto assign_ub_finish;
            unbounded_ubvar = y;
          } else {
            Wt ymax(ntow::convert(*(y_ub.number()), overflow));
            if (overflow) {
              continue;
            }
            exp_ub += ymax * coeff;
            ub_terms.push_back({y, ymax});
          }
        }
      }

      if (unbounded_ubvar) {
        // There is exactly one unbounded variable.
        ub.push_back({*unbounded_ubvar, exp_ub});
      } else {
        for (auto p : ub_terms) {
          ub.push_back({p.first, exp_ub - p.second});
        }
      }
    }
  assign_ub_finish :

  {
    boost::optional<variable_t> unbounded_lbvar;
    bool overflow;

    Wt exp_lb(ntow::convert(exp.constant(), overflow));
    if (overflow) {
      return;
    }
    std::vector<std::pair<variable_t, Wt>> lb_terms;
    for (auto p : exp) {
      Wt coeff(ntow::convert(p.first, overflow));
      if (overflow) {
        continue;
      }
      if (coeff < Wt(0)) {
        // Again, can't do anything with negative coefficients.
        bound_t y_ub = operator[](p.second).ub();
        if (y_ub.is_infinite())
          goto assign_lb_finish;
        exp_lb += (ntow::convert(*(y_ub.number()), overflow)) * coeff;
        if (overflow) {
          continue;
        }
      } else {
        variable_t y(p.second);
        bound_t y_lb = operator[](y).lb();
        if (y_lb.is_infinite()) {
          if (unbounded_lbvar || coeff != Wt(1))
            goto assign_lb_finish;
          unbounded_lbvar = y;
        } else {
          Wt ymin(ntow::convert(*(y_lb.number()), overflow));
          if (overflow) {
            continue;
          }
          exp_lb += ymin * coeff;
          lb_terms.push_back({y, ymin});
        }
      }
    }

    if (unbounded_lbvar) {
      lb.push_back({*unbounded_lbvar, exp_lb});
    } else {
      for (auto p : lb_terms) {
        lb.push_back({p.first, exp_lb - p.second});
      }
    }
  }
  assign_lb_finish:
    return;
  }

  // GKG: I suspect there're some sign/bound direction errors in the
  // following.
  void diffcsts_of_lin_leq(const linear_expression_t &exp,
                           std::vector<diffcst_t> &csts,
                           std::vector<std::pair<variable_t, Wt>> &lbs,
                           std::vector<std::pair<variable_t, Wt>> &ubs) {
    // Process upper bounds.
    Wt unbounded_lbcoeff;
    Wt unbounded_ubcoeff;
    boost::optional<variable_t> unbounded_lbvar;
    boost::optional<variable_t> unbounded_ubvar;
    bool underflow, overflow;

    Wt exp_ub = -(ntow::convert(exp.constant(), overflow));
    if (overflow) {
      return;
    }

    // temporary hack
    ntow::convert(exp.constant() - 1, underflow);
    if (underflow) {
      // We don't like MIN either because the code will compute
      // minus MIN and it will silently overflow.
      return;
    }

    std::vector<std::pair<std::pair<Wt, variable_t>, Wt>> pos_terms, neg_terms;
    for (auto p : exp) {
      Wt coeff(ntow::convert(p.first, overflow));
      if (overflow) {
        continue;
      }
      if (coeff > Wt(0)) {
        variable_t y(p.second);
        bound_t y_lb = operator[](y).lb();
        if (y_lb.is_infinite()) {
          if (unbounded_lbvar)
            goto diffcst_finish;
          unbounded_lbvar = y;
          unbounded_lbcoeff = coeff;
        } else {
          Wt ymin(ntow::convert(*(y_lb.number()), overflow));
          if (overflow) {
            continue;
          }
          // Coeff is negative, so it's still add
          exp_ub -= ymin * coeff;
          pos_terms.push_back({{coeff, y}, ymin});
        }
      } else {
        variable_t y(p.second);
        bound_t y_ub = operator[](y).ub();
        if (y_ub.is_infinite()) {
          if (unbounded_ubvar)
            goto diffcst_finish;
          unbounded_ubvar = y;
          unbounded_ubcoeff = -coeff;
        } else {
          Wt ymax(ntow::convert(*(y_ub.number()), overflow));
          if (overflow) {
            continue;
          }
          exp_ub -= ymax * coeff;
          neg_terms.push_back({{-coeff, y}, ymax});
        }
      }
    }

    if (unbounded_lbvar) {
      variable_t x(*unbounded_lbvar);
      if (unbounded_ubvar) {
        if (unbounded_lbcoeff != Wt(1) || unbounded_ubcoeff != Wt(1))
          goto diffcst_finish;
        variable_t y(*unbounded_ubvar);
        csts.push_back({{x, y}, exp_ub});
      } else {
        if (unbounded_lbcoeff == Wt(1)) {
          for (auto p : neg_terms)
            csts.push_back({{x, p.first.second}, exp_ub - p.second});
        }
        // Add bounds for x
        ubs.push_back({x, exp_ub / unbounded_lbcoeff});
      }
    } else {
      if (unbounded_ubvar) {
        variable_t y(*unbounded_ubvar);
        if (unbounded_ubcoeff == Wt(1)) {
          for (auto p : pos_terms)
            csts.push_back({{p.first.second, y}, exp_ub + p.second});
        }
        // Bounds for y
        lbs.push_back({y, -exp_ub / unbounded_ubcoeff});
      } else {
        for (auto pl : neg_terms) {
          for (auto pu : pos_terms) {
            csts.push_back({{pu.first.second, pl.first.second},
                            exp_ub - pl.second + pu.second});
          }
        }

        for (auto pl : neg_terms) {
          lbs.push_back(
              {pl.first.second, -exp_ub / pl.first.first + pl.second});
        }
        for (auto pu : pos_terms) {
          ubs.push_back({pu.first.second, exp_ub / pu.first.first + pu.second});
        }
      }
    }
  diffcst_finish:
    return;
  }

  bool add_linear_leq(const linear_expression_t &exp) {
    CRAB_LOG("zones-sparse", linear_expression_t exp_tmp(exp);
             crab::outs() << "Adding: " << exp_tmp << "<= 0"
                          << "\n");
    std::vector<std::pair<variable_t, Wt>> lbs, ubs;
    std::vector<diffcst_t> csts;
    diffcsts_of_lin_leq(exp, csts, lbs, ubs);

    assert(check_potential(g, potential));

    Wt_min min_op;

    edge_vector es;
    for (auto p : lbs) {
      es.push_back({{get_vert(p.first), 0}, -p.second});
    }
    for (auto p : ubs) {
      es.push_back({{0, get_vert(p.first)}, p.second});
    }
    for (auto diff : csts) {
      CRAB_LOG("zones-sparse", crab::outs() << diff.first.first << "-"
                                            << diff.first.second
                                            << "<=" << diff.second << "\n";);
      es.push_back({{get_vert(diff.first.second), get_vert(diff.first.first)},
                    diff.second});
    }

    for (auto edge : es) {
      // CRAB_LOG("zones-sparse",
      // crab::outs() << diff.first.first<< "-"<< diff.first.second<< "<="
      //              << diff.second<<"\n";);

      vert_id src = edge.first.first;
      vert_id dest = edge.first.second;
      g.update_edge(src, edge.second, dest, min_op);
      if (!repair_potential(src, dest)) {
        set_to_bottom();
        return false;
      }
      assert(check_potential(g, potential));

      close_over_edge(src, dest);
      assert(check_potential(g, potential));
    }

    assert(check_potential(g, potential));
    return true;
  }

  // x != n
  void add_univar_disequation(variable_t x, number_t n) {
    bool overflow;
    interval_t i = get_interval(x);
    interval_t new_i = linear_interval_solver_impl::trim_interval<interval_t>(
        i, interval_t(n));
    if (new_i.is_bottom()) {
      set_to_bottom();
    } else if (!new_i.is_top() && (new_i <= i)) {
      vert_id v = get_vert(x);
      Wt_min min_op;
      typename graph_t::mut_val_ref_t w;
      if (new_i.lb().is_finite()) {
        // strenghten lb
        Wt lb_val = ntow::convert(-(*(new_i.lb().number())), overflow);
        if (overflow) {
          return;
        }
        if (g.lookup(v, 0, &w) && lb_val < w) {
          g.set_edge(v, lb_val, 0);
          if (!repair_potential(v, 0)) {
            set_to_bottom();
            return;
          }
          assert(check_potential(g, potential));
          // Update other bounds
          for (auto e : g.e_preds(v)) {
            if (e.vert == 0)
              continue;
            g.update_edge(e.vert, e.val + lb_val, 0, min_op);
            if (!repair_potential(e.vert, 0)) {
              set_to_bottom();
              return;
            }
            assert(check_potential(g, potential));
          }
        }
      }
      if (new_i.ub().is_finite()) {
        // strengthen ub
        Wt ub_val = ntow::convert(*(new_i.ub().number()), overflow);
        if (overflow) {
          return;
        }
        if (g.lookup(0, v, &w) && (ub_val < w)) {
          g.set_edge(0, ub_val, v);
          if (!repair_potential(0, v)) {
            set_to_bottom();
            return;
          }
          assert(check_potential(g, potential));
          // Update other bounds
          for (auto e : g.e_succs(v)) {
            if (e.vert == 0)
              continue;
            g.update_edge(0, e.val + ub_val, e.vert, min_op);
            if (!repair_potential(0, e.vert)) {
              set_to_bottom();
              return;
            }
            assert(check_potential(g, potential));
          }
        }
      }
    }
  }

  void add_disequation(linear_expression_t e) {
    // XXX: similar precision as the interval domain
    for (typename linear_expression_t::iterator it = e.begin(); it != e.end();
         ++it) {
      variable_t pivot = it->second;
      interval_t i = compute_residual(e, pivot) / interval_t(it->first);
      if (auto k = i.singleton()) {
        add_univar_disequation(pivot, *k);
      }
    }
  }

  // Restore potential after an edge addition
  bool repair_potential(vert_id src, vert_id dest) {
    return GrOps::repair_potential(g, potential, src, dest);
  }

  // Restore closure after a single edge addition
  void close_over_edge(vert_id ii, vert_id jj) {
    Wt_min min_op;

    Wt c = g.edge_val(ii, jj);

    typename graph_t::mut_val_ref_t w;

    // There may be a cheaper way to do this.
    // GKG: Now implemented.
    std::vector<std::pair<vert_id, Wt>> src_dec;

    for (auto edge : g.e_preds(ii)) {
      vert_id se = edge.vert;
      Wt w_si = edge.val;
      Wt wt_sij = w_si + c;

      assert(g.succs(se).begin() != g.succs(se).end());
      if (se != jj) {
        if (g.lookup(se, jj, &w)) {
          if (w.get() <= wt_sij)
            continue;
          w = wt_sij;
        } else {
          g.add_edge(se, wt_sij, jj);
        }
        // assert(potential[se] + g.edge_val(se, jj) - potential[jj] >= Wt(0));
        src_dec.push_back({se, w_si});

        /*
         for(auto edge : g.e_succs(jj))
         {
           vert_id de = edge.vert;
           if(se != de)
           {
             Wt wt_sijd = wt_sij + edge.val;
             if(g.lookup(se, de, &w))
             {
               if((*w) <= wt_sijd)
                 continue;
               (*w) = wt_sijd;
             } else {
               g.add_edge(se, wt_sijd, de);
             }
           }
         }
         */
      }
    }

    std::vector<std::pair<vert_id, Wt>> dest_dec;
    for (auto edge : g.e_succs(jj)) {
      vert_id de = edge.vert;
      Wt w_jd = edge.val;
      Wt wt_ijd = w_jd + c;
      if (de != ii) {
        if (g.lookup(ii, de, &w)) {
          if (w.get() <= wt_ijd)
            continue;
          w = wt_ijd;
        } else {
          g.add_edge(ii, wt_ijd, de);
        }
        // assert(potential[ii] + g.edge_val(ii, de) - potential[de] >= Wt(0));
        // dest_dec.push_back(std::make_pair(de, edge.val));
        dest_dec.push_back({de, w_jd});
      }
    }
    // Look at (src, dest) pairs with updated edges.
    for (auto s_p : src_dec) {
      vert_id se = s_p.first;
      Wt wt_sij = c + s_p.second;
      for (auto d_p : dest_dec) {
        vert_id de = d_p.first;
        Wt wt_sijd = wt_sij + d_p.second;
        if (g.lookup(se, de, &w)) {
          if (w.get() <= wt_sijd)
            continue;
          w = wt_sijd;
        } else {
          g.add_edge(se, wt_sijd, de);
        }
        //  assert(potential[se] + g.edge_val(se, de) - potential[de] >= Wt(0));
      }
    }
    // Closure is now updated.
  }

  // Restore closure after a variable assignment
  // Assumption: x = f(y_1, ..., y_n) cannot induce non-trivial
  // relations between (y_i, y_j)
  /*
  bool close_after_assign(vert_id v)
  {
    // Run Dijkstra's forward to collect successors of v,
    // and backward to collect predecessors
    edge_vector delta;
    if(!GrOps::close_after_assign(g, potential, v, delta))
      return false;
    GrOps::apply_delta(g, delta);
    return true;
  }

  bool closure(void)
  {
    // Full Johnson-style all-pairs shortest path
    CRAB_ERROR("SparseWtGraph::closure not yet implemented.");
  }
  */
  
public:
  SparseDBM_(bool is_bottom = false) : _is_bottom(is_bottom) {
    g.growTo(1); // Allocate the zero vector
    potential.push_back(Wt(0));
    rev_map.push_back(boost::none);
  }

  // FIXME: Rewrite to avoid copying if o is _|_
  SparseDBM_(const DBM_t &o)
      : vert_map(o.vert_map), rev_map(o.rev_map), g(o.g),
        potential(o.potential), unstable(o.unstable), _is_bottom(false) {

    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");

    if (o._is_bottom)
      set_to_bottom();

    if (!_is_bottom)
      assert(g.size() > 0);
  }

  SparseDBM_(DBM_t &&o)
      : vert_map(std::move(o.vert_map)), rev_map(std::move(o.rev_map)),
        g(std::move(o.g)), potential(std::move(o.potential)),
        unstable(std::move(o.unstable)), _is_bottom(o._is_bottom) {}

  // Magical rvalue ownership stuff for efficient initialization
  SparseDBM_(vert_map_t &&_vert_map, rev_map_t &&_rev_map, graph_t &&_g,
             std::vector<Wt> &&_potential, vert_set_t &&_unstable)
      : vert_map(std::move(_vert_map)), rev_map(std::move(_rev_map)),
        g(std::move(_g)), potential(std::move(_potential)),
        unstable(std::move(_unstable)), _is_bottom(false) {

    CRAB_LOG("zones-sparse-size", auto p = size();
             crab::outs() << "#nodes = " << p.first << " #edges=" << p.second
                          << "\n";);

    assert(g.size() > 0);
  }

  SparseDBM_ &operator=(const SparseDBM_ &o) {
    crab::CrabStats::count(getDomainName() + ".count.copy");
    crab::ScopedCrabStats __st__(getDomainName() + ".copy");

    if (this != &o) {
      if (o._is_bottom)
        set_to_bottom();
      else {
        _is_bottom = false;
        vert_map = o.vert_map;
        rev_map = o.rev_map;
        g = o.g;
        potential = o.potential;
        unstable = o.unstable;
        assert(g.size() > 0);
      }
    }
    return *this;
  }

  SparseDBM_ &operator=(SparseDBM_ &&o) {
    if (o._is_bottom) {
      set_to_bottom();
    } else {
      _is_bottom = false;
      vert_map = std::move(o.vert_map);
      rev_map = std::move(o.rev_map);
      g = std::move(o.g);
      potential = std::move(o.potential);
      unstable = std::move(o.unstable);
    }
    return *this;
  }

  void set_to_top() {
    SparseDBM_ abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    vert_map.clear();
    rev_map.clear();
    g.clear();
    potential.clear();
    unstable.clear();
    _is_bottom = true;
  }

  // void set_to_bottom() {
  // 	SparseDBM_ abs(true);
  // 	std::swap(*this, abs);
  // }

  bool is_bottom() {
    // if(!_is_bottom && g.has_negative_cycle())
    // _is_bottom = true;
    return _is_bottom;
  }

  bool is_top() {
    if (_is_bottom)
      return false;
    return g.is_empty();
  }

  bool operator<=(DBM_t o) {
    crab::CrabStats::count(getDomainName() + ".count.leq");
    crab::ScopedCrabStats __st__(getDomainName() + ".leq");

    // cover all trivial cases to avoid allocating a dbm matrix
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

      // CRAB_LOG("zones-sparse",
      //          crab::outs() << "operator<=: "<< *this<< "<=?"<< o << "\n");

      if (vert_map.size() < o.vert_map.size())
        return false;

      typename graph_t::mut_val_ref_t wx;

      // Set up a mapping from o to this.
      std::vector<unsigned int> vert_renaming(o.g.size(), -1);
      vert_renaming[0] = 0;
      for (auto p : o.vert_map) {
        auto it = vert_map.find(p.first);
        // We can't have this <= o if we're missing some
        // vertex.
        if (it == vert_map.end())
          return false;
        vert_renaming[p.second] = (*it).second;
      }

      assert(g.size() > 0);
      // GrPerm g_perm(vert_renaming, g);

      for (vert_id ox : o.g.verts()) {
        assert(vert_renaming[ox] != -1);
        vert_id x = vert_renaming[ox];
        for (auto edge : o.g.e_succs(ox)) {
          vert_id oy = edge.vert;
          assert(vert_renaming[ox] != -1);
          vert_id y = vert_renaming[oy];
          Wt ow = edge.val;

          if (!g.lookup(x, y, &wx) || (ow < wx))
            return false;
        }
      }
      return true;
    }
  }

  // FIXME: can be done more efficient
  void operator|=(DBM_t o) { *this = *this | o; }

  DBM_t operator|(DBM_t o) {
    crab::CrabStats::count(getDomainName() + ".count.join");
    crab::ScopedCrabStats __st__(getDomainName() + ".join");

    if (is_bottom() || o.is_top())
      return o;
    else if (is_top() || o.is_bottom())
      return *this;
    else {
      CRAB_LOG("zones-sparse", crab::outs() << "Before join:\n"
                                            << "DBM 1\n"
                                            << *this << "\n"
                                            << "DBM 2\n"
                                            << o << "\n");

      normalize();
      o.normalize();

      assert(check_potential(g, potential));
      assert(check_potential(o.g, o.potential));

      // Figure out the common renaming, initializing the
      // resulting potentials as we go.
      std::vector<vert_id> perm_x;
      std::vector<vert_id> perm_y;
      std::vector<variable_t> perm_inv;

      std::vector<Wt> pot_rx;
      std::vector<Wt> pot_ry;
      vert_map_t out_vmap;
      rev_map_t out_revmap;
      // Add the zero vertex
      assert(potential.size() > 0);
      pot_rx.push_back(0);
      pot_ry.push_back(0);
      perm_x.push_back(0);
      perm_y.push_back(0);
      out_revmap.push_back(boost::none);

      for (auto p : vert_map) {
        auto it = o.vert_map.find(p.first);
        // Variable exists in both
        if (it != o.vert_map.end()) {
          out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
          out_revmap.push_back(p.first);

          pot_rx.push_back(potential[p.second] - potential[0]);
          pot_ry.push_back(o.potential[(*it).second] - o.potential[0]);
          perm_inv.push_back(p.first);
          perm_x.push_back(p.second);
          perm_y.push_back((*it).second);
        }
      }
      // unsigned int sz = perm_x.size();

      // Build the permuted view of x and y.
      assert(g.size() > 0);
      GrPerm gx(perm_x, g);
      assert(o.g.size() > 0);
      GrPerm gy(perm_y, o.g);

      // We now have the relevant set of relations. Because g_rx and g_ry are
      // closed, the result is also closed.
      Wt_min min_op;
      graph_t join_g(GrOps::join(gx, gy));

      // Now garbage collect any unused vertices
      for (vert_id v : join_g.verts()) {
        if (v == 0)
          continue;
        if (join_g.succs(v).size() == 0 && join_g.preds(v).size() == 0) {
          join_g.forget(v);
          if (out_revmap[v]) {
            out_vmap.erase(*(out_revmap[v]));
            out_revmap[v] = boost::none;
          }
        }
      }

      DBM_t res(std::move(out_vmap), std::move(out_revmap), std::move(join_g),
                std::move(pot_rx), vert_set_t());
      CRAB_LOG("zones-sparse", crab::outs() << "Result join:\n"
                                            << res << "\n";);

      return res;
    }
  }

  DBM_t operator||(DBM_t o) {
    crab::CrabStats::count(getDomainName() + ".count.widening");
    crab::ScopedCrabStats __st__(getDomainName() + ".widening");

    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else {
      CRAB_LOG("zones-sparse",
               DBM_t left(*this); // to avoid closure on left operand
               crab::outs() << "Before widening:\n"
                            << "DBM 1\n"
                            << left << "\n"
                            << "DBM 2\n"
                            << o << "\n";);
      o.normalize();

      // Figure out the common renaming
      std::vector<vert_id> perm_x;
      std::vector<vert_id> perm_y;
      vert_map_t out_vmap;
      rev_map_t out_revmap;
      std::vector<Wt> widen_pot;
      vert_set_t widen_unstable(unstable);

      assert(potential.size() > 0);
      widen_pot.push_back(Wt(0));
      perm_x.push_back(0);
      perm_y.push_back(0);
      out_revmap.push_back(boost::none);
      for (auto p : vert_map) {
        auto it = o.vert_map.find(p.first);
        // Variable exists in both
        if (it != o.vert_map.end()) {
          out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
          out_revmap.push_back(p.first);

          widen_pot.push_back(potential[p.second] - potential[0]);
          perm_x.push_back(p.second);
          perm_y.push_back((*it).second);
        }
      }

      // Build the permuted view of x and y.
      assert(g.size() > 0);
      GrPerm gx(perm_x, g);
      assert(o.g.size() > 0);
      GrPerm gy(perm_y, o.g);

      // Now perform the widening
      std::vector<vert_id> destabilized;
      graph_t widen_g(GrOps::widen(gx, gy, destabilized));
      for (vert_id v : destabilized)
        widen_unstable.insert(v);

      DBM_t res(std::move(out_vmap), std::move(out_revmap), std::move(widen_g),
                std::move(widen_pot), std::move(widen_unstable));

      CRAB_LOG("zones-sparse", DBM_t res_copy(res); crab::outs()
                                                    << "Result widening:\n"
                                                    << res_copy << "\n";);
      return res;
    }
  }

  DBM_t operator&(DBM_t o) {
    crab::CrabStats::count(getDomainName() + ".count.meet");
    crab::ScopedCrabStats __st__(getDomainName() + ".meet");

    if (is_bottom() || o.is_bottom())
      return DBM_t::bottom();
    else if (is_top())
      return o;
    else if (o.is_top())
      return *this;
    else {
      CRAB_LOG("zones-sparse", crab::outs() << "Before meet:\n"
                                            << "DBM 1\n"
                                            << *this << "\n"
                                            << "DBM 2\n"
                                            << o << "\n";);
      normalize();
      o.normalize();

      // We map vertices in the left operand onto a contiguous range.
      // This will often be the identity map, but there might be gaps.
      vert_map_t meet_verts;
      rev_map_t meet_rev;

      std::vector<vert_id> perm_x;
      std::vector<vert_id> perm_y;
      std::vector<Wt> meet_pi;
      perm_x.push_back(0);
      perm_y.push_back(0);
      meet_pi.push_back(Wt(0));
      meet_rev.push_back(boost::none);
      for (auto p : vert_map) {
        vert_id vv = perm_x.size();
        meet_verts.insert(vmap_elt_t(p.first, vv));
        meet_rev.push_back(p.first);

        perm_x.push_back(p.second);
        perm_y.push_back(-1);
        meet_pi.push_back(potential[p.second] - potential[0]);
      }

      // Add missing mappings from the right operand.
      for (auto p : o.vert_map) {
        auto it = meet_verts.find(p.first);

        if (it == meet_verts.end()) {
          vert_id vv = perm_y.size();
          meet_rev.push_back(p.first);

          perm_y.push_back(p.second);
          perm_x.push_back(-1);
          meet_pi.push_back(o.potential[p.second] - o.potential[0]);
          meet_verts.insert(vmap_elt_t(p.first, vv));
        } else {
          perm_y[(*it).second] = p.second;
        }
      }

      // Build the permuted view of x and y.
      assert(g.size() > 0);
      GrPerm gx(perm_x, g);
      assert(o.g.size() > 0);
      GrPerm gy(perm_y, o.g);

      // Compute the syntactic meet of the permuted graphs.
      bool is_closed;
      graph_t meet_g(GrOps::meet(gx, gy, is_closed));

      // Compute updated potentials on the zero-enriched graph
      // std::vector<Wt> meet_pi(meet_g.size());
      // We've warm-started pi with the operand potentials
      if (!GrOps::select_potentials(meet_g, meet_pi)) {
        // Potentials cannot be selected -- state is infeasible.
        return DBM_t::bottom();
      }

      if (!is_closed) {
        edge_vector delta;
        if (Params::chrome_dijkstra)
          GrOps::close_after_meet(meet_g, meet_pi, gx, gy, delta);
        else
          GrOps::close_johnson(meet_g, meet_pi, delta);

        GrOps::apply_delta(meet_g, delta);
      }
      assert(check_potential(meet_g, meet_pi));
      DBM_t res(std::move(meet_verts), std::move(meet_rev), std::move(meet_g),
                std::move(meet_pi), vert_set_t());
      CRAB_LOG("zones-sparse", crab::outs() << "Result meet:\n"
                                            << res << "\n";);
      return res;
    }
  }

  DBM_t operator&&(DBM_t o) {
    crab::CrabStats::count(getDomainName() + ".count.narrowing");
    crab::ScopedCrabStats __st__(getDomainName() + ".narrowing");

    if (is_bottom() || o.is_bottom())
      return DBM_t::bottom();
    else if (is_top())
      return o;
    else {
      CRAB_LOG("zones-sparse", crab::outs() << "Before narrowing:\n"
                                            << "DBM 1\n"
                                            << *this << "\n"
                                            << "DBM 2\n"
                                            << o << "\n";);

      // FIXME: Implement properly
      // Narrowing as a no-op should be sound.
      normalize();
      DBM_t res(*this);

      CRAB_LOG("zones-sparse", crab::outs() << "Result narrowing:\n"
                                            << res << "\n";);
      return res;
    }
  }

  DBM_t widening_thresholds(DBM_t o,
                            const iterators::thresholds<number_t> &ts) {
    // TODO: use thresholds
    return (*this || o);
  }

  void normalize() {
// dbm_canonical(_dbm);
// Always maintained in normal form, except for widening
#ifdef SDBM_NO_NORMALIZE
    return;
#endif

    if (unstable.size() == 0)
      return;

    edge_vector delta;

    if (Params::widen_restabilize)
      GrOps::close_after_widen(g, potential, vert_set_wrap_t(unstable), delta);
    else
      GrOps::close_johnson(g, potential, delta);

    GrOps::apply_delta(g, delta);

    unstable.clear();
  }

  void minimize() {}

  void operator-=(variable_t v) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");

    if (is_bottom())
      return;
    normalize();

    auto it = vert_map.find(v);
    if (it != vert_map.end()) {
      CRAB_LOG("zones-sparse", crab::outs() << "Before forget " << it->second
                                            << ": " << g << "\n";);
      g.forget(it->second);
      CRAB_LOG("zones-sparse", crab::outs() << "After: " << g << "\n";);

      rev_map[it->second] = boost::none;
      vert_map.erase(v);
    }
  }

  // Assumption: state is currently feasible.
  void assign(variable_t x, linear_expression_t e) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    if (is_bottom()) {
      return;
    }

    normalize();

    assert(check_potential(g, potential));

    // If it's a constant, just assign the interval.
    if (e.is_constant()) {
      set(x, e.constant());
    } else {
      interval_t x_int = eval_interval(e);

      boost::optional<Wt> lb_w, ub_w;
      bool overflow;
      if (x_int.lb().is_finite()) {
        lb_w = ntow::convert(-(*(x_int.lb().number())), overflow);
        if (overflow) {
          operator-=(x);
          CRAB_LOG("zones-sparse", crab::outs()
                                       << "---" << x << ":=" << e << "\n"
                                       << *this << "\n");
          return;
        }
      }
      if (x_int.ub().is_finite()) {
        ub_w = ntow::convert(*(x_int.ub().number()), overflow);
        if (overflow) {
          operator-=(x);
          CRAB_LOG("zones-sparse", crab::outs()
                                       << "---" << x << ":=" << e << "\n"
                                       << *this << "\n");
          return;
        }
      }

      std::vector<std::pair<variable_t, Wt>> diffs_lb, diffs_ub;
      // Construct difference constraints from the assignment
      diffcsts_of_assign(x, e, diffs_lb, diffs_ub);
      if (diffs_lb.size() > 0 || diffs_ub.size() > 0) {
        if (Params::special_assign) {
          bool overflow;
          Wt e_val = eval_expression(e, overflow);
          if (overflow) {
            operator-=(x);
            return;
          }

          // Allocate a new vertex for x
          vert_id v = g.new_vertex();
          assert(v <= rev_map.size());
          if (v == rev_map.size()) {
            rev_map.push_back(x);
            potential.push_back(potential[0] + e_val);
          } else {
            potential[v] = potential[0] + e_val;
            rev_map[v] = x;
          }

          edge_vector delta;
          for (auto diff : diffs_lb) {
            delta.push_back({{v, get_vert(diff.first)}, -diff.second});
          }

          for (auto diff : diffs_ub) {
            delta.push_back({{get_vert(diff.first), v}, diff.second});
          }

          if (lb_w) {
            delta.push_back({{v, 0}, *lb_w});
          }

          if (ub_w) {
            delta.push_back({{0, v}, *ub_w});
          }

          GrOps::apply_delta(g, delta);
          delta.clear();
          GrOps::close_after_assign(g, potential, v, delta);
          GrOps::apply_delta(g, delta);

          // Clear the old x vertex
          operator-=(x);
          vert_map.insert(vmap_elt_t(x, v));
        } else {
          vert_id v = g.new_vertex();
          assert(v <= rev_map.size());
          if (v == rev_map.size()) {
            rev_map.push_back(x);
            potential.push_back(Wt(0));
          } else {
            assert(!rev_map[v]);
            potential[v] = Wt(0);
            rev_map[v] = x;
          }
          Wt_min min_op;
          edge_vector cst_edges;

          if (lb_w) {
            cst_edges.push_back({{v, 0}, *lb_w});
          }
          if (ub_w) {
            cst_edges.push_back({{0, v}, *ub_w});
          }

          for (auto diff : diffs_lb) {
            cst_edges.push_back({{v, get_vert(diff.first)}, -diff.second});
          }

          for (auto diff : diffs_ub) {
            cst_edges.push_back({{get_vert(diff.first), v}, diff.second});
          }

          for (auto diff : cst_edges) {
            vert_id src = diff.first.first;
            vert_id dest = diff.first.second;
            g.update_edge(src, diff.second, dest, min_op);
            if (!repair_potential(src, dest)) {
              assert(0 && "Unreachable");
              set_to_bottom();
            }
            assert(check_potential(g, potential));
            close_over_edge(src, dest);
            assert(check_potential(g, potential));
          }
          // Clear the old x vertex
          operator-=(x);
          vert_map.insert(vmap_elt_t(x, v));
        }
        assert(check_potential(g, potential));
      } else {
        set(x, x_int);
      }
      // CRAB_WARN("DBM only supports a cst or var on the rhs of assignment");
      // this->operator-=(x);
    }

    // g.check_adjs();

    assert(check_potential(g, potential));
    CRAB_LOG("zones-sparse", crab::outs() << "---" << x << ":=" << e << "\n"
                                          << *this << "\n";);
  }

  void apply(ikos::operation_t op, variable_t x, variable_t y, variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (is_bottom()) {
      return;
    }

    normalize();

    switch (op) {
    case ikos::OP_ADDITION:
      assign(x, y + z);
      break;
    case ikos::OP_SUBTRACTION:
      assign(x, y - z);
      break;
    // For the rest of operations, we fall back on intervals.
    case ikos::OP_MULTIPLICATION:
      set(x, get_interval(y) * get_interval(z));
      break;
    case ikos::OP_SDIV:
      set(x, get_interval(y) / get_interval(z));
      break;
    case ikos::OP_UDIV:
      set(x, get_interval(y).UDiv(get_interval(z)));
      break;
    case ikos::OP_SREM:
      set(x, get_interval(y).SRem(get_interval(z)));
      break;
    case ikos::OP_UREM:
      set(x, get_interval(y).URem(get_interval(z)));
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }
    CRAB_LOG("zones-sparse", crab::outs()
                                 << "---" << x << ":=" << y << op << z << "\n"
                                 << *this << "\n";);
  }

  void apply(ikos::operation_t op, variable_t x, variable_t y, number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (is_bottom()) {
      return;
    }

    normalize();

    switch (op) {
    case ikos::OP_ADDITION:
      assign(x, y + k);
      break;
    case ikos::OP_SUBTRACTION:
      assign(x, y - k);
      break;
    // For the rest of operations, we fall back on intervals.
    case ikos::OP_MULTIPLICATION:
      set(x, get_interval(y) * interval_t(k));
      break;
    case ikos::OP_SDIV:
      set(x, get_interval(y) / interval_t(k));
      break;
    case ikos::OP_UDIV:
      set(x, get_interval(y).UDiv(interval_t(k)));
      break;
    case ikos::OP_SREM:
      set(x, get_interval(y).SRem(interval_t(k)));
      break;
    case ikos::OP_UREM:
      set(x, get_interval(y).URem(interval_t(k)));
      break;
    default:
      CRAB_ERROR("Operation ", op, " not supported");
    }

    CRAB_LOG("zones-sparse", crab::outs()
                                 << "---" << x << ":=" << y << op << k << "\n"
                                 << *this << "\n";);
  }

  void operator+=(linear_constraint_t cst) {
    crab::CrabStats::count(getDomainName() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(getDomainName() + ".add_constraints");

    // XXX: we do nothing with unsigned linear inequalities
    if (cst.is_inequality() && cst.is_unsigned()) {
      CRAB_WARN("unsigned inequality ", cst, " skipped by split_dbm domain");
      return;
    }

    if (is_bottom())
      return;

    normalize();

    if (cst.is_tautology())
      return;

    // g.check_adjs();

    if (cst.is_contradiction()) {
      set_to_bottom();
      return;
    }

    if (cst.is_inequality()) {
      if (!add_linear_leq(cst.expression())) {
        set_to_bottom();
      }
      // g.check_adjs();
      CRAB_LOG("zones-sparse", crab::outs() << "--- " << cst << "\n"
                                            << *this << "\n";);
      return;
    }

    if (cst.is_strict_inequality()) {
      // We try to convert a strict to non-strict.
      auto nc = linear_constraint_impl::strict_to_non_strict_inequality(cst);
      if (nc.is_inequality()) {
        // here we succeed
        if (!add_linear_leq(nc.expression())) {
          set_to_bottom();
        }
        CRAB_LOG("zones-split", crab::outs() << "--- " << cst << "\n"
                                             << *this << "\n");
        return;
      }
    }

    if (cst.is_equality()) {
      linear_expression_t exp = cst.expression();
      if (!add_linear_leq(exp) || !add_linear_leq(-exp)) {
        CRAB_LOG("zones-sparse", crab::outs() << " ~~> _|_"
                                              << "\n";);
        set_to_bottom();
      }
      // g.check_adjs();
      CRAB_LOG("zones-sparse", crab::outs() << "--- " << cst << "\n"
                                            << *this << "\n";);
      return;
    }

    if (cst.is_disequation()) {
      add_disequation(cst.expression());
      return;
    }

    CRAB_WARN("Unhandled constraint ", cst, " by split_dbm");
    CRAB_LOG("zones-sparse", crab::outs() << "---" << cst << "\n"
                                          << *this << "\n";);
    return;
  }

  void operator+=(linear_constraint_system_t csts) {
    if (is_bottom())
      return;

    for (auto cst : csts) {
      operator+=(cst);
    }
  }

  interval_t operator[](variable_t x) {
    crab::CrabStats::count(getDomainName() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(getDomainName() + ".to_intervals");

    // if (is_top()) return interval_t::top();
    if (is_bottom()) {
      return interval_t::bottom();
    } else {
      variable_t vx(x);
      return get_interval(vert_map, g, vx);
    }
  }

  void set(variable_t x, interval_t intv) {
    crab::CrabStats::count(getDomainName() + ".count.assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".assign");

    if (is_bottom()) {
      return;
    }

    if (intv.is_bottom()) {
      set_to_bottom();
      return;
    }

    this->operator-=(x);

    if (intv.is_top()) {
      return;
    }

    vert_id v = get_vert(x);
    bool overflow;
    if (intv.ub().is_finite()) {
      Wt ub = ntow::convert(*(intv.ub().number()), overflow);
      if (overflow) {
        return;
      }
      potential[v] = potential[0] + ub;
      g.set_edge(0, ub, v);
      close_over_edge(0, v);
    }
    if (intv.lb().is_finite()) {
      Wt lb = ntow::convert(*(intv.lb().number()), overflow);
      if (overflow) {
        return;
      }
      potential[v] = potential[0] + lb;
      g.set_edge(v, -lb, 0);
      close_over_edge(v, 0);
    }
  }

  // backward arithmetic operators
  void backward_assign(variable_t x, linear_expression_t e, DBM_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_assign");

    crab::domains::BackwardAssignOps<DBM_t>::assign(*this, x, e, inv);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, number_t z,
                      DBM_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    crab::domains::BackwardAssignOps<DBM_t>::apply(*this, op, x, y, z, inv);
  }

  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      DBM_t inv) {
    crab::CrabStats::count(getDomainName() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".backward_apply");

    crab::domains::BackwardAssignOps<DBM_t>::apply(*this, op, x, y, z, inv);
  }

  // cast operators

  void apply(int_conv_operation_t /*op*/, variable_t dst, variable_t src) {
    // since reasoning about infinite precision we simply assign and
    // ignore the widths.
    assign(dst, src);
  }

  // bitwise operators

  void apply(ikos::bitwise_operation_t op, variable_t x, variable_t y,
             variable_t z) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (is_bottom())
      return;
    normalize();

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi = operator[](z);
    interval_t xi = interval_t::bottom();
    switch (op) {
    case ikos::OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case ikos::OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case ikos::OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case ikos::OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case ikos::OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case ikos::OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("DBM: unreachable");
    }
    set(x, xi);
  }

  void apply(ikos::bitwise_operation_t op, variable_t x, variable_t y,
             number_t k) {
    crab::CrabStats::count(getDomainName() + ".count.apply");
    crab::ScopedCrabStats __st__(getDomainName() + ".apply");

    if (is_bottom())
      return;
    normalize();

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi(k);
    interval_t xi = interval_t::bottom();

    switch (op) {
    case ikos::OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case ikos::OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case ikos::OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case ikos::OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case ikos::OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    case ikos::OP_ASHR: {
      xi = yi.AShr(zi);
      break;
    }
    default:
      CRAB_ERROR("DBM: unreachable");
    }
    set(x, xi);
  }

  /*
     Begin unimplemented operations

     SparseDBM implements only standard abstract operations of a
     numerical domain.  The implementation of boolean, array, or
     pointer operations is empty because they should never be
     called.
  */

  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                DBM_t invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                DBM_t invariant) {}
  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z, DBM_t invariant) {
  }
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t v) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           DBM_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           DBM_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, DBM_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, DBM_t invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v, DBM_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v, DBM_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs, DBM_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_store(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    crab::CrabStats::count(getDomainName() + ".count.rename");
    crab::ScopedCrabStats __st__(getDomainName() + ".rename");

    if (is_top() || is_bottom())
      return;

    // renaming vert_map by creating a new vert_map since we are
    // modifying the keys.
    // rev_map is modified in-place since we only modify values.
    CRAB_LOG("zones-sparse", crab::outs() << "Renaming {";
             for (auto v
                  : from) crab::outs()
             << v << ";";
             crab::outs() << "} with "; for (auto v
                                             : to) crab::outs()
                                        << v << ";";
             crab::outs() << "}:\n"; crab::outs() << *this << "\n";);


    for (unsigned i=0, sz=from.size(); i<sz; ++i) {
      variable_t v = from[i];
      variable_t new_v = to[i];
      if (v == new_v) { // nothing to rename
        continue;
      }

      { auto it = vert_map.find(new_v);
	if (it != vert_map.end()) {
	  CRAB_ERROR(getDomainName() + "::rename assumes that ", new_v, " does not exist");	  
	}
      }

      auto it = vert_map.find(v);
      if (it != vert_map.end()) {
	vert_id dim = it->second;
	vert_map.erase(it);
	vert_map.insert(vmap_elt_t(new_v, dim));
        rev_map[dim] = new_v;	
      }
    }

    CRAB_LOG("zones-sparse", crab::outs() << "RESULT=" << *this << "\n");
  }

  void forget(const variable_vector_t &variables) {
    crab::CrabStats::count(getDomainName() + ".count.forget");
    crab::ScopedCrabStats __st__(getDomainName() + ".forget");

    if (is_bottom() || is_top())
      return;

    for (auto v : variables) {
      auto it = vert_map.find(v);
      if (it != vert_map.end()) {
        operator-=(v);
      }
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

    normalize();

    std::vector<bool> save(rev_map.size(), false);
    for (auto x : variables) {
      auto it = vert_map.find(x);
      if (it != vert_map.end()) {
        save[(*it).second] = true;
      }
    }

    for (vert_id v = 0; v < rev_map.size(); v++) {
      if (!save[v] && rev_map[v])
        operator-=((*rev_map[v]));
    }
  }

  void expand(variable_t x, variable_t y) {
    crab::CrabStats::count(getDomainName() + ".count.expand");
    crab::ScopedCrabStats __st__(getDomainName() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    CRAB_LOG("zones-sparse", crab::outs() << "Before expand " << x << " into "
                                          << y << ":\n"
                                          << *this << "\n");

    auto it = vert_map.find(y);
    if (it != vert_map.end()) {
      CRAB_ERROR("sparse_dbm expand operation failed because y already exists");
    }

    vert_id ii = get_vert(x);
    vert_id jj = get_vert(y);

    for (auto edge : g.e_preds(ii)) {
      g.add_edge(edge.vert, edge.val, jj);
    }

    for (auto edge : g.e_succs(ii)) {
      g.add_edge(jj, edge.val, edge.vert);
    }

    potential[jj] = potential[ii];

    CRAB_LOG("zones-sparse", crab::outs() << "After expand " << x << " into "
                                          << y << ":\n"
                                          << *this << "\n");
  }

  void extract(const variable_t &x, linear_constraint_system_t &csts,
               bool only_equalities) {
    crab::CrabStats::count(getDomainName() + ".count.extract");
    crab::ScopedCrabStats __st__(getDomainName() + ".extract");

    normalize();
    if (is_bottom()) {
      return;
    }

    auto it = vert_map.find(x);
    if (it != vert_map.end()) {
      vert_id s = (*it).second;
      if (rev_map[s]) {
        variable_t vs = *rev_map[s];
        SubGraph<graph_t> g_excl(g, 0);
        for (vert_id d : g_excl.verts()) {
          if (rev_map[d]) {
            variable_t vd = *rev_map[d];
            // We give priority to equalities since some domains
            // might not understand inequalities
            if (g_excl.elem(s, d) && g_excl.elem(d, s) &&
                g_excl.edge_val(s, d) == Wt(0) &&
                g_excl.edge_val(d, s) == Wt(0)) {
              linear_constraint_t cst(linear_expression_t(vs) == vd);
              csts += cst;
            } else {
              if (!only_equalities && g_excl.elem(s, d)) {
                linear_constraint_t cst(vd - vs <=
                                        number_t(g_excl.edge_val(s, d)));
                csts += cst;
              }
              if (!only_equalities && g_excl.elem(d, s)) {
                linear_constraint_t cst(vs - vd <=
                                        number_t(g_excl.edge_val(d, s)));
                csts += cst;
              }
            }
          }
        }
      }
    }
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
			  DBM_t invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", getDomainName());    
  }
  /* end intrinsics operations */
  
  // Output function
  void write(crab_os &o) {
    crab::CrabStats::count(getDomainName() + ".count.write");
    crab::ScopedCrabStats __st__(getDomainName() + ".write");

    normalize();

    if (is_bottom()) {
      o << "_|_";
      return;
    } else if (is_top()) {
      o << "{}";
      return;
    } else {
      // Intervals
      bool first = true;
      o << "{";
      // Extract all the edges
      SubGraph<graph_t> g_excl(g, 0);
      for (vert_id v : g_excl.verts()) {
        if (!rev_map[v])
          continue;
        if (!g.elem(0, v) && !g.elem(v, 0))
          continue;
        interval_t v_out = interval_t(g.elem(v, 0) ? -number_t(g.edge_val(v, 0))
                                                   : bound_t::minus_infinity(),
                                      g.elem(0, v) ? number_t(g.edge_val(0, v))
                                                   : bound_t::plus_infinity());

        if (first)
          first = false;
        else
          o << ", ";
        o << *(rev_map[v]) << " -> " << v_out;
      }

      for (vert_id s : g_excl.verts()) {
        if (!rev_map[s])
          continue;
        variable_t vs = *rev_map[s];
        for (vert_id d : g_excl.succs(s)) {
          if (!rev_map[d])
            continue;
          variable_t vd = *rev_map[d];
          if (first)
            first = false;
          else
            o << ", ";
          o << vd << "-" << vs << "<=" << g_excl.edge_val(s, d);
        }
      }
      o << "}";

      // linear_constraint_system_t inv = to_linear_constraint_system();
      // o << inv;
    }
  }

  linear_constraint_system_t to_linear_constraint_system() {
    crab::CrabStats::count(getDomainName() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(getDomainName() +
                                 ".to_linear_constraint_system");

    normalize();

    linear_constraint_system_t csts;

    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    // Extract all the edges

    SubGraph<graph_t> g_excl(g, 0);

    for (vert_id v : g_excl.verts()) {
      if (!rev_map[v])
        continue;
      if (g.elem(v, 0))
        csts += linear_constraint_t(linear_expression_t(*rev_map[v]) >=
                                    -number_t(g.edge_val(v, 0)));
      if (g.elem(0, v))
        csts += linear_constraint_t(linear_expression_t(*rev_map[v]) <=
                                    number_t(g.edge_val(0, v)));
    }

    for (vert_id s : g_excl.verts()) {
      if (!rev_map[s])
        continue;
      variable_t vs = *rev_map[s];
      for (vert_id d : g_excl.succs(s)) {
        if (!rev_map[d])
          continue;
        variable_t vd = *rev_map[d];
        csts += linear_constraint_t(vd - vs <= number_t(g_excl.edge_val(s, d)));
      }
    }

    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  // return number of vertices and edges
  std::pair<std::size_t, std::size_t> size() const {
    return {g.size(), g.num_edges()};
  }

  static std::string getDomainName() { return "SparseDBM"; }

}; // class SparseDBM_

#if 1
template <class Number, class VariableName,
          class Params = DBM_impl::DefaultParams<Number>>
using SparseDBM = SparseDBM_<Number, VariableName, Params>;
#else

template <typename Number, typename VariableName, typename Params>
struct abstract_domain_traits<SparseDBM_<Number, VariableName, Params>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

// Quick wrapper which uses shared references with copy-on-write.
template <class Number, class VariableName,
          class Params = DBM_impl::DefaultParams<Number>>
class SparseDBM final
    : public abstract_domain<SparseDBM<Number, VariableName, Params>> {
  typedef SparseDBM<Number, VariableName, Params> DBM_t;
  typedef abstract_domain<DBM_t> abstract_domain_t;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::pointer_constraint_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  typedef Number number_t;
  typedef VariableName varname_t;
  typedef typename linear_constraint_t::kind_t constraint_kind_t;
  typedef ikos::interval<number_t> interval_t;

public:
  typedef SparseDBM_<number_t, varname_t, Params> dbm_impl_t;
  typedef std::shared_ptr<dbm_impl_t> dbm_ref_t;

  SparseDBM(dbm_ref_t _ref) : norm_ref(_ref) {}

  SparseDBM(dbm_ref_t _base, dbm_ref_t _norm)
      : base_ref(_base), norm_ref(_norm) {}

  DBM_t create(dbm_impl_t &&t) {
    return std::make_shared<dbm_impl_t>(std::move(t));
  }

  DBM_t create_base(dbm_impl_t &&t) {
    dbm_ref_t base = std::make_shared<dbm_impl_t>(t);
    dbm_ref_t norm = std::make_shared<dbm_impl_t>(std::move(t));
    return DBM_t(base, norm);
  }

  void lock(void) {
    // Allocate a fresh copy.
    if (!norm_ref.unique())
      norm_ref = std::make_shared<dbm_impl_t>(*norm_ref);
    base_ref.reset();
  }

public:
  void set_to_top() {
    SparseDBM abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() {
    SparseDBM abs(true);
    std::swap(*this, abs);
  }

  SparseDBM(bool is_bottom = false)
      : norm_ref(std::make_shared<dbm_impl_t>(is_bottom)) {}

  SparseDBM(const DBM_t &o) : base_ref(o.base_ref), norm_ref(o.norm_ref) {}

  SparseDBM &operator=(const DBM_t &o) {
    if (this != &o) {
      base_ref = o.base_ref;
      norm_ref = o.norm_ref;
    }
    return *this;
  }

  dbm_impl_t &base(void) {
    if (base_ref)
      return *base_ref;
    else
      return *norm_ref;
  }
  dbm_impl_t &norm(void) { return *norm_ref; }

  bool is_bottom() { return norm().is_bottom(); }
  bool is_top() { return norm().is_top(); }
  bool operator<=(DBM_t o) { return norm() <= o.norm(); }
  void operator|=(DBM_t o) {
    lock();
    norm() |= o.norm();
  }
  DBM_t operator|(DBM_t o) { return create(norm() | o.norm()); }
  DBM_t operator||(DBM_t o) { return create_base(base() || o.norm()); }
  // DBM_t operator||(DBM_t o) { return create(norm() || o.norm()); }
  DBM_t operator&(DBM_t o) { return create(norm() & o.norm()); }
  DBM_t operator&&(DBM_t o) { return create(norm() && o.norm()); }
  DBM_t widening_thresholds(DBM_t o,
                            const iterators::thresholds<number_t> &ts) {
    return create_base(base().widening_thresholds(o.norm(), ts));
  }

  void normalize() {
    lock();
    norm().normalize();
  }
  void minimize() {}

  void operator+=(linear_constraint_system_t csts) {
    lock();
    norm() += csts;
  }
  void operator-=(variable_t v) {
    lock();
    norm() -= v;
  }
  interval_t operator[](variable_t x) { return norm()[x]; }
  void set(variable_t x, interval_t intv) {
    lock();
    norm().set(x, intv);
  }

  void assign(variable_t x, linear_expression_t e) {
    lock();
    norm().assign(x, e);
  }
  void apply(ikos::operation_t op, variable_t x, variable_t y, number_t k) {
    lock();
    norm().apply(op, x, y, k);
  }
  void apply(ikos::operation_t op, variable_t x, variable_t y, variable_t z) {
    lock();
    norm().apply(op, x, y, z);
  }
  void apply(int_conv_operation_t op, variable_t dst, variable_t src) {
    lock();
    norm().apply(op, dst, src);
  }
  void backward_assign(variable_t x, linear_expression_t e, DBM_t invariant) {
    lock();
    norm().backward_assign(x, e, invariant.norm());
  }
  void backward_apply(operation_t op, variable_t x, variable_t y, number_t k,
                      DBM_t invariant) {
    lock();
    norm().backward_apply(op, x, y, k, invariant.norm());
  }
  void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
                      DBM_t invariant) {
    lock();
    norm().backward_apply(op, x, y, z, invariant.norm());
  }
  void apply(ikos::bitwise_operation_t op, variable_t x, variable_t y,
             number_t k) {
    lock();
    norm().apply(op, x, y, k);
  }
  void apply(ikos::bitwise_operation_t op, variable_t x, variable_t y,
             variable_t z) {
    lock();
    norm().apply(op, x, y, z);
  }

  /* Begin unimplemented operations */
  // boolean operations
  void assign_bool_cst(variable_t lhs, linear_constraint_t rhs) {}
  void assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs) {}
  void apply_binary_bool(bool_operation_t op, variable_t x, variable_t y,
                         variable_t z) {}
  void assume_bool(variable_t v, bool is_negated) {}
  // backward boolean operations
  void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                DBM_t invariant) {}
  void backward_assign_bool_var(variable_t lhs, variable_t rhs, bool is_not_rhs,
                                DBM_t invariant) {}
  void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                  variable_t y, variable_t z, DBM_t invariant) {
  }
  // array operations
  void array_init(variable_t a, linear_expression_t elem_size,
                  linear_expression_t lb_idx, linear_expression_t ub_idx,
                  linear_expression_t val) {}
  void array_load(variable_t lhs, variable_t a, linear_expression_t elem_size,
                  linear_expression_t i) {}
  void array_store(variable_t a, linear_expression_t elem_size,
                   linear_expression_t i, linear_expression_t v,
                   bool is_strong_update) {}
  void array_store(variable_t a_new, variable_t a_old,
                   linear_expression_t elem_size, linear_expression_t i,
                   linear_expression_t v, bool is_strong_update) {}
  void array_store_range(variable_t a, linear_expression_t elem_size,
                         linear_expression_t i, linear_expression_t j,
                         linear_expression_t v) {}
  void array_store_range(variable_t a_new, variable_t a_old,
                         linear_expression_t elem_size, linear_expression_t i,
                         linear_expression_t j, linear_expression_t v) {}
  void array_assign(variable_t lhs, variable_t rhs) {}
  // backward array operations
  void backward_array_init(variable_t a, linear_expression_t elem_size,
                           linear_expression_t lb_idx,
                           linear_expression_t ub_idx, linear_expression_t val,
                           DBM_t invariant) {}
  void backward_array_load(variable_t lhs, variable_t a,
                           linear_expression_t elem_size, linear_expression_t i,
                           DBM_t invariant) {}
  void backward_array_store(variable_t a, linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, DBM_t invariant) {}
  void backward_array_store(variable_t a_new, variable_t a_old,
                            linear_expression_t elem_size,
                            linear_expression_t i, linear_expression_t v,
                            bool is_strong_update, DBM_t invariant) {}
  void backward_array_store_range(variable_t a, linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v, DBM_t invariant) {}
  void backward_array_store_range(variable_t a_new, variable_t a_old,
                                  linear_expression_t elem_size,
                                  linear_expression_t i, linear_expression_t j,
                                  linear_expression_t v, DBM_t invariant) {}
  void backward_array_assign(variable_t lhs, variable_t rhs, DBM_t invariant) {}
  // pointer operations
  void pointer_load(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_store(variable_t lhs, variable_t rhs, linear_expression_t elem_size) {}
  void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) {}
  void pointer_mk_obj(variable_t lhs, ikos::index_t address) {}
  void pointer_function(variable_t lhs, varname_t func) {}
  void pointer_mk_null(variable_t lhs) {}
  void pointer_assume(pointer_constraint_t cst) {}
  void pointer_assert(pointer_constraint_t cst) {}
  /* End unimplemented operations */

  void rename(const variable_vector_t &from, const variable_vector_t &to) {
    lock();
    norm().rename(from, to);
  }

  void forget(const variable_vector_t &variables) {
    lock();
    norm().forget(variables);
  }

  void expand(variable_t x, variable_t y) {
    lock();
    norm().expand(x, y);
  }

  void project(const variable_vector_t &variables) {
    lock();
    norm().project(variables);
  }

  void extract(const variable_t &x, linear_constraint_system_t &csts,
               bool only_equalities) {
    norm().extract(x, csts, only_equalities);
  }

  /* begin intrinsics operations */  
  void intrinsic(std::string name,
		 const variable_vector_t &inputs,
		 const variable_vector_t &outputs) override {
    lock();
    norm().intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
			  const variable_vector_t &inputs,
			  const variable_vector_t &outputs,
			  DBM_t invariant) override {
    lock();
    norm().backward_intrinsic(name, inputs, outputs, invariant);
  }
  /* end intrinsics operations */
  
  void write(crab_os &o) { norm().write(o); }

  linear_constraint_system_t to_linear_constraint_system() {
    return norm().to_linear_constraint_system();
  }
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    return norm().to_disjunctive_linear_constraint_system();
  }

  static std::string getDomainName() { return dbm_impl_t::getDomainName(); }

  std::pair<std::size_t, std::size_t> size() const { return norm().size(); }

protected:
  dbm_ref_t base_ref;
  dbm_ref_t norm_ref;
};
#endif

template <typename Number, typename VariableName, typename Params>
struct abstract_domain_traits<SparseDBM<Number, VariableName, Params>> {
  typedef Number number_t;
  typedef VariableName varname_t;
};

template <typename Number, typename VariableName, typename Params>
class reduced_domain_traits<SparseDBM<Number, VariableName, Params>> {
public:
  typedef SparseDBM<Number, VariableName, Params> sdbm_domain_t;
  typedef typename sdbm_domain_t::variable_t variable_t;
  typedef typename sdbm_domain_t::linear_constraint_system_t
      linear_constraint_system_t;

  static void extract(sdbm_domain_t &dom, const variable_t &x,
                      linear_constraint_system_t &csts, bool only_equalities) {
    dom.extract(x, csts, only_equalities);
  }
};

} // namespace domains

} // namespace crab

#pragma GCC diagnostic pop
