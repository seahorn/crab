/*******************************************************************************
 *
 * Difference Bound Matrix domain based on the paper "Exploiting
 * Sparsity in Difference-Bound Matrices" by Gange, Navas, Schachte,
 * Sondergaard, and Stuckey published in SAS'16.
 *
 * A re-engineered implementation of the Difference Bound Matrix
 * domain, which maintains bounds and relations separately.
 *
 * Closure operations based on the paper "Fast and Flexible Difference
 * Constraint Propagation for DPLL(T)" by Cotton and Maler.
 *
 * Author: Graeme Gange (gkgange@unimelb.edu.au)
 *
 * Contributors: Jorge A. Navas (jorge.navas@sri.com)
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/graphs/graph_config.hpp>
#include <crab/domains/graphs/graph_ops.hpp>
#include <crab/domains/graphs/graph_views.hpp>
#include <crab/domains/interval.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <boost/optional.hpp>
#include <unordered_set>

#define JOIN_CLOSE_AFTER_MEET
//#define CHECK_POTENTIAL
//#define SDBM_NO_NORMALIZE
#define USE_FLAT_MAP

#ifdef USE_FLAT_MAP
#include <boost/container/flat_map.hpp>
#else
// Operations like rename are much faster using unordered_map
#include <unordered_map>
#endif

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {

namespace domains {

template <class Number, class VariableName,
          class Params = DBM_impl::DefaultParams<Number>>
class split_dbm_domain final
    : public abstract_domain_api<
          split_dbm_domain<Number, VariableName, Params>> {
  using DBM_t = split_dbm_domain<Number, VariableName, Params>;
  using abstract_domain_t = abstract_domain_api<DBM_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;    
  using number_t = Number;
  using varname_t = VariableName;
  using constraint_kind_t = typename linear_constraint_t::kind_t;

private:
  using bound_t = ikos::bound<number_t>;
  using Wt = typename Params::Wt;
  using graph_t = typename Params::graph_t;
  using ntow = DBM_impl::NtoW<number_t, Wt>;
  using vert_id = typename graph_t::vert_id;
  using wt_ref_t = typename graph_t::wt_ref_t;
#ifdef USE_FLAT_MAP
  using vert_map_t = boost::container::flat_map<variable_t, vert_id>;
#else
  using vert_map_t = std::unordered_map<variable_t, vert_id>;
#endif
  using vmap_elt_t = typename vert_map_t::value_type;
  using rev_map_t = std::vector<boost::optional<variable_t>>;
  using GrOps = GraphOps<graph_t>;
  using GrPerm = GraphPerm<graph_t>;
  using edge_vector = typename GrOps::edge_vector;
  // < <x, y>, k> == x - y <= k.
  using diffcst_t = std::pair<std::pair<variable_t, variable_t>, Wt>;
  using vert_set_t = std::unordered_set<vert_id>;

protected:
  //================
  // Domain data
  //================
  // GKG: ranges are now maintained in the graph
  vert_map_t vert_map; // Mapping from variables to vertices
  rev_map_t rev_map;
  graph_t g;                 // The underlying relation graph
  std::vector<Wt> potential; // Stored potential for the vertex
  vert_set_t unstable;
  bool _is_bottom;

  class Wt_max {
  public:
    Wt_max() {}
    Wt apply(const Wt &x, const Wt &y) { return std::max(x, y); }
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
    vert_map.insert(vmap_elt_t(v, vert));
    // Initialize
    assert(vert <= rev_map.size());
    if (vert < rev_map.size()) {
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

  boost::optional<vert_id> get_vert(const variable_t &v) const {
    auto it = vert_map.find(v);
    if (it != vert_map.end()) {
      return (*it).second;
    } else {
      return boost::none;
    } 
  }
  
  vert_id get_vert(graph_t &g, vert_map_t &vmap, rev_map_t &rmap,
                   std::vector<Wt> &pot, variable_t v) {
    auto it = vmap.find(v);
    if (it != vmap.end())
      return (*it).second;

    vert_id vert(g.new_vertex());
    vmap.insert(vmap_elt_t(v, vert));
    // Initialize
    assert(vert <= rmap.size());
    if (vert < rmap.size()) {
      pot[vert] = Wt(0);
      rmap[vert] = v;
    } else {
      pot.push_back(Wt(0));
      rmap.push_back(v);
    }
    vmap.insert(vmap_elt_t(v, vert));

    return vert;
  }

  template <class G, class P>
  inline void check_potential(const G &g, const P &p, unsigned line) const {
#ifdef CHECK_POTENTIAL
    for (vert_id v : g.verts()) {
      for (vert_id d : g.succs(v)) {
        if (p[v] + g.edge_val(v, d) - p[d] < Wt(0)) {
          CRAB_ERROR("Invalid potential at line ", line, ":", "pot[", v,
                     "]=", p[v], " ", "pot[", d, "]=", p[d], " ", "edge(", v,
                     ",", d, ")=", g.edge_val(v, d));
        }
      }
    }
#endif
  }

  class vert_set_wrap_t {
  public:
    vert_set_wrap_t(const vert_set_t &_vs) : vs(_vs) {}

    bool operator[](vert_id v) const { return vs.find(v) != vs.end(); }
    const vert_set_t &vs;
  };

  // Evaluate the potential value of a variable.
  Wt pot_value(const variable_t &v) {
    auto it = vert_map.find(v);
    if (it != vert_map.end())
      return potential[(*it).second];
    return ((Wt)0);
  }

  // Evaluate an expression under the chosen potentials
  Wt eval_expression(const linear_expression_t &e, bool &overflow) {
    overflow = false;
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

  interval_t eval_interval(const linear_expression_t &e) {
    interval_t r = e.constant();
    for (auto p : e)
      r += p.first * operator[](p.second);
    return r;
  }

  interval_t compute_residual(const linear_expression_t &e, variable_t pivot) {
    interval_t residual(-e.constant());
    for (auto kv : e) {
      const variable_t &v = kv.second;
      if (v.index() != pivot.index()) {
        residual = residual - (interval_t(kv.first) * this->operator[](v));
      }
    }
    return residual;
  }

  /**
   *  Turn an assignment into a set of difference constraints.
   *
   *  Given v := a*x + b*y + k, where a,b >= 0, we generate the
   *  difference constraints:
   *
   *  if extract_upper_bounds
   *     v - x <= ub((a-1)*x + b*y + k)
   *     v - y <= ub(a*x + (b-1)*y + k)
   *  else
   *     x - v <= lb((a-1)*x + b*y + k)
   *     y - v <= lb(a*x + (b-1)*y + k)
   **/
  void diffcsts_of_assign(const variable_t &x, const linear_expression_t &exp,
                          /* if true then process the upper
                             bounds, else the lower bounds */
                          bool extract_upper_bounds,
                          /* foreach {v, k} \in diff_csts we have
                             the difference constraint v - k <= k */
                          std::vector<std::pair<variable_t, Wt>> &diff_csts) {

    boost::optional<variable_t> unbounded_var;
    std::vector<std::pair<variable_t, Wt>> terms;
    bool overflow;

    Wt residual(ntow::convert(exp.constant(), overflow));
    if (overflow) {
      return;
    }

    for (auto p : exp) {
      Wt coeff(ntow::convert(p.first, overflow));
      if (overflow) {
        continue;
      }

      variable_t y(p.second);
      if (coeff < Wt(0)) {
        // Can't do anything with negative coefficients.
        bound_t y_val =
            (extract_upper_bounds ? operator[](y).lb() : operator[](y).ub());

        if (y_val.is_infinite()) {
          return;
        }
        residual += ntow::convert(*(y_val.number()), overflow) * coeff;
        if (overflow) {
          continue;
        }

      } else {
        bound_t y_val =
            (extract_upper_bounds ? operator[](y).ub() : operator[](y).lb());

        if (y_val.is_infinite()) {
          if (unbounded_var || coeff != Wt(1)) {
            return;
          }
          unbounded_var = y;
        } else {
          Wt ymax(ntow::convert(*(y_val.number()), overflow));
          if (overflow) {
            continue;
          }
          residual += ymax * coeff;
          terms.push_back({y, ymax});
        }
      }
    }

    if (unbounded_var) {
      // There is exactly one unbounded variable with unit
      // coefficient
      diff_csts.push_back({*unbounded_var, residual});
    } else {
      for (auto p : terms) {
        diff_csts.push_back({p.first, residual - p.second});
      }
    }
  }

  // Turn an assignment into a set of difference constraints.
  void diffcsts_of_assign(const variable_t &x, const linear_expression_t &exp,
                          std::vector<std::pair<variable_t, Wt>> &lb,
                          std::vector<std::pair<variable_t, Wt>> &ub) {
    diffcsts_of_assign(x, exp, true, ub);
    diffcsts_of_assign(x, exp, false, lb);
  }

  /**
   * Turn a linear inequality into a set of difference
   * constraints.
   **/
  void diffcsts_of_lin_leq(const linear_expression_t &exp,
                           /* difference contraints */
                           std::vector<diffcst_t> &csts,
                           /* x >= lb for each {x,lb} in lbs */
                           std::vector<std::pair<variable_t, Wt>> &lbs,
                           /* x <= ub for each {x,ub} in ubs */
                           std::vector<std::pair<variable_t, Wt>> &ubs) const {

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
        bound_t y_lb = at(y).lb();
        if (y_lb.is_infinite()) {
          if (unbounded_lbvar) {
            return;
          }
          unbounded_lbvar = y;
          unbounded_lbcoeff = coeff;
        } else {
          Wt ymin(ntow::convert(*(y_lb.number()), overflow));
          if (overflow) {
            continue;
          }
          exp_ub -= ymin * coeff;
          pos_terms.push_back({{coeff, y}, ymin});
        }
      } else {
        variable_t y(p.second);
        bound_t y_ub = at(y).ub();
        if (y_ub.is_infinite()) {
          if (unbounded_ubvar) {
            return;
          }
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
        if (unbounded_lbcoeff != Wt(1) || unbounded_ubcoeff != Wt(1)) {
          return;
        }
        variable_t y(*unbounded_ubvar);
        csts.push_back({{x, y}, exp_ub});
      } else {
        if (unbounded_lbcoeff == Wt(1)) {
          for (auto p : neg_terms) {
            csts.push_back({{x, p.first.second}, exp_ub - p.second});
          }
        }
        // Add bounds for x
        ubs.push_back({x, exp_ub / unbounded_lbcoeff});
      }
    } else {
      if (unbounded_ubvar) {
        variable_t y(*unbounded_ubvar);
        if (unbounded_ubcoeff == Wt(1)) {
          for (auto p : pos_terms) {
            csts.push_back({{p.first.second, y}, exp_ub + p.second});
          }
        }
        // Add bounds for y
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
  }

  bool add_linear_leq(const linear_expression_t &exp) {
    CRAB_LOG("zones-split", linear_expression_t exp_tmp(exp);
             crab::outs() << "Adding: " << exp_tmp << "<= 0"
                          << "\n");
    std::vector<std::pair<variable_t, Wt>> lbs, ubs;
    std::vector<diffcst_t> csts;
    diffcsts_of_lin_leq(exp, csts, lbs, ubs);

    check_potential(g, potential, __LINE__);

    Wt_min min_op;
    wt_ref_t w;
    for (auto p : lbs) {
      CRAB_LOG("zones-split", crab::outs()
                                  << p.first << ">=" << p.second << "\n");
      variable_t x(p.first);
      vert_id v = get_vert(p.first);
      if (g.lookup(v, 0, w) && w.get() <= -p.second)
        continue;
      g.set_edge(v, -p.second, 0);

      if (!repair_potential(v, 0)) {
        set_to_bottom();
        return false;
      }
      check_potential(g, potential, __LINE__);
      // Compute other updated bounds
      if (crab_domain_params_man::get().zones_close_bounds_inline()) {
        for (auto e : g.e_preds(v)) {
          if (e.vert == 0)
            continue;
          g.update_edge(e.vert, e.val - p.second, 0, min_op);

          if (!repair_potential(e.vert, 0)) {
            set_to_bottom();
            return false;
          }
          check_potential(g, potential, __LINE__);
        }
      }
    }
    for (auto p : ubs) {
      CRAB_LOG("zones-split", crab::outs()
                                  << p.first << "<=" << p.second << "\n");
      variable_t x(p.first);
      vert_id v = get_vert(p.first);
      if (g.lookup(0, v, w) && w.get() <= p.second)
        continue;
      g.set_edge(0, p.second, v);
      if (!repair_potential(0, v)) {
        set_to_bottom();
        return false;
      }
      check_potential(g, potential, __LINE__);

      if (crab_domain_params_man::get().zones_close_bounds_inline()) {
        for (auto e : g.e_succs(v)) {
          if (e.vert == 0)
            continue;
          g.update_edge(0, e.val + p.second, e.vert, min_op);
          if (!repair_potential(0, e.vert)) {
            set_to_bottom();
            return false;
          }
          check_potential(g, potential, __LINE__);
        }
      }
    }

    for (auto diff : csts) {
      CRAB_LOG("zones-split", crab::outs() << diff.first.first << "-"
                                           << diff.first.second
                                           << "<=" << diff.second << "\n");

      vert_id src = get_vert(diff.first.second);
      vert_id dest = get_vert(diff.first.first);

      // Check if the edge (src,dest) via bounds already exists
      wt_ref_t w1, w2;
      if (g.lookup(src, 0, w1) && g.lookup(0, dest, w2) &&
          (w1.get() + w2.get()) <= diff.second) {
        continue;
      }

      g.update_edge(src, diff.second, dest, min_op);
      if (!repair_potential(src, dest)) {
        set_to_bottom();
        return false;
      }
      check_potential(g, potential, __LINE__);
      close_over_edge(src, dest);
      check_potential(g, potential, __LINE__);
    }
    // Collect bounds
    // GKG: Now done in close_over_edge

    if (!crab_domain_params_man::get().zones_close_bounds_inline()) {
      edge_vector delta;
      GrOps::close_after_assign(g, potential, 0, delta);
      GrOps::apply_delta(g, delta);
    }

    check_potential(g, potential, __LINE__);
    // CRAB_WARN("split_dbm_domain::add_linear_leq not yet implemented.");
    return true;
  }

  // x != n
  bool add_univar_disequation(const variable_t &x, number_t n) {
    bool overflow;
    interval_t i = get_interval(x);
    interval_t ni(n);
    interval_t new_i =
        ikos::linear_interval_solver_impl::trim_interval<interval_t>(
            i, ni);
    if (new_i.is_bottom()) {
      set_to_bottom();
      return false;
    } else if (!new_i.is_top() && (new_i <= i)) {
      vert_id v = get_vert(x);
      wt_ref_t w; 
      Wt_min min_op;
      if (new_i.lb().is_finite()) {
        // strenghten lb
        Wt lb_val = ntow::convert(-(*(new_i.lb().number())), overflow);
        if (overflow) {
          return true;
        }

        if (g.lookup(v, 0, w) && lb_val < w.get()) {
          g.set_edge(v, lb_val, 0);
          if (!repair_potential(v, 0)) {
            set_to_bottom();
            return false;
          }
          check_potential(g, potential, __LINE__);
          // Update other bounds
          for (auto e : g.e_preds(v)) {
            if (e.vert == 0)
              continue;
            g.update_edge(e.vert, e.val + lb_val, 0, min_op);
            if (!repair_potential(e.vert, 0)) {
              set_to_bottom();
              return false;
            }
            check_potential(g, potential, __LINE__);
          }
        }
      }
      if (new_i.ub().is_finite()) {
        // strengthen ub
        Wt ub_val = ntow::convert(*(new_i.ub().number()), overflow);
        if (overflow) {
          return true;
        }

        if (g.lookup(0, v, w) && (ub_val < w.get())) {
          g.set_edge(0, ub_val, v);
          if (!repair_potential(0, v)) {
            set_to_bottom();
            return false;
          }
          check_potential(g, potential, __LINE__);
          // Update other bounds
          for (auto e : g.e_succs(v)) {
            if (e.vert == 0)
              continue;
            g.update_edge(0, e.val + ub_val, e.vert, min_op);
            if (!repair_potential(0, e.vert)) {
              set_to_bottom();
              return false;
            }
            check_potential(g, potential, __LINE__);
          }
        }
      }
    }
    return true;
  }

  void add_disequation(const linear_expression_t &e) {
    // XXX: similar precision as the interval domain

    for (auto kv : e) {
      const variable_t &pivot = kv.second;
      interval_t i = compute_residual(e, pivot) / interval_t(kv.first);
      if (auto k = i.singleton()) {
        if (!add_univar_disequation(pivot, *k)) {
          // set_to_bottom() was already called
          return;
        }
      }
    }

    /*
    // Can only exploit \sum_i c_i x_i \neq k if:
    // (1) exactly one x_i is unfixed
    // (2) lb(x_i) or ub(x_i) = k - \sum_i' c_i' x_i'
    Wt k = exp.constant();
    auto it = exp.begin();
    for(; it != exp.end(); ++it)
    {
      if(!var_is_fixed((*it).second))
        break;
      k -= (*it).first*get_value((*it).second);
    }

    // All variables are fixed
    if(it == exp.end())
    {
      if(k == Wt(0))
        set_to_bottom();
      return;
    }

    // Found one unfixed variable; collect the rest.
    Wt ucoeff = (*it).first;
    varname_t uvar((*it).second;
    interval_t u_int = get_interval(ranges, uvar);
    // We need at least one side of u to be finite.
    if(u_int.lb().is_infinite() && u_int.ub().is_infinite())
      return;

    for(++it; it != exp.end(); ++it)
    {
      // Two unfixed variables; nothing we can do.
      if(!var_is_fixed((*it).second))
        return;
      k -= (*it).first*get_value((*it).second);
    }
    */
  }

  interval_t get_interval(const variable_t &x) const {
    return get_interval(vert_map, g, x);
  }

  interval_t get_interval(const vert_map_t &m, const graph_t &g,
                          const variable_t &x) const {
    auto it = m.find(x);
    if (it == m.end()) {
      return interval_t::top();
    }
    vert_id v = (*it).second;
    interval_t x_out = interval_t(
        g.elem(v, 0) ? -number_t(g.edge_val(v, 0)) : bound_t::minus_infinity(),
        g.elem(0, v) ? number_t(g.edge_val(0, v)) : bound_t::plus_infinity());
    return x_out;
  }

  // Resore potential after an edge addition
  bool repair_potential(vert_id src, vert_id dest) {
    return GrOps::repair_potential(g, potential, src, dest);
  }

  // Restore closure after a single edge addition
  void close_over_edge(vert_id ii, vert_id jj) {
    assert(ii != 0 && jj != 0);

    Wt_min min_op;
    wt_ref_t w;
    SubGraph<graph_t> g_excl(g, 0);
    Wt c = g_excl.edge_val(ii, jj);


    if (crab_domain_params_man::get().zones_close_bounds_inline()) {
      if (g.lookup(0, ii, w))
        g.update_edge(0, w.get() + c, jj, min_op);
      if (g.lookup(jj, 0, w))
        g.update_edge(ii, w.get() + c, 0, min_op);
    }

    // There may be a cheaper way to do this.
    // GKG: Now implemented.
    std::vector<std::pair<vert_id, Wt>> src_dec;
    // we add in delta so that we don't invalidate graph iterators
    edge_vector delta;
    for (auto edge : g_excl.e_preds(ii)) {
      vert_id se = edge.vert;
      Wt wt_sij = edge.val + c;
      assert(g_excl.succs(se).begin() != g_excl.succs(se).end());
      if (se != jj) {
        if (g_excl.lookup(se, jj, w)) {
          if (w.get() <= wt_sij) {
            continue;
	  }
	  // REVISIT(PERFORMANCE): extra call to lookup
	  g.set_edge(se, wt_sij, jj);
        } else {
          delta.push_back({{se, jj}, wt_sij});
        }
        src_dec.push_back(std::make_pair(se, edge.val));
        if (crab_domain_params_man::get().zones_close_bounds_inline()) {
          if (g.lookup(0, se, w))
            g.update_edge(0, w.get() + wt_sij, jj, min_op);
          if (g.lookup(jj, 0, w))
            g.update_edge(se, w.get() + wt_sij, 0, min_op);
        }
      }
    }
    GrOps::apply_delta(g, delta);
    delta.clear();
    std::vector<std::pair<vert_id, Wt>> dest_dec;
    for (auto edge : g_excl.e_succs(jj)) {
      vert_id de = edge.vert;
      Wt wt_ijd = edge.val + c;
      if (de != ii) {
        if (g_excl.lookup(ii, de, w)) {
          if (w.get() <= wt_ijd) {
            continue;
	  }
	  // REVISIT(PERFORMANCE): extra call to lookup
	  g.set_edge(ii, wt_ijd, de);
        } else {
          delta.push_back({{ii, de}, {wt_ijd}});
        }
        dest_dec.push_back(std::make_pair(de, edge.val));
        if (crab_domain_params_man::get().zones_close_bounds_inline()) {
          if (g.lookup(0, ii, w))
            g.update_edge(0, w.get() + wt_ijd, de, min_op);
          if (g.lookup(de, 0, w))
            g.update_edge(ii, w.get() + wt_ijd, 0, min_op);
        }
      }
    }
    GrOps::apply_delta(g, delta);

    for (auto s_p : src_dec) {
      vert_id se = s_p.first;
      Wt wt_sij = c + s_p.second;
      for (auto d_p : dest_dec) {
        vert_id de = d_p.first;
        Wt wt_sijd = wt_sij + d_p.second;
        if (g.lookup(se, de, w)) {
          if (w.get() <= wt_sijd) {
            continue;
	  }
	  // REVISIT(PERFORMANCE): extra call to lookup
	  g.set_edge(se, wt_sijd, de);
        } else {
          g.add_edge(se, wt_sijd, de);
        }
        if (crab_domain_params_man::get().zones_close_bounds_inline()) {
          if (g.lookup(0, se, w))
            g.update_edge(0, w.get() + wt_sijd, de, min_op);
          if (g.lookup(de, 0, w))
            g.update_edge(se, w.get() + wt_sijd, 0, min_op);
        }
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

  // return true if edge from x to y with weight k is unsatisfiable
  bool is_unsat_edge(vert_id x, vert_id y, Wt k) const {
    wt_ref_t w;
    if (g.lookup(y, x, w)) {
      return ((w.get() + k) < Wt(0));
    } else {
      interval_t intv_x = interval_t::top();
      interval_t intv_y = interval_t::top();
      if (g.elem(0, x) || g.elem(x, 0)) {
        intv_x = interval_t(g.elem(x, 0) ? -number_t(g.edge_val(x, 0))
                                         : bound_t::minus_infinity(),
                            g.elem(0, x) ? number_t(g.edge_val(0, x))
                                         : bound_t::plus_infinity());
      }
      if (g.elem(0, y) || g.elem(y, 0)) {
        intv_y = interval_t(g.elem(y, 0) ? -number_t(g.edge_val(y, 0))
                                         : bound_t::minus_infinity(),
                            g.elem(0, y) ? number_t(g.edge_val(0, y))
                                         : bound_t::plus_infinity());
      }
      if (intv_x.is_top() || intv_y.is_top()) {
        return false;
      } else {
        return (!((intv_y - intv_x).lb() <= (number_t)k));
      }
    }
  }

  // return true iff cst is unsatisfiable without modifying the DBM
  bool is_unsat(const linear_constraint_t &cst) const {
    if (is_bottom() || cst.is_contradiction()) {
      return true;
    }

    if (is_top() || cst.is_tautology()) {
      return false;
    }

    std::vector<std::pair<variable_t, Wt>> lbs, ubs;
    std::vector<diffcst_t> diffcsts;

    if (cst.is_inequality()) {
      linear_expression_t exp = cst.expression();
      diffcsts_of_lin_leq(exp, diffcsts, lbs, ubs);
    } else if (cst.is_strict_inequality()) {
      auto nc = ikos::linear_constraint_impl::strict_to_non_strict_inequality(cst);
      if (nc.is_inequality()) {
        linear_expression_t exp = nc.expression();
        diffcsts_of_lin_leq(exp, diffcsts, lbs, ubs);
      } else {
        // we couldn't convert the strict into a non-strict
        return false;
      }
    } else if (cst.is_equality()) {
      linear_expression_t exp = cst.expression();
      diffcsts_of_lin_leq(exp, diffcsts, lbs, ubs);
      diffcsts_of_lin_leq(-exp, diffcsts, lbs, ubs);
    } else if (cst.is_disequation()) {
      CRAB_WARN("disequalities ", cst, " not implemented by ", domain_name(), "::is_unsat");      
      return false;
    } else {
      return false;
    }

    // check difference constraints
    for (auto diffcst : diffcsts) {
      variable_t x = diffcst.first.first;
      variable_t y = diffcst.first.second;
      Wt k = diffcst.second;

      auto vy = get_vert(y);
      auto vx = get_vert(x);
      if (vx && vy && is_unsat_edge(*vy, *vx, k)) {
        return true;
      }
    }

    // check interval constraints
    for (auto ub : ubs) {
      auto vx = get_vert(ub.first);
      if (vx && is_unsat_edge(0, *vx , ub.second)) {
        return true;
      }
    }
    for (auto lb : lbs) {
      auto vx = get_vert(lb.first);
      if (vx && is_unsat_edge(*vx, 0, -lb.second)) {
        return true;
      }
    }

    return false;
  }
  
  // Join of gx and gy.
  graph_t join(GrPerm &gx, GrPerm &gy, unsigned sz, std::vector<Wt> &pot_rx,
               std::vector<Wt> &pot_ry) const {

    // Compute the deferred relations
    graph_t g_ix_ry;
    wt_ref_t ws, wd;    
    g_ix_ry.growTo(sz);
    SubGraph<GrPerm> gy_excl(gy, 0);
    for (vert_id s : gy_excl.verts()) {
      for (vert_id d : gy_excl.succs(s)) {
        if (gx.lookup(s, 0, ws) && gx.lookup(0, d, wd)) {
          g_ix_ry.add_edge(s, ws.get() + wd.get(), d);
        }
      }
    }
    // Apply the deferred relations, and re-close.
    bool is_closed;
    graph_t g_rx(GrOps::meet(gx, g_ix_ry, is_closed));
    check_potential(g_rx, pot_rx, __LINE__);
#ifdef JOIN_CLOSE_AFTER_MEET
    // Conjecture: g_rx is closed
    if (!is_closed) {
      edge_vector delta;
      SubGraph<graph_t> g_rx_excl(g_rx, 0);
      GrOps::close_after_meet(g_rx_excl, pot_rx, gx, g_ix_ry, delta);
      GrOps::apply_delta(g_rx, delta);
    }
#endif

    graph_t g_rx_iy;
    g_rx_iy.growTo(sz);
    SubGraph<GrPerm> gx_excl(gx, 0);
    for (vert_id s : gx_excl.verts()) {
      for (vert_id d : gx_excl.succs(s)) {
        // Assumption: gx.mem(s, d) -> gx.edge_val(s, d) <=
        //             ranges[var(s)].ub() - ranges[var(d)].lb()
        // That is, if the relation exists, it's at least as strong as the
        // bounds.
        if (gy.lookup(s, 0, ws) && gy.lookup(0, d, wd))
          g_rx_iy.add_edge(s, ws.get() + wd.get(), d);
      }
    }
    // Similarly, should use a SubGraph view.
    graph_t g_ry(GrOps::meet(gy, g_rx_iy, is_closed));
    check_potential(g_ry, pot_ry, __LINE__);
#ifdef JOIN_CLOSE_AFTER_MEET
    // Conjecture: g_ry is closed
    if (!is_closed) {
      edge_vector delta;
      SubGraph<graph_t> g_ry_excl(g_ry, 0);
      GrOps::close_after_meet(g_ry_excl, pot_ry, gy, g_rx_iy, delta);
      GrOps::apply_delta(g_ry, delta);
    }
#endif

    // We now have the relevant set of relations. Because g_rx
    // and g_ry are closed, the result is also closed.
    Wt_min min_op;
    graph_t join_g(GrOps::join(g_rx, g_ry));

    // Now reapply the missing independent relations.
    // Need to derive vert_ids from lb_up/lb_down, and make sure the vertices
    // exist
    std::vector<vert_id> lb_up;
    std::vector<vert_id> lb_down;
    std::vector<vert_id> ub_up;
    std::vector<vert_id> ub_down;

    wt_ref_t wx, wy;
    for (vert_id v : gx_excl.verts()) {
      if (gx.lookup(0, v, wx) && gy.lookup(0, v, wy)) {
        if (wx.get() < wy.get())
          ub_up.push_back(v);
        if (wy.get() < wx.get())
          ub_down.push_back(v);
      }
      if (gx.lookup(v, 0, wx) && gy.lookup(v, 0, wy)) {
        if (wx.get() < wy.get())
          lb_down.push_back(v);
        if (wy.get() < wx.get())
          lb_up.push_back(v);
      }
    }

    for (vert_id s : lb_up) {
      Wt dx_s = gx.edge_val(s, 0);
      Wt dy_s = gy.edge_val(s, 0);
      for (vert_id d : ub_up) {
        if (s == d)
          continue;
        join_g.update_edge(
            s, std::max(dx_s + gx.edge_val(0, d), dy_s + gy.edge_val(0, d)), d,
            min_op);
      }
    }

    for (vert_id s : lb_down) {
      Wt dx_s = gx.edge_val(s, 0);
      Wt dy_s = gy.edge_val(s, 0);
      for (vert_id d : ub_down) {
        if (s == d)
          continue;
        join_g.update_edge(
            s, std::max(dx_s + gx.edge_val(0, d), dy_s + gy.edge_val(0, d)), d,
            min_op);
      }
    }
    return join_g;
  }

  template <class G1, class G2>
  graph_t split_widen(G1 &l, G2 &r, std::vector<vert_id> &unstable) const {
    assert(l.size() == r.size());
    size_t sz = l.size();
    graph_t g;
    g.growTo(sz);

    Wt_min min_op;
    wt_ref_t wx, wy;

    auto update_edge_widen_g = [&g](vert_id src, vert_id dst, Wt val) {
      wt_ref_t wz;
      if (g.lookup(src, dst, wz)) {
        if (wz.get() > val) {
          g.set_edge(src, val, dst);
          return true;
        }
      } else {
        g.add_edge(src, val, dst);
        return true;
      }
      return false;
    };

    /**
     * Check for stable implicit relationships in r
     **/
    for (auto edge_pred : r.e_preds(0)) {
      vert_id s = edge_pred.vert;
      for (auto edge_succ : r.e_succs(0)) {
        vert_id d = edge_succ.vert;
        if (s == d)
          continue;
        /* for each edge(s,0,d) in r check if exists edge(s,d) in l */
        if (l.lookup(s, d, wx) &&
            ((edge_pred.val + edge_succ.val) <= wx.get())) {
          bool res = update_edge_widen_g(s, d, wx.get());
          if (res) {
            CRAB_LOG("zones-split-widening", auto vs = rev_map[s];
                     auto vd = rev_map[d];
                     crab::outs() << "Widening 1: added " << *vd << "-" << *vs
		                  << "<=" << wx.get() << "\n";);
          }
        }
      }
    }

    /**
     * Check for stable explicit relationships in r
     **/
    for (vert_id s : r.verts()) {
      for (auto e : r.e_succs(s)) {
        vert_id d = e.vert;
        /* for each edge(s,d) in r check if exists edge(s,d) in l */
        if (l.lookup(s, d, wx) && e.val <= wx.get()) {
          bool res = update_edge_widen_g(s, d, wx.get());
          if (res) {
            CRAB_LOG(
                "zones-split-widening", auto vs = rev_map[s];
                auto vd = rev_map[d]; if (s == 0 && d != 0) {
                  crab::outs() << "Widening 2: added " << *vd
                               << "<=" << wx.get() << "\n";
                } else if (s != 0 && d == 0) {
                  crab::outs() << "Widening 2: added "
                               << "-" << *vs << "<=" << wx.get() << "\n";
                } else {
                  crab::outs() << "Widening 2: added " << *vd << "-" << *vs
                               << "<=" << wx.get() << "\n";
                });
          }
        }
      }

      // Check if this vertex is stable
      for (vert_id d : l.succs(s)) {
        if (!g.elem(s, d)) {
          unstable.push_back(s);
          CRAB_LOG("zones-split-widening",
                   if (s == 0) {
                     crab::outs() << "Widening 5: added v0"
                                  << " in the normalization queue\n";
                   } else {
                     auto vs = rev_map[s];
                     crab::outs() << "Widening 5: added " << *vs
                                  << " in the normalization queue\n";
                   });
          break;
        }
      }
    }

    // for(vert_id s: r.verts()) {
    //   for(vert_id d : l.succs(s)) {
    //     if(!g.elem(s, d)) {
    //       unstable.push_back(s);
    //       break;
    //     }
    //   }
    // }
    return g;
  }

  bool need_normalization() const {
#ifdef SDBM_NO_NORMALIZE
    return false;
#endif
    return unstable.size() > 0;
  }

  // dbm is already normalized
  linear_constraint_system_t
  to_linear_constraint_system(const DBM_t &dbm) const {
    linear_constraint_system_t csts;

    if (dbm.is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    // Extract all the edges
    SubGraph<graph_t> g_excl(const_cast<graph_t &>(dbm.g), 0);
    for (vert_id v : g_excl.verts()) {
      if (!dbm.rev_map[v])
        continue;
      if (dbm.g.elem(v, 0)) {
        variable_t vv = *dbm.rev_map[v];
        Wt c = dbm.g.edge_val(v, 0);
        csts += linear_constraint_t(linear_expression_t(vv) >= -number_t(c));
      }
      if (dbm.g.elem(0, v)) {
        variable_t vv = *dbm.rev_map[v];
        Wt c = dbm.g.edge_val(0, v);
        csts += linear_constraint_t(linear_expression_t(vv) <= number_t(c));
      }
    }

    for (vert_id s : g_excl.verts()) {
      if (!dbm.rev_map[s])
        continue;
      variable_t vs = *dbm.rev_map[s];
      for (vert_id d : g_excl.succs(s)) {
        if (!dbm.rev_map[d])
          continue;
        variable_t vd = *dbm.rev_map[d];
        csts += linear_constraint_t(vd - vs <= number_t(g_excl.edge_val(s, d)));
      }
    }
    return csts;
  }

  // dbm is already normalized
  void write(crab_os &o, const DBM_t &dbm) const {
#if 0
    o << "edges={";
    for(vert_id v : dbm.g.verts()) {
      for(vert_id d : dbm.g.succs(v)) {
	if(!dbm.rev_map[v] || !dbm.rev_map[d]) {
	  CRAB_WARN("Edge incident to un-mapped vertex.");
	  continue;
	}
	variable_t vv = *dbm.rev_map[v];
	variable_t vd = *dbm.rev_map[d];
	o << "(" << vv << "," << vd << ":"
	  << dbm.g.edge_val(v,d) << ")";
      }
    }
    o << "}";
    crab::outs() << "rev_map={";
    for(unsigned i=0, e = dbm.rev_map.size(); i!=e; i++) {
      if (dbm.rev_map[i]) {
	variable_t vi = *dbm.rev_map[i];
	crab::outs() << vi << "(" << i << ");";
      }
    }
    crab::outs() << "}\n";
#endif
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
      SubGraph<graph_t> g_excl(const_cast<graph_t &>(dbm.g), 0);
      for (vert_id v : g_excl.verts()) {
        if (!dbm.rev_map[v])
          continue;
        if (!dbm.g.elem(0, v) && !dbm.g.elem(v, 0))
          continue;
        interval_t v_out =
            interval_t(dbm.g.elem(v, 0) ? -number_t(dbm.g.edge_val(v, 0))
                                        : bound_t::minus_infinity(),
                       dbm.g.elem(0, v) ? number_t(dbm.g.edge_val(0, v))
                                        : bound_t::plus_infinity());
        if (first)
          first = false;
        else
          o << ", ";
        variable_t vv = *dbm.rev_map[v];
        o << vv << " -> " << v_out;
      }

      for (vert_id s : g_excl.verts()) {
        if (!dbm.rev_map[s])
          continue;
        variable_t vs = *dbm.rev_map[s];
        for (vert_id d : g_excl.succs(s)) {
          if (!dbm.rev_map[d])
            continue;
          variable_t vd = *dbm.rev_map[d];

          if (first)
            first = false;
          else
            o << ", ";
          o << vd << "-" << vs << "<=" << g_excl.edge_val(s, d);
        }
      }
      o << "}";
    }
  }

  split_dbm_domain(vert_map_t &&_vert_map, rev_map_t &&_rev_map, graph_t &&_g,
                   std::vector<Wt> &&_potential, vert_set_t &&_unstable)
      : vert_map(std::move(_vert_map)), rev_map(std::move(_rev_map)),
        g(std::move(_g)), potential(std::move(_potential)),
        unstable(std::move(_unstable)), _is_bottom(false) {

    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");

    if (is_top()) {
      // Garbage collection from unconstrained variables in vert_map
      // and rev_map.
      set_to_top();
    }

    CRAB_LOG("zones-split-size",
	     auto p = size();
	     print_dbm_size(p.first, p.second));
  }

  void print_dbm_size(unsigned nodes, unsigned edges) {
    if (nodes > 1) {
      crab::outs() << "#nodes=" << nodes << " "
		   << "#edges=" << edges << " "
		   << "#max-edges=" << nodes*nodes << " "
		   << "(" << ((float)edges / (float)(nodes*nodes))*100 << ")"
		   << "\n";
    }
  }

public:
  split_dbm_domain(bool is_bottom = false) : _is_bottom(is_bottom) {
    g.growTo(1); // Allocate the zero vector
    potential.push_back(Wt(0));
    rev_map.push_back(boost::none);
  }

  // FIXME: Rewrite to avoid copying if o is _|_
  split_dbm_domain(const DBM_t &o)
      : vert_map(o.vert_map), rev_map(o.rev_map), g(o.g),
        potential(o.potential), unstable(o.unstable), _is_bottom(false) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");

    CRAB_LOG("zones-split-size",
	     auto p = size();
	     print_dbm_size(p.first, p.second));
    
    if (o._is_bottom)
      set_to_bottom();

    if (!_is_bottom)
      assert(g.size() > 0);
  }

  split_dbm_domain(DBM_t &&o)
      : vert_map(std::move(o.vert_map)), rev_map(std::move(o.rev_map)),
        g(std::move(o.g)), potential(std::move(o.potential)),
        unstable(std::move(o.unstable)), _is_bottom(o._is_bottom) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  split_dbm_domain &operator=(const split_dbm_domain &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");

    if (this != &o) {
      if (o._is_bottom) {
        set_to_bottom();
      } else {
        _is_bottom = false;
        vert_map = o.vert_map;
        rev_map = o.rev_map;
        g = o.g;
        potential = o.potential;
        unstable = o.unstable;
        assert(g.size() > 0);
      }
    }

    CRAB_LOG("zones-split-size",
	     auto p = size();
	     print_dbm_size(p.first, p.second));
    
    return *this;
  }

  split_dbm_domain &operator=(split_dbm_domain &&o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");

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

  DBM_t make_top() const override { return DBM_t(false); }

  DBM_t make_bottom() const override { return DBM_t(true); }

  void set_to_top() override {
    split_dbm_domain abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    vert_map.clear();
    rev_map.clear();
    g.clear();
    potential.clear();
    unstable.clear();
    _is_bottom = true;
  }

  bool is_bottom() const override { return _is_bottom; }

  bool is_top() const override {
    if (_is_bottom)
      return false;
    return g.is_empty();
  }

  bool operator<=(const DBM_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

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
      auto leq_op = [](const DBM_t& left, const DBM_t& right) -> bool {
	// left is normalized but right doesn't need to.

	wt_ref_t wx, wy;
	
	if (left.vert_map.size() < right.vert_map.size()) {
	  return false;
	}
	
	// Set up a mapping from o to this.
	std::vector<unsigned int> vert_renaming(right.g.size(), -1);
	vert_renaming[0] = 0;
	for (auto p : right.vert_map) {
	  if (right.g.succs(p.second).size() == 0 &&
	      right.g.preds(p.second).size() == 0)
	    continue;
	  
	  auto it = left.vert_map.find(p.first);
	  // We can't have this <= o if we're missing some
	  // vertex.
	  if (it == left.vert_map.end())
	    return false;
	  vert_renaming[p.second] = (*it).second;
	  // vert_renaming[(*it).second] = p.second;
	}
	
	assert(left.g.size() > 0);
	// GrPerm g_perm(vert_renaming, g);
	
	for (vert_id ox : right.g.verts()) {
	  if (right.g.succs(ox).size() == 0)
	    continue;
	  
	  assert(vert_renaming[ox] != -1);
	  vert_id x = vert_renaming[ox];
	  for (auto edge : right.g.e_succs(ox)) {
	    vert_id oy = edge.vert;
	    assert(vert_renaming[oy] != -1);
	    vert_id y = vert_renaming[oy];
	    Wt ow = edge.val;
	    if (left.g.lookup(x, y, wx) && (wx.get() <= ow))
            continue;
	    if (!left.g.lookup(x, 0, wx) || !left.g.lookup(0, y, wy))
	      return false;
	    if (!(wx.get() + wy.get() <= ow))
	      return false;
        }
	}
	return true;
      };
      
      if (need_normalization()) {
	DBM_t left(*this);
	left.normalize();
	return leq_op(left, o);
      } else {
	return leq_op(*this, o);
      }
    }
  }

  void operator|=(const DBM_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    CRAB_LOG("zones-split", crab::outs() << "Before join:\n"
                                         << "DBM 1\n"
                                         << *this << "\n"
                                         << "DBM 2\n"
                                         << o << "\n");

    if (is_bottom()) {
      *this = o;
    } else if (o.is_top()) {
      set_to_top();
    } else if (is_top() || o.is_bottom()) {
      // do nothing
    } else {

      auto join_op = [this](DBM_t &left, const DBM_t& right)  {
	// Both left and right are normalized
	
	check_potential(left.g, left.potential, __LINE__);
	check_potential(right.g, right.potential, __LINE__);
	
	// Figure out the common renaming, initializing the
	// resulting potentials as we go.
	std::vector<vert_id> perm_x, perm_y;
	std::vector<variable_t> perm_inv;
	
	std::vector<Wt> pot_rx, pot_ry;
	vert_map_t out_vmap;
	rev_map_t out_revmap;
	// Add the zero vertex
	assert(left.potential.size() > 0);
	pot_rx.push_back(0);
	pot_ry.push_back(0);
	perm_x.push_back(0);
	perm_y.push_back(0);
	out_revmap.push_back(boost::none);
	
	for (auto p : left.vert_map) {
	  auto it = right.vert_map.find(p.first);
	  // Variable exists in both
	  if (it != right.vert_map.end()) {
	    out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
	    out_revmap.push_back(p.first);
	    pot_rx.push_back(left.potential[p.second] - left.potential[0]);
	    // XXX JNL: check this out
	    // pot_ry.push_back(right.potential[p.second] - right.potential[0]);
	    pot_ry.push_back(right.potential[(*it).second] - right.potential[0]);
	    perm_inv.push_back(p.first);
	    perm_x.push_back(p.second);
	    perm_y.push_back((*it).second);
	  }
	}
	unsigned int sz = perm_x.size();
	
	// Build the permuted view of x and y.
	assert(left.g.size() > 0);
	GrPerm gx(perm_x, left.g);
	assert(right.g.size() > 0);
	GrPerm gy(perm_y, right.g);
	
	graph_t join_g = join(gx, gy, sz, pot_rx, pot_ry);
	// Conjecture: join_g remains closed.
	
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
	
	left.vert_map = std::move(out_vmap);
	left.rev_map = std::move(out_revmap);
	left.g = std::move(join_g);
	left.potential = std::move(pot_rx);
	left.unstable.clear();
	left._is_bottom = false;
	CRAB_LOG("zones-split", crab::outs() << "Result join:\n" << left << "\n");
    };

      DBM_t &left = *this;
      left.normalize();
      if (o.need_normalization()) {
	DBM_t right(o);
	right.normalize();
	join_op(left, right);
      } else {
	join_op(left, o);
      } 
    }
  }

  DBM_t operator|(const DBM_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom()) {
      return o;
    } else if (o.is_top() || is_top()) {
      DBM_t res;
      return res;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      CRAB_LOG("zones-split", crab::outs() << "Before join:\n"
                                           << "DBM 1\n"
                                           << *this << "\n"
                                           << "DBM 2\n"
                                           << o << "\n");

      auto join_op = [this](const DBM_t &left, const DBM_t& right) -> DBM_t {
	// Both left and right are normalized
	
	check_potential(left.g, left.potential, __LINE__);
	check_potential(right.g, right.potential, __LINE__);
	
	// Figure out the common renaming, initializing the
	// resulting potentials as we go.
	std::vector<vert_id> perm_x, perm_y;
	std::vector<variable_t> perm_inv;
	
	std::vector<Wt> pot_rx, pot_ry;
	vert_map_t out_vmap;
	rev_map_t out_revmap;
	// Add the zero vertex
	assert(left.potential.size() > 0);
	pot_rx.push_back(0);
	pot_ry.push_back(0);
	perm_x.push_back(0);
	perm_y.push_back(0);
	out_revmap.push_back(boost::none);
	
	for (auto p : left.vert_map) {
	  auto it = right.vert_map.find(p.first);
	  // Variable exists in both
	  if (it != right.vert_map.end()) {
	    out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
	    out_revmap.push_back(p.first);
	    
	    pot_rx.push_back(left.potential[p.second] - left.potential[0]);
	    // XXX JNL: check this out
	    // pot_ry.push_back(right.potential[p.second] - right.potential[0]);
	    pot_ry.push_back(right.potential[(*it).second] - right.potential[0]);
	    perm_inv.push_back(p.first);
	    perm_x.push_back(p.second);
	    perm_y.push_back((*it).second);
        }
	}
	unsigned int sz = perm_x.size();
	
	// Build the permuted view of x and y.
	assert(left.g.size() > 0);
	GrPerm gx(perm_x, left.g);
	assert(right.g.size() > 0);
	GrPerm gy(perm_y, right.g);
	
	graph_t join_g = join(gx, gy, sz, pot_rx, pot_ry);
	// Conjecture: join_g remains closed.
	
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
	
	// DBM_t res(join_range, out_vmap, out_revmap, join_g, join_pot);
	DBM_t res(std::move(out_vmap), std::move(out_revmap), std::move(join_g),
		  std::move(pot_rx), vert_set_t());
	// join_g.check_adjs();
	CRAB_LOG("zones-split", crab::outs() << "Result join:\n" << res << "\n");
	return res;
      };

      if (need_normalization() && o.need_normalization()) {
	DBM_t left(*this);
	DBM_t right(o);
	left.normalize();
	right.normalize();
	return join_op(left, right);
      } else if (need_normalization()) {
	DBM_t left(*this);
	const DBM_t &right = o;
	left.normalize();
	return join_op(left, right);
      } else if (o.need_normalization()) {
	const DBM_t &left = *this;
	DBM_t right(o);
	right.normalize();
	return join_op(left, right);
      } else {
	return join_op(*this, o);
      }
    }
  }

  DBM_t operator||(const DBM_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else {
      CRAB_LOG("zones-split",
               DBM_t left(*this); // to avoid closure on left operand
               crab::outs() << "Before widening:\n"
                            << "DBM 1\n"
                            << left << "\n"
                            << "DBM 2\n"
                            << o << "\n");

      auto widen_op = [this](const DBM_t &left, const DBM_t &right) -> DBM_t {
	// Only right is normalized
	
	// Figure out the common renaming
	std::vector<vert_id> perm_x, perm_y;
	vert_map_t out_vmap;
	rev_map_t out_revmap;
	std::vector<Wt> widen_pot;
	vert_set_t widen_unstable(left.unstable);
	
	assert(left.potential.size() > 0);
	widen_pot.push_back(Wt(0));
	perm_x.push_back(0);
	perm_y.push_back(0);
	out_revmap.push_back(boost::none);
	for (auto p : left.vert_map) {
	  auto it = right.vert_map.find(p.first);
	  // Variable exists in both
	  if (it != right.vert_map.end()) {
	    out_vmap.insert(vmap_elt_t(p.first, perm_x.size()));
	    out_revmap.push_back(p.first);
	    
	    widen_pot.push_back(left.potential[p.second] - left.potential[0]);
	    perm_x.push_back(p.second);
	    perm_y.push_back((*it).second);
	  }
	}
	
	// Build the permuted view of x and y.
	assert(left.g.size() > 0);
	GrPerm gx(perm_x, left.g);
	assert(right.g.size() > 0);
	GrPerm gy(perm_y, right.g);
	
	// Now perform the widening
	std::vector<vert_id> destabilized;
	graph_t widen_g(split_widen(gx, gy, destabilized));
	for (vert_id v : destabilized)
	  widen_unstable.insert(v);
	
	DBM_t res(std::move(out_vmap), std::move(out_revmap), std::move(widen_g),
		  std::move(widen_pot), std::move(widen_unstable));
	
	CRAB_LOG("zones-split",crab::outs() << "Result widening:\n" << res << "\n");
	return res;
      };

      // Do not normalize left operand
      const DBM_t &left = *this;
      if (o.need_normalization()) {
	DBM_t right(o);
	right.normalize();
	return widen_op(left, right);
      } else {
	return widen_op(left, o);
      }
    }
  }

  DBM_t widening_thresholds(const DBM_t &o,
			    const thresholds<number_t> &ts) const override {
    // TODO: use thresholds
    return (*this || o);
  }

  void operator&=(const DBM_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      // do nothing
    } else if (is_top() || o.is_bottom()) {
      *this = o;
    } else {
      CRAB_LOG("zones-split", crab::outs() << "Before meet:\n"
                                           << "DBM 1\n"
                                           << *this << "\n"
                                           << "DBM 2\n"
                                           << o << "\n");

      auto meet_op = [this](DBM_t &left, const DBM_t & right) {
	// Both left and right are normalized
	
	check_potential(left.g, left.potential, __LINE__);
	check_potential(right.g, right.potential, __LINE__);

	// We map vertices in the left operand onto a contiguous range.
	// This will often be the identity map, but there might be gaps.
	vert_map_t meet_verts;
	rev_map_t meet_rev;
	
	std::vector<vert_id> perm_x, perm_y;
	std::vector<Wt> meet_pi;
	perm_x.push_back(0);
	perm_y.push_back(0);
	meet_pi.push_back(Wt(0));
	meet_rev.push_back(boost::none);
	for (auto p : left.vert_map) {
	  vert_id vv = perm_x.size();
	  meet_verts.insert(vmap_elt_t(p.first, vv));
	  meet_rev.push_back(p.first);
	  
	  perm_x.push_back(p.second);
	  perm_y.push_back(-1);
	  meet_pi.push_back(left.potential[p.second] - left.potential[0]);
	}
	
	// Add missing mappings from the right operand.
	for (auto p : right.vert_map) {
	  auto it = meet_verts.find(p.first);
	  
	  if (it == meet_verts.end()) {
	    vert_id vv = perm_y.size();
	    meet_rev.push_back(p.first);
	    
	    perm_y.push_back(p.second);
	    perm_x.push_back(-1);
	    meet_pi.push_back(right.potential[p.second] - right.potential[0]);
	    meet_verts.insert(vmap_elt_t(p.first, vv));
	  } else {
	    perm_y[(*it).second] = p.second;
	  }
	}
	
	// Build the permuted view of x and y.
	assert(left.g.size() > 0);
	GrPerm gx(perm_x, left.g);
	assert(right.g.size() > 0);
	GrPerm gy(perm_y, right.g);
	
	// Compute the syntactic meet of the permuted graphs.
	bool is_closed;
	graph_t meet_g(GrOps::meet(gx, gy, is_closed));
	
	// Compute updated potentials on the zero-enriched graph
	// vector<Wt> meet_pi(meet_g.size());
	// We've warm-started pi with the operand potentials
	if (!GrOps::select_potentials(meet_g, meet_pi)) {
	  // Potentials cannot be selected -- state is infeasible.
	  set_to_bottom();
	  return;
	}
	
	if (!is_closed) {
	  edge_vector delta;
	  SubGraph<graph_t> meet_g_excl(meet_g, 0);
	  // GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
	  
	  if (crab_domain_params_man::get().zones_chrome_dijkstra())
	    GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
	  else
	    GrOps::close_johnson(meet_g_excl, meet_pi, delta);
	  
	  GrOps::apply_delta(meet_g, delta);
	  
	  // Recover updated LBs and UBs.
	  if (crab_domain_params_man::get().zones_close_bounds_inline()) {
	    Wt_min min_op;
	    for (auto e : delta) {
	      if (meet_g.elem(0, e.first.first))
		meet_g.update_edge(0,
				   meet_g.edge_val(0, e.first.first) + e.second,
				   e.first.second, min_op);
	      if (meet_g.elem(e.first.second, 0))
		meet_g.update_edge(e.first.first,
				   meet_g.edge_val(e.first.second, 0) + e.second,
				   0, min_op);
	    }
	  } else {
	    delta.clear();
	    GrOps::close_after_assign(meet_g, meet_pi, 0, delta);
	    GrOps::apply_delta(meet_g, delta);
	  }
	}
      
	check_potential(meet_g, meet_pi, __LINE__);

      
	left.vert_map = std::move(meet_verts);
	left.rev_map = std::move(meet_rev);
	left.g = std::move(meet_g);
	left.potential = std::move(meet_pi);
	left.unstable.clear();
	left._is_bottom = false;

	CRAB_LOG("zones-split", crab::outs() << "Result meet:\n" << left << "\n");
      };

      DBM_t &left = *this;
      left.normalize();

      if (o.need_normalization()) {
	DBM_t right(o);
	right.normalize();
	meet_op(left, right);
      } else {
	meet_op(left, o);
      }
    }
  }

  
  DBM_t operator&(const DBM_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top())
      return *this;
    else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      CRAB_LOG("zones-split", crab::outs() << "Before meet:\n"
                                           << "DBM 1\n"
                                           << *this << "\n"
                                           << "DBM 2\n"
                                           << o << "\n");

      auto meet_op = [this](const DBM_t &left, const DBM_t &right) -> DBM_t {
	// Both left and right are normalized 
	check_potential(left.g, left.potential, __LINE__);
	check_potential(right.g, right.potential, __LINE__);
	
	// We map vertices in the left operand onto a contiguous range.
	// This will often be the identity map, but there might be gaps.
	vert_map_t meet_verts;
	rev_map_t meet_rev;
	
	std::vector<vert_id> perm_x, perm_y;
	std::vector<Wt> meet_pi;
	perm_x.push_back(0);
	perm_y.push_back(0);
	meet_pi.push_back(Wt(0));
	meet_rev.push_back(boost::none);
	for (auto p : left.vert_map) {
	  vert_id vv = perm_x.size();
	  meet_verts.insert(vmap_elt_t(p.first, vv));
	  meet_rev.push_back(p.first);
	  
	  perm_x.push_back(p.second);
	  perm_y.push_back(-1);
	  meet_pi.push_back(left.potential[p.second] - left.potential[0]);
	}
	
	// Add missing mappings from the right operand.
	for (auto p : right.vert_map) {
	  auto it = meet_verts.find(p.first);
	  
	  if (it == meet_verts.end()) {
	    vert_id vv = perm_y.size();
	    meet_rev.push_back(p.first);
	    
	    perm_y.push_back(p.second);
	    perm_x.push_back(-1);
	    meet_pi.push_back(right.potential[p.second] - right.potential[0]);
	    meet_verts.insert(vmap_elt_t(p.first, vv));
	  } else {
	    perm_y[(*it).second] = p.second;
	  }
	}
	
	// Build the permuted view of x and y.
	assert(left.g.size() > 0);
	GrPerm gx(perm_x, left.g);
	assert(right.g.size() > 0);
	GrPerm gy(perm_y, right.g);
	
	// Compute the syntactic meet of the permuted graphs.
	bool is_closed;
	graph_t meet_g(GrOps::meet(gx, gy, is_closed));
	
	// Compute updated potentials on the zero-enriched graph
	// vector<Wt> meet_pi(meet_g.size());
	// We've warm-started pi with the operand potentials
	if (!GrOps::select_potentials(meet_g, meet_pi)) {
	  // Potentials cannot be selected -- state is infeasible.
	  DBM_t res;
	  res.set_to_bottom();
	  return res;
	}
	
	if (!is_closed) {
	  edge_vector delta;
	  SubGraph<graph_t> meet_g_excl(meet_g, 0);
	  // GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
	  
	  if (crab_domain_params_man::get().zones_chrome_dijkstra())
	    GrOps::close_after_meet(meet_g_excl, meet_pi, gx, gy, delta);
	  else
	    GrOps::close_johnson(meet_g_excl, meet_pi, delta);
	  
	  GrOps::apply_delta(meet_g, delta);
	  
	  // Recover updated LBs and UBs.
	  if (crab_domain_params_man::get().zones_close_bounds_inline()) {
	    Wt_min min_op;
	    for (auto e : delta) {
	      if (meet_g.elem(0, e.first.first))
		meet_g.update_edge(0,
				   meet_g.edge_val(0, e.first.first) + e.second,
                                 e.first.second, min_op);
	      if (meet_g.elem(e.first.second, 0))
		meet_g.update_edge(e.first.first,
				   meet_g.edge_val(e.first.second, 0) + e.second,
				   0, min_op);
	    }
	  } else {
	    delta.clear();
	    GrOps::close_after_assign(meet_g, meet_pi, 0, delta);
	    GrOps::apply_delta(meet_g, delta);
	  }
	}
	check_potential(meet_g, meet_pi, __LINE__);
	DBM_t res(std::move(meet_verts), std::move(meet_rev), std::move(meet_g),
		  std::move(meet_pi), vert_set_t());
	CRAB_LOG("zones-split", crab::outs() << "Result meet:\n" << res << "\n");
	return res;
      };

      if (need_normalization() && o.need_normalization()) {
	DBM_t left(*this);
	DBM_t right(o);
	left.normalize();
	right.normalize();
	return meet_op(left, right);
      } else if (need_normalization()) {
	DBM_t left(*this);
	const DBM_t &right = o;
	left.normalize();
	return meet_op(left, right);
      } else if (o.need_normalization()) {
	const DBM_t &left = *this;
	DBM_t right(o);
	right.normalize();
	return meet_op(left, right);
      } else {
	return meet_op(*this, o);
      }
    }
  }

  DBM_t operator&&(const DBM_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_top())
      return *this;
    else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      // TODO: Implement properly
#if 1
      // Narrowing implemented as meet might not terminate.
      // Make sure that there is always a maximum bound for narrowing iterations.
      return *this & o;
#else
      // Narrowing as a no-op: sound and it will terminate 
      CRAB_LOG("zones-split", crab::outs() << "Before narrowing:\n"
	       << "DBM 1\n" << *this << "\n" << "DBM 2\n" << o << "\n");
      
      if (need_normalization()) {
	DBM_t res(*this);
	res.normalize();      
	CRAB_LOG("zones-split", crab::outs() << "Result narrowing:\n" << res << "\n");
	return res;
      } else {
	CRAB_LOG("zones-split", crab::outs() << "Result narrowing:\n" << *this << "\n");
	return *this;
      }
#endif 	       
    }
  }

  void normalize() override {
    crab::CrabStats::count(domain_name() + ".count.normalize");
    crab::ScopedCrabStats __st__(domain_name() + ".normalize");

    // Always maintained in normal form, except for widening
    if (!need_normalization()) {
      return;
    }

    edge_vector delta;
    // GrOps::close_after_widen(g, potential, vert_set_wrap_t(unstable), delta);
    // GKG: Check
    SubGraph<graph_t> g_excl(g, 0);
    if (crab_domain_params_man::get().zones_widen_restabilize())
      GrOps::close_after_widen(g_excl, potential, vert_set_wrap_t(unstable),
                               delta);
    else
      GrOps::close_johnson(g_excl, potential, delta);
    // Retrive variable bounds
    GrOps::close_after_assign(g, potential, 0, delta);

    GrOps::apply_delta(g, delta);

    unstable.clear();
  }

  void minimize() override {}

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom())
      return;
    normalize();

    auto it = vert_map.find(v);
    if (it != vert_map.end()) {
      CRAB_LOG("zones-split", crab::outs() << "Before forget " << it->second
                                           << ": " << g << "\n");
      g.forget(it->second);
      CRAB_LOG("zones-split", crab::outs() << "After: " << g << "\n");
      rev_map[it->second] = boost::none;
      vert_map.erase(v);
    }
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (is_bottom()) {
      return;
    }

    CRAB_LOG("zones-split", crab::outs() << "Before assign: " << *this << "\n");
    CRAB_LOG("zones-split", crab::outs() << x << ":=" << e << "\n");
    normalize();

    check_potential(g, potential, __LINE__);

    interval_t x_int = eval_interval(e);

    boost::optional<Wt> lb_w, ub_w;
    bool overflow;
    if (x_int.lb().is_finite()) {
      lb_w = ntow::convert(-(*(x_int.lb().number())), overflow);
      if (overflow) {
        operator-=(x);
        CRAB_LOG("zones-split", crab::outs() << "---" << x << ":=" << e << "\n"
                                             << *this << "\n");
        return;
      }
    }
    if (x_int.ub().is_finite()) {
      ub_w = ntow::convert(*(x_int.ub().number()), overflow);
      if (overflow) {
        operator-=(x);
        CRAB_LOG("zones-split", crab::outs() << "---" << x << ":=" << e << "\n"
                                             << *this << "\n");
        return;
      }
    }

    bool is_rhs_constant = false;
    // If it's a constant, just assign the interval.
    if (!crab_domain_params_man::get().zones_close_bounds_inline()) {
      // JN: it seems that we can only do this if
      // close_bounds_inline is disabled. Otherwise, the meet
      // operator misses some non-redundant edges. Need to
      // investigate more this.
      if (boost::optional<number_t> x_n = x_int.singleton()) {
        set(x, *x_n);
        is_rhs_constant = true;
      }
    } else {
      if (e.is_constant()) {
        set(x, e.constant());
        is_rhs_constant = true;
      }
    }

    if (!is_rhs_constant) {
      std::vector<std::pair<variable_t, Wt>> diffs_lb, diffs_ub;
      // Construct difference constraints from the assignment
      diffcsts_of_assign(x, e, diffs_lb, diffs_ub);
      if (diffs_lb.size() > 0 || diffs_ub.size() > 0) {
        if (crab_domain_params_man::get().zones_special_assign()) {
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

          // apply_delta should be safe here, as x has no edges in G.
          GrOps::apply_delta(g, delta);
          delta.clear();
          SubGraph<graph_t> g_excl(g, 0);
          GrOps::close_after_assign(g_excl, potential, v, delta);
          GrOps::apply_delta(g, delta);

          Wt_min min_op;
          if (lb_w) {
            g.update_edge(v, *lb_w, 0, min_op);
          }
          if (ub_w) {
            g.update_edge(0, *ub_w, v, min_op);
          }
          // Clear the old x vertex
          operator-=(x);
          vert_map.insert(vmap_elt_t(x, v));
        } else {
          // Assignment as a sequence of edge additions.
          vert_id v = g.new_vertex();
          assert(v <= rev_map.size());
          if (v == rev_map.size()) {
            rev_map.push_back(x);
            potential.push_back(Wt(0));
          } else {
            potential[v] = Wt(0);
            rev_map[v] = x;
          }
          Wt_min min_op;
          edge_vector cst_edges;

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
            check_potential(g, potential, __LINE__);
            close_over_edge(src, dest);
            check_potential(g, potential, __LINE__);
          }

          if (lb_w) {
            g.update_edge(v, *lb_w, 0, min_op);
          }
          if (ub_w) {
            g.update_edge(0, *ub_w, v, min_op);
          }

          // Clear the old x vertex
          operator-=(x);
          vert_map.insert(vmap_elt_t(x, v));
        }
      } else {
        set(x, x_int);
      }
    }

    // CRAB_WARN("DBM only supports a cst or var on the rhs of assignment");
    // this->operator-=(x);
    // g.check_adjs();

    check_potential(g, potential, __LINE__);
    CRAB_LOG("zones-split", crab::outs() << "---" << x << ":=" << e << "\n"
                                         << *this << "\n");
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom()) {
      return;
    }

    normalize();

    switch (op) {
    case OP_ADDITION:
      assign(x, y + z);
      return;
    case OP_SUBTRACTION:
      assign(x, y - z);
      return;
    // For the rest of operations, we fall back on intervals.
    case OP_MULTIPLICATION:
      set(x, get_interval(y) * get_interval(z));
      break;
    case OP_SDIV:
      set(x, get_interval(y) / get_interval(z));
      break;
    case OP_UDIV:
      set(x, get_interval(y).UDiv(get_interval(z)));
      break;
    case OP_SREM:
      set(x, get_interval(y).SRem(get_interval(z)));
      break;
    default:
      //case OP_UREM:
      set(x, get_interval(y).URem(get_interval(z)));
      break;
    }

    CRAB_LOG("zones-split", crab::outs()
                                << "---" << x << ":=" << y << op << z << "\n"
                                << *this << "\n");
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom()) {
      return;
    }

    normalize();

    switch (op) {
    case OP_ADDITION:
      assign(x, y + k);
      return;
    case OP_SUBTRACTION:
      assign(x, y - k);
      return;
    case OP_MULTIPLICATION:
      assign(x, k * y);
      return;
    // For the rest of operations, we fall back on intervals.
    case OP_SDIV:
      set(x, get_interval(y) / interval_t(k));
      break;
    case OP_UDIV:
      set(x, get_interval(y).UDiv(interval_t(k)));
      break;
    case OP_SREM:
      set(x, get_interval(y).SRem(interval_t(k)));
      break;
    default:
      //case OP_UREM:
      set(x, get_interval(y).URem(interval_t(k)));
      break;
    }

    CRAB_LOG("zones-split", crab::outs()
                                << "---" << x << ":=" << y << op << k << "\n"
                                << *this << "\n");
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const DBM_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    crab::domains::BackwardAssignOps<DBM_t>::assign(*this, x, e, inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const DBM_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<DBM_t>::apply(*this, op, x, y, z, inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const DBM_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<DBM_t>::apply(*this, op, x, y, z, inv);
  }

  void operator+=(const linear_constraint_t &cst) {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (cst.is_tautology()) {
      return;
    }

    if (cst.is_contradiction()) {
      set_to_bottom();
      return;
    }
    
    if (is_bottom()) {
      return;
    }
    
    normalize();
    // g.check_adjs();
    if (cst.is_inequality()) {
      if (!add_linear_leq(cst.expression())) {
        set_to_bottom();
      }
    } else if (cst.is_strict_inequality()) {
      // We try to convert a strict to non-strict.
      auto nc = ikos::linear_constraint_impl::strict_to_non_strict_inequality(cst);
      if (nc.is_inequality()) {
        // here we succeed
        if (!add_linear_leq(nc.expression())) {
          set_to_bottom();
        }
      }
    } else if (cst.is_equality()) {
      const linear_expression_t &exp = cst.expression();
      if (!add_linear_leq(exp) || !add_linear_leq(-exp)) {
        set_to_bottom();
      }
      // g.check_adjs();
    } else if (cst.is_disequation()) {
      // We handle here the case x !=y by converting the disequation
      // into a strict inequality if possible.
      linear_constraint_system_t csts;
      constraint_simp_domain_traits<DBM_t>::lower_disequality(*this, cst, csts);
      for (auto const& c: csts) {
	// We try to convert a strict inequality into non-strict one
	auto nc = ikos::linear_constraint_impl::strict_to_non_strict_inequality(c);
	if (nc.is_inequality()) {
	  // here we succeed
	  if (!add_linear_leq(nc.expression())) {
	    set_to_bottom();
	  }
	}
      }
      
      if (!is_bottom()) {
	// We handle here the case x != c
	add_disequation(cst.expression());
      }
    }

    CRAB_LOG("zones-split", crab::outs() << "---" << cst << "\n"
	     << *this << "\n");
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    if (is_bottom())
      return;

    for (auto cst : csts) {
      operator+=(cst);
    }
  }

  virtual bool entails(const linear_constraint_t &cst) const override {	
    if (is_bottom()) {							
      return true;							
    } if (cst.is_tautology()) {						
      return true;							
    } if (cst.is_contradiction()) {					
      return false;							
    }

    bool res;
    if (cst.is_disequation()) {
      // |= c1.x1 + ... + cn.xn != k is iff 
      // (1) |= c1.x1 + ... + cn.xn < k OR
      // (2) |= c1.x1 + ... + cn.xn > k
      linear_constraint_t pob1(cst.expression(), linear_constraint_t::kind_t::STRICT_INEQUALITY);
      res = is_unsat(pob1.negate());
      if (!res) {
	linear_constraint_t pob2(cst.expression() * number_t(-1), linear_constraint_t::kind_t::STRICT_INEQUALITY);
	res = is_unsat(pob2.negate());
      }
    } else if (cst.is_equality()) {
      // |= c1.x1 + ... + cn.xn == k is iff 
      // (1) |= c1.x1 + ... + cn.xn <= k AND
      // (2) |= c1.x1 + ... + cn.xn >= k
      linear_constraint_t pob1(cst.expression(), linear_constraint_t::kind_t::INEQUALITY);
      res = is_unsat(pob1.negate());
      if (res) {
	linear_constraint_t pob2(cst.expression() * number_t(-1), linear_constraint_t::kind_t::INEQUALITY);
	res = is_unsat(pob2.negate());
      }
    } else {
      // cst is an inequality
      res = is_unsat(cst.negate());
    }
    
    return res;
  }
  
  interval_t operator[](const variable_t &x) override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");
    normalize();
    return (is_bottom() ? interval_t::bottom() :
	    get_interval(vert_map, g, x));    
  }

  interval_t at(const variable_t &x) const override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");
    return (is_bottom() ? interval_t::bottom() :
	    get_interval(vert_map, g, x));
  }

  void set(const variable_t &x, interval_t intv) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (is_bottom())
      return;

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
    }
    if (intv.lb().is_finite()) {
      Wt lb = ntow::convert(*(intv.lb().number()), overflow);
      if (overflow) {
        return;
      }
      potential[v] = potential[0] + lb;
      g.set_edge(v, -lb, 0);
    }
  }

  // int_cast_operators_api
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    int_cast_domain_traits<DBM_t>::apply(*this, op, dst, src);    
  }

  // bitwise_operators_api
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;
    normalize();

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi = operator[](z);
    interval_t xi = interval_t::bottom();
    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    default:
      //case OP_ASHR: 
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;
    normalize();

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi(k);
    interval_t xi = interval_t::bottom();

    switch (op) {
    case OP_AND: {
      xi = yi.And(zi);
      break;
    }
    case OP_OR: {
      xi = yi.Or(zi);
      break;
    }
    case OP_XOR: {
      xi = yi.Xor(zi);
      break;
    }
    case OP_SHL: {
      xi = yi.Shl(zi);
      break;
    }
    case OP_LSHR: {
      xi = yi.LShr(zi);
      break;
    }
    default:
      //case OP_ASHR: 
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }

  /// split_dbm_domain implements only standard abstract operations
  /// of a numerical domain so it is intended to be used as a leaf
  /// domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(DBM_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(DBM_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(DBM_t)
  DEFAULT_SELECT(DBM_t)
  DEFAULT_WEAK_ASSIGN(DBM_t)
    
  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

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
      if (it != vert_map.end())
        save[(*it).second] = true;
    }

    for (vert_id v = 0; v < rev_map.size(); v++) {
      if (!save[v] && rev_map[v]) {
        variable_t vv = *rev_map[v];
        operator-=(vv);
      }
    }
  }

  void forget(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top()) {
      return;
    }

    for (auto v : variables) {
      operator-=(v);
    }
  }

  void expand(const variable_t &x, const variable_t &y) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    CRAB_LOG("zones-split", crab::outs() << "Before expand " << x << " into "
                                         << y << ":\n"
                                         << *this << "\n");

    auto it = vert_map.find(y);
    if (it != vert_map.end()) {
      CRAB_ERROR("split_dbm expand operation failed because y already exists");
    }

    vert_id ii = get_vert(x);
    vert_id jj = get_vert(y);
    edge_vector delta;
    for (auto edge : g.e_preds(ii)) {
      delta.push_back({{edge.vert, jj}, edge.val});
    }

    for (auto edge : g.e_succs(ii)) {
      delta.push_back({{jj, edge.vert}, edge.val});
    }
    GrOps::apply_delta(g, delta);

    potential[jj] = potential[ii];

    CRAB_LOG("zones-split", crab::outs() << "After expand " << x << " into "
                                         << y << ":\n"
                                         << *this << "\n");
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_top() || is_bottom())
      return;

    CRAB_LOG("zones-split", crab::outs() << "Renaming {";
             for (auto v
                  : from) crab::outs()
             << v << ";";
             crab::outs() << "} with "; for (auto v
                                             : to) crab::outs()
                                        << v << ";";
             crab::outs() << "}:\n"; crab::outs() << *this << "\n";);

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      variable_t v = from[i];
      variable_t new_v = to[i];
      if (v == new_v) { // nothing to rename
        continue;
      }

      { // We do garbage collection of unconstrained variables only
        // after joins so it's possible to find new_v but we are ok as
        // long as it's unconstrained.
        auto it = vert_map.find(new_v);
        if (it != vert_map.end()) {
          vert_id dim = it->second;
          if (g.succs(dim).size() != 0 || g.preds(dim).size() != 0) {
            CRAB_ERROR(domain_name() + "::rename assumes that ", new_v,
                       " does not exist");
          }
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

    CRAB_LOG("zones-split", crab::outs() << "RESULT=" << *this << "\n");
  }

  void extract(const variable_t &x, linear_constraint_system_t &csts,
               bool only_equalities) {
    crab::CrabStats::count(domain_name() + ".count.extract");
    crab::ScopedCrabStats __st__(domain_name() + ".extract");

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

  // intrinsics operations
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const DBM_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  // Output function
  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    // linear_constraint_system_t inv = to_linear_constraint_system();
    // o << inv;

    if (need_normalization()) {
      DBM_t tmp(*this);
      tmp.normalize();
      write(o, tmp);
    } else {
      write(o, *this);
    }
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() + ".count.to_linear_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".to_linear_constraints");

    if (need_normalization()) {
      DBM_t tmp(*this);
      tmp.normalize();
      return to_linear_constraint_system(tmp);
    } else {
      return to_linear_constraint_system(*this);
    }
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
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

  std::string domain_name() const override { return "SplitDBM"; }
}; // class split_dbm_domain

template <typename Number, typename VariableName, typename Params>
struct abstract_domain_traits<split_dbm_domain<Number, VariableName, Params>> {
  using number_t = Number;
  using varname_t = VariableName;
};

template <typename Number, typename VariableName, typename Params>
class reduced_domain_traits<split_dbm_domain<Number, VariableName, Params>> {
public:
  using sdbm_domain_t = split_dbm_domain<Number, VariableName, Params>;
  using variable_t = typename sdbm_domain_t::variable_t;
  using linear_constraint_system_t =
      typename sdbm_domain_t::linear_constraint_system_t;

  static void extract(sdbm_domain_t &dom, const variable_t &x,
                      linear_constraint_system_t &csts, bool only_equalities) {
    dom.extract(x, csts, only_equalities);
  }
};
} // namespace domains
} // namespace crab

#pragma GCC diagnostic pop
