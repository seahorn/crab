/*******************************************************************************
 * Implementation of Octagons based on Split Normal Form (SNF).
 *
 * Based on the draft "A Fresh Look at Zones and Octagons" by Gange,
 * Ma, Navas, Schachte, Sondergaard, and Stuckey.
 *
 * A graph representing octagon constraints is in SNF if
 *
 * 1) the graph explicitly stores the strongest binary constraints
 *    implied by other binary constraints but not necessarily by bounds.
 * 2) the graph explicitly stores the strongest variable bounds.
 *
 * - Closure and coherence are preserved incrementally while adding
 *   constraints and assignments.
 *
 *   coherence means that the DBM is kind of symmetric:
 *   forall i,j in DBM :: mi,j = mj',i'
 *   where i' =  i+1 if i is even or i-1 if odd
 *
 * - The strengthening step means that an edge can be improved via
 *   bounds: forall i,j in DBM :: mi,j = min (mi,j, (mi,i'+ mj',j)/2)
 *
 *   Very importantly, unlike traditional octagon domain, we don't
 *   need to perform explicitly this step since the DBM is in SNF.
 *   Instead, we adapt all the transfer functions to consider that a
 *   binary constraint can be either explicit or implicit via two
 *   unary constraints.
 *
 * Todos:
 *
 * - Integer tightening is implemented partially: potential function
 *   is not refined (i.e., compute-potential-integer Fig 32 is not
 *   implemented).
 * - Linear division can be improved.
 * - Meet operator needs more testing.
 *******************************************************************************/

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

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

#define JOIN_CLOSE_AFTER_MEET
//#define CHECK_POTENTIAL
#define INTEGER_TIGHTENING
#define USE_FLAT_MAP

#ifdef USE_FLAT_MAP
#include <boost/container/flat_map.hpp>
#else
// Operations like rename are much faster using unordered_map
#include <unordered_map>
#endif

namespace crab {
namespace domains {
namespace split_octagons_impl {
/*
   Special view that skips bounds from a DBM containing octagon constraints.
   Given pos(x) (neg(x)) its successor/predecessor neg(x) (pos(x)) is skipped.
*/
template <class G> class SplitOctGraph {
public:
  using vert_id = typename G::vert_id;
  using Wt = typename G::Wt;
  // These defined below
  // using const_pred_range = ...
  // using const_succ_range = ...
  // using const_e_pred_range = ...
  // using const_e_succ_range = ...
  using wt_ref_t = typename G::wt_ref_t;

  SplitOctGraph(const G &_g) : g(_g) {}

  bool elem(vert_id x, vert_id y) const {
    return ((x / 2 != y / 2) && g.elem(x, y));
  }

  bool lookup(vert_id x, vert_id y, wt_ref_t &w) const {
    return ((x / 2 != y / 2) && g.lookup(x, y, w));
  }

  Wt edge_val(vert_id x, vert_id y) const { return g.edge_val(x, y); }

  // Precondition: elem(x, y) is true.
  Wt operator()(vert_id x, vert_id y) const { return g(x, y); }

  // Number of allocated vertices
  int size(void) const { return g.size(); }

  class vert_iterator {
  public:
    vert_iterator(const typename G::vert_iterator &_iG) : iG(_iG) {}

    // Skipping of v_ex is done entirely by !=.
    // So we _MUST_ test it != verts.end() before dereferencing.
    vert_id operator*(void) { return *iG; }
    vert_iterator operator++(void) {
      ++iG;
      return *this;
    }
    bool operator!=(const vert_iterator &o) { return iG != o.iG; }

    typename G::vert_iterator iG;
  };

  class vert_range {
  public:
    vert_range(const typename G::vert_range &_rG) : rG(_rG) {}

    vert_iterator begin(void) const { return vert_iterator(rG.begin()); }
    vert_iterator end(void) const { return vert_iterator(rG.end()); }

    typename G::vert_range rG;
  };

  vert_range verts(void) const { return vert_range(g.verts()); }

  /* exclude _v_ex from the view */
  template <class It> class adj_iterator {
  public:
    adj_iterator(const It &_iG, vert_id _v_ex) : iG(_iG), v_ex(_v_ex) {}
    vert_id operator*(void)const { return *iG; }
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

  /* exclude _v_ex from the view */
  template <class It> class e_adj_iterator {
  public:
    using edge_wrapper = typename It::edge_wrapper;

    e_adj_iterator(const It &_iG, vert_id _v_ex) : iG(_iG), v_ex(_v_ex) {}
    edge_wrapper operator*(void)const { return *iG; }
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

  /* exclude _v_ex from the view */
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

  // If v is v+ (v-) then v- (v+) does not appear as a successor
  const_succ_range succs(vert_id v) const {
    // assert(v != v_ex);
    vert_id v_opp;
    if (v % 2 == 0)
      v_opp = v + 1;
    else
      v_opp = v - 1;
    return const_succ_range(g.succs(v), v_opp);
  }

  // If v is v+ (v-) then v- (v+) does not appear as a predecessor
  const_pred_range preds(vert_id v) const {
    // assert(v != v_ex);
    vert_id v_opp;
    if (v % 2 == 0)
      v_opp = v + 1;
    else
      v_opp = v - 1;
    return const_pred_range(g.preds(v), v_opp);
  }

  const_e_succ_range e_succs(vert_id v) const {
    vert_id v_opp;
    if (v % 2 == 0)
      v_opp = v + 1;
    else
      v_opp = v - 1;
    return const_e_succ_range(g.e_succs(v), v_opp);
  }
  const_e_pred_range e_preds(vert_id v) const {
    vert_id v_opp;
    if (v % 2 == 0)
      v_opp = v + 1;
    else
      v_opp = v - 1;
    return const_e_pred_range(g.e_preds(v), v_opp);
  }

  const G &g;
};
} // end namespace split_octagons_impl

template <class Number, class VariableName,
          class Params = DBM_impl::DefaultParams<Number>>
class split_oct_domain final
    : public abstract_domain_api<
          split_oct_domain<Number, VariableName, Params>> {
  using split_oct_domain_t = split_oct_domain<Number, VariableName, Params>;
  using abstract_domain_t = abstract_domain_api<split_oct_domain_t>;

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
  using wt_ref_t = typename graph_t::wt_ref_t;  
  using ntow = DBM_impl::NtoW<number_t, Wt>;
  using vert_id = typename graph_t::vert_id;
  // (variable: (vert_id of pos, vert_id of neg))
#ifdef USE_FLAT_MAP
  using vert_map_t =
      boost::container::flat_map<variable_t, std::pair<vert_id, vert_id>>;
#else
  using vert_map_t =
      std::unordered_map<variable_t, std::pair<vert_id, vert_id>>;
#endif
  using vmap_elt_t = typename vert_map_t::value_type;
  using rev_map_t = std::vector<boost::optional<variable_t>>;
  using GrOps = GraphOps<graph_t>;
  using GrPerm = GraphPerm<graph_t>;
  using edge_vector = typename GrOps::edge_vector;
  using diffcst_t = std::pair<std::pair<variable_t, variable_t>, Wt>;
  using vert_set_t = std::unordered_set<vert_id>;

  // Map variables to graph nodes
  vert_map_t m_vert_map;
  rev_map_t m_rev_vert_map;
  /*
     Octagonal constraints are encoded as difference constraints in
     the graph.  For each program variable v, we have two graph nodes
     pos(v) and neg(v). Then, we encode constraints as follows:

     x - y <= k  is encoded as pos(x)-pos(y) <= k and neg(y)-neg(x) <= k
     x + y <= k  is encoded as pos(x)-neg(y) <= k and pos(y)-neg(x) <= k
     -x -y <= k  is encoded as neg(x)-pos(y) <= k and neg(y)-pos(x) <= k
         x <= k  is encoded as pos(x)-neg(x) <= 2*k
        -x <= k  is encoded as neg(x)-pos(x) <= 2*k
   */
  graph_t m_graph;
  // Potential funtions: model of octagon constraints
  std::vector<Wt> m_potential;
  // Graph nodes that need to be closed after widening
  vert_set_t m_unstable;
  bool m_is_bottom;

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

  // get vert_id of both v+ and v- after creating vertices in graph
  // return vert_id of positive node, negative node is positive + 1
  vert_id get_vert(variable_t v) {
    auto it = m_vert_map.find(v);
    if (it != m_vert_map.end())
      return (*it).second.first;

    vert_id vert_pos(m_graph.new_vertex());
    vert_id vert_neg(m_graph.new_vertex());
    if (vert_pos > vert_neg) {
      vert_id tmp = vert_pos;
      vert_pos = vert_neg;
      vert_neg = tmp;
    }
    m_vert_map.insert(vmap_elt_t(v, {vert_pos, vert_neg}));
    assert(vert_pos <= m_rev_vert_map.size());
    assert(vert_neg <= m_rev_vert_map.size() + 1);

    if (vert_pos < m_rev_vert_map.size()) {
      m_potential[vert_pos] = Wt(0);
      m_rev_vert_map[vert_pos] = v;
    } else {
      m_potential.push_back(Wt(0));
      m_rev_vert_map.push_back(v);
    }

    if (vert_neg < m_rev_vert_map.size()) {
      m_potential[vert_neg] = Wt(0);
      m_rev_vert_map[vert_neg] = v;
    } else {
      m_potential.push_back(Wt(0));
      m_rev_vert_map.push_back(v);
    }

    return vert_pos;
  }

  class vert_set_wrap_t {
  public:
    const vert_set_t &vs;
    vert_set_wrap_t(const vert_set_t &_vs) : vs(_vs) {}
    bool operator[](vert_id v) const { return vs.find(v) != vs.end(); }
  };

  template <class G, class P>
  inline bool check_potential(G &g, P &p, unsigned line) const {
#ifdef CHECK_POTENTIAL
    crab::ScopedCrabStats __st__(domain_name() + ".check_potential");
    for (vert_id v : g.verts()) {
      for (vert_id d : g.succs(v)) {
        if (p[v] + g.edge_val(v, d) - p[d] < Wt(0)) {
          g.write(crab::outs());
          crab::outs() << "\n";
          CRAB_ERROR("Invalid potential at line ", line, ":", "pot[", v,
                     "]=", p[v], " ", "pot[", d, "]=", p[d], " ", "edge(", v,
                     ",", d, ")=", g.edge_val(v, d));
          return false;
        }
      }
    }
#endif
    return true;
  }

  // get potential of a variable
  Wt pot_value(const variable_t &v) {
    auto it = m_vert_map.find(v);
    if (it != m_vert_map.end()) {
      vert_id v = (*it).second.first;
      return (m_potential[v] - m_potential[v + 1]) / (Wt)2;
    }
    return ((Wt)0);
  }

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
      v += (pot_value(p.second)) * coef;
    }
    return v;
  }

  interval_t eval_interval(linear_expression_t e) {
    interval_t r = e.constant();
    for (auto p : e)
      r += p.first * operator[](p.second);
    return r;
  }

  // Update a set of edges while repairing potential
  // JN: not sure repair_potential is needed.
  bool update_delta(graph_t &g, std::vector<Wt> &pot,
                    edge_vector &delta) const {
    Wt_min min_op;
    for (auto edge : delta) {
      vert_id src = edge.first.first;
      vert_id dst = edge.first.second;
      Wt w = edge.second;
      g.update_edge(src, w, dst, min_op);
      if (!GrOps::repair_potential(g, pot, src, dst)) {
        return false;
      }
    }
    check_potential(g, pot, __LINE__);
    return true;
  }

  // JN: apply_delta takes as input edges that we know they already
  // have more precise weights than existing ones. If that's the case
  // this function is unnecessary.
  template <class G> void update_delta(G &g, edge_vector &delta) const {
    wt_ref_t w;
    for (std::pair<std::pair<vert_id, vert_id>, Wt> &e : delta) {
      if (!g.lookup(e.first.first, e.first.second, w) || e.second < w.get()) {
        g.set_edge(e.first.first, Wt(e.second), e.first.second);
      }
    }
  }

  /**
   *  Turn an assignment into a set of octagon constraints.
   *
   *  Given v := a*x + b*y + E + k, we generate the octagon
   *  constraints:
   *
   *  if extract_upper_bounds
   *     v - x <= ub((a-1)*x + b*y + E + k) if a >= 0
   *     v + x <= ub((a+1)*x + b*y + E + k) if a < 0
   *     v - y <= ub(a*x + (b-1)*y + E + k) if b >= 0
   *     v + y <= ub(a*x + (b+1)*y + E + k) if b < 0
   *  else
   *     v - x >= lb((a-1)*x + b*y + E + k) if a >= 0
   *     v + x >= lb((a+1)*x + b*y + E + k) if a < 0
   *     v - y >= lb(a*x + (b-1)*y + E + k) if b >= 0
   *     v + y >= lb(a*x + (b+1)*y + E + k) if b < 0
   **/
  void oct_csts_of_assign(const variable_t &x, const linear_expression_t &exp,
                          bool extract_upper_bounds,
                          /* foreach {v, k} \in diff_csts we have
                             x - v <= k (if extract_upper_bounds), or
                             x - v >= k (if not extract_upper_bounds) */
                          std::vector<std::pair<variable_t, Wt>> &diff_csts,
                          /* foreach {v, k} \in oct_csts we have
                             x + v <= k (if extract_upper_bounds), or
                             x + v >= k (if not extract_upper_bounds) */
                          std::vector<std::pair<variable_t, Wt>> &oct_csts) {

    crab::ScopedCrabStats __st__(domain_name() + ".oct_csts_of_assign");
    boost::optional<variable_t> unbounded_pos_var, unbounded_neg_var;
    std::vector<std::pair<variable_t, Wt>> diff_terms, oct_terms;
    bool overflow;

    Wt residual(ntow::convert(exp.constant(), overflow));
    if (overflow) {
      return;
    }
    for (auto p : exp) {
      Wt coeff(ntow::convert(p.first, overflow));
      variable_t y(p.second);
      if (overflow) {
        continue;
      }

      if (coeff < Wt(0)) {
        bound_t y_val =
            (extract_upper_bounds ? operator[](y).lb() : operator[](y).ub());

        if (y_val.is_infinite()) {
          if (unbounded_pos_var || unbounded_neg_var || coeff != Wt(-1)) {
            return;
          }
          unbounded_neg_var = y;
        } else {
          Wt ymax(ntow::convert(*(y_val.number()), overflow));
          if (overflow) {
            continue;
          }
          residual += ymax * coeff;
          oct_terms.push_back({y, ymax});
        }
      } else {
        bound_t y_val =
            (extract_upper_bounds ? operator[](y).ub() : operator[](y).lb());

        if (y_val.is_infinite()) {
          if (unbounded_pos_var || unbounded_neg_var || coeff != Wt(1)) {
            return;
          }
          unbounded_pos_var = y;
        } else {
          Wt ymax(ntow::convert(*(y_val.number()), overflow));
          if (overflow) {
            continue;
          }
          residual += ymax * coeff;
          diff_terms.push_back({y, ymax});
        }
      }
    }

    if (unbounded_pos_var) {
      // There is exactly one unbounded positive variable with
      // unit coefficient
      diff_csts.push_back({*unbounded_pos_var, residual});
      CRAB_LOG("octagon-assign",
               if (extract_upper_bounds) {
                 crab::outs() << x << "-" << *unbounded_pos_var
                              << "<=" << residual << ";";
               } else {
                 crab::outs() << "-" << x << "+" << *unbounded_pos_var
                              << "<=" << -residual << ";";
               });
    } else if (unbounded_neg_var) {
      // There is exactly one unbounded negative variable with
      // unit coefficient.
      oct_csts.push_back({*unbounded_neg_var, residual});
      CRAB_LOG("octagon-assign",
               if (extract_upper_bounds) {
                 crab::outs() << x << "+" << *unbounded_neg_var
                              << "<=" << residual << ";";
               } else {
                 crab::outs() << "-" << x << "-" << *unbounded_neg_var
                              << "<=" << -residual << ";";
               });
    } else {
      for (auto p : diff_terms) {
        diff_csts.push_back({p.first, residual - p.second});
        CRAB_LOG("octagon-assign",
                 if (extract_upper_bounds) {
                   crab::outs() << x << "-" << p.first
                                << "<=" << (residual - p.second) << ";";
                 } else {
                   crab::outs() << "-" << x << "+" << p.first
                                << "<=" << -(residual - p.second) << ";";
                 });
      }
      for (auto p : oct_terms) {
        oct_csts.push_back({p.first, residual + p.second});
        CRAB_LOG("octagon-assign",
                 if (extract_upper_bounds) {
                   crab::outs() << x << "+" << p.first
                                << "<=" << (residual + p.second) << ";";
                 } else {
                   crab::outs() << "-" << x << "-" << p.first
                                << "<=" << -(residual + p.second) << ";";
                 });
      }
    }
  }

  // Turn an assignment into a set of octagon constraints.
  void oct_csts_of_assign(const variable_t &x, const linear_expression_t &exp,
                          std::vector<std::pair<variable_t, Wt>> &diff_lb,
                          std::vector<std::pair<variable_t, Wt>> &oct_lb,
                          std::vector<std::pair<variable_t, Wt>> &diff_ub,
                          std::vector<std::pair<variable_t, Wt>> &oct_ub) {
    oct_csts_of_assign(x, exp, true, diff_ub, oct_ub);
    oct_csts_of_assign(x, exp, false, diff_lb, oct_lb);
  }

  // Turn a linear inequality into a set of octagon constraints
  void oct_csts_of_lin_leq(const linear_expression_t &exp,
                           std::vector<diffcst_t> &diff_csts,
                           std::vector<diffcst_t> &pos_oct_csts,
                           std::vector<diffcst_t> &neg_oct_csts,
                           std::vector<std::pair<variable_t, Wt>> &lbs,
                           std::vector<std::pair<variable_t, Wt>> &ubs) {

    crab::ScopedCrabStats __st__(domain_name() + ".oct_csts_of_lin_leq");

    /*
      consider x + y - z <= k and assume lb(x),lb(y) and ub(z) are finite bounds

        pos_terms = [(1,x,lb(x)),(1,y,lb(y))]
        neg_terms = [(1,z,ub(z))]
        exp_ub = maximize(k - x - y + z) = k -lb(x) -lb(y) + ub(z)

        octagon constraints =
            x + y <= exp_ub + lb(x) + lb(y)
            x - z <= exp_ub + lb(x) - ub(z)
            y - z <= exp_ub + lb(y) - ub(z)

      If lower/upper bounds are unknown then we can still generate
      octagon constraints if the number of variables with unknown
      bounds is not greater than 2.
     */

    // {{v,coeff}, polarity} where polarity is true iff the term is positive
    // Invariant: unbounded_var2 != None => unbounded_var1 != None
    boost::optional<std::pair<std::pair<variable_t, Wt>, bool>> unbounded_var1,
        unbounded_var2;
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

    // crab::outs() << "Extracting octagon constraints from " << exp << "\n";

    std::vector<std::pair<std::pair<Wt, variable_t>, Wt>> pos_terms, neg_terms;
    for (auto p : exp) {
      Wt coeff(ntow::convert(p.first, overflow));
      if (overflow) {
        continue;
      }
      if (coeff > Wt(0)) {
        variable_t y(p.second);
        bound_t y_lb = operator[](y).lb();
        // crab::outs() << "###" << y << ">=" << y_lb << "\n";
        if (y_lb.is_infinite()) {
          // positive term with unknown lower bound
          if (!unbounded_var1) {
            unbounded_var1 = std::make_pair(std::make_pair(y, coeff), true);
          } else if (!unbounded_var2) {
            unbounded_var2 = std::make_pair(std::make_pair(y, coeff), true);
          } else {
            return;
          }
        } else {
          Wt ymin(ntow::convert(*(y_lb.number()), overflow));
          if (overflow) {
            continue;
          }
          // Coeff is negative, so it's still add
          exp_ub -= ymin * coeff;
          pos_terms.push_back({{coeff, y}, ymin});
          // crab::outs() << "###" << "Added pos term " << "(" << coeff << ","
          // << y << ")  ymin=" << ymin << "\n";
        }
      } else {
        variable_t y(p.second);
        bound_t y_ub = operator[](y).ub();
        // crab::outs() << "###" << y << "<=" << y_ub << "\n";
        if (y_ub.is_infinite()) {
          // negative term with unknown upper-bound
          if (!unbounded_var1) {
            unbounded_var1 = std::make_pair(std::make_pair(y, -coeff), false);
          } else if (!unbounded_var2) {
            unbounded_var2 = std::make_pair(std::make_pair(y, -coeff), false);
          } else {
            return;
          }
        } else {
          Wt ymax(ntow::convert(*(y_ub.number()), overflow));
          if (overflow) {
            continue;
          }
          exp_ub -= ymax * coeff;
          neg_terms.push_back({{-coeff, y}, ymax});
          // crab::outs() << "###" << "Added negative term " << "(" << -coeff <<
          // "," << y
          // 	       << ")  ymax=" << ymax << "\n";
        }
      }
    }

    if (unbounded_var1) {
      variable_t x((*unbounded_var1).first.first);
      Wt coeff_x = (*unbounded_var1).first.second;
      bool polarity_x = (*unbounded_var1).second;
      if (unbounded_var2) {
        variable_t y((*unbounded_var2).first.first);
        Wt coeff_y = (*unbounded_var2).first.second;
        bool polarity_y = (*unbounded_var2).second;
        if (coeff_x != Wt(1) || coeff_y != Wt(1)) {
          return;
        }
        //   pos and pos
        if (polarity_x && polarity_y) {
          pos_oct_csts.push_back({{x, y}, exp_ub});
        } // pos and neg
        else if (polarity_x && !polarity_y) {
          diff_csts.push_back({{x, y}, exp_ub});
        } // neg and pos
        else if (!polarity_x && polarity_y) {
          diff_csts.push_back({{y, x}, exp_ub});
        } // neg and neg
        else {
          assert(!polarity_x && !polarity_y);
          neg_oct_csts.push_back({{x, y}, exp_ub});
        }
      } else {
        if (coeff_x == Wt(1)) {
          if (polarity_x) {
            for (auto p : neg_terms) {
              diff_csts.push_back({{x, p.first.second}, exp_ub - p.second});
            }
            for (auto p : pos_terms) {
              pos_oct_csts.push_back({{x, p.first.second}, exp_ub + p.second});
            }
          } else {
            for (auto p : pos_terms) {
              diff_csts.push_back({{p.first.second, x}, exp_ub + p.second});
            }
            for (auto p : neg_terms) {
              neg_oct_csts.push_back({{x, p.first.second}, exp_ub - p.second});
            }
          }
        }
        // Bounds for x
        if (polarity_x) {
          ubs.push_back({x, exp_ub / coeff_x});
        } else {
          lbs.push_back({x, -exp_ub / coeff_x});
        }
      }
    } else {
      assert(!unbounded_var1 && !unbounded_var2);

      for (auto pl : neg_terms) {
        for (auto pu : pos_terms) {
          diff_csts.push_back({{pu.first.second, pl.first.second},
                               exp_ub - pl.second + pu.second});
        }
      }

      for (unsigned i = 0, e = pos_terms.size(); i < e; ++i) {
        for (unsigned j = i; j < e; ++j) {
          auto pu1 = pos_terms[i];
          auto pu2 = pos_terms[j];
          if (pu1.first.second == pu2.first.second)
            continue;
          pos_oct_csts.push_back({{pu1.first.second, pu2.first.second},
                                  exp_ub + pu1.second + pu2.second});
        }
      }

      for (unsigned i = 0, e = neg_terms.size(); i < e; ++i) {
        for (unsigned j = i; j < e; ++j) {
          auto pu1 = neg_terms[i];
          auto pu2 = neg_terms[j];
          if (pu1.first.second == pu2.first.second)
            continue;
          neg_oct_csts.push_back({{pu1.first.second, pu2.first.second},
                                  exp_ub - pu1.second - pu2.second});
        }
      }

      for (auto pl : neg_terms) {
        lbs.push_back({pl.first.second, -exp_ub / pl.first.first + pl.second});
      }

      for (auto pu : pos_terms) {
        ubs.push_back({pu.first.second, exp_ub / pu.first.first + pu.second});
      }
    }
  }

  bool repair_potential(vert_id src, vert_id dest) {
    return GrOps::repair_potential(m_graph, m_potential, src, dest);
  }

  /*
   * Update other bounds after a lower-bound interval constraint has
   * been added.
   *                                k
   * We just added the edge pos(v) ----> neg(v):
   *
   * (1) This is the first rule of inference:
   *     -v <= -k
   *   y +v <= k'
   *  -------------
   *      y <= k'-k
   *
   * There are two cases to implement because there are two ways of
   * representing y+v <= k':
   *
   * CASE A
   *          k            k'             k+k'
   * pos(v)------>neg(v) -----> pos(y) <- - - - neg(y)
   *
   * CASE B
   *        k'            k
   * neg(y)----->pos(v)------>neg(v)       pos(y)
   *   |                                     ^
   *   |                                     |
   *    - - - - - - - - - - - - - - - - - - -
   *
   *
   * (2) This is the second rule of inference:
   *     -v <= -k
   *   v -y <= k'
   *  ------------
   *     -y <= k'-k
   *
   * Again, there are two cases to implement because there are two
   * ways of representing v-y <= k':
   *
   * CASE C:
   *            k            k'                 k+k'
   *   pos(v) -----> neg(v) -----> neg(y) <- - - - - - pos(y)
   *
   * CASE D
   *           k+k'           k'           k
   *   neg(y)<- - - pos(y) -----> pos(v) -----> neg(v)
   *
   *
   */
  bool update_bounds_lb(vert_id v, Wt k) {
    assert(v % 2 == 0);
    // we add in delta so that we don't invalidate graph iterators
    edge_vector delta_lb, delta_ub;

    for (auto e : m_graph.e_succs(v + 1)) {
      if (e.vert % 2 == 0) {
        // CASE A
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-lb-a) add edge " << e.vert + 1 << " --> "
                     << e.vert << " with " << ((Wt)2 * e.val) + k << "\n";);
        delta_ub.push_back({{e.vert + 1, e.vert}, ((Wt)2 * e.val) + k});
      } else {
        // CASE C
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-lb-c) add edge " << e.vert - 1 << " --> "
                     << e.vert << " with " << ((Wt)2 * e.val) + k << "\n";);
        delta_lb.push_back({{e.vert - 1, e.vert}, ((Wt)2 * e.val) + k});
      }
    }

    for (auto e : m_graph.e_preds(v)) {
      if (e.vert % 2 != 0) {
        // CASE B
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-lb-b) add edge: " << e.vert << " --> "
                     << e.vert - 1 << " with " << ((Wt)2 * e.val) + k << "\n";);
        delta_ub.push_back({{e.vert, e.vert - 1}, ((Wt)2 * e.val) + k});
      } else {
        // CASE D
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-lb-d) add edge: " << e.vert << " --> "
                     << e.vert + 1 << " with " << ((Wt)2 * e.val) + k << "\n";);
        delta_lb.push_back({{e.vert, e.vert + 1}, ((Wt)2 * e.val) + k});
      }
    }

    if (!update_delta(m_graph, m_potential, delta_lb)) {
      set_to_bottom();
      return false;
    }

    if (!update_delta(m_graph, m_potential, delta_ub)) {
      set_to_bottom();
      return false;
    }

    return true;
  }

  /*
   * Update other bounds after an upper-bound interval constraint has
   * been added.
   *                                k
   * We just added the edge neg(v) ----> pos(v):

   * (1) This is the first rule of inference:
   *      v <= k
   *   y -v <= k'
   *  -------------
   *   y    <= k'+k
   *
   * There are two cases to implement because there are two ways of
   * representing y-v <= k':
   *
   * CASE A
   *          k            k'             k+k'
   * neg(v)------>pos(v) -----> pos(y) <- - - - neg(y)
   *
   * CASE B
   *        k'            k
   * neg(y)----->neg(v)------>pos(v)       pos(y)
   *   |                                     ^
   *   |          k+k'                       |
   *    - - - - - - - - - - - - - - - - - - -
   *
   * (2) This is the second rule of inference:
   *      v <= k
   *  -v -y <= k'
   *  ------------
   *     -y <= k'+k
   *
   * Again, there are two cases to implement because there are two
   * ways of representing -v-y <= k':
   *
   * CASE C:
   *            k            k'                 k+k'
   *   neg(v) -----> pos(v) -----> neg(y) <- - - - - - pos(y)
   *
   * CASE D
   *           k+k'           k'           k
   *   neg(y)<- - - pos(y) -----> neg(v) -----> pos(v)
   *
  */
  bool update_bounds_ub(vert_id v, Wt k) {
    assert(v % 2 == 0);
    // we add in delta so that we don't invalidate graph iterators
    edge_vector delta_lb, delta_ub;

    for (auto e : m_graph.e_succs(v)) {
      if (e.vert % 2 == 0) {
        // CASE A
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-ub-a) add edge " << e.vert + 1 << " --> "
                     << e.vert << " with " << ((Wt)2 * e.val) + k << "\n";);
        delta_ub.push_back({{e.vert + 1, e.vert}, (Wt)2 * e.val + k});
      } else {
        // CASE C
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-ub-c) add edge " << e.vert - 1 << " --> "
                     << e.vert << " with " << ((Wt)2 * e.val) + k << "\n";);
        delta_lb.push_back({{e.vert - 1, e.vert}, (Wt)2 * e.val + k});
      }
    }
    check_potential(m_graph, m_potential, __LINE__);

    for (auto e : m_graph.e_preds(v + 1)) {
      if (e.vert % 2 != 0) {
        // CASE B
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-ub-b) add edge: " << e.vert << " --> "
                     << e.vert - 1 << " with " << (Wt)2 * e.val + k << "\n";);
        delta_ub.push_back({{e.vert, e.vert - 1}, (Wt)2 * e.val + k});
      } else {
        // CASE D
        CRAB_LOG("octagon-bounds",
                 crab::outs()
                     << "\t"
                     << "(close-bounds-ub-d) add edge: " << e.vert << " --> "
                     << e.vert + 1 << " with " << (Wt)2 * e.val + k << "\n";);
        delta_lb.push_back({{e.vert, e.vert + 1}, (Wt)2 * e.val + k});
      }
    }

    if (!update_delta(m_graph, m_potential, delta_lb)) {
      set_to_bottom();
      return false;
    }

    if (!update_delta(m_graph, m_potential, delta_ub)) {
      set_to_bottom();
      return false;
    }

    // assert(check_potential(m_graph, m_potential, __LINE__));
    return true;
  }

  // Update bounds after adding an edge from i to j with weight c.
  void update_bounds(graph_t &g, vert_id i, vert_id j, Wt c,
                     edge_vector &delta) {
    CRAB_LOG("octagon-bounds",
             crab::outs() << "Updating bounds after adding an edge from " << i
                          << " to " << j << " with k=" << c << "...\n";);

    wt_ref_t w;
    if (i % 2 == 0 && j % 2 == 0) {
      /* as inference rule:
         j - i <= c
         j + i <= w
         -------------
         j     <=  (c+w)/2
       */
      if (g.lookup(j + 1, i, w)) {
        Wt k = (w.get() + c); // don't multiply by 2
        if (!g.lookup(j + 1, j, w) || k < w.get()) {
          delta.push_back({{j + 1, j}, k}); // ub of j
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 1)] added edge " << j + 1
                       << " --> " << j << " with " << k << "\n";);
        }
      }
      /* as inference rule:
          j - i <= c
         -j - i <= w
         -------------
              -i <= (c+w)/2
       */
      if (g.lookup(j, i + 1, w)) {
        Wt k = (w.get() + c); // don't multiply by -2
        if (!g.lookup(i, i + 1, w) || k < w.get()) {
          delta.push_back({{i, i + 1}, k}); // lb of i
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 2)] added edge " << i
                       << " --> " << i + 1 << " with " << k << "\n";);
        }
      }
    } else if (i % 2 == 0 && j % 2 != 0) {
      /* as inference rule:
         -j - i <= c
         -j + i <= w
         -------------
          -j    <=  (c+w)/2
       */
      if (g.lookup(j - 1, i, w)) {
        Wt k = (w.get() + c); // don't multiply by -2
        if (!g.lookup(j - 1, j, w) || k < w.get()) {
          delta.push_back({{j - 1, j}, k}); // lb of j
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 3)] added edge " << j - 1
                       << " --> " << j << " with " << k << "\n";);
        }
      }
      /* as inference rule:
         -j - i <= c
         +j - i <= w
         -------------
            -i  <=  (c+w)/2
       */
      if (g.lookup(j, i + 1, w)) {
        Wt k = (w.get() + c); // don't multiply by -2
        if (!g.lookup(i, i + 1, w) || k < w.get()) {
          delta.push_back({{i, i + 1}, k}); // lb of i
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 4)] added edge " << i
                       << " --> " << i + 1 << " with " << k << "\n";);
        }
      }
    } else if (i % 2 != 0 && j % 2 == 0) {
      /* as inference rule:
          j + i <= c
         +j - i <= w
         -------------
          j      <= (c+w)/2
       */
      if (g.lookup(j + 1, i, w)) {
        Wt k = (w.get() + c); // don't multiply by 2
        if (!g.lookup(j + 1, j, w) || k < w.get()) {
          delta.push_back({{j + 1, j}, k}); // ub of j
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 5)] added edge " << j + 1
                       << " --> " << j << " with " << k << "\n";);
        }
      }
      /* as inference rule:
          j + i <= c
         -j + i <= w
         -------------
               i <= (c+w)/2
       */
      if (g.lookup(j, i - 1, w)) {
        Wt k = (w.get() + c); // don't multiply by 2
        if (!g.lookup(i, i - 1, w) || k < w.get()) {
          delta.push_back({{i, i - 1}, k}); // ub of i
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 6)] added edge " << i
                       << " --> " << i - 1 << " with " << k << "\n";);
        }
      }
    } else if (i % 2 != 0 && j % 2 != 0) {
      /* as inference rule:
         -j + i <= c
         -j - i <= w
         -------------
          -j    <= (c+w)/2
       */
      if (g.lookup(j - 1, i, w)) {
        Wt k = (w.get() + c); // don't multiply by -2
        if (!g.lookup(j - 1, j, w) || k < w.get()) {
          delta.push_back({{j - 1, j}, k}); // lb of j
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 7)] added edge " << j - 1
                       << " --> " << j << " with " << k << "\n";);
        }
      }
      /* as inference rule:
         -j + i <= c
          j + i <= w
         -------------
               i <= (c+w)/2
       */
      if (g.lookup(j, i - 1, w)) {
        Wt k = (w.get() + c); // don't multiply by 2
        if (!g.lookup(i, i - 1, w) || k < w.get()) {
          delta.push_back({{i, i - 1}, k}); // ub of i
          CRAB_LOG("octagon-bounds",
                   crab::outs()
                       << "[close-edge (update_bounds 8)] added edge " << i
                       << " --> " << i - 1 << " with " << k << "\n";);
        }
      }
    }
  }

  bool add_linear_leq(const linear_expression_t &exp) {
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".add_constraints.add_linear_leq");
    normalize();
    CRAB_LOG("octagon-add", linear_expression_t exp_tmp(exp);
             crab::outs() << "Begin adding: " << exp_tmp << "<= 0"
                          << " to:\n"
                          << *this << "\n"
                          << "graph=" << m_graph << "\n";
             for (vert_id v
                  : m_graph.verts()) {
               if (auto v_name = m_rev_vert_map[v]) {
                 crab::outs() << "\tpot[" << *v_name << "=" << v
                              << "]=" << m_potential[v] << "\n";
               } else {
                 crab::outs()
                     << "\tpot[" << v << "]=" << m_potential[v] << "\n";
               }
             });

    // x <= k constraints
    std::vector<std::pair<variable_t, Wt>> lbs, ubs;
    //  x-y constraints
    std::vector<diffcst_t> csts;
    //  x+y <=k constraints
    std::vector<diffcst_t> pos_csts;
    // -x-y <=k constraints
    std::vector<diffcst_t> neg_csts;
    oct_csts_of_lin_leq(exp, csts, pos_csts, neg_csts, lbs, ubs);

    check_potential(m_graph, m_potential, __LINE__);
    Wt_min min_op;
    wt_ref_t w;
    for (auto p : lbs) {
      crab::ScopedCrabStats __st__(domain_name() +
                                   ".add_constraints.add_lb_csts");
      CRAB_LOG("octagon-add",
               crab::outs() << "\t" << p.first << " >= " << p.second << "\n";);

      variable_t x(p.first);
      vert_id v = get_vert(p.first);
      if (m_graph.lookup(v, v + 1, w) && w.get() <= (Wt)-2 * p.second)
        continue;
      m_graph.set_edge(v, (Wt)-2 * p.second, v + 1);
      CRAB_LOG("octagon-add", crab::outs()
                                  << "\t1 Add edge: " << v << " --> " << v + 1
                                  << " with " << (Wt)-2 * p.second << "\n";);

      if (!repair_potential(v, v + 1)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 1\n";);
        return false;
      }

      check_potential(m_graph, m_potential, __LINE__);

      // Update other bounds
      if (!update_bounds_lb(v, (Wt)-2 * p.second)) {
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 2\n";);
        return false;
      }
    }

    for (auto p : ubs) {
      crab::ScopedCrabStats __st__(domain_name() +
                                   ".add_constraints.add_ub_csts");
      CRAB_LOG("octagon-add",
               crab::outs() << "\t" << p.first << "<=" << p.second << "\n");
      variable_t x(p.first);
      vert_id v = get_vert(p.first);
      if (m_graph.lookup(v + 1, v, w) && w.get() <= (Wt)2 * p.second)
        continue;
      m_graph.set_edge(v + 1, (Wt)2 * p.second, v);
      CRAB_LOG("octagon-add", crab::outs()
                                  << "\t3 Add edge: " << v + 1 << " --> " << v
                                  << " with " << (Wt)2 * p.second << "\n";);
      if (!repair_potential(v + 1, v)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 4\n";);
        return false;
      }
      check_potential(m_graph, m_potential, __LINE__);

      // Update other bounds
      if (!update_bounds_ub(v, (Wt)2 * p.second)) {
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 5\n";);
        return false;
      }
    }

    for (auto diff : csts) {
      crab::ScopedCrabStats __st__(domain_name() +
                                   ".add_constraints.add_diff_csts");
      CRAB_LOG("octagon-add", crab::outs() << "\t" << diff.first.first << "-"
                                           << diff.first.second
                                           << "<=" << diff.second << "\n");

      vert_id src = get_vert(diff.first.second);
      vert_id dest = get_vert(diff.first.first);
      Wt w = diff.second;

      /*
        Check if the binary constraint already exists via bounds:
           x - y <= k is skipped if x <= k1 and -y <= k2 and k1+k2 <= k
       */
      wt_ref_t w1, w2;
      if (m_graph.lookup(dest + 1, dest, w1) &&
          m_graph.lookup(src, src + 1, w2) &&
          ((w1.get() + w2.get()) / (Wt)2) <= w) {
        continue;
      }

      m_graph.update_edge(src, w, dest, min_op);
      CRAB_LOG("octagon-add", crab::outs() << "\t5 Add edge: " << src << " --> "
                                           << dest << " with " << w << "\n";);
      if (!repair_potential(src, dest)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 7\n";);
        return false;
      }
      check_potential(m_graph, m_potential, __LINE__);
      close_over_edge(src, dest);
      check_potential(m_graph, m_potential, __LINE__);

      /* begin preserve coherence */
      m_graph.update_edge(dest + 1, w, src + 1, min_op);
      CRAB_LOG("octagon-add", crab::outs()
                                  << "\t6 Add edge: " << dest + 1 << " --> "
                                  << src + 1 << " with " << w << "\n";);
      if (!repair_potential(dest + 1, src + 1)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 8\n";);
        return false;
      }
      check_potential(m_graph, m_potential, __LINE__);
      close_over_edge(dest + 1, src + 1);
      check_potential(m_graph, m_potential, __LINE__);
      /* end preserve coherence*/
    }

    for (auto diff : pos_csts) {
      crab::ScopedCrabStats __st__(domain_name() +
                                   ".add_constraints.add_pos_csts");
      CRAB_LOG("octagon-add", crab::outs() << "\t" << diff.first.first << "+"
                                           << diff.first.second
                                           << "<=" << diff.second << "\n");

      vert_id src = get_vert(diff.first.second);
      vert_id dest = get_vert(diff.first.first);
      Wt w = diff.second;

      /*
        Check if the binary constraint already exists via bounds:
           x + y <= k is skipped if x <= k1 and y <= k2 and k1+k2 <= k
       */
      wt_ref_t w1, w2;
      if (m_graph.lookup(dest + 1, dest, w1) &&
          m_graph.lookup(src + 1, src, w2) &&
          ((w1.get() + w2.get()) / (Wt)2) <= w) {
        continue;
      }

      vert_id x = src + 1;
      vert_id y = dest;
      m_graph.update_edge(x, w, y, min_op);
      CRAB_LOG("octagon-add", crab::outs() << "\t7 Add edge: " << x << " --> "
                                           << y << " with " << w << "\n";);
      if (!repair_potential(x, y)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 9\n";);
        return false;
      }
      check_potential(m_graph, m_potential, __LINE__);
      close_over_edge(x, y);
      check_potential(m_graph, m_potential, __LINE__);

      /* begin preserve coherence */
      x = dest + 1;
      y = src;
      m_graph.update_edge(x, w, y, min_op);
      CRAB_LOG("octagon-add", crab::outs() << "\t8 Add edge: " << x << " --> "
                                           << y << " with " << w << "\n";);

      if (!repair_potential(x, y)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 10\n";);
        return false;
      }
      check_potential(m_graph, m_potential, __LINE__);
      close_over_edge(x, y);
      check_potential(m_graph, m_potential, __LINE__);
      /* end preserve coherence */
    }

    for (auto diff : neg_csts) {
      crab::ScopedCrabStats __st__(domain_name() +
                                   ".add_constraints.add_neg_csts");
      CRAB_LOG("octagon-add", crab::outs() << "\t-" << diff.first.first << "-"
                                           << diff.first.second
                                           << "<=" << diff.second << "\n");
      check_potential(m_graph, m_potential, __LINE__);
      vert_id src = get_vert(diff.first.second);
      vert_id dest = get_vert(diff.first.first);
      Wt w = diff.second;

      /*
        Check if the binary constraint already exists via bounds:
          -x - y <= k is skipped if -x <= k1 and -y <= k2 and k1+k2 <= k
       */
      wt_ref_t w1, w2;
      if (m_graph.lookup(dest, dest + 1, w1) &&
          m_graph.lookup(src, src + 1, w2) &&
          ((w1.get() + w2.get()) / (Wt)2) <= w) {
        continue;
      }

      vert_id x = src;
      vert_id y = dest + 1;
      m_graph.update_edge(x, w, y, min_op);
      CRAB_LOG("octagon-add", crab::outs() << "\t9 Add edge: " << x << " --> "
                                           << y << " with " << w << "\n";);
      if (!repair_potential(x, y)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 11\n";);
        return false;
      }
      check_potential(m_graph, m_potential, __LINE__);
      close_over_edge(x, y);
      check_potential(m_graph, m_potential, __LINE__);

      /* begin preserve coherence */
      x = dest;
      y = src + 1;
      m_graph.update_edge(x, w, y, min_op);
      CRAB_LOG("octagon-add", crab::outs() << "\t10 Add edge: " << x << " --> "
                                           << y << " with " << w << "\n";);
      if (!repair_potential(x, y)) {
        set_to_bottom();
        CRAB_LOG("octagon-add", crab::outs() << "\tDetected bottom 12\n";);
        return false;
      }
      check_potential(m_graph, m_potential, __LINE__);
      close_over_edge(x, y);
      check_potential(m_graph, m_potential, __LINE__);
      /* end preserve coherence */
    }

    integer_tightening();

    check_potential(m_graph, m_potential, __LINE__);
    CRAB_LOG("octagon-add", crab::outs() << "End adding: " << *this << "\n"
                                         << "with graph" << m_graph << "\n";
             for (vert_id v
                  : m_graph.verts()) {
               crab::outs() << "\tpot[" << v << "]=" << m_potential[v] << "\n";
             });
    return true;
  }

  bool add_univar_disequation(const variable_t &x, number_t n) {
    bool overflow;
    interval_t i = get_interval(x);
    interval_t ni(n);
    interval_t new_i =
        ikos::linear_interval_solver_impl::trim_interval<interval_t>(
            i, ni);
    CRAB_LOG("octagon", crab::outs()
                            << "Adding disequation: " << x << "!=" << n << "\n"
                            << new_i << "\n");

    Wt_min min_op;
    if (new_i.is_bottom()) {
      set_to_bottom();
      return false;
    } else if (!new_i.is_top() && (new_i <= i)) {
      vert_id v = get_vert(x);
      wt_ref_t w;
      if (new_i.lb().is_finite()) {
        Wt lb_val =
            ntow::convert(-number_t(2) * (*(new_i.lb().number())), overflow);
        if (overflow) {
          return true;
        }
        if (m_graph.lookup(v, v + 1, w) && lb_val < w.get()) {
          m_graph.set_edge(v, lb_val, v + 1);
          if (!repair_potential(v, v + 1)) {
            set_to_bottom();
            return false;
          }
          check_potential(m_graph, m_potential, __LINE__);
          if (!update_bounds_lb(v, lb_val)) {
            return false;
          }
        }
      }

      if (new_i.ub().is_finite()) {
        Wt ub_val =
            ntow::convert(number_t(2) * (*(new_i.ub().number())), overflow);
        if (overflow) {
          return true;
        }
        if (m_graph.lookup(v + 1, v, w) && (ub_val < w.get())) {
          m_graph.set_edge(v + 1, ub_val, v);
          if (!repair_potential(v + 1, v)) {
            set_to_bottom();
            return false;
          }
          check_potential(m_graph, m_potential, __LINE__);
          if (!update_bounds_ub(v, ub_val)) {
            return false;
          }
        }
      }
    }

    integer_tightening();
    return true;
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

  void add_disequation(const linear_expression_t &e) {
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".add_constraints.add_disequation");
    for (auto kv : e) {
      variable_t pivot = kv.second;
      interval_t i = compute_residual(e, pivot) / interval_t(kv.first);
      if (auto k = i.singleton()) {
        if (!add_univar_disequation(pivot, *k)) {
          // already set to bottom
          return;
        }
      }
    }
    return;
  }

  interval_t get_interval(const variable_t &x) const {
    return get_interval(m_vert_map, m_graph, x);
  }

  interval_t get_interval(const vert_map_t &m, const graph_t &r, const variable_t &x) const {
                          
    auto it = m.find(x);
    if (it == m.end())
      return interval_t::top();
    vert_id v = (*it).second.first;
    interval_t x_out =
        interval_t(r.elem(v, v + 1) ? -number_t(r.edge_val(v, v + 1)) / 2
                                    : bound_t::minus_infinity(),
                   r.elem(v + 1, v) ? number_t(r.edge_val(v + 1, v)) / 2
                                    : bound_t::plus_infinity());
    return x_out;
  }

  // FIXME HACK: big assumption that Wt can be casted to float. This
  // is true if Wt is long but not if, for instance, we use bignums.
  void integer_tightening() {
#ifdef INTEGER_TIGHTENING
    for (vert_id v : m_graph.verts()) {
      wt_ref_t w;
      if (v % 2 == 0 && (m_graph.lookup(v, v + 1, w)) && (w.get() & 1)) {
        // weight is odd so we need to tighten it
        CRAB_LOG("octagon-integer", crab::outs() << "Tightening edge " << v
                                                 << " --> " << v + 1 << " from "
                                                 << w.get() << " to ";);
	// REVISIT(PERFORMANCE): extra call to lookup
	Wt tightened_w = 2 * (Wt)std::floor((float)w.get() / 2);
	m_graph.set_edge(v, tightened_w, v + 1);
        CRAB_LOG("octagon-integer", crab::outs() << tightened_w << "\n";);
      }
      if (v % 2 != 0 && (m_graph.lookup(v, v - 1, w)) && (w.get() & 1)) {
        // weight is odd so we need to tighten it
        CRAB_LOG("octagon-integer", crab::outs() << "Tightening edge " << v
                                                 << " --> " << v + 1 << " from "
                                                 << w.get() << " to ";);
	// REVISIT(PERFORMANCE): extra call to lookup
	Wt tightened_w = 2 * (Wt)std::floor((float)w.get() / 2);
	m_graph.set_edge(v, tightened_w, v - 1);
        CRAB_LOG("octagon-integer", crab::outs() << tightened_w << "\n";);
      }
    }
#endif
  }

  void close_over_edge(vert_id ii, vert_id jj) {
    assert(ii / 2 != jj / 2);

    wt_ref_t w;
    Wt c = m_graph.edge_val(ii, jj);

    split_octagons_impl::SplitOctGraph<graph_t> g_oct(m_graph);
    // we add in delta so that we don't invalidate graph iterators
    edge_vector delta;

    update_bounds(m_graph, ii, jj, c, delta);
    GrOps::apply_delta(m_graph, delta);
    delta.clear();
    CRAB_LOG("octagon-close-edge", crab::outs()
                                       << "After adding " << ii << "-->" << jj
                                       << "w=" << c << ": ";
             m_graph.write(crab::outs()); crab::outs() << "\n";);

    std::vector<std::pair<vert_id, Wt>> src_dec;
    for (auto edge : g_oct.e_preds(ii)) {
      vert_id se = edge.vert;
      Wt wt_sij = edge.val + c;
      assert(g_oct.succs(se).begin() != g_oct.succs(se).end());
      if (se != jj) {
        if (/*g_oct*/ m_graph.lookup(se, jj, w)) {
          if (w.get() <= wt_sij)
            continue;
          CRAB_LOG("octagon-close-edge",
                   crab::outs()
                       << "[close-edge 1] improved edge " << se << "-->" << jj
                       << " from " << w.get() << " to " << wt_sij << "\n";);
	  // REVISIT(PERFORMANCE): extra call to lookup
	  m_graph.set_edge(se, wt_sij, jj);
        } else {
          CRAB_LOG("octagon-close-edge",
                   crab::outs() << "[close-edge 1] added edge " << se << "-->"
                                << jj << " with w= " << wt_sij << "\n";);
          delta.push_back({{se, jj}, wt_sij});
        }
        src_dec.push_back({se, edge.val});
        update_bounds(m_graph, se, jj, wt_sij, delta);
      }
    }
    GrOps::apply_delta(m_graph, delta);
    delta.clear();
    CRAB_LOG("octagon-close-edge",
             crab::outs() << "After closing all predecessors of " << ii << ": ";
             m_graph.write(crab::outs()); crab::outs() << "\n";);

    std::vector<std::pair<vert_id, Wt>> dest_dec;
    for (auto edge : g_oct.e_succs(jj)) {
      vert_id de = edge.vert;
      Wt wt_ijd = edge.val + c;
      if (de != ii) {
        if (/*g_oct*/ m_graph.lookup(ii, de, w)) {
          if (w.get() <= wt_ijd)
            continue;
          CRAB_LOG("octagon-close-edge",
                   crab::outs()
                       << "[close-edge 2] improved edge " << ii << "-->" << de
                       << " from " << w.get() << " to " << wt_ijd << "\n";);
	  // REVISIT(PERFORMANCE): extra call to lookup
	  m_graph.set_edge(ii, wt_ijd, de);
        } else {
          delta.push_back({{ii, de}, wt_ijd});
          CRAB_LOG("octagon-close-edge",
                   crab::outs() << "[close-edge 2] added edge " << ii << "-->"
                                << de << " with w= " << wt_ijd << "\n";);
        }
        dest_dec.push_back({de, edge.val});
        update_bounds(m_graph, ii, de, wt_ijd, delta);
      }
    }

    GrOps::apply_delta(m_graph, delta);
    delta.clear();
    CRAB_LOG("octagon-close-edge",
             crab::outs() << "After closing all successors of " << jj << ": ";
             m_graph.write(crab::outs()); crab::outs() << "\n";);

    for (auto s_p : src_dec) {
      vert_id se = s_p.first;
      Wt wt_sij = c + s_p.second;
      for (auto d_p : dest_dec) {
        vert_id de = d_p.first;
        if (se == de)
          continue;
        Wt wt_sijd = wt_sij + d_p.second;
        if (m_graph.lookup(se, de, w)) {
          if (w.get() <= wt_sijd)
            continue;
          CRAB_LOG("octagon-close-edge",
                   crab::outs()
                       << "[close-edge 3] improved edge " << se << "-->" << de
                       << " from " << w.get() << " to " << wt_sijd << "\n";);
	  // REVISIT(PERFORMANCE): extra call to lookup
	  m_graph.set_edge(se, wt_sijd, de);
        } else {
          // m_graph.add_edge(se, wt_sijd, de);
          delta.push_back({{se, de}, wt_sijd});
          CRAB_LOG("octagon-close-edge",
                   crab::outs() << "[close-edge 3] added edge " << se << "-->"
                                << de << " with w= " << wt_sijd << "\n";);
        }
        update_bounds(m_graph, se, de, wt_sijd, delta);
      }
    }
    GrOps::apply_delta(m_graph, delta);
    delta.clear();
    CRAB_LOG("octagon-close-edge",
             crab::outs() << "After closing paths between predecessors of "
                          << ii << " and successors of " << jj << ": ";
             m_graph.write(crab::outs()); crab::outs() << "\n";);
  }

  // helper for the join operation.
  template <typename Gr> graph_t split_rels(Gr &gx, Gr &gy, unsigned sz) const {
    crab::CrabStats::count(domain_name() + ".count.join.split_rels");
    crab::ScopedCrabStats __st__(domain_name() + ".join.split_rels");

    // Compute the deferred relations for gx
    graph_t g_ix_ry;
    g_ix_ry.growTo(sz);

    split_octagons_impl::SplitOctGraph<GrPerm> gy_oct(gy);
    for (vert_id s : gy_oct.verts()) {
      for (vert_id d : gy_oct.succs(s)) {
        Wt_min min_op;
        wt_ref_t ws, wd;
        bool added = false;
        if (s % 2 == 0 && d % 2 == 0) { // both pos
          // if  d-s <= k in gy and -2s <= k1  and 2d <= k2 in gx then
          //     add d-s <= k1+k2/2
          if (gx.lookup(s, s + 1, ws) && gx.lookup(d + 1, d, wd)) {
            added = true;
          }
        } else if (s % 2 == 0 && d % 2 != 0) { // s pos
          // if -d-s <= k in gy and -2s <= k1 and -2d <= k2 in gx then
          //    add -d-s <= k1+k2/2
          if (gx.lookup(s, s + 1, ws) && gx.lookup(d - 1, d, wd)) {
            added = true;
          }
        } else if (s % 2 != 0 && d % 2 == 0) { // d pos
          // if d+s <= k in gy and 2s <= k1 and 2d <= k2 in gx then
          //    add d+s <= k1+k2/2
          if (gx.lookup(s, s - 1, ws) && gx.lookup(d + 1, d, wd)) {
            added = true;
          }
        } else if (s % 2 != 0 && d % 2 != 0) { // both neg
          // if -d+s <= k in gy and 2s <= k1 and -2d <= k2 in gx then
          //    add -d+s <= k1+k2/2
          if (gx.lookup(s, s - 1, ws) && gx.lookup(d - 1, d, wd)) {
            added = true;
          }
        }
        if (added) {
          g_ix_ry.update_edge(s, (Wt)(ws.get() + wd.get()) / (Wt)2, d, min_op);
        }
      }
    }
    return g_ix_ry;
  }

  /**
   * Compute stable relations by considering explicit edges in the
   * left operand and implicit one in the right operand.
   *
   * g, l, and r are not modified so they should be "const".
   **/
  template <typename Gr>
  void split_widen_rels(graph_t &g, Gr &l, Gr &r, edge_vector &delta) const {

    Wt_min min_op;
    Wt e_val;
    wt_ref_t ws, wd, wy;

    // Skip bounds on the non-closed left
    split_octagons_impl::SplitOctGraph<GrPerm> l_oct(l);

    for (vert_id s : l_oct.verts()) {
      for (vert_id d : l_oct.succs(s)) {
        e_val = l.edge_val(s, d);
        if (s % 2 == 0 && d % 2 == 0) { // both pos
          // if d-s <= k in l and -2s <= k1 and 2d <= k2 in r then
          //    add d-s <= e_val if k1+k2/2 <= e_val
          if (r.lookup(s, s + 1, ws) && r.lookup(d + 1, d, wd) &&
              ((Wt)(ws.get() + wd.get()) / (Wt)2 <= e_val) &&
              (!g.lookup(s, d, wy) || (e_val < wy.get()))) {
            /* begin preserve coherence */
            if (g.elem(d + 1, s + 1)) {
              continue;
            }
            /* end preserve coherence */
            CRAB_LOG("octagon-widening", auto vs = m_rev_vert_map[s];
                     auto vd = m_rev_vert_map[d];
                     crab::outs() << "Widening implicit added " << *vd << "-"
                                  << *vs << "<=" << e_val << "  " << s << "-->"
                                  << d << "=" << e_val << "\n";);
            // g.update_edge(s, e_val, d, min_op);
            delta.push_back({{s, d}, e_val});
          }
        } else if (s % 2 == 0 && d % 2 != 0) { // s pos
          // if -d-s <= k in l and -2s <= k1 and -2d <= k2 in r then
          //    add -d-s <= e_val if k1+k2/2 <= e_val
          if (r.lookup(s, s + 1, ws) && r.lookup(d - 1, d, wd) &&
              ((Wt)(ws.get() + wd.get()) / (Wt)2 <= e_val) &&
              (!g.lookup(s, d, wy) || (e_val < wy.get()))) {
            /* begin preserve coherence */
            if (g.elem(d - 1, s + 1)) {
              continue;
            }
            /* end preserve coherence */
            CRAB_LOG("octagon-widening", auto vs = m_rev_vert_map[s];
                     auto vd = m_rev_vert_map[d];
                     crab::outs()
                     << "Widening implicit added "
                     << "-" << *vd << "-" << *vs << "<=" << e_val << "  " << s
                     << "-->" << d << "=" << e_val << "\n";);
            // g.update_edge(s, e_val, d, min_op);
            delta.push_back({{s, d}, e_val});
          }
        } else if (s % 2 != 0 && d % 2 == 0) { // d pos
          // if d+s <= k in l and 2s <= k1 and 2d <= k2 in r then
          //    add d+s <= e_val if k1+k2/2  <= e_val
          if (r.lookup(s, s - 1, ws) && r.lookup(d + 1, d, wd) &&
              ((Wt)(ws.get() + wd.get()) / (Wt)2 <= e_val) &&
              (!g.lookup(s, d, wy) || (e_val < wy.get()))) {
            /* begin preserve coherence */
            if (g.elem(d + 1, s - 1)) {
              continue;
            }
            /* end preserve coherence */
            CRAB_LOG("octagon-widening", auto vs = m_rev_vert_map[s];
                     auto vd = m_rev_vert_map[d];
                     crab::outs() << "Widening implicit added " << *vd << "+"
                                  << *vs << "<=" << e_val << "  " << s << "-->"
                                  << d << "=" << e_val << "\n";);
            // g.update_edge(s, e_val, d, min_op);
            delta.push_back({{s, d}, e_val});
          }
        } else if (s % 2 != 0 && d % 2 != 0) { // both neg
          // if -d+s <= k in l and 2s <= k1 and -2d <= k2 in r then
          //    add -d+s <= e_val if k1+k2/2  <= e_val
          if (r.lookup(s, s - 1, ws) && r.lookup(d - 1, d, wd) &&
              ((Wt)(ws.get() + wd.get()) / (Wt)2 <= e_val) &&
              (!g.lookup(s, d, wy) || (e_val < wy.get()))) {
            /* begin preserve coherence */
            if (g.elem(d - 1, s - 1)) {
              continue;
            }
            /* end preserve coherence */
            CRAB_LOG("octagon-widening", auto vs = m_rev_vert_map[s];
                     auto vd = m_rev_vert_map[d];
                     crab::outs()
                     << "Widening implicit added "
                     << "-" << *vd << "+" << *vs << "<=" << e_val << "  " << s
                     << "-->" << d << "=" << e_val << "\n";);
            // g.update_edge(s, e_val, d, min_op);
            delta.push_back({{s, d}, e_val});
          }
        }
      }
    }
  }

  template <class G1, class G2>
  graph_t split_widen(G1 &l, G2 &r, vert_set_t &unstable) const {
    assert(l.size() == r.size());
    size_t sz = l.size();
    graph_t g;
    g.growTo(sz);

    Wt_min min_op;
    wt_ref_t wx;

    /**
     * Add stable relationships which are explicit both in l and r
     **/
    for (vert_id s : r.verts()) {
      for (auto e : r.e_succs(s)) {
        vert_id d = e.vert;
        /* check if s->d is explicit in l */
        if (l.lookup(s, d, wx) && e.val <= wx.get()) {
          g.add_edge(s, wx.get(), d);
          CRAB_LOG(
              "octagon-widening", auto vs = m_rev_vert_map[s];
              auto vd = m_rev_vert_map[d]; if (s % 2 == 0 && s + 1 == d) {
                crab::outs() << "Widening added "
                             << "-" << *vd << "<=" << wx.get() / (Wt)2 << "   "
                             << s << "-->" << d << "=" << wx.get() << "\n";
              } else if (s % 2 != 0 && s == d + 1) {
                crab::outs()
                    << "Widening added " << *vs << "<=" << wx.get() / (Wt)2
                    << "   " << s << "-->" << d << "=" << wx.get() << "\n";
              } else if (s % 2 == 0 && d % 2 == 0) {
                crab::outs() << "Widening added " << *vd << "-" << *vs
                             << "<=" << wx.get() << "   " << s << "-->" << d
                             << "=" << wx.get() << "\n";
              } else if (s % 2 == 0 && d % 2 != 0) {
                crab::outs()
                    << "Widening added "
                    << "-" << *vd << "-" << *vs << "<=" << wx.get() << " " << s
                    << "-->" << d << "=" << wx.get() << "\n";
              } else if (s % 2 != 0 && d % 2 == 0) {
                crab::outs() << "Widening added " << *vd << "+" << *vs
                             << "<=" << wx.get() << "   " << s << "-->" << d
                             << "=" << wx.get() << "\n";
              } else {
                assert(s % 2 != 0 && d % 2 != 0);
                crab::outs()
                    << "Widening added "
                    << "-" << *vd << "+" << *vs << "<=" << wx.get() << "   "
                    << s << "-->" << d << "=" << wx.get() << "\n";
              });
        }
      }
    }

    /**
     * Add stable relationships which are implicit in r and explicit
     * in l
     **/
    edge_vector delta;
    split_widen_rels(g, l, r, delta);
    update_delta(g, delta);
    delta.clear();

    /**
     * Remember vertices that need to be re-closed
     **/
    for (vert_id s : l.verts()) {
      for (vert_id d : l.succs(s)) {
        if (!g.elem(s, d)) {
          unstable.insert(s);
          break;
        }
      }
    }

    /**
     * Add stable relationships which are explicit in r and implicit in l.
     * Special care to ensure termination.
     **/
    split_widen_rels(g, r, l, delta);

    while (!delta.empty()) {
      auto p = delta.back();
      delta.pop_back();
      vert_id s = p.first.first;
      vert_id d = p.first.second;
      Wt w = p.second;
      // Important for termination: only if either src or dst is
      // unstable we add the binary stable edge.  In that way, we
      // ensure that we can only add the binary stable edge once.
      if (unstable.count(s) <= 0 && unstable.count(d) <= 0) {
        CRAB_LOG(
            "octagon-widening", auto vs = m_rev_vert_map[s];
            auto vd = m_rev_vert_map[d]; if (s % 2 == 0 && s + 1 == d) {
              crab::outs() << "Widening discarded "
                           << "-" << *vd << "<=" << wx.get() / (Wt)2 << "   "
                           << s << "-->" << d << "=" << wx.get() << "\n";
            } else if (s % 2 != 0 && s == d + 1) {
              crab::outs() << "Widening discarded " << *vs
                           << "<=" << wx.get() / (Wt)2 << "   " << s << "-->"
                           << d << "=" << wx.get() << "\n";
            } else if (s % 2 == 0 && d % 2 == 0) {
              crab::outs() << "Widening discarded " << *vd << "-" << *vs
                           << "<=" << wx.get() << "   " << s << "-->" << d
                           << "=" << wx.get() << "\n";
            } else if (s % 2 == 0 && d % 2 != 0) {
              crab::outs() << "Widening discarded "
                           << "-" << *vd << "-" << *vs << "<=" << wx.get()
                           << " " << s << "-->" << d << "=" << wx.get() << "\n";
            } else if (s % 2 != 0 && d % 2 == 0) {
              crab::outs() << "Widening discarded " << *vd << "+" << *vs
                           << "<=" << wx.get() << "   " << s << "-->" << d
                           << "=" << wx.get() << "\n";
            } else {
              assert(s % 2 != 0 && d % 2 != 0);
              crab::outs() << "Widening discarded "
                           << "-" << *vd << "+" << *vs << "<=" << wx.get()
                           << "   " << s << "-->" << d << "=" << wx.get()
                           << "\n";
            });
        continue;
      }
      g.update_edge(s, w, d, min_op);
    }

    // TODOXXX: misisng implicit in l and implicit in r.
    return g;
  }

  bool need_normalization() const {
    return !m_unstable.empty();
  }
  
  split_oct_domain(vert_map_t &&_vert_map, rev_map_t &&_rev_map,
                   graph_t &&_graph, std::vector<Wt> &&_potential,
                   vert_set_t &&_unstable)
      : m_vert_map(std::move(_vert_map)), m_rev_vert_map(std::move(_rev_map)),
        m_graph(std::move(_graph)), m_potential(std::move(_potential)),
        m_unstable(std::move(_unstable)), m_is_bottom(false) {
    assert(m_graph.size() >= 0);
  }

public:
  split_oct_domain(bool is_bottom = false) : m_is_bottom(is_bottom) {}

  split_oct_domain(const split_oct_domain_t &o)
      : m_vert_map(o.m_vert_map), m_rev_vert_map(o.m_rev_vert_map),
        m_graph(o.m_graph), m_potential(o.m_potential),
        m_unstable(o.m_unstable), m_is_bottom(false) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (o.m_is_bottom)
      set_to_bottom();
    if (!m_is_bottom)
      assert(m_graph.size() >= 0);
  }

  split_oct_domain(split_oct_domain_t &&o)
      : m_vert_map(std::move(o.m_vert_map)),
        m_rev_vert_map(std::move(o.m_rev_vert_map)),
        m_graph(std::move(o.m_graph)), m_potential(std::move(o.m_potential)),
        m_unstable(std::move(o.m_unstable)), m_is_bottom(o.m_is_bottom) {}

  split_oct_domain &operator=(const split_oct_domain &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");

    if (this != &o) {
      if (o.m_is_bottom) {
        set_to_bottom();
      } else {
        m_is_bottom = false;
        m_vert_map = o.m_vert_map;
        m_rev_vert_map = o.m_rev_vert_map;
        m_graph = o.m_graph;
        m_potential = o.m_potential;
        m_unstable = o.m_unstable;
        assert(m_graph.size() >= 0);
      }
    }
    return *this;
  }

  split_oct_domain &operator=(split_oct_domain &&o) {
    if (o.m_is_bottom) {
      set_to_bottom();
    } else {
      m_is_bottom = false;
      m_vert_map = std::move(o.m_vert_map);
      m_rev_vert_map = std::move(o.m_rev_vert_map);
      m_graph = std::move(o.m_graph);
      m_potential = std::move(o.m_potential);
      m_unstable = std::move(o.m_unstable);
    }
    return *this;
  }

  split_oct_domain_t make_top() const override {
    return split_oct_domain_t(false);
  }

  split_oct_domain_t make_bottom() const override {
    return split_oct_domain_t(true);
  }

  void set_to_top() override {
    split_oct_domain_t tmp(false);
    std::swap(*this, tmp);
  }

  void set_to_bottom() override {
    m_vert_map.clear();
    m_rev_vert_map.clear();
    m_graph.clear();
    m_potential.clear();
    m_unstable.clear();
    m_is_bottom = true;
  }

  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override { return (!is_bottom() && m_graph.is_empty()); }

  bool operator<=(const split_oct_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    if (is_bottom())
      return true;
    else if (o.is_bottom())
      return false;
    else if (o.is_top())
      return true;
    else if (is_top())
      return false;
    else {

      auto leq_op =
	[](const split_oct_domain_t&left, const split_oct_domain_t& right) -> bool {
	// left operand is normalized, right doesn't need to.

	wt_ref_t wx, wy, wz;
	
	if (left.m_vert_map.size() < right.m_vert_map.size()) {
	  return false;
	}


	std::vector<unsigned int> vert_renaming(right.m_graph.size(), -1);
	for (auto p : right.m_vert_map) {
	  auto it = left.m_vert_map.find(p.first);
	  if (it == left.m_vert_map.end()) {
	    return false;
	  }
	  if (right.m_graph.succs(p.second.first).size() == 0 &&
	      right.m_graph.succs(p.second.second).size() == 0 &&
	      right.m_graph.preds(p.second.first).size() == 0 &&
	      right.m_graph.preds(p.second.second).size() == 0)
	    continue;
	  
	  vert_renaming[p.second.first] = (*it).second.first;
	  vert_renaming[p.second.second] = (*it).second.second;
	}
	
	assert(left.m_graph.size() >= 0);
	for (vert_id ox : right.m_graph.verts()) {
	  if (right.m_graph.succs(ox).size() == 0)
	    continue;
	  
	  assert(vert_renaming[ox] != -1);
	  vert_id x = vert_renaming[ox];
	  for (auto edge : right.m_graph.e_succs(ox)) {
	    vert_id oy = edge.vert;
	    assert(vert_renaming[oy] != -1);
	    if (ox == oy)
	      continue;
	    vert_id y = vert_renaming[oy];
	    Wt ow = edge.val;
	    
	    // explicit edge on the left operand
	    if (left.m_graph.lookup(x, y, wx) && (wx.get() <= ow))
	      continue;
	    
          // implicit edge on the left operand
	    if (x % 2 == 0 && y % 2 == 0) { // both pos
	      // y-x <= k in o: check if -2x <= k1 and 2y <= k2 in this and
	      //                         k1+k2/2 <= k
	      if (!left.m_graph.lookup(x, x + 1, wx) ||
		  !left.m_graph.lookup(y + 1, y, wy)) {
		return false;
	      }
	      if (!((Wt)(wx.get() + wy.get()) / (Wt)2 <= ow)) {
		return false;
	      }
	      
          } else if (x % 2 == 0 && y % 2 != 0) { // x pos
	      // -y-x <= k in o: check if -2x <= k1 and -2y <= k2 in this and
	      //                         k1+k2/2 <= k
	      if (!left.m_graph.lookup(x, x + 1, wx) ||
		  !left.m_graph.lookup(y - 1, y, wy)) {
		return false;
	      }
	      if (!((Wt)(wx.get() + wy.get()) / (Wt)2 <= ow)) {
		return false;
	      }
	    } else if (x % 2 != 0 && y % 2 == 0) { // y pos
	      //  y+x <= k in o: check if  2x <= k1 and 2y <= k2 in this and
	      //                         k1+k2/2 <= k
	      if (!left.m_graph.lookup(x, x - 1, wx) ||
		  !left.m_graph.lookup(y + 1, y, wy)) {
		return false;
	      }
	      if (!((Wt)(wx.get() + wy.get()) / (Wt)2 <= ow)) {
		return false;
	      }
	    } else if (x % 2 != 0 && y % 2 != 0) { // both neg
	      //  -y+x <= k in o: check if  2x <= k1 and -2y <= k2 in this and
	      //                         k1+k2/2 <= k
	      if (!left.m_graph.lookup(x, x - 1, wx) ||
		  !left.m_graph.lookup(y - 1, y, wy)) {
              return false;
	      }
	      if (!((Wt)(wx.get() + wy.get()) / (Wt)2 <= ow)) {
		return false;
	      }
	    }
	  }
	}
	return true;
      };

      if (need_normalization()) {
	split_oct_domain_t left(*this);
	left.normalize();
	return leq_op(left, o);
      } else {
	return leq_op(*this, o);
      }
    }
  }

  void operator|=(const split_oct_domain_t &o) override { *this = *this | o; }

  split_oct_domain_t operator|(const split_oct_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || o.is_top())
      return o;
    else if (is_top() || o.is_bottom())
      return *this;
    else {
      CRAB_LOG("octagon", crab::outs() << "Before join:\n"
                                       << "DBM 1\n"
                                       << *this << "\n"
                                       << m_graph << "\n"
                                       << "DBM 2\n"
                                       << o << "\n"
	                               << o.m_graph << "\n");

      auto join_op =
	[this](const split_oct_domain_t &left,  const split_oct_domain_t &right) -> split_oct_domain_t {
	  // Both left and right are normalized
	  
	  CRAB_LOG("octagon-join", crab::outs() << "rev_map 1={";
		   for (unsigned i = 0, e = left.m_rev_vert_map.size(); i != e; i++) {
		     if (left.m_rev_vert_map[i])
		       crab::outs() << *(left.m_rev_vert_map[i]) << "(" << i << ");";
		   } crab::outs()
		   << "}\n";
		   crab::outs() << "rev_map 2={"; for (unsigned i = 0,
							 e = right.m_rev_vert_map.size();
						       i != e; i++) {
		     if (right.m_rev_vert_map[i])
		       crab::outs() << *(right.m_rev_vert_map[i]) << "(" << i << ");";
		   } crab::outs() << "}\n";);
	  

	  check_potential(left.m_graph, left.m_potential, __LINE__);
	  check_potential(right.m_graph, right.m_potential, __LINE__);
	  
	  std::vector<vert_id> perm_x, perm_y;
	  std::vector<Wt> pot_rx, pot_ry;
	  vert_map_t out_vmap;
	  rev_map_t out_revmap;
	  
	  for (auto p : left.m_vert_map) {
	    auto it = right.m_vert_map.find(p.first);
	    if (it != right.m_vert_map.end()) {
	      out_vmap.insert(vmap_elt_t(p.first, {perm_x.size(), perm_x.size() + 1}));
	      out_revmap.push_back(p.first);
	      out_revmap.push_back(p.first);
	      pot_rx.push_back(left.m_potential[p.second.first]);
	      pot_rx.push_back(left.m_potential[p.second.second]);
	      pot_ry.push_back(right.m_potential[(*it).second.first]);
	      pot_ry.push_back(right.m_potential[(*it).second.second]);
	      perm_x.push_back(p.second.first);
	      perm_x.push_back(p.second.second);
	      perm_y.push_back((*it).second.first);
	      perm_y.push_back((*it).second.second);
	    }
	  }
	  
	  unsigned int sz = perm_x.size();
	  
	  // Build the permuted view of x and y.
	  assert(left.m_graph.size() > 0);
	  GrPerm gx(perm_x, left.m_graph);
	  assert(right.m_graph.size() > 0);
	  GrPerm gy(perm_y, right.m_graph);
	  
	  // Compute the deferred relations for gy
	  graph_t g_ix_ry = split_rels(gx, gy, sz);
	  // Apply the deferred relations, and re-close.
	  bool is_closed;
	  graph_t g_rx(GrOps::meet(gx, g_ix_ry, is_closed));
#ifdef JOIN_CLOSE_AFTER_MEET
	  // Conjecture: g_rx is closed
	  if (!is_closed) {
	    edge_vector delta;
	    split_octagons_impl::SplitOctGraph<graph_t> g_rx_oct(g_rx);
	    GrOps::close_after_meet(g_rx_oct, pot_rx, gx, g_ix_ry, delta);
	    update_delta(g_rx, delta);
	  }
#endif
	  // Compute the deferred relations for gx
	  graph_t g_rx_iy = split_rels(gy, gx, sz);
	  // Apply the deferred relations, and re-close.
	  graph_t g_ry(GrOps::meet(gy, g_rx_iy, is_closed));
#ifdef JOIN_CLOSE_AFTER_MEET
	  // Conjecture: g_ry is closed
	  if (!is_closed) {
	    edge_vector delta;
	    split_octagons_impl::SplitOctGraph<graph_t> g_ry_oct(g_ry);
	    GrOps::close_after_meet(g_ry_oct, pot_ry, gy, g_rx_iy, delta);
	    update_delta(g_ry, delta);
	  }
#endif
	  CRAB_LOG("octagon-join", crab::outs() << "rev_map={";
		   for (unsigned i = 0, e = out_revmap.size(); i != e; i++) {
		     if (out_revmap[i])
		       crab::outs() << *(out_revmap[i]) << "(" << i << ");";
		   } crab::outs()
		   << "}\n";);
	  
	  CRAB_LOG("octagon-join", crab::outs()
		   << "\tBefore joined:\n\tg_rx:" << g_rx
		   << "\n\tg_ry:" << g_ry << "\n");

	  // We now have the relevant set of relations. Because g_rx and
	  // g_ry are closed, the result is also closed.
	  graph_t join_g(GrOps::join(g_rx, g_ry));
	  
	  CRAB_LOG("octagon-join", crab::outs() << "Joined graph:\n"
		   << join_g << "\n");

	  /*
	   * Infer relationships from bounds
	   */

	  //         a                                                c
	  // pos(x) --> neg(x) in one graph and the other has pos(x) ---> neg(x)
	  // such that a < c
	  std::vector<vert_id> lb_left;
	  //         a                                                c
	  // pos(x) --> neg(x) in one graph and the other has pos(x) ---> neg(x)
	  // such that a > c
	  std::vector<vert_id> lb_right;
	  //         a                                                c
	  // neg(x) --> pos(x) in one graph and the other has neg(x) ---> pos(x)
	  // such that a < c
	  std::vector<vert_id> ub_left;
	  //         a                                                c
	  // neg(x) --> pos(x) in one graph and the other has neg(x) ---> pos(x)
	  // such that a > c
	  std::vector<vert_id> ub_right;

	  wt_ref_t wx;
	  wt_ref_t wy;

	  for (vert_id v : gx.verts()) {
	    if (v % 2 != 0)
	      continue;
	    if (gx.lookup(v + 1, v, wx) && gy.lookup(v + 1, v, wy)) {
	      if (wx.get() < wy.get())
		ub_left.push_back(v);
	      if (wy.get() < wx.get())
		ub_right.push_back(v);
	    }
	    if (gx.lookup(v, v + 1, wx) && gy.lookup(v, v + 1, wy)) {
	      if (wx.get() < wy.get())
		lb_left.push_back(v);
	      if (wy.get() < wx.get())
		lb_right.push_back(v);
	    }
	  }

	  CRAB_LOG("octagon-join",
		   crab::outs() << "lb_right={";
		   for (vert_id s : lb_right) {
		     crab::outs() << s << ";";
		   } crab::outs() << "}\n";
		   crab::outs() << "lb_left={";
		   for (vert_id s : lb_left) {
		     crab::outs() << s << ";";
		   } crab::outs() << "}\n";
		   crab::outs() << "ub_right={";
		   for (vert_id s : ub_right) {
		     crab::outs() << s << ";";
		   } crab::outs() << "}\n";
		   crab::outs() << "ub_left={";
		   for (vert_id s : ub_left) {
		     crab::outs() << s << ";";
		   } crab::outs() << "}\n";);
	  
	  Wt_min min_op;
	  for (vert_id s : lb_left) {
	    Wt dx_s = gx.edge_val(s, s + 1) / (Wt)2;
	    Wt dy_s = gy.edge_val(s, s + 1) / (Wt)2;
            for (vert_id d : ub_right) { // CASE 1
            /*
                     a                   b
            pos(x) ---> neg(x), neg(y) ---> pos(y)  in G1  // {-x <= a/2, y <=
            b/2} c                   d pos(x) ---> neg(x), neg(y) ---> pos(y) in
            G2  // {-x <= c/2, y <= d/2}

            ===> -x+y <= max(a+b, c+d)  if max(a+b, c+d) < max(a,c) + max(b,d)
            edge from pos(y) to pos(x) or from neg(x) to neg(y)
            */
	      if (s == d)
		continue;
	      join_g.update_edge(s,
				 std::max(dx_s + gx.edge_val(d + 1, d) / (Wt)2,
					  dy_s + gy.edge_val(d + 1, d) / (Wt)2),
				 d, min_op);
	      join_g.update_edge(d + 1,
				 std::max(dx_s + gx.edge_val(d + 1, d) / (Wt)2,
					  dy_s + gy.edge_val(d + 1, d) / (Wt)2),
				 s + 1, min_op);
	    }
	    for (vert_id d : lb_right) { // CASE 2
	    /*
                    a                   b
            pos(x) ---> neg(x), pos(y) ---> neg(y)  in G1 // {-x <= a/2, -y <=
            b/2} c                   d pos(x) ---> neg(x), pos(y) ---> neg(y) in
            G2 // {-x <= c/2, -y <= d/2}

            ===>  -x-y <= max(a+b, c+d) if max(a+b, c+d) < max(a,c) + max(b,d)
            edge from neg(x) to pos(y) or from neg(y) to pos(x)
          */
	      if (s == d)
		continue;
	      join_g.update_edge(d,
				 std::max(dx_s + gx.edge_val(d, d + 1) / (Wt)2,
					  dy_s + gy.edge_val(d, d + 1) / (Wt)2),
				 s + 1, min_op);
	      join_g.update_edge(s,
				 std::max(dx_s + gx.edge_val(d, d + 1) / (Wt)2,
					  dy_s + gy.edge_val(d, d + 1) / (Wt)2),
				 d + 1, min_op);
	    }
      }
	  
	  for (vert_id s : ub_left) {
	    Wt dx_s = gx.edge_val(s + 1, s) / (Wt)2;
	    Wt dy_s = gy.edge_val(s + 1, s) / (Wt)2;
	    for (vert_id d : lb_right) { // CASE 3
            /*
                    a                   b
            neg(x) ---> pos(x), pos(y) ---> neg(y)  in G1 // {x <= a/2, -y <=
            b/2} c                   d neg(x) ---> pos(x), pos(y) ---> neg(y) in
            G2 // {x <= c/2, -y <= d/2}

            ===>  x-y <= max(a+b, c+d) if max(a+b, c+d) < max(a,c) + max(b,d)
            edge from pos(x) to pos(y) or from neg(y) to neg(x)
           */
	      if (s == d)
		continue;
	      join_g.update_edge(d,
				 std::max(dx_s + gx.edge_val(d, d + 1) / (Wt)2,
					  dy_s + gy.edge_val(d, d + 1) / (Wt)2),
				 s, min_op);
	      join_g.update_edge(s + 1,
				 std::max(dx_s + gx.edge_val(d, d + 1) / (Wt)2,
					  dy_s + gy.edge_val(d, d + 1) / (Wt)2),
				 d + 1, min_op);
	    }
	    for (vert_id d : ub_right) { // CASE 4
            /*
                    a                   b
            neg(x) ---> pos(x), neg(y) ---> pos(y)  in G1 // {x <= a/2, y <=
            b/2} c                   d neg(x) ---> pos(x), neg(y) ---> pos(y) in
            G2 // {x <= c/2, y <= d/2}

            ===>  x+y <= max(a+b, c+d) if max(a+b, c+d) < max(a,c) + max(b,d)
            edge from pos(x) to neg(y) or from pos(y) to neg(x)
            */
	      if (s == d)
		continue;
	      join_g.update_edge(s + 1,
				 std::max(dx_s + gx.edge_val(d + 1, d) / (Wt)2,
					  dy_s + gy.edge_val(d + 1, d) / (Wt)2),
				 d, min_op);
	      join_g.update_edge(d + 1,
				 std::max(dx_s + gx.edge_val(d + 1, d) / (Wt)2,
					  dy_s + gy.edge_val(d + 1, d) / (Wt)2),
				 s, min_op);
	    }
	  }
	  // Conjecture: join_g remains closed.
	  // Now garbage collect any unused vertices
	  CRAB_LOG("octagon-join", crab::outs()
		   << "Joined graph (with implicit edges):\n"
		   << join_g << "\n");
	  for (vert_id v : join_g.verts()) {
	    if (v % 2 != 0)
	      continue;
	    if (join_g.succs(v).size() == 0 && join_g.preds(v).size() == 0 &&
		join_g.succs(v + 1).size() == 0 &&
		join_g.preds(v + 1).size() == 0) {
	      join_g.forget(v);
	      join_g.forget(v + 1);
	      if (out_revmap[v]) {
		out_vmap.erase(*(out_revmap[v]));
		out_revmap[v] = boost::none;
		out_revmap[v + 1] = boost::none;
	      }
	    }
	  }
	  split_oct_domain_t res(std::move(out_vmap), std::move(out_revmap),
				 std::move(join_g), std::move(pot_rx),
				 vert_set_t());
	  // join_g.check_adjs();
	  CRAB_LOG("octagon", crab::outs() << "Result join:\n" << res << "\n");
	  CRAB_LOG("octagon-join", crab::outs() << join_g << "\n";);
	  check_potential(res.m_graph, res.m_potential, __LINE__);
	  return res;
	};

      if (need_normalization() && o.need_normalization()) {
	split_oct_domain_t left(*this);
	split_oct_domain_t right(o);
	left.normalize();
	right.normalize();
	return join_op(left, right);
      } else if (need_normalization()) {
	split_oct_domain_t left(*this);
	const split_oct_domain_t &right = o;
	left.normalize();
	return join_op(left, right);
      } else if (o.need_normalization()) {
	const split_oct_domain_t &left = *this;
	split_oct_domain_t right(o);
	right.normalize();
	return join_op(left, right);
      } else {
	return join_op(*this, o);
      }
    }
  }

  // widening
  split_oct_domain_t operator||(const split_oct_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else {
      CRAB_LOG("octagon",
               crab::outs() << "Before widening:\n";
               crab::outs() << "DBM 1\n"
                            << *this << "\n";
               crab::outs() << "DBM 2\n"
                            << o << "\n";);

      auto widen_op =
	[this](const split_oct_domain_t &left, const split_oct_domain_t &right) -> split_oct_domain_t {
	  // Figure out the common renaming
	  std::vector<vert_id> perm_x, perm_y;
	  vert_map_t out_vmap;
	  rev_map_t out_revmap;
	  std::vector<Wt> widen_pot;
	  vert_set_t widen_unstable(left.m_unstable);
	  
	  assert(left.m_potential.size() > 0);
	  for (auto p : left.m_vert_map) {
	    auto it = right.m_vert_map.find(p.first);
	    // Variable exists in both
	    if (it != right.m_vert_map.end()) {
	      out_vmap.insert(vmap_elt_t(p.first, {perm_x.size(), perm_x.size() + 1}));
	      out_revmap.push_back(p.first);
	      out_revmap.push_back(p.first);
	      widen_pot.push_back(left.m_potential[p.second.first]);
	      widen_pot.push_back(left.m_potential[p.second.second]);
	      perm_x.push_back(p.second.first);
	      perm_x.push_back(p.second.second);
	      perm_y.push_back((*it).second.first);
	      perm_y.push_back((*it).second.second);
	    }
	  }
	  // Build the permuted view of x and y.
	  assert(left.m_graph.size() > 0);
	  GrPerm gx(perm_x, left.m_graph);
	  assert(right.m_graph.size() > 0);
	  GrPerm gy(perm_y, right.m_graph);
	  
	  CRAB_LOG("octagon-widening", crab::outs()
                                       << "== After permutations == \n";
               crab::outs() << "DBM 1:\n"; gx.write(crab::outs());
               crab::outs() << "\n"; crab::outs() << "DBM 2:\n";
               gy.write(crab::outs()); crab::outs() << "\n";
               crab::outs() << "rev_map={"; for (unsigned i = 0,
                                                 e = out_revmap.size();
                                                 i != e; i++) {
                 if (out_revmap[i])
                   crab::outs() << *(out_revmap[i]) << "(" << i << ");";
               } crab::outs() << "}\n";);

	  // Now perform the widening
	  graph_t widen_g(split_widen(gx, gy, widen_unstable));
	  
	  CRAB_LOG("octagon-widening", crab::outs()
		   << "Unstable list after point-wise widening{";
          for (vert_id v : widen_unstable) {
		 
            crab::outs() << *out_revmap[v] << ((v % 2 == 0) ? "+" : "-") << ";";
          } crab::outs()
          << "}\n");

	  split_oct_domain_t res(std::move(out_vmap), std::move(out_revmap),
				 std::move(widen_g), std::move(widen_pot),
				 std::move(widen_unstable));
	  
	  CRAB_LOG("octagon-widening", crab::outs() << "Result widening:\n"
		   << res.m_graph << "\n";);
	  CRAB_LOG("octagon",
		   crab::outs() << "Result widening:\n" << res << "\n";);
	  return res;
	};
      
      // Do not normalize left operand
      const split_oct_domain_t &left = *this;
      if (o.need_normalization()) {
	split_oct_domain_t right(o);
	right.normalize();
	return widen_op(left, right);
      } else {
	return widen_op(left, o);
      }
      
    }
  }

  // meet
  split_oct_domain_t operator&(const split_oct_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top())
      return *this;
    else if (is_top() || o.is_bottom()) {
      return o;
    } else {

      CRAB_LOG("octagon", crab::outs() << "Before meet:\n"
                                       << "DBM 1\n"
                                       << *this << "\n"
                                       << "DBM 2\n"
                                       << o << "\n");

      auto meet_op =
	[this](const split_oct_domain_t &left, const split_oct_domain_t &right) -> split_oct_domain_t {
	  vert_map_t meet_verts;
	  rev_map_t meet_rev;
	  std::vector<vert_id> perm_x, perm_y;
	  std::vector<Wt> meet_pi;
	  
	  for (auto p : left.m_vert_map) {
	    vert_id vv = perm_x.size();
	    meet_verts.insert(vmap_elt_t(p.first, {vv, vv + 1}));
	    meet_rev.push_back(p.first);
	    meet_rev.push_back(p.first);
	    
	    perm_x.push_back(p.second.first);
	    perm_x.push_back(p.second.second);
	    perm_y.push_back(-1);
	    perm_y.push_back(-1);
	    meet_pi.push_back(left.m_potential[p.second.first]);
	    meet_pi.push_back(left.m_potential[p.second.second]);
	  }
	  
	  // Add missing mappings from the right operand.
	  for (auto p : right.m_vert_map) {
	    auto it = meet_verts.find(p.first);
	    
	    if (it == meet_verts.end()) {
	      vert_id vv = perm_y.size();
	      meet_rev.push_back(p.first);
	      meet_rev.push_back(p.first);
	      
	      perm_y.push_back(p.second.first);
	      perm_y.push_back(p.second.second);
	      perm_x.push_back(-1);
	      perm_x.push_back(-1);
	      meet_pi.push_back(right.m_potential[p.second.first]);
	      meet_pi.push_back(right.m_potential[p.second.second]);
	      meet_verts.insert(vmap_elt_t(p.first, {vv, vv + 1}));
	    } else {
	      perm_y[(*it).second.first] = p.second.first;
	      perm_y[(*it).second.second] = p.second.second;
	    }
	  }
	  
	  // Build the permuted view of x and y.
	  assert(left.m_graph.size() > 0);
	  GrPerm gx(perm_x, left.m_graph);
	  assert(right.m_graph.size() > 0);
	  GrPerm gy(perm_y, right.m_graph);
	  
	  // Compute the syntactic meet of the permuted graphs.
	  bool is_closed;
	  graph_t meet_g(GrOps::meet(gx, gy, is_closed));

	  // Compute updated potentials on the zero-enriched graph
	  // vector<Wt> meet_pi(meet_g.size());
	  // We've warm-started pi with the operand potentials
	  if (!GrOps::select_potentials(meet_g, meet_pi)) {
	    // Potentials cannot be selected -- state is infeasible.
	    return make_bottom();
	  }

	  if (!is_closed) {
	    // JN: this code needs to be tested
	    edge_vector delta;
	    split_octagons_impl::SplitOctGraph<graph_t> meet_g_oct(meet_g);
	    
	    if (crab_domain_params_man::get().oct_chrome_dijkstra()) {
	      GrOps::close_after_meet(meet_g_oct, meet_pi, gx, gy, delta);
	    } else {
	      GrOps::close_johnson(meet_g_oct, meet_pi, delta);
	    }
	    
	    // JN: we should be fine calling GrOps::apply_delta
	    update_delta(meet_g, delta);

	    // Wt_min min_op;
	    // for(auto e : delta) {
	    //   if (e.first.first % 2 == 2) {
	    //     if(meet_g.elem(e.first.first+1, e.first.first))
	    // 	meet_g.update_edge(e.first.first+1,
	    // 			   meet_g.edge_val(e.first.first+1,
	    // e.first.first)
	    //  			   (Wt) 2*e.second, 			   e.first.second,
	    //  min_op);
	    //     if(meet_g.elem(e.first.second, e.first.first+1))
	    // 	meet_g.update_edge(e.first.first,
	    // 			   meet_g.edge_val(e.first.second,
	    // e.first.first+1)
	    //  			   (Wt) 2*e.second, 			   e.first.first+1,
	    //  min_op);
	    //   } else if (e.first.first % 2 == 2) {
	    //     if(meet_g.elem(e.first.first, e.first.first+1))
	    // 	meet_g.update_edge(e.first.first,
	    // 			   meet_g.edge_val(e.first.first,
	    // e.first.first+1)
	    //  			   (Wt) 2*e.second, 			   e.first.second,
	    //  min_op);
	    //     if(meet_g.elem(e.first.second, e.first.first+1))
	    // 	meet_g.update_edge(e.first.first+1,
	    // 			   meet_g.edge_val(e.first.second+1,
	    // e.first.first)
	    //  			   (Wt) 2*e.second, 			   e.first.first,
	    //  min_op);
	    //   }
	    // }
	  }

	  check_potential(meet_g, meet_pi, __LINE__);
	  split_oct_domain_t res(std::move(meet_verts), std::move(meet_rev),
				 std::move(meet_g), std::move(meet_pi),
				 vert_set_t());
	  
	  CRAB_LOG("octagon", crab::outs() << "Result meet:\n" << res << "\n");
	  return res;
	};

      if (need_normalization() && o.need_normalization()) {
	split_oct_domain_t left(*this);
	split_oct_domain_t right(o);
	left.normalize();
	right.normalize();
	return meet_op(left, right);
      } else if (need_normalization()) {
	split_oct_domain_t left(*this);
	const split_oct_domain_t &right = o;
	left.normalize();
	return meet_op(left, right);
      } else if (o.need_normalization()) {
	const split_oct_domain_t &left = *this;
	split_oct_domain_t right(o);
	right.normalize();
	return meet_op(left, right);
      } else {
	return meet_op(*this, o);
      }

      
    }
  }

  void operator&=(const split_oct_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      // do nothing
    } else if (is_top() || o.is_bottom()) {
      *this = o;
    } else {

      CRAB_LOG("octagon", crab::outs() << "Before meet:\n"
                                       << "DBM 1\n"
                                       << *this << "\n"
                                       << "DBM 2\n"
                                       << o << "\n");

      auto meet_op = [this](split_oct_domain_t &left, const split_oct_domain_t &right)  {
	vert_map_t meet_verts;
	rev_map_t meet_rev;
	std::vector<vert_id> perm_x, perm_y;
	std::vector<Wt> meet_pi;

	for (auto p : left.m_vert_map) {
	  vert_id vv = perm_x.size();
	  meet_verts.insert(vmap_elt_t(p.first, {vv, vv + 1}));
	  meet_rev.push_back(p.first);
	  meet_rev.push_back(p.first);
	  
	  perm_x.push_back(p.second.first);
	  perm_x.push_back(p.second.second);
	  perm_y.push_back(-1);
	  perm_y.push_back(-1);
	  meet_pi.push_back(left.m_potential[p.second.first]);
	  meet_pi.push_back(left.m_potential[p.second.second]);
	}
	
	// Add missing mappings from the right operand.
	for (auto p : right.m_vert_map) {
	  auto it = meet_verts.find(p.first);
	  
	  if (it == meet_verts.end()) {
	    vert_id vv = perm_y.size();
	    meet_rev.push_back(p.first);
	    meet_rev.push_back(p.first);
	    
	    perm_y.push_back(p.second.first);
	    perm_y.push_back(p.second.second);
	    perm_x.push_back(-1);
	    perm_x.push_back(-1);
	    meet_pi.push_back(right.m_potential[p.second.first]);
	    meet_pi.push_back(right.m_potential[p.second.second]);
	    meet_verts.insert(vmap_elt_t(p.first, {vv, vv + 1}));
	  } else {
	    perm_y[(*it).second.first] = p.second.first;
	    perm_y[(*it).second.second] = p.second.second;
	  }
	}
	
	// Build the permuted view of x and y.
	assert(left.m_graph.size() > 0);
	GrPerm gx(perm_x, left.m_graph);
	assert(right.m_graph.size() > 0);
	GrPerm gy(perm_y, right.m_graph);
	
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
	  split_octagons_impl::SplitOctGraph<graph_t> meet_g_oct(meet_g);
	  if (crab_domain_params_man::get().oct_chrome_dijkstra()) {
	    GrOps::close_after_meet(meet_g_oct, meet_pi, gx, gy, delta);
	  } else {
	    GrOps::close_johnson(meet_g_oct, meet_pi, delta);
	  }
	  update_delta(meet_g, delta);
	}
	
	check_potential(meet_g, meet_pi, __LINE__);
	
	left.m_vert_map = std::move(meet_verts);
	left.m_rev_vert_map = std::move(meet_rev);
	left.m_graph = std::move(meet_g);
	left.m_potential = std::move(meet_pi);
	left.m_unstable.clear();
	left.m_is_bottom = false;
	
	CRAB_LOG("octagon", crab::outs() << "Result meet:\n" << left << "\n");
      };

      split_oct_domain_t &left = *this;
      left.normalize();

      if (o.need_normalization()) {
	split_oct_domain_t right(o);
	right.normalize();
	meet_op(left, right);
      } else {
	meet_op(left, o);
      }
     
    }
  }
  
  split_oct_domain_t widening_thresholds(
      const split_oct_domain_t &o,
      const thresholds<number_t> &ts) const override {
    return (*this || o);
  }

  split_oct_domain_t operator&&(const split_oct_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_bottom())
      return *this;
    else if (is_top())
      return o;
    else {
      CRAB_LOG("octagon", crab::outs() << "Before narrowing:\n"
                                       << "DBM 1\n"
                                       << *this << "\n"
                                       << "DBM 2\n"
                                       << o << "\n");

      // TODO: Implement properly
#if 1
      // Narrowing implemented as meet might not terminate.
      // Make sure that there is always a maximum bound for narrowing iterations.
      return *this & o;
#else
      // Narrowing as a no-op: sound and it will terminate
      if (need_normalization()) {
	split_oct_domain_t res(*this);
	res.normalize();      
	CRAB_LOG("octagon", crab::outs() << "Result narrowing:\n" << res << "\n");
	return res;
      } else {
	CRAB_LOG("octagon", crab::outs() << "Result narrowing:\n" << *this << "\n");
	return *this;
      }
#endif      
    }
  }

  void minimize() override {}

  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom()) {
      return;
    }

    normalize();
    auto it = m_vert_map.find(v);
    if (it != m_vert_map.end()) {
      m_graph.forget((*it).second.first);
      m_graph.forget((*it).second.second);
      m_rev_vert_map[(*it).second.first] = boost::none;
      m_rev_vert_map[(*it).second.second] = boost::none;
      m_vert_map.erase(v);
    }
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
    } else if (cst.is_disequation()) {
      // We handle here the case x !=y by converting the disequation
      // into a strict inequality if possible.
      linear_constraint_system_t csts;
      constraint_simp_domain_traits<split_oct_domain_t>::lower_disequality(*this, cst, csts);
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

    CRAB_LOG("octagon", crab::outs() << "---" << cst << "\n" << *this << "\n");
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    if (is_bottom())
      return;
    for (auto cst : csts) {
      operator+=(cst);
    }
  }

  DEFAULT_ENTAILS(split_oct_domain_t)
  
  interval_t operator[](const variable_t &x) override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");
    // Needed for accuracy
    normalize();

    if (this->is_bottom()) {
      return interval_t::bottom();
    } else {
      return get_interval(m_vert_map, m_graph, x);
    }
  }

  interval_t at(const variable_t &x) const override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");
    
    return (is_bottom() ? interval_t::bottom() :
	    get_interval(m_vert_map, m_graph, x));
  }
  
  void normalize(void) override {
    crab::CrabStats::count(domain_name() + ".count.normalize");
    crab::ScopedCrabStats __st__(domain_name() + ".normalize");

    /// The DBM is always in split normal form except after widening

    if (m_unstable.empty()) {
      return;
    }

    CRAB_LOG(
        "octagon-normalize", crab::outs()
                                 << "Normalization...\nUnstable list {";
        for (vert_id v
             : m_unstable) {
          crab::outs() << *m_rev_vert_map[v] << ((v % 2 == 0) ? "+" : "-")
                       << ";";
        } crab::outs()
        << "}\n";);

    edge_vector delta;
    split_octagons_impl::SplitOctGraph<graph_t> g_oct(m_graph);
    if (crab_domain_params_man::get().oct_widen_restabilize()) {
      GrOps::close_after_widen(g_oct, m_potential, vert_set_wrap_t(m_unstable),
                               delta);
      m_unstable.clear();
    } else {
      GrOps::close_johnson(g_oct, m_potential, delta);
    }

    CRAB_LOG(
        "octagon-normalize", crab::outs() << "Edges after close_after_widen {";
        for (auto edge
             : delta) {
          vert_id src = edge.first.first;
          vert_id dst = edge.first.second;
          Wt w = edge.second;
          crab::outs() << *m_rev_vert_map[src] << ((src % 2 == 0) ? "+" : "-")
                       << "-->" << *m_rev_vert_map[dst]
                       << ((dst % 2 == 0) ? "+" : "-") << " w=" << w << ";";
        } crab::outs()
        << "}\n";);

    // JN: we should be fine calling GrOps::apply_delta
    update_delta(m_graph, delta);

    // Recover updated bounds
    for (auto edge : delta) {
      vert_id src = edge.first.first;
      vert_id dst = edge.first.second;
      Wt w = edge.second;
      edge_vector bounds_delta;
      update_bounds(m_graph, src, dst, w, bounds_delta);
      GrOps::apply_delta(m_graph, bounds_delta);
    }
  }

  // set a variable to an interval
  void set(const variable_t &x, interval_t intv) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");
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
      Wt ub = ntow::convert(number_t(2) * (*(intv.ub().number())), overflow);
      if (overflow) {
        return;
      }
      m_potential[v] = ub / (Wt)2;
      m_potential[v + 1] = -ub / (Wt)2;
      m_graph.set_edge(v + 1, ub, v);
    }
    if (intv.lb().is_finite()) {
      Wt lb = ntow::convert(-number_t(2) * (*(intv.lb().number())), overflow);
      if (overflow) {
        return;
      }
      m_potential[v] = lb / (Wt)-2;
      m_potential[v + 1] = lb / (Wt)2;
      m_graph.set_edge(v, lb, v + 1);
    }
  }

  // assign an exact expression to a variable
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (is_bottom())
      return;

    CRAB_LOG("octagon", crab::outs() << "--- " << x << ":=" << e << "\n";);

    normalize();
    check_potential(m_graph, m_potential, __LINE__);

    if (e.is_constant()) {
      set(x, e.constant());
    } else {
      interval_t x_int = eval_interval(e);

      boost::optional<Wt> lb_w, ub_w;
      bool overflow;
      if (x_int.lb().is_finite()) {
        lb_w = ntow::convert(-number_t(2) * (*(x_int.lb().number())), overflow);
        if (overflow) {
          operator-=(x);
          CRAB_LOG("octagon", crab::outs() << "LB overflow=" << x_int << "\n";);
          CRAB_LOG("octagon", crab::outs() << "---" << x << ":=" << e << "\n"
                                           << *this << "\n");
          return;
        }
      }
      if (x_int.ub().is_finite()) {
        ub_w = ntow::convert(number_t(2) * (*(x_int.ub().number())), overflow);
        if (overflow) {
          operator-=(x);
          CRAB_LOG("octagon", crab::outs() << "UB overflow=" << x_int << "\n";);
          CRAB_LOG("octagon", crab::outs() << "---" << x << ":=" << e << "\n"
                                           << *this << "\n");
          return;
        }
      }

      if (boost::optional<number_t> x_n = x_int.singleton()) {
        set(x, *x_n);
      } else {
        std::vector<std::pair<variable_t, Wt>> diffs_lb, diffs_ub;
        std::vector<std::pair<variable_t, Wt>> octs_lb, octs_ub;
        CRAB_LOG(
            "octagon-assign", crab::outs() << *this << "\n";
            crab::outs() << m_graph << "\n"; for (auto kv
                                                  : m_vert_map) {
              crab::outs() << "\t" << kv.first << "=(" << kv.second.first << ","
                           << kv.second.second << ")\n";
            } crab::outs() << "Octagon constraints from "
                           << x << ":=" << e << "={";);
        oct_csts_of_assign(x, e, diffs_lb, octs_lb, diffs_ub, octs_ub);
        CRAB_LOG("octagon-assign", crab::outs() << "}\n";);

        if (diffs_lb.size() > 0 || diffs_ub.size() > 0 || octs_lb.size() > 0 ||
            octs_ub.size() > 0) {

          /********************************************************/
          /**     Assignment as a sequence of edge additions     **/
          /********************************************************/

          vert_id v = m_graph.new_vertex();
          assert(v <= m_rev_vert_map.size());
          vert_id w = m_graph.new_vertex();

          if (w < v) {
            std::swap(w, v);
          }

          if (v == m_rev_vert_map.size()) {
            m_rev_vert_map.push_back(x);
            m_potential.push_back(Wt(0));
          } else {
            m_potential[v] = Wt(0);
            m_rev_vert_map[v] = x;
          }
          if (w == m_rev_vert_map.size()) {
            m_rev_vert_map.push_back(x);
            m_potential.push_back(Wt(0));
          } else {
            m_potential[w] = Wt(0);
            m_rev_vert_map[w] = x;
          }

          edge_vector cst_edges;
          for (auto diff : diffs_lb) {
            variable_t y = diff.first;
            Wt k = diff.second;
            // x - y >= k ---> pos(Vy) - pos(Vx) <= -k --> edge(pos(Vx), -k,
            // pos(Vy))
            cst_edges.push_back({{v, get_vert(y)}, -k});
            // to preserve coherence
            cst_edges.push_back({{get_vert(y) + 1, v + 1}, -k});
          }

          for (auto diff : diffs_ub) {
            variable_t y = diff.first;
            Wt k = diff.second;
            // x - y <= k ---> pos(Vx) - pos(Vy) <= k --> edge(pos(Vy), k,
            // pos(Vx))
            cst_edges.push_back({{get_vert(y), v}, k});
            // to preserve coherence
            cst_edges.push_back({{v + 1, get_vert(y) + 1}, k});
          }

          for (auto oct : octs_lb) {
            variable_t y = oct.first;
            Wt k = oct.second;
            // x + y >= k --> neg(Vx) - pos(Vy) <= -k --> edge(pos(Vy), -k,
            // neg(Vx))
            cst_edges.push_back({{get_vert(y), w}, -k});
            // to preserve coherence
            cst_edges.push_back({{w - 1, get_vert(y) + 1}, -k});
          }

          for (auto oct : octs_ub) {
            variable_t y = oct.first;
            Wt k = oct.second;
            // x + y <= k --> pos(Vx) - neg(Vy) <= k --> edge(neg(Vy), k,
            // pos(Vx))
            cst_edges.push_back({{get_vert(y) + 1, v}, k});
            // to preserve coherence
            cst_edges.push_back({{v + 1, get_vert(y)}, k});
          }

          Wt_min min_op;
          for (auto diff : cst_edges) {
            vert_id src = diff.first.first;
            vert_id dest = diff.first.second;
            m_graph.update_edge(src, diff.second, dest, min_op);
            CRAB_LOG("octagon-assign",
                     crab::outs() << "\tAdded edge " << src << "-->" << dest
                                  << " w=" << diff.second << "\n";);
            if (!repair_potential(src, dest)) {
              CRAB_ERROR("bottom when split_oct::assign");
              // set_to_bottom();
            }
            close_over_edge(src, dest);
          }

          check_potential(m_graph, m_potential, __LINE__);

          if (lb_w) {
            CRAB_LOG("octagon-assign",
                     crab::outs() << "\tAdded edge (interval) " << v << "-->"
                                  << w << " w=" << *lb_w << "\n";);
            m_graph.update_edge(v, *lb_w, w, min_op);
            if (!update_bounds_lb(v, *lb_w)) {
              CRAB_ERROR("bottom when split_oct::update_bounds_lb");
            }
          }

          if (ub_w) {
            CRAB_LOG("octagon-assign",
                     crab::outs() << "\tAdded edge (interval) " << w << "-->"
                                  << v << " w=" << *ub_w << "\n";);
            m_graph.update_edge(w, *ub_w, v, min_op);
            if (!update_bounds_ub(v, *ub_w)) {
              CRAB_ERROR("bottom when split_oct::update_bounds_ub");
            }
          }
          check_potential(m_graph, m_potential, __LINE__);
          // Clear the old x vertex
          operator-=(x);
          assert(!is_bottom());
          m_vert_map.insert(vmap_elt_t(x, {v, w}));
        }
      }
    }

    check_potential(m_graph, m_potential, __LINE__);
    integer_tightening();
    check_potential(m_graph, m_potential, __LINE__);
    CRAB_LOG("octagon", crab::outs() << "---" << x << ":=" << e << "\n"
                                     << *this << "\n";);
  }

  void forget(const variable_vector_t &vars) override {
    if (is_bottom() || is_top())
      return;

    for (auto v : vars) {
      auto it = m_vert_map.find(v);
      if (it != m_vert_map.end()) {
        operator-=(v);
      }
    }
  }

  void project(const variable_vector_t &vars) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");
    if (is_bottom() || is_top())
      return;

    normalize();
    std::vector<bool> save(m_rev_vert_map.size(), false);
    for (auto x : vars) {
      auto it = m_vert_map.find(x);
      if (it != m_vert_map.end()) {
        save[(*it).second.first] = true;
        save[(*it).second.second] = true;
      }
    }
    for (vert_id v = 0; v < m_rev_vert_map.size(); v++) {
      if (!save[v] && m_rev_vert_map[v])
        operator-=((*m_rev_vert_map[v]));
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;

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
      // TODO: we can improve here by linearizing one of the
      // operands and then calling assign.
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
      // case OP_UREM:
      set(x, get_interval(y).URem(get_interval(z)));
      break;
    }
    CRAB_LOG("octagon", crab::outs()
                            << "---" << x << ":=" << y << op << z << "\n"
                            << *this << "\n");
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;

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
    CRAB_LOG("octagon", crab::outs()
                            << "---" << x << ":=" << y << op << k << "\n"
                            << *this << "\n");
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    int_cast_domain_traits<split_oct_domain_t>::apply(*this, op, dst, src);    
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;

    normalize();
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
      //  case OP_ASHR: 
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    // Convert to intervals and perform the operation
    if (is_bottom())
      return;

    normalize();
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

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const split_oct_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");
    if (is_bottom())
      return;
    normalize();
    crab::domains::BackwardAssignOps<split_oct_domain_t>::assign(*this, x, e,
                                                                 inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const split_oct_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    if (is_bottom())
      return;
    normalize();
    crab::domains::BackwardAssignOps<split_oct_domain_t>::apply(*this, op, x, y,
                                                                z, inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const split_oct_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    if (is_bottom())
      return;
    normalize();
    crab::domains::BackwardAssignOps<split_oct_domain_t>::apply(*this, op, x, y,
                                                                z, inv);
  }

  void expand(const variable_t &x, const variable_t &y) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top())
      return;

    normalize();

    CRAB_LOG(
        "octagon", crab::outs()
                       << "Before expand " << x << " into " << y << ":\n"
                       << *this << "\n"
                       << m_graph << "\n"
                       << "rev_map ={";
        for (unsigned i = 0, e = m_rev_vert_map.size(); i != e; i++) {
          if (m_rev_vert_map[i])
            crab::outs() << *(m_rev_vert_map[i]) << "(" << i << ");";
        } crab::outs()
        << "}\n";);

    auto it = m_vert_map.find(variable_t(y));
    if (it != m_vert_map.end()) {
      CRAB_ERROR("split_dbm expand operation failed because y already exists");
    }

    vert_id ii = get_vert(x);
    vert_id jj = get_vert(y);

    edge_vector delta;
    for (auto edge : m_graph.e_preds(ii)) {
      delta.push_back(
          {{edge.vert == ii + 1 ? jj + 1 : edge.vert, jj}, edge.val});
    }

    for (auto edge : m_graph.e_succs(ii)) {
      delta.push_back(
          {{jj, edge.vert == ii + 1 ? jj + 1 : edge.vert}, edge.val});
    }

    for (auto edge : m_graph.e_preds(ii + 1)) {
      delta.push_back({{edge.vert == ii ? jj : edge.vert, jj + 1}, edge.val});
    }

    for (auto edge : m_graph.e_succs(ii + 1)) {
      delta.push_back({{jj + 1, edge.vert == ii ? jj : edge.vert}, edge.val});
    }
    GrOps::apply_delta(m_graph, delta);

    m_potential[jj] = m_potential[ii];
    m_potential[jj + 1] = m_potential[ii + 1];

    CRAB_LOG(
        "octagon", crab::outs()
                       << "After expand " << x << " into " << y << ":\n"
                       << *this << "\n"
                       << m_graph << "\n"
                       << "rev_map ={";
        for (unsigned i = 0, e = m_rev_vert_map.size(); i != e; i++) {
          if (m_rev_vert_map[i])
            crab::outs() << *(m_rev_vert_map[i]) << "(" << i << ");";
        } crab::outs()
        << "}\n";);
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_top() || is_bottom())
      return;

    normalize();

    CRAB_LOG("octagon", crab::outs() << "Replacing {"; for (auto v
                                                            : from) crab::outs()
                                                       << v << ";";
             crab::outs() << "} with "; for (auto v
                                             : to) crab::outs()
                                        << v << ";";
             crab::outs() << "}:\n"; crab::outs() << *this << "\n";);

    vert_map_t new_vert_map;
    for (auto kv : m_vert_map) {
      ptrdiff_t pos = std::distance(
          from.begin(), std::find(from.begin(), from.end(), kv.first));
      if (pos < from.size()) {
        variable_t new_v(to[pos]);
        new_vert_map.insert(vmap_elt_t(new_v, kv.second));
        m_rev_vert_map[kv.second.first] = new_v;
        m_rev_vert_map[kv.second.second] = new_v;
      } else {
        new_vert_map.insert(kv);
      }
    }
    std::swap(m_vert_map, new_vert_map);
    CRAB_LOG("octagon", crab::outs() << "RESULT=" << *this << "\n");
  }

  /// split_oct_domain_t implements only standard abstract operations
  /// of a numerical domain so it is intended to be used as a leaf
  /// domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(split_oct_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(split_oct_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(split_oct_domain_t)
  DEFAULT_SELECT(split_oct_domain_t)
  DEFAULT_WEAK_ASSIGN(split_oct_domain_t)  

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const split_oct_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

#if 1
    /// The normalization that happens in linear_constraint_t
    /// removes redundant constraints (see below).
    auto csts = to_linear_constraint_system().normalize();
    o << csts;
#else
    /// This code can print the same constraint twice. E.g., the
    /// constraint x - y <= k can be printed if we have edges from
    /// pos(y) to pos(x) and from neg(x) to neg(y).
    wt_ref_t w;
    if (is_bottom()) {
      o << "_|_";
      return;
    } else if (is_top()) {
      o << "{}";
      return;
    } else {
      split_oct_domain_t tmp(*this);
      tmp.normalize();
      bool first = true;
      o << "{";
      for (vert_id v : tmp.m_graph.verts()) {
        if (v % 2 != 0)
          continue;
        if (!tmp.m_rev_vert_map[v])
          continue;
        if (!tmp.m_graph.elem(v, v + 1) && !tmp.m_graph.elem(v + 1, v))
          continue;
        interval_t v_out =
            interval_t(tmp.m_graph.elem(v, v + 1)
                           ? -number_t(tmp.m_graph.edge_val(v, v + 1)) / 2
                           : bound_t::minus_infinity(),
                       tmp.m_graph.elem(v + 1, v)
                           ? number_t(tmp.m_graph.edge_val(v + 1, v)) / 2
                           : bound_t::plus_infinity());

        if (first)
          first = false;
        else
          o << ", ";
        o << *(tmp.m_rev_vert_map[v]) << " -> " << v_out;
      }
      split_octagons_impl::SplitOctGraph<graph_t> g_oct(tmp.m_graph);
      for (vert_id s : tmp.m_graph.verts()) {
        if (!tmp.m_rev_vert_map[s])
          continue;
        variable_t vs = *(tmp.m_rev_vert_map[s]);
        for (vert_id d : g_oct.succs(s)) {
          if (!tmp.m_rev_vert_map[d])
            continue;
          variable_t vd = *(tmp.m_rev_vert_map[d]);
          if (first)
            first = false;
          else
            o << ", ";

          if (s % 2 == 0 && d % 2 == 0) {
            o << vd << "-" << vs << "<=" << g_oct.edge_val(s, d);
          } else if (s % 2 != 0 && d % 2 == 0) {
            o << vd << "+" << vs << "<=" << g_oct.edge_val(s, d);
          } else if (s % 2 == 0 && d % 2 != 0) {
            o << "-" << vd << "-" << vs << "<=" << g_oct.edge_val(s, d);
          } else if (s % 2 != 0 && d % 2 != 0) {
            o << "-" << vd << "+" << vs << "<=" << g_oct.edge_val(s, d);
          }
        }
      }
      o << '}';
    }
#endif
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;
    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
      return csts;
    }

    split_oct_domain_t tmp(*this);
    tmp.normalize();

    split_octagons_impl::SplitOctGraph<graph_t> g_oct(tmp.m_graph);
    for (vert_id v : g_oct.verts()) {
      if (v % 2 != 0)
        continue;
      if (!tmp.m_rev_vert_map[v])
        continue;
      if (tmp.m_graph.elem(v, v + 1)) {
        csts += linear_constraint_t(*tmp.m_rev_vert_map[v] >=
                                    -(tmp.m_graph.edge_val(v, v + 1) / (Wt)2));
      }
      if (tmp.m_graph.elem(v + 1, v)) {
        csts += linear_constraint_t(*tmp.m_rev_vert_map[v] <=
                                    (tmp.m_graph.edge_val(v + 1, v) / (Wt)2));
      }
    }

    for (vert_id s : g_oct.verts()) {
      if (!tmp.m_rev_vert_map[s])
        continue;
      variable_t vs = *tmp.m_rev_vert_map[s];
      linear_expression_t s_exp(vs);
      for (vert_id d : g_oct.succs(s)) {
        if (!tmp.m_rev_vert_map[d])
          continue;
        variable_t vd = *tmp.m_rev_vert_map[d];
        linear_expression_t d_exp(vd);
        auto w = g_oct.edge_val(s, d);

        if (s % 2 == 0 && d % 2 == 0) {
          csts += linear_constraint_t(d_exp - s_exp <= w);
        } else if (s % 2 != 0 && d % 2 == 0) {
          csts += linear_constraint_t(d_exp + s_exp <= w);
        } else if (s % 2 == 0 && d % 2 != 0) {
          csts += linear_constraint_t(-d_exp - s_exp <= w);
        } else if (s % 2 != 0 && d % 2 != 0) {
          csts += linear_constraint_t(-d_exp + s_exp <= w);
        }
      }
    }
    return csts;
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

  std::string domain_name() const override { return "SplitOctagons"; }

}; // end class split_oct_domain

template <typename Number, typename VariableName, typename Params>
struct abstract_domain_traits<split_oct_domain<Number, VariableName, Params>> {
  using number_t = Number;
  using varname_t = VariableName;
};

} // end namespace domains
} // end namespace crab

#pragma GCC diagnostic pop
