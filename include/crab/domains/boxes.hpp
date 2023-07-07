#pragma once

/**************************************************************************
 * A wrapper for a disjunctive domain of intervals from "Boxes: A
 * Symbolic Abstract Domain of Boxes" by A. Gurfinkel and S. Chaki
 * published in SAS'10.
 **************************************************************************/
#include <crab/config.h>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#ifndef HAVE_LDD
/*
 * Dummy implementation if ldd not found
 */
#include <crab/domains/dummy_abstract_domain.hpp>
namespace crab {
namespace domains {
template <typename N, typename V>
class boxes_domain final
    : public dummy_abstract_domain<boxes_domain<N, V>> {
public:
  std::string not_implemented_msg() const override {
    return "No LDD. Run cmake with -DCRAB_USE_LDD=ON";
  }
};
} // namespace domains
} // namespace crab
#else
/*
 *  Real implementation starts here
 */

#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/ldd/ldd.hpp>
#include <crab/domains/ldd/ldd_print.hpp>

#include <algorithm>
#include <boost/bimap.hpp>
#include <boost/optional.hpp>

namespace crab {

namespace domains {

using namespace crab::domains::ldd;

/*
 * The wrapper has two global datastructures:
 * 1) a ldd manager and
 * 2) a map from VariableName to ldd dimension.
 *
 * FIXME: since the ldd manager is shared we need to fix a
 * single size for all ldds. Since ldds are sparse we can fix a
 * size big enough for our programs.
 *
 * FIXME: Ldd_TermReplace seems to leak memory sometimes.
 */
template <typename Number, typename VariableName>
class boxes_domain final
    : public abstract_domain_api<boxes_domain<Number, VariableName>> {

  using interval_domain_t = ikos::interval_domain<Number, VariableName>;
  using boxes_domain_t = boxes_domain<Number, VariableName>;
  using abstract_domain_t = abstract_domain_api<boxes_domain_t>;

public:
  using number_t = Number;
  using varname_t = VariableName;
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
  
private:
  using kind_t = typename linear_constraint_t::kind_t;

  // --- map from crab variable index to ldd term index
  using var_map_t = boost::bimap<variable_t, int>;
  using binding_t = typename var_map_t::value_type;

  LddNodePtr m_ldd;
  static LddManager *s_ldd_man;
  static var_map_t s_var_map;

  // -- bool reasoning is mostly based on disjunctions so for
  //    efficiency we might want to disable it if precision
  //    gains do not pay off.
  const bool m_bool_reasoning = true;

  static LddManager *get_ldd_man() {
    if (!s_ldd_man) {
      DdManager *cudd = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
      unsigned ldd_size = crab_domain_params_man::get().boxes_ldd_size();
      bool dynamic_reordering = crab_domain_params_man::get().boxes_dynamic_reordering();
      theory_t *theory = ldd::create_box_theory<number_t>(ldd_size);
      CRAB_LOG("boxes", crab::outs() << "Created a ldd of size "
                                     << ldd_size << "\n";);
      s_ldd_man = Ldd_Init(cudd, theory);
      if (dynamic_reordering) {
	Cudd_AutodynEnable(cudd, CUDD_REORDER_GROUP_SIFT);
      }
      Ldd_SanityCheck(get_ldd_man());
    }
    return s_ldd_man;
  }

  theory_t *get_theory() { return Ldd_GetTheory(get_ldd_man()); }

  const theory_t *get_theory() const { return Ldd_GetTheory(get_ldd_man()); }

  static int num_of_vars() { return s_var_map.left.size(); }

  int get_var_dim(variable_t v) const {
    auto it = s_var_map.left.find(v);
    if (it != s_var_map.left.end()) {
      return it->second;
    } else {
      // Reserved dim 0 for SPECIAL use (as a temporary var)
      unsigned int id = s_var_map.size() + 1;
      unsigned ldd_size = crab_domain_params_man::get().boxes_ldd_size();
      if (id >= ldd_size) {
        CRAB_ERROR("The Ldd size of ", ldd_size, " needs to be larger");
      }
      s_var_map.insert(binding_t(v, id));
      return id;
    }
  }

  inline constant_t mk_cst(ikos::z_number k) {
    ikos::q_number qk(k);
    return (constant_t)tvpi_create_cst(qk.get_mpq_t());
  }

  inline constant_t mk_cst(ikos::q_number k) {
    return (constant_t)tvpi_create_cst(k.get_mpq_t());
  }

  // convex approximation
  void convex_approx() {
    m_ldd = lddPtr(get_ldd_man(), Ldd_TermMinmaxApprox(get_ldd_man(), &*m_ldd));
  }

  LddNodePtr convex_approx(LddNodePtr ldd) const {
    return lddPtr(get_ldd_man(), Ldd_TermMinmaxApprox(get_ldd_man(), &*ldd));
  }

  LddNodePtr project(variable_t v) const {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");
    
    std::vector<int> qvars;
    // num_of_vars is shared by all ldd's
    qvars.reserve(num_of_vars() - 1);
    for (auto p : s_var_map.left) {
      if (!(p.first == v)) {
        qvars.push_back(get_var_dim(p.first));
      }
    }
    return lddPtr(get_ldd_man(), Ldd_MvExistAbstract(get_ldd_man(), &*m_ldd,
						     &qvars[0], qvars.size()));
  }

  LddNodePtr join(LddNodePtr v1, LddNodePtr v2) const {
    return lddPtr(get_ldd_man(), Ldd_Or(get_ldd_man(), &*v1, &*v2));
  }

  /** return term for variable v, neg for negation of variable */
  linterm_t term_from_var(variable_t v, bool neg = false) {
    int dim = get_var_dim(v);
    int sgn = neg ? -1 : 1;
    linterm_t term =
        Ldd_GetTheory(get_ldd_man())->create_linterm_sparse_si(&dim, &sgn, 1);
    return term;
  }

  /** return term from SPECIAL variable $0 */
  linterm_t term_from_special_var(bool neg = false) {
    int dim = 0; // reserved for SPECIAL variable
    int sgn = neg ? -1 : 1;
    return Ldd_GetTheory(get_ldd_man())
        ->create_linterm_sparse_si(&dim, &sgn, 1);
  }

  void copy_term(variable_t v, boost::optional<variable_t> x) {
    if (is_top() || is_bottom()) {
      return;
    }

    this->operator-=(v); // remove v before assigning new term

    linterm_t lhs = term_from_var(v);
    linterm_t rhs = (x ? term_from_var(*x) : term_from_special_var());
    
    m_ldd =
        lddPtr(get_ldd_man(), Ldd_TermCopy(get_ldd_man(), &(*m_ldd), lhs, rhs));
    Ldd_GetTheory(get_ldd_man())->destroy_term(lhs);
    Ldd_GetTheory(get_ldd_man())->destroy_term(rhs);

    if (!x) {
      // XXX: if we copy the SPECIAL var to v we forget the SPECIAL
      //      var after we copy.
      crab::CrabStats::count(domain_name() + ".count.forget");
      crab::ScopedCrabStats __st__(domain_name() + ".forget");
      int dim = 0; // SPECIAL variable is dim 0
      m_ldd = lddPtr(get_ldd_man(),
                     Ldd_ExistsAbstract(get_ldd_man(), &*m_ldd, dim));
    }
  }
 
  
  /**
   **  All the expressiveness about numerical operations with boxes
   **  is limited to:
   **     v := a * x + [k.lb(),k.ub()],
   **     where a is a constant, k an interval and x variable
   **  Each numerical operation is reduced to this form.
   */
  void apply_ldd(const variable_t &v, const variable_t &x, number_t a, interval_t k) {
    if (is_top() || is_bottom())
      return;

    linterm_t t = term_from_var(v);
    linterm_t r = term_from_var(x);

    constant_t c = mk_cst(a);

    constant_t kmin = NULL;
    constant_t kmax = NULL;
    if (k.lb().is_finite())
      kmin = mk_cst(*(k.lb().number()));
    if (k.ub().is_finite())
      kmax = mk_cst(*(k.ub().number()));

    m_ldd = lddPtr(get_ldd_man(), Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t,
                                                  r, c, kmin, kmax));

    Ldd_GetTheory(get_ldd_man())->destroy_term(t);
    Ldd_GetTheory(get_ldd_man())->destroy_term(r);
    Ldd_GetTheory(get_ldd_man())->destroy_cst(c);
    if (kmin)
      Ldd_GetTheory(get_ldd_man())->destroy_cst(kmin);
    if (kmax)
      Ldd_GetTheory(get_ldd_man())->destroy_cst(kmax);
  }

  /* Special case of apply_ldd : v := [lb,ub] */
  LddNodePtr apply_interval(const variable_t &v, interval_t ival) {
    assert(!is_bottom());

    constant_t kmin = NULL, kmax = NULL;
    if (boost::optional<number_t> l = ival.lb().number())
      kmin = mk_cst(*l);
    if (boost::optional<number_t> u = ival.ub().number())
      kmax = mk_cst(*u);

    linterm_t t = term_from_var(v);
    LddNodePtr res =
        lddPtr(get_ldd_man(), Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t, NULL,
                                              NULL, kmin, kmax));
    if (kmin)
      get_theory()->destroy_cst(kmin);
    if (kmax)
      get_theory()->destroy_cst(kmax);
    get_theory()->destroy_term(t);
    return res;
  }

  void num_from_ldd_cst(constant_t cst, ikos::z_number &res) const {
    // XXX We know that the theory is tvpi, use its method direclty.
    ikos::q_number q;
    tvpi_cst_set_mpq(q.get_mpq_t(), (tvpi_cst_t)cst);
    res = q.round_to_lower();
  }

  void num_from_ldd_cst(constant_t cst, ikos::q_number &res) const {
    // XXX We know that the theory is tvpi, use its method direclty.
    tvpi_cst_set_mpq(res.get_mpq_t(), (tvpi_cst_t)cst);
  }

  linear_expression_t expr_from_ldd_term(linterm_t term) const {
    linear_expression_t e(0);
    for (size_t i = 0; i < (size_t)get_theory()->term_size(term); i++) {
      number_t k(0); // any value
      num_from_ldd_cst(get_theory()->term_get_coeff(term, i), k);
      variable_t v(getVarName(get_theory()->term_get_var(term, i)));
      e = e + (k * linear_expression_t(v));
    }
    return e;
  }

  ikos::linear_constraint<ikos::z_number, varname_t> cst_from_ldd_strict_cons(
      const ikos::linear_expression<ikos::z_number, varname_t> &l,
      const ikos::linear_expression<ikos::z_number, varname_t> &r) const {
    // l < r <-> l+1 <= r
    ikos::linear_expression<ikos::z_number, varname_t> e = l + 1;
    e = e - r;
    return ikos::linear_constraint<ikos::z_number, varname_t>(
        e, linear_constraint_t::INEQUALITY);
  }

  ikos::linear_constraint<ikos::q_number, varname_t> cst_from_ldd_strict_cons(
      const ikos::linear_expression<ikos::q_number, varname_t> &l,
      const ikos::linear_expression<ikos::q_number, varname_t> &r) const {
    return ikos::linear_constraint<ikos::q_number, varname_t>(
        l - r, linear_constraint_t::STRICT_INEQUALITY);
  }

  linear_constraint_t cst_from_ldd_cons(lincons_t lincons) const {

    linear_expression_t lhs =
        expr_from_ldd_term(get_theory()->get_term(lincons));
    number_t rhs(0); // any value
    num_from_ldd_cst(get_theory()->get_constant(lincons), rhs);

    if (get_theory()->is_strict(lincons)) {
      // lhs < rhs
      return cst_from_ldd_strict_cons(lhs, rhs);
    } else {
      // lhs <= rhs
      linear_expression_t e = lhs - rhs;
      return linear_constraint_t(e, linear_constraint_t::INEQUALITY);
    }
  }

  LddNodePtr make_unit_constraint(number_t coef, variable_t var, kind_t kind, number_t k) {
    assert(coef == 1 || coef == -1);    
    linterm_t term = term_from_var(var, (coef == 1 ? false : true));
    return make_unit_constraint(term, kind, k);
  }

  LddNodePtr make_unit_constraint(number_t coef, kind_t kind, number_t k) {
    assert(coef == 1 || coef == -1);    
    linterm_t term = term_from_special_var((coef == 1 ? false : true));
    return make_unit_constraint(term, kind, k);
  }
  
  // Create a new ldd representing the constraint e
  // where e can be one of
  //    term <= k | term < k | term == k | term != k
  //    where term can be a variable either with positive or negative sign
  //          k is an integer constant
  LddNodePtr make_unit_constraint(linterm_t term, kind_t kind, number_t k) {
    switch(kind) {
    case kind_t::EQUALITY: {
      constant_t c = mk_cst(k);
      // x>=k            
      lincons_t cons1 = get_theory()->create_cons(term, 0 /*non-strict*/, c);
      // x<=k      
      lincons_t cons2 = get_theory()->create_cons(get_theory()->negate_term(term), 0, get_theory()->negate_cst(c));
      LddNodePtr n1 = lddPtr(get_ldd_man(), get_theory()->to_ldd(get_ldd_man(), cons1));
      LddNodePtr n2 = lddPtr(get_ldd_man(), get_theory()->to_ldd(get_ldd_man(), cons2));
      LddNodePtr n = lddPtr(get_ldd_man(), Ldd_And(get_ldd_man(), &*n1, &*n2));
      get_theory()->destroy_lincons(cons1);
      get_theory()->destroy_lincons(cons2);                  
      return n;
    }
    case kind_t::INEQUALITY: {
      // x <= k
      constant_t c = mk_cst(k);
      lincons_t cons = get_theory()->create_cons(term, 0 /*non-strict*/, c);
      LddNodePtr n = lddPtr(get_ldd_man(), get_theory()->to_ldd(get_ldd_man(), cons));
      get_theory()->destroy_lincons(cons);
      return n;
    }
    case kind_t::STRICT_INEQUALITY: {
      // x < k
      constant_t c = mk_cst(k);
      lincons_t cons = get_theory()->create_cons(term, 1 /*strict*/, c);
      LddNodePtr n =
	lddPtr(get_ldd_man(), get_theory()->to_ldd(get_ldd_man(), cons));
      get_theory()->destroy_lincons(cons);
      return n;
    }
    default: {
      assert(kind == kind_t::DISEQUATION);
      // x >= k+1
      auto cons2 = make_unit_constraint(get_theory()->negate_term(term), kind_t::INEQUALITY, -(k+1));
      // x <= k-1
      auto cons1 = make_unit_constraint(term, kind_t::INEQUALITY, k-1);
      return lddPtr(get_ldd_man(), Ldd_Or(get_ldd_man(), &*cons1, &*cons2));
    }
    }
  }

  void apply_unit_constraint(number_t coef, variable_t x, kind_t kind, number_t k) {
    LddNodePtr n = make_unit_constraint(coef, x, kind, k);
    m_ldd = lddPtr(get_ldd_man(), Ldd_And(get_ldd_man(), &*m_ldd, &*n));
  }

  // Extract interval constraints from linear expression exp. This
  // operation uses "at" which loses precision. It should be
  // only used with non-unit expressions.
  void unitcsts_of_exp(const linear_expression_t &exp,
                       std::vector<std::pair<variable_t, number_t>> &lbs,
                       std::vector<std::pair<variable_t, number_t>> &ubs) const {

    number_t unbounded_lbcoeff;
    number_t unbounded_ubcoeff;
    boost::optional<variable_t> unbounded_lbvar;
    boost::optional<variable_t> unbounded_ubvar;
    number_t exp_ub = -(exp.constant());
    std::vector<std::pair<std::pair<number_t, variable_t>, number_t>> pos_terms;
    std::vector<std::pair<std::pair<number_t, variable_t>, number_t>> neg_terms;
    for (auto p : exp) {
      number_t coeff(p.first);
      if (coeff > number_t(0)) {
        variable_t y(p.second);
        // evaluate the variable in the domain
        auto y_lb = at(y).lb();
        if (y_lb.is_infinite()) {
          if (unbounded_lbvar)
            return;
          unbounded_lbvar = y;
          unbounded_lbcoeff = coeff;
        } else {
          number_t ymin(*(y_lb.number()));
          exp_ub -= ymin * coeff;
          pos_terms.push_back({{coeff, y}, ymin});
        }
      } else {
        variable_t y(p.second);
        // evaluate the variable in the domain
        auto y_ub = at(y).ub();
        if (y_ub.is_infinite()) {
          if (unbounded_ubvar)
            return;
          unbounded_ubvar = y;
          unbounded_ubcoeff = -(coeff);
        } else {
          number_t ymax(*(y_ub.number()));
          exp_ub -= ymax * coeff;
          neg_terms.push_back({{-coeff, y}, ymax});
        }
      }
    }

    if (unbounded_lbvar) {
      if (!unbounded_ubvar) {
        // Add bounds for x
        variable_t x(*unbounded_lbvar);
        ubs.push_back({x, exp_ub / unbounded_lbcoeff});
      }
    } else {
      if (unbounded_ubvar) {
        // Bounds for y
        variable_t y(*unbounded_ubvar);
        lbs.push_back({y, -exp_ub / unbounded_ubcoeff});
      } else {
        for (auto pl : neg_terms)
          lbs.push_back(
              {pl.first.second, -exp_ub / pl.first.first + pl.second});
        for (auto pu : pos_terms)
          ubs.push_back({pu.first.second, exp_ub / pu.first.first + pu.second});
      }
    }
  }

  interval_t eval_interval(const linear_expression_t &e) const {
    interval_t r = e.constant();
    for (auto p : e)
      r += p.first * at(p.second);
    return r;
  }

  // Given a constraint a1*x1 + ... + an*xn <= k and pivot xi,
  // it computes the interval:
  //  k - intv (a1*x1 + ... + ai-1*xi-1 + ai+1*xi+1 + ... an*xn)
  interval_t compute_residual(const linear_expression_t &e, variable_t pivot) const {
    interval_t residual(-e.constant());
    for (typename linear_expression_t::const_iterator it = e.begin();
         it != e.end(); ++it) {
      variable_t v = it->second;
      if (!(v == pivot)) {
        residual = residual - (interval_t(it->first) * at(v));
      }
    }
    return residual;
  }

  bool add_linear_leq(const linear_expression_t &e) {
    std::vector<std::pair<variable_t, number_t>> lbs;
    std::vector<std::pair<variable_t, number_t>> ubs;

    unitcsts_of_exp(e, lbs, ubs);
    CRAB_LOG("boxes-extract-leq",
             for (auto p
                  : lbs) {
               crab::outs() << "\t" << p.first << ">=" << p.second << "\n";
             } for (auto p
                    : ubs) {
               crab::outs() << "\t" << p.first << "<=" << p.second << "\n";
             });

    for (auto p : lbs) {
      apply_unit_constraint(-number_t(1) /*coef*/, p.first, kind_t::INEQUALITY,
                            -p.second);
      if (is_bottom()) {
        return false;
      }
    }
    for (auto p : ubs) {
      apply_unit_constraint(number_t(1) /*coef*/, p.first, kind_t::INEQUALITY,
                            p.second);
      if (is_bottom()) {
        return false;
      }
    }
    return true;
  }

  bool add_linear_diseq(const linear_expression_t &e) {
    assert(!is_bottom());

    for (typename linear_expression_t::const_iterator it = e.begin();
         it != e.end(); ++it) {
      variable_t v = it->second;
      interval_t res_i = compute_residual(e, v) / interval_t(it->first);
      boost::optional<number_t> k = res_i.singleton();
      if (k) {
        linear_constraint_t cst(v != *k);
        operator+=(cst);
        if (is_bottom())
          return false;
      } else {
        // XXX: we lose precision here because res_i is an
        // interval which already over-approximates the possible
        // values of the residual.
      }
    }
    return true;
  }

  boxes_domain(LddNodePtr ldd) : m_ldd(ldd) {
    if (crab_domain_params_man::get().boxes_convexify_threshold() > 0) {
      // XXX: the value of convexify_threshold is quite arbitrary. A
      // good value seems around 1000000
      unsigned threshold = num_of_vars() *
	crab_domain_params_man::get().boxes_convexify_threshold();
      // unsigned num_paths = Ldd_PathSize(NULL, &*m_ldd);
      unsigned num_paths = (unsigned)Cudd_CountPath(&*m_ldd);
      if (threshold > 0 && num_paths > threshold) {
        convex_approx();
        CRAB_WARN("ldd size was too large: ", num_paths, ". Made ldd convex.");
      }
    }
  }

  interval_domain_t to_intervals(LddNodePtr &ldd) const {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");

    interval_domain_t res; // top

    if (&*ldd == Ldd_GetFalse(get_ldd_man())) {
      res.set_to_bottom();
      return res;
    }
    if (&*ldd == Ldd_GetTrue(get_ldd_man())) {
      return res;
    }

    LddNodePtr ldd_copy(ldd);
    ldd_copy = convex_approx(ldd_copy);

    auto disjs = to_disjunctive_linear_constraint_system(ldd_copy);
    if (disjs.is_true()) {
      return res;
    } else if (disjs.is_false()) {
      res.set_to_bottom();
      return res;
    } else {
      if (disjs.size() != 1) {
        CRAB_ERROR("Boxes::to_intervals: it should not be disjunctive ", disjs);
      }
      for (auto c : *(disjs.begin())) {
        res += c;
      }
      return res;
    }
  }

  void to_disjunctive_linear_constraint_system_aux(
      LddManager *ldd, LddNode *n, disjunctive_linear_constraint_system_t &e,
      std::vector<int> &list) const {
    /**
     * the latest negative constraint to be printed.
     * is NULL if there isn't one
     */
    lincons_t negc = nullptr;

    DdNode *N = Cudd_Regular(n);

    if (cuddIsConstant(N)) {
      /* n == N here implies that n is one */
      if (n == N) {

        linear_constraint_system_t csts;

        /* for each level */
        for (int i = 0; i < ldd->cudd->size; i++) {
          /* let p be the index of level i */
          int p = ldd->cudd->invperm[i];
          /* let v be the value of p */
          int v = list[p];

          /* skip don't care */
          if (v == 2)
            continue;

          if (v == 0 && ldd->ddVars[p] != nullptr) {
            lincons_t c;
            c = ldd->theory->negate_cons(ldd->ddVars[p]);

            if (negc != nullptr) {
              /* consider negative constraint if it is not implied
                 by c
              */
              if (!ldd->theory->is_stronger_cons(c, negc)) {
                csts += cst_from_ldd_cons(negc);
              }
              ldd->theory->destroy_lincons(negc);
            }

            /* store the current constraint to be considered later */
            negc = c;
            continue;
          }

          /* if there is a negative constraint waiting to be
             considered, conjoin it now
          */
          if (negc != nullptr) {
            csts += cst_from_ldd_cons(negc);
            ldd->theory->destroy_lincons(negc);
            negc = nullptr;
          }

          /* if v is not a don't care but p does not correspond to a
           * constraint, consider it as a Boolean variable */
          if (v != 2 && ldd->ddVars[p] == nullptr) {
            // XXX: I don't know what to do here
            // fprintf(stderr, "%sb%d",(v == 0 ? "!" : " "), p);
          }
          /* v is true */
          else if (v == 1) {
            csts += cst_from_ldd_cons(ldd->ddVars[p]);
          }
        } // end for

        /* if there is a constraint waiting to be considered, do it
           now */
        if (negc != nullptr) {
          csts += cst_from_ldd_cons(negc);
          ldd->theory->destroy_lincons(negc);
          negc = nullptr;
        }

        if (csts.is_true()) {
          // FIXME: if csts is true then e should become true
        } else {
          e += csts;
        }
      }
    } else {
      DdNode *Nv = Cudd_NotCond(cuddT(N), N != n);
      DdNode *Nnv = Cudd_NotCond(cuddE(N), N != n);
      int index = N->index;
      list[index] = 0;
      to_disjunctive_linear_constraint_system_aux(ldd, Nnv, e, list);
      list[index] = 1;
      to_disjunctive_linear_constraint_system_aux(ldd, Nv, e, list);
      list[index] = 2;
    }
    return;
  }

  inline LddNodePtr mk_false(boost::optional<variable_t> x) {
    return (x ? make_unit_constraint(number_t(1), variable_t(*x),
				     linear_constraint_t::EQUALITY, number_t(0)):
	    make_unit_constraint(number_t(1), 
				 linear_constraint_t::EQUALITY, number_t(0)));
  }

  inline linear_constraint_t mk_false_cst(const variable_t &x) {
    return linear_constraint_t(x == 0);
  }

  // Important: Crab follows LLVM semantics where TRUE is 1. Note that
  // in C, TRUE is any value different from 0. Following LLVM
  // semantics avoids many disjunctions in the boolean operations.
  inline LddNodePtr mk_true(boost::optional<variable_t> x) {
    return (x ? make_unit_constraint(number_t(1), variable_t(*x),
				     linear_constraint_t::EQUALITY, number_t(1)):
	        make_unit_constraint(number_t(1),
				     linear_constraint_t::EQUALITY, number_t(1)));
    // return (x ? make_unit_constraint(number_t(1), variable_t(*x),
    // 				     linear_constraint_t::DISEQUATION, number_t(0)):
    // 	        make_unit_constraint(number_t(1),
    // 				     linear_constraint_t::DISEQUATION, number_t(0)));
  }

  inline linear_constraint_t mk_true_cst(const variable_t &x) {
    return linear_constraint_t(x == 1);
    //return linear_constraint_t(x != 0);
  }
  

  // Return a node of the form ITE(y Op z, x=TRUE, x=FALSE)
  LddNodePtr gen_binary_bool(bool_operation_t op, boost::optional<variable_t> x,
                             const variable_t &y, const variable_t &z) {
    switch (op) {
    case OP_BAND: {
      LddNodePtr c = lddPtr(get_ldd_man(),
                            Ldd_And(get_ldd_man(), &*mk_true(y), &*mk_true(z)));
      return lddPtr(get_ldd_man(),
                    Ldd_Ite(get_ldd_man(), &*c, &*mk_true(x), &*mk_false(x)));
    }
    case OP_BOR: {
      LddNodePtr c = lddPtr(get_ldd_man(),
                            Ldd_Or(get_ldd_man(), &*mk_true(y), &*mk_true(z)));
      return lddPtr(get_ldd_man(),
                    Ldd_Ite(get_ldd_man(), &*c, &*mk_true(x), &*mk_false(x)));
    }
    case OP_BXOR: {
      LddNodePtr c = lddPtr(get_ldd_man(),
                            Ldd_Xor(get_ldd_man(), &*mk_true(y), &*mk_true(z)));
      return lddPtr(get_ldd_man(),
                    Ldd_Ite(get_ldd_man(), &*c, &*mk_true(x), &*mk_false(x)));

    }
    }
  }
  
  variable_t getVarName(int v) const {
    auto it = s_var_map.right.find(v);
    if (it != s_var_map.right.end())
      return it->second;
    else {
      CRAB_ERROR("Index ", v, " cannot be mapped back to a variable name");
    }
  }
  
public:
  static void clear_global_state() { s_var_map.clear(); }

  boxes_domain() : m_ldd(lddPtr(get_ldd_man(), Ldd_GetTrue(get_ldd_man()))) {}

  ~boxes_domain() {
    // DdManager *cudd = nullptr;
    // theory_t *theory = nullptr;
    // if (s_ldd_man)  {
    //   cudd = Ldd_GetCudd(s_ldd_man);
    //   theory = Ldd_GetTheory(s_ldd_man);
    //   Ldd_Quit(s_ldd_man);
    // }
    // if (theory) tvpi_destroy_theory(theory);
    // if (cudd) Cudd_Quit(cudd);
  }

  boxes_domain_t make_top() const override {
    boxes_domain_t out(lddPtr(get_ldd_man(), Ldd_GetTrue(get_ldd_man())));
    return out;
  }

  boxes_domain_t make_bottom() const override {
    boxes_domain_t out(lddPtr(get_ldd_man(), Ldd_GetFalse(get_ldd_man())));
    return out;
  }

  void set_to_top() override {
    boxes_domain_t abs(lddPtr(get_ldd_man(), Ldd_GetTrue(get_ldd_man())));
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    boxes_domain_t abs(lddPtr(get_ldd_man(), Ldd_GetFalse(get_ldd_man())));
    std::swap(*this, abs);
  }

  boxes_domain(const boxes_domain_t &other)
      : // m_ldd(other.m_ldd)
        m_ldd(lddPtr(get_ldd_man(), &(*other.m_ldd))) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  boxes_domain(boxes_domain_t &&other) : m_ldd(std::move(other.m_ldd)) {}

  boxes_domain_t &operator=(const boxes_domain_t &other) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &other) {
      // m_ldd = other.m_ldd;
      m_ldd = lddPtr(get_ldd_man(), &(*other.m_ldd));
    }
    return *this;
  }

  bool is_bottom() const override {
    return &*m_ldd == Ldd_GetFalse(get_ldd_man());
  }

  bool is_top() const override { return &*m_ldd == Ldd_GetTrue(get_ldd_man()); }

  bool operator<=(const boxes_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    bool res = Ldd_TermLeq(get_ldd_man(), &(*m_ldd), &(*other.m_ldd));

    CRAB_LOG("boxes", boxes_domain_t left(*this); boxes_domain_t right(other);
             crab::outs() << "Check if " << left << " <= " << right << " ---> "
                          << res << "\n";);
    return res;
  }

  void operator|=(const boxes_domain_t &other) override {
    *this = *this | other;
  }

  boxes_domain_t operator|(const boxes_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    return boxes_domain_t(join(m_ldd, other.m_ldd));
  }

  void operator&=(const boxes_domain_t &other) override {
    *this = *this & other;
  }
  
  boxes_domain_t operator&(const boxes_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    return boxes_domain_t(
        lddPtr(get_ldd_man(), Ldd_And(get_ldd_man(), &*m_ldd, &*other.m_ldd)));
  }

  boxes_domain_t operator||(const boxes_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // It is not necessarily true that the new value is bigger
    // than the old value so we apply
    // widen(old, new) = widen(old,(join(old,new)))
    LddNodePtr v = join(m_ldd, other.m_ldd);
    LddNodePtr w =
        lddPtr(get_ldd_man(), Ldd_BoxWiden2(get_ldd_man(), &*m_ldd, &*v));
    boxes_domain_t res(w);
    CRAB_LOG("boxes", crab::outs() << "Performed widening \n"
                                   << "**" << *this << "\n"
                                   << "** " << other << "\n"
                                   << "= " << res << "\n";);
    return res;
  }

  boxes_domain_t widening_thresholds(
      const boxes_domain_t &other,
      const thresholds<number_t> & /*ts*/) const override {
    // CRAB_WARN(" boxes widening operator with thresholds not implemented");
    return (*this || other);
  }

  boxes_domain_t operator&&(const boxes_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    boxes_domain_t res(*this & other);
    // CRAB_WARN(" boxes narrowing operator replaced with meet");
    CRAB_LOG("boxes", crab::outs() << "Performed narrowing \n"
                                   << "**" << *this << "\n"
                                   << "** " << other << "\n"
                                   << "= " << res << "\n";);
    return res;
  }

  void operator-=(const variable_t &var) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top())
      return;

    int id = get_var_dim(var);
    m_ldd =
        lddPtr(get_ldd_man(), Ldd_ExistsAbstract(get_ldd_man(), &*m_ldd, id));
  }

  void forget(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top())
      return;

    std::vector<int> qvars;
    qvars.reserve(variables.size());
    for (variable_t v : variables) {
      qvars.push_back(get_var_dim(v));
    }

    m_ldd = lddPtr(get_ldd_man(), Ldd_MvExistAbstract(get_ldd_man(), &*m_ldd,
                                                      &qvars[0], qvars.size()));
  }

  void project(const variable_vector_t &variables) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top())
      return;

    if (variables.empty()) {
      set_to_top();
      return;
    }

    std::set<variable_t> s1, s2;
    variable_vector_t s3;
    for (auto p : s_var_map.left)
      s1.insert(p.first);
    s2.insert(variables.begin(), variables.end());
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(s3));
    forget(s3);
  }

  void expand(const variable_t &v, const variable_t &new_v) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_top() || is_bottom())
      return;

    // new_v should be completely unconstrained
    this->operator-=(new_v);

    linterm_t lnew_v = term_from_var(new_v);
    linterm_t lv = term_from_var(v);

    m_ldd = lddPtr(get_ldd_man(),
                   Ldd_TermCopy(get_ldd_man(), &(*m_ldd), lnew_v, lv));
    Ldd_GetTheory(get_ldd_man())->destroy_term(lnew_v);
    Ldd_GetTheory(get_ldd_man())->destroy_term(lv);
  }

  void operator+=(const linear_constraint_t &cst) {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    CRAB_LOG("boxes", crab::outs() << "--- assume(" << cst << ") --> ";);

    if (is_bottom() || cst.is_tautology()) {
      return;
    }

    if (cst.is_contradiction()) {
      set_to_bottom();
      return;
    }

    linear_expression_t exp = cst.expression();
    unsigned int size = exp.size();
    if (size == 1) {
      // Handle precisely unit constraints
      auto it = exp.begin();
      number_t cx = it->first;
      if (cx == number_t(1) || cx == number_t(-1)) {
        number_t k = -exp.constant();
        variable_t x = it->second;
        apply_unit_constraint(cx, x, cst.kind(), k);
      } else {
        CRAB_WARN("non-unit coefficients not implemented in boxes");
      }
    }
    // Lose precision with non-unit constraints
    else if (size >= 2) {
      if (cst.is_disequation()) {
        if (!add_linear_diseq(cst.expression())) {
          CRAB_LOG("boxes", crab::outs() << "_|_\n";);
          return;
        }
      } else if (cst.is_inequality()) {
        if (!add_linear_leq(cst.expression())) {
          CRAB_LOG("boxes", crab::outs() << "_|_\n";);
          return;
        }
      } else if (cst.is_strict_inequality()) {
        // We try to convert a strict to non-strict.
        auto nc =
            ikos::linear_constraint_impl::strict_to_non_strict_inequality(cst);
        if (nc.is_inequality()) {
          // here we succeed
          if (!add_linear_leq(nc.expression())) {
            CRAB_LOG("boxes", crab::outs() << "_|_\n";);
            return;
          }
        }
      } else if (cst.is_equality()) {
        linear_expression_t exp = cst.expression();
        // split equality into two inequalities
        if (!add_linear_leq(exp) || !add_linear_leq(-exp)) {
          CRAB_LOG("boxes", crab::outs() << "_|_\n";);
          return;
        }
      }
    }

    CRAB_LOG("boxes", crab::outs() << *this << "\n";);
  }

  DEFAULT_ENTAILS(boxes_domain_t)
  
  void normalize() override {}

  void minimize() override {}

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const boxes_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    assert(from.size() == to.size());

    if (is_top() || is_bottom())
      return;

    CRAB_WARN(domain_name(), "::rename not implemented");
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    if (!is_bottom()) {
      for (auto cst : csts) {
	operator+=(cst);
      }
    }
  }

  void set(const variable_t &v, interval_t ival) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      m_ldd = apply_interval(v, ival);
    }
  }

  virtual interval_t operator[](const variable_t &v) override {
    return at(v);
  }

  virtual interval_t at(const variable_t &v) const override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");

    if (is_bottom()) {
      return interval_t::bottom();
    }
    if (is_top()) {
      return interval_t::top();
    }

    LddNodePtr ldd = project(v);
    interval_domain_t intv = to_intervals(ldd);
    interval_t i = intv[v];
    CRAB_LOG("boxes-project", crab::outs() << "Before projecting on " << v
                                           << ": " << *this << "\n"
                                           << "Projection " << i << "\n";);
    return i;
  }
  
  // x := e
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (is_bottom()) {
      return;
    }

    if (e.is_constant()) {
      // Handle precisely if e is constant
      constant_t c = mk_cst(e.constant());
      linterm_t t = term_from_var(x);
      m_ldd = lddPtr(get_ldd_man(), Ldd_TermReplace(get_ldd_man(), &(*m_ldd), t,
                                                    NULL, NULL, c, c));
      get_theory()->destroy_cst(c);
      get_theory()->destroy_term(t);
    } else if (e.size() == 1) {
      // Handle precise if e is a*x + k
      auto it = e.begin();
      apply_ldd(x, it->second, it->first, e.constant());
    } else {
      crab::ScopedCrabStats __st__(domain_name() + ".assign.converting_to_intervals");
      // Lose precision in the general case
      // XXX: we can probably do a bit better here
      m_ldd = apply_interval(x, eval_interval(e));
    }

    CRAB_LOG("boxes", crab::outs() << "--- " << x << ":=" << e << "\n"
                                   << *this << "\n";);
  }

  // x := y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom()) {
      return;
    }

    if (op >= OP_ADDITION && op <= OP_SUBTRACTION /*OP_MULTIPLICATION*/) {
      // FIXME: we skip multiplication because it seems imprecise.  It
      // produces imprecise lower bounds in unittests-boxes.
      switch (op) {
      case OP_ADDITION:
        apply_ldd(x, y, 1, k);
        break;
      case OP_SUBTRACTION:
        apply_ldd(x, y, 1, -k);
        break;
      case OP_MULTIPLICATION:
        apply_ldd(x, y, k, number_t(0));
        break;
      default:
        CRAB_ERROR("Unexpected operator ", op);
      }
      CRAB_LOG("boxes", crab::outs() << "--- " << x << " := " << y << " " << op
                                     << " " << k << "\n"
                                     << *this << "\n";);
    } else {
      // Convert to intervals
      interval_t yi(this->operator[](y));
      interval_t zi(k);
      interval_t xi(interval_t::bottom());
      switch (op) {
      case OP_MULTIPLICATION:
        xi = yi * zi;
        break;
      case OP_SDIV:
        xi = yi / zi;
        break;
      case OP_UDIV:
        xi = yi.UDiv(zi);
        break;
      case OP_SREM:
        xi = yi.SRem(zi);
        break;
      case OP_UREM:
        xi = yi.URem(zi);
        break;
      default:
        CRAB_ERROR("Unexpected operator ", op);
      }
      set(x, xi);
    }
  }

  // x := y op z
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom()) {
      return;
    }

    // --- if z is a singleton we do not lose precision
    interval_t zi = this->operator[](z);
    if (auto k = zi.singleton()) {
      apply(op, x, y, *k);
      return;
    }

    // --- if y or z is top then we give up
    interval_t yi = this->operator[](y);
    if (yi.is_top() || zi.is_top()) {
      this->operator-=(x);
      goto apply_end;
    }

    // --- if y is a singleton we do not lose precision
    if (auto k = yi.singleton()) {
      if (op == OP_ADDITION || op == OP_MULTIPLICATION) {
        apply(op, x, z, *k);
        return;
      } else if (op == OP_SUBTRACTION) {
        // x = yi - z <--> x = -z + [yi.lb,yi.ub]
        apply_ldd(x, z, -1, *k);
        goto apply_end;
      }
    }

    // XXX: ideally we would like to extract all the interval
    // constraints from x:= y op z, forget x and, then add all
    // the interval constraints. However, we would need to
    // project on y and z while keeping all the disjunctive
    // information about these variables. Note that the method
    // operator[] loses all the disjunctive information.

    // -- We need to abstract either y or z and lose some
    // precision. We abstract the one with the smaller interval.
    if (op == OP_ADDITION) {
      if (yi <= zi) { // abstract y
        // x = yi + z <--> x = z + yi
        apply_ldd(x, z, 1, yi);
      } else { // abstract z
        // x = y + zi
        apply_ldd(x, y, 1, zi);
      }
      goto apply_end;
    }

    // -- We need to abstract either y or z and lose some
    // precision. We abstract the one with the smaller interval.
    if (op == OP_SUBTRACTION) {
      if (yi <= zi) { // abstract y
        // x = yi - z <-->  x = -z + yi
        apply_ldd(x, z, -1, yi);
      } else { // abstract z
        // x = y - zi
        apply_ldd(x, y, 1, zi * number_t(-1));
      }
      goto apply_end;
    }

    // --- we go to intervals and lose precision
    switch (op) {
    case OP_MULTIPLICATION:
      set(x, yi * zi);
      break;
    case OP_SDIV:
      set(x, yi / zi);
      break;
    case OP_UDIV:
      set(x, yi.UDiv(zi));
      break;
    case OP_SREM:
      set(x, yi.SRem(zi));
      break;
    case OP_UREM:
      set(x, yi.URem(zi));
      break;
    default:
      CRAB_ERROR("Unexpected operator ", op);
    }

  apply_end:
    CRAB_LOG("boxes", crab::outs() << "--- " << x << " := " << y << " " << op
                                   << " " << z << "\n"
                                   << *this << "\n";);
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const boxes_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    CRAB_LOG("boxes", crab::outs() << "Backward " << x << ":=" << e
                                   << "\n\tPOST=" << *this << "\n");
    BackwardAssignOps<boxes_domain_t>::assign(*this, x, e, invariant);
    CRAB_LOG("boxes", crab::outs() << "\tPRE=" << *this << "\n");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const boxes_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    CRAB_LOG("boxes", crab::outs() << "Backward " << x << ":=" << y << op << z
                                   << "\n\tPOST=" << *this << "\n");
    BackwardAssignOps<boxes_domain_t>::apply(*this, op, x, y, z, invariant);
    CRAB_LOG("boxes", crab::outs() << "\tPRE=" << *this << "\n");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const boxes_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    CRAB_LOG("boxes", crab::outs() << "Backward " << x << ":=" << y << op << z
                                   << "\n\tPOST=" << *this << "\n");
    BackwardAssignOps<boxes_domain_t>::apply(*this, op, x, y, z, invariant);
    CRAB_LOG("boxes", crab::outs() << "\tPRE=" << *this << "\n");
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    // since reasoning about infinite precision we simply assign and
    // ignore the widths.
    assign(dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;

    // Convert to intervals and perform the operation
    interval_t yi = this->operator[](y);
    interval_t zi = this->operator[](z);
    interval_t xi = interval_t::bottom();
    switch (op) {
    case OP_AND:
      xi = yi.And(zi);
      break;
    case OP_OR:
      xi = yi.Or(zi);
      break;
    case OP_XOR:
      xi = yi.Xor(zi);
      break;
    case OP_SHL:
      xi = yi.Shl(zi);
      break;
    case OP_LSHR:
      xi = yi.LShr(zi);
      break;
    case OP_ASHR:
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

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi(k);
    interval_t xi = interval_t::bottom();
    switch (op) {
    case OP_AND:
      xi = yi.And(zi);
      break;
    case OP_OR:
      xi = yi.Or(zi);
      break;
    case OP_XOR:
      xi = yi.Xor(zi);
      break;
    case OP_SHL:
      xi = yi.Shl(zi);
      break;
    case OP_LSHR:
      xi = yi.LShr(zi);
      break;
    case OP_ASHR:
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }

  ////////
  //// boolean operations
  ////////
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &cst) override {
    if (!m_bool_reasoning) {
      return;
    }
    
    crab::CrabStats::count(domain_name() + ".count.assign_bool_cst");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst");

    auto is_interval_cst = [](const linear_constraint_t &c) {
    	   return (c.expression().size() == 1 &&
    		   (c.expression().begin()->first == number_t(1) ||
    		    c.expression().begin()->first == number_t(-1))); };
    
    if (is_bottom()) {
      return;
    } 

    // Important to forget first lhs
    this->operator-=(lhs);
    
    if (cst.is_tautology()) {
      this->operator+=(mk_true_cst(lhs));
    } else if (cst.is_contradiction()) {
      this->operator+=(mk_false_cst(lhs));
    } else {

      if (is_interval_cst(cst)) {
	crab::CrabStats::count(domain_name() + ".count.assign_bool_cst.interval");
	crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst.interval");
      	/// Interval constraint
      	linear_expression_t exp = cst.expression();
        number_t cx = exp.begin()->first;
      	variable_t vx = exp.begin()->second;
      	number_t k = -exp.constant();
      	// m_ldd &= ite (cst, lhs = TRUE, lhs = FALSE);
      	LddNodePtr ldd_cst = make_unit_constraint(cx, vx, cst.kind(), k);
      	LddNodePtr ldd_rhs = lddPtr(get_ldd_man(), Ldd_Ite(get_ldd_man(), &*(ldd_cst),
      							   &*mk_true(lhs), &*mk_false(lhs)));
      	m_ldd = lddPtr(get_ldd_man(), Ldd_And(get_ldd_man(), &*m_ldd, &*ldd_rhs));
      } else {
      	/// Non-interval constraint
	crab::CrabStats::count(domain_name() + ".count.assign_bool_cst.non-interval");
	crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_cst.non-interval");
	// Note that we could project away irrelevant variables before
	// performing the join but this would lose precision.
        boxes_domain_t tt(*this);
        boxes_domain_t ff(*this);
        tt += cst;
        tt += mk_true_cst(lhs);
        ff += cst.negate();
        ff += mk_false_cst(lhs);
        if (::crab::CrabSanityCheckFlag) {
          if (tt.is_bottom() && ff.is_bottom()) {
            CRAB_ERROR("Something went wrong in assign_bool_cst ", lhs,
                       " := ", cst, " with ", *this);
          }
        }
        *this = (tt | ff);
      }
    }
      
    CRAB_LOG("boxes", crab::outs() << lhs << ":= "
                                   << "(" << cst << ")\n"
                                   << *this << "\n");
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &cst) override {
    if (is_bottom())
      return;

    this->operator-=(lhs);
    CRAB_WARN(
        "boxes boolean assignment of reference constraints not implemented");
  }

  void assign_bool_var(const variable_t &x, const variable_t &y,
                       bool is_not_y) override {
    if (!m_bool_reasoning)
      return;
    
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (is_bottom()) {
      return;
    }
    
    if (is_not_y) {
      // forget left-hand side
      operator-=(x);

      if (x == y) {
	// TODO
	CRAB_WARN("boxes skipped ", x, " not(", y, ")");	
	return;
      }
      
      // Note that we could project away irrelevant variables before
      // performing the join but this would lose precision.
      boxes_domain_t inv1(*this);
      boxes_domain_t inv2(*this);
      inv1 += mk_true_cst(y);
      inv1 += mk_false_cst(x);
      inv2 += mk_false_cst(y);
      inv2 += mk_true_cst(x);
      *this = (inv1 | inv2);
    } else {
      apply_ldd(x, y, number_t(1), number_t(0));      
    }

    CRAB_LOG("boxes", crab::outs() << x << ":=";
             if (is_not_y) crab::outs() << "not(" << y << ")";
             else crab::outs() << y; crab::outs() << "\n"
                                                  << *this << "\n");
  }

  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {

    if (!m_bool_reasoning)
      return;

    crab::CrabStats::count(domain_name() + ".count.apply_bin_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_bin_bool");

    if (is_bottom()) {
      return;
    }
    
    // XXX: if *lhs is null then it represents the SPECIAL
    // variable $0.
    boost::optional<variable_t> lhs;

    if (!(x == y) && !(x == z)) {
      lhs = boost::optional<variable_t>(x);
      // XXX: x does not appear on the rhs so we can remove it
      // without losing precision.
      this->operator-=(x);
    }

    auto add_assumptions = [&y,&z, this]() {
			     // the Ldd does not know that y and z can
			     // only be either 0 or 1.  The Ite term
			     // needs to negate the values of x and
			     // y. So this negation can add a
			     // disjunction to express != 0 or != 1.
			     this->operator+=(y >= 0);
			     this->operator+=(y <= 1);
			     this->operator+=(z >= 0);
			     this->operator+=(z <= 1);
			   };


    m_ldd = lddPtr(get_ldd_man(),
		   Ldd_And(get_ldd_man(), &*m_ldd,
			   &*gen_binary_bool(op, lhs, y, z)));
    add_assumptions();

    if ((x == y) || (x == z)) {
      // XXX: if we are here we added ite(y op z, $0 >= 1, $0 <= 0);
      // so we still need to assign $0 to x:
      copy_term(x, boost::optional<variable_t>());
    }

    CRAB_LOG("boxes", crab::outs()
                          << x << ":= " << y << " " << op << " " << z << "\n"
                          << *this << "\n");
  }

  void assume_bool(const variable_t &x, bool is_negated) override {
    if (!m_bool_reasoning)
      return;

    crab::CrabStats::count(domain_name() + ".count.assume_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".assume_bool");

    if (is_bottom()) {
      return;
    }
    
    m_ldd = lddPtr(get_ldd_man(),
                   Ldd_And(get_ldd_man(), &*m_ldd,
                           ((is_negated) ? &*mk_false(x) : &*mk_true(x))));

    CRAB_LOG("boxes",
             if (!is_negated) {
               crab::outs() << "--- bool_assume(" << x << ")"
                            << "\n"
                            << *this << "\n";
             } else {
               crab::outs() << "--- bool_assume(not(" << x << "))"
                            << "\n"
                            << *this << "\n";
             });
  }

  // Backward boolean operations
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const boxes_domain_t &inv) override {
    if (is_bottom()) {
      return;
    }

    this->operator-=(lhs);
    CRAB_WARN("boxes backward boolean assignment not implemented");
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const boxes_domain_t &inv) override {
    if (is_bottom()) {
      return;
    }

    this->operator-=(lhs);
    CRAB_WARN("boxes backward boolean assignment of reference constraints not "
              "implemented");
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const boxes_domain_t &inv) override {
    if (is_bottom()) {
      return;
    }

    this->operator-=(lhs);
    CRAB_WARN("boxes backward boolean assignment not implemented");
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const boxes_domain_t &inv) override {
    if (is_bottom()) {
      return;
    }

    this->operator-=(x);
    CRAB_WARN("boxes backward boolean apply not implemented");
  }

  /// Boxes domain implements only standard abstract operations of a
  /// numerical domain plus boolean operations. It is intended to be
  /// used as a leaf domain in the hierarchy of domains.
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(boxes_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(boxes_domain_t)

  DEFAULT_SELECT(boxes_domain_t)
  DEFAULT_SELECT_BOOL(boxes_domain_t)
  DEFAULT_WEAK_ASSIGN(boxes_domain_t)
  DEFAULT_WEAK_BOOL_ASSIGN(boxes_domain_t)
  
  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;

    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
    } else if (is_top()) {
      csts += linear_constraint_t::get_true();
    } else {
      // --- produce convex approximation
      LddNodePtr ldd(m_ldd);
      ldd = convex_approx(ldd);
      // --- extract linear inequalities from the convex ldd
      auto disjs = to_disjunctive_linear_constraint_system(ldd);
      if (disjs.is_false()) {
        csts += linear_constraint_t::get_false();
      } else if (disjs.is_true()) {
        csts += linear_constraint_t::get_true();
      } else {
        if (disjs.size() != 1) {
          CRAB_ERROR("Boxes::to_linear_constraint_system: it should not be "
                     "disjunctive ",
                     disjs);
        }
        for (auto c : *(disjs.begin())) {
          csts += c;
        }
      }
    }

    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system(const LddNodePtr &ldd) const {
    LddManager *ldd_man = get_ldd_man();

    std::vector<int> list;
    list.reserve(ldd_man->cudd->size);
    for (int i = 0; i < ldd_man->cudd->size; i++)
      list.push_back(2);

    disjunctive_linear_constraint_system_t r;
    to_disjunctive_linear_constraint_system_aux(ldd_man, ldd.get(), r, list);
    return r;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    return to_disjunctive_linear_constraint_system(m_ldd);
  }

  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    if (is_top()) {
      o << "{}";
    } else if (is_bottom()) {
      o << "_|_";
    } else {
#if 0
      auto r = to_linear_constraint_system();
#else      
      auto r = to_disjunctive_linear_constraint_system();
#endif
      o << r;
    }
  }

  std::string domain_name() const override { return "Boxes"; }
};

template <typename N, typename V>
LddManager *boxes_domain<N, V>::s_ldd_man = nullptr;

template <typename N, typename V>
typename boxes_domain<N, V>::var_map_t boxes_domain<N, V>::s_var_map;

template <typename Number, typename VariableName>
class special_domain_traits<boxes_domain<Number, VariableName>> {
public:
  static void clear_global_state(void) {
    boxes_domain<Number, VariableName>::clear_global_state();
  }
};

// template <typename Number, typename VariableName>
// class checker_domain_traits<boxes_domain<Number, VariableName>> {
// public:
//   using this_type = boxes_domain<Number, VariableName>;
//   using linear_constraint_t = typename this_type::linear_constraint_t;
//   using disjunctive_linear_constraint_system_t =
//       typename this_type::disjunctive_linear_constraint_system_t;
//   static bool entail(this_type &lhs,
//                      const disjunctive_linear_constraint_system_t &rhs) {
//     // -- trivial cases first
//     if (rhs.is_false()) {
//       return false;
//     } else if (rhs.is_true()) {
//       return true;
//     } else if (lhs.is_bottom()) {
//       return true;
//     } else if (lhs.is_top()) {
//       return false;
//     }
//     this_type inv;
//     inv.set_to_bottom();
//     for (auto const &csts : rhs) {
//       this_type conj;
//       conj += csts;
//       inv |= conj;
//     }
//     return (lhs & inv.complement()).is_bottom();
//   }
//   static bool entail(const disjunctive_linear_constraint_system_t &lhs,
//                      this_type &rhs) {
//     // -- trivial cases first
//     if (rhs.is_bottom()) {
//       return false;
//     } else if (rhs.is_top()) {
//       return true;
//     } else if (lhs.is_false()) {
//       return true;
//     } else if (lhs.is_true()) {
//       return false;
//     }
//     this_type inv;
//     inv.set_to_bottom();
//     for (auto const &csts : lhs) {
//       this_type conj;
//       conj += csts;
//       inv |= conj;
//     }
//     return (inv & rhs.complement()).is_bottom();
//   }
//   static bool intersect(this_type &inv, const linear_constraint_t &cst) {
//     // default code

//     // -- trivial cases first
//     if (inv.is_bottom() || cst.is_contradiction())
//       return false;
//     if (inv.is_top() || cst.is_tautology())
//       return true;

//     this_type cst_inv;
//     cst_inv += cst;
//     return (!(cst_inv & inv).is_bottom());
//   }
// };
  
} // namespace domains
} // namespace crab
#endif /* HAVE_LDD */

namespace crab {
namespace domains {
template <typename Number, typename VariableName>
struct abstract_domain_traits<boxes_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
};
} // namespace domains
} // namespace crab
