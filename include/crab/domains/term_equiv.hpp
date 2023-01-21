/*******************************************************************************
 * Abstract domain based on the paper "An Abstract Domain of
 * Uninterpreted Functions" by Gange, Navas, Schachte, Sondergaard,
 * and Stuckey published in VMCAI'16.
 *
 * Author: Graeme Gange (gkgange@unimelb.edu.au)
 *
 * Contributors: Jorge A. Navas (jorge.navas@sri.com)
 *
 * The "term" domain is a customized reduction between a Herbrand
 * domain which maps variables to Herbrand terms
 * (https://en.wikipedia.org/wiki/Herbrand_structure) and a numerical
 * domain. The Herbrand domain infers equalities between variables
 * (typically not inferable by common numerical domains) which are
 * used to improve the numerical domain.
 * 
 * Note about how to read printed invariants.
 *
 * The term domain will print invariants in this form:
 * 
 * {x -> t1[_x1], y -> t0[_x0], z -> t0[_x0]}{_x0 -> [5, 10]; _x1 -> [10, 20]}
 *
 * The 1st component {x -> t1[_x1], y -> t0[_x0], z -> t0[_x0]}
 * contains the term equivalences while the 2nd component is the
 * numerical domain. For the 1st component, each program variable
 * (x,y,z) is mapped to t[v] where t represents some term (the shape
 * of the term is not printed) and v is a ghost variable that can
 * appear on the 2nd component (the numerical abstract state).  The
 * relevant information in the 1st compoment is when two variables are
 * mapped to the same term (e.g., y and z in our example), meaning
 * that they are equal. So, reading our sample invariant, we know that
 * y and z are equal and their values range in [5,10]. About x the
 * only thing we know is that its value ranges in [10,20].
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/domains/term/inverse.hpp>
#include <crab/domains/term/simplify.hpp>
#include <crab/domains/term/term_expr.hpp>
#include <crab/domains/term/term_operators.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>
#include <crab/types/varname_factory.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/optional.hpp>

#define USE_TERM_INTERVAL_NORMALIZER

#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wsign-compare"

namespace crab {
namespace domains {

namespace term {
template <class Num, class VName, class Abs> class TDomInfo {
public:
  using Number = Num;
  using VariableName = VName;
  using variable_t = variable<Num, VName>;
  using Alloc = crab::var_factory_impl::str_var_alloc_col;
  using domain_t = Abs;

  static_assert(
      std::is_same<typename Abs::varname_t, Alloc::varname_t>::value,
      "Base abstract domain must use str_var_alloc_col as varname factory");
};
} // namespace term

template <class Info, class Abs> class TermNormalizer;

template <typename Info>
class term_domain final : public abstract_domain_api<term_domain<Info>> {
  friend class TermNormalizer<Info, typename Info::domain_t>;

  // Number and VariableName can be different from
  // dom_t::number_t and dom_t::varname_t although currently
  // Number and dom_t::number_t must be the same type.
  using Number = typename Info::Number;
  using VariableName = typename Info::VariableName;
  using dom_t = typename Info::domain_t;
  using dom_var_t = typename dom_t::variable_t;
  using dom_var_alloc_t = typename Info::Alloc;
  using dom_varname_t = typename dom_var_alloc_t::varname_t;
  using domvar_set_t = ikos::patricia_tree_set<dom_var_t>;
  using bound_t = ikos::bound<Number>;

  using term_domain_t = term_domain<Info>;
  using abstract_domain_t = abstract_domain_api<term_domain_t>;

public:
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using number_t = Number;
  using varname_t = VariableName;

private:
  using ttbl_t = term::term_table<number_t, term::term_operator_t>;
  using term_id_t = typename ttbl_t::term_id_t;
  using term_t = typename ttbl_t::term_t;

  // WARNING: assumes the underlying domain uses the same number type.
  using dom_number = typename Info::Number;
  using dom_lincst_t = typename dom_t::linear_constraint_t;
  using dom_linsys_t = typename dom_t::linear_constraint_system_t;
  using dom_linexp_t = typename dom_t::linear_expression_t;
  using linterm_t = typename linear_expression_t::component_t;
  using term_map_t = boost::container::flat_map<term_id_t, dom_var_t>;
  using rev_map_t = boost::container::flat_map<dom_var_t, variable_t>;
  using term_set_t = boost::container::flat_set<term_id_t>;
  using var_map_t = boost::container::flat_map<variable_t, term_id_t>;
  // the reverse of var_map: from term_id_t to a set of variable_t
  using var_set_t = std::set<variable_t>;
  using rev_var_map_t = boost::container::flat_map<term_id_t, var_set_t>;
  using simplifier_t = term::NumSimplifier<number_t>;

  bool _is_bottom;
  // Uses a single state of the underlying domain.
  ttbl_t _ttbl;
  dom_t _impl;
  dom_var_alloc_t _alloc;
  var_map_t _var_map;
  rev_var_map_t _rev_var_map; // to extract equalities efficiently
  term_map_t _term_map;
  term_set_t changed_terms;

  term_domain(bool is_top) : _is_bottom(!is_top) {}

  term_domain(dom_var_alloc_t &&alloc, var_map_t &&vm, rev_var_map_t &&rvm,
              ttbl_t &&tbl, term_map_t &&tmap, dom_t &&impl)
      : _is_bottom((impl.is_bottom()) ? true : false), _ttbl(std::move(tbl)),
        _impl(std::move(impl)), _alloc(std::move(alloc)),
        _var_map(std::move(vm)), _rev_var_map(std::move(rvm)),
        _term_map(std::move(tmap)) {
    check_terms(__LINE__);
  }

  // x = y op [lb,ub]
  term_id_t term_of_itv(bound_t lb, bound_t ub) {
    boost::optional<number_t> n_lb = lb.number();
    boost::optional<number_t> n_ub = ub.number();

    if (n_lb && n_ub && (*n_lb == *n_ub))
      return term_of_const(*n_lb);

    term_id_t t_itv = _ttbl.fresh_var();
    dom_var_t dom_itv = domvar_of_term(t_itv);
    _impl.set(dom_itv, interval_t(lb, ub));
    return t_itv;
  }

  void apply(dom_t &dom, term::term_operator_t op, dom_var_t x, dom_var_t y,
             dom_var_t z) {
    if (auto top = term::conv2arith(op)) {
      dom.apply(*top, x, y, z);
    } else if (auto top = term::conv2bitwise(op)) {
      dom.apply(*top, x, y, z);
    } else if (op == term::TERM_OP_FUNCTION) {
      // uninterpreted function: do nothing in the underlying
      // numerical domain.
    } else {
      CRAB_ERROR("unsupported binary operator ", op);
    }
  }

  // Apply a given functor in the underlying domain.
  // GKG: Looks the current implementation could actually
  // lose information; as it's not taking the meet with
  // the current value.
  void eval_ftor(dom_t &dom, ttbl_t &tbl, term_id_t t) {
    // Get the term info.
    term_t *t_ptr = tbl.get_term_ptr(t);

    // Only apply functors.
    if (t_ptr->kind() == term::TERM_APP) {
      term::term_operator_t op = term::term_ftor(t_ptr);

      const std::vector<term_id_t> &args(term::term_args(t_ptr));
      assert(args.size() == 2);
      apply(dom, op, domvar_of_term(t), domvar_of_term(args[0]),
            domvar_of_term(args[1]));
    }
  }

  void eval_ftor_down(dom_t &dom, ttbl_t &tbl, term_id_t t) {
    // Get the term info.
    term_t *t_ptr = tbl.get_term_ptr(t);

    // Only apply functors.
    if (t_ptr->kind() == term::TERM_APP) {
      term::term_operator_t op = term::term_ftor(t_ptr);
      const std::vector<term_id_t> &args(term::term_args(t_ptr));
      assert(args.size() == 2);

      if (boost::optional<arith_operation_t> arith_op = term::conv2arith(op)) {
        term::InverseOps<dom_number, dom_var_t, dom_t>::apply(
            dom, *arith_op, domvar_of_term(t), domvar_of_term(args[0]),
            domvar_of_term(args[1]));
      }
    }
  }

  dom_t eval_ftor_copy(dom_t &dom, ttbl_t &tbl, term_id_t t) {
    dom_t ret = dom;
    eval_ftor(ret, tbl, t);
    return ret;
  }

  void check_terms(int line) const {
    CRAB_LOG(
        "terms-check-terms",
        for (auto const &p
             : _var_map) {
          if (!(p.second < _ttbl.size())) {
            CRAB_ERROR("term_equiv.hpp at line=", line, ": ",
                       "term id is not the table term");
          }
        }

        for (auto kv
             : _rev_var_map) {
          for (auto v : kv.second) {
            auto it = _var_map.find(v);
            if (it->second != kv.first) {
              CRAB_ERROR("term_equiv.hpp at line=", line, ": ", v,
                         " is mapped to t", it->second,
                         " but the reverse map says that should be t",
                         kv.first);
            }
          }
        });
  }

  void deref(term_id_t t) {
    std::vector<term_id_t> forgotten;
    _ttbl.deref(t, forgotten);
    for (term_id_t f : forgotten) {
      typename term_map_t::iterator it(_term_map.find(f));
      if (it != _term_map.end()) {
        _impl -= (*it).second;
        _term_map.erase(it);
      }
    }
  }

  /* Begin manipulate the reverse variable map */
  void add_rev_var_map(rev_var_map_t &rvmap, term_id_t t, variable_t v) const {
    auto it = rvmap.find(t);
    if (it != rvmap.end()) {
      it->second.insert(v);
    } else {
      var_set_t varset;
      varset.insert(v);
      rvmap.insert(std::make_pair(t, varset));
    }
  }

  void remove_rev_var_map(term_id_t t, const variable_t &v) {
    auto it = _rev_var_map.find(t);
    if (it != _rev_var_map.end()) {
      it->second.erase(v);
      if (it->second.empty()) {
        _rev_var_map.erase(it);
      }
    }
  }
  /* End manipulate the reverse variable map */

  void rebind_var(const variable_t &x, term_id_t tx) {
    _ttbl.add_ref(tx);

    auto it(_var_map.find(x));
    if (it != _var_map.end()) {
      remove_rev_var_map((*it).second, x);
      deref((*it).second);
      _var_map.erase(it);
    }
    _var_map.insert(std::make_pair(x, tx));
    add_rev_var_map(_rev_var_map, tx, x);
  }

  // Build the tree for a linexpr, and ensure that
  // values for the subterms are sustained.

  term_id_t term_of_const(const number_t &n) {
    dom_number dom_n(n);
    boost::optional<term_id_t> opt_n(_ttbl.find_const(dom_n));
    if (opt_n) {
      return *opt_n;
    } else {
      term_id_t term_n(_ttbl.make_const(dom_n));
      dom_var_t v = domvar_of_term(term_n);

      dom_linexp_t exp(n);
      _impl.assign(v, exp);
      return term_n;
    }
  }

  term_id_t term_of_var(variable_t v, var_map_t &var_map,
                        rev_var_map_t &rvar_map, ttbl_t &ttbl) {
    auto it(var_map.find(v));
    if (it != var_map.end()) {
      // assert((*it).first == v);
      assert(ttbl.size() > (*it).second);
      return (*it).second;
    } else {
      // Allocate a fresh term
      term_id_t id(ttbl.fresh_var());
      var_map[v] = id;
      add_rev_var_map(rvar_map, id, v);
      ttbl.add_ref(id);
      return id;
    }
  }

  term_id_t term_of_var(variable_t v) {
    return term_of_var(v, _var_map, _rev_var_map, _ttbl);
  }

  term_id_t term_of_linterm(linterm_t term) {
    if (term.first == 1) {
      return term_of_var(term.second);
    } else {
      return build_term(OP_MULTIPLICATION, term_of_const(term.first),
                        term_of_var(term.second));
    }
  }

  // OpTy = [arith_operation_t | bitwise_operation_t]
  template <typename OpTy>
  term_id_t build_term(OpTy op, term_id_t ty, term_id_t tz) {
    // Check if the term already exists
    term::term_operator_t binop = term::conv2termop(op);
    boost::optional<term_id_t> eopt(_ttbl.find_ftor(binop, ty, tz));
    if (eopt) {
      return *eopt;
    } else {
      // Create the term
      term_id_t tx = _ttbl.apply_ftor(binop, ty, tz);
      dom_var_t v(domvar_of_term(tx));
      dom_var_t y(domvar_of_term(ty));
      dom_var_t z(domvar_of_term(tz));

      // Set evaluation
      CRAB_LOG("term", crab::outs() << "Prev: " << _impl << "\n");
      _impl.apply(op, v, y, z);

      CRAB_LOG("term", crab::outs() << "Should have " << v << " := " << y << op
                                    << z << "\n";
               crab::outs() << _impl << "\n";);
      return tx;
    }
  }

  term_id_t build_function(term_id_t ty, term_id_t tz) {
    term::term_operator_t op = term::TERM_OP_FUNCTION;
    // Check if the term already exists
    boost::optional<term_id_t> eopt(_ttbl.find_ftor(op, ty, tz));
    if (eopt) {
      return *eopt;
    } else {
      // Create the term
      term_id_t tx = _ttbl.apply_ftor(op, ty, tz);
      return tx;
    }
  }

  term_id_t build_linexpr(const linear_expression_t &e) {
    number_t cst = e.constant();
    typename linear_expression_t::const_iterator it(e.begin());
    if (it == e.end())
      return term_of_const(cst);

    term_id_t t;
    if (cst == 0) {
      t = term_of_linterm(*it);
      ++it;
    } else {
      t = term_of_const(cst);
    }
    for (; it != e.end(); ++it) {
      t = build_term(OP_ADDITION, t, term_of_linterm(*it));
    }

    CRAB_LOG("term", crab::outs() << "Should have " << domvar_of_term(t)
                                  << " := " << e << "\n"
                                  << _impl << "\n");
    return t;
  }

  dom_var_t domvar_of_term(term_id_t id) {
    auto it = _term_map.find(id);
    if (it != _term_map.end()) {
      return (*it).second;
    } else {
      // Allocate a fresh variable
      dom_var_t dvar(_alloc.next());
      _term_map.insert(std::make_pair(id, dvar));
      return dvar;
    }
  }

  // used by to_linear_constraint_system()
  boost::optional<dom_var_t> domvar_of_term(term_id_t id) const {
    auto it = _term_map.find(id);
    if (it != _term_map.end()) {
      return (*it).second;
    } else {
      return boost::none;
    }
  }

  dom_var_t domvar_of_var(variable_t v) {
    return domvar_of_term(term_of_var(v));
  }

  // Remap a linear constraint to the domain.
  dom_linexp_t rename_linear_expr(linear_expression_t exp) {
    number_t cst(exp.constant());
    dom_linexp_t dom_exp(cst);
    for (auto const &v : exp.variables()) {
      dom_exp = dom_exp + exp[v] * domvar_of_var(v);
    }
    return dom_exp;
  }

  dom_lincst_t rename_linear_cst(linear_constraint_t cst) {
    return dom_lincst_t(rename_linear_expr(cst.expression()),
                        (typename dom_lincst_t::kind_t)cst.kind());
  }

  // Assumption: vars(exp) subseteq keys(map)
  // XXX JNL: exp can have variables that are not in rev_map (e.g.,
  // some generated by build_linexpr).
  boost::optional<linear_expression_t>
  rename_linear_expr_rev(dom_linexp_t exp, rev_map_t rev_map) const {
    number_t cst(exp.constant());
    linear_expression_t rev_exp(cst);
    for (auto const &v : exp.variables()) {
      auto it = rev_map.find(v);
      if (it != rev_map.end()) {
        variable_t v_out((*it).second);
        rev_exp = rev_exp + exp[v] * v_out;
      } else
        return boost::optional<linear_expression_t>();
    }
    return rev_exp;
  }

  boost::optional<linear_constraint_t>
  rename_linear_cst_rev(dom_lincst_t cst, rev_map_t rev_map) const {
    boost::optional<linear_expression_t> e =
        rename_linear_expr_rev(cst.expression(), rev_map);
    if (e)
      return linear_constraint_t(
          *e, (typename linear_constraint_t::kind_t)cst.kind());
    else
      return boost::optional<linear_constraint_t>();
  }

  boost::optional<std::pair<variable_t, variable_t>>
  get_eq_or_diseq(linear_constraint_t cst) {
    if (cst.is_equality() || cst.is_disequation()) {
      if (cst.size() == 2 && cst.constant() == 0) {
        auto it = cst.begin();
        auto nx = it->first;
        auto vx = it->second;
        ++it;
        assert(it != cst.end());
        auto ny = it->first;
        auto vy = it->second;
        if (nx == (ny * -1)) {
          return std::make_pair(vx, vy);
        }
      }
    }
    return boost::optional<std::pair<variable_t, variable_t>>();
  }

  struct WidenOp {
    dom_t apply(dom_t before, dom_t after) { return before || after; }
  };

  template <typename Thresholds> struct WidenWithThresholdsOp {
    const Thresholds &m_ts;

    WidenWithThresholdsOp(const Thresholds &ts) : m_ts(ts) {}

    dom_t apply(dom_t before, dom_t after) {
      return before.widening_thresholds(after, m_ts);
    }
  };

  template <typename WidenOp>
  term_domain_t widening(const term_domain_t &o, WidenOp widen_op) const {

    // The left operand of the widenning cannot be closed, otherwise
    // termination is not ensured. However, if the right operand is
    // close precision may be improved.
    term_domain_t right(o);
    right.normalize();
    if (is_bottom()) {
      return right;
    } else if (right.is_bottom()) {
      return *this;
    } else {
      term_domain_t left(*this);

      // First, we need to compute the new term table.
      ttbl_t out_tbl;
      // Mapping of (term, term) pairs to terms in the join state
      typename ttbl_t::gener_map_t gener_map;

      var_map_t out_vmap;
      rev_var_map_t out_rvmap;
      dom_var_alloc_t palloc(left._alloc, right._alloc);

      // For each program variable in state, compute a generalization
      for (auto p : left._var_map) {
        const variable_t &v = p.first;
        // term_id_t tx(term_of_var(v));
        term_id_t tx = p.second;
        term_id_t ty(right.term_of_var(v));

        term_id_t tz =
            left._ttbl.generalize(right._ttbl, tx, ty, out_tbl, gener_map);
        out_vmap[v] = tz;
        add_rev_var_map(out_rvmap, tz, v);
      }

      // Rename the common terms together
      dom_t x_impl(std::move(left._impl));
      dom_t y_impl(std::move(right._impl));

      // Perform the mapping
      term_map_t out_map;
      std::vector<dom_var_t> out_varnames;
      for (auto p : gener_map) {
        auto txy = p.first;
        term_id_t tz = p.second;
        dom_var_t vt(palloc.next());
        out_map.insert(std::make_pair(tz, vt));

        dom_var_t vx = left.domvar_of_term(txy.first);
        dom_var_t vy = right.domvar_of_term(txy.second);

        out_varnames.push_back(vt);

        x_impl.assign(vt, vx);
        y_impl.assign(vt, vy);
      }

      x_impl.project(out_varnames);
      y_impl.project(out_varnames);

      dom_t x_widen_y = widen_op.apply(x_impl, y_impl);

      for (auto p : out_vmap)
        out_tbl.add_ref(p.second);

      CRAB_LOG("term", crab::outs()
                           << "============ WIDENING ==================";
               crab::outs() << x_impl << "\n~~~~~~~~~~~~~~~~";
               crab::outs() << y_impl << "\n----------------";
               crab::outs() << x_widen_y << "\n================"
                            << "\n");

      term_domain_t res(std::move(palloc), std::move(out_vmap),
                        std::move(out_rvmap), std::move(out_tbl),
                        std::move(out_map), std::move(x_widen_y));
      return res;
    }
  }

  // Choose one non-var term from the equivalence class
  // associated with t.
  template <typename Range>
  boost::optional<term_id_t> choose_non_var(ttbl_t &ttbl,
                                            const Range &terms) const {
    std::vector<term_id_t> non_var_terms(terms.size());
    auto it = std::copy_if(terms.begin(), terms.end(), non_var_terms.begin(),
                           [&ttbl](term_id_t t) {
                             term_t *t_ptr = ttbl.get_term_ptr(t);
                             return (t_ptr && t_ptr->kind() == term::TERM_APP);
                           });
    non_var_terms.resize(std::distance(non_var_terms.begin(), it));
    if (non_var_terms.empty()) {
      return boost::optional<term_id_t>();
    } else {
      // TODO: the heuristics as described in the VMCAI'16 paper that
      // chooses the one that has more references each class.  For
      // now, we just make sure that we always pick the same term in a
      // deterministic way.
      std::sort(non_var_terms.begin(), non_var_terms.end(),
		[&ttbl](const term_id_t &t1, const term_id_t &t2) {
		  term_t *t1_ptr = ttbl.get_term_ptr(t1);
		  assert(t1_ptr);
		  term_t *t2_ptr = ttbl.get_term_ptr(t2);
		  assert(t2_ptr);
		  return *t1_ptr < *t2_ptr;
		});
      return *(non_var_terms.begin());
    }
  }

  term_id_t build_dag_term(ttbl_t &ttbl, int t,
                           term::congruence_closure_solver<ttbl_t> &solver,
                           ttbl_t &out_ttbl, std::vector<int> &stack,
                           std::map<int, term_id_t> &cache) const {

    // already processed
    auto it = cache.find(t);
    if (it != cache.end()) {
      CRAB_LOG("terms-meet", crab::outs() << "build_dag_term. Found in cache: ";
               crab::outs() << "t" << t << " --> "
                            << "t" << it->second << "\n";);
      return it->second;
    }

    // break the cycle with a fresh variable
    if (std::find(stack.begin(), stack.end(), t) != stack.end()) {
      term_id_t v = out_ttbl.fresh_var();
      CRAB_LOG("terms-meet", crab::outs() << "build_dag_term. Detected cycle: ";
               crab::outs() << "t" << t << " --> "
                            << "t" << v << "\n";);
      return v;
    }

    stack.push_back(t);
    auto membs = solver.get_members(t);
    boost::optional<term_id_t> f = choose_non_var(ttbl, membs);

    if (!f) {
      // no concrete definition exists return a fresh variable
      term_id_t v = out_ttbl.fresh_var();
      auto res = cache.insert(std::make_pair(t, v));
      stack.pop_back();
      CRAB_LOG("terms-meet", crab::outs()
                                 << "build_dag_term. No concrete definition: ";
               crab::outs() << "t" << t << " --> "
                            << "t" << (res.first)->second << "\n";);
      return (res.first)->second;
    } else {
      // traverse recursively the term
      term_t *f_ptr = ttbl.get_term_ptr(*f);
      CRAB_LOG("terms-meet",
               crab::outs()
                   << "build_dag_term. Traversing recursively the term "
                   << "t" << *f << ":";
               crab::outs() << *f_ptr << "\n";);
      const std::vector<term_id_t> &args(term::term_args(f_ptr));
      std::vector<term_id_t> res_args;
      res_args.reserve(args.size());
      for (term_id_t c : args) {
        res_args.push_back(build_dag_term(ttbl, solver.get_class(c), solver,
                                          out_ttbl, stack, cache));
      }
      auto res = cache.insert(std::make_pair(
          t, out_ttbl.apply_ftor(term::term_ftor(f_ptr), res_args)));
      stack.pop_back();
      CRAB_LOG("terms-meet", crab::outs()
                                 << "build_dag_term. Finished recursive case: ";
               crab::outs() << "t" << t << " --> "
                            << "t" << (res.first)->second << "\n";);
      return (res.first)->second;
    }
  }

public:
  term_domain_t make_top() const override { return term_domain_t(true); }

  term_domain_t make_bottom() const override { return term_domain_t(false); }

  void set_to_top() override {
    term_domain abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    term_domain abs(false);
    std::swap(*this, abs);
  }

  term_domain() : _is_bottom(false) {}

  term_domain(const term_domain_t &o)
      : _is_bottom(o._is_bottom), _ttbl(o._ttbl), _impl(o._impl),
        _alloc(o._alloc), _var_map(o._var_map), _rev_var_map(o._rev_var_map),
        _term_map(o._term_map), changed_terms(o.changed_terms) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    check_terms(__LINE__);
  }

  term_domain_t &operator=(const term_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");

    o.check_terms(__LINE__);
    if (this != &o) {
      _is_bottom = o._is_bottom;
      _ttbl = o._ttbl;
      _impl = o._impl;
      _alloc = o._alloc;
      _var_map = o._var_map;
      _rev_var_map = o._rev_var_map;
      _term_map = o._term_map;
      changed_terms = o.changed_terms;
    }
    check_terms(__LINE__);
    return *this;
  }

  bool is_bottom() const override { return _is_bottom; }

  bool is_top() const override { return !_var_map.size() && !is_bottom(); }

  bool is_normalized() {
    return changed_terms.size() == 0;
    // return _is_normalized(_impl);
  }

  // Lattice operations
  bool operator<=(const term_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    // Require normalization of the first argument
    term_domain_t left(*this);
    left.normalize();

    if (left.is_bottom()) {
      return true;
    } else if (o.is_bottom()) {
      return false;
    } else {
      term_domain_t right(o);
      typename ttbl_t::term_map_t gen_map;
      dom_var_alloc_t palloc(left._alloc, right._alloc);

      // Build up the mapping of right onto left, variable by variable.
      // Assumption: the set of variables in left & right are common.
      for (auto p : left._var_map) {
        if (!left._ttbl.map_leq(right._ttbl, left.term_of_var(p.first),
                                right.term_of_var(p.first), gen_map))
          return false;
      }
      // We now have a mapping of reachable y-terms to x-terms.
      // Create copies of left._impl and right._impl with a common
      // variable set.
      dom_t x_impl(std::move(left._impl));
      dom_t y_impl(std::move(right._impl));

      // Perform the mapping
      std::vector<dom_var_t> out_varnames;
      for (auto p : gen_map) {
        // dom_var_t vt = _alloc.next();
        dom_var_t vt(palloc.next());
        dom_var_t vx = left.domvar_of_term(p.second);
        dom_var_t vy = right.domvar_of_term(p.first);

        out_varnames.push_back(vt);

        x_impl.assign(vt, vx);
        y_impl.assign(vt, vy);
      }

      x_impl.project(out_varnames);
      y_impl.project(out_varnames);

      return x_impl <= y_impl;
    }
  }

  void operator|=(const term_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    // Requires normalization of both operands
    term_domain_t right(o);
    normalize();
    right.normalize();

    if (is_bottom() || right.is_top()) {
      *this = right;
    } else if (right.is_bottom() || is_top()) {
      return;
    } else {
      // First, we need to compute the new term table.
      ttbl_t out_tbl;

      // Mapping of (term, term) pairs to terms in the join state
      typename ttbl_t::gener_map_t gener_map;

      var_map_t out_vmap;
      rev_var_map_t out_rvmap;
      dom_var_alloc_t palloc(_alloc, right._alloc);

      // For each program variable in state, compute a generalization
      for (auto p : _var_map) {
        const variable_t &v = p.first;
        // term_id_t tx(term_of_var(v));
        term_id_t tx = p.second;
        term_id_t ty(right.term_of_var(v));

        term_id_t tz =
            _ttbl.generalize(right._ttbl, tx, ty, out_tbl, gener_map);
        assert(tz < out_tbl.size());
        out_vmap[v] = tz;
        add_rev_var_map(out_rvmap, tz, v);
      }

      // Rename the common terms together
      // Perform the mapping
      term_map_t out_map;
      std::vector<dom_var_t> out_varnames;
      for (auto p : gener_map) {
        auto txy = p.first;
        term_id_t tz = p.second;
        dom_var_t vt(palloc.next());
        out_map.insert(std::make_pair(tz, vt));

        dom_var_t vx = domvar_of_term(txy.first);
        dom_var_t vy = right.domvar_of_term(txy.second);

        out_varnames.push_back(vt);

        _impl.assign(vt, vx);
        right._impl.assign(vt, vy);
      }

      _impl.project(out_varnames);
      right._impl.project(out_varnames);

      _impl |= right._impl;

      for (auto p : out_vmap)
        out_tbl.add_ref(p.second);

      _alloc = std::move(palloc);
      _var_map = std::move(out_vmap);
      _rev_var_map = std::move(out_rvmap);
      _ttbl = std::move(out_tbl);
      _term_map = std::move(out_map);
      _is_bottom = (_impl.is_bottom() ? true : false);
    }
  }

  term_domain_t operator|(const term_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    term_domain_t left(*this);
    term_domain_t right(o);

    // Requires normalization of both operands
    left.normalize();
    right.normalize();

    if (left.is_bottom() || right.is_top()) {
      return right;
    } else if (right.is_bottom() || left.is_top()) {
      return left;
    } else {

      // First, we need to compute the new term table.
      ttbl_t out_tbl;
      // Mapping of (term, term) pairs to terms in the join state
      typename ttbl_t::gener_map_t gener_map;

      var_map_t out_vmap;
      rev_var_map_t out_rvmap;
      dom_var_alloc_t palloc(left._alloc, right._alloc);

      // For each program variable in state, compute a generalization
      for (auto p : left._var_map) {
        const variable_t &v = p.first;
        // term_id_t tx(term_of_var(v));
        term_id_t tx = p.second;
        term_id_t ty = right.term_of_var(v);

        term_id_t tz =
            left._ttbl.generalize(right._ttbl, tx, ty, out_tbl, gener_map);
        assert(tz < out_tbl.size());
        out_vmap[v] = tz;
        add_rev_var_map(out_rvmap, tz, v);
      }

      // Rename the common terms together
      dom_t x_impl(std::move(left._impl));
      dom_t y_impl(std::move(right._impl));

      // Perform the mapping
      term_map_t out_map;
      std::vector<dom_var_t> out_varnames;
      for (auto p : gener_map) {
        auto txy = p.first;
        term_id_t tz = p.second;
        dom_var_t vt(palloc.next());
        out_map.insert(std::make_pair(tz, vt));

        dom_var_t vx = left.domvar_of_term(txy.first);
        dom_var_t vy = right.domvar_of_term(txy.second);

        out_varnames.push_back(vt);

        x_impl.assign(vt, vx);
        y_impl.assign(vt, vy);
      }

      CRAB_LOG("term", crab::outs() << "============ JOIN =================="
                                    << *this << "\n"
                                    << "~~~~~~~~~~~~~~~~" << o << "\n"
                                    << "----------------"
                                    << "ren_0(x) = " << x_impl
                                    << "ren_0(y) = " << y_impl << "\n");

      x_impl.project(out_varnames);
      y_impl.project(out_varnames);

      dom_t x_join_y = x_impl | y_impl;

      for (auto p : out_vmap)
        out_tbl.add_ref(p.second);

      term_domain_t res(std::move(palloc), std::move(out_vmap),
                        std::move(out_rvmap), std::move(out_tbl),
                        std::move(out_map), std::move(x_join_y));

      CRAB_LOG("term", crab::outs() << "After elimination:\n" << res << "\n");

      return res;
    }
  }

  term_domain_t operator||(const term_domain_t &other) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    WidenOp op;
    return widening(other, op);
  }

  term_domain_t widening_thresholds(
      const term_domain_t &other,
      const thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    WidenWithThresholdsOp<thresholds<number_t>> op(ts);
    return widening(other, op);
  }

  // Meet
  term_domain_t operator&(const term_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    // Does not require normalization of any of the two operands
    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      ttbl_t out_ttbl(_ttbl);
      std::map<term_id_t, term_id_t> copy_map;
      // bring all terms to one ttbl
      for (auto p : o._var_map) {
        // const variable_t &v = p.first;
        // term_id_t tx(o.term_of_var(v));
        term_id_t tx = p.second;
        out_ttbl.copy_term(o._ttbl, tx, copy_map);
      }

      // build unifications between terms from this and o
      std::vector<std::pair<term_id_t, term_id_t>> eqs;
      for (auto p : _var_map) {
        variable_t v(p.first);
        auto it = o._var_map.find(v);
        if (it != o._var_map.end()) {
          // term_id_t tx(term_of_var(v));
          term_id_t tx = p.second;
          eqs.push_back(std::make_pair(tx, copy_map[it->second]));
        }
      }

      // compute equivalence classes
      term::congruence_closure_solver<ttbl_t> solver(out_ttbl);
      solver.run(eqs);

      std::vector<int> stack;
      std::map<int, term_id_t> cache;
      var_map_t out_vmap;
      rev_var_map_t out_rvmap;
      // new map from variable to an acyclic term
      for (auto p : _var_map) {
        const variable_t &v = p.first;
        // term_id_t t_old(term_of_var(v));
        term_id_t t_old = p.second;
        term_id_t t_new = build_dag_term(out_ttbl, solver.get_class(t_old),
                                         solver, out_ttbl, stack, cache);
        out_vmap[v] = t_new;
        add_rev_var_map(out_rvmap, t_new, v);
      }
      for (auto p : o._var_map) {
        variable_t v(p.first);
        if (out_vmap.find(v) != out_vmap.end())
          continue;
        term_id_t t_old(copy_map[p.second /*o.term_of_var(v)*/]);
        term_id_t t_new = build_dag_term(out_ttbl, solver.get_class(t_old),
                                         solver, out_ttbl, stack, cache);
        out_vmap[v] = t_new;
        add_rev_var_map(out_rvmap, t_new, v);
      }

      for (auto p : out_vmap)
        out_ttbl.add_ref(p.second);

      // Rename the base domains
      dom_var_alloc_t palloc(_alloc, o._alloc);
      dom_t x_impl(_impl);
      dom_t y_impl(o._impl);
      term_map_t out_map; // map term to dom var
      std::vector<dom_var_t> out_varnames;
      for (auto p : out_vmap) {
        const variable_t &v = p.first;
        term_id_t t_new = p.second;
        dom_var_t vt(palloc.next());
        out_map.insert(std::make_pair(t_new, vt));
        // renaming this's base domain
        auto xit = _var_map.find(v);
        if (xit != _var_map.end()) {
          if (boost::optional<dom_var_t> vx = domvar_of_term(xit->second)) {
            x_impl.assign(vt, *vx);
          }
        }
        // renaming o's base domain
        auto yit = o._var_map.find(v);
        if (yit != o._var_map.end()) {
          if (boost::optional<dom_var_t> vy = o.domvar_of_term(yit->second)) {
            y_impl.assign(vt, *vy);
          }
        }
        out_varnames.push_back(vt);
      }

      x_impl.project(out_varnames);
      y_impl.project(out_varnames);

      dom_t x_meet_y = x_impl & y_impl;

      term_domain_t res(std::move(palloc), std::move(out_vmap),
                        std::move(out_rvmap), std::move(out_ttbl),
                        std::move(out_map), std::move(x_meet_y));

      CRAB_LOG("term", crab::outs() << "============ MEET ==================";
               crab::outs() << *this << "\n----------------";
               crab::outs() << o << "\n----------------";
               crab::outs() << res << "\n================"
                            << "\n");
      return res;
    }
  }

  void operator&=(const term_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    // TODO: improve by avoiding the copy of the left operand
    *this = *this & o;
  }
  
  // Narrowing
  term_domain_t operator&&(const term_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    CRAB_WARN("Term narrowing operator replaced with meet");
    return *this & o;
  }

  // Remove a variable from the scope
  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    auto it(_var_map.find(v));
    if (it != _var_map.end()) {
      term_id_t t = (*it).second;
      _var_map.erase(it);
      remove_rev_var_map(t, v);
      deref(t);
    }
    CRAB_LOG("term", crab::outs()
                         << "After removing " << v << ": " << *this << "\n";);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (this->is_bottom()) {
      return;
    } else {
      // dom_linexp_t dom_e(rename_linear_expr(e));
      term_id_t tx(build_linexpr(e));
      rebind_var(x, tx);

      check_terms(__LINE__);

      CRAB_LOG("term", crab::outs() << "*** Assign " << x << ":=" << e << ":"
                                    << *this << "\n");
      return;
    }
  }

  // Apply operations to variables.

  // x = y op z
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    check_terms(__LINE__);
    if (this->is_bottom()) {
      return;
    } else {
      term_id_t tx(build_term(op, term_of_var(y), term_of_var(z)));
      rebind_var(x, tx);
    }
    check_terms(__LINE__);
    CRAB_LOG("term", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << z << ":" << *this << "\n");
  }

  // x = y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (this->is_bottom()) {
      return;
    } else {
      term_id_t tx(build_term(op, term_of_var(y), term_of_const(k)));
      rebind_var(x, tx);
    }
    check_terms(__LINE__);
    CRAB_LOG("term", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << k << ":" << *this << "\n");
    return;
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const term_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    crab::domains::BackwardAssignOps<term_domain_t>::assign(*this, x, e, inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const term_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<term_domain_t>::apply(*this, op, x, y, z,
                                                           inv);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const term_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    crab::domains::BackwardAssignOps<term_domain_t>::apply(*this, op, x, y, z,
                                                           inv);
  }

  void operator+=(const linear_constraint_t &cst) {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    CRAB_LOG("term", crab::outs() << "*** Before assume " << cst << ":" << *this
                                  << "\n");

    using pair_var_t = std::pair<variable_t, variable_t>;

    if (boost::optional<pair_var_t> eq = get_eq_or_diseq(cst)) {
      term_id_t tx(term_of_var((*eq).first));
      term_id_t ty(term_of_var((*eq).second));
      if (cst.is_disequation()) {
        if (tx == ty) {
          set_to_bottom();
          CRAB_LOG("term", crab::outs() << "*** After assume " << cst << ":"
                                        << *this << "\n");
          return;
        }
      } else {
        // not bother if they are already equal
        if (tx == ty)
          return;

        // congruence closure to compute equivalence classes
        term::congruence_closure_solver<ttbl_t> solver(_ttbl);
        std::vector<std::pair<term_id_t, term_id_t>> eqs = {
            std::make_pair(tx, ty)};
        solver.run(eqs);

        std::vector<int> stack;
        std::map<int, term_id_t> cache;
        dom_t x_impl(_impl);
        std::vector<dom_var_t> out_varnames;
        // new map from variable to an acyclic term
        // and also renaming of the base domain
        for (auto p : _var_map) {
          const variable_t &v = p.first;
          term_id_t t_old(term_of_var(v));
          term_id_t t_new = build_dag_term(_ttbl, solver.get_class(t_old),
                                           solver, _ttbl, stack, cache);

          dom_var_t vt = domvar_of_term(t_new);
          dom_var_t vx = domvar_of_term(t_old);

          out_varnames.push_back(vt);
          x_impl.assign(vt, vx);

          rebind_var(v, t_new);
        }
        x_impl.project(out_varnames);
        std::swap(_impl, x_impl);
      }
    }

    dom_lincst_t cst_rn(rename_linear_cst(cst));
    _impl += cst_rn;
    // Possibly tightened some variable in cst
    for (auto const &v : cst.variables()) {
      CRAB_LOG("term-normalization",
               crab::outs() << "Added to the normalization queue "
                            << "t" << term_of_var(v) << "["
                            << domvar_of_term(term_of_var(v)) << "] "
                            << "from variable " << v << "\n";);
      changed_terms.insert(term_of_var(v));
    }

    // Probably doesn't need to done so eagerly.
    normalize();

    CRAB_LOG("term", crab::outs()
                         << "*** After assume " << cst << ":" << *this << "\n");
    return;
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    for (auto cst : csts) {
      this->operator+=(cst);
    }
  }

  DEFAULT_ENTAILS(term_domain_t)
  
  /*
  // If the children of t have changed, see if re-applying
  // the definition of t tightens the domain.
  void tighten_term(term_id_t t)
  {
  dom_t tight = _impl&eval_ftor_copy(_impl, _ttbl, t);

  if(!(_impl <= tight))
  {
  // Applying the functor has changed the domain
  _impl = tight;
  for(term_id_t p : _ttbl.parents(t))
  tighten_term(p);
  }
  check_terms(__LINE__);
  }
  */

  interval_t operator[](const variable_t &x) override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");

    // Needed for accuracy
    normalize();

    if (is_bottom()) {
      return interval_t::bottom();
    } else {
      auto it = _var_map.find(x);
      if (it == _var_map.end()) {
	return interval_t::top();
      } else {
	dom_var_t dom_x = domvar_of_term(it->second);
	return _impl[dom_x];
      }
    }
  }

  interval_t at(const variable_t &x) const override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");

    if (is_bottom()) {
      return interval_t::bottom();
    } else {
      auto it = _var_map.find(x);
      if (it == _var_map.end()) {
	return interval_t::top();
      } else {
	if (boost::optional<dom_var_t> dom_x = domvar_of_term(it->second)) {
	  return _impl.at(*dom_x);
	} else {
	  return interval_t::top();
	}
      }
    }
  }  

  void set(const variable_t &x, interval_t intv) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    rebind_var(x, term_of_itv(intv.lb(), intv.ub()));
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    int_cast_domain_traits<term_domain_t>::apply(*this, op, dst, src);    
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (this->is_bottom()) {
      return;
    } else {
      term_id_t tx(build_term(op, term_of_var(y), term_of_var(z)));
      rebind_var(x, tx);
    }
    check_terms(__LINE__);
    CRAB_LOG("term", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << z << ":" << *this << "\n");
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (this->is_bottom()) {
      return;
    } else {
      term_id_t tx(build_term(op, term_of_var(y), term_of_const(k)));
      rebind_var(x, tx);
    }

    check_terms(__LINE__);
    CRAB_LOG("term", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << k << ":" << *this << "\n");
    return;
  }

  /* Array operations */

  virtual void array_init(const variable_t & /*a*/,
                          const linear_expression_t & /*elem_size*/,
                          const linear_expression_t & /*lb_idx*/,
                          const linear_expression_t & /*ub_idx*/,
                          const linear_expression_t & /*val*/) override {
    // TODO: perform a loop of array stores if [lb_idx, ub_idx]
    //       is finite.
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t & /*elem_size*/,
                          const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.array_read");
    crab::ScopedCrabStats __st__(domain_name() + ".array_read");

    if (this->is_bottom()) {
      return;
    } else {
      /**
       *  We treat the array load as an uninterpreted function
       *  lhs := array_load(a, i) -->  lhs := f(a,i)
       */
      term_id_t t_uf(build_function(term_of_var(a), build_linexpr(i)));
      rebind_var(lhs, t_uf);
    }
    check_terms(__LINE__);
    CRAB_LOG("term", crab::outs() << lhs << ":=" << a << "[" << i << "]  -- "
                                  << *this << "\n";);
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t & /*elem_size*/,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool /*is_strong_update*/) override {
    crab::CrabStats::count(domain_name() + ".count.array_store");
    crab::ScopedCrabStats __st__(domain_name() + ".array_store");

    if (this->is_bottom()) {
      return;
    } else {
      auto &vfac = const_cast<varname_t *>(&(a.name()))->get_var_factory();
      /**
       *  We treat the array store as an uninterpreted function
       *  array_store(a, i, val) -->  tmp := f(a,i); assume(tmp == val);
       */
       /// -- tmp := f(a,i)
      term_id_t t_uf(build_function(term_of_var(a), build_linexpr(i)));
      // FIXME:
      // 1. we are creating a fresh varname each time we are called.
      // 2. we are creating an untyped variable. This is ok if the
      // underlying domain is, for instance, intervals but it will
      // fail if we use a more type-dependent numerical domain such as
      // wrapped intervals.
      variable_t tmp(vfac.get());
      // forget tmp
      this->operator-=(tmp);
      // forget the old value for t_uf, otherwise we can get
      // incorrectly bottom when we add the constraint val == tmp.
      _impl -= domvar_of_term(t_uf);
      rebind_var(tmp, t_uf);
      /// -- assume(tmp == val)
      this->operator+=(val == tmp);
    }
    check_terms(__LINE__);
    CRAB_LOG("term", crab::outs() << a << "[" << i << "]:=" << val << " -- "
                                  << *this << "\n";);
  }

  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &i,
                                 const linear_expression_t &j,
                                 const linear_expression_t &v) override {
    // do nothing
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {
    // do nothing
  }

  // backward array operations
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const term_domain_t &invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const term_domain_t &invariant) override {
    *this -= lhs;
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const term_domain_t &invariant) override {}
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const term_domain_t &invariant) override {}
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const term_domain_t &invariant) override {}

  /// term_domain implements standard abstract operations of a
  /// numerical domain and some array operations. It is intended to be
  /// used as a leaf domain in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(term_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(term_domain_t)
  DEFAULT_SELECT(term_domain_t)
  DEFAULT_WEAK_ASSIGN(term_domain_t)

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_top() || is_bottom())
      return;

    CRAB_LOG("term", crab::outs() << "Renaming {"; for (auto v
                                                        : from) crab::outs()
                                                   << v << ";";
             crab::outs() << "} with "; for (auto v
                                             : to) crab::outs()
                                        << v << ";";
             crab::outs() << "}:\n"; crab::outs() << *this << "\n";);

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      const variable_t &v = from[i];
      const variable_t &new_v = to[i];
      if (v == new_v) { // nothing to rename
        continue;
      }

      {
        auto it = _var_map.find(new_v);
        if (it != _var_map.end()) {
          CRAB_ERROR(domain_name() + "::rename assumes that ", new_v,
                     " does not exist");
        }
      }

      auto it = _var_map.find(v);
      if (it != _var_map.end()) {
        term_id_t id = it->second;
        _var_map.erase(it);
        _var_map.insert(std::make_pair(new_v, id));
        remove_rev_var_map(id, v);
        add_rev_var_map(_rev_var_map, id, new_v);
      }
    }

    CRAB_LOG("term", crab::outs() << "RESULT=" << *this << "\n");
  }

  // extract operation is used during reduction with other domains.
  void extract(const variable_t &x, linear_constraint_system_t &csts,
               bool only_equalities /*unused*/) {
    crab::CrabStats::count(domain_name() + ".count.extract");
    crab::ScopedCrabStats __st__(domain_name() + ".extract");

    if (!is_normalized())
      normalize();

    if (is_bottom()) {
      return;
    }

    // TODO: make user parameter
    const unsigned max_eq = 10;
    // We limit the number of equalities via
    // constants. Otherwise, the number of equalities can be too
    // large (e.g., with domains like array_expansion) and it
    // would make the reduction very slow.

    // Extract equalities
    auto it = _var_map.find(x);
    if (it != _var_map.end()) {
      term_id_t tx = it->second;
      auto tx_ptr = _ttbl.get_term_ptr(tx);

      bool active_threshold = false;
      if (tx_ptr->kind() == term::TERM_CONST) {
        active_threshold = true;
      }
      auto &varset = _rev_var_map[tx];
      unsigned num_eq = 0;
      for (auto var : varset) {
        if (active_threshold && num_eq > max_eq) {
          return;
        }
        if (var.index() != x.index()) {
          num_eq++;
          linear_constraint_t cst(linear_expression_t(x) ==
                                  linear_expression_t(var));
          CRAB_LOG("terms", crab::outs() << "Extracting " << cst << "\n";);
          csts += cst;
        }
      }
    }
  }

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top())
      return;

    for (auto v : variables) {
      *this -= v;
    }
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
    for (auto p : _var_map)
      s1.insert(p.first);
    s2.insert(variables.begin(), variables.end());
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(s3));
    forget(s3);
  }

  void expand(const variable_t &x, const variable_t &y) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top()) {
      return;
    }

    linear_expression_t e(x);
    term_id_t tx(build_linexpr(e));
    rebind_var(y, tx);

    check_terms(__LINE__);
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const term_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  // Propagate information from tightened terms to
  // parents/children.
  void normalize() override {
    crab::CrabStats::count(domain_name() + ".count.normalize");
    crab::ScopedCrabStats __st__(domain_name() + ".normalize");
    TermNormalizer<Info, typename Info::domain_t>::normalize(*this);
  }

  void minimize() override {}

  /// Simplify the term associated with x by given the standard
  /// arithmetic meaning to the functors
  bool simplify(const variable_t &x) {
    crab::CrabStats::count(domain_name() + ".count.simplify");
    crab::ScopedCrabStats __st__(domain_name() + ".simplify");

    auto it = _var_map.find(x);
    if (it != _var_map.end()) {
      term_id_t t = it->second;
      simplifier_t simp(_ttbl);
      auto nt = simp.simplify_term(t);
      if (nt) {
        rebind_var(x, *nt);
        return true;
      }
    }
    return false;
  }

  // Output function
  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    // XXX: we can avoid the copy as we do in
    // to_linear_constraint_system but we don't bother for
    // pretty-printing.
    term_domain_t tmp(*this);

    // Normalization is not enforced in order to maintain accuracy
    // but we force it to display all the relationships.
    tmp.normalize();

    if (tmp.is_bottom()) {
      o << "_|_";
      return;
    }
    if (tmp._var_map.empty()) {
      o << "{}";
      return;
    }

    bool first = true;
    o << "{";
    for (auto p : tmp._var_map) {
      if (first)
        first = false;
      else
        o << ", ";
      o << p.first << " -> t" << p.second << "[" << tmp.domvar_of_term(p.second)
        << "]";
    }
    o << "}";

    // print underlying domain
    o << tmp._impl;

    CRAB_LOG("terms-print-ttbl",
             /// For debugging purposes
             o << " ttbl={" << tmp._ttbl << "}\n";);
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    // Collect the visible terms
    rev_map_t rev_map;
    std::vector<std::pair<variable_t, variable_t>> equivs;
    for (auto p : _var_map) {
      boost::optional<dom_var_t> dv_opt = domvar_of_term(p.second);
      if (!dv_opt)
        continue;

      dom_var_t dv = *dv_opt;
      auto it = rev_map.find(dv);
      if (it == rev_map.end()) {
        // The term has not yet been seen.
        rev_map.insert(std::make_pair(dv, p.first));
      } else {
        // The term is already mapped to (*it).second,
        // so add an equivalence.
        equivs.push_back(std::make_pair((*it).second, p.first));
      }
    }

    // Create a copy of _impl with only visible variables.
    dom_t d_vis(_impl);
    for (auto p : _term_map) {
      dom_var_t dv = p.second;
      if (rev_map.find(dv) == rev_map.end())
        d_vis -= dv;
    }

    // Now build and rename the constraint system, plus equivalences.
    dom_linsys_t dom_sys(d_vis.to_linear_constraint_system());

    linear_constraint_system_t out_sys;
    for (dom_lincst_t cst : dom_sys) {
      auto out_cst = rename_linear_cst_rev(cst, rev_map);
      if (!out_cst)
        continue;
      out_sys += *out_cst;
    }

    for (auto p : equivs) {
      CRAB_LOG("term", crab::outs() << "Added equivalence " << p.first << "="
                                    << p.second << "\n");
      out_sys += (p.first - p.second == 0);
    }

    // Now rename it back into the external scope.
    return out_sys;
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

  std::string domain_name() const override {
    std::string name("Term(" + _impl.domain_name() + ")");
    return name;
  }

}; // class term_domain

// Propagate information from tightened terms to
// parents/children.
template <class Info, class Abs> class TermNormalizer {
public:
  using term_domain_t = typename term_domain<Info>::term_domain_t;
  using term_id_t = typename term_domain_t::term_id_t;
  using dom_t = Abs;
  using ttbl_t = typename term_domain_t::ttbl_t;
  using term_t = typename ttbl_t::term_t;
  using term_set_t = boost::container::flat_set<term_id_t>;

  static void queue_push(ttbl_t &tbl,
                         std::vector<std::vector<term_id_t>> &queue,
                         term_id_t t) {
    int d = tbl.depth(t);
    if (d == 0) {
      // if depth 0 then t is a free variable: nothing to propagate.
      return;
    }
    while (queue.size() <= d) {
      queue.push_back(std::vector<term_id_t>());
    }
    queue[d].push_back(t);
  }

  static void normalize(term_domain_t &abs) {
    // First propagate down, then up.
    std::vector<std::vector<term_id_t>> queue;

    ttbl_t &ttbl(abs._ttbl);
    dom_t &impl = abs._impl;

    for (term_id_t t : abs.changed_terms) {
      queue_push(ttbl, queue, t);
    }

    dom_t d_prime = impl;
    // Propagate information to children.
    // Don't need to propagate level 0, since it's for free variables
    for (int d = queue.size() - 1; d > 0; d--) {
      for (term_id_t t : queue[d]) {
        term_t *t_ptr = ttbl.get_term_ptr(t);
        if (t_ptr->kind() != term::TERM_APP)
          continue;

        CRAB_LOG("term-normalization", crab::outs() << "Propagate to children ";
                 t_ptr->write(crab::outs()); crab::outs() << "\n";
                 crab::outs() << "\tTerm table: "; ttbl.write(crab::outs());
                 crab::outs() << "\n";);

        abs.eval_ftor_down(d_prime, ttbl, t);
        if (!(abs._impl <= d_prime)) {
          impl = d_prime;

          CRAB_LOG(
              "term-normalization",
              crab::outs() << "\trefinement done: enqueue children.\n";
              if (impl.is_bottom()) { crab::outs() << "\tfound bottom\n"; });

          // Enqueue the args.
          term_t *t_ptr = ttbl.get_term_ptr(t);
          const std::vector<term_id_t> &args(term::term_args(t_ptr));
          for (term_id_t c : args) {
            if (abs.changed_terms.find(c) == abs.changed_terms.end()) {
              abs.changed_terms.insert(c);
              queue[ttbl.depth(c)].push_back(c);
            }
          }
        } else {
          CRAB_LOG("term-normalization", crab::outs()
                                             << "\tno refinement done.\n";);
        }
      }
    }

    // Collect the parents of changed terms.
    term_set_t up_terms;
    std::vector<std::vector<term_id_t>> up_queue;
    for (term_id_t t : abs.changed_terms) {
      for (term_id_t p : ttbl.parents(t)) {
        if (up_terms.find(p) == up_terms.end()) {
          up_terms.insert(p);
          CRAB_LOG("term-normalization",
                   crab::outs()
                       << "t" << p << "[" << abs.domvar_of_term(p) << "]"
                       << " is a parent of "
                       << "t" << t << "[" << abs.domvar_of_term(t) << "]\n";);
          queue_push(ttbl, up_queue, p);
        }
      }
    }

    // Now propagate up, level by level.
    // This may miss inferences; for example with [[x = y - z]]
    // information about y can propagate to z.
    assert(up_queue.size() == 0 || up_queue[0].size() == 0);
    for (int d = 1; d < up_queue.size(); d++) {
      // up_queue[d] shouldn't change.
      for (term_id_t t : up_queue[d]) {
        abs.eval_ftor(d_prime, ttbl, t);
        CRAB_LOG("term-normalization", term_t *t_ptr = ttbl.get_term_ptr(t);
                 crab::outs() << "Propagate to parent ";
                 t_ptr->write(crab::outs()); crab::outs() << "\n";
                 crab::outs() << "\tTerm table: "; ttbl.write(crab::outs());
                 crab::outs() << "\n";);

        if (!(impl <= d_prime)) {
          CRAB_LOG("term-normalization",
                   crab::outs() << "Before up propagation: " << impl << "\n";
                   crab::outs()
                   << "After up propagation : " << d_prime << "\n";);

          // We need to do a meet here, as
          // impl and F(stmt)impl may be
          // incomparable
          impl = impl & d_prime;
          // impl = d_prime; // Old code

          CRAB_LOG(
              "term-normalization",
              crab::outs() << "\trefinement done: enqueue parents.\n";
              if (impl.is_bottom()) { crab::outs() << "\tfound bottom\n"; });

          for (term_id_t p : ttbl.parents(t)) {
            if (up_terms.find(p) == up_terms.end()) {
              up_terms.insert(p);
              queue_push(ttbl, up_queue, p);
            }
          }
        } else {
          CRAB_LOG("term-normalization", crab::outs()
                                             << "\tno refinement done.\n";);
        }
      }
    }

    abs.changed_terms.clear();

    if (abs._impl.is_bottom())
      abs.set_to_bottom();
  }
};

// Specialized implementation for interval domain.
// GKG: Should modify to work with any independent attribute domain.
#ifdef USE_TERM_INTERVAL_NORMALIZER
template <class Info, class Num, class Var>
class TermNormalizer<Info, ikos::interval_domain<Num, Var>> {
public:
  using term_domain_t = typename term_domain<Info>::term_domain_t;
  using term_id_t = typename term_domain_t::term_id_t;
  using dom_t = ikos::interval_domain<Num, Var>;
  using var_t = typename term_domain_t::dom_var_t;
  using ttbl_t = typename term_domain_t::ttbl_t;
  using term_t = typename ttbl_t::term_t;
  using term_set_t = boost::container::flat_set<term_id_t>;

  using interval_t = typename dom_t::interval_t;

  static void queue_push(ttbl_t &tbl,
                         std::vector<std::vector<term_id_t>> &queue,
                         term_id_t t) {
    int d = tbl.depth(t);
    if (d == 0) {
      // if depth 0 then t is a free variable: nothing to
      // propagate.
      return;
    }
    while (queue.size() <= d) {
      queue.push_back(std::vector<term_id_t>());
    }
    queue[d].push_back(t);
  }

  static void normalize(term_domain_t &abs) {
    // First propagate down, then up.
    std::vector<std::vector<term_id_t>> queue;

    ttbl_t &ttbl(abs._ttbl);
    dom_t &impl = abs._impl;
    if (impl.is_bottom()) {
      abs.set_to_bottom();
      return;
    }

    for (term_id_t t : abs.changed_terms) {
      queue_push(ttbl, queue, t);
    }

    // Propagate information to children.
    // Don't need to propagate level 0, since it's for free variables
    for (int d = queue.size() - 1; d > 0; d--) {
      for (term_id_t t : queue[d]) {
        term_t *t_ptr = ttbl.get_term_ptr(t);
        if (t_ptr->kind() != term::TERM_APP)
          continue;

        const std::vector<term_id_t> &args(term::term_args(t_ptr));
        std::vector<interval_t> arg_intervals;
        for (term_id_t c : args)
          arg_intervals.push_back(impl[abs.domvar_of_term(c)]);
        abs.eval_ftor_down(impl, ttbl, t);

        // Enqueue the args
        for (size_t ci = 0; ci < args.size(); ci++) {
          term_id_t c(args[ci]);
          var_t v = abs.domvar_of_term(c);
          interval_t v_upd(impl[v]);
          if (!(arg_intervals[ci] <= v_upd)) {
            impl.set(v, arg_intervals[ci] & v_upd);
            if (abs.changed_terms.find(c) == abs.changed_terms.end()) {
              abs.changed_terms.insert(c);
              queue[ttbl.depth(c)].push_back(c);
            }
          }
        }
      }
    }

    // Collect the parents of changed terms.
    term_set_t up_terms;
    std::vector<std::vector<term_id_t>> up_queue;
    for (term_id_t t : abs.changed_terms) {
      for (term_id_t p : ttbl.parents(t)) {
        if (up_terms.find(p) == up_terms.end()) {
          up_terms.insert(p);
          queue_push(ttbl, up_queue, p);
        }
      }
    }

    // Now propagate up, level by level.
    assert(up_queue.size() == 0 || up_queue[0].size() == 0);
    for (int d = 1; d < up_queue.size(); d++) {
      // up_queue[d] shouldn't change.
      for (term_id_t t : up_queue[d]) {
        var_t v = abs.domvar_of_term(t);
        interval_t v_interval = impl[v];

        abs.eval_ftor(impl, ttbl, t);

        interval_t v_upd = impl[v];
        if (!(v_interval <= v_upd)) {
          impl.set(v, v_interval & v_upd);
          for (term_id_t p : ttbl.parents(t)) {
            if (up_terms.find(p) == up_terms.end()) {
              up_terms.insert(p);
              queue_push(ttbl, up_queue, p);
            }
          }
        }
      }
    }

    abs.changed_terms.clear();

    if (impl.is_bottom())
      abs.set_to_bottom();
  }
};
#endif

template <typename Info> struct abstract_domain_traits<term_domain<Info>> {
  using number_t = typename Info::Number;
  using varname_t = typename Info::VariableName;
}; // end term_domain

template <typename Info> class reduced_domain_traits<term_domain<Info>> {
public:
  using term_domain_t = term_domain<Info>;
  using variable_t = typename term_domain_t::variable_t;
  using linear_constraint_system_t =
      typename term_domain_t::linear_constraint_system_t;

  static void extract(term_domain_t &dom, const variable_t &x,
                      linear_constraint_system_t &csts, bool only_equalities) {
    dom.extract(x, csts, only_equalities);
  }
};

} // namespace domains
} // namespace crab

#pragma GCC diagnostic pop
