/*******************************************************************************
 * Abstract domain described in Section 3 in the paper "An Abstract
 * Domain of Uninterpreted Functions" by Gange, Navas, Schachte,
 * Sondergaard, and Stuckey published in VMCAI'16.
 *
 * Each program variable is mapped to a Herbrand term
 * (https://en.wikipedia.org/wiki/Herbrand_structure). The join is
 * anti-unification and the meet is a pseudo-meet based on the
 * classical congruence closure algorithm. The domain is suitable to
 * infer equalities between variables.
 * 
 * TODO: conceptually uf_domain is a subset of term_domain because
 *       uf_domain implements just the Herbrand domain. We need either
 *       modify term_domain so that we can enable/disable the
 *       reduction with the numerical domain or factorize common code
 *       between uf_domain and term_domain to avoid duplication.
 ******************************************************************************/

#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/term/term_expr.hpp>
#include <crab/domains/term/term_operators.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <algorithm>
#include <map>
#include <set>
#include <utility>
#include <vector>

#include <boost/container/flat_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/iterator/transform_iterator.hpp>
#include <boost/optional.hpp>
#include <boost/utility/string_ref.hpp>

namespace crab {
namespace domains {

template <typename Number, typename VariableName>
class uf_domain final
    : public abstract_domain_api<uf_domain<Number, VariableName>> {

  using uf_domain_t = uf_domain<Number, VariableName>;
  using abstract_domain_t = abstract_domain_api<uf_domain_t>;

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
public:
  using term_id_t = typename ttbl_t::term_id_t;
  using term_t = typename ttbl_t::term_t;
  using var_set_t = std::set<variable_t>;  
private:  
  using const_term_t = typename ttbl_t::const_term_t;
  using var_term_t = typename ttbl_t::var_term_t;
  using ftor_term_t = typename ttbl_t::ftor_term_t;

  using var_map_t = boost::container::flat_map<variable_t, term_id_t>;

  struct get_term_transform :
    public std::unary_function<typename var_map_t::value_type,
			       std::pair<variable_t, const term_t*>> {
    const ttbl_t &m_ttbl;
    get_term_transform(const ttbl_t& ttbl): m_ttbl(ttbl) {}
    std::pair<variable_t, const term_t*>
    operator()(const typename var_map_t::value_type &p) const {
      return {p.first, m_ttbl.get_term_ptr(p.second)};
    }
  };
public:  
  using const_iterator =
    boost::transform_iterator<get_term_transform, typename var_map_t::const_iterator>;
private:  
  using rev_var_map_t = boost::container::flat_map<term_id_t, var_set_t>;
  using linterm_t = typename linear_expression_t::component_t;

  bool m_is_bottom;
  ttbl_t m_ttbl;
  var_map_t m_var_map;
  rev_var_map_t m_rev_var_map;

  uf_domain(bool is_top) : m_is_bottom(!is_top) {}

  uf_domain(ttbl_t &&tbl, var_map_t &&vm, rev_var_map_t &&rvm)
      : m_is_bottom(false), m_ttbl(std::move(tbl)), m_var_map(std::move(vm)),
        m_rev_var_map(std::move(rvm)) {
    check_terms(__LINE__);
  }

  void check_terms(int line) const {
    CRAB_LOG(
        "uf-check-terms",
        for (auto const &p
             : m_var_map) {
          if (!(p.second < m_ttbl.size())) {
            CRAB_ERROR("term_equiv.hpp at line=", line, ": ",
                       "term id is not the table term");
          }
        }

        for (auto kv
             : m_rev_var_map) {
          for (auto v : kv.second) {
            auto it = m_var_map.find(v);
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
    std::vector<term_id_t> forgotten /*unused*/;
    m_ttbl.deref(t, forgotten);
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
    auto it = m_rev_var_map.find(t);
    if (it != m_rev_var_map.end()) {
      it->second.erase(v);
      if (it->second.empty()) {
        m_rev_var_map.erase(it);
      }
    }
  }
  /* End manipulate the reverse variable map */

   
  void rebind_var(const variable_t &x, term_id_t tx) {
    m_ttbl.add_ref(tx);

    auto it(m_var_map.find(x));
    if (it != m_var_map.end()) {
      remove_rev_var_map((*it).second, x);
      deref((*it).second);
      m_var_map.erase(it);
    }
    m_var_map.insert(std::make_pair(x, tx));
    add_rev_var_map(m_rev_var_map, tx, x);
  }

  term_id_t term_of_const(const number_t &n) {
    boost::optional<term_id_t> opt_tn(m_ttbl.find_const(n));
    if (opt_tn) {
      return *opt_tn;
    } else {
      return m_ttbl.make_const(n);
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
    return term_of_var(v, m_var_map, m_rev_var_map, m_ttbl);
  }

  term_id_t term_of_linterm(linterm_t term) {
    if (term.first == 1) {
      return term_of_var(term.second);
    } else {
      return build_term(term::conv2termop(OP_MULTIPLICATION),
                        term_of_const(term.first), term_of_var(term.second));
    }
  }

  term_id_t build_term(term::term_operator_t op, term_id_t tx) {
    std::vector<term_id_t> ids = {tx};
    return build_term(op, ids);
  }
  
  term_id_t build_term(term::term_operator_t op, term_id_t tx, term_id_t ty) {
    std::vector<term_id_t> ids = {tx,ty};
    return build_term(op, ids);
  }
      
  term_id_t build_term(term::term_operator_t op, const std::vector<term_id_t> &ids) {
    boost::optional<term_id_t> eopt(m_ttbl.find_ftor(op, ids));
    if (eopt) {
      return *eopt;
    } else {
      term_id_t tx = m_ttbl.apply_ftor(op, ids);
      return tx;
    }
  }

  term_id_t build_linexpr(const linear_expression_t &e) {
    number_t cst = e.constant();
    typename linear_expression_t::const_iterator it(e.begin());
    if (it == e.end()) {
      return term_of_const(cst);
    }

    term_id_t t;
    if (cst == 0) {
      t = term_of_linterm(*it);
      ++it;
    } else {
      t = term_of_const(cst);
    }
    for (; it != e.end(); ++it) {
      t = build_term(term::conv2termop(OP_ADDITION), t, term_of_linterm(*it));
    }

    return t;
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

  // helper for pseudo-meet: choose one non-var term from the
  // equivalence class associated with t.
  template <typename Range>
  boost::optional<term_id_t> choose_non_var(ttbl_t &ttbl,
                                            const Range &terms) const {
    std::vector<term_id_t> non_var_terms(terms.size());
    auto it = std::copy_if(terms.begin(), terms.end(), non_var_terms.begin(),
                           [&ttbl](term_id_t t) {
                             term_t *t_ptr = ttbl.get_term_ptr(t);
                             return (t_ptr && (t_ptr->kind() == term::TERM_APP));
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


  // helper for pseudo-meet
  term_id_t build_dag_term(ttbl_t &ttbl, int t,
                           term::congruence_closure_solver<ttbl_t> &solver,
                           ttbl_t &out_ttbl, std::vector<int> &stack,
                           std::map<int, term_id_t> &cache) const {

    // already processed
    auto it = cache.find(t);
    if (it != cache.end()) {
      CRAB_LOG("uf-meet", crab::outs() << "build_dag_term. Found in cache: ";
               crab::outs() << "t" << t << " --> "
                            << "t" << it->second << "\n";);
      return it->second;
    }

    // break the cycle with a fresh variable
    if (std::find(stack.begin(), stack.end(), t) != stack.end()) {
      term_id_t v = out_ttbl.fresh_var();
      CRAB_LOG("uf-meet", crab::outs() << "build_dag_term. Detected cycle: ";
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
      CRAB_LOG("uf-meet", crab::outs()
                              << "build_dag_term. No concrete definition: ";
               crab::outs() << "t" << t << " --> "
                            << "t" << (res.first)->second << "\n";);
      return (res.first)->second;
    } else {
      term_t *f_ptr = ttbl.get_term_ptr(*f);
      // traverse recursively the term
      CRAB_LOG("uf-meet",
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
      auto res = cache.insert(std::make_pair(t, out_ttbl.apply_ftor(term::term_ftor(f_ptr), res_args)));
      stack.pop_back();
      CRAB_LOG("uf-meet",
	       crab::outs() << "build_dag_term. Finished recursive case: "
		              << "t" << t << " --> "
		              << "t" << (res.first)->second << "\n";);
      return (res.first)->second;
    }
  }

  void print_term(const term_t &t, crab_os &o) const {
    if (t.kind() == term::TERM_CONST) {
      const const_term_t *ct = static_cast<const const_term_t *>(&t);
      o << ct->val;
    } else if (t.kind() == term::TERM_VAR) {
      const var_term_t *vt = static_cast<const var_term_t *>(&t);
      o << "$VAR_" << vt->var;
    } else {
      assert(t.kind() == term::TERM_APP && "term should be a function");
      const ftor_term_t *ft = static_cast<const ftor_term_t *>(&t);
      o << ft->ftor;
      if (!ft->args.empty()) {
	o << "(";
	for (unsigned i = 0, sz = ft->args.size(); i < sz;) {
	  print_term(get_term(ft->args[i]), o);
	  ++i;
	  if (i < sz) {
	    o << ",";
	  }
	}
	o << ")";
      }
    }
  }

public:
  uf_domain_t make_top() const override { return uf_domain_t(true); }

  uf_domain_t make_bottom() const override { return uf_domain_t(false); }

  void set_to_top() override {
    uf_domain abs(true);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    uf_domain abs(false);
    std::swap(*this, abs);
  }

  uf_domain() : m_is_bottom(false) {}

  uf_domain(const uf_domain_t &o)
      : m_is_bottom(o.m_is_bottom), m_ttbl(o.m_ttbl), m_var_map(o.m_var_map),
        m_rev_var_map(o.m_rev_var_map) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    check_terms(__LINE__);
  }

  uf_domain_t &operator=(const uf_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");

    o.check_terms(__LINE__);
    if (this != &o) {
      m_is_bottom = o.m_is_bottom;
      m_ttbl = o.m_ttbl;
      m_var_map = o.m_var_map;
      m_rev_var_map = o.m_rev_var_map;
    }
    check_terms(__LINE__);
    return *this;
  }
  
  bool is_bottom() const override { return m_is_bottom; }

  bool is_top() const override { return !m_var_map.size() && !is_bottom(); }

  // Lattice operations
  bool operator<=(const uf_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    if (is_bottom()) {
      return true;
    } else if (o.is_bottom()) {
      return false;
    } else {
      // FIXME: avoid this copy
      uf_domain_t left(*this);
      uf_domain_t right(o);
      typename ttbl_t::term_map_t gen_map /*unused*/;

      // Build up the mapping of right onto left, variable by variable.
      // Assumption: the set of variables in left & right are common.
      for (auto p : left.m_var_map) {
        if (!left.m_ttbl.map_leq(right.m_ttbl, left.term_of_var(p.first),
                                 right.term_of_var(p.first), gen_map))
          return false;
      }
      return true;
    }
  }

  void operator|=(const uf_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || o.is_top()) {
      *this = o;
    } else if (o.is_bottom() || is_top()) {
      return;
    } else {
      // FIXME: avoid this copy
      uf_domain_t right(o);
      ttbl_t out_tbl;
      var_map_t out_vmap;
      rev_var_map_t out_rvmap;
      typename ttbl_t::gener_map_t gener_map /*unused*/;

      for (auto p : m_var_map) {
        const variable_t &v = p.first;
        term_id_t tx = p.second;
        term_id_t ty(right.term_of_var(v));
        term_id_t tz =
            m_ttbl.generalize(right.m_ttbl, tx, ty, out_tbl, gener_map);
        assert(tz < out_tbl.size());
        out_vmap[v] = tz;
        add_rev_var_map(out_rvmap, tz, v);
      }

      for (auto p : out_vmap) {
        out_tbl.add_ref(p.second);
      }

      m_is_bottom = false;
      std::swap(m_ttbl, out_tbl);
      std::swap(m_var_map, out_vmap);
      std::swap(m_rev_var_map, out_rvmap);
    }
  }

  uf_domain_t operator|(const uf_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom() || o.is_top()) {
      return o;
    } else if (o.is_bottom() || is_top()) {
      return *this;
    } else {
      // FIXME: avoid this copy
      uf_domain_t left(*this);
      uf_domain_t right(o);
      ttbl_t out_tbl;
      typename ttbl_t::gener_map_t gener_map /*unused*/;
      var_map_t out_vmap;
      rev_var_map_t out_rvmap;

      // For each program variable in state, compute a generalization.
      for (auto p : left.m_var_map) {
        const variable_t &v = p.first;
        term_id_t tx = p.second;
        term_id_t ty = right.term_of_var(v);
        term_id_t tz =
            left.m_ttbl.generalize(right.m_ttbl, tx, ty, out_tbl, gener_map);
        assert(tz < out_tbl.size());
        out_vmap[v] = tz;
        add_rev_var_map(out_rvmap, tz, v);
      }

      for (auto p : out_vmap) {
        out_tbl.add_ref(p.second);
      }

      uf_domain_t res(std::move(out_tbl), std::move(out_vmap),
                      std::move(out_rvmap));

      CRAB_LOG("uf", crab::outs() << "============ JOIN ==================";
               crab::outs() << *this << "\n----------------";
               crab::outs() << o << "\n----------------";
               crab::outs() << res << "\n================"
                            << "\n");

      return res;
    }
  }

  uf_domain_t operator||(const uf_domain_t &other) const override {
    return *this | other;
  }

  uf_domain_t widening_thresholds(
      const uf_domain_t &other,
      const thresholds<number_t> &ts) const override {
    return *this | other;
  }

  /**
   * Meet
   *
   * Since the domain is purely syntactic two different terms can
   * represent the same concrete expression. For instance,
   * '+'(TERM(y), TERM(z)) and '+'(TERM(z), TERM(y)) are equivalent
   * assuming that '+' is the arithmetic addition over integers since
   * that "+' is commutative. However, as Herbrand terms the
   * unification of '+'(TERM(y), TERM(z)) and '+'(TERM(z), TERM(y))
   * would fail because they are syntactically different. The other
   * tricky thing with using unification (most general unifier) as
   * meet is that the mgu of two finite terms might not be another
   * finite term but instead by a rational term which is beyond what
   * the domain can express. As a result, to ensure both soundness and
   * termination, the meet is quite imprecise. The algorithm for the
   * meet is defined in Fig 6, Sec 3 in the VMCAI'16 paper.
   * 
   * The meet described in the paper (and implemented here) is a bit
   * more precise but in the worst case does the following simple
   * thing:
   *
   *    meet({x -> TERM1}, {x-> TERM2}) = {x -> TERM1}
   *
   *    assuming TERM1 << TERM2 for some term ordering <<.
   *
   * This is even the case if TERM1 and TERM2 are two different
   * uninterpreted **symbols**.  That is,
   * 
   *    meet({x -> a}, {x-> b}) = {x -> a} assuming a << b.
   **/
  uf_domain_t operator&(const uf_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      ttbl_t out_ttbl(m_ttbl);
      std::map<term_id_t, term_id_t> copy_map;
      std::vector<int> stack;
      std::map<int, term_id_t> cache;
      var_map_t out_vmap;
      rev_var_map_t out_rvmap;

      // bring all terms to one ttbl
      for (auto p : o.m_var_map) {
        term_id_t tx = p.second;
        out_ttbl.copy_term(o.m_ttbl, tx, copy_map);
      }

      // build unifications between terms from left and right
      std::vector<std::pair<term_id_t, term_id_t>> eqs;
      for (auto p : m_var_map) {
        variable_t v(p.first);
        auto it = o.m_var_map.find(v);
        if (it != o.m_var_map.end()) {
          term_id_t tx = p.second;
          eqs.push_back(std::make_pair(tx, copy_map[it->second]));
        }
      }

      // compute equivalence classes
      term::congruence_closure_solver<ttbl_t> solver(out_ttbl);
      solver.run(eqs);

      // new map from variable to an acyclic term
      for (auto p : m_var_map) {
        const variable_t &v = p.first;
        term_id_t t_old = p.second;
        term_id_t t_new;
	if (o.m_var_map.find(v) == o.m_var_map.end()) {
	  // only on left
	  t_new = t_old;
	} else {
	  // both operands
	  t_new = build_dag_term(out_ttbl, solver.get_class(t_old),
				 solver, out_ttbl, stack, cache);
	}
        out_vmap[v] = t_new;
        add_rev_var_map(out_rvmap, t_new, v);
      }
      for (auto p : o.m_var_map) {
        const variable_t &v = p.first;
        if (out_vmap.find(v) == out_vmap.end()) {
	  // only on right operand
	  term_id_t t_new(copy_map[p.second]);
	  out_vmap[v] = t_new;
	  add_rev_var_map(out_rvmap, t_new, v);
	}
      }

      for (auto p : out_vmap) {
        out_ttbl.add_ref(p.second);
      }

      uf_domain_t res(std::move(out_ttbl), std::move(out_vmap),
                      std::move(out_rvmap));

      CRAB_LOG("uf", crab::outs() << "============ MEET ==================";
               crab::outs() << *this << "\n----------------";
               crab::outs() << o << "\n----------------";
               crab::outs() << res << "\n================"
                            << "\n");
      return res;
    }
  }

  void operator&=(const uf_domain_t &other) override {
    // TODO: avoid the copy of the left operand.
    *this = *this & other;
  }
  
  uf_domain_t operator&&(const uf_domain_t &o) const override {
    return *this & o;
  }

  // Remove a variable from the scope
  void operator-=(const variable_t &v) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    auto it(m_var_map.find(v));
    if (it != m_var_map.end()) {
      term_id_t t = (*it).second;
      m_var_map.erase(it);
      remove_rev_var_map(t, v);
      deref(t);
    }
    CRAB_LOG("uf", crab::outs()
                       << "After removing " << v << ": " << *this << "\n";);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (!is_bottom()) {
      term_id_t tx(build_linexpr(e));
      rebind_var(x, tx);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << "*** Assign " << x << ":=" << e << ":"
                                  << *this << "\n");
    }
  }

  DEFAULT_WEAK_ASSIGN(uf_domain_t)
  
  /* Begin uf-domain API */

  // Precondition:
  // - value must be > term_operator_t::first_nonreserved_value()
  static term::term_operator_t make_uf(uint32_t value) {
    auto op = term::term_operator_t::make_operator(value);
    return op;
  }

  // Assign uninterpreted symbol to x
  void set(const variable_t &x, term::term_operator_t symbol) {
    if (!is_bottom()) {
      std::vector<term_id_t> targs;
      term_id_t tx = build_term(symbol, targs);
      rebind_var(x, tx);
      check_terms(__LINE__);
    }    
  }

  // Assign uninterpreted function "functor(args)" to x
  void set(const variable_t &x, term::term_operator_t functor, const variable_vector_t &args) {
    if (!is_bottom()) {
      std::vector<term_id_t> targs;
      for (auto const&v: args) {
	targs.push_back(term_of_var(v));
      }
      term_id_t tx = build_term(functor, targs);
      rebind_var(x, tx);
      check_terms(__LINE__);      
    }
  }

  // API to traverse the term associated to a variable.
  // Pre-condition: t must be a term id generated m_ttbl.
  // FIXME: probably ugly to expose term_id_t but it's needed right now.
  const term_t &get_term(term_id_t t) const {
    const term_t *ptr_t = m_ttbl.get_term_ptr(t);
    assert(ptr_t);
    return *ptr_t;
  }

  const term_t *get_term(const variable_t &x) const {
    auto it = m_var_map.find(x);
    if (it != m_var_map.end()) {
      return &(get_term(it->second));
    } else {
      return nullptr;
    } 
  }

  // Iterate over the map from program variables to terms. This map is
  // called eta in the VMCAI'16 paper, start of section 3.
  const_iterator begin() const {
    return boost::make_transform_iterator(m_var_map.begin(),
					  get_term_transform(m_ttbl));
  }
  const_iterator end() const {
    return boost::make_transform_iterator(m_var_map.end(),
					  get_term_transform(m_ttbl));
  }
    
  // Return program variables that are mapped to a given term.
  //
  // Pre-condition: t must be generated by m_ttbl.
  // Post-condition: if the return value is not nullptr then the set
  // is not empty.
  const var_set_t *get_variables(const term_t *t) const {
    // it's fine to remove constness here because find_term won't
    // modify it.
    typename ttbl_t::term_ref_t t_ref(const_cast<term_t*>(t));
    if (boost::optional<term_id_t> t_id = m_ttbl.find_term(t_ref)) {
      auto it = m_rev_var_map.find(*t_id);
      if (it != m_rev_var_map.end()) {	
	return (it->second.empty() ? nullptr: &(it->second));
      }
    }        
    return nullptr;
  }
  /* End uf-domain API */

  
  // Apply operations to variables.

  // x = y op z
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");
    check_terms(__LINE__);

    if (!is_bottom()) {
      term_id_t tx(
          build_term(term::conv2termop(op), term_of_var(y), term_of_var(z)));
      rebind_var(x, tx);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << z << ":" << *this << "\n");
    }
  }

  // x = y op k
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      term_id_t tx(
          build_term(term::conv2termop(op), term_of_var(y), term_of_const(k)));
      rebind_var(x, tx);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << k << ":" << *this << "\n");
    }
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const uf_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");
    if (!is_bottom()) {
      CRAB_WARN("backward_assign not implemented by ", domain_name());
    }
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const uf_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    if (!is_bottom()) {
      CRAB_WARN("backward_apply not implemented by ", domain_name());
    }
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const uf_domain_t &inv) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");
    if (!is_bottom()) {
      CRAB_WARN("backward_apply not implemented by ", domain_name());
    }
  }

  void operator+=(const linear_constraint_t &cst) {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    CRAB_LOG("uf", crab::outs()
                       << "*** Before assume " << cst << ":" << *this << "\n");

    if (is_bottom()) {
      return;
    }

    using pair_var_t = std::pair<variable_t, variable_t>;

    if (boost::optional<pair_var_t> eq = get_eq_or_diseq(cst)) {
      term_id_t tx(term_of_var((*eq).first));
      term_id_t ty(term_of_var((*eq).second));
      if (cst.is_disequation()) {
        if (tx == ty) {
          set_to_bottom();
          CRAB_LOG("uf", crab::outs() << "*** After assume " << cst << ":"
                                      << *this << "\n");
          return;
        }
      } else {
        // not bother if they are already equal
        if (tx == ty) {
          return;
        }

        std::vector<int> stack;
        std::map<int, term_id_t> cache;
        // congruence closure to compute equivalence classes
        term::congruence_closure_solver<ttbl_t> solver(m_ttbl);
        std::vector<std::pair<term_id_t, term_id_t>> eqs = {{tx, ty}};
        solver.run(eqs);

        // new map from variable to an acyclic term
        for (auto p : m_var_map) {
          const variable_t &v = p.first;
          term_id_t t_old(term_of_var(v));
          term_id_t t_new = build_dag_term(m_ttbl, solver.get_class(t_old),
                                           solver, m_ttbl, stack, cache);
          rebind_var(v, t_new);
        }
      }
    }

    CRAB_LOG("uf", crab::outs()
                       << "*** After assume " << cst << ":" << *this << "\n");
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    for (auto cst : csts) {
      this->operator+=(cst);
    }
  }

  DEFAULT_ENTAILS(uf_domain_t)
  
  interval_t operator[](const variable_t &x) override {
    return at(x);
  }

  interval_t at(const variable_t &x) const override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");
    if (is_bottom()) {
      return interval_t::bottom();
    } else {
      return interval_t::top();
    }
  }  

  void apply(int_conv_operation_t /*op*/, const variable_t &dst,
             const variable_t &src) override {
    // since reasoning about infinite precision we simply assign and
    // ignore the widths.  Note that dst can be a boolean and src and
    // integer, or viceversa. 
    assign(dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      term_id_t tx(
          build_term(term::conv2termop(op), term_of_var(y), term_of_var(z)));
      rebind_var(x, tx);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << z << ":" << *this << "\n");
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (!is_bottom()) {
      term_id_t tx(
          build_term(term::conv2termop(op), term_of_var(y), term_of_const(k)));
      rebind_var(x, tx);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << k << ":" << *this << "\n");
    }
  }

  /* Array operations */

  virtual void array_init(const variable_t & /*a*/,
                          const linear_expression_t & /*elem_size*/,
                          const linear_expression_t & /*lb_idx*/,
                          const linear_expression_t & /*ub_idx*/,
                          const linear_expression_t & /*val*/) override {
    // do nothing
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t & /*elem_size*/,
                          const linear_expression_t &i) override {
    crab::CrabStats::count(domain_name() + ".count.array_read");
    crab::ScopedCrabStats __st__(domain_name() + ".array_read");

    if (!is_bottom()) {
      /**
       *  We treat the array load as an uninterpreted function
       *  lhs := array_load(a, i) -->  lhs := f(a,i)
       */
      term_id_t t_uf(
          build_term(term::TERM_OP_FUNCTION, term_of_var(a), build_linexpr(i)));
      rebind_var(lhs, t_uf);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << lhs << ":=" << a << "[" << i << "]  -- "
                                  << *this << "\n";);
    }
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t & /*elem_size*/,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool /*is_strong_update*/) override {
    // do nothing
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
                           const uf_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const uf_domain_t &invariant) override {
    *this -= lhs;
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const uf_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const uf_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const uf_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  DEFAULT_SELECT(uf_domain_t)

  // boolean operators
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    // TODO
    operator-=(lhs);
  }

  virtual void assign_bool_ref_cst(const variable_t &lhs,
                                   const reference_constraint_t &rhs) override {
    // TODO
    operator-=(lhs);
  }

  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    crab::CrabStats::count(domain_name() + ".count.assign_bool_var");
    crab::ScopedCrabStats __st__(domain_name() + ".assign_bool_var");

    if (!is_bottom()) {
      check_terms(__LINE__);
      if (is_not_rhs) {
        term_id_t tx(build_term(term::TERM_OP_NOT, term_of_var(rhs)));
        rebind_var(lhs, tx);
      } else {
        term_id_t tx(term_of_var(rhs));
        rebind_var(lhs, tx);
      }
      check_terms(__LINE__);

      CRAB_LOG(
          "uf", crab::outs() << "*** " << lhs << ":="; if (is_not_rhs) {
            crab::outs() << "not(" << rhs << ")";
          } else { crab::outs() << rhs; } crab::outs() << ":"
                                                       << *this << "\n");
    }
  }


  DEFAULT_WEAK_BOOL_ASSIGN(uf_domain_t)
  
  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y,
                                 const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply_binary_bool");
    crab::ScopedCrabStats __st__(domain_name() + ".apply_binary_bool");

    if (!is_bottom()) {
      check_terms(__LINE__);
      term_id_t tx(
          build_term(term::conv2termop(op), term_of_var(y), term_of_var(z)));
      rebind_var(x, tx);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << "*** " << x << ":=" << y << " " << op
                                  << " " << z << ":" << *this << "\n");
    }
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    // do nothing
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {
    operator-=(lhs);
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(const variable_t &lhs,
                                        const linear_constraint_t &rhs,
                                        const uf_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_assign_bool_ref_cst(const variable_t &lhs,
                                            const reference_constraint_t &rhs,
                                            const uf_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_assign_bool_var(const variable_t &lhs,
                                        const variable_t &rhs, bool is_not_rhs,
                                        const uf_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_apply_binary_bool(bool_operation_t op,
                                          const variable_t &x,
                                          const variable_t &y,
                                          const variable_t &z,
                                          const uf_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  // Region operations
  virtual void region_init(const variable_t &reg) override {
    // do nothing
  }

  virtual void region_copy(const variable_t &lhs_reg,
                           const variable_t &rhs_reg) override {
    // do nothing
  }

  virtual void region_cast(const variable_t &src_reg,
                           const variable_t &dst_reg) override {
    // do nothing
  }

  virtual void ref_make(const variable_t &ref, const variable_t &reg,
                        const variable_or_constant_t &size,
                        const allocation_site &as) override {
    // do nothing
  }

  virtual void ref_free(const variable_t &reg, const variable_t &ref) override {
    // do nothing
  }

  virtual void ref_load(const variable_t &ref, const variable_t &reg,
                        const variable_t &res) override {
    crab::CrabStats::count(domain_name() + ".count.ref_load");
    crab::ScopedCrabStats __st__(domain_name() + ".ref_load");

    if (!is_bottom()) {
      /**
       *  We treat the load as an uninterpreted function:
       *  res := ref_load(reg, ref) -->  res := f(reg,ref)
       */
      term_id_t t_uf(build_term(term::TERM_OP_FUNCTION, term_of_var(reg),
                                build_linexpr(ref)));
      rebind_var(res, t_uf);
      check_terms(__LINE__);
      CRAB_LOG("uf", crab::outs() << res << ":=ref_load(" << reg << "," << ref
                                  << ")  -- " << *this << "\n";);
    }
  }

  virtual void ref_store(const variable_t &ref, const variable_t &reg,
                         const variable_or_constant_t &val) override {
    // do nothing
  }

  virtual void ref_gep(const variable_t &ref1, const variable_t &reg1,
                       const variable_t &ref2, const variable_t &reg2,
                       const linear_expression_t &offset) override {
    // do nothing
  }

  virtual void ref_assume(const reference_constraint_t &cst) override {
    // do nothing
  }

  void ref_to_int(const variable_t &reg, const variable_t &ref_var,
                  const variable_t &int_var) override {
    // do nothing
  }

  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref_var) override {
    // do nothing
  }

  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    // do nothing
  }

  boolean_value is_null_ref(const variable_t &ref) override {
    // do nothing
    return boolean_value();
  }
  bool
  get_allocation_sites(const variable_t &ref,
                       std::vector<allocation_site> &alloc_sites) override {
    // do nothing
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &tags) override {
    // do nothing
    return false;
  }

  // Miscellaneous
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_top() || is_bottom()) {
      return;
    }

    CRAB_LOG(
        "uf", crab::outs() << "Renaming {"; for (auto v
                                                 : from) {
          crab::outs() << v << ";";
        } crab::outs() << "} with ";
        for (auto v
             : to) { crab::outs() << v << ";"; } crab::outs()
        << "}:\n";
        crab::outs() << *this << "\n";);

    auto error_if_found = [this](const variable_t &v) {
      auto it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        CRAB_ERROR(domain_name() + "::rename assumes that ", v,
                   " does not exist");
      }
    };

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      const variable_t &v = from[i];
      const variable_t &new_v = to[i];
      if (v == new_v) { // nothing to rename
        continue;
      }

      error_if_found(new_v);

      auto it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        term_id_t id = it->second;
        m_var_map.erase(it);
        m_var_map.insert(std::make_pair(new_v, id));
        remove_rev_var_map(id, v);
        add_rev_var_map(m_rev_var_map, id, new_v);
      }
    }
    CRAB_LOG("uf", crab::outs() << "RESULT=" << *this << "\n");
  }

  void forget(const variable_vector_t &variables) override {
    if (is_bottom() || is_top()) {
      return;
    }

    for (auto v : variables) {
      *this -= v;
    }
  }

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

    std::set<variable_t> s1, s2;
    variable_vector_t s3;
    for (auto p : m_var_map) {
      s1.insert(p.first);
    }
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
                          const uf_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void normalize() override {}
  void minimize() override {}

  // Output function
  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

    if (is_bottom()) {
      o << "_|_";
      return;
    }
    if (m_var_map.empty()) {
      o << "{}";
      return;
    }

    bool first = true;
    o << "{";
    for (auto p : m_var_map) {
      if (first) {
        first = false;
      } else {
        o << ", ";
      }
      o << p.first << " -> ";
      print_term(get_term(p.second), o);
    }
    o << "}";

    CRAB_LOG("ufo-print-ttbl",
             /// For debugging purposes
             o << " ttbl={" << m_ttbl << "}\n";);
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t out_csts;
    if (is_bottom()) {
      out_csts += linear_constraint_t::get_false();
    } else if (!is_top()) {
      // Extract equalities

      // Seen equalities to avoid adding twice the same.
      std::set<std::pair<variable_t, variable_t>> seen;
      for (auto &kv : m_var_map) {
        const variable_t &x = kv.first;
        term_id_t tx = kv.second;
        auto it = m_rev_var_map.find(tx);
        if (it == m_rev_var_map.end()) {
          // this shouldn't happen
          continue;
        }
        for (auto var : it->second) {
          if (var.index() != x.index()) {
            if (seen.count(std::make_pair(var, x)) <= 0) {
              seen.insert(std::make_pair(x, var));
              out_csts += linear_constraint_t(linear_expression_t(x) ==
                                              linear_expression_t(var));
            }
          }
        }
      }
    }
    return out_csts;
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

  std::string domain_name() const override { return "UFDomain"; }
}; // class uf_domain

template <typename Number, typename VariableName>
struct abstract_domain_traits<uf_domain<Number, VariableName>> {
  using number_t = Number;
  using varname_t = VariableName;
}; // end uf_domain

} // namespace domains
} // namespace crab
