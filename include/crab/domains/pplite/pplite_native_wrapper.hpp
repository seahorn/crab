#pragma once

#include <crab/config.h>

#ifdef HAVE_PPLITE

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

#include <algorithm>
#include <cassert>
#include <iostream>
#include <memory>
#include <string>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include <pplite/pplite.hh>

// If needed, reset NDEBUG and include (again) cassert
#if defined(_DEBUG) && defined(NDEBUG)
# undef NDEBUG
# include <cassert>
#endif

namespace crab {
namespace domains {

/*
  This enumeration needs to be kept in sync with
    enum pplite::dynamic::Abs_Poly::Kind
  The listed enum values are up-to-date with PPLite 0.11.
  Note: even values are for domains; odd values for _STATS variants.
*/
enum pplite_domain_id_t : unsigned int {
  POLY, POLY_STATS,
  B_POLY, B_POLY_STATS,
  F_POLY, F_POLY_STATS,
  U_POLY, U_POLY_STATS,
  UF_POLY, UF_POLY_STATS,
  P_SET, P_SET_STATS,
  FP_SET, FP_SET_STATS
};

} // namespace domains
} // namespace crab

#define PPLITE_DOMAIN_SCOPED_STATS(NAME) \
  CRAB_DOMAIN_SCOPED_STATS(this, NAME, 0)
#define PPLITE_DOMAIN_SCOPED_STATS_ASSIGN_CTOR(NAME) \
  CRAB_DOMAIN_SCOPED_STATS(&o, NAME, 0)

#define CHECK_INV(obj) \
  assert((obj).check_inv());


#define NOISY_DEBUG false

#if NOISY_DEBUG

#define PRE()                                  \
  crab::outs() << "pre " << __func__ << " : "; \
  CHECK_INV(*this);                            \
  dump();

#define PRE1(title)              \
  crab::outs() << title << "\n"; \
  PRE();

#define PRE2(title, code)        \
  crab::outs() << title << "\n"; \
  do { code } while (false);     \
  PRE();

#define POST()                                  \
  crab::outs() << "post " << __func__ << " : "; \
  CHECK_INV(*this);                             \
  dump();

#define POST1(title)             \
  crab::outs() << title << "\n"; \
  POST();

#define POST2(title, code) \
  crab::outs() << title << "\n"; \
  do { code } while (false);     \
  POST();

#define DUMP_RESULT(map, poly)     \
  crab::outs() << "\n  result = "; \
  dump(map, poly, crab::outs());

#else /* !NOISY_DEBUG */

#define PRE() CHECK_INV(*this)
#define PRE1(title) PRE()
#define PRE2(title, code) PRE()

#define POST() CHECK_INV(*this)
#define POST1(title) POST()
#define POST2(title, code) POST()

#define DUMP_RESULT(map, poly)

#endif /* !NOISY_DEBUG */

/**
 * If template parameter Number is ikos::q_number then the PPLite
 * domain will assume variables can take rational values.
 * Otherwise, all variables are assumed to only take integral values.
 **/

namespace crab {
namespace domains {

inline crab::crab_os&
operator<<(crab::crab_os& os, const pplite::dynamic::Dyn_Poly& p) {
  // EZ FIXME: how to extract the standard ostream from crab_os?
  if (&os == &crab::outs())
    p.print(std::cout);
  else if (&os == &crab::errs())
    p.print(std::cerr);
  else {
    // os is based on a string stream
    std::ostringstream ss;
    p.print(ss);
    os << ss.str();
  }
  return os;
}

template <typename Number, typename VariableName, pplite_domain_id_t K,
          class Params = PPLiteDefaultParams<Number>>
class pplite_domain final
  : public abstract_domain_api<pplite_domain<Number, VariableName, K, Params>> {
public:
  using number_t = Number;
  using varname_t = VariableName;
  using params_t = Params;

private:
  using pplite_domain_t = pplite_domain<number_t, varname_t, K, params_t>;
  using abstract_domain_t = abstract_domain_api<pplite_domain_t>;

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

private:
  using dim_type = pplite::dim_type;
  using dim_range = pplite::dim_range;
  using Affine_Expr = pplite::Affine_Expr;
  using Con = pplite::Con;
  using Cons = pplite::Cons;
  using Dims = pplite::Dims;
  using Dyn_Poly = pplite::dynamic::Dyn_Poly;
  using Dyn_Poly_Kind = pplite::dynamic::Abs_Poly::Kind;
  using Index_Set = pplite::Index_Set;
  using Integer = pplite::Integer;
  using Itv = pplite::Itv;
  using Linear_Expr = pplite::Linear_Expr;
  using Rational = pplite::Rational;
  using Spec_Elem = pplite::Spec_Elem;
  using Topol = pplite::Topol;
  using Var = pplite::Var;

  static constexpr pplite_domain_id_t domain_id = K;
  static constexpr auto poly_kind = static_cast<Dyn_Poly_Kind>(domain_id);

  // Trick to force library init code before building objects.
  struct init_check_t {
    init_check_t() {
      static bool is_initialized = false;
      if (not is_initialized) {
        bool requires_stats = (domain_id % 2 == 1);
        // Note: if (requires_stats == false), we do NOT reset
        // PPLite's noisy stats to false, because it is a library-wide
        // setting. This way, stats will be printed (at program exit)
        // if there *exists* a template instance requiring them,
        // no matter the order of template instantiation.
        if (requires_stats) {
          pplite::set_noisy_stats(requires_stats);
        }
        is_initialized = true;
      }
    }
  }; // struct init_check_t

  class var_map_t {
  public:
    using value_type = std::pair<const variable_t, dim_type>;
    using opt_dim_type = boost::optional<dim_type>;
    using inverse_t = std::vector<variable_t>;

  private:
    using map_t = std::unordered_map<variable_t, dim_type>;
    map_t map_to_dim;
    // optional: (re-) built on demand when needed
    std::unique_ptr<inverse_t> map_to_var;

  public:
    var_map_t() = default;
    var_map_t(var_map_t&& var_map) = default;
    var_map_t& operator=(var_map_t&& var_map) = default;
    ~var_map_t() = default;

    // define copy ctor and assignment (due to std::unique_ptr)
    var_map_t(const var_map_t& other)
      : map_to_dim(other.map_to_dim),
        map_to_var() {
      if (other.map_to_var)
        map_to_var.reset(new inverse_t(*other.map_to_var));
    }
    var_map_t& operator=(const var_map_t& other) {
      if (this != &other) {
        map_to_dim = other.map_to_dim;
        if (other.map_to_var)
          map_to_var.reset(new inverse_t(*other.map_to_var));
        else
          map_to_var.reset();
      }
      return *this;
    }

    using iterator = typename map_t::iterator;
    using const_iterator = typename map_t::const_iterator;

    iterator begin() { return map_to_dim.begin(); }
    iterator end() { return map_to_dim.end(); }

    const_iterator begin() const { return map_to_dim.begin(); }
    const_iterator end() const { return map_to_dim.end(); }
    const_iterator cbegin() const { return map_to_dim.cbegin(); }
    const_iterator cend() const { return map_to_dim.cend(); }

    opt_dim_type get_var_dim(const variable_t& var) const {
      auto it = map_to_dim.find(var);
      return (it != map_to_dim.end())
        ? opt_dim_type(it->second)
        : opt_dim_type();
    }

    bool has_variable(const variable_t& var) const {
      return map_to_dim.find(var) != map_to_dim.end();
    }

    iterator find(const variable_t& var) {
      return map_to_dim.find(var);
    }

    void insert(const variable_t& var, dim_type dim) {
      insert(value_type(var, dim));
    }

    void insert(value_type&& cp) {
      drop_inverse();
      map_to_dim.insert(std::move(cp));
    }

    void erase(const variable_t& key) {
      drop_inverse();
      map_to_dim.erase(key);
    }

    void erase(const_iterator it) {
      drop_inverse();
      map_to_dim.erase(it);
    }

    bool set(const variable_t& key, dim_type val) {
      auto it = map_to_dim.find(key);
      if (it == map_to_dim.end())
        return false;
      drop_inverse();
      it->second = val;
      return true;
    }

    size_t size() const {
      return map_to_dim.size();
    }

    const inverse_t& get_inverse() const {
      maybe_regen_inverse();
      return *map_to_var;
    }

    void swap(var_map_t& o) {
      map_to_dim.swap(o.map_to_dim);
      map_to_var.swap(o.map_to_var);
    }

    void print(crab::crab_os& os) const {
      for (const auto& p : map_to_dim)
        os << p.first << " --> " << p.second << "\n";
    }

  private:
    void drop_inverse() {
      map_to_var.reset();
    }

    void maybe_regen_inverse() const {
      if (map_to_var == nullptr) {
        auto& x = const_cast<var_map_t&>(*this);
        x.regen_inverse();
      }
    }

    void regen_inverse() {
      auto sz = map_to_dim.size();
      if (sz == 0) {
        map_to_var.reset(new inverse_t());
        return;
      }
      // needed because variable_t has no default ctor
      const variable_t& var = map_to_dim.begin()->first;
      map_to_var.reset(new inverse_t(sz, var));
      for (const auto& e : map_to_dim)
        (*map_to_var)[e.second] = e.first;
    }

  }; // class var_map_t

  using poly_t = Dyn_Poly;

  using interval_domain_t = ikos::interval_domain<number_t, varname_t>;
  using bound_t = ikos::bound<number_t>;
  using binding_t = typename var_map_t::value_type;
  using opt_dim_type = typename var_map_t::opt_dim_type;

  init_check_t m_init_check;
  var_map_t m_var_map;
  poly_t m_poly;

  static bool is_real() {
    return std::is_same<number_t, ikos::q_number>::value;
  }

  static bool is_integer() {
    return not is_real();
  }

  static size_t get_dims(const poly_t& ph) { return ph.space_dim(); }
  size_t get_dims() const { return get_dims(m_poly); }

  /**** Conversion helpers: from Crab to PPLite ****/

  static Integer to_Integer(const ikos::z_number& z) {
    // EZ CHECKME: it is unclear why z_number::get_mpz_t()
    // is non-const (and its use deprecated).
    auto& zz = const_cast<ikos::z_number&>(z);
    return Integer(zz.get_mpz_t());
  }
  // Note: correct rounding is up to the caller
  static Integer to_Integer(const ikos::q_number& q) {
    assert(q.denominator() == 1);
    return to_Integer(q.numerator());
  }

  static Rational to_Rational(const ikos::z_number& z) {
    assert(is_integer());
    return Rational(to_Integer(z));
  }

  static Rational to_Rational(const ikos::q_number& q) {
    assert(is_real());
    return Rational(to_Integer(q.numerator()), to_Integer(q.denominator()));
  }

  // If v is in the map then it maps v to a dimension, otherwise null
  static opt_dim_type
  get_var_dim(const var_map_t& m, const variable_t& v) {
    return m.get_var_dim(v);
  }
  opt_dim_type get_var_dim(const variable_t& v) const {
    return m_var_map.get_var_dim(v);
  }
  dim_type get_or_insert_var_dim(const variable_t& v) const {
    assert(m_var_map.size() == get_dims());
    if (auto opt_dim = get_var_dim(v))
      return *opt_dim;
    dim_type dim = m_var_map.size();
    // Semantically const
    auto& obj = const_cast<pplite_domain_t&>(*this);
    obj.m_var_map.insert(v, dim);
    obj.m_poly.add_space_dims(1);
    return dim;
  }

  Var to_Var(const variable_t& x) const {
    return Var(get_or_insert_var_dim(x));
  }

  using z_linexpr = ikos::linear_expression<ikos::z_number, varname_t>;
  using q_linexpr = ikos::linear_expression<ikos::q_number, varname_t>;

  Affine_Expr
  to_Affine_Expr(const z_linexpr& z_le) const {
    assert(is_integer());
    Affine_Expr res { to_Integer(z_le.constant()) };
    for (const auto& p : z_le) {
      auto z = to_Integer(p.first);
      auto v = to_Var(p.second);
      add_mul_assign(res.expr, z, v);
    }
    return res;
  }

  Affine_Expr
  to_Affine_Expr(const q_linexpr& q_le) const {
    Integer den = 1;
    return to_Affine_Expr_with_den(q_le, den);
  }

  Affine_Expr
  to_Affine_Expr_with_den(const z_linexpr& z_le, Integer& den) const {
    assert(den == 1);
    return to_Affine_Expr(z_le);
  }

  Affine_Expr
  to_Affine_Expr_with_den(const q_linexpr& q_le, Integer& den) const {
    assert(is_real());
    // compute into den the lcm of denominators
    den = to_Integer(q_le.constant().denominator());
    for (const auto& p : q_le) {
      const auto& q = p.first;
      const auto& q_den = q.denominator();
      if (q_den != 1)
        lcm_assign(den, den, to_Integer(q_den));
    }

    const bool need_adjustment = (den != 1);
    Integer adjust_factor;

    auto convert_with_lcm = [&](const ikos::q_number& q) -> Integer {
      auto res = to_Integer(q.numerator());
      if (need_adjustment) {
        auto q_den = to_Integer(q.denominator());
        exact_div_assign(adjust_factor, den, q_den);
        res *= adjust_factor;
      }
      return res;
    };

    Affine_Expr res;
    res.inhomo = convert_with_lcm(q_le.constant());
    for (const auto& p : q_le) {
      if (p.first == 0)
        continue;
      auto z = convert_with_lcm(p.first);
      auto v = to_Var(p.second);
      add_mul_assign(res.expr, z, v);
    }
    return res;
  }

  Con to_Con(const linear_constraint_t& cst) const {
    auto to_Type = [](typename linear_constraint_t::kind_t kind) {
      assert(kind != linear_constraint_t::DISEQUATION);
      switch (kind) {
      case linear_constraint_t::EQUALITY:
        return Con::EQUALITY;
      case linear_constraint_t::INEQUALITY:
        return Con::NONSTRICT_INEQUALITY;
      case linear_constraint_t::STRICT_INEQUALITY:
        return  Con::STRICT_INEQUALITY;
      case linear_constraint_t::DISEQUATION:
        CRAB_ERROR("conversion to Con not supported for disequations");
      }
    };

    // Note: no need to keep denominator (even when is_real() holds).
    auto ae = to_Affine_Expr(cst.expression());
    auto type = to_Type(cst.kind());
    if (type != Con::EQUALITY) {
      // Negation required (PPLite uses >= when Crab uses <=).
      pplite::neg_assign(ae);
    }
    return Con(std::move(ae), type);
  }


  /**** Conversion helpers: from PPLite to Crab ****/

  static void convert_Integer(const Integer& z_src, ikos::z_number& z_dst) {
    mpz_class tmp = z_src;
    mpz_swap(tmp.get_mpz_t(), z_dst.get_mpz_t());
  }
  static void convert_Integer(const Integer& z_src, ikos::q_number& q_dst) {
    ikos::z_number z_dst;
    convert_Integer(z_src, z_dst);
    q_dst = z_dst;
  }

  // Note: proper rounding is up to the caller
  static void convert_Rational(const Rational& q_src, ikos::z_number& z_dst) {
    assert(q_src.get_den().is_one());
    convert_Integer(q_src.get_num(), z_dst);
  }
  static void convert_Rational(const Rational& q_src, ikos::q_number& q_dst) {
    mpq_class tmp = q_src;
    q_dst = ikos::q_number::from_mpq_t(tmp.get_mpq_t());
  }

  static number_t to_number(const Integer& z) {
    number_t res;
    convert_Integer(z, res);
    return res;
  }
  static number_t to_number(const Rational& q) {
    number_t res;
    convert_Rational(q, res);
    return res;
  }

  // Note: integral refinement of itv (when safe) is up to the caller.
  static interval_t to_interval(const Itv& itv) {
    if (itv.is_empty())
      return interval_t::bottom();
    bound_t lb = itv.has_lb()
      ? bound_t(to_number(itv.lb))
      : bound_t::minus_infinity();
    bound_t ub = itv.has_ub()
      ? bound_t(to_number(itv.ub))
      : bound_t::plus_infinity();
    return interval_t(lb, ub);
  }

  linear_expression_t
  to_linear_expression(const Linear_Expr& expr,
                       const Integer& inhomo) const {
    linear_expression_t res = to_number(inhomo);
    if (expr.space_dim() == 0)
      return res;
    const auto& vars = m_var_map.get_inverse();
    for (auto d : dim_range(expr)) {
      if (not expr[d].is_zero())
        res = res + to_number(expr[d]) * vars[d];
    }
    return res;
  }

  linear_constraint_t
  to_linear_constraint(const Con& con) const {
    auto to_kind = [](typename Con::Type type) {
      switch (type) {
      case Con::EQUALITY:
        return linear_constraint_t::EQUALITY;
      case Con::NONSTRICT_INEQUALITY:
        return linear_constraint_t::INEQUALITY;
      case Con::STRICT_INEQUALITY:
        return  linear_constraint_t::STRICT_INEQUALITY;
      }
    };
    auto expr = to_linear_expression(con.linear_expr(), con.inhomo_term());
    auto kind = to_kind(con.type());
    if (kind != linear_constraint_t::EQUALITY) {
      // Negation required (PPLite uses >= when Crab uses <=).
      expr = -expr;
    }
    return linear_constraint_t(expr, kind);
  }

  /**** End of conversion helpers ****/

  // Checks the class invariant:
  // m_var_map and m_poly should be consistent.
  bool check_inv() const {
    auto& os = crab::errs();
    dim_type size = m_var_map.size();
    dim_type dim = m_poly.space_dim();
    if (size != dim) {
      os << "pplite_domain: var map size and poly dim mismatch\n";
      return false;
    }
    if (size == 0)
      return true;
    Index_Set dims;
    for (const auto& p : m_var_map)
      dims.set(p.second);
    if (dims.size() != size) {
      os << "pplite_domain: var map has index bigger than poly dim\n";
      return false;
    }
    if (dims.last() != size - 1) {
      os << "pplite_domain: var map is not injective\n";
      return false;
    }
    // all checks passed
    return true;
  }

  static void
  remove_dimensions(const Index_Set& to_remove,
                    var_map_t& m, poly_t& p) {
    assert(m.size() == get_dims(p));
    if (to_remove.empty())
      return;
    // Note: copy of inverse map is meant!
    // (m.set() and m.erase() invalidate inverse map)
    auto inverse = m.get_inverse();
    auto removed = 0;
    for (auto d : dim_range(p)) {
      const auto& var = inverse[d];
      if (to_remove.test(d)) {
        m.erase(var);
        ++removed;
      } else if (removed > 0) {
        m.set(var, d - removed);
      }
    }
    p.remove_space_dims(to_remove);
    assert(m.size() == get_dims(p));
  }

  static void
  remove_some_unconstrained(const Index_Set& uncon,
                            const var_map_t& m, const poly_t& p) {
    if (uncon.empty())
      return;
    // Here we remove some unconstrained space dims:
    // we modify m and p, but we are semantically const.
    auto& mm = const_cast<var_map_t&>(m);
    auto& pp = const_cast<poly_t&>(p);
    remove_dimensions(uncon, mm, pp);
  }

  void remove_unconstrained() const {
    remove_some_unconstrained(m_poly.get_unconstrained(),
                              m_var_map, m_poly);
  }

  static void
  merge_var_map(var_map_t& m_x, const var_map_t& m_y) {
    for (const auto& p : m_y) {
      if (not m_x.has_variable(p.first))
        m_x.insert(p.first, m_x.size());
    }
  }

  static void
  extend_to_var_map(const var_map_t& m, poly_t& p) {
    assert(m.size() >= get_dims(p));
    dim_type new_dims = m.size() - get_dims(p);
    p.add_space_dims(new_dims);
  }

  static void
  adapt_to_var_map(const var_map_t& new_map,
                   const var_map_t& old_map, poly_t& p) {
    const dim_type new_dim = new_map.size();
    const dim_type old_dim = old_map.size();
    assert(new_dim >= old_dim);
    assert(old_dim == p.space_dim());
    // add missing space dims
    p.add_space_dims(new_dim - old_dim);
    // permute p dimensions according to new_map
    Dims perm(new_dim);
    const auto& new_inverse = new_map.get_inverse();
    auto add_d = old_dim;
    for (auto d = 0; d != new_dim; ++d) {
      const auto& v = new_inverse[d];
      if (auto old_d = old_map.get_var_dim(v))
        perm[*old_d] = d;
      else {
        perm[add_d] = d;
        ++add_d;
      }
    }
    p.map_space_dims(perm);
  }

  // Dirty trick to possibly avoid copy of polyhedron
  // (if no new dim has to be added, then adaptation can be done
  // on the input polyhedron, because it is semantically const).
  // We return a pair containing a reference and a smart pointer;
  // the caller should use the reference; the smart pointer is
  // only meant to ensure proper clean up (caller should not touch it).
  static std::pair<poly_t&, std::unique_ptr<poly_t>>
  maybe_copy_and_adapt(const var_map_t& new_map,
                       const var_map_t& old_map,
                       const poly_t& old_poly) {
    const bool needs_copy = (new_map.size() > old_map.size());
    if (needs_copy) {
      std::unique_ptr<poly_t> ptr(new poly_t(old_poly));
      auto& pp = *ptr;
      adapt_to_var_map(new_map, old_map, pp);
      return { pp, std::move(ptr) };
    } else {
      // semantically const: we only permute dims
      std::unique_ptr<poly_t> ptr = nullptr;
      auto& pp = const_cast<poly_t&>(old_poly);
      adapt_to_var_map(new_map, old_map, pp);
      auto& mm = const_cast<var_map_t&>(old_map);
      mm = new_map;
      return { pp, std::move(ptr) };
    }
  }

  static bool
  merge_var_map_leq(var_map_t& m_x,
                    const var_map_t& m_y, const poly_t& p_y) {
    // Here we are checking if p_x <= p_y; if p_y constrains a variable
    // which is not in m_x (hence, unconstrained in p_x), then
    // containment does not hold and we eagerly return false.
    for (const auto& p : m_y) {
      if (not m_x.has_variable(p.first)) {
        if (p_y.constrains(Var(p.second)))
          return false;
        m_x.insert(p.first, m_x.size());
      }
    }
    return true;
  }

  // Helper for debugging purposes
  static void
  dump(const var_map_t& m, const poly_t& p,
       crab::crab_os& os = crab::outs()) {
    const char* spacer = "  ";
    auto sdim = p.space_dim();
    os << "\n";
    os << spacer << "m_poly space dim = " << sdim << "\n";
    const auto& names = m.get_inverse();
    os << spacer << "m_var_map = [ ";
    for (auto i : pplite::index_range(names)) {
      os << i << " -> " << names[i] << "; ";
    }
    os << "]\n";
    auto var_printer = [&names](std::ostream& s, Var v) {
      crab::crab_os c_os(&s);
      c_os << names[v.id()];
    };
    Var::output_function.set(var_printer);
    os << spacer << "m_poly = { " << p << " }\n";
  }

  void dump(crab::crab_os& os = crab::outs()) const {
    dump(m_var_map, m_poly, os);
  }

  // FIXME: let it produce pplite::Cons
  void
  inequalities_from_disequation(const variable_t& x, const number_t& n,
                                linear_constraint_system_t& out) const {
    interval_t xi = at(x);
    interval_t ni(n);
    interval_t new_xi = ikos::linear_interval_solver_impl
      ::trim_interval<interval_t>(xi, ni);
    if (new_xi.is_bottom()) {
      out += linear_constraint_t::get_false();
    } else if (!new_xi.is_top() && (new_xi <= xi)) {
      if (new_xi.lb().is_finite()) {
        // strenghten lb
        out += linear_constraint_t(x >= *(new_xi.lb().number()));
      }
      if (new_xi.ub().is_finite()) {
        // strenghten ub
        out += linear_constraint_t(x <= *(new_xi.ub().number()));
      }
    }
  }

  interval_t
  compute_residual(const linear_expression_t &e,
                   const variable_t& pivot) const {
    interval_t residual(-e.constant());
    for (auto kv : e) {
      const variable_t& v = kv.second;
      if (v.index() != pivot.index()) {
        residual = residual - (interval_t(kv.first) * at(v));
      }
    }
    return residual;
  }

  // FIXME: let it produce pplite::Cons
  void
  inequalities_from_disequation(const linear_expression_t &e,
                                linear_constraint_system_t &o) const {
    for (auto kv : e) {
      const variable_t& pivot = kv.second;
      interval_t i = compute_residual(e, pivot) / interval_t(kv.first);
      if (auto k = i.singleton()) {
        inequalities_from_disequation(pivot, *k, o);
      }
    }
  }

private:
  pplite_domain(var_map_t&& m, poly_t&& p, bool compact = true)
    : m_var_map(std::move(m)), m_poly(std::move(p)) {
    if (compact)
      remove_unconstrained();
    assert(check_inv());
  }

public:
  explicit pplite_domain(bool isBot = false)
    : m_var_map(),
      m_poly(isBot ? Spec_Elem::EMPTY : Spec_Elem::UNIVERSE,
             0, Topol::CLOSED, poly_kind) {
    assert(check_inv());
  }

  pplite_domain(pplite_domain_t&&) = default;
  pplite_domain_t& operator=(pplite_domain_t&&) = default;
  ~pplite_domain() = default;

  pplite_domain(const pplite_domain_t& o)
    : m_var_map(o.m_var_map), m_poly(o.m_poly) {
    PPLITE_DOMAIN_SCOPED_STATS(".copy");
  }

  pplite_domain_t& operator=(const pplite_domain_t& o) {
    PPLITE_DOMAIN_SCOPED_STATS_ASSIGN_CTOR(".copy");
    if (this != &o) {
      m_var_map = o.m_var_map;
      m_poly = o.m_poly;
    }
    return *this;
  }

  pplite_domain_t make_top() const override {
    pplite_domain_t res(false);
    return res;
  }

  pplite_domain_t make_bottom() const override {
    pplite_domain_t res(true);
    return res;
  }

  void set_to_top() override {
    m_poly.set_universe();
  }

  void set_to_bottom() override {
    m_poly.set_empty();
  }

  bool is_bottom() const override {
    return m_poly.is_empty();
  }

  bool is_top() const override {
    return m_poly.is_universe();
  }

  bool operator<=(const pplite_domain_t &o) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".leq");
    PRE();

    const auto& x = *this;
    const auto& y = o;

    if (x.is_bottom())
      return true;
    if (y.is_bottom())
      return false;
    if (y.is_top())
      return true;
    if (x.is_top())
      return false;

    auto m = x.m_var_map;
    if (not merge_var_map_leq(m, y.m_var_map, y.m_poly))
      return false;

    auto xx = maybe_copy_and_adapt(m, x.m_var_map, x.m_poly);
    auto yy = maybe_copy_and_adapt(m, y.m_var_map, y.m_poly);
    return yy.first.contains(xx.first);
  }

  void operator|=(const pplite_domain_t& o) override {
    PPLITE_DOMAIN_SCOPED_STATS(".join");
    PRE();

    auto& x = *this;
    const auto& y = o;
    if (x.is_top() || y.is_bottom())
      return;
    if (x.is_bottom()) {
      x = y;
      return;
    }
    if (y.is_top()) {
      x.set_to_top();
      return;
    }

    // Force compactness, as it improves efficiency.
    x.remove_unconstrained();
    y.remove_unconstrained();

    merge_var_map(x.m_var_map, y.m_var_map);
    extend_to_var_map(x.m_var_map, x.m_poly);
    auto yy = maybe_copy_and_adapt(x.m_var_map, y.m_var_map, y.m_poly);
    x.m_poly.join_assign(yy.first);
    POST();
  }

  pplite_domain_t operator|(const pplite_domain_t& o) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".join");
    PRE();

    const auto& x = *this;
    const auto& y = o;
    if (x.is_bottom() || y.is_top())
      return y;
    if (y.is_bottom() || x.is_top())
      return x;

    // Compactness improves efficiency.
    x.remove_unconstrained();
    y.remove_unconstrained();

    var_map_t m = x.m_var_map;
    merge_var_map(m, y.m_var_map);
    poly_t p = x.m_poly;
    extend_to_var_map(m, p);
    auto yy = maybe_copy_and_adapt(m, y.m_var_map, y.m_poly);
    p.join_assign(yy.first);
    auto res = pplite_domain_t(std::move(m), std::move(p));
    CHECK_INV(res);
    return res;
  }

  void operator&=(const pplite_domain_t& o) override {
    PPLITE_DOMAIN_SCOPED_STATS(".meet");
    PRE();

    auto& x = *this;
    const auto& y = o;
    if (x.is_bottom() || y.is_top())
      return;
    if (x.is_top() || y.is_bottom()) {
      x = y;
      return;
    }

    // Compactness improves efficiency.
    x.remove_unconstrained();
    y.remove_unconstrained();

    merge_var_map(x.m_var_map, y.m_var_map);
    extend_to_var_map(x.m_var_map, x.m_poly);
    auto yy = maybe_copy_and_adapt(x.m_var_map, y.m_var_map, y.m_poly);
    x.m_poly.intersection_assign(yy.first);
    POST();
  }

  pplite_domain_t operator&(const pplite_domain_t& o) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".meet");
    PRE();

    const auto& x = *this;
    const auto& y = o;
    if (x.is_bottom() || y.is_top())
      return x;
    if (y.is_bottom() || x.is_top())
      return y;

    // Compactness improves efficiency.
    x.remove_unconstrained();
    y.remove_unconstrained();

    var_map_t m = x.m_var_map;
    merge_var_map(m, y.m_var_map);
    poly_t p = x.m_poly;
    extend_to_var_map(m, p);
    auto yy = maybe_copy_and_adapt(m, y.m_var_map, y.m_poly);
    p.intersection_assign(yy.first);
    DUMP_RESULT(m, p);
    return pplite_domain_t(std::move(m), std::move(p));
  }

private:
  // Helper for widening
  pplite_domain_t widening_aux(const pplite_domain_t &o) const {
    PRE();

    const auto& prev = *this;
    const auto& curr = o;
    if (prev.is_bottom() || curr.is_top())
      return curr;
    if (curr.is_bottom() || prev.is_top())
      return prev;

    // Compactness improves efficiency.
    prev.remove_unconstrained();
    curr.remove_unconstrained();

    // In PPLite widening, the roles of prev and curr are swapped
    // (i.e., this is curr, rather than prev);
    // so we copy and extend curr and we adapt prev.
    var_map_t m = curr.m_var_map;
    merge_var_map(m, prev.m_var_map);
    poly_t p = curr.m_poly;
    extend_to_var_map(m, p);
    auto prev_adapted = maybe_copy_and_adapt(m, prev.m_var_map, prev.m_poly);

    // Retrieve widening spec and impl from domain parameters.
    using namespace pplite;
#if 0 // TO BE ACTIVATED WHEN SUPPORTING PARAMETERS
    const auto& params = crab_domain_params_man::get();
    auto w_spec = params.pplite_safe_widening()
      ? Widen_Spec::SAFE : Widen_Spec::RISKY;
    auto w_impl = static_cast<Widen_Impl>(params.pplite_widening_impl());
#else // USE CONTANTS FOR THE TIME BEING
    auto w_spec = Widen_Spec::RISKY;
    auto w_impl = Widen_Impl::BOXED_H79;
#endif
    // Risky widening specification requires join.
    if (w_spec == Widen_Spec::RISKY) {
      p.join_assign(prev_adapted.first);
    }
    // Now we can apply widening
    p.widening_assign(prev_adapted.first, w_impl, w_spec);

    POST();
    return pplite_domain_t(std::move(m), std::move(p));
  }

public:
  pplite_domain_t operator||(const pplite_domain_t &o) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".widening");
    return widening_aux(o);
  }

  pplite_domain_t
  widening_thresholds(const pplite_domain_t& o,
                      const thresholds<number_t>&) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".widening");
    // EZ: CHECK ME
    //////
    // We cannot refine the result of widening with
    // widening w/ thresholds over intervals because it might
    // cause non-termination.
    /////
    // This causes a loss of precision in a couple of tests:
    // - tests/domains/test2-rat.cc
    // - tests/domains/test3-rat.cc
    /////
    return widening_aux(o);
  }

  pplite_domain_t operator&&(const pplite_domain_t& o) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".narrowing");
    PRE();

    const auto& x = *this;
    const auto& y = o;
    if (x.is_bottom() || y.is_top())
      return x;
    if (y.is_bottom() || x.is_top())
      return y;

    // Compactness improves efficiency.
    x.remove_unconstrained();
    y.remove_unconstrained();

    var_map_t m = x.m_var_map;
    merge_var_map(m, y.m_var_map);
    poly_t p = x.m_poly;
    extend_to_var_map(m, p);
    auto yy = maybe_copy_and_adapt(m, y.m_var_map, y.m_poly);
    p.intersection_assign(yy.first);
    DUMP_RESULT(m, p);
    return pplite_domain_t(std::move(m), std::move(p));
  }

  void forget(const variable_vector_t& vars) override {
    PPLITE_DOMAIN_SCOPED_STATS(".forget");

    PRE2("forget",
        crab::outs()<<"   vars=[";
        for (const auto& v : vars)
          crab::outs()<< v<<", ";
        crab::outs()<<"]";
    );

    Index_Set to_remove;
    for (const auto v : vars) {
      if (auto d = m_var_map.get_var_dim(v))
        to_remove.set(*d);
    }
    remove_dimensions(to_remove, m_var_map, m_poly);
    POST();
  }

  void operator-=(const variable_t& var) override {
    PRE2("operator-=",
        crab::outs()<<"   var= "<<var;
    );

    if (auto dim = get_var_dim(var)) {
      Var v = Var(*dim);
      m_poly.remove_space_dim(v);
      m_var_map.erase(var);
      for (auto& p : m_var_map) {
        if (p.second > *dim)
          --p.second;
      }
    }
    POST();
  }

  // remove all variables except vars
  void project(const variable_vector_t& vars) override {
    PPLITE_DOMAIN_SCOPED_STATS(".project");

    PRE2("project",
        crab::outs()<<"   vars=[";
        for (const auto& v : vars)
          crab::outs()<< v<<", ";
        crab::outs()<<"]";
    );

    if (is_bottom() || is_top())
      return;

    Index_Set to_keep;
    for (const auto& v : vars) {
      if (auto d = get_var_dim(v))
        to_keep.set(*d);
    }
    if (to_keep.empty()) {
      set_to_top();
      return;
    }

    Index_Set to_remove;
    for (const auto& p : m_var_map) {
      if (not to_keep.test(p.second))
        to_remove.set(p.second);
    }
    remove_dimensions(to_remove, m_var_map, m_poly);
    POST();
  }

  interval_t operator[](const variable_t& v) override {
    return at(v);
  }

  interval_t at(const variable_t& v) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".to_intervals");

    PRE2("at(var)", crab::outs()<<"   var= "<<v; );

    if (is_bottom())
      return interval_t::bottom();
    if (is_top() || not m_var_map.has_variable(v))
      return interval_t::top();

    auto var = Var(*get_var_dim(v));
    Itv itv = m_poly.get_bounds(var);
    auto res = (is_integer() && itv.refine_as_integral())
      ? interval_t::bottom()
      : to_interval(itv);
    POST2("at(var)",
          crab::outs()<<"with var="<<v<<"="<<res<<"\n";
    );
    return res;
  }

  void operator+=(const linear_constraint_system_t& _csts) override {
    PPLITE_DOMAIN_SCOPED_STATS(".add_constraints");

    // analyze separately disequalities.
    auto preprocess = [this](const linear_constraint_t& c,
                             linear_constraint_system_t& csts) {
      switch (c.kind()) {
      case linear_constraint_t::INEQUALITY:
      case linear_constraint_t::EQUALITY:
        csts += c;
        break;
      case linear_constraint_t::STRICT_INEQUALITY:
        { // Try to convert a strict to non-strict.
          using namespace ikos::linear_constraint_impl;
          csts += strict_to_non_strict_inequality(c);
        }
        break;
      case linear_constraint_t::DISEQUATION:
        { // Maybe convert a disequation into a strict inequality
          // EZ FIXME: CHECKME: if is_integer() holds, these should be
          // refined as non-strict inequalities? (see above)
          using ctraits = constraint_simp_domain_traits<pplite_domain_t>;
          ctraits::lower_disequality(*this, c, csts);
          // Maybe convert disequation into conjunctive inequalities
          inequalities_from_disequation(c.expression(), csts);
        }
        break;
      }
    };

    PRE2("operator+=", _csts.write(crab::outs()); );

    if (is_bottom() || _csts.is_true())
      return;

    if (_csts.is_false()) {
      set_to_bottom();
      return;
    }

    linear_constraint_system_t csts;
    for (const auto& c : _csts)
      preprocess(c, csts);
    if (csts.is_true())
      return;
    if (csts.is_false()) {
      // csts can be false after breaking disequalities into inequalities
      set_to_bottom();
      return;
    }

    Cons cs;
    cs.reserve(csts.size());
    for (const auto& c : csts)
      cs.push_back(to_Con(c));
    m_poly.add_cons(cs);
    POST();
  }

  // EZ FIXME: improvable
  DEFAULT_ENTAILS(pplite_domain_t)

  void assign(const variable_t& x, const linear_expression_t &e) override {
    PPLITE_DOMAIN_SCOPED_STATS(".assign");

    PRE2("assign(x = expr)",
        crab::outs()<<x<<" = ";
        e.write(crab::outs());
    );

    if (is_bottom())
      return;

    auto x_var = to_Var(x);
    Integer den = Integer::one();
    auto ae = to_Affine_Expr_with_den(e, den);
    m_poly.affine_image(x_var, ae.expr, ae.inhomo, den);
    POST1("assign(x = expr)");
  }

private:
  // A few helpers for implementing the apply methods.

  void set(const variable_t& x, const interval_t& ival) {
    if (ival.is_top()) {
      *this -= x;
      return;
    }

    opt_dim_type x_opt_dim = m_var_map.get_var_dim(x);
    if (x_opt_dim)
      m_poly.unconstrain(Var(*x_opt_dim));
    else {
      dim_type x_dim = m_var_map.size();
      m_var_map.insert(x, x_dim);
      m_poly.add_space_dims(1);
      x_opt_dim = x_dim;
    }
    assert(x_opt_dim);
    auto x_var = Var(*x_opt_dim);

    if (auto k = ival.singleton()) {
      // add constraint x == k
      auto r = to_Rational(*k);
      Con c(r.get_den() * x_var == r.get_num());
      m_poly.add_con(c);
      assert(check_inv());
      return;
    }

    auto lb = ival.lb();
    if (lb.is_finite()) {
      // add constraint x >= lb
      auto r = to_Rational(*lb.number());
      Con c(r.get_den() * x_var >= r.get_num());
      m_poly.add_con(c);
    }
    auto ub = ival.ub();
    if (ub.is_finite()) {
      // add constraint x <= ub
      auto r = to_Rational(*ub.number());
      Con c(r.get_den() * x_var <= r.get_num());
      m_poly.add_con(c);
    }
    assert(check_inv());
  }

  void
  apply_impl_on_intervals(arith_operation_t op, const variable_t& x,
                          const interval_t& yi, const interval_t& zi) {
    interval_t xi = interval_t::top();
    switch (op) {
    case OP_ADDITION:
      xi = yi + zi;
      break;
    case OP_SUBTRACTION:
      xi = yi - zi;
      break;
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
      CRAB_ERROR("apply: operator not supported");
    }
    set(x, xi);
  }

  void
  apply_impl_on_intervals(bitwise_operation_t op, const variable_t& x,
                          const interval_t& yi, const interval_t& zi) {
    PRE();
    interval_t xi = interval_t::top();
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
    POST();
  }

  // Special handling for integer division by a constant.
  void
  apply_impl_integer_sdiv(const variable_t& x,
                          const variable_t& y, const number_t& z) {
    assert(is_integer());
    // Check for division by zero
    if (z == 0) {
      set_to_bottom();
      return;
    }
    // Check for "fake" divisions
    if (z == 1 || z == -1) {
      auto x_var = to_Var(x);
      auto y_var = to_Var(y);
      auto le = Linear_Expr(y_var);
      if (z == -1)
        pplite::neg_assign(le);
      m_poly.affine_image(x_var, le);
      POST();
      return;
    }

    // It is a proper integer division.
    assert(z <= -2 || z >= 2);
    /* Rationale: let delta = abs(z) - 1.

       If z is positive we have:
       let z' = z - 1; (note: 0 < z' < z; also, delta = z')
       y/z - 1 < x < y/z + 1
         ==>
       y/z - z'/z <= x <= y/z + z'/z
         ==> (multiplying by pos z keeps relops)
       y - z' <= z*x <= y + z'
         ==>
       -z' <= z*x - y <= z'
         ==>
       -delta <= z*x - y <= delta

       If z is negative we have:
       let z' = z + 1;
       (note: z < z' < 0; also, delta = -z')
       y/z - 1 < x < y/z + 1
         ==>
       y/z - z'/z <= x <= y/z + z'/z
         ==> (multiplying by neg z changes relops)
       y - z' >= z*x >= y + z'
         ==>
       y + z' <= z*x <= y - z'
         ==>
       z' <= z*x - y <= -z'
         ==>
       -delta <= z*x - y <= delta

       Hence we can add the same constraints in both cases.

       NOTE: precision further improved when y is known
       to be positive or negative, or if it is known to be a singleton.
    */
    auto z_int = to_Integer(z);
    auto ub = pplite::abs(z_int) - 1;
    auto lb = ub;
    pplite::neg_assign(lb);

    // refine lb/ub if y is known to be positive/negative
    if (auto y_opt_dim = m_var_map.get_var_dim(y)) {
      Var y_var(*y_opt_dim);
      auto y_itv = m_poly.get_bounds(y_var);
      y_itv.refine_as_integral();
      if (y_itv.is_empty()) {
        // may happen due to integral refinement
        set_to_bottom();
        return;
      }
      if (y_itv.is_singleton()) {
        // y is constant: compute exact integer division
        auto res = y_itv.ub.get_num() / z_int;
        auto x_var = to_Var(x);
        m_poly.unconstrain(x_var);
        m_poly.add_con(x_var == res);
        POST();
        return;
      }
      if (y_itv.has_lb() && y_itv.lb >= Rational::zero()) {
        // y is non-negative: refine lb to 0
        lb = 0;
      } else if (y_itv.has_ub() && y_itv.ub <= Rational::zero()) {
        // y is non-positive: refine ub to 0
        ub = 0;
      }
    }

    if (x == y) {
      // Adding a new dimension for (updated) x value
      auto y_var = to_Var(y);
      const auto x_old_id = y_var.id();
      const auto x_new_id = m_poly.space_dim();
      Var x_var(x_new_id);
      Cons cs {
        lb <= z_int * x_var - y_var,
        z_int * x_var - y_var <= ub
      };
      m_poly.add_space_dims(1);
      m_poly.add_cons(cs);
      // Now remove old dimension for x (i.e., y)
      m_poly.remove_space_dim(y_var);
      // Adjust map to compensate index shift
      for (auto& p : m_var_map) {
        if (p.second >= x_old_id) {
          if (p.second > x_old_id)
            --p.second;
          else
            p.second = x_new_id - 1;
        }
      }
    } else {
      // Here x and y are different
      auto x_var = to_Var(x);
      auto y_var = to_Var(y);
      Cons cs {
        lb <= z_int * x_var - y_var,
        z_int * x_var - y_var <= ub
      };
      m_poly.unconstrain(x_var);
      m_poly.add_cons(cs);
    }
    POST();
  }

  void
  apply_impl(arith_operation_t op, const variable_t& x,
             const variable_t& y, number_t z) {
    if (is_bottom())
      return;

    if ((op >= OP_ADDITION && op <= OP_MULTIPLICATION)
        ||
        (op == OP_SDIV && is_real())) {
      // Handle as affine image
      auto x_var = to_Var(x);
      auto y_var = to_Var(y);

      auto z_rat = to_Rational(z);
      auto z_num = z_rat.get_num();
      auto z_den = z_rat.get_den();

      auto lexpr = Linear_Expr(y_var);
      auto inhomo = Integer::zero();
      auto den = Integer::one();

      switch (op) {
      case OP_ADDITION:
        lexpr *= z_den;
        inhomo = z_num;
        den = z_den;
        break;
      case OP_SUBTRACTION:
        lexpr *= z_den;
        inhomo -= z_num;
        den = z_den;
        break;
      case OP_MULTIPLICATION:
        lexpr *= z_num;
        den = z_den;
        break;
      case OP_SDIV:
        assert(is_real() && not z_num.is_zero());
        if (z == 0) {
          set_to_bottom();
          return;
        }
        lexpr *= z_den;
        den = z_num;
        break;
      default:
        CRAB_ERROR("arithmetic operator not supported");
      }
      m_poly.affine_image(x_var, lexpr, inhomo, den);
      return;
    }

    // Special handling for integer division.
    if (op == OP_SDIV) {
      apply_impl_integer_sdiv(x, y, z);
      return;
    }

    // Not a linear expression: use intervals
    interval_t yi = at(y);
    interval_t zi(z);
    apply_impl_on_intervals(op, x, yi, zi);
  }

  void
  apply_impl(arith_operation_t op, const variable_t& x,
             const variable_t& y, const variable_t& z) {
    if (is_bottom())
      return;

    if (op == OP_ADDITION || op == OP_SUBTRACTION) {
      auto x_var = to_Var(x);
      auto y_var = to_Var(y);
      auto z_var = to_Var(z);
      auto lexpr = Linear_Expr(y_var);
      if (op == OP_ADDITION)
        lexpr += z_var;
      else
        lexpr -= z_var;
      m_poly.affine_image(x_var, lexpr);
      return;
    }

    if (op == OP_MULTIPLICATION) {
      interval_t yi = at(y);
      if (auto yval = yi.singleton()) {
        apply_impl(op, x, z, *yval);
        return;
      }
      interval_t zi = at(z);
      if (auto zval = zi.singleton()) {
        apply_impl(op, x, y, *zval);
        return;
      }
      apply_impl_on_intervals(op, x, yi, zi);
      return;
    }

    if (op == OP_SDIV) {
      interval_t zi = at(z);
      if (auto zval = zi.singleton()) {
        apply_impl(op, x, y, *zval);
        return;
      }
      interval_t yi = at(y);
      apply_impl_on_intervals(op, x, yi, zi);
      return;
    }

    // other operators
    interval_t yi = at(y);
    interval_t zi = at(z);
    apply_impl_on_intervals(op, x, yi, zi);
  }

public:
  void apply(arith_operation_t op, const variable_t& x,
             const variable_t& y, number_t z) override {
    PPLITE_DOMAIN_SCOPED_STATS(".apply");

    PRE2("apply(x = y [op] z)",
         crab::outs() << x << " = " << y << " " << op << " " << z;
    );
    apply_impl(op, x, y, z);
    POST1("apply(x = y [op] z)");
  }

  void apply(arith_operation_t op, const variable_t& x,
             const variable_t& y, const variable_t& z) override {
    PPLITE_DOMAIN_SCOPED_STATS(".apply");

    PRE2("apply(x = y [op] z)",
         crab::outs() << x << " = " << y << " " << op << " " << z;
    );
    apply_impl(op, x, y, z);
    POST1("apply(x = y [op] z)");
  }

  void apply(int_conv_operation_t op, const variable_t& dst,
             const variable_t& src) override {
    int_cast_domain_traits<pplite_domain_t>::apply(*this, op, dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t& x,
             const variable_t& y, const variable_t& z) override {
    PPLITE_DOMAIN_SCOPED_STATS(".apply");
    PRE2("apply(x = y [op] z)",
         crab::outs() << x << " = " << y << " " << op << " " << z;
    );
    // Convert to intervals and perform the operation
    interval_t yi = at(y);
    interval_t zi = at(z);
    apply_impl_on_intervals(op, x, yi, zi);
    POST1("apply(x = y [op] z)");
  }

  void apply(bitwise_operation_t op, const variable_t& x,
             const variable_t& y, number_t k) override {
    PPLITE_DOMAIN_SCOPED_STATS(".apply");
    PRE2("apply(x = y [op] z)",
         crab::outs() << x << " = " << y << " " << op << " " << k;
    );
    // Convert to intervals and perform the operation
    interval_t yi = at(y);
    interval_t zi(k);
    apply_impl_on_intervals(op, x, yi, zi);
    POST1("apply(x = y [op] k)");
  }

private:
  pplite_domain_t split_on_constraint(const linear_constraint_t& c) {
    // If it is a disequation, negate it before calling split operator;
    // in this case we will later swap *this and result.
    const bool need_swap = c.is_disequation();
    const auto& split_cond = need_swap ? c.negate() : c;
    auto con = to_Con(split_cond);

    // NOTE: splits are only generated for integral expressions.
    auto m_res = m_var_map;
    auto p_res = m_poly.integral_split(con);
    auto res = pplite_domain_t(std::move(m_res), std::move(p_res), false);
    if (need_swap) {
      using std::swap;
      swap(m_var_map, res.m_var_map);
      swap(m_poly, res.m_poly);
    }
    return res;
  }

public:
  void select(const variable_t& lhs,
              const linear_constraint_t& cond,
              const linear_expression_t& e1,
              const linear_expression_t& e2) override {
    PPLITE_DOMAIN_SCOPED_STATS(".select");
    if (is_bottom())
      return;
    auto& inv1 = *this;
    auto inv2 = inv1.split_on_constraint(cond);
    inv1.assign(lhs, e1);
    inv2.assign(lhs, e2);
    inv1 |= inv2;
  }

  DEFAULT_WEAK_ASSIGN(pplite_domain_t)

  void backward_assign(const variable_t& x, const linear_expression_t& e,
                       const pplite_domain_t& invariant) override {
    PPLITE_DOMAIN_SCOPED_STATS(".backward_assign");

    if (is_bottom())
      return;
    const auto& y = invariant;
    if (y.is_bottom()) {
      set_to_bottom();
      return;
    }

    auto x_var = to_Var(x);
    Integer den = Integer::one();
    auto ae = to_Affine_Expr_with_den(e, den);

    merge_var_map(m_var_map, y.m_var_map);
    extend_to_var_map(m_var_map, m_poly);
    auto yy = maybe_copy_and_adapt(m_var_map, y.m_var_map, y.m_poly);

    m_poly.affine_preimage(x_var, ae.expr, ae.inhomo, den);
    m_poly.intersection_assign(yy.first);

    CRAB_LOG("pplite", crab::outs() << "--- " << x << " :=_bwd " << e << " --> "
             << *this << "\n";);
  }

  void backward_apply(arith_operation_t op, const variable_t& x,
                      const variable_t& y, number_t z,
                      const pplite_domain_t &invariant) override {
    PPLITE_DOMAIN_SCOPED_STATS(".backward_apply");

    if (is_bottom())
      return;

    if (invariant.is_bottom()) {
      set_to_bottom();
      return;
    }

    if ((op >= OP_ADDITION && op <= OP_MULTIPLICATION)
        ||
        (op >= OP_SDIV && is_real() && z != 0)) {
      auto x_var = to_Var(x);
      auto y_var = to_Var(y);

      auto z_rat = to_Rational(z);
      auto z_num = z_rat.get_num();
      auto z_den = z_rat.get_den();

      auto lexpr = Linear_Expr(y_var);
      auto inhomo = Integer::zero();
      auto den = Integer::one();

      switch (op) {
      case OP_ADDITION:
        lexpr *= z_den;
        inhomo = z_num;
        den = z_den;
        break;
      case OP_SUBTRACTION:
        lexpr *= z_den;
        inhomo -= z_num;
        den = z_den;
        break;
      case OP_MULTIPLICATION:
        lexpr *= z_num;
        den = z_den;
        break;
      case OP_SDIV:
        assert(is_real() && z != 0);
        lexpr *= z_den;
        den = z_num;
        break;
      default:
        CRAB_ERROR("arithmetic operator not supported");
      }

      merge_var_map(m_var_map, invariant.m_var_map);
      extend_to_var_map(m_var_map, m_poly);
      auto ii = maybe_copy_and_adapt(m_var_map,
                                     invariant.m_var_map, invariant.m_poly);
      m_poly.affine_preimage(x_var, lexpr, inhomo, den);
      m_poly.intersection_assign(ii.first);

      CRAB_LOG("pplite", crab::outs() << "--- " << x << " :=_bwd " << y << op
               << z << " --> " << *this << "\n";);
      return;
    }

    // Handle "fake" integer sdiv.
    if (is_integer() && op == OP_SDIV && (z == 1 || z == -1)) {
      auto x_var = to_Var(x);
      auto y_var = to_Var(y);
      auto le = Linear_Expr(y_var);
      if (z == -1)
        pplite::neg_assign(le);
      m_poly.affine_preimage(x_var, le);
      POST();
      return;
    }

    // EZ FIXME: deal with integral division as in the forward case?

    CRAB_WARN("backwards x = y ", op, " z not implemented");
    this->operator-=(x);
  }

  void backward_apply(arith_operation_t op, const variable_t& x,
                      const variable_t& y, const variable_t& z,
                      const pplite_domain_t &invariant) override {
    PPLITE_DOMAIN_SCOPED_STATS(".backward_apply");

    if (is_bottom())
      return;
    if (invariant.is_bottom()) {
      set_to_bottom();
      return;
    }

    if (op == OP_ADDITION || op == OP_SUBTRACTION) {
      auto x_var = to_Var(x);
      auto y_var = to_Var(y);
      auto z_var = to_Var(z);
      auto lexpr = Linear_Expr(y_var);
      if (op == OP_ADDITION)
        lexpr += z_var;
      else
        lexpr -= z_var;

      merge_var_map(m_var_map, invariant.m_var_map);
      extend_to_var_map(m_var_map, m_poly);
      auto ii = maybe_copy_and_adapt(m_var_map,
                                     invariant.m_var_map, invariant.m_poly);
      m_poly.affine_preimage(x_var, lexpr);
      m_poly.intersection_assign(ii.first);
      CRAB_LOG("pplite", crab::outs() << "--- " << x << " :=_bwd " << y << op
               << z << " --> " << *this << "\n";);
      return;
    }

    // EZ FIXME: consider mult and sdiv.

    CRAB_WARN("backwards x = y ", op, " z not implemented");
    this->operator-=(x);
  }

  /// PPLite domains implement only standard abstract operations of a
  /// numerical domain so it is intended to be used as a leaf domain
  /// in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(pplite_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(pplite_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(pplite_domain_t)

  void callee_entry(const callsite_info<variable_t> &callsite,
		    const pplite_domain_t &caller) override {
    PPLITE_DOMAIN_SCOPED_STATS(".callee_entry");
    inter_abstract_operations<pplite_domain_t, Params::implement_inter_transformers>::
      callee_entry(callsite, caller, *this);
  }

  void caller_continuation(const callsite_info<variable_t> &callsite,
			   const pplite_domain_t &callee) override {
    PPLITE_DOMAIN_SCOPED_STATS(".caller_cont");
    inter_abstract_operations<pplite_domain_t, Params::implement_inter_transformers>::
      caller_continuation(callsite, callee, *this);
  }

  interval_domain_t to_interval_domain() const {
    PPLITE_DOMAIN_SCOPED_STATS(".to_interval_domain");

    interval_domain_t res;
    if (is_bottom()) {
      res.set_to_bottom();
      return res;
    }
    if (is_top())
      return res;

    pplite::BBox bbox = m_poly.get_bounding_box();
    for (const auto& p : m_var_map) {
      const auto& itv = bbox.get_bounds(p.second);
      if (is_integer() && itv.refine_as_integral()) {
        res.set_to_bottom();
        return res;
      }
      res.set(p.first, to_interval(itv));
    }
    return res;
  }

  linear_constraint_system_t
  to_linear_constraint_system() const override {
    PPLITE_DOMAIN_SCOPED_STATS(".to_linear_constraint_system");

    linear_constraint_system_t csts;
    if (is_bottom())
      csts += linear_constraint_t::get_false();
    else if (is_top())
      csts += linear_constraint_t::get_true();
    else {
      // normalize();
      for (const auto& c : m_poly.cons())
        csts += to_linear_constraint(c);
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    // Ad-hoc handling for powersets (which are disjunctive).
    if (m_poly.is_disjunctive()) {
      if (is_bottom())
        return disjunctive_linear_constraint_system_t(true /*is_false*/);
      else if (is_top())
        return disjunctive_linear_constraint_system_t(false /*is_false*/);

      // normalize();
      auto num_disj = m_poly.num_disjuncts();
      auto res = disjunctive_linear_constraint_system_t(false /*is_false*/);
      for (auto n = 0; n < num_disj; ++n) {
        linear_constraint_system_t csts;
        for (const auto& c : m_poly.disjunct_cons(n))
          csts += to_linear_constraint(c);
        res += csts;
      }
      return res;
    }

    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false())
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    else if (lin_csts.is_true())
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    else
      return disjunctive_linear_constraint_system_t(lin_csts);
  }

  void rename(const variable_vector_t& from,
              const variable_vector_t& to) override {
    PPLITE_DOMAIN_SCOPED_STATS(".rename");

    PRE2(__func__,
        crab::outs()<<"   from=[";
        for (const auto& v : from)
          crab::outs()<< v<<", ";
        crab::outs()<<"]\n";
        crab::outs()<<"   to=[";
        for (const auto& v : to)
          crab::outs()<< v<<", ";
        crab::outs()<<"]\n";
    );

    if (is_top() || is_bottom())
      return;

    CRAB_LOG("pplite",
            crab::outs() << "Renaming {"; for (auto v : from) crab::outs() << v << ";";
            crab::outs() << "} with "; for (auto v : to) crab::outs() << v << ";";
            crab::outs() << "}:\n";
            crab::outs() << *this << "\n";);

    if (from.size() == 0 || to.size() == 0)
      return;

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      const variable_t& v = from[i];
      const variable_t& new_v = to[i];
      if (v == new_v) // nothing to rename
        continue;

      // We do garbage collection of unconstrained variables only
      // after joins so it's possible to find new_v but we are ok as
      // long as it's unconstrained.

      auto it = m_var_map.find(new_v);
      if (it != m_var_map.end()) {
        dim_type dim = it->second;
        if (m_poly.constrains(Var(dim)))
          CRAB_ERROR(domain_name() + "::rename assumes that ",
                     new_v, " does not exist");
      }

      it = m_var_map.find(v);
      if (it != m_var_map.end()) {
        dim_type dim = it->second;
        m_var_map.erase(it);
        m_var_map.insert(new_v, dim);
      }
    }
    POST();
    CRAB_LOG("pplite", crab::outs() << "RESULT=" << *this << "\n");
  }

  void expand(const variable_t& x, const variable_t& dup) override {
    PPLITE_DOMAIN_SCOPED_STATS(".expand");

    PRE2("expand",
        crab::outs()<<x<<" into "<<dup;
    );

    if (is_bottom() || is_top())
      return;

    if (m_var_map.has_variable(dup)) {
      CRAB_ERROR("expand second parameter ", dup,
                 " cannot be already a variable in the domain element ", *this);
    }

    if (auto x_dim = get_var_dim(x)) {
      // map dup as the last space dimension
      m_var_map.insert(dup, m_poly.space_dim());
      m_poly.expand_space_dim(Var(*x_dim), 1);
    } else {
      // x is unconstrained, hence dup too: nothing to do
    }
    POST();
  }

  void normalize() override {
    PPLITE_DOMAIN_SCOPED_STATS(".normalize");
    PRE();
    m_poly.minimize();
    remove_unconstrained();
    POST();
  }

  // reduce the size of the internal representation
  void minimize() override {
    PPLITE_DOMAIN_SCOPED_STATS(".minimize");
    PRE();
    m_poly.minimize();
    remove_unconstrained();
    POST();
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const pplite_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void write(crab_os &os) const override {
    PPLITE_DOMAIN_SCOPED_STATS(".write");
    if (is_bottom())
      os << "_|_";
    else if (is_top())
      os << "{}";
    // EZ FIXME: avoid conversion and directly print from PPLite?
    else if (m_poly.is_disjunctive())
      os << to_disjunctive_linear_constraint_system();
    else
      os << to_linear_constraint_system();
  }

  std::string domain_name() const override {
    return abs_poly_kind_to_name(m_poly.poly_kind());
  }

}; // class pplite_domain

} // namespace domains
} // namespace crab

#endif // HAVE_PPLITE
