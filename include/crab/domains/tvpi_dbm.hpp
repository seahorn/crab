/*******************************************************************************
 *
 * Template DBM (tDBM) domain based on the paper "Template DBM: A New
 * Weakly Relational Domain for Efficient Memory-Access Validation" by
 * Su, Navas, and Gurfinkel published in VSTTE'25.
 *
 * Author: Yusen Su (yusen.su@uwaterloo.ca)
 ******************************************************************************/

#pragma once

#include <algorithm>
#include <boost/optional.hpp>
#include <functional>
#include <limits>
#include <string>
#include <type_traits>
#include <unordered_set>
#include <utility>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/inter_abstract_operations.hpp>
#include <crab/domains/tvpi/ghost_variable_manager.hpp>
#include <crab/numbers/bignums.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

namespace crab {
namespace domains {

#define TVPI_DBM_DOMAIN_SCOPED_STATS(NAME)                                     \
  CRAB_DOMAIN_SCOPED_STATS(this, NAME, 1)
#define TVPI_DBM_DOMAIN_COUNT_STATS(NAME) CRAB_DOMAIN_COUNT_STATS(NAME, 0)

/// Domain policy.  @c normalize = 1 runs full DbmTvpiSaturation
/// (@c normalize()) after every lattice/transfer operation (eager, always
/// closed); @c normalize = 0 relies on the incremental saturation hooks plus
/// explicit @c normalize() calls (lazy).
class TVPIDBMDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
  enum { normalize = 0 };
};

class TVPIDBMNormalizeParams {
public:
  enum { implement_inter_transformers = 0 };
  enum { normalize = 1 };
};

/**
 * @brief Template DBM (tDBM) abstract domain parameterized over any
 *        @c abstract_domain_api base domain.
 *
 * Given a coefficient template @c T = {1, a, b, ...} (positive integers,
 * 1 ∈ T), tDBM represents
 * @verbatim
 *   I_T  =  { a*x - b*y <= c  |  x, y ∈ Vars,  a, b ∈ T,  c ∈ ℤ }
 *   Zones (T = {1})  ⊂  tDBM (I_T)  ⊂  TVPI  ⊂  Polyhedra
 * @endverbatim
 *
 * Encoding: ghost variable @c ghost(v,N) ≙ @c N*v for each @c N ∈ T, so
 * @c a*x - b*y <= c is the unit-coefficient difference constraint
 * @c ghost(x,a) - ghost(y,b) <= c, representable by any DBM-like @p BaseDom.
 *
 * Saturation makes implied constraints explicit via the **Resultant** rule
 * (Fourier-Motzkin for two-variable constraints):
 * @verbatim
 *   a*x - b*y <= c   ∧   d*y - e*z <= f
 *   ──────────────────────────────────────────   g = gcd(b,d), λ₁=d/g, λ₂=b/g
 *   (λ₁*a)*x  -  (λ₂*e)*z  <=  λ₁*c + λ₂*f
 * @endverbatim
 * **DbmTvpiSaturation** alternates TvpiReduce (b ≠ d, this class) with
 * DbmClosure (b = d, performed by @p BaseDom whenever a constraint is
 * added).  @ref normalize runs the full algorithm;
 * @ref incremental_tvpi_reduce restores closure per added constraint.
 *
 * Original variables and their ghosts live in one base value @c m_absval;
 * @c m_vars records which variables are originals.
 *
 * @tparam BaseDom  @c abstract_domain_api domain that also provides
 *                  @c difference_bound(x,y): the tightest known @c c with
 *                  @c y - x <= c.
 * @tparam Params   @c TVPIDBMDefaultParams (no @c normalize()) or
 *                  @c TVPIDBMNormalizeParams (eager @c normalize()).
 */
template <typename BaseDom, typename Params = TVPIDBMDefaultParams>
class tvpi_dbm_domain
    : public abstract_domain_api<tvpi_dbm_domain<BaseDom, Params>> {
public:
  using tvpi_dbm_domain_t = tvpi_dbm_domain<BaseDom, Params>;
  using abstract_domain_api_t = abstract_domain_api<tvpi_dbm_domain_t>;
  using typename abstract_domain_api_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_api_t::interval_t;
  using typename abstract_domain_api_t::linear_constraint_system_t;
  using typename abstract_domain_api_t::linear_constraint_t;
  using typename abstract_domain_api_t::linear_expression_t;
  using typename abstract_domain_api_t::number_t;
  using typename abstract_domain_api_t::reference_constraint_t;
  using typename abstract_domain_api_t::variable_or_constant_t;
  using typename abstract_domain_api_t::variable_or_constant_vector_t;
  using typename abstract_domain_api_t::variable_t;
  using typename abstract_domain_api_t::variable_vector_t;
  using typename abstract_domain_api_t::varname_t;
  static_assert(std::is_same<typename abstract_domain_api_t::number_t,
                             ikos::z_number>::value,
                "abstract_domain_api_t::number_t must be the ikos::z_number");

private:
  /// Detects  D::difference_bound(x, y) -> boost::optional<number_t>.
  template <typename D, typename = void>
  struct has_difference_bound : std::false_type {};
  template <typename D>
  struct has_difference_bound<
      D, typename std::enable_if<std::is_convertible<
             decltype(std::declval<const D &>().difference_bound(
                 std::declval<const variable_t &>(),
                 std::declval<const variable_t &>())),
             boost::optional<number_t>>::value>::type> : std::true_type {};

public:
  static_assert(has_difference_bound<BaseDom>::value,
                "BaseDom must provide difference_bound(x, y) -> "
                "boost::optional<number_t>: the tightest known c with "
                "y - x <= c");

private:
  using base_domain_t = BaseDom;
  using variable_set_t = std::unordered_set<variable_t>;
  using bound_t = ikos::bound<number_t>;
  using ghost_man_t = tvpi_utils::ghost_variable_manager<variable_t>;

  /// Single base domain holding both original program variables and ghosts.
  base_domain_t m_absval;
  /// Tracks only the *original* (non-ghost) variables currently in scope.
  variable_set_t m_vars;
  /// Guards against re-entrant calls to incremental_tvpi_reduce from within
  /// tvpi_reduce() or another incremental pass.
  bool m_reduce_active = false;

  // ============================================================
  // Numeric utilities
  // ============================================================

  /// Euclidean gcd over nonnegative number_t.
  static number_t gcd_num(number_t a, number_t b) {
    while (b != 0) {
      number_t t = a % b;
      a = b;
      b = std::move(t);
    }
    return a;
  }

  /// Floor division for a positive divisor: z_number's operator/ truncates
  /// toward zero, but sound-and-tight bound scaling needs ⌊a/b⌋.
  static number_t fdiv(const number_t &a, const number_t &b) {
    number_t q = a / b;
    if (a % b != 0 && a < 0) {
      q = q - 1;
    }
    return q;
  }

  // ============================================================
  // Ghost variable helpers
  // ============================================================

  /**
   * @brief Return (creating if necessary) the ghost variable for @p v at
   *        coefficient @p coefficient, i.e., the auxiliary dimension that
   *        represents @c coefficient*v.
   *
   * Coefficient 1 is the identity: @c get_ghost_var(v, 1) == v.
   */
  variable_t get_ghost_var(const variable_t &v, const number_t &coefficient) {
    if (coefficient <= 0) {
      CRAB_ERROR("Coefficient must be > 0");
    } else if (coefficient == 1) {
      return v;
    }

    return ghost_man_t::get_or_insert(v, coefficient);
  }

  /**
   * @brief Non-creating @ref get_ghost_var: @c boost::none if the ghost does
   *        not already exist.
   *
   * Any coefficient of the global template set counts as existing
   * (materialized on demand).
   */
  boost::optional<variable_t>
  find_ghost_var(const variable_t &v, const number_t &coefficient) const {
    if (coefficient <= 0) {
      return boost::none;
    } else if (coefficient == 1) {
      return v;
    }
    // The global template set is unsigned; anything outside its range is
    // simply not in the template.
    if (coefficient.fits_int64() && static_cast<int64_t>(coefficient) <=
                                        std::numeric_limits<unsigned>::max()) {
      const auto &coeffs = crab_domain_params_man::get().coefficients();
      const unsigned c =
          static_cast<unsigned>(static_cast<int64_t>(coefficient));
      if (std::find(coeffs.begin(), coeffs.end(), c) != coeffs.end()) {
        return ghost_man_t::get_or_insert(v, coefficient);
      }
    }
    return boost::none;
  }

  /// Drop every ghost of @p v from the base domain.  Used when @p v is
  /// redefined and ghost(v,c) = c*v cannot be re-established.
  void forget_ghost_vars(const variable_t &v) {
    for (unsigned c : crab_domain_params_man::get().coefficients()) {
      if (c <= 1)
        continue; // never drop v itself
      if (auto gv = find_ghost_var(v, c)) {
        m_absval -= *gv;
      }
    }
  }

  // ============================================================
  // Linear constraint rewriting helpers
  // ============================================================

  /// Canonical form of a constraint: divide through by g = gcd of the
  /// |coefficients| (not the constant).  Over the integers
  ///   Σ ai·xi ≤ c   ⟺   Σ (ai/g)·xi ≤ ⌊c/g⌋
  /// so every TVPI constraint on a pair of variables uses one canonical
  /// coefficient pair, and constraints with equal coefficients always meet
  /// on the same ghost variables.  An equality whose constant is not
  /// divisible by g has no integer solution: the false constraint is
  /// returned.
  linear_constraint_t
  canonicalize_constraint(const linear_constraint_t &cst) const {
    if (!cst.is_inequality() && !cst.is_equality()) {
      return cst;
    }
    const linear_expression_t &e = cst.expression();
    number_t g(0);
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const number_t &coeff = (*it).first;
      number_t a = coeff < 0 ? -coeff : coeff;
      g = g == 0 ? a : gcd_num(g, a);
      if (g == 1) {
        return cst;
      }
    }
    if (g <= 1) {
      return cst;
    }
    const number_t c = -e.constant();
    if (cst.is_equality() && c % g != 0) {
      return linear_constraint_t::get_false();
    }
    linear_expression_t e2;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      e2 = e2 + ((*it).first / g) * (*it).second;
    }
    const number_t c2 = cst.is_equality() ? c / g : fdiv(c, g);
    return linear_constraint_t(e2 - c2, cst.kind());
  }

  /// GCD of all |coefficients| and the |constant| of @p cst (0 if empty);
  /// used by the Scaling rule before ghost-variable rewriting.
  number_t constraint_gcd(const linear_constraint_t &cst) const {
    const linear_expression_t &e = cst.expression();
    number_t d(0);
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const number_t &coeff = (*it).first;
      number_t abs_coeff = coeff < 0 ? -coeff : coeff;
      d = d == number_t(0) ? abs_coeff : gcd_num(d, abs_coeff);
      if (d == number_t(1))
        break;
    }
    const number_t k = e.constant();
    number_t abs_k = k < 0 ? -k : k;
    d = d == number_t(0) ? abs_k : gcd_num(d, abs_k);
    return d;
  }

  /**
   * @brief Rewrite a linear expression by replacing each term @c ci*xi with
   *        the ghost variable @c ghost(xi, |ci|), producing a unit-coefficient
   *        expression understood by the base domain.
   *
   * For positive coefficient @c ci: replace with @c +ghost(xi, ci).
   * For negative coefficient @c ci: replace with @c -ghost(xi, |ci|).
   * In fixed-coefficient mode returns @p e unchanged if any required ghost is
   * not in the tracked set.
   */
  linear_expression_t rewrite_linear_expression(const linear_expression_t &e) {
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      if (coeff == 0) {
        continue;
      } else if (coeff > 0) {
        if (find_ghost_var(v, coeff) == boost::none) { // give up
          return e;
        }
        res = res + get_ghost_var(v, coeff);
      } else {                                          // coeff < 0
        if (find_ghost_var(v, -coeff) == boost::none) { // give up
          return e;
        }
        res = res - get_ghost_var(v, -coeff);
      }
    }
    res = res + e.constant();
    return res;
  }

  /**
   * @brief Rewrite @c c1*x1 + ... + cn*xn <= k into ghost form.  With
   *        @c d = constraint_gcd(cst), which folds in @c |k| so @c k/d is
   *        exact:
   * @verbatim
   *   d > 1:  Σ sign(ci)*ghost(xi, |ci|/d)  <=  k/d    (Scaling rule)
   *   else :  Σ sign(ci)*ghost(xi, |ci|)    <=  k
   * @endverbatim
   * Since @c ghost(x,m) = m*x, both forms denote exactly @p cst.
   *
   * Ghost lookup goes through the @p resolve functor (the only difference
   * between the create and query entry points); returns @c none if @p resolve
   * fails on any term.  @c const — ghost creation, if any, is @p resolve's.
   */
  template <typename ResolveGhost>
  boost::optional<linear_constraint_t>
  rewrite_linear_constraint_impl(const linear_constraint_t &cst,
                                 ResolveGhost resolve) const {
    const linear_expression_t e = cst.expression();
    const number_t &k = e.constant();
    const number_t d = constraint_gcd(cst);
    // d > 1 => divide through (Scaling rule); else map each ci*xi directly.
    const bool scale = (d != number_t(1) && d != number_t(0));
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      const number_t abs_coeff = coeff < 0 ? -coeff : coeff;
      const bool neg = coeff < 0;
      const number_t ghost_coeff = scale ? abs_coeff / d : abs_coeff;
      boost::optional<variable_t> gv = resolve(v, ghost_coeff);
      if (!gv) { // required ghost unavailable: let the caller decide
        return boost::none;
      }
      res = neg ? res - (*gv) : res + (*gv);
    }
    res = res + (scale ? k / d : k);
    return linear_constraint_t(res, cst.kind());
  }

  /// Mutating rewrite entry point (used by @c operator+=): creates the ghosts
  /// it needs.  Returns @p cst unchanged when rewriting is impossible — the
  /// @c ecst.equal(cst) signal (fixed mode, ghost outside the template set).
  linear_constraint_t
  rewrite_linear_constraint(const linear_constraint_t &cst) {
    auto ecst = rewrite_linear_constraint_impl(
        cst,
        [this](const variable_t &v,
               const number_t &coeff) -> boost::optional<variable_t> {
          if (find_ghost_var(v, coeff) == boost::none) { // give up
            return boost::none;
          }
          return get_ghost_var(v, coeff);
        });
    return ecst ? *ecst : cst;
  }

  /// Query rewrite entry point (used by @ref entails): never creates ghosts;
  /// @c none if a required ghost does not exist.
  boost::optional<linear_constraint_t>
  try_rewrite_linear_constraint(const linear_constraint_t &cst) const {
    return rewrite_linear_constraint_impl(
        cst,
        [this](const variable_t &v,
               const number_t &coeff) -> boost::optional<variable_t> {
          return find_ghost_var(v, coeff);
        });
  }

  /// Rewrite @p e scaled by @p coefficient: each @c ci*xi maps to
  /// @c ghost(xi, ci*coefficient); @c none if such a ghost does not exist.
  /// Used to propagate @c x := e into ghost space: @c ghost(x,c) :=
  /// rewrite(e,c).
  boost::optional<linear_expression_t>
  try_rewrite_linear_expression(const linear_expression_t &e,
                                const number_t &coefficient) const {
    if (e.is_constant()) {
      // Return c*k so callers set ghost(x,c) = c*k (not forget it)
      return linear_expression_t(e.constant() * coefficient);
    }
    linear_expression_t res;
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coeff = (*it).first;
      if (coeff == 0) {
        continue;
      }
      number_t abs_coeff = coeff < 0 ? -coeff : coeff;
      bool neg = coeff < 0;
      number_t new_coeff = abs_coeff * coefficient;
      auto gv = find_ghost_var(v, new_coeff);
      if (gv == boost::none) {
        return boost::none;
      }
      if (neg) {
        res = res - (*gv);
      } else {
        res = res + (*gv);
      }
    }
    const number_t k = e.constant();
    if (k != number_t(0)) {
      res = res + k * coefficient;
    }
    return res;
  }

  /**
   * @brief Add the ghost-space equality implied by the assignment @c x := e.
   *
   * Rewrites @p e term by term into ghost form (each @c ci*xi becomes
   * @c ±ghost(xi,|ci|), the constant carries over) and adds @c x == rewrite(e)
   * — e.g. @c x := 2*y + z  adds  @c x == ghost(y,2) + z, which the base
   * domain stores as difference constraints.  If a required coefficient is
   * outside the template the rewrite returns @p e unchanged and nothing is
   * added.
   */
  void rewrite_assign(const variable_t &x, const linear_expression_t &e) {
    auto e1 = rewrite_linear_expression(e);
    if (!e1.equal(e)) {
      CRAB_LOG("tvpi-dbm-assign", crab::outs() << "processing rewritten " << x
                                               << " := " << e1 << "\n");
      m_absval += linear_constraint_t(linear_expression_t(x) - e1,
                                      linear_constraint_t::EQUALITY);
    } else {
      CRAB_LOG("tvpi-dbm-assign",
               crab::outs() << "cannot rewrite: " << x << " := " << e << "\n");
    }
  }

  /**
   * @brief Propagate an arithmetic operation with a scalar RHS into ghost
   *        space for a fixed @p coefficient.
   *
   * Given @c x := y op z (scalar @p z), adds the ghost-scaled version
   * @c ghost(x,c) := ghost(y,c) op z' where @c z' depends on the operation:
   *  - ADD/SUB: @c z' = z * coefficient  (scale the addend)
   *  - MUL/DIV: @c z' = z                (the multiplier stays)
   *
   * If @c ghost(x,c) exists but cannot be re-established (ghost(y,c)
   * absent, or the op does not distribute over scaling), it is dropped.
   */
  void rewrite_apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, number_t z,
                     const number_t &coefficient) {
    assert(coefficient > 1);

    auto gx = find_ghost_var(x, coefficient);
    if (gx == boost::none) {
      return; // no ghost(x,c): nothing to maintain
    }
    auto gy = find_ghost_var(y, coefficient);
    switch (op) {
    case OP_MULTIPLICATION:
    case OP_SDIV:
    case OP_UDIV: // c*(y op z) == (c*y) op z for scalar z
      if (gy != boost::none) {
        m_absval.apply(op, *gx, *gy, z);
        return;
      }
      break;
    case OP_ADDITION:
    case OP_SUBTRACTION: // c*(y +/- z) == c*y +/- c*z
      if (gy != boost::none) {
        m_absval.apply(op, *gx, *gy, z * coefficient);
        return;
      }
      break;
    default: // c*(y % z) etc. is not expressible in ghost space
      break;
    }
    // ghost(x,c) could not be re-established: drop it (else it stays stale,
    // still denoting c*x_old).
    m_absval -= *gx;
  }

  /**
   * @brief Propagate an arithmetic operation with a variable RHS into ghost
   *        space for a fixed @p coefficient.
   *
   * Adds @c ghost(x,c) := ghost(y,c) op ghost(z,c) — valid only for the
   * linear ops (ADD/SUB).  For everything else, or when an input ghost is
   * absent, @c ghost(x,c) is DROPPED instead of left stale.
   */
  void rewrite_apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, const variable_t &z,
                     const number_t &coefficient) {
    assert(coefficient > 1);
    auto gx = find_ghost_var(x, coefficient);
    if (gx == boost::none) {
      return; // no ghost(x,c): nothing to maintain
    }
    // Only the LINEAR ops distribute over scaling:
    //   c*(y ± z) = c*y ± c*z,   but  c*(y*z) != (c*y)*(c*z)
    //   and (c*y)/(c*z) = y/z != c*(y/z).
    if (op == OP_ADDITION || op == OP_SUBTRACTION) {
      auto gy = find_ghost_var(y, coefficient);
      auto gz = find_ghost_var(z, coefficient);
      if (gy != boost::none && gz != boost::none) {
        m_absval.apply(op, *gx, *gy, *gz);
        return;
      }
    }
    // ghost(x,c) could not be re-established: drop it.
    m_absval -= *gx;
  }

  // ============================================================
  // DBM edge-weight query
  // ============================================================

  /// Tightest known @c c with @c y - x <= c, from
  /// @c BaseDom::difference_bound (the only BaseDom method used beyond
  /// @c abstract_domain_api).  @c boost::none if unbounded.
  boost::optional<number_t> difference_bound(const variable_t &x,
                                             const variable_t &y) {
    return m_absval.difference_bound(x, y);
  }

  // ============================================================
  // Primitive TVPI constraint insertion
  // ============================================================

  /// Insert a derived Resultant output  na*lv - nb*rv <= nc, dispatching on
  /// its shape: contradiction check, bound, unit difference, or ghost edge
  /// (the last only when both coefficients are in the template).
  /// @return @c true iff redundant (nothing tightened).
  bool add_derived_constraint(const number_t &na, const variable_t &lv,
                              const number_t &nb, const variable_t &rv,
                              const number_t &nc) {
    if (na == 0 && nb == 0) {
      if (nc < 0) {
        set_to_bottom();
      }
      return true;
    } else if (na == 1 && nb == 1) {
      return add_tvpi_constraint(lv, rv, nc);
    } else if (na == 1 && nb == 0) {
      return add_ub_constraint(lv, nc);
    } else if (na == 0 && nb == 1) {
      return add_lb_constraint(rv, nc);
    } else {
      if (find_ghost_var(lv, na) && find_ghost_var(rv, nb)) {
        return add_tvpi_constraint(get_ghost_var(lv, na), get_ghost_var(rv, nb),
                                   nc);
      }
      return true; // result coefficient not in template — skip
    }
  }

  /// v₀ bound patterns for  a*x - b*y <= c :
  ///   b≠1:  y <= ub(y)  →  upper bound on x;
  ///   a≠1:  x >= lb(x)  →  lower bound on y.
  /// @return @c true iff a bound was tightened.
  bool apply_v0_patterns(const number_t &a, const variable_t &x,
                         const number_t &b, const variable_t &y,
                         const number_t &c) {
    bool changed = false;
    if (b != 1) {
      if (auto ub = m_absval.at(y).ub().number()) {
        auto ret = resultant(a, x, -b, y, c, number_t(1), number_t(-1),
                             boost::none, *ub);
        if (!add_ub_constraint(x, std::get<2>(ret)))
          changed = true;
      }
    }
    if (a != 1) {
      if (auto lb = m_absval.at(x).lb().number()) {
        auto ret = resultant(-b, y, a, x, c, number_t(-1), number_t(0),
                             boost::none, (-*lb));
        if (!add_lb_constraint(y, std::get<2>(ret)))
          changed = true;
      }
    }
    return changed;
  }

  /// Add @c ax - by <= c (@p ax, @p by ghosts, hence unit-difference in the
  /// base domain) unless the existing bound is at least as tight.
  /// @return @c true iff redundant (nothing added).
  bool add_tvpi_constraint(const variable_t &ax, const variable_t &by,
                           const number_t &c) {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".tvpi");
    // add ax - by <= c
    auto oldcopt = difference_bound(by, ax);
    if (!oldcopt || c < *oldcopt) { // if new bound is tighter
      m_absval += (ax - by <= c);
      return false;
    }
    return true;
  }

  /**
   * @brief Tighten the upper bound of @p x to @p ub if strictly tighter.
   * @return @c true if the bound was already at least as tight.
   */
  bool add_ub_constraint(const variable_t &x, const number_t &ub) {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".ub");
    // add x <= ub
    bound_t x_ub = m_absval.at(x).ub();
    bound_t new_ub = bound_t(ub);
    if (new_ub < x_ub) { // if new bound is tighter
      m_absval += (x <= ub);
      return false;
    }
    return true;
  }

  /**
   * @brief Tighten the lower bound of @p x to @p lb (i.e., add @c -x <= -lb)
   *        if strictly tighter.
   * @return @c true if the bound was already at least as tight.
   */
  bool add_lb_constraint(const variable_t &x, const number_t &lb) {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".lb");
    // add -x <= lb
    bound_t x_lb = m_absval.at(x).lb();
    bound_t new_lb = bound_t(-lb);
    if (new_lb > x_lb) { // if new bound is tighter
      m_absval += (-x <= lb);
      return false;
    }
    return true;
  }

  // ============================================================
  // Saturation math: Scaling and Resultant rules
  // ============================================================

  /**
   * @brief **Resultant** rule: eliminate the shared @p y.  Signed
   *        coefficients make one method cover every orientation:
   * @verbatim
   *   cx*x + cyl*y <= c   ∧   cyr*y + cz*z <= f     (sign cyl ≠ sign cyr)
   *   ─────────────────────────────────────────     g = gcd(|cyl|,|cyr|),
   *   (λ₁*cx)*x + (λ₂*cz)*z  <=  λ₁*c + λ₂*f        λ₁=|cyr|/g, λ₂=|cyl|/g
   * @endverbatim
   * then normalized by the Scaling rule: divide by g' = gcd of the result
   * coefficients, flooring the bound (the LHS is a multiple of g', so
   * LHS <= c ⟺ LHS/g' <= ⌊c/g'⌋).  E.g.
   * @c ax-by / @c dy-ez is @c cx=+a, cyl=-b, cyr=+d, cz=-e.
   *
   * Special cases: @p z == @c none → the z term drops (v₀ bound on x);
   * @p x == @p z → both terms fold into one signed bound.
   *
   * @return @c {new_a, new_b, new_c} for @c new_a*x - new_b*z <= new_c.
   *         @c new_c is orientation-invariant — bound-only callers read it
   *         alone.
   */
  std::tuple<number_t, number_t, number_t>
  resultant(const number_t &cx, const variable_t &x, const number_t &cyl,
            const variable_t &y, const number_t &c, const number_t &cyr,
            const number_t &cz, const boost::optional<variable_t> &z,
            const number_t &f) const {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".resultant");
    // Eliminate y from  cx*x + cyl*y <= c  and  cyr*y + cz*z <= f  (cyl, cyr
    // opposite signs).  Scale each side by λ so the y coefficients match.
    number_t abs_cyl = cyl < 0 ? -cyl : cyl;
    number_t abs_cyr = cyr < 0 ? -cyr : cyr;
    number_t gcd = gcd_num(abs_cyl, abs_cyr);
    number_t lambda1 = abs_cyr / gcd;
    number_t lambda2 = abs_cyl / gcd;
    // Signed coefficients of the derived constraint  A*x + B*z <= c_p.
    number_t A = lambda1 * cx;
    number_t B = lambda2 * cz;
    number_t c_p = c * lambda1 + f * lambda2;
    CRAB_LOG("tvpi-dbm-resultant", crab::outs() << A << x << " + " << B
                                                << (z ? z->name().str() : "v0")
                                                << " <= " << c_p << "\n");

    if (z == boost::none || x == *z) {
      // z term drops (v0) or folds into x: single-variable bound a_p*x <= c_p.
      number_t a_p = (z == boost::none) ? A : A + B;
      if (a_p == 0) {
        return {number_t(0), number_t(0), c_p};
      }
      number_t abs_a = a_p < 0 ? -a_p : a_p;
      bool neg = a_p < 0;
      // Scaling rule:  a_p*x <= c_p  ⇒  x <= ⌊c_p/|a_p|⌋.
      // neg => the surviving variable is on the negative side (-x <= new_c).
      return {neg ? number_t(0) : number_t(1), neg ? number_t(1) : number_t(0),
              fdiv(c_p, abs_a)};
    } else {
      // Genuine difference constraint: A and B have opposite signs.  Scaling
      // rule: divide by gcd(|A|,|B|); new_c = ⌊c_p/gcd⌋ is independent of the
      // orientation.
      number_t absA = A < 0 ? -A : A;
      number_t absB = B < 0 ? -B : B;
      number_t gcd2 = gcd_num(absA, absB);
      return {absA / gcd2, absB / gcd2, fdiv(c_p, gcd2)};
    }
  }

  tvpi_dbm_domain(base_domain_t &&absval, variable_set_t &&vars)
      : m_absval(std::move(absval)), m_vars(std::move(vars)) {}

public:
  DEFAULT_SELECT(tvpi_dbm_domain_t)
  BOOL_OPERATIONS_NOT_IMPLEMENTED(tvpi_dbm_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(tvpi_dbm_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(tvpi_dbm_domain_t)

  tvpi_dbm_domain() {}

  tvpi_dbm_domain(const tvpi_dbm_domain_t &o) = default;
  tvpi_dbm_domain(tvpi_dbm_domain_t &&o) = default;
  tvpi_dbm_domain_t &operator=(const tvpi_dbm_domain_t &o) = default;
  tvpi_dbm_domain_t &operator=(tvpi_dbm_domain_t &&o) = default;

  bool is_asc_phase() const override { return m_absval.is_asc_phase(); }

  void set_phase(bool is_ascending) override {
    m_absval.set_phase(is_ascending);
  }

  void set_to_top() override { m_absval.set_to_top(); }

  void set_to_bottom() override { m_absval.set_to_bottom(); }

  tvpi_dbm_domain_t make_bottom() const override {
    tvpi_dbm_domain_t res;
    res.set_to_bottom();
    return res;
  }

  tvpi_dbm_domain_t make_top() const override {
    tvpi_dbm_domain_t res;
    return res;
  }

  bool is_bottom() const override {
    bool res = m_absval.is_bottom();
    return res;
  }

  bool is_top() const override {
    bool res = m_absval.is_top();
    return res;
  }

  // ============================================================
  // Saturation: TvpiReduce and incremental saturation
  // ============================================================

  /**
   * @brief **TvpiReduce** — step 1 of DbmTvpiSaturation.
   *        Resultant over all quadruples @c (a,b,d,e) ∈ T⁴ with @c b ≠ d:
   * @verbatim
   *   ax - by ≤ c    ∧    dy - ez ≤ f    (b ≠ d)
   *   ────────────────────────────────────────────
   *   (λ₁a)x - (λ₂e)z ≤ λ₁c + λ₂f        g=gcd(b,d), λ₁=d/g, λ₂=b/g
   * @endverbatim
   * @c b = d is exact DBM transitivity — BaseDom's closure handles it after
   * each insertion.  Runs @c ⌈log₂(N)⌉-1 rounds (Nelson 1978) with early
   * fixpoint exit; the symmetric iteration over @c (x,y,a,b) covers both
   * orientations.
   *
   * v₀ patterns:  @c b≠1: ax-by≤c + y≤ub → ub(x);
   *               @c a≠1: ax-by≤c + x≥lb → lb(y).
   *
   * @note Only called from @ref normalize().
   */
  void tvpi_reduce() {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".reduce");
    if (is_bottom() || m_vars.size() < 2)
      return;
    m_reduce_active = true;

    CRAB_LOG("tvpi-dbm-reduce", crab::outs()
                                    << "Before tvpi_reduce: " << *this << "\n");

    auto vars = variable_vector_t(m_vars.begin(), m_vars.end());
    const int n = (int)vars.size();

    // ⌈log₂(N)⌉ - 1 rounds, minimum 1.  Early exit on fixpoint.
    auto ceil_log2 = [](int k) -> int {
      if (k <= 1)
        return 0;
      int r = 0;
      k -= 1;
      while (k > 0) {
        r++;
        k >>= 1;
      }
      return r;
    };
    const int rounds = std::max(1, ceil_log2(n) - 1);

    auto for_each_coeff = [&](const std::function<void(const number_t &)> &f) {
      f(number_t(1));
      for (unsigned c : crab_domain_params_man::get().coefficients()) {
        f(number_t(c));
      }
    };

    for (int round = 0; round < rounds; ++round) {
      if (is_bottom())
        return;
      bool changed = false;

      for (const auto &x : vars) {
        for_each_coeff([&](const number_t &a) {
          auto gax = get_ghost_var(x, a);
          for (const auto &y : vars) {
            if (x == y)
              continue;
            for_each_coeff([&](const number_t &b) {
              auto gby = get_ghost_var(y, b);
              auto copt = difference_bound(gby, gax); // ax - by <= c
              if (!copt)
                return;

              // General Resultant: combine ax-by<=c with dy-ez<=f for all
              // (d,e) ∈ T² with d ≠ b (b==d handled by BaseDom DbmClosure).
              for (const auto &z : vars) {
                if (z == y)
                  continue; // d·y−e·y is a scaled bound: v₀'s case

                for_each_coeff([&](const number_t &d) {
                  if (d == b)
                    return; // b==d: BaseDom's closure handles it
                  for_each_coeff([&](const number_t &e) {
                    auto gdy = get_ghost_var(y, d);
                    auto gez = get_ghost_var(z, e);
                    auto fopt = difference_bound(gez, gdy); // dy - ez <= f
                    if (!fopt)
                      return;
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs()
                                 << "Resultant(" << a << x << "-" << b << y
                                 << "<=" << *copt << ", " << d << y << "-" << e
                                 << z << "<=" << *fopt << ")\n");
                    auto ret = resultant(a, x, -b, y, *copt, d, -e, z, *fopt);
                    bool skip = add_derived_constraint(std::get<0>(ret), x,
                                                       std::get<1>(ret), z,
                                                       std::get<2>(ret));
                    if (!skip)
                      changed = true;
                    CRAB_LOG("tvpi-dbm-reduce2",
                             crab::outs()
                                 << "=>>>" << std::get<0>(ret) << x << "-"
                                 << std::get<1>(ret) << z
                                 << "<=" << std::get<2>(ret)
                                 << (skip ? ", skip" : ", added") << "\n");
                  });
                });
              }

              // v₀ patterns (bounds as virtual v₀ edges).
              if (apply_v0_patterns(a, x, b, y, *copt))
                changed = true;
            });
          }
        });
      }

      if (!changed)
        break; // fixpoint: no new constraints in this round
    }

    CRAB_LOG("tvpi-dbm-reduce", crab::outs()
                                    << "After tvpi_reduce: " << *this << "\n");
    m_reduce_active = false;
  }

  /**
   * @brief Restore closure after adding ONE constraint  a·x − b·y ≤ c  to a
   *        closed tDBM.  A new consequence combines the new constraint with
   *        at most one existing constraint on x and one on y, giving three
   *        derivation cases:
   * @verbatim
   *   Purple:  x–z  :=  (a·x − b·y ≤ c) ∘ (d·y − e·z ≤ f)      eliminate y
   *   Orange:  w–y  :=  (p·w − q·x ≤ m) ∘ (a·x − b·y ≤ c)      eliminate x
   *   Green:   w–z  :=  Purple applied to each Orange result   (worklist W_y)
   * @endverbatim
   * An elimination with equal coefficients (d = b, resp. q = a) is
   * transitivity, performed by the base domain when the constraint is added;
   * otherwise the Resultant rule applies (@ref resultant).  Bounds combine
   * through the virtual variable v₀ (@ref apply_v0_patterns).
   *
   * This is TvpiIncrementalSaturation as a two-phase split: phase 1 = the
   * Resultant cases above; phase 2 = the base domain's closure at each
   * insertion.
   * @note Caller manages @c m_reduce_active (must be @c false on entry).
   */
  void incremental_tvpi_reduce(const number_t &a, const variable_t &x,
                               const number_t &b, const variable_t &y,
                               const number_t &c) {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".incr_reduce");
    if (is_bottom() || m_vars.size() < 2)
      return;

    CRAB_LOG("tvpi-dbm-incr", crab::outs()
                                  << "incr_tvpi_reduce: " << a << x << " - "
                                  << b << y << " <= " << c << "\n");

    auto vars = variable_vector_t(m_vars.begin(), m_vars.end());

    auto for_each_coeff = [&](const std::function<void(const number_t &)> &f) {
      f(number_t(1));
      for (unsigned c : crab_domain_params_man::get().coefficients()) {
        f(number_t(c));
      }
    };

    // W_y: Orange results  na·w − nb·y ≤ nc  (y not yet eliminated).
    std::vector<std::tuple<number_t, variable_t, number_t, number_t>> wy;

    // Eliminates y between  a·x − b·y ≤ c  and every  d·y − e·z ≤ f  known
    // to the base domain.  Implements Purple; Green reuses it on W_y.
    auto eliminate_y = [&](const number_t &a, const variable_t &x,
                           const number_t &b, const variable_t &y,
                           const number_t &c) {
      for (const auto &z : vars) {
        if (z == y)
          continue; // d·y−e·y is a scaled bound: v₀'s territory
        for_each_coeff([&](const number_t &d) {
          if (d == b)
            return; // aligned: base closure's job
          for_each_coeff([&](const number_t &e) {
            if (is_bottom())
              return;
            auto gdy = get_ghost_var(y, d);
            auto gez = get_ghost_var(z, e);
            auto fopt = difference_bound(gez, gdy); // dy - ez <= f
            if (!fopt)
              return;
            CRAB_LOG("tvpi-dbm-incr",
                     crab::outs() << "  purple: " << a << x << "-" << b << y
                                  << "<=" << c << " + " << d << y << "-" << e
                                  << z << "<=" << *fopt << "\n");
            auto ret = resultant(a, x, -b, y, c, d, -e, z, *fopt);
            add_derived_constraint(std::get<0>(ret), x, std::get<1>(ret), z,
                                   std::get<2>(ret));
          });
        });
      }
    };

    // Eliminates x between every  p·w − q·x ≤ m  known to the base domain
    // and  a·x − b·y ≤ c.  Implements Orange; two-variable results go to
    // the Green worklist W_y.
    auto eliminate_x = [&](const number_t &a, const variable_t &x,
                           const number_t &b, const variable_t &y,
                           const number_t &c) {
      for (const auto &w : vars) {
        if (w == x)
          continue; // p·x−q·x is a scaled bound: v₀'s territory
        for_each_coeff([&](const number_t &p) {
          for_each_coeff([&](const number_t &q) {
            if (is_bottom())
              return;
            auto gw = get_ghost_var(w, p);
            auto gx = get_ghost_var(x, q);
            auto mopt = difference_bound(gx, gw); // pw - qx <= m
            if (!mopt)
              return;
            if (q == a) {
              // Aligned: closure derives  p·w − b·y ≤ m+c; record for Green.
              CRAB_LOG("tvpi-dbm-incr",
                       crab::outs() << "  orange(seed): " << p << w << "-" << b
                                    << y << "<=" << (*mopt + c) << "\n");
              wy.emplace_back(p, w, b, *mopt + c);
              return;
            }
            CRAB_LOG("tvpi-dbm-incr",
                     crab::outs() << "  orange: " << p << w << "-" << q << x
                                  << "<=" << *mopt << " + " << a << x << "-"
                                  << b << y << "<=" << c << "\n");
            auto ret = resultant(p, w, -q, x, *mopt, a, -b, y, c);
            add_derived_constraint(std::get<0>(ret), w, std::get<1>(ret), y,
                                   std::get<2>(ret));
            if (std::get<0>(ret) >= 1 && std::get<1>(ret) >= 1)
              wy.emplace_back(std::get<0>(ret), w, std::get<1>(ret),
                              std::get<2>(ret));
          });
        });
      }
    };

    // Case 1 (Purple): the x–z family, plus the bounds the new constraint
    // implies.
    eliminate_y(a, x, b, y, c);
    apply_v0_patterns(a, x, b, y, c);

    // Case 2 (Orange): the w–y family; fills W_y.
    eliminate_x(a, x, b, y, c);

    // Case 3 (Green): the w–z family — eliminate the remaining y of each
    // W_y constraint.  A single pass; results are not fed back.
    for (const auto &e : wy) {
      if (is_bottom())
        break;
      eliminate_y(std::get<0>(e), std::get<1>(e), std::get<2>(e), y,
                  std::get<3>(e));
      apply_v0_patterns(std::get<0>(e), std::get<1>(e), std::get<2>(e), y,
                        std::get<3>(e));
    }
  }

  /// Re-entrancy-guarded wrapper around @ref incremental_tvpi_reduce:
  /// no-op while another reduce is running.
  void incremental_reduce(const number_t &a, const variable_t &x,
                          const number_t &b, const variable_t &y,
                          const number_t &c) {
    if (m_reduce_active) {
      return;
    }
    m_reduce_active = true;
    incremental_tvpi_reduce(a, x, b, y, c);
    m_reduce_active = false;
  }

  /// Run incremental saturation if @p cst is TVPI over original variables:
  /// @c a*x - b*y <= k  (strict integer inequality: k-1; equality: both
  /// orientations).  Bounds, sums @c a*x + b*y <= k, disequations and
  /// constraints with more than two variables are not reduced here; they
  /// are handled by the base domain and by @ref tvpi_reduce.
  void incremental_reduce_constraint(const linear_constraint_t &cst) {
    if (m_reduce_active || is_bottom()) {
      return;
    }
    if (!cst.is_inequality() && !cst.is_strict_inequality() &&
        !cst.is_equality()) {
      return;
    }
    // cst is  e ⋈ 0; require e = c1*v1 + c2*v2 + k with c1, c2 of opposite
    // sign, i.e., the difference form  a*x - b*y <= -k.
    const linear_expression_t &e = cst.expression();
    auto it = e.begin(), et = e.end();
    if (it == et) {
      return;
    }
    number_t c1 = (*it).first;
    variable_t v1 = (*it).second;
    ++it;
    if (it == et) {
      // Single-variable bound: not reduced incrementally (bounds are the
      // most frequent constraint; a per-bound alignment sweep costs
      // O(|T|²·N)).  @ref tvpi_reduce derives its consequences.
      return;
    }
    number_t c2 = (*it).first;
    variable_t v2 = (*it).second;
    ++it;
    if (it != et || c1 == 0 || c2 == 0 || (c1 < 0) == (c2 < 0)) {
      return;
    }
    if (c1 < 0) { // orient as  c1*v1 - |c2|*v2
      std::swap(c1, c2);
      std::swap(v1, v2);
    }
    c2 = -c2;
    number_t k = -e.constant();
    if (cst.is_strict_inequality()) { // integers:  e < 0  ⟺  e <= -1
      k = k - 1;
    }
    incremental_reduce(c1, v1, c2, v2, k);
    if (cst.is_equality() && !is_bottom()) {
      incremental_reduce(c2, v2, c1, v1, -k);
    }
  }

  // ============================================================
  // Lattice operations
  // ============================================================

  /// Abstract inclusion @c *this ⊑ other, decided by the base domain.
  /// Precise on closed operands (the eager variant keeps values closed).
  bool operator<=(const tvpi_dbm_domain_t &other) const override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".leq");
    if (is_bottom() || other.is_top()) {
      return true;
    } else if (is_top() || other.is_bottom()) {
      return false;
    } else {
      CRAB_LOG("tvpi-dbm-leq", crab::outs() << "[leq]\n"
                                            << *this << "\n<=\n"
                                            << other << "\n");
      return m_absval <= other.m_absval;
    }
  }

  /**
   * @brief In-place abstract join @c *this = *this ⊔ other.
   *
   * Computes the base-domain join, then calls @c normalize() when
   * @c Params::normalize == 1.  @c m_vars is set to the union of both
   * operands' tracked variables.
   */
  void operator|=(const tvpi_dbm_domain_t &other) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".self_join");
    if (is_bottom() || other.is_top()) {
      *this = other;
    } else if (other.is_bottom() || is_top()) {
      // do nothing
    } else {
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[join]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
      m_absval |= other.m_absval;
      m_vars.insert(other.m_vars.begin(), other.m_vars.end());
      if constexpr (Params::normalize)
        normalize();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[res]\n" << *this << "\n");
    }
  }

  /// Abstract join @c *this ⊔ other, returning a fresh domain.
  tvpi_dbm_domain_t operator|(const tvpi_dbm_domain_t &other) const override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".join");
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[join]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
      tvpi_dbm_domain_t res = *this;
      res.m_absval |= other.m_absval;
      res.m_vars.insert(other.m_vars.begin(), other.m_vars.end());
      if constexpr (Params::normalize)
        res.normalize();
      CRAB_LOG("tvpi-dbm-join", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  /**
   * @brief In-place abstract meet @c *this = *this ⊓ other.
   *
   * Computes the base-domain meet; sets to bottom if the result is detected
   * as unsatisfiable.  Calls @c normalize() when @c Params::normalize == 1.
   */
  void operator&=(const tvpi_dbm_domain_t &other) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".self_meet");
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (other.is_bottom() || is_top()) {
      *this = other;
    } else {
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[meet]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
      m_absval &= other.m_absval;
      m_vars.insert(other.m_vars.begin(), other.m_vars.end());
      if (m_absval.is_bottom()) {
        set_to_bottom();
      } else if constexpr (Params::normalize) {
        normalize();
      }
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[res]\n" << *this << "\n");
    }
  }

  /// Abstract meet @c *this ⊓ other, returning a fresh domain.
  tvpi_dbm_domain_t operator&(const tvpi_dbm_domain_t &other) const override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".meet");
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[meet]\n"
                                             << *this << "\nwith\n"
                                             << other << "\n");
      tvpi_dbm_domain_t res = *this;
      res.m_absval &= other.m_absval;
      res.m_vars.insert(other.m_vars.begin(), other.m_vars.end());
      if (res.m_absval.is_bottom()) {
        res.set_to_bottom();
      } else if constexpr (Params::normalize) {
        res.normalize();
      }
      CRAB_LOG("tvpi-dbm-meet", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  /// Abstract widening @c *this ▽ other.
  tvpi_dbm_domain_t operator||(const tvpi_dbm_domain_t &other) const override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".widen");
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[widen]\n"
                                              << *this << "\nwith\n"
                                              << other << "\n");
      base_domain_t out_absval = m_absval || other.m_absval;
      variable_set_t out_vars = m_vars;
      out_vars.insert(other.m_vars.begin(), other.m_vars.end());
      tvpi_dbm_domain_t res(std::move(out_absval), std::move(out_vars));
      if constexpr (Params::normalize)
        res.normalize();
      CRAB_LOG("tvpi-dbm-widen", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  /// Widening with thresholds.
  tvpi_dbm_domain_t
  widening_thresholds(const tvpi_dbm_domain_t &other,
                      const thresholds<number_t> &ts) const override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".widen.thresholds");
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      base_domain_t out_absval =
          m_absval.widening_thresholds(other.m_absval, ts);
      variable_set_t out_vars = m_vars;
      out_vars.insert(other.m_vars.begin(), other.m_vars.end());
      tvpi_dbm_domain_t res(std::move(out_absval), std::move(out_vars));
      if constexpr (Params::normalize)
        res.normalize();
      return res;
    }
  }

  /// Abstract narrowing @c *this △ other.
  tvpi_dbm_domain_t operator&&(const tvpi_dbm_domain_t &other) const override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".narrow");
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      base_domain_t out_absval = m_absval && other.m_absval;
      variable_set_t out_vars = m_vars;
      out_vars.insert(other.m_vars.begin(), other.m_vars.end());
      tvpi_dbm_domain_t res(std::move(out_absval), std::move(out_vars));
      if constexpr (Params::normalize)
        res.normalize();
      CRAB_LOG("tvpi-dbm-narrow", crab::outs() << "[res]\n" << res << "\n");
      return res;
    }
  }

  // ============================================================
  // Transfer functions
  // ============================================================

  /**
   * @brief Assume @p csts.  Per constraint:
   *  1. add it to @c m_absval (BaseDom closure runs on insertion);
   *  2. add each c-scaled ghost version, @c c ∈ T
   *     (e.g. @c i-len<=-1 → @c ghost(i,c)-ghost(len,c)<=-c);
   *  3. add the GCD-normalized ghost form (Scaling rule);
   *  4. lazy variant: incremental saturation on the new TVPI edge.
   * Eager variant runs @c normalize() once at the end instead of step 4.
   */
  void operator+=(const linear_constraint_system_t &csts) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".+=");
    CRAB_LOG("tvpi-dbm-+=",
             crab::outs() << "Before assume(" << csts << ")=" << *this << "\n");
    if (!is_bottom()) {
      for (auto const &raw_cst : csts) {
        if (raw_cst.is_contradiction()) {
          set_to_bottom();
          break;
        }

        if (raw_cst.is_tautology()) {
          continue;
        }

        // Canonical (lowest coefficient terms) form; may detect an
        // integer-infeasible equality outright.
        const linear_constraint_t cst = canonicalize_constraint(raw_cst);
        if (cst.is_contradiction()) {
          set_to_bottom();
          break;
        }

        CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                    << "processing original: " << cst << "\n");

        for (auto const &v : cst.variables()) {
          m_vars.insert(v);
        }

        m_absval += cst;
        if (m_absval.is_bottom()) {
          set_to_bottom();
          break;
        }

        // Add coefficient-scaled ghost versions: for each tracking
        // coefficient c of the global template, add the
        // c-scaled constraint using ghost variables, e.g. "i-x<=-1" with c=4
        // adds "ghost(i,4)-ghost(x,4)<=-4".
        {
          // Scale by the global template list in both modes.
          const auto &coeffs = crab_domain_params_man::get().coefficients();
          for (unsigned c : coeffs) {
            auto scaled_expr =
                try_rewrite_linear_expression(cst.expression(), c);
            if (scaled_expr) {
              linear_constraint_t scaled_cst(*scaled_expr, cst.kind());
              CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                          << "process_tvpi[*" << c
                                          << "]: " << scaled_cst << "\n");
              m_absval += scaled_cst;
              if (m_absval.is_bottom()) {
                set_to_bottom();
                break;
              }
            } else {
              CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                          << "process_tvpi[*" << c
                                          << "]: SKIP (ghost not found for: "
                                          << cst.expression() << ")\n");
            }
          }
          if (is_bottom())
            break;
        }

        auto ecst = rewrite_linear_constraint(cst);
        if (ecst.equal(cst)) {
          CRAB_LOG("tvpi-dbm-+=", crab::outs()
                                      << "cannot rewrite: " << cst << "\n");
        } else {
          CRAB_LOG("tvpi-dbm-+=",
                   crab::outs() << "processing rewritten " << ecst << "\n");
          m_absval += ecst;
          if (m_absval.is_bottom()) {
            set_to_bottom();
            break;
          }
        }

        // Lazy variant: restore closure incrementally for the newly added
        // TVPI constraint (the eager variant runs the full reduce below).
        if constexpr (!Params::normalize) {
          incremental_reduce_constraint(cst);
          if (is_bottom()) {
            break;
          }
        }

      } // end for
      if constexpr (Params::normalize)
        normalize();
    }

    CRAB_LOG("tvpi-dbm-+=",
             crab::outs() << "After assume(" << csts << ")=" << *this << "\n");
  }

  /**
   * @brief Does the domain imply @p cst?  Tries, in cost order:
   *  1. direct @c m_absval.entails;
   *  2. same on the ghost-rewritten form of @p cst;
   *  3. refutation: assume the integer negation on a copy (@c += runs
   *     incremental saturation); bottom ⇒ entailed.  Equalities split into
   *     their two inequalities;
   *  4. full @ref normalize on a copy, then (1)+(2) again — materialises
   *     ghost-chain facts the cheaper steps may miss.
   */
  bool entails(const linear_constraint_t &cst) const override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".entails");
    if (is_bottom()) {
      return true;
    } else if (cst.is_tautology()) {
      return true;
    } else if (cst.is_contradiction()) {
      return false;
    } else if (m_absval.entails(cst)) {
      return true;
    } else {
      auto ecst = try_rewrite_linear_constraint(cst);
      if (ecst && !(*ecst).equal(cst)) {
        if (m_absval.entails((*ecst))) {
          return true;
        }
      }

      // Refutation: bottom on γ(this) ∧ ¬cst implies γ(this) ⊆ cst,
      // regardless of closure.
      auto refuted = [this](const linear_constraint_t &neg_cst) {
        tvpi_dbm_domain_t tmp = *this;
        tmp += neg_cst;
        return tmp.is_bottom();
      };
      if (cst.is_inequality() || cst.is_strict_inequality()) {
        if (refuted(cst.negate())) {
          return true;
        }
      } else if (cst.is_equality()) {
        // e == 0  ⟺  e <= 0 ∧ -e <= 0: refute each side's negation.
        linear_constraint_t leq(cst.expression(),
                                linear_constraint_t::INEQUALITY);
        linear_constraint_t geq(-(cst.expression()),
                                linear_constraint_t::INEQUALITY);
        if (refuted(leq.negate()) && refuted(geq.negate())) {
          return true;
        }
      }

      tvpi_dbm_domain_t tmp = *this;
      tmp.normalize();
      CRAB_LOG("tvpi-dbm-entails", crab::outs()
                                       << "normalized: " << tmp << "\n");

      if (tmp.m_absval.entails(cst)) {
        return true;
      }

      if (ecst && !(*ecst).equal(cst)) {
        if (tmp.m_absval.entails(*ecst)) {
          return true;
        }
      }

      return false;
    }
  }

  /**
   * @brief Shared body of @ref assign and @ref weak_assign: @c x := e, then
   * re-establish the ghost invariant per template coefficient @c c —
   * @c ghost(x,c) := c-scaled(e) when expressible, else that ghost is
   * forgotten.  Strong assignments additionally add the structural ghost
   * equality (@ref rewrite_assign, e.g. @c x = b*y → @c x == ghost(y,b));
   * a weak assignment must not (it could constrain @c x beyond weak
   * semantics), and it weak-assigns the ghosts instead.
   */
  void assign_impl(const variable_t &x, const linear_expression_t &e,
                   bool weak) {
    if (is_bottom())
      return;
    CRAB_LOG("tvpi-dbm-assign",
             crab::outs() << "Before " << (weak ? "weak " : "") << "assign("
                          << x << " := " << e << ")=" << *this << "\n");
    m_vars.insert(x);
    for (auto const &v : e.variables()) {
      m_vars.insert(v);
    }
    if (weak) {
      m_absval.weak_assign(x, e);
    } else {
      m_absval.assign(x, e);
      rewrite_assign(x, e);
    }
    // Ascending order required: when e mentions x itself, ghost(x,c) :=
    // c-scaled(e) reads ghost(x, |b|*c) with |b| >= 2, so larger
    // coefficients must still hold their PRE-state.  Sortedness is an
    // invariant of fixed_tvpi_domain_params.
    const auto &coeffs = crab_domain_params_man::get().coefficients();
    assert(std::is_sorted(coeffs.begin(), coeffs.end()));
    for (auto c : coeffs) {
      auto cx = get_ghost_var(x, c);
      if (auto ce = try_rewrite_linear_expression(e, c)) {
        CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                        << "processing rewritten " << cx
                                        << " := " << *ce << "\n");
        if (weak) {
          m_absval.weak_assign(cx, *ce);
        } else {
          m_absval.assign(cx, *ce);
        }
      } else {
        m_absval -= cx;
        CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                        << "cannot rewrite: " << c << " * " << x
                                        << " := " << c << " * "
                                        << "(" << e << ")"
                                        << "\n");
      }
    }
    CRAB_LOG("tvpi-dbm-assign", crab::outs()
                                    << "After assign=" << *this << "\n");
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".assign");
    assign_impl(x, e, false /*weak*/);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".weak_assign");
    assign_impl(x, e, true /*weak*/);
  }

  /**
   * @brief @c x := y op z (scalar @p z).  Precision special cases:
   *  - division by ±1 → routed through @ref assign (keeps @c x - y = 0);
   *  - @c x = y*z, |z|>1 → adds @c x == ±ghost(y,|z|);
   *  - @c x = y/z, |z|>1 → assigns @c ghost(x,|z|) := ±y;
   * then @ref rewrite_apply per tracked coefficient, plus incremental
   * saturation for the new TVPI equalities.
   */
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".apply");
    if (!is_bottom()) {
      m_vars.insert(x);
      m_vars.insert(y);
      // Division by ±1 is exact: route through assign so that sdbm adds
      // difference constraints (x - y = 0) rather than losing them.
      if ((op == OP_SDIV || op == OP_UDIV) &&
          (z == number_t(1) || z == number_t(-1))) {
        linear_expression_t rhs = (z == number_t(1)) ? linear_expression_t(y)
                                                     : linear_expression_t(-y);
        assign(x, rhs);
        return;
      }
      m_absval.apply(op, x, y, z);
      bool neg = (z < 0);
      number_t sign_one = neg ? number_t(-1) : number_t(1);
      number_t z_abs = neg ? -z : z;

      // rewrite("x := y op z")
      switch (op) {
      case OP_ADDITION:
      case OP_SUBTRACTION:
        break;
      case OP_MULTIPLICATION: // x := y * z
        if (z_abs > number_t(1)) {
          // TVPI form x == ±z*y over original variables, passed to
          // incremental saturation (a negative sign gives a sum, which the
          // hook rejects).
          const linear_constraint_t tvpi_cst(linear_expression_t(x) -
                                                 (sign_one * z_abs) * y,
                                             linear_constraint_t::EQUALITY);
          // Add constraint x == ±ghost(y, z_abs) without overwriting x's
          // precise value
          if (auto gv_opt = find_ghost_var(y, z_abs)) {
            linear_expression_t ghost_e = sign_one > 0
                                              ? linear_expression_t(*gv_opt)
                                              : linear_expression_t(-(*gv_opt));
            m_absval += linear_constraint_t(linear_expression_t(x) - ghost_e,
                                            linear_constraint_t::EQUALITY);
            incremental_reduce_constraint(tvpi_cst);
          }
        }
        break;
      case OP_SDIV:                // x := y /s z
      case OP_UDIV:                // x := y /u z
        if (z_abs > number_t(1)) { // "zx := y for z > 1"
          // Mirror of the ghost update below:  z*x == ±y.
          const linear_constraint_t tvpi_cst(z_abs * x - sign_one * y,
                                             linear_constraint_t::EQUALITY);
          if (find_ghost_var(x, z_abs)) {
            m_absval.apply(op, get_ghost_var(x, z_abs), y, sign_one);
            incremental_reduce_constraint(tvpi_cst);
          }
        }
        break;
      case OP_SREM: // x := y %/s z
      case OP_UREM: // x := y %/u z
        // This can be rewritten into b - c * floor(b / c) if
        // b / c in low level domain is equivalent to floor(b / c).
        break;
      default:
        break;
      }

      auto &coeffs = crab_domain_params_man::get().coefficients();
      for (auto coefficient : coeffs) {
        rewrite_apply(op, x, y, z, coefficient);
      }
    }
  }

  /// Reduce @c x := y op z to scalar form when one operand is a singleton
  /// interval:  z = c  gives  x := y op c  for ANY op;  y = c  gives
  /// x := z op c  only for commutative ops (c - z and c / z are not
  /// expressible as z op c — those stay with the base domain).
  void eval_apply(arith_operation_t op, const variable_t &x,
                  const variable_t &y, const variable_t &z) {
    if (auto c = m_absval[z].singleton()) { // x := y op c
      apply(op, x, y, *c);
    } else if (auto c = m_absval[y].singleton()) {
      switch (op) {
      case OP_ADDITION:
      case OP_MULTIPLICATION: // commutative: x := c op z == z op c
        apply(op, x, z, *c);
        break;
      default: // x := c - z, c / z, c % z: not expressible as z op c
        break;
      }
    }
  }

  /**
   * @brief Apply an arithmetic operation with a variable RHS: @c x := y op z.
   *
   * First delegates to the base domain, then calls @ref eval_apply to attempt
   * reduction to a scalar form.  Finally calls @ref rewrite_apply for each
   * tracked coefficient.
   */
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".applyv");
    if (!is_bottom()) {
      m_vars.insert(x);
      m_vars.insert(y);
      m_vars.insert(z);
      m_absval.apply(op, x, y, z);
      eval_apply(op, x, y, z);
      auto &coeffs = crab_domain_params_man::get().coefficients();
      for (auto coefficient : coeffs) {
        rewrite_apply(op, x, y, z, coefficient);
      }
    }
  }

  /// Integer conversion (truncation / zero-/sign-extension) @c dst := op(src).
  /// Truncation is not linear, so @c dst's ghosts cannot be re-established —
  /// they are dropped.
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom()) {
      m_vars.insert(dst);
      m_vars.insert(src);
      m_absval.apply(op, dst, src);
      forget_ghost_vars(dst);
    }
  }

  /// Bitwise operation @c x := y op z (variable RHS).  Not linear: @c x's
  /// ghosts are dropped rather than left stale.
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      m_vars.insert(x);
      m_vars.insert(y);
      m_vars.insert(z);
      m_absval.apply(op, x, y, z);
      forget_ghost_vars(x);
    }
  }

  /// Bitwise operation @c x := y op z (scalar RHS).  Not linear: @c x's
  /// ghosts are dropped rather than left stale.
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    if (!is_bottom()) {
      m_vars.insert(x);
      m_vars.insert(y);
      m_absval.apply(op, x, y, z);
      forget_ghost_vars(x);
    }
  }

  // ============================================================
  // Backward and inter-procedural operations
  // ============================================================

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_assign not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void callee_entry(const callsite_info<variable_t> &callsite,
                    const tvpi_dbm_domain_t &caller) override {
    inter_abstract_operations<
        tvpi_dbm_domain_t,
        Params::implement_inter_transformers>::callee_entry(callsite, caller,
                                                            *this);
  }

  void caller_continuation(const callsite_info<variable_t> &callsite,
                           const tvpi_dbm_domain_t &callee) override {
    inter_abstract_operations<
        tvpi_dbm_domain_t,
        Params::implement_inter_transformers>::caller_continuation(callsite,
                                                                   callee,
                                                                   *this);
  }

  // ============================================================
  // Conversions, projection, renaming
  // ============================================================

  /// Constraints over original variables only.  @c m_absval holds both
  /// originals and ghosts; every ghost term is translated back through the
  /// ghost invariant  ghost(v,k) = k*v, e.g. the base constraint
  /// ghost(x,3) - y <= 4 is returned as  3*x - y <= 4.  Constraints
  /// mentioning a dimension that is neither tracked nor a template ghost
  /// are dropped (conservative).
  linear_constraint_system_t to_linear_constraint_system() const override {

    if (is_bottom()) {
      return linear_constraint_system_t(linear_constraint_t::get_false());
    }

    if (is_top()) {
      return linear_constraint_system_t(linear_constraint_t::get_true());
    }

    // ghost variable -> (original variable, coefficient), for every tracked
    // original and template coefficient.
    std::unordered_map<variable_t, std::pair<variable_t, number_t>> g2o;
    for (const auto &v : m_vars) {
      for (unsigned k : crab_domain_params_man::get().coefficients()) {
        if (k <= 1)
          continue;
        if (auto gv = find_ghost_var(v, number_t(k))) {
          g2o.emplace(*gv, std::make_pair(v, number_t(k)));
        }
      }
    }

    linear_constraint_system_t res;
    for (auto const &cst : m_absval.to_linear_constraint_system()) {
      const linear_expression_t &e = cst.expression();
      linear_expression_t e2(e.constant());
      bool ok = true;
      for (auto it = e.begin(), et = e.end(); it != et; ++it) {
        const number_t &coef = (*it).first;
        const variable_t &var = (*it).second;
        if (m_vars.find(var) != m_vars.end()) {
          e2 = e2 + coef * var;
        } else {
          auto g = g2o.find(var);
          if (g == g2o.end()) { // unknown dimension: drop conservatively
            ok = false;
            break;
          }
          // ghost(v,k) = k*v:  coef * ghost(v,k)  ==  (coef*k) * v.
          e2 = e2 + (coef * g->second.second) * g->second.first;
        }
      }
      if (ok) {
        res += linear_constraint_t(e2, cst.kind());
      }
    }
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_WARN(domain_name(),
              "::to_disjunctive_linear_constraint_system not implemented");
    disjunctive_linear_constraint_system_t res;
    return res;
  }

  /// Forget @p var and all its ghosts.  The eager variant normalizes first
  /// so constraints implied through @p var are explicit before it is dropped.
  void operator-=(const variable_t &var) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".-=");
    if (!(is_bottom() || is_top())) {
      if constexpr (Params::normalize)
        normalize();
      m_absval -= var;
      m_vars.erase(var);
      forget_ghost_vars(var);
    }
  }

  /// Interval of @p v.
  interval_t operator[](const variable_t &v) override {
    if (is_bottom()) {
      return interval_t::bottom();
    }
    return m_absval[v];
  }

  /// Interval of @p v (const, no side effects).
  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    }
    return m_absval.at(v);
  }

  /// Batched @ref operator-=.
  void forget(const variable_vector_t &variables) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".forget");
    for (auto const &v : variables) {
      *this -= v;
    }
  }

  /// Project onto @p variables plus their ghosts; everything else is dropped
  /// from @c m_absval, @c m_vars and the coefficient map.
  void project(const variable_vector_t &variables) override {
    TVPI_DBM_DOMAIN_SCOPED_STATS(".project");
    if (!is_bottom()) {
      variable_vector_t allvars(variables);
      for (auto const &v : variables) {
        auto &coeffs = crab_domain_params_man::get().coefficients();
        for (auto coefficient : coeffs) {
          variable_t gv = get_ghost_var(v, coefficient);
          allvars.push_back(gv);
        }
      }
      m_absval.project(allvars);
      m_vars = variable_set_t(variables.begin(), variables.end());
    }
  }

  /// Rename @c from[i] → to[i] pairwise, extending both lists with the
  /// corresponding ghost variables.
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (!is_bottom()) {
      variable_vector_t extd_from(from);
      variable_vector_t extd_to(to);
      for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
        const variable_t &f = from[i];
        const variable_t &t = to[i];
        if (m_vars.find(f) != m_vars.end()) {
          m_vars.erase(f);
          m_vars.insert(t);
        }
        auto &coeffs = crab_domain_params_man::get().coefficients();
        for (auto coefficient : coeffs) {
          extd_from.push_back(get_ghost_var(f, coefficient));
          extd_to.push_back(get_ghost_var(t, coefficient));
        }
      }
      m_absval.rename(extd_from, extd_to);
    }
  }

  /// Duplicate @p var (and its ghosts) as @p new_var — inter-procedural use.
  void expand(const variable_t &var, const variable_t &new_var) override {
    if (is_bottom() || is_top()) {
      return;
    }

    if (m_vars.find(var) != m_vars.end()) {
      m_vars.insert(new_var);
    }
    m_absval.expand(var, new_var);
    auto &coeffs = crab_domain_params_man::get().coefficients();
    for (auto coefficient : coeffs) {
      variable_t gv = get_ghost_var(var, coefficient);
      variable_t gnv = get_ghost_var(new_var, coefficient);
      m_absval.expand(gv, gnv);
    }
  }

  // ============================================================
  // Normalization and printing
  // ============================================================

  /**
   * @brief Normalize to the closed tDBM form:
   *  1. @c m_absval.normalize(): close the base DBM, so all constraint
   *     pairs with equal coefficients on the shared variable are combined;
   *  2. @ref tvpi_reduce(): apply the Resultant rule to the remaining,
   *     coefficient-misaligned pairs; every derived constraint is closed by
   *     the base domain as it is added.
   *
   * This realizes DbmTvpiSaturation, which alternates TvpiReduce and
   * DbmClosure for ⌈log₂(N)⌉-1 rounds: the rounds run inside
   * @ref tvpi_reduce, and DbmClosure happens at every insertion rather than
   * once per round.  Called automatically after each operation in the eager
   * variant; explicitly by the client otherwise.
   */
  void normalize() override {
    if (is_bottom() || is_top())
      return;
    m_absval.normalize();
    if (!is_bottom())
      tvpi_reduce();
  }
  void minimize() override {}

  void intrinsic(std::string /*name*/,
                 const variable_or_constant_vector_t & /*inputs*/,
                 const variable_vector_t & /*outputs*/) override {}

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const tvpi_dbm_domain_t &invariant) override {
    CRAB_WARN(domain_name(), "::backward_intrinsic for ", name,
              " not implemented");
  }

  /// Print as @c {coeffs:<set>, absval:<base-state>} (ghosts visible; use
  /// @ref to_linear_constraint_system for originals only).
  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      o << "{";
      {
        o << "coeffs:[";
        const auto &coeffs = crab_domain_params_man::get().coefficients();
        for (auto it = coeffs.begin(); it != coeffs.end(); ++it) {
          if (it != coeffs.begin()) {
            o << ",";
          }
          o << (*it);
        }
        o << "], ";
      }
      o << "absval:";
      m_absval.write(o);
      o << "}";
    }
  }

  friend crab_os &operator<<(crab_os &o, const tvpi_dbm_domain_t &val) {
    val.write(o);
    return o;
  }

  /// Return the domain name string, e.g. @c "TVPI(SplitDBM)".
  std::string domain_name() const override {
    base_domain_t absval;
    std::string base_name = absval.domain_name();
    const char *prefix = "TVPI";
    std::string name;
    name.reserve(base_name.size() + 8);
    name.append(prefix);
    name.append("(");
    name.append(base_name);
    name.append(")");
    return name;
  }
};

template <typename BaseDom, typename Params>
struct abstract_domain_traits<tvpi_dbm_domain<BaseDom, Params>> {
  using number_t = typename BaseDom::number_t;
  using varname_t = typename BaseDom::varname_t;
};

} // end namespace domains
} // end namespace crab
