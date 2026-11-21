#pragma once

#include <crab/domains/boolean.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <vector>

namespace crab {
namespace domains {
namespace fun_array_impl {

/**
 * Paper's B(E) (POPL'11, §8): a set of expressions that are all equal in
 * the concrete.
 *
 * Use of Curiously Recurring Template Pattern (CRTP): Bound is the
 * implementation, Expr its expression domain (a segment_expression_api),
 * ScalarDom the scalar domain consulted by the comparisons.
 *
 * Required constructors of Bound (constructors cannot be virtual):
 *
 *   Bound();                       // empty set
 *   explicit Bound(const Expr &e); // singleton
 **/
template <class Bound, class Expr, class ScalarDom> class segment_bound_api {
public:
  using expression_t = Expr;
  using variable_t = typename ScalarDom::variable_t;
  using variable_vector_t = typename ScalarDom::variable_vector_t;
  using linear_expression_t = typename ScalarDom::linear_expression_t;
  using number_t = typename ScalarDom::number_t;

  virtual ~segment_bound_api() {}

  // The empty set denotes unreachability
  virtual bool is_empty() const = 0;
  virtual bool contains(const Expr &e) const = 0;
  // The expressions, sorted by Expr::operator<
  virtual std::vector<Expr> expressions() const = 0;
  // Variables mentioned by the expressions (used by purification)
  virtual std::vector<variable_t> variables() const = 0;

  // Used by Clousot's unification heuristics
  // (SegmentLimit.TryFindConstant)
  virtual boost::optional<number_t> find_constant() const = 0;
  // (SegmentLimit.KeepOnlyConstantExpressionIfAny)
  virtual Bound keep_only_constants() const = 0;
  // Some e in o has e + 1 in this bound (SegmentLimit.IsSuccessorOf)
  virtual bool is_successor_of(const Bound &o) const = 0;

  // Intersection
  virtual Bound operator&(const Bound &o) const = 0;
  // Union
  virtual Bound operator|(const Bound &o) const = 0;
  // Difference
  virtual Bound operator-(const Bound &o) const = 0;
  // Subset
  virtual bool operator<=(const Bound &o) const = 0;
  // Every e becomes e + k
  virtual Bound shifted(number_t k) const = 0;

  virtual void add(const Expr &e) = 0;
  // Drop the expressions mentioning x
  virtual void forget(const variable_t &x) = 0;
  // x := e inside every expression; drop those not in normal form after it
  virtual void substitute(const variable_t &x,
                          const linear_expression_t &e) = 0;
  virtual void rename(const variable_vector_t &from,
                      const variable_vector_t &to) = 0;

  // Three-valued comparisons of whole bounds (§8.2). E.g., B < B' is true if
  // some e in B and e' in B' satisfy e < e', false if some pair satisfies
  // e >= e', unknown otherwise.
  virtual boolean_value is_eq(const Bound &o, const ScalarDom &s) const = 0;
  virtual boolean_value is_lt(const Bound &o, const ScalarDom &s) const = 0;
  virtual boolean_value is_le(const Bound &o, const ScalarDom &s) const = 0;

  virtual void write(crab_os &o) const = 0;
};

} // end namespace fun_array_impl
} // end namespace domains
} // end namespace crab
