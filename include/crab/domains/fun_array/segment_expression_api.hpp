#pragma once

#include <crab/domains/boolean.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>

namespace crab {
namespace domains {
namespace fun_array_impl {

/**
 * Paper's E (POPL'11, §7): one bound expression in normal form.
 *
 * Use of Curiously Recurring Template Pattern (CRTP): Expr is the
 * implementation, ScalarDom the scalar domain consulted by the
 * comparisons.
 *
 * Required static member of Expr (static members cannot be virtual):
 *
 *   static boost::optional<Expr> make(const linear_expression_t &e);
 *
 * It returns e in normal form, or none if e cannot be put in normal form.
 **/
template <class Expr, class ScalarDom> class segment_expression_api {
public:
  using variable_t = typename ScalarDom::variable_t;
  using variable_vector_t = typename ScalarDom::variable_vector_t;
  using linear_expression_t = typename ScalarDom::linear_expression_t;
  using number_t = typename ScalarDom::number_t;

  virtual ~segment_expression_api() {}

  // For queries to the scalar domain and for exporting constraints
  virtual linear_expression_t to_linear_expression() const = 0;

  // Shape queries used by unification heuristics and purification
  virtual bool is_constant() const = 0;
  // Defined only if is_constant()
  virtual boost::optional<number_t> constant() const = 0;
  // The x in x + k
  virtual boost::optional<variable_t> variable() const = 0;

  virtual bool mentions(const variable_t &x) const = 0;
  // x := e inside this expression; none if the result is not in normal form
  virtual boost::optional<Expr>
  substitute(const variable_t &x, const linear_expression_t &e) const = 0;
  virtual Expr rename(const variable_vector_t &from,
                      const variable_vector_t &to) const = 0;
  // this + k
  virtual Expr shifted(number_t k) const = 0;

  // Three-valued comparisons: syntactic first, then the scalar domain (§7.4)
  virtual boolean_value is_eq(const Expr &o, const ScalarDom &s) const = 0;
  virtual boolean_value is_lt(const Expr &o, const ScalarDom &s) const = 0;
  virtual boolean_value is_le(const Expr &o, const ScalarDom &s) const = 0;

  // Syntactic equality and a total order, so bounds can store expressions
  // in sorted containers
  virtual bool operator==(const Expr &o) const = 0;
  virtual bool operator<(const Expr &o) const = 0;

  virtual void write(crab_os &o) const = 0;
};

} // end namespace fun_array_impl
} // end namespace domains
} // end namespace crab
