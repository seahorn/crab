#pragma once

#include <crab/domains/fun_array/segment_bound_api.hpp>
#include <crab/domains/fun_array/segment_expression_api.hpp>

#include <boost/optional.hpp>
#include <set>
#include <type_traits>
#include <utility>
#include <vector>

namespace crab {
namespace domains {
namespace fun_array_impl {

/**
 * Default bound (POPL'11, §8): a finite set of expressions, all equal in
 * the concrete. Written once for any expression domain Expr.
 **/
template <class Expr, class ScalarDom>
class expression_set_bound final
    : public segment_bound_api<expression_set_bound<Expr, ScalarDom>, Expr,
                               ScalarDom> {
public:
  using expression_set_bound_t = expression_set_bound<Expr, ScalarDom>;
  using expression_t = Expr;
  using variable_t = typename ScalarDom::variable_t;
  using variable_vector_t = typename ScalarDom::variable_vector_t;
  using linear_expression_t = typename ScalarDom::linear_expression_t;
  using number_t = typename ScalarDom::number_t;

  static_assert(
      std::is_base_of<segment_expression_api<Expr, ScalarDom>, Expr>::value,
      "Expr must derive from segment_expression_api");
  static_assert(
      std::is_same<decltype(Expr::make(
                       std::declval<const linear_expression_t &>())),
                   boost::optional<Expr>>::value,
      "Expr must have static boost::optional<Expr> "
      "make(const linear_expression_t &)");

private:
  // Ordered by Expr::operator<
  std::set<Expr> m_exprs;

  explicit expression_set_bound(std::set<Expr> &&exprs)
      : m_exprs(std::move(exprs)) {}

  // Some pair (e, e') with e in this and e' in o satisfies cmp
  template <class Cmp>
  bool some_pair(const expression_set_bound_t &o, Cmp cmp) const {
    for (auto const &e : m_exprs) {
      for (auto const &f : o.m_exprs) {
        if (cmp(e, f)) {
          return true;
        }
      }
    }
    return false;
  }

public:
  expression_set_bound() {}

  explicit expression_set_bound(const Expr &e) { m_exprs.insert(e); }

  bool is_empty() const override { return m_exprs.empty(); }

  bool contains(const Expr &e) const override {
    return m_exprs.count(e) > 0;
  }

  std::vector<Expr> expressions() const override {
    return std::vector<Expr>(m_exprs.begin(), m_exprs.end());
  }

  std::vector<variable_t> variables() const override {
    std::set<variable_t> vars;
    for (auto const &e : m_exprs) {
      if (auto x = e.variable()) {
        vars.insert(*x);
      }
    }
    return std::vector<variable_t>(vars.begin(), vars.end());
  }

  boost::optional<number_t> find_constant() const override {
    for (auto const &e : m_exprs) {
      if (auto k = e.constant()) {
        return k;
      }
    }
    return boost::none;
  }

  expression_set_bound_t keep_only_constants() const override {
    std::set<Expr> res;
    for (auto const &e : m_exprs) {
      if (e.is_constant()) {
        res.insert(e);
      }
    }
    return expression_set_bound_t(std::move(res));
  }

  bool is_successor_of(const expression_set_bound_t &o) const override {
    for (auto const &e : o.m_exprs) {
      if (contains(e.shifted(number_t(1)))) {
        return true;
      }
    }
    return false;
  }

  expression_set_bound_t
  operator&(const expression_set_bound_t &o) const override {
    std::set<Expr> res;
    for (auto const &e : m_exprs) {
      if (o.contains(e)) {
        res.insert(e);
      }
    }
    return expression_set_bound_t(std::move(res));
  }

  expression_set_bound_t
  operator|(const expression_set_bound_t &o) const override {
    std::set<Expr> res(m_exprs);
    res.insert(o.m_exprs.begin(), o.m_exprs.end());
    return expression_set_bound_t(std::move(res));
  }

  expression_set_bound_t
  operator-(const expression_set_bound_t &o) const override {
    std::set<Expr> res;
    for (auto const &e : m_exprs) {
      if (!o.contains(e)) {
        res.insert(e);
      }
    }
    return expression_set_bound_t(std::move(res));
  }

  bool operator<=(const expression_set_bound_t &o) const override {
    for (auto const &e : m_exprs) {
      if (!o.contains(e)) {
        return false;
      }
    }
    return true;
  }

  expression_set_bound_t shifted(number_t k) const override {
    std::set<Expr> res;
    for (auto const &e : m_exprs) {
      res.insert(e.shifted(k));
    }
    return expression_set_bound_t(std::move(res));
  }

  void add(const Expr &e) override { m_exprs.insert(e); }

  void forget(const variable_t &x) override {
    for (auto it = m_exprs.begin(); it != m_exprs.end();) {
      if (it->mentions(x)) {
        it = m_exprs.erase(it);
      } else {
        ++it;
      }
    }
  }

  void substitute(const variable_t &x, const linear_expression_t &e) override {
    std::set<Expr> res;
    for (auto const &f : m_exprs) {
      if (auto g = f.substitute(x, e)) {
        res.insert(*g);
      }
    }
    m_exprs = std::move(res);
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    std::set<Expr> res;
    for (auto const &e : m_exprs) {
      res.insert(e.rename(from, to));
    }
    m_exprs = std::move(res);
  }

  // A true pair decides the comparison before a false one: both exist only
  // if s is bottom.
  boolean_value is_eq(const expression_set_bound_t &o,
                      const ScalarDom &s) const override {
    if (some_pair(o, [&s](const Expr &e, const Expr &f) {
          return e.is_eq(f, s).is_true();
        })) {
      return boolean_value::get_true();
    }
    if (some_pair(o, [&s](const Expr &e, const Expr &f) {
          return e.is_eq(f, s).is_false();
        })) {
      return boolean_value::get_false();
    }
    return boolean_value::top();
  }

  boolean_value is_lt(const expression_set_bound_t &o,
                      const ScalarDom &s) const override {
    if (some_pair(o, [&s](const Expr &e, const Expr &f) {
          return e.is_lt(f, s).is_true();
        })) {
      return boolean_value::get_true();
    }
    if (some_pair(o, [&s](const Expr &e, const Expr &f) {
          return e.is_lt(f, s).is_false();
        })) {
      return boolean_value::get_false();
    }
    return boolean_value::top();
  }

  boolean_value is_le(const expression_set_bound_t &o,
                      const ScalarDom &s) const override {
    if (some_pair(o, [&s](const Expr &e, const Expr &f) {
          return e.is_le(f, s).is_true();
        })) {
      return boolean_value::get_true();
    }
    if (some_pair(o, [&s](const Expr &e, const Expr &f) {
          return e.is_le(f, s).is_false();
        })) {
      return boolean_value::get_false();
    }
    return boolean_value::top();
  }

  // {e1, e2, ...}
  void write(crab_os &o) const override {
    o << "{";
    for (auto it = m_exprs.begin(); it != m_exprs.end(); ++it) {
      if (it != m_exprs.begin()) {
        o << ", ";
      }
      o << *it;
    }
    o << "}";
  }
};

template <class Expr, class ScalarDom>
inline crab_os &operator<<(crab_os &o,
                           const expression_set_bound<Expr, ScalarDom> &b) {
  b.write(o);
  return o;
}

} // end namespace fun_array_impl
} // end namespace domains
} // end namespace crab
