#pragma once

#include <crab/domains/fun_array/segment_expression_api.hpp>
#include <crab/numbers/bignums.hpp>

#include <boost/optional.hpp>
#include <cassert>
#include <type_traits>

namespace crab {
namespace domains {
namespace fun_array_impl {

/**
 * Default bound expression (POPL'11, §7.2): x + k or a constant k.
 *
 * Indices are integers: e < e' is checked as e - e' <= -1.
 **/
template <class ScalarDom>
class offset_expression final
    : public segment_expression_api<offset_expression<ScalarDom>, ScalarDom> {
public:
  using offset_expression_t = offset_expression<ScalarDom>;
  using variable_t = typename ScalarDom::variable_t;
  using variable_vector_t = typename ScalarDom::variable_vector_t;
  using linear_expression_t = typename ScalarDom::linear_expression_t;
  using linear_constraint_t = typename ScalarDom::linear_constraint_t;
  using number_t = typename ScalarDom::number_t;

  static_assert(std::is_same<number_t, ikos::z_number>::value,
                "offset_expression requires integer numbers");

private:
  boost::optional<variable_t> m_var;
  number_t m_k;

  offset_expression(boost::optional<variable_t> x, number_t k)
      : m_var(x), m_k(k) {}

  // Same variable, or both constants: the comparison is decided by the
  // offsets alone.
  bool same_base(const offset_expression_t &o) const { return m_var == o.m_var; }

public:
  static boost::optional<offset_expression_t>
  make(const linear_expression_t &e) {
    if (e.is_constant()) {
      return offset_expression_t(boost::none, e.constant());
    }
    if (e.size() == 1) {
      auto it = e.begin();
      if (it->first == number_t(1)) {
        return offset_expression_t(it->second, e.constant());
      }
    }
    return boost::none;
  }

  linear_expression_t to_linear_expression() const override {
    if (m_var) {
      return linear_expression_t(*m_var) + m_k;
    } else {
      return linear_expression_t(m_k);
    }
  }

  bool is_constant() const override { return !m_var; }

  boost::optional<number_t> constant() const override {
    if (m_var) {
      return boost::none;
    }
    return m_k;
  }

  boost::optional<variable_t> variable() const override { return m_var; }

  bool mentions(const variable_t &x) const override {
    return m_var && *m_var == x;
  }

  boost::optional<offset_expression_t>
  substitute(const variable_t &x, const linear_expression_t &e) const override {
    if (!mentions(x)) {
      return *this;
    }
    return make(e + m_k);
  }

  offset_expression_t rename(const variable_vector_t &from,
                             const variable_vector_t &to) const override {
    assert(from.size() == to.size());
    if (m_var) {
      for (unsigned i = 0, n = from.size(); i < n; ++i) {
        if (from[i] == *m_var) {
          return offset_expression_t(to[i], m_k);
        }
      }
    }
    return *this;
  }

  offset_expression_t shifted(number_t k) const override {
    return offset_expression_t(m_var, m_k + k);
  }

  boolean_value is_eq(const offset_expression_t &o,
                      const ScalarDom &s) const override {
    if (same_base(o)) {
      return m_k == o.m_k ? boolean_value::get_true()
                          : boolean_value::get_false();
    }
    linear_expression_t e = to_linear_expression();
    linear_expression_t f = o.to_linear_expression();
    if (s.entails(linear_constraint_t(e - f, linear_constraint_t::EQUALITY))) {
      return boolean_value::get_true();
    }
    if (s.entails(e - f <= number_t(-1)) || s.entails(f - e <= number_t(-1))) {
      return boolean_value::get_false();
    }
    return boolean_value::top();
  }

  boolean_value is_lt(const offset_expression_t &o,
                      const ScalarDom &s) const override {
    if (same_base(o)) {
      return m_k < o.m_k ? boolean_value::get_true()
                         : boolean_value::get_false();
    }
    linear_expression_t e = to_linear_expression();
    linear_expression_t f = o.to_linear_expression();
    if (s.entails(e - f <= number_t(-1))) {
      return boolean_value::get_true();
    }
    if (s.entails(f - e <= number_t(0))) {
      return boolean_value::get_false();
    }
    return boolean_value::top();
  }

  boolean_value is_le(const offset_expression_t &o,
                      const ScalarDom &s) const override {
    if (same_base(o)) {
      return m_k <= o.m_k ? boolean_value::get_true()
                          : boolean_value::get_false();
    }
    linear_expression_t e = to_linear_expression();
    linear_expression_t f = o.to_linear_expression();
    if (s.entails(e - f <= number_t(0))) {
      return boolean_value::get_true();
    }
    if (s.entails(f - e <= number_t(-1))) {
      return boolean_value::get_false();
    }
    return boolean_value::top();
  }

  bool operator==(const offset_expression_t &o) const override {
    return m_var == o.m_var && m_k == o.m_k;
  }

  // Constants first, then by variable, then by offset
  bool operator<(const offset_expression_t &o) const override {
    if (m_var != o.m_var) {
      if (!m_var || !o.m_var) {
        return !m_var;
      }
      return *m_var < *o.m_var;
    }
    return m_k < o.m_k;
  }

  void write(crab_os &o) const override {
    if (!m_var) {
      o << m_k;
    } else if (m_k == 0) {
      o << *m_var;
    } else if (m_k > 0) {
      o << *m_var << "+" << m_k;
    } else {
      o << *m_var << m_k;
    }
  }
};

template <class ScalarDom>
inline crab_os &operator<<(crab_os &o, const offset_expression<ScalarDom> &e) {
  e.write(o);
  return o;
}

} // end namespace fun_array_impl
} // end namespace domains
} // end namespace crab
