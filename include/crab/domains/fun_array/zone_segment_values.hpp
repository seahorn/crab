#pragma once

#include <crab/domains/fun_array/segment_value_api.hpp>
#include <crab/domains/split_dbm.hpp>
#include <crab/types/variable.hpp>

#include <algorithm>
#include <cassert>
#include <utility>

namespace crab {
namespace domains {
namespace fun_array_impl {

/**
 * Segment values: universally quantified invariants over one index.
 *
 * The value of a segment [l, h) is a zones state (split_dbm, in split
 * normal form) that stands for
 *
 *   forall j in [l, h). phi(A[j], B[j], ..., j, scalars)
 *
 * over
 *  - the array variables of the group: in phi, A stands for A[j], the
 *    element at index j;
 *  - the index placeholder (m_index), which stands for j itself;
 *  - scalars brought in by loads, stores and singleton.
 *
 * For instance:
 *   A == 0      forall j in [l, h). A[j] == 0
 *   A == B      forall j in [l, h). A[j] == B[j]
 *   A == j      forall j in [l, h). A[j] == j   (needs the index placeholder)
 *   A <= t      forall j in [l, h). A[j] <= t, for a scalar t
 *
 * Only one index is quantified: invariants relating two elements of the
 * same segment, such as sortedness (forall j1 < j2. A[j1] <= A[j2]), need a
 * different value domain.
 *
 * Array variables never reach the scalar domain, so using them in phi
 * cannot clash with scalar constraints.
 **/
template <class ScalarDom>
class zone_segment_values final
    : public segment_value_api<zone_segment_values<ScalarDom>, ScalarDom> {
public:
  using zone_segment_values_t = zone_segment_values<ScalarDom>;
  using number_t = typename ScalarDom::number_t;
  using varname_t = typename ScalarDom::varname_t;
  using variable_t = typename ScalarDom::variable_t;
  using variable_vector_t = typename ScalarDom::variable_vector_t;
  using linear_expression_t = typename ScalarDom::linear_expression_t;
  using linear_constraint_t = typename ScalarDom::linear_constraint_t;
  using linear_constraint_system_t =
      typename ScalarDom::linear_constraint_system_t;
  using zones_t = split_dbm_domain<number_t, varname_t>;

private:
  // The quantified index j. Only needed to relate an element to its own
  // index (e.g., A[j] == j); invariants such as A[j] == 0 or A[j] == B[j]
  // do not mention it. It is a ghost variable named after the group's
  // representative array, so every segment of the group shares it.
  variable_t m_index;
  zones_t m_state;

  zone_segment_values(variable_t index, zones_t state)
      : m_index(index), m_state(std::move(state)) {}

  zone_segment_values_t with_state(const zone_segment_values_t &o,
                                   zones_t state) const {
    assert(m_index == o.m_index);
    return zone_segment_values_t(m_index, std::move(state));
  }

  // The scalars tracked by m_state: neither arrays nor the index
  variable_vector_t scalars(const zones_t &state) const {
    variable_vector_t res;
    for (auto const &v : state.tracked_variables()) {
      if (!v.get_type().is_array() && !(v == m_index)) {
        res.push_back(v);
      }
    }
    return res;
  }

  // What s knows about vars
  static linear_constraint_system_t project(const ScalarDom &s,
                                            const variable_vector_t &vars) {
    ScalarDom t(s);
    t.project(vars);
    return t.to_linear_constraint_system();
  }

public:
  // The index placeholder of a group, named after its representative array
  static variable_t make_index_placeholder(const variable_t &representative) {
    auto &vfac =
        const_cast<varname_t *>(&(representative.name()))->get_var_factory();
    return variable_t(vfac.get_or_insert_varname(representative.name(), "$j"),
                      variable_type(MATH_INT_TYPE));
  }

  // Top
  explicit zone_segment_values(const variable_t &index) : m_index(index) {}

  // Top, for a group whose representative array is a
  static zone_segment_values_t make_top_for(const variable_t &a) {
    return zone_segment_values_t(make_index_placeholder(a));
  }

  zone_segment_values_t make_bottom() const override {
    return zone_segment_values_t(m_index, zones_t(true));
  }
  zone_segment_values_t make_top() const override {
    return zone_segment_values_t(m_index);
  }
  void set_to_bottom() override { m_state.set_to_bottom(); }
  void set_to_top() override { m_state.set_to_top(); }
  bool is_bottom() const override { return m_state.is_bottom(); }
  bool is_top() const override { return m_state.is_top(); }

  bool operator<=(const zone_segment_values_t &o) const override {
    return m_state <= o.m_state;
  }
  zone_segment_values_t
  operator|(const zone_segment_values_t &o) const override {
    return with_state(o, m_state | o.m_state);
  }
  zone_segment_values_t
  operator&(const zone_segment_values_t &o) const override {
    return with_state(o, m_state & o.m_state);
  }
  zone_segment_values_t
  operator||(const zone_segment_values_t &o) const override {
    return with_state(o, m_state || o.m_state);
  }
  zone_segment_values_t
  operator&&(const zone_segment_values_t &o) const override {
    return with_state(o, m_state && o.m_state);
  }

  // What holds for every index of [l,h) holds for every index of its parts
  std::pair<zone_segment_values_t, zone_segment_values_t>
  split() const override {
    return {*this, *this};
  }
  zone_segment_values_t
  merge(const zone_segment_values_t &right) const override {
    return *this | right;
  }
  zone_segment_values_t
  singleton(const linear_expression_t &i) const override {
    zone_segment_values_t r(*this);
    r.m_state += linear_constraint_t(linear_expression_t(m_index) - i,
                                     linear_constraint_t::EQUALITY);
    return r;
  }

  void store(const variable_t &a, const linear_expression_t &val,
             const ScalarDom &s) override {
    m_state.assign(a, val);
    variable_vector_t vars(val.variables_begin(), val.variables_end());
    if (!vars.empty()) {
      m_state += project(s, vars);
    }
  }

  void load(const variable_t &a, const variable_t &lhs,
            ScalarDom &s) override {
    zones_t t(m_state);
    t.assign(lhs, a);
    // Projecting on the scalars drops the arrays and the index
    t.project(scalars(t));
    s += t.to_linear_constraint_system();

    m_state += linear_constraint_t(linear_expression_t(a) - lhs,
                                   linear_constraint_t::EQUALITY);
  }

  bool mentions(const variable_t &x) const override {
    return m_state.is_tracked(x);
  }
  void forget_array(const variable_t &a) override { m_state -= a; }

  void forget(const variable_vector_t &scalars) override {
    m_state.forget(scalars);
  }
  void project(const variable_vector_t &keep) override {
    variable_vector_t vars;
    for (auto const &v : m_state.tracked_variables()) {
      if (v.get_type().is_array() || v == m_index ||
          std::find(keep.begin(), keep.end(), v) != keep.end()) {
        vars.push_back(v);
      }
    }
    m_state.project(vars);
  }
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    m_state.rename(from, to);
  }
  void reduce(const ScalarDom &s) override {
    variable_vector_t vars = scalars(m_state);
    if (!vars.empty()) {
      m_state += project(s, vars);
    }
  }

  void write(crab_os &o) const override { m_state.write(o); }
};

template <class ScalarDom>
inline crab_os &operator<<(crab_os &o,
                           const zone_segment_values<ScalarDom> &v) {
  v.write(o);
  return o;
}

} // end namespace fun_array_impl
} // end namespace domains
} // end namespace crab
