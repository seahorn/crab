#pragma once

#include <boost/optional.hpp>
#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/boolean.hpp>

#include <memory>

namespace crab {
namespace domains {
namespace object_domain_impl {

/* Copy-on-write wrapper around a concrete abstract domain: it holds one
   std::shared_ptr<Domain> and deep-copies it only when a mutating operation
   is requested on a shared value. Used by odi_map_domain for the summary,
   cache and field-equality components of an abstract object.
   Not to be confused with crab::domains::abstract_domain_ref
   (generic_abstract_domain.hpp), which wraps the *type-erased* domain and
   keeps a separate normalized value for the asc/desc phase mechanism. */
template <typename Variable, typename Domain>
class abstract_domain_ref final
    : abstract_domain_api<abstract_domain_ref<Variable, Domain>> {
public:
  using number_t = typename Variable::number_t;
  using varname_t = typename Variable::varname_t;
  using linear_expression_t = ikos::linear_expression<number_t, varname_t>;
  using linear_constraint_t = ikos::linear_constraint<number_t, varname_t>;
  using linear_constraint_system_t =
      ikos::linear_constraint_system<number_t, varname_t>;
  using disjunctive_linear_constraint_system_t =
      ikos::disjunctive_linear_constraint_system<number_t, varname_t>;
  using variable_t = variable<number_t, varname_t>;
  using variable_or_constant_t = variable_or_constant<number_t, varname_t>;
  using variable_vector_t = std::vector<variable_t>;
  using variable_or_constant_vector_t = std::vector<variable_or_constant_t>;
  using reference_constraint_t = reference_constraint<number_t, varname_t>;
  using interval_t = ikos::interval<number_t>;

private:
  using this_domain_t = abstract_domain_ref<Variable, Domain>;
  using abstract_domain_t = Domain;
  using abstract_domain_ptr_t = std::shared_ptr<abstract_domain_t>;

  abstract_domain_ptr_t m_base_ref;

  abstract_domain_ref(abstract_domain_ptr_t ref) : m_base_ref(ref) {}

  this_domain_t create(abstract_domain_t &&abs) const {
    return this_domain_t(std::make_shared<abstract_domain_t>(std::move(abs)));
  }

  void detach(void) { m_base_ref.reset(new abstract_domain_t(*m_base_ref)); }

  abstract_domain_t &detach_and_get_absval() {
    // NOTE: use_count()==1 instead of unique(): the latter is deprecated in
    // C++17 and removed in C++20.
    if (m_base_ref.use_count() != 1) {
      detach();
    }
    return *m_base_ref;
  }

  const abstract_domain_t &get_absval() const { return *m_base_ref; }

public:
  abstract_domain_ref()
      : m_base_ref(std::make_shared<abstract_domain_t>(abstract_domain_t())) {}

  ~abstract_domain_ref() = default;

  abstract_domain_ref(const this_domain_t &o) = default;

  this_domain_t &operator=(const this_domain_t &o) = default;

  abstract_domain_ref(this_domain_t &&o) = default;

  this_domain_t &operator=(this_domain_t &&o) = default;

  bool is_asc_phase() const override { return get_absval().is_asc_phase(); }

  void set_phase(bool is_asc) override {
    detach_and_get_absval().set_phase(is_asc);
  }

  this_domain_t make_top() const override {
    return create(get_absval().make_top());
  }

  this_domain_t make_bottom() const override {
    return create(get_absval().make_bottom());
  }

  // set_to_top/set_to_bottom never read the current value, so a shared
  // value is replaced with a fresh one instead of detach()'s deep copy
  // (the clone would be overwritten immediately).
  void set_to_top() override {
    if (m_base_ref.use_count() == 1) {
      m_base_ref->set_to_top();
    } else {
      m_base_ref = std::make_shared<abstract_domain_t>(m_base_ref->make_top());
    }
  }

  void set_to_bottom() override {
    if (m_base_ref.use_count() == 1) {
      m_base_ref->set_to_bottom();
    } else {
      m_base_ref =
          std::make_shared<abstract_domain_t>(m_base_ref->make_bottom());
    }
  }

  bool is_bottom() const override { return get_absval().is_bottom(); }

  bool is_top() const override { return get_absval().is_top(); }

  bool operator<=(const this_domain_t &o) const override {
    if (same_absval(o)) { // same allocation: trivially included
      return true;
    }
    return get_absval() <= o.get_absval();
  }

  void operator|=(const this_domain_t &o) override {
    if (same_absval(o)) { // x op x == x for all the lattice operations
      return;
    }
    detach_and_get_absval() |= o.get_absval();
  }

  this_domain_t operator|(const this_domain_t &o) const override {
    if (same_absval(o)) {
      return *this;
    }
    return create(get_absval() | o.get_absval());
  }

  void operator&=(const this_domain_t &o) override {
    if (same_absval(o)) {
      return;
    }
    detach_and_get_absval() &= o.get_absval();
  }

  this_domain_t operator&(const this_domain_t &o) const override {
    if (same_absval(o)) {
      return *this;
    }
    return create(get_absval() & o.get_absval());
  }

  this_domain_t operator||(const this_domain_t &o) const override {
    if (same_absval(o)) {
      return *this;
    }
    return create(get_absval() || o.get_absval());
  }

  this_domain_t operator&&(const this_domain_t &o) const override {
    if (same_absval(o)) {
      return *this;
    }
    return create(get_absval() && o.get_absval());
  }

  this_domain_t
  widening_thresholds(const this_domain_t &o,
                      const thresholds<number_t> &ts) const override {
    return create(get_absval().widening_thresholds(o.get_absval(), ts));
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    detach_and_get_absval().apply(op, x, y, z);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    detach_and_get_absval().apply(op, x, y, k);
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    detach_and_get_absval().assign(x, e);
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    detach_and_get_absval().weak_assign(x, e);
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    detach_and_get_absval().operator+=(csts);
  }

  bool entails(const linear_constraint_t &cst) const override {
    return get_absval().entails(cst);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {

    detach_and_get_absval().apply(op, x, y, z);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    detach_and_get_absval().apply(op, x, y, k);
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    detach_and_get_absval().apply(op, dst, src);
  }

  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    detach_and_get_absval().select(lhs, cond, e1, e2);
  }

  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    detach_and_get_absval().assign_bool_cst(lhs, rhs);
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {
    detach_and_get_absval().assign_bool_ref_cst(lhs, rhs);
  }

  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    detach_and_get_absval().assign_bool_var(lhs, rhs, is_not_rhs);
  }

  void weak_assign_bool_cst(const variable_t &lhs,
                            const linear_constraint_t &rhs) override {
    detach_and_get_absval().weak_assign_bool_cst(lhs, rhs);
  }

  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                            bool is_not_rhs) override {
    detach_and_get_absval().weak_assign_bool_var(lhs, rhs, is_not_rhs);
  }

  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {

    detach_and_get_absval().apply_binary_bool(op, x, y, z);
  }

  void assume_bool(const variable_t &v, bool is_negated) override {
    detach_and_get_absval().assume_bool(v, is_negated);
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {
    detach_and_get_absval().select_bool(lhs, cond, b1, b2);
  }

  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {
    detach_and_get_absval().array_init(a, elem_size, lb_idx, ub_idx, val);
  }

  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    detach_and_get_absval().array_load(lhs, a, elem_size, i);
  }

  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &val,
                   bool is_strong_update) override {
    detach_and_get_absval().array_store(a, elem_size, i, val, is_strong_update);
  }

  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &val) override {
    detach_and_get_absval().array_store_range(a, elem_size, i, j, val);
  }

  void array_assign(const variable_t &a, const variable_t &b) override {
    detach_and_get_absval().array_assign(a, b);
  }

  void region_init(const variable_t &reg) override {
    detach_and_get_absval().region_init(reg);
  }

  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {
    detach_and_get_absval().region_copy(lhs_reg, rhs_reg);
  }

  void region_cast(const variable_t &src_reg,
                   const variable_t &dst_reg) override {
    detach_and_get_absval().region_cast(src_reg, dst_reg);
  }

  void ref_make(const variable_t &ref, const variable_t &reg,
                const variable_or_constant_t &size,
                const allocation_site &as) override {
    detach_and_get_absval().ref_make(ref, reg, size, as);
  }

  void ref_free(const variable_t &reg, const variable_t &ref) override {
    detach_and_get_absval().ref_free(reg, ref);
  }

  void ref_load(const variable_t &ref, const variable_t &reg,
                const variable_t &res) override {
    detach_and_get_absval().ref_load(ref, reg, res);
  }

  void ref_store(const variable_t &ref, const variable_t &reg,
                 const variable_or_constant_t &val) override {
    detach_and_get_absval().ref_store(ref, reg, val);
  }

  void ref_gep(const variable_t &ref1, const variable_t &reg1,
               const variable_t &ref2, const variable_t &reg2,
               const linear_expression_t &offset) override {
    detach_and_get_absval().ref_gep(ref1, reg1, ref2, reg2, offset);
  }

  void ref_assume(const reference_constraint_t &cst) override {
    detach_and_get_absval().ref_assume(cst);
  }

  void ref_to_int(const variable_t &reg, const variable_t &ref_var,
                  const variable_t &int_var) override {
    detach_and_get_absval().ref_to_int(reg, ref_var, int_var);
  }

  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref_var) override {
    detach_and_get_absval().int_to_ref(int_var, reg, ref_var);
  }

  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    detach_and_get_absval().select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2,
                                       rgn2);
  }

  boolean_value is_null_ref(const variable_t &ref) override {
    return detach_and_get_absval().is_null_ref(ref);
  }

  bool
  get_allocation_sites(const variable_t &ref,
                       std::vector<allocation_site> &alloc_sites) override {
    return detach_and_get_absval().get_allocation_sites(ref, alloc_sites);
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
                std::vector<uint64_t> &tags) override {
    return detach_and_get_absval().get_tags(rgn, ref, tags);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const this_domain_t &invariant) override {
    detach_and_get_absval().backward_apply(op, x, y, z, invariant.get_absval());
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const this_domain_t &invariant) override {
    detach_and_get_absval().backward_apply(op, x, y, k, invariant.get_absval());
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const this_domain_t &invariant) override {
    detach_and_get_absval().backward_assign(x, e, invariant.get_absval());
  }

  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const this_domain_t &invariant) override {
    detach_and_get_absval().backward_assign_bool_cst(lhs, rhs,
                                                     invariant.get_absval());
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const this_domain_t &invariant) override {
    detach_and_get_absval().backward_assign_bool_ref_cst(
        lhs, rhs, invariant.get_absval());
  }

  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const this_domain_t &invariant) override {
    detach_and_get_absval().backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                                     invariant.get_absval());
  }

  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const this_domain_t &invariant) override {
    detach_and_get_absval().backward_apply_binary_bool(op, x, y, z,
                                                       invariant.get_absval());
  }

  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const this_domain_t &invariant) override {
    detach_and_get_absval().backward_array_init(a, elem_size, lb_idx, ub_idx,
                                                val, invariant.get_absval());
  }

  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const this_domain_t &invariant) override {
    detach_and_get_absval().backward_array_load(lhs, a, elem_size, i,
                                                invariant.get_absval());
  }

  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const this_domain_t &invariant) override {
    detach_and_get_absval().backward_array_store(
        a, elem_size, i, v, is_strong_update, invariant.get_absval());
  }

  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const this_domain_t &invariant) override {
    detach_and_get_absval().backward_array_store_range(a, elem_size, i, j, v,
                                                       invariant.get_absval());
  }

  void backward_array_assign(const variable_t &a, const variable_t &b,
                             const this_domain_t &invariant) override {
    detach_and_get_absval().backward_array_assign(a, b, invariant.get_absval());
  }

  void callee_entry(const callsite_info<variable_t> &callsite,
                    const this_domain_t &caller) override {
    detach_and_get_absval().callee_entry(callsite, caller.get_absval());
  }

  void caller_continuation(const callsite_info<variable_t> &callsite,
                           const this_domain_t &callee) override {
    detach_and_get_absval().caller_continuation(callsite, callee.get_absval());
  }

  void operator-=(const variable_t &v) override {
    detach_and_get_absval().operator-=(v);
  }

  interval_t operator[](const variable_t &v) override {
    return detach_and_get_absval().operator[](v);
  }

  interval_t at(const variable_t &v) const override {
    return get_absval().at(v);
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    return get_absval().to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    return get_absval().to_disjunctive_linear_constraint_system();
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    detach_and_get_absval().rename(from, to);
  }

  void normalize() override { detach_and_get_absval().normalize(); }

  void minimize() override { detach_and_get_absval().minimize(); }

  void forget(const variable_vector_t &variables) override {
    detach_and_get_absval().forget(variables);
  }

  void project(const variable_vector_t &variables) override {
    detach_and_get_absval().project(variables);
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    detach_and_get_absval().expand(var, new_var);
  }

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    detach_and_get_absval().intrinsic(name, inputs, outputs);
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const this_domain_t &invariant) override {
    detach_and_get_absval().backward_intrinsic(name, inputs, outputs,
                                               invariant.get_absval());
  }

  std::string domain_name() const override {
    return get_absval().domain_name();
  }

  void write(crab::crab_os &o) const override { get_absval().write(o); }

  /// @brief true iff both wrappers share the same underlying allocation --
  /// exact value identity at O(1), used for leaf-equality in the odi map
  bool same_absval(const this_domain_t &o) const {
    return m_base_ref == o.m_base_ref;
  }

  // Dereference operator
  const abstract_domain_t &operator*() const { return get_absval(); }
  abstract_domain_t &operator*() { return detach_and_get_absval(); }
};

} // end namespace object_domain_impl

template <typename Variable, typename Domain>
struct abstract_domain_traits<object_domain_impl::abstract_domain_ref<Variable, Domain>> {
  using number_t = typename Variable::number_t;
  using varname_t = typename Variable::varname_t;
};
} // end namespace domains
} // end namespace crab