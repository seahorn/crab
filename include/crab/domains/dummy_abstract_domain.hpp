#pragma once

#include <crab/support/debug.hpp>
#include <string>

namespace crab {
namespace domains {
/* Dummy abstract domain that aborts on all operations */
template <class Derived>
class dummy_abstract_domain : public abstract_domain_api<Derived> {
protected:
  virtual std::string not_implemented_msg() const = 0;

public:
  using this_type = dummy_abstract_domain<Derived>;
  using abstract_domain_api_t = abstract_domain_api<Derived>;
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

public:
  dummy_abstract_domain() {}
  void set_to_top() override { CRAB_ERROR(not_implemented_msg()); }
  void set_to_bottom() override { CRAB_ERROR(not_implemented_msg()); }
  Derived make_bottom() const override { CRAB_ERROR(not_implemented_msg()); }
  Derived make_top() const override { CRAB_ERROR(not_implemented_msg()); }
  bool is_bottom() const override { CRAB_ERROR(not_implemented_msg()); }
  bool is_top() const override { CRAB_ERROR(not_implemented_msg()); }
  bool operator<=(const Derived &other) const override {
    CRAB_ERROR(not_implemented_msg());
  }
  void operator|=(const Derived &other) override {
    CRAB_ERROR(not_implemented_msg());
  }
  Derived operator|(const Derived &other) const override {
    CRAB_ERROR(not_implemented_msg());
  }
  void operator&=(const Derived &other) override {
    CRAB_ERROR(not_implemented_msg());
  }  
  Derived operator&(const Derived &other) const override {
    CRAB_ERROR(not_implemented_msg());
  }
  Derived operator||(const Derived &other) const override {
    CRAB_ERROR(not_implemented_msg());
  }
  Derived widening_thresholds(
      const Derived &e,
      const thresholds<number_t> &ts) const override {
    CRAB_ERROR(not_implemented_msg());
  }
  Derived operator&&(const Derived &other) const override {
    CRAB_ERROR(not_implemented_msg());
  }
  void operator-=(const variable_t &var) override {
    CRAB_ERROR(not_implemented_msg());
  }
  interval_t operator[](const variable_t &v) override {
    CRAB_ERROR(not_implemented_msg());
  }
  interval_t at(const variable_t &v) const override {
    CRAB_ERROR(not_implemented_msg());
  }  
  void operator+=(const linear_constraint_system_t &csts) override {
    CRAB_ERROR(not_implemented_msg());
  }
  bool entails(const linear_constraint_t &cst) const override {
    CRAB_ERROR(not_implemented_msg());
  }
  void assign(const variable_t &x, const linear_expression_t &e) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    CRAB_ERROR(not_implemented_msg());
  }  
  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void weak_assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    CRAB_ERROR(not_implemented_msg());
  }  
  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void assume_bool(const variable_t &v, bool is_negated) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                                  const variable_t &y, const variable_t &z,
                                  const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &v) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void region_init(const variable_t &reg) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void region_cast(const variable_t &src, const variable_t &dst) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void ref_make(const variable_t &ref, const variable_t &reg,
                const variable_or_constant_t &size,
                const allocation_site &as) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void ref_free(const variable_t &reg, const variable_t &ref) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void ref_load(const variable_t &ref, const variable_t &reg,
                const variable_t &res) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void ref_store(const variable_t &ref, const variable_t &reg,
                 const variable_or_constant_t &val) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void ref_gep(const variable_t &ref1, const variable_t &reg1,
               const variable_t &ref2, const variable_t &reg2,
               const linear_expression_t &offset) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void ref_assume(const reference_constraint_t &cst) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void ref_to_int(const variable_t &reg, const variable_t &ref,
                  const variable_t &int_var) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    CRAB_ERROR(not_implemented_msg());
  }
  boolean_value is_null_ref(const variable_t &ref) override {
    CRAB_ERROR("mierda");
    CRAB_ERROR(not_implemented_msg());
  }
  bool
  get_allocation_sites(const variable_t &ref,
                       std::vector<allocation_site> &alloc_sites) override {
    CRAB_ERROR(not_implemented_msg());
  }
  bool get_tags(const variable_t &rng, const variable_t &ref,
                std::vector<uint64_t> &out) override {
    CRAB_ERROR(not_implemented_msg());
  }
  linear_constraint_system_t to_linear_constraint_system() const override {
    CRAB_ERROR(not_implemented_msg());
  }
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_ERROR(not_implemented_msg());
  }
  void forget(const variable_vector_t &variables) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void project(const variable_vector_t &variables) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void expand(const variable_t &var, const variable_t &new_var) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void normalize() override { CRAB_ERROR(not_implemented_msg()); }
  void minimize() override { CRAB_ERROR(not_implemented_msg()); }
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const Derived &invariant) override {
    CRAB_ERROR(not_implemented_msg());
  }
  void write(crab_os &o) const override { CRAB_ERROR(not_implemented_msg()); }
  std::string domain_name() const override {
    CRAB_ERROR(not_implemented_msg());
  }
};

} // end namespace domains
} // end namespace crab
