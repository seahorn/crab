#pragma once

#include <crab/domains/abstract_domain.hpp>

#include <memory>

namespace crab {
namespace domains {

/*
 * Use of the Type-Erasure idiom to encapsulate an arbitrary abstract
 * domain.
 *
 * Important notes:
 *
 * 1. abstract_domain can only be the __root__ of the
 * hierarchy of abstract domains. That is, __don't__ create a reduced
 * product of abstract_domain with something else, or
 * __don't__ pass abstract_domain as parameter to an array
 * domain.
 *
 * 2. when calling binary operations make sure that the two operands
 * contain the same abstract domain. Otherwise, you will get a runtime
 * error for sure. Using the type-erasure idiom does not allow the
 * compiler to check for this kind of errors anymore.
 */

template <typename Variable>
class abstract_domain final
    : abstract_domain_api<abstract_domain<Variable>> {
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
  using variable_vector_t = std::vector<variable_t>;
  using reference_constraint_t = reference_constraint<number_t, varname_t>;
  using interval_t = ikos::interval<number_t>;

private:
  class abstract_domain_concept {
  public:
    abstract_domain_concept() = default;
    virtual ~abstract_domain_concept() = default;
    abstract_domain_concept(const abstract_domain_concept &) = delete;
    abstract_domain_concept(abstract_domain_concept &&) = delete;
    abstract_domain_concept &
    operator=(const abstract_domain_concept &) = delete;
    abstract_domain_concept &operator=(abstract_domain_concept &&) = delete;
    virtual std::unique_ptr<abstract_domain_concept> clone() const = 0;
    virtual std::unique_ptr<abstract_domain_concept> make_top() const = 0;
    virtual std::unique_ptr<abstract_domain_concept> make_bottom() const = 0;
    virtual void set_to_top() = 0;
    virtual void set_to_bottom() = 0;
    virtual bool is_bottom() const = 0;
    virtual bool is_top() const = 0;
    virtual bool operator<=(const abstract_domain_concept &abs) const = 0;
    virtual std::unique_ptr<abstract_domain_concept>
    operator|(const abstract_domain_concept &abs) const = 0;
    virtual void operator|=(const abstract_domain_concept &abs) = 0;
    virtual std::unique_ptr<abstract_domain_concept>
    operator&(const abstract_domain_concept &abs) const = 0;
    virtual std::unique_ptr<abstract_domain_concept>
    operator||(const abstract_domain_concept &abs) const = 0;
    virtual std::unique_ptr<abstract_domain_concept>
    operator&&(const abstract_domain_concept &abs) const = 0;
    virtual std::unique_ptr<abstract_domain_concept> widening_thresholds(
        const abstract_domain_concept &abs,
        const crab::iterators::thresholds<number_t> &ts) const = 0;
    virtual void apply(arith_operation_t op, const variable_t &x,
                       const variable_t &y, const variable_t &z) = 0;
    virtual void apply(arith_operation_t op, const variable_t &x,
                       const variable_t &y, number_t k) = 0;
    virtual void assign(const variable_t &x, const linear_expression_t &e) = 0;
    virtual void operator+=(const linear_constraint_system_t &csts) = 0;
    virtual void apply(bitwise_operation_t op, const variable_t &x,
                       const variable_t &y, const variable_t &z) = 0;
    virtual void apply(bitwise_operation_t op, const variable_t &x,
                       const variable_t &y, number_t k) = 0;
    virtual void apply(int_conv_operation_t op, const variable_t &dst,
                       const variable_t &src) = 0;
    virtual void assign_bool_cst(const variable_t &lhs,
                                 const linear_constraint_t &rhs) = 0;
    virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                 bool is_not_rhs) = 0;
    virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                   const variable_t &y,
                                   const variable_t &z) = 0;
    virtual void assume_bool(const variable_t &v, bool is_negated) = 0;
    virtual void array_init(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &lb_idx,
                            const linear_expression_t &ub_idx,
                            const linear_expression_t &val) = 0;
    virtual void array_load(const variable_t &lhs, const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i) = 0;
    virtual void array_store(const variable_t &a,
                             const linear_expression_t &elem_size,
                             const linear_expression_t &i,
                             const linear_expression_t &val,
                             bool is_strong_update) = 0;
    virtual void array_store_range(const variable_t &a,
                                   const linear_expression_t &elem_size,
                                   const linear_expression_t &i,
                                   const linear_expression_t &j,
                                   const linear_expression_t &val) = 0;
    virtual void array_assign(const variable_t &a, const variable_t &b) = 0;
    virtual void region_init(const memory_region &reg) = 0;
    virtual void ref_make(const variable_t &ref, const memory_region &reg) = 0;
    virtual void ref_load(const variable_t &ref, const memory_region &reg,
                          const variable_t &res) = 0;
    virtual void ref_store(const variable_t &ref, const memory_region &reg,
                           const linear_expression_t &val) = 0;
    virtual void ref_gep(const variable_t &ref1, const memory_region &reg1,
                         const variable_t &ref2, const memory_region &reg2,
                         const linear_expression_t &offset) = 0;
    virtual void ref_load_from_array(const variable_t &lhs,
                                     const variable_t &ref,
                                     const memory_region &region,
                                     const linear_expression_t &index,
                                     const linear_expression_t &elem_size) = 0;
    virtual void ref_store_to_array(const variable_t &ref,
                                    const memory_region &region,
                                    const linear_expression_t &index,
                                    const linear_expression_t &elem_size,
                                    const linear_expression_t &val) = 0;
    virtual void ref_assume(const reference_constraint_t &cst) = 0;
    virtual void backward_apply(arith_operation_t op, const variable_t &x,
                                const variable_t &y, const variable_t &z,
                                const abstract_domain_concept &invariant) = 0;
    virtual void backward_apply(arith_operation_t op, const variable_t &x,
                                const variable_t &y, number_t k,
                                const abstract_domain_concept &invariant) = 0;
    virtual void backward_assign(const variable_t &x,
                                 const linear_expression_t &e,
                                 const abstract_domain_concept &invariant) = 0;
    virtual void
    backward_assign_bool_cst(const variable_t &lhs,
                             const linear_constraint_t &rhs,
                             const abstract_domain_concept &invariant) = 0;
    virtual void
    backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                             bool is_not_rhs,
                             const abstract_domain_concept &invariant) = 0;
    virtual void
    backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                               const variable_t &y, const variable_t &z,
                               const abstract_domain_concept &invariant) = 0;
    virtual void backward_array_init(
        const variable_t &a, const linear_expression_t &elem_size,
        const linear_expression_t &lb_idx, const linear_expression_t &ub_idx,
        const linear_expression_t &val,
        const abstract_domain_concept &invariant) = 0;
    virtual void
    backward_array_load(const variable_t &lhs, const variable_t &a,
                        const linear_expression_t &elem_size,
                        const linear_expression_t &i,
                        const abstract_domain_concept &invariant) = 0;
    virtual void backward_array_store(
        const variable_t &a, const linear_expression_t &elem_size,
        const linear_expression_t &i, const linear_expression_t &v,
        bool is_strong_update, const abstract_domain_concept &invariant) = 0;
    virtual void backward_array_store_range(
        const variable_t &a, const linear_expression_t &elem_size,
        const linear_expression_t &i, const linear_expression_t &j,
        const linear_expression_t &v,
        const abstract_domain_concept &invariant) = 0;
    virtual void
    backward_array_assign(const variable_t &a, const variable_t &b,
                          const abstract_domain_concept &invariant) = 0;
    virtual void operator-=(const variable_t &v) = 0;
    virtual interval_t operator[](const variable_t &v) = 0;
    virtual linear_constraint_system_t to_linear_constraint_system() const = 0;
    virtual disjunctive_linear_constraint_system_t
    to_disjunctive_linear_constraint_system() const = 0;
    virtual void rename(const variable_vector_t &from,
                        const variable_vector_t &to) = 0;
    virtual void normalize() = 0;
    virtual void minimize() = 0;
    virtual void forget(const variable_vector_t &variables) = 0;
    virtual void project(const variable_vector_t &variables) = 0;
    virtual void expand(const variable_t &var, const variable_t &new_var) = 0;
    virtual void intrinsic(std::string name, const variable_vector_t &inputs,
                           const variable_vector_t &outputs) = 0;
    virtual void
    backward_intrinsic(std::string name, const variable_vector_t &inputs,
                       const variable_vector_t &outputs,
                       const abstract_domain_concept &invariant) = 0;
    virtual std::string domain_name() const = 0;
    virtual void write(crab::crab_os &o) const = 0;
  }; // end class abstract_domain_concept

  template <typename Domain>
  class abstract_domain_model final : public abstract_domain_concept {
  private:
    Domain m_inv;

  public:
    explicit abstract_domain_model(Domain inv) : m_inv(std::move(inv)) {}

    std::unique_ptr<abstract_domain_concept> clone() const override {
      // it would be nice to have std::make_unique
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv));
      return std::move(res);
    }

    std::unique_ptr<abstract_domain_concept> make_top() const override {
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv.make_top()));
      return std::move(res);
    }

    std::unique_ptr<abstract_domain_concept> make_bottom() const override {
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv.make_bottom()));
      return std::move(res);
    }

    void set_to_top() override { m_inv.set_to_top(); }
    void set_to_bottom() override { m_inv.set_to_bottom(); }
    bool is_bottom() const override { return m_inv.is_bottom(); }
    bool is_top() const override { return m_inv.is_top(); }

    // unsafe: if the underlying domain in abs is not Domain then it will crash
    bool operator<=(const abstract_domain_concept &abs) const override {
      return m_inv.operator<=(
          static_cast<const abstract_domain_model *>(&abs)->m_inv);
    }

    // unsafe: if the underlying domain in abs is not Domain then it will crash
    void operator|=(const abstract_domain_concept &abs) override {
      m_inv |= static_cast<const abstract_domain_model *>(&abs)->m_inv;
    }

    // unsafe: if the underlying domain in abs is not Domain then it will crash
    std::unique_ptr<abstract_domain_concept>
    operator|(const abstract_domain_concept &abs) const override {
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv.operator|(
              static_cast<const abstract_domain_model *>(&abs)->m_inv)));
      return std::move(res);
    }

    // unsafe: if the underlying domain in abs is not Domain then it will crash
    std::unique_ptr<abstract_domain_concept>
    operator&(const abstract_domain_concept &abs) const override {
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv.operator&(
              static_cast<const abstract_domain_model *>(&abs)->m_inv)));
      return std::move(res);
    }

    // unsafe: if the underlying domain in abs is not Domain then it will crash
    std::unique_ptr<abstract_domain_concept>
    operator||(const abstract_domain_concept &abs) const override {
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv.operator||(
              static_cast<const abstract_domain_model *>(&abs)->m_inv)));
      return std::move(res);
    }

    // unsafe: if the underlying domain in abs is not Domain then it will crash
    std::unique_ptr<abstract_domain_concept>
    operator&&(const abstract_domain_concept &abs) const override {
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv.operator&&(
              static_cast<const abstract_domain_model *>(&abs)->m_inv)));
      return std::move(res);
    }

    // unsafe: if the underlying domain in abs is not Domain then it will crash
    std::unique_ptr<abstract_domain_concept> widening_thresholds(
        const abstract_domain_concept &abs,
        const crab::iterators::thresholds<number_t> &ts) const override {
      std::unique_ptr<abstract_domain_concept> res(
          new abstract_domain_model(m_inv.widening_thresholds(
              static_cast<const abstract_domain_model *>(&abs)->m_inv, ts)));
      return std::move(res);
    }
    void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
               const variable_t &z) override {
      m_inv.apply(op, x, y, z);
    }
    void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
               number_t k) override {
      m_inv.apply(op, x, y, k);
    }
    void assign(const variable_t &x, const linear_expression_t &e) override {
      m_inv.assign(x, e);
    }
    void operator+=(const linear_constraint_system_t &csts) override {
      m_inv.operator+=(csts);
    }
    void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
               const variable_t &z) override {
      m_inv.apply(op, x, y, z);
    }
    void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
               number_t k) override {
      m_inv.apply(op, x, y, k);
    }
    void apply(int_conv_operation_t op, const variable_t &dst,
               const variable_t &src) override {
      m_inv.apply(op, dst, src);
    }
    void assign_bool_cst(const variable_t &lhs,
                         const linear_constraint_t &rhs) override {
      m_inv.assign_bool_cst(lhs, rhs);
    }
    void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                         bool is_not_rhs) override {
      m_inv.assign_bool_var(lhs, rhs, is_not_rhs);
    }
    void apply_binary_bool(bool_operation_t op, const variable_t &x,
                           const variable_t &y, const variable_t &z) override {
      m_inv.apply_binary_bool(op, x, y, z);
    }
    void assume_bool(const variable_t &v, bool is_negated) override {
      m_inv.assume_bool(v, is_negated);
    }
    void array_init(const variable_t &a, const linear_expression_t &elem_size,
                    const linear_expression_t &lb_idx,
                    const linear_expression_t &ub_idx,
                    const linear_expression_t &val) override {
      m_inv.array_init(a, elem_size, lb_idx, ub_idx, val);
    }
    void array_load(const variable_t &lhs, const variable_t &a,
                    const linear_expression_t &elem_size,
                    const linear_expression_t &i) override {
      m_inv.array_load(lhs, a, elem_size, i);
    }
    void array_store(const variable_t &a, const linear_expression_t &elem_size,
                     const linear_expression_t &i,
                     const linear_expression_t &val,
                     bool is_strong_update) override {
      m_inv.array_store(a, elem_size, i, val, is_strong_update);
    }
    void array_store_range(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const linear_expression_t &j,
                           const linear_expression_t &val) override {
      m_inv.array_store_range(a, elem_size, i, j, val);
    }
    void array_assign(const variable_t &a, const variable_t &b) override {
      m_inv.array_assign(a, b);
    }
    void region_init(const memory_region &reg) override {
      m_inv.region_init(reg);
    }
    void ref_make(const variable_t &ref, const memory_region &reg) override {
      m_inv.ref_make(ref, reg);
    }
    void ref_load(const variable_t &ref, const memory_region &reg,
                  const variable_t &res) override {
      m_inv.ref_load(ref, reg, res);
    }
    void ref_store(const variable_t &ref, const memory_region &reg,
                   const linear_expression_t &val) override {
      m_inv.ref_store(ref, reg, val);
    }
    void ref_gep(const variable_t &ref1, const memory_region &reg1,
                 const variable_t &ref2, const memory_region &reg2,
                 const linear_expression_t &offset) override {
      m_inv.ref_gep(ref1, reg1, ref2, reg2, offset);
    }
    void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
                             const memory_region &region,
                             const linear_expression_t &index,
                             const linear_expression_t &elem_size) override {
      m_inv.ref_load_from_array(lhs, ref, region, index, elem_size);
    }
    void ref_store_to_array(const variable_t &ref, const memory_region &region,
                            const linear_expression_t &index,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &val) override {
      m_inv.ref_store_to_array(ref, region, index, elem_size, val);
    }
    void ref_assume(const reference_constraint_t &cst) override {
      m_inv.ref_assume(cst);
    }

    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_apply(arith_operation_t op, const variable_t &x,
                        const variable_t &y, const variable_t &z,
                        const abstract_domain_concept &invariant) override {
      m_inv.backward_apply(
          op, x, y, z,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_apply(arith_operation_t op, const variable_t &x,
                        const variable_t &y, number_t k,
                        const abstract_domain_concept &invariant) override {
      m_inv.backward_apply(
          op, x, y, k,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_assign(const variable_t &x, const linear_expression_t &e,
                         const abstract_domain_concept &invariant) override {
      m_inv.backward_assign(
          x, e, static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_assign_bool_cst(
        const variable_t &lhs, const linear_constraint_t &rhs,
        const abstract_domain_concept &invariant) override {
      m_inv.backward_assign_bool_cst(
          lhs, rhs,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_assign_bool_var(
        const variable_t &lhs, const variable_t &rhs, bool is_not_rhs,
        const abstract_domain_concept &invariant) override {
      m_inv.backward_assign_bool_var(
          lhs, rhs, is_not_rhs,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_apply_binary_bool(
        bool_operation_t op, const variable_t &x, const variable_t &y,
        const variable_t &z,
        const abstract_domain_concept &invariant) override {
      m_inv.backward_apply_binary_bool(
          op, x, y, z,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_array_init(
        const variable_t &a, const linear_expression_t &elem_size,
        const linear_expression_t &lb_idx, const linear_expression_t &ub_idx,
        const linear_expression_t &val,
        const abstract_domain_concept &invariant) override {
      m_inv.backward_array_init(
          a, elem_size, lb_idx, ub_idx, val,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void
    backward_array_load(const variable_t &lhs, const variable_t &a,
                        const linear_expression_t &elem_size,
                        const linear_expression_t &i,
                        const abstract_domain_concept &invariant) override {
      m_inv.backward_array_load(
          lhs, a, elem_size, i,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void
    backward_array_store(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &v, bool is_strong_update,
                         const abstract_domain_concept &invariant) override {
      m_inv.backward_array_store(
          a, elem_size, i, v, is_strong_update,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_array_store_range(
        const variable_t &a, const linear_expression_t &elem_size,
        const linear_expression_t &i, const linear_expression_t &j,
        const linear_expression_t &v,
        const abstract_domain_concept &invariant) override {
      m_inv.backward_array_store_range(
          a, elem_size, i, j, v,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void
    backward_array_assign(const variable_t &a, const variable_t &b,
                          const abstract_domain_concept &invariant) override {
      m_inv.backward_array_assign(
          a, b, static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    void operator-=(const variable_t &v) override { m_inv.operator-=(v); }
    interval_t operator[](const variable_t &v) override {
      return m_inv.operator[](v);
    }
    linear_constraint_system_t to_linear_constraint_system() const override {
      return m_inv.to_linear_constraint_system();
    }
    disjunctive_linear_constraint_system_t
    to_disjunctive_linear_constraint_system() const override {
      return m_inv.to_disjunctive_linear_constraint_system();
    }
    void rename(const variable_vector_t &from,
                const variable_vector_t &to) override {
      m_inv.rename(from, to);
    }
    void normalize() override { m_inv.normalize(); }
    void minimize() override { m_inv.minimize(); }
    void forget(const variable_vector_t &variables) override {
      m_inv.forget(variables);
    }
    void project(const variable_vector_t &variables) override {
      m_inv.project(variables);
    }
    void expand(const variable_t &var, const variable_t &new_var) override {
      m_inv.expand(var, new_var);
    }
    void intrinsic(std::string name, const variable_vector_t &inputs,
                   const variable_vector_t &outputs) override {
      m_inv.intrinsic(name, inputs, outputs);
    }
    // unsafe: if the underlying domain in invariant is not Domain then it will
    // crash
    void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                            const variable_vector_t &outputs,
                            const abstract_domain_concept &invariant) override {
      m_inv.backward_intrinsic(
          name, inputs, outputs,
          static_cast<const abstract_domain_model *>(&invariant)->m_inv);
    }
    std::string domain_name() const override { return m_inv.domain_name(); }
    void write(crab::crab_os &o) const override { m_inv.write(o); }
  }; // end class abstract_domain_model

  std::unique_ptr<abstract_domain_concept> m_concept;

  explicit abstract_domain(
      std::unique_ptr<abstract_domain_concept> concept)
      : m_concept(std::move(concept)) {}

public:
  /**===================================================**/
  /**                    External API                   **/
  /**===================================================**/

  template <typename Domain>
  abstract_domain(Domain inv)
      : m_concept(new abstract_domain_model<Domain>(std::move(inv))) {}

  ~abstract_domain() = default;

  abstract_domain(const abstract_domain &o)
      : m_concept(o.m_concept->clone()) {}

  abstract_domain &operator=(const abstract_domain &o) {
    if (this != &o) {
      m_concept = o.m_concept->clone();
    }
    return *this;
  }

  abstract_domain(abstract_domain &&o) = default;

  abstract_domain &operator=(abstract_domain &&o) = default;

  abstract_domain make_top() const override {
    return abstract_domain(std::move(m_concept->make_top()));
  }

  abstract_domain make_bottom() const override {
    return abstract_domain(std::move(m_concept->make_bottom()));
  }

  void set_to_top() override { m_concept->set_to_top(); }
  void set_to_bottom() override { m_concept->set_to_bottom(); }
  bool is_bottom() const override { return m_concept->is_bottom(); }
  bool is_top() const override { return m_concept->is_top(); }
  bool operator<=(const abstract_domain &abs) const override {
    return m_concept->operator<=(*(abs.m_concept));
  }
  void operator|=(const abstract_domain &abs) override {
    m_concept->operator|=(*(abs.m_concept));
  }
  abstract_domain
  operator|(const abstract_domain &abs) const override {
    return abstract_domain(
        std::move(m_concept->operator|(*(abs.m_concept))));
  }
  abstract_domain
  operator&(const abstract_domain &abs) const override {
    return abstract_domain(
        std::move(m_concept->operator&(*(abs.m_concept))));
  }
  abstract_domain
  operator||(const abstract_domain &abs) const override {
    return abstract_domain(
        std::move(m_concept->operator||(*(abs.m_concept))));
  }
  abstract_domain
  operator&&(const abstract_domain &abs) const override {
    return abstract_domain(
        std::move(m_concept->operator&&(*(abs.m_concept))));
  }
  abstract_domain widening_thresholds(
      const abstract_domain &abs,
      const crab::iterators::thresholds<number_t> &ts) const override {
    return abstract_domain(
        std::move(m_concept->widening_thresholds(*(abs.m_concept), ts)));
  }
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_concept->apply(op, x, y, z);
  }
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_concept->apply(op, x, y, k);
  }
  void assign(const variable_t &x, const linear_expression_t &e) override {
    m_concept->assign(x, e);
  }
  void operator+=(const linear_constraint_system_t &csts) override {
    m_concept->operator+=(csts);
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    m_concept->apply(op, x, y, z);
  }
  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    m_concept->apply(op, x, y, k);
  }
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    m_concept->apply(op, dst, src);
  }
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    m_concept->assign_bool_cst(lhs, rhs);
  }
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    m_concept->assign_bool_var(lhs, rhs, is_not_rhs);
  }
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    m_concept->apply_binary_bool(op, x, y, z);
  }
  void assume_bool(const variable_t &v, bool is_negated) override {
    m_concept->assume_bool(v, is_negated);
  }
  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {
    m_concept->array_init(a, elem_size, lb_idx, ub_idx, val);
  }
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    m_concept->array_load(lhs, a, elem_size, i);
  }
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &val,
                   bool is_strong_update) override {
    m_concept->array_store(a, elem_size, i, val, is_strong_update);
  }
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &val) override {
    m_concept->array_store_range(a, elem_size, i, j, val);
  }
  void array_assign(const variable_t &a, const variable_t &b) override {
    m_concept->array_assign(a, b);
  }
  void region_init(const memory_region &reg) override {
    m_concept->region_init(reg);
  }
  void ref_make(const variable_t &ref, const memory_region &reg) override {
    m_concept->ref_make(ref, reg);
  }
  void ref_load(const variable_t &ref, const memory_region &reg,
                const variable_t &res) override {
    m_concept->ref_load(ref, reg, res);
  }
  void ref_store(const variable_t &ref, const memory_region &reg,
                 const linear_expression_t &val) override {
    m_concept->ref_store(ref, reg, val);
  }
  void ref_gep(const variable_t &ref1, const memory_region &reg1,
               const variable_t &ref2, const memory_region &reg2,
               const linear_expression_t &offset) override {
    m_concept->ref_gep(ref1, reg1, ref2, reg2, offset);
  }
  void ref_load_from_array(const variable_t &lhs, const variable_t &ref,
                           const memory_region &region,
                           const linear_expression_t &index,
                           const linear_expression_t &elem_size) override {
    m_concept->ref_load_from_array(lhs, ref, region, index, elem_size);
  }
  void ref_store_to_array(const variable_t &ref, const memory_region &region,
                          const linear_expression_t &index,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &val) override {
    m_concept->ref_store_to_array(ref, region, index, elem_size, val);
  }
  void ref_assume(const reference_constraint_t &cst) override {
    m_concept->ref_assume(cst);
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const abstract_domain &invariant) override {
    m_concept->backward_apply(op, x, y, z, *invariant.m_concept);
  }
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const abstract_domain &invariant) override {
    m_concept->backward_apply(op, x, y, k, *invariant.m_concept);
  }
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const abstract_domain &invariant) override {
    m_concept->backward_assign(x, e, *invariant.m_concept);
  }
  void
  backward_assign_bool_cst(const variable_t &lhs,
                           const linear_constraint_t &rhs,
                           const abstract_domain &invariant) override {
    m_concept->backward_assign_bool_cst(lhs, rhs, *invariant.m_concept);
  }
  void
  backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                           bool is_not_rhs,
                           const abstract_domain &invariant) override {
    m_concept->backward_assign_bool_var(lhs, rhs, is_not_rhs,
                                        *invariant.m_concept);
  }
  void backward_apply_binary_bool(
      bool_operation_t op, const variable_t &x, const variable_t &y,
      const variable_t &z, const abstract_domain &invariant) override {
    m_concept->backward_apply_binary_bool(op, x, y, z, *invariant.m_concept);
  }
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const abstract_domain &invariant) override {
    m_concept->backward_array_init(a, elem_size, lb_idx, ub_idx, val,
                                   *invariant.m_concept);
  }
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const abstract_domain &invariant) override {
    m_concept->backward_array_load(lhs, a, elem_size, i, *invariant.m_concept);
  }
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const abstract_domain &invariant) override {
    m_concept->backward_array_store(a, elem_size, i, v, is_strong_update,
                                    *invariant.m_concept);
  }
  void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &i, const linear_expression_t &j,
      const linear_expression_t &v,
      const abstract_domain &invariant) override {
    m_concept->backward_array_store_range(a, elem_size, i, j, v,
                                          *invariant.m_concept);
  }
  void
  backward_array_assign(const variable_t &a, const variable_t &b,
                        const abstract_domain &invariant) override {
    m_concept->backward_array_assign(a, b, *invariant.m_concept);
  }
  void operator-=(const variable_t &v) override { m_concept->operator-=(v); }
  interval_t operator[](const variable_t &v) override {
    return m_concept->operator[](v);
  }
  linear_constraint_system_t to_linear_constraint_system() const override {
    return m_concept->to_linear_constraint_system();
  }
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    return m_concept->to_disjunctive_linear_constraint_system();
  }
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    m_concept->rename(from, to);
  }
  void normalize() override { m_concept->normalize(); }
  void minimize() override { m_concept->minimize(); }
  void forget(const variable_vector_t &variables) override {
    m_concept->forget(variables);
  }
  void project(const variable_vector_t &variables) override {
    m_concept->project(variables);
  }
  void expand(const variable_t &var, const variable_t &new_var) override {
    m_concept->expand(var, new_var);
  }
  void intrinsic(std::string name, const variable_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    m_concept->intrinsic(name, inputs, outputs);
  }
  void backward_intrinsic(std::string name, const variable_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const abstract_domain &invariant) override {
    m_concept->backward_intrinsic(name, inputs, outputs, *invariant.m_concept);
  }
  std::string domain_name() const override { return m_concept->domain_name(); }
  void write(crab::crab_os &o) const override { m_concept->write(o); }
  friend crab::crab_os &
  operator<<(crab::crab_os &o, const abstract_domain<Variable> &dom) {
    dom.write(o);
    return o;
  }
};


  
template <typename Variable>
struct abstract_domain_traits<abstract_domain<Variable>> {
  using number_t = typename Variable::number_t;
  using varname_t = typename Variable::varname_t;
};

} // end namespace domains
} // end namespace crab
