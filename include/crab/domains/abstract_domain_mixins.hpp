#pragma once

/**
 * Mixins providing empty or default implementations of abstract
 * operations.
 *
 * A domain lists the mixins it wants in its base clause:
 *
 *   class my_domain final
 *     : public abstract_domain_base<my_domain, leaf_numerical_domain> { ... };
 *
 * The mixins are chained linearly, so there is exactly one
 * abstract_domain_api<Dom> subobject and an unambiguous override order.
 *
 * A domain that does not list a mixin leaves the corresponding pure
 * virtuals unimplemented and fails to compile. This is intentional: it
 * keeps visible which operations a domain does not implement and which
 * ones it could implement more efficiently.
 *
 * NOTE: the list is outermost-first. The FIRST mixin sits closest to the
 * domain, the LAST closest to abstract_domain_api. If two mixins define the
 * same operation the first one silently wins, so every mixin declares a
 * "provided_ops" list and abstract_domain_base static_asserts that no
 * operation is claimed twice.
 *
 * WHAT "NOT IMPLEMENTED" MEANS
 *
 * The *_operations_not_implemented mixins do NOT mean that the operation
 * must not be called: none of them raises an error. It means the domain
 * does not care about that family of operations. Transformers are no-ops
 * that leave the state unchanged, and queries answer conservatively:
 *
 *   entails               -> false (cannot prove anything)
 *   at, operator[]        -> top
 *   is_null_ref           -> top
 *   get_allocation_sites  -> false
 *   get_tags              -> false
 *   array_load            -> forgets its left-hand side
 *
 * The queries are safe by construction, but ignoring a transformer is
 * only sound if the domain really tracks nothing that the operation
 * would have had to update. A domain that does track such state and
 * silently ignores the update is unsound, and nothing here diagnoses it.
 * It is up to the caller -- typically the enclosing functor or product
 * domain -- to ensure a domain only receives operations it can safely
 * ignore.
 */

#include <crab/domains/abstract_domain.hpp>
#include <crab/support/stats.hpp>
#include <crab/domains/inter_abstract_operations.hpp>

#include <boost/optional.hpp>
#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace crab {
namespace domains {

/*===================================================================*/
/* Which abstract operation a mixin provides                          */
/*===================================================================*/
enum class domain_op {
  numerical_core,
  bool_core,
  array_core,
  region_core,
  select,
  select_bool,
  select_ref,
  entails,
  weak_assign,
  weak_assign_bool,
  inter_ops
};

template <domain_op... Ops> struct op_list {};

/** Concatenate several op_lists into one **/
template <class... Lists> struct op_concat;

template <> struct op_concat<> { using type = op_list<>; };

template <domain_op... A> struct op_concat<op_list<A...>> {
  using type = op_list<A...>;
};

template <domain_op... A, domain_op... B, class... Rest>
struct op_concat<op_list<A...>, op_list<B...>, Rest...>
    : op_concat<op_list<A..., B...>, Rest...> {};

/** True iff no operation appears twice **/
template <domain_op... Ops> constexpr bool no_duplicate_ops(op_list<Ops...>) {
  constexpr std::size_t n = sizeof...(Ops);
  const domain_op ops[n == 0 ? 1 : n] = {Ops...};
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      if (ops[i] == ops[j]) {
        return false;
      }
    }
  }
  return true;
}

/*===================================================================*/
/* Common base of every mixin: CRTP downcast + the domain typedefs.   */
/* Written once here rather than in each mixin.                       */
/*===================================================================*/
template <class Dom, class Base> class domain_mixin : public Base {
protected:
  const Dom &self() const { return static_cast<const Dom &>(*this); }
  Dom &self() { return static_cast<Dom &>(*this); }

public:
  using api_t = abstract_domain_api<Dom>;
  using number_t = typename api_t::number_t;
  using varname_t = typename api_t::varname_t;
  using variable_t = typename api_t::variable_t;
  using variable_vector_t = typename api_t::variable_vector_t;
  using variable_or_constant_t = typename api_t::variable_or_constant_t;
  using linear_expression_t = typename api_t::linear_expression_t;
  using linear_constraint_t = typename api_t::linear_constraint_t;
  using linear_constraint_system_t = typename api_t::linear_constraint_system_t;
  using reference_constraint_t = typename api_t::reference_constraint_t;
  using interval_t = typename api_t::interval_t;
};

/*===================================================================*/
/* Build the mixin chain                                              */
/*===================================================================*/
template <class Dom, template <class, class> class... Mixins>
struct domain_base_chain;

template <class Dom> struct domain_base_chain<Dom> {
  using type = abstract_domain_api<Dom>;
  using provided_ops = op_list<>;
};

template <class Dom, template <class, class> class M,
          template <class, class> class... Rest>
struct domain_base_chain<Dom, M, Rest...> {
private:
  using tail = domain_base_chain<Dom, Rest...>;

public:
  using type = M<Dom, typename tail::type>;
  using provided_ops =
      typename op_concat<typename M<Dom, typename tail::type>::provided_ops,
                         typename tail::provided_ops>::type;
};

template <class Dom, template <class, class> class... Mixins>
class abstract_domain_base_impl {
  using chain = domain_base_chain<Dom, Mixins...>;
  static_assert(no_duplicate_ops(typename chain::provided_ops{}),
                "two mixins provide the same abstract operation: the first "
                "one in the list would silently win");

public:
  using type = typename chain::type;
};

template <class Dom, template <class, class> class... Mixins>
using abstract_domain_base =
    typename abstract_domain_base_impl<Dom, Mixins...>::type;

/*===================================================================*/
/* Numerical operations not implemented.                              */
/*                                                                    */
/* Note that this also provides select, entails and weak_assign:      */
/* a domain without numerical operations has nothing for the          */
/* corresponding default_* mixins to build on, so listing both is a   */
/* mistake and provided_ops makes it a compile error.                 */
/*===================================================================*/
template <class Dom, class Base>
class numerical_operations_not_implemented : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops =
      op_list<domain_op::numerical_core, domain_op::select,
              domain_op::entails, domain_op::weak_assign>;
  using typename base_t::interval_t;
  using typename base_t::linear_constraint_system_t;
  using typename base_t::linear_constraint_t;
  using typename base_t::linear_expression_t;
  using typename base_t::number_t;
  using typename base_t::variable_t;

  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {}
  void apply(crab::domains::arith_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {}
  void assign(const variable_t &x, const linear_expression_t &e) override {}
  void weak_assign(const variable_t &x,
                   const linear_expression_t &e) override {}
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {}
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const Dom &invariant) override {}
  void backward_apply(crab::domains::arith_operation_t op,
                      const variable_t &x, const variable_t &y, number_t z,
                      const Dom &invariant) override {}
  void backward_apply(crab::domains::arith_operation_t op,
                      const variable_t &x, const variable_t &y,
                      const variable_t &z, const Dom &invariant) override {}
  void operator+=(const linear_constraint_system_t &csts) override {}
  bool entails(const linear_constraint_t &rhs) const override { return false; }
  interval_t operator[](const variable_t &x) override { return this->at(x); }
  interval_t at(const variable_t &x) const override {
    return interval_t::top();
  }
  void apply(crab::domains::int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {}
  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {}
  void apply(crab::domains::bitwise_operation_t op, const variable_t &x,
             const variable_t &y, number_t z) override {}
};

/*===================================================================*/
/* Boolean operations not implemented                                 */
/*===================================================================*/
template <class Dom, class Base>
class bool_operations_not_implemented : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::bool_core>;
  using typename base_t::linear_constraint_t;
  using typename base_t::reference_constraint_t;
  using typename base_t::variable_t;

  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {}
  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {}
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {}
  void weak_assign_bool_cst(const variable_t &lhs,
                            const linear_constraint_t &rhs) override {}
  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                            bool is_not_rhs) override {}
  void apply_binary_bool(crab::domains::bool_operation_t op,
                         const variable_t &x, const variable_t &y,
                         const variable_t &z) override {}
  void assume_bool(const variable_t &v, bool is_negated) override {}
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {}
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const Dom &invariant) override {}
  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const Dom &invariant) override {}
  void backward_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                                bool is_not_rhs,
                                const Dom &invariant) override {}
  void backward_apply_binary_bool(crab::domains::bool_operation_t op,
                                  const variable_t &x, const variable_t &y,
                                  const variable_t &z,
                                  const Dom &invariant) override {}
};

/*===================================================================*/
/* Array operations not implemented                                   */
/*===================================================================*/
template <class Dom, class Base>
class array_operations_not_implemented : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::array_core>;
  using typename base_t::linear_expression_t;
  using typename base_t::variable_t;

  void array_init(const variable_t &a, const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {}
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    this->operator-=(lhs);
  }
  void array_store(const variable_t &a, const linear_expression_t &elem_size,
                   const linear_expression_t &i, const linear_expression_t &v,
                   bool is_strong_update) override {}
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &v) override {}
  void array_assign(const variable_t &lhs, const variable_t &rhs) override {}
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const Dom &invariant) override {}
  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const Dom &invariant) override {}
  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v, bool is_strong_update,
                            const Dom &invariant) override {}
  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const Dom &invariant) override {}
  void backward_array_assign(const variable_t &lhs, const variable_t &rhs,
                             const Dom &invariant) override {}
};

/*===================================================================*/
/* Region and reference operations not implemented                    */
/*===================================================================*/
template <class Dom, class Base>
class region_and_reference_operations_not_implemented
    : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::region_core>;
  using typename base_t::linear_expression_t;
  using typename base_t::reference_constraint_t;
  using typename base_t::variable_or_constant_t;
  using typename base_t::variable_t;

  void region_init(const variable_t &reg) override {}
  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {}
  void region_cast(const variable_t &src_reg,
                   const variable_t &dst_reg) override {}
  void ref_make(const variable_t &ref, const variable_t &reg,
                const variable_or_constant_t &size,
                const crab::allocation_site &as) override {}
  void ref_free(const variable_t &reg, const variable_t &ref) override {}
  void ref_load(const variable_t &ref, const variable_t &reg,
                const variable_t &res) override {}
  void ref_store(const variable_t &ref, const variable_t &reg,
                 const variable_or_constant_t &val) override {}
  void ref_gep(const variable_t &ref1, const variable_t &reg1,
               const variable_t &ref2, const variable_t &reg2,
               const linear_expression_t &offset) override {}
  void ref_assume(const reference_constraint_t &cst) override {}
  void ref_to_int(const variable_t &reg, const variable_t &ref,
                  const variable_t &int_var) override {}
  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref) override {}
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {}
  crab::domains::boolean_value is_null_ref(const variable_t &ref) override {
    return crab::domains::boolean_value::top();
  }
  bool get_allocation_sites(const variable_t &ref,
                            std::vector<crab::allocation_site> &out) override {
    return false;
  }
  bool get_tags(const variable_t &rng, const variable_t &ref,
                std::vector<uint64_t> &out) override {
    return false;
  }
};

/*===================================================================*/
/* Default select: case split on the condition and join.              */
/* A domain that can build the restricted states more directly        */
/* should implement select natively.                                  */
/*===================================================================*/
template <class Dom, class Base>
class default_select : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::select>;
  using typename base_t::linear_constraint_t;
  using typename base_t::linear_expression_t;
  using typename base_t::variable_t;

  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    CRAB_DOMAIN_SCOPED_STATS(this, ".select", 0);
    if (!this->is_bottom()) {
      Dom inv1(this->self());
      inv1 += cond;
      if (inv1.is_bottom()) {
        this->assign(lhs, e2);
        return;
      }
      Dom inv2(this->self());
      inv2 += cond.negate();
      if (inv2.is_bottom()) {
        this->assign(lhs, e1);
        return;
      }
      inv1.assign(lhs, e1);
      inv2.assign(lhs, e2);
      this->self() = inv1 | inv2;
    }
  }
};

/*===================================================================*/
/* Default select_bool: case split on the boolean condition and join. */
/*===================================================================*/
template <class Dom, class Base>
class default_select_bool : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::select_bool>;
  using typename base_t::variable_t;

  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {
    CRAB_DOMAIN_SCOPED_STATS(this, ".select_bool", 0);
    if (!this->is_bottom()) {
      const bool negate = true;
      Dom inv1(this->self());
      inv1.assume_bool(cond, !negate);
      if (inv1.is_bottom()) {
        this->assign_bool_var(lhs, b2, !negate);
        return;
      }
      Dom inv2(this->self());
      inv2.assume_bool(cond, negate);
      if (inv2.is_bottom()) {
        this->assign_bool_var(lhs, b1, !negate);
        return;
      }
      inv1.assign_bool_var(lhs, b1, !negate);
      inv2.assign_bool_var(lhs, b2, !negate);
      this->self() = inv1 | inv2;
    }
  }
};

/*===================================================================*/
/* Default select_ref: case split on the boolean condition, evaluate  */
/* each reference (null, or a gep at offset zero) and join.           */
/*===================================================================*/
template <class Dom, class Base>
class default_select_ref : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::select_ref>;
  using typename base_t::linear_expression_t;
  using typename base_t::number_t;
  using typename base_t::reference_constraint_t;
  using typename base_t::variable_or_constant_t;
  using typename base_t::variable_t;

  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond, const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    CRAB_DOMAIN_SCOPED_STATS(this, ".select_ref", 0);
    if (!this->is_bottom()) {
      auto eval_true_value = [&lhs_ref, &lhs_rgn, &ref1, &rgn1](Dom &out) {
        if (ref1.is_reference_null()) {
          out -= lhs_ref;
          out.ref_assume(reference_constraint_t::mk_null(lhs_ref));
        } else {
          assert(ref1.is_variable());
          assert(rgn1);
          linear_expression_t zero_offset(number_t(0));
          out.ref_gep(ref1.get_variable(), *rgn1, lhs_ref, lhs_rgn,
                      zero_offset);
        }
      };
      auto eval_false_value = [&lhs_ref, &lhs_rgn, &ref2, &rgn2](Dom &out) {
        if (ref2.is_reference_null()) {
          out -= lhs_ref;
          out.ref_assume(reference_constraint_t::mk_null(lhs_ref));
        } else {
          assert(ref2.is_variable());
          assert(rgn2);
          linear_expression_t zero_offset(number_t(0));
          out.ref_gep(ref2.get_variable(), *rgn2, lhs_ref, lhs_rgn,
                      zero_offset);
        }
      };
      const bool negate = true;
      Dom inv1(this->self());
      inv1.assume_bool(cond, !negate);
      if (inv1.is_bottom()) {
        eval_false_value(this->self());
        return;
      }
      Dom inv2(this->self());
      inv2.assume_bool(cond, negate);
      if (inv2.is_bottom()) {
        eval_true_value(this->self());
        return;
      }
      eval_true_value(inv1);
      eval_false_value(inv2);
      this->self() = inv1 | inv2;
    }
  }
};

/*===================================================================*/
/* Default entails: check that the negation of the constraint is      */
/* unsatisfiable. An equality is split into two inequalities.         */
/*===================================================================*/
template <class Dom, class Base>
class default_entails : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::entails>;
  using typename base_t::linear_constraint_system_t;
  using typename base_t::linear_constraint_t;
  using typename base_t::number_t;

  bool entails(const linear_constraint_t &cst) const override {
    if (this->is_bottom()) {
      return true;
    }
    if (cst.is_tautology()) {
      return true;
    }
    if (cst.is_contradiction()) {
      return false;
    }
    auto entailmentFn = [this](const linear_constraint_t &c) -> bool {
      Dom dom(this->self());
      linear_constraint_t neg_c = c.negate();
      dom += neg_c;
      return dom.is_bottom();
    };
    if (cst.is_equality()) {
      linear_constraint_system_t inequalities;
      inequalities += linear_constraint_t(cst.expression(),
                                          linear_constraint_t::INEQUALITY);
      inequalities += linear_constraint_t(cst.expression() * number_t(-1),
                                          linear_constraint_t::INEQUALITY);
      return std::all_of(inequalities.begin(), inequalities.end(),
                         entailmentFn);
    } else {
      return entailmentFn(cst);
    }
  }
};

/*===================================================================*/
/* Default weak assignment: assign on a copy and join it back.        */
/*===================================================================*/
template <class Dom, class Base>
class default_weak_assign : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::weak_assign>;
  using typename base_t::linear_expression_t;
  using typename base_t::variable_t;

  void weak_assign(const variable_t &x,
                   const linear_expression_t &e) override {
    if (!this->is_bottom()) {
      Dom other(this->self());
      other.assign(x, e);
      this->self() |= other;
    }
  }
};

/*===================================================================*/
/* Default weak boolean assignment: as default_weak_assign, for the   */
/* two boolean assignments.                                           */
/*===================================================================*/
template <class Dom, class Base>
class default_weak_bool_assign : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

public:
  using provided_ops = op_list<domain_op::weak_assign_bool>;
  using typename base_t::linear_constraint_t;
  using typename base_t::variable_t;

  void weak_assign_bool_cst(const variable_t &lhs,
                            const linear_constraint_t &rhs) override {
    if (!this->is_bottom()) {
      Dom other(this->self());
      other.assign_bool_cst(lhs, rhs);
      this->self() |= other;
    }
  }

  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                            bool is_not_rhs) override {
    if (!this->is_bottom()) {
      Dom other(this->self());
      other.assign_bool_var(lhs, rhs, is_not_rhs);
      this->self() |= other;
    }
  }
};

/*===================================================================*/
/* Default inter-procedural operations: forward to                    */
/* inter_abstract_operations, which is where the real implementation  */
/* lives. Whether that implementation does anything is decided by     */
/* Dom's own parameters, so the domain must expose them:              */
/*                                                                    */
/*   using params_t = Params;                                         */
/*                                                                    */
/* A domain whose call transfer functions are not a plain forward     */
/* (region_domain, the type-erasure wrappers, dummy_abstract_domain)  */
/* simply does not list this mixin and writes them by hand.           */
/*===================================================================*/
template <class Dom, class Base>
class default_inter_operations : public domain_mixin<Dom, Base> {
  using base_t = domain_mixin<Dom, Base>;

  // NOTE: Dom is incomplete while its own base clause is being formed,
  // so params_t can only be named from a member function body, which is
  // instantiated later.
  template <class D>
  using inter_ops_t =
      inter_abstract_operations<D, D::params_t::implement_inter_transformers>;

public:
  using provided_ops = op_list<domain_op::inter_ops>;
  using typename base_t::variable_t;

  void callee_entry(const callsite_info<variable_t> &callsite,
                    const Dom &caller) override {
    CRAB_DOMAIN_SCOPED_STATS(this, ".callee_entry", 0);
    inter_ops_t<Dom>::callee_entry(callsite, caller, this->self());
  }

  void caller_continuation(const callsite_info<variable_t> &callsite,
                           const Dom &callee) override {
    CRAB_DOMAIN_SCOPED_STATS(this, ".caller_cont", 0);
    inter_ops_t<Dom>::caller_continuation(callsite, callee, this->self());
  }
};

/*===================================================================*/
/* Convenience for a domain that implements only the standard         */
/* numerical operations, i.e. a leaf domain in the hierarchy of       */
/* domains: no boolean, array, region or reference operations.        */
/*                                                                    */
/* provided_ops must list the union of what the composed mixins       */
/* provide: each one hides the provided_ops of the next, so it cannot */
/* be inherited. Keep it in sync when changing the composition.       */
/*===================================================================*/
template <class Dom, class Base>
class leaf_numerical_domain
    : public bool_operations_not_implemented<
          Dom, array_operations_not_implemented<
                   Dom, region_and_reference_operations_not_implemented<
                            Dom, Base>>> {
public:
  using provided_ops = op_list<domain_op::bool_core, domain_op::array_core,
                               domain_op::region_core>;
};

} // namespace domains
} // namespace crab
