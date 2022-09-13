#pragma once

#include <crab/domains/abstract_domain_operators.hpp>
#include <crab/domains/boolean.hpp>
#include <crab/fixpoint/thresholds.hpp>
#include <crab/types/linear_constraints.hpp>
#include <crab/types/reference_constraints.hpp>
#include <crab/types/variable.hpp>
#include <crab/types/tag.hpp>

#include <vector>

namespace crab {

using allocation_site = crab::tag;
  
namespace domains {

template <class Dom> struct abstract_domain_traits;
template <class Number, class VariableName> class abstract_domain_results_api;
  
/**
 * All abstract domains must derive from the abstract_domain_api class
 * and expose publicly all its public typedef's.
 *
 * Use of Curiously Recurring Template Pattern (CRTP).
 *
 * This is a sample of how to implement a new abstract domain:
 *
 * template<typename Number, typename VariableName>
 * class my_new_domain final: public
 *     abstract_domain_api<my_new_domain<Number,VariableName>> {
 *     ...
 *     bool is_bottom() const override {...}
 *     bool is_top() const override {...}
 *     ...
 * };
 *
 *
 * template<typename Number, typename VariableName>
 * struct abstract_domain_traits<my_new_domain<Number,VariableName>> {
 *   using number_t = Number;
 *   using varname_t = VariableName;
 * };
 *
 * The API is divided into:
 *
 * (1) Lattice operations
 * (2) (forward and backward) Numerical operations
 * (3) (forward and backward) Boolean operations
 * (4) (only forward) Region and reference operations
 * (5) (forward and backward) Array operations
 * 
 * Where forward (backward) means forward (backward) semantics. The
 * abstract_domain_api API doesn't provide backward versions for (4)
 * but it should. 
 **/
template <class Dom> class abstract_domain_api:
  public abstract_domain_results_api<typename abstract_domain_traits<Dom>::number_t,
				     typename abstract_domain_traits<Dom>::varname_t> {
public:
  using number_t = typename abstract_domain_traits<Dom>::number_t;
  using varname_t = typename abstract_domain_traits<Dom>::varname_t;

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

  abstract_domain_api() = default;
  virtual ~abstract_domain_api() = default;

  /**************************** Lattice operations ****************************/

  // Return a top abstract value
  virtual Dom make_bottom() const = 0;
  // Return a bottom abstract value
  virtual Dom make_top() const = 0;
  // set *this to top
  virtual void set_to_top() = 0;
  // set *this to bottom
  virtual void set_to_bottom() = 0;
  // return true if the abstract state is bottom
  virtual bool is_bottom() const = 0;
  // return true if the abstract state is top
  virtual bool is_top() const = 0;

  // Inclusion operator: return true if *this is equal or more precise than abs
  virtual bool operator<=(const Dom &abs) const = 0;
  // Join operator: return join(*this, abs)
  virtual Dom operator|(const Dom &abs) const = 0;
  // *this = join(*this, abs)
  virtual void operator|=(const Dom &abs) = 0;
  // *this = meet(*this, abs)
  virtual void operator&=(const Dom &abs) = 0;
  // Meet operator: return meet(*this, abs)
  virtual Dom operator&(const Dom &abs) const = 0;
  // Widening operator: return widening(*this, abs)
  virtual Dom operator||(const Dom &abs) const = 0;
  // Narrowing operator: return narrowing(*this, abs)
  virtual Dom operator&&(const Dom &abs) const = 0;
  // Widening with thresholds: return widening(*this, abs) using thresholds ts
  virtual Dom widening_thresholds(
      const Dom &abs,
      const thresholds<number_t> &ts) const = 0;

  /**************************** Numerical operations *************************/
  // x := y op z
  virtual void apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, const variable_t &z) = 0;
  // x := y op k
  virtual void apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, number_t k) = 0;
  // x := e
  virtual void assign(const variable_t &x, const linear_expression_t &e) = 0;
  // join(*this, copy_of_this(x := e))
  virtual void weak_assign(const variable_t &x, const linear_expression_t &e) = 0;
  // add all constraints \in csts
  virtual void operator+=(const linear_constraint_system_t &csts) = 0;
  // return true only if *this entails rhs
  virtual bool entails(const linear_constraint_t &rhs) const = 0;
  // x := y op z
  virtual void apply(bitwise_operation_t op, const variable_t &x,
                     const variable_t &y, const variable_t &z) = 0;
  // x := y op k
  virtual void apply(bitwise_operation_t op, const variable_t &x,
                     const variable_t &y, number_t k) = 0;
  // dst := src
  virtual void apply(int_conv_operation_t op, const variable_t &dst,
                     const variable_t &src) = 0;

  // if(cond) lhs := e1 else lhs := e2
  virtual void select(const variable_t &lhs, const linear_constraint_t &cond,
		      const linear_expression_t &e1,  const linear_expression_t &e2) = 0;
  
  /**************************** Boolean operations ****************************/
  // lhs := rhs
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) = 0;  
  virtual void assign_bool_ref_cst(const variable_t &lhs,
                                   const reference_constraint_t &rhs) = 0; 
  // lhs := not(rhs) if is_not_rhs
  // lhs := rhs      otherwise
  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) = 0;

  // join(*this, copy_of_this(lhs:=rhs))
  virtual void weak_assign_bool_cst(const variable_t &lhs,
				    const linear_constraint_t &rhs) = 0;
  // join(*this, copy_of_this(lhs:=not(rhs))) if is_not_rhs
  // join(*this, copy_of_this(lhs:=rhs))      otherwise
  virtual void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
				    bool is_not_rhs) = 0;
  // x := y op z
  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y, const variable_t &z) = 0;
  // assume(not(v)) if is_negated
  // assume(v)      otherwise
  virtual void assume_bool(const variable_t &v, bool is_negated) = 0;

  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  virtual void select_bool(const variable_t &lhs, const variable_t &cond,
			   const variable_t &b1, const variable_t &b2) = 0;
  
  /**************************** Array operations *****************************/
  // make a fresh array with contents a[j] initialized to val such that
  // j \in [lb_idx,ub_idx] and j % elem_size == val.
  // elem_size is in bytes.
  virtual void array_init(const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &lb_idx,
                          const linear_expression_t &ub_idx,
                          const linear_expression_t &val) = 0;
  // lhs := a[i] where elem_size is in bytes
  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &i) = 0;
  // a[i] := val where elem_size is in bytes
  virtual void array_store(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const linear_expression_t &val,
                           bool is_strong_update) = 0;
  // forall i<=k<j and k % elem_size == 0 :: a[k] := val.
  // elem_size is in bytes
  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &i,
                                 const linear_expression_t &j,
                                 const linear_expression_t &val) = 0;
  // forall i :: a[i] := b[i]
  virtual void array_assign(const variable_t &a, const variable_t &b) = 0;

  /***************** Regions and reference operations *****************/
  // Initialize a region 
  virtual void region_init(const variable_t &reg) = 0;
  // Make a copy of a region
  virtual void region_copy(const variable_t &lhs_reg,
                           const variable_t &rhs_reg) = 0;
  // Cast between regions of different types
  virtual void region_cast(const variable_t &src_reg,
                           const variable_t &dst_reg) = 0;
  // Create a new reference ref associated with as within region reg 
  virtual void ref_make(const variable_t &ref, const variable_t &reg,
			/* size of the allocation in bytes */
			const variable_or_constant_t &size,
			/* identifier for the allocation site */
			const allocation_site &as) = 0;
  // Remove a reference ref within region reg 
  virtual void ref_free(const variable_t &reg, const variable_t &ref) = 0;
  // Read the content of reference ref within reg. The content is
  // stored in res.
  virtual void ref_load(const variable_t &ref, const variable_t &reg,
                        const variable_t &res) = 0;
  // Write the content of val to the address pointed by ref in region
  // reg.
  virtual void ref_store(const variable_t &ref, const variable_t &reg,
                         const variable_or_constant_t &val) = 0;
  // Create a new reference ref2 to region reg2.
  // The reference ref2 is created by adding offset to ref1.
  virtual void ref_gep(const variable_t &ref1, const variable_t &reg1,
                       const variable_t &ref2, const variable_t &reg2,
                       const linear_expression_t &offset) = 0; 
  // Add constraints between references
  virtual void ref_assume(const reference_constraint_t &cst) = 0;
  // Convert a reference to an integer variable
  virtual void ref_to_int(const variable_t &reg, const variable_t &ref,
                          const variable_t &int_var) = 0;
  // Convert an integer variable to a reference
  virtual void int_to_ref(const variable_t &int_var, const variable_t &reg,
                          const variable_t &ref) = 0;
  // if (cond) ref_gep(ref1, rgn1, lhs_ref, lhs_rgn, 0) else
  //           ref_gep(ref2, rgn2, lhs_ref, lhs_rgn, 0)
  // cond is a boolean variable
  virtual void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
			  const variable_t &cond,
			  const variable_or_constant_t &ref1,
			  const boost::optional<variable_t> &rgn1,
			  const variable_or_constant_t &ref2,
			  const boost::optional<variable_t> &rgn2) = 0;
  
  /**************************** Backward numerical operations ***************/
  // x = y op z
  // Substitute x with y op z in the abstract value
  // The result is meet with invariant.
  virtual void backward_apply(arith_operation_t op, const variable_t &x,
                              const variable_t &y, const variable_t &z,
                              const Dom &invariant) = 0;
  // x = y op k
  // Substitute x with y op k in the abstract value
  // The result is meet with invariant.
  virtual void backward_apply(arith_operation_t op, const variable_t &x,
                              const variable_t &y, number_t k,
                              const Dom &invariant) = 0;
  // x = e
  // Substitute x with e in the abstract value
  // The result is meet with invariant.
  virtual void backward_assign(const variable_t &x,
                               const linear_expression_t &e,
                               const Dom &invariant) = 0;

  /**************************** Backward boolean operations ******************/
  virtual void backward_assign_bool_cst(const variable_t &lhs,
                                        const linear_constraint_t &rhs,
                                        const Dom &invariant) = 0;
  virtual void backward_assign_bool_ref_cst(const variable_t &lhs,
                                            const reference_constraint_t &rhs,
                                            const Dom &invariant) = 0;
  virtual void backward_assign_bool_var(const variable_t &lhs,
                                        const variable_t &rhs, bool is_not_rhs,
                                        const Dom &invariant) = 0;
  virtual void backward_apply_binary_bool(bool_operation_t op,
                                          const variable_t &x,
                                          const variable_t &y,
                                          const variable_t &z,
                                          const Dom &invariant) = 0;

  /**************************** Backward array operations ******************/
  virtual void backward_array_init(const variable_t &a,
                                   const linear_expression_t &elem_size,
                                   const linear_expression_t &lb_idx,
                                   const linear_expression_t &ub_idx,
                                   const linear_expression_t &val,
                                   const Dom &invariant) = 0;
  virtual void backward_array_load(const variable_t &lhs, const variable_t &a,
                                   const linear_expression_t &elem_size,
                                   const linear_expression_t &i,
                                   const Dom &invariant) = 0;
  virtual void backward_array_store(const variable_t &a,
                                    const linear_expression_t &elem_size,
                                    const linear_expression_t &i,
                                    const linear_expression_t &v,
                                    bool is_strong_update,
                                    const Dom &invariant) = 0;
  virtual void backward_array_store_range(const variable_t &a,
                                          const linear_expression_t &elem_size,
                                          const linear_expression_t &i,
                                          const linear_expression_t &j,
                                          const linear_expression_t &v,
                                          const Dom &invariant) = 0;
  virtual void backward_array_assign(const variable_t &a, const variable_t &b,
                                     const Dom &invariant) = 0;

  /**************************** Miscellaneous operations ****************/
  // Forget v
  virtual void operator-=(const variable_t &v) = 0;

  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain. Calling this method might trigger
  // normalization if the underlying domain requires so.
  virtual interval_t operator[](const variable_t &v) = 0;

  // Similar to operator[] but it doesn't modify the internal state.
  virtual interval_t at(const variable_t &v) const = 0;
  
  // Convert the abstract state into a conjunction of linear constraints.
  virtual linear_constraint_system_t to_linear_constraint_system() const = 0;

  // Convert the abstract state into a disjunction of conjunction
  // of linear constraints.
  virtual disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const = 0;

  // Rename in the abstract state the variables "from" with those from to.
  //
  // If any variable from "to" exists already in the abstract state
  // then an error will be raised. This might be a bit restrictive and
  // it can be relaxed if needed in the future.
  virtual void rename(const variable_vector_t &from,
                      const variable_vector_t &to) = 0;

  // Normalize the abstract domain if such notion exists.
  virtual void normalize() = 0;

  // Reduce the size of the abstract domain representation.
  virtual void minimize() = 0;

  // Forget variables form the abstract domain
  virtual void forget(const variable_vector_t &variables) = 0;

  // Project the abstract domain onto variables (dual to forget)
  virtual void project(const variable_vector_t &variables) = 0;

  // Make a new copy of var without relating var with new_var
  virtual void expand(const variable_t &var, const variable_t &new_var) = 0;

  // Function whose semantics is defined by the particular abstract
  // domain
  virtual void intrinsic(std::string name,
			 const variable_or_constant_vector_t &inputs,
                         const variable_vector_t &outputs) = 0;

  virtual void backward_intrinsic(std::string name,
                                  const variable_or_constant_vector_t &inputs,
                                  const variable_vector_t &outputs,
                                  const Dom &invariant) = 0;

  // Print the internal state of the abstract domain
  virtual void write(crab::crab_os &o) const = 0;

  // Return a string the abstract domain name
  virtual std::string domain_name(void) const = 0;

  friend crab::crab_os &operator<<(crab::crab_os &o,
                                   const abstract_domain_api<Dom> &dom) {
    dom.write(o);
    return o;
  }
};

/* 
 * Extend abstract_domain_api to answer specialized queries which are
 * not possible by abstract_domain_api methods such as operator[],
 * to_linear_constraint_system and
 * to_disjunctive_linear_constraint_system.
 */
template<class Number, class VariableName>  
class abstract_domain_results_api {
public:
  using number_t = Number;
  using variable_t = variable<Number, VariableName>;

  virtual ~abstract_domain_results_api() {}
  
  // Return a 3-valued boolean about whether the reference variable
  // ref is null.
  virtual boolean_value is_null_ref(const variable_t &ref) = 0;

  // If return true then out contains all the *possible* allocation
  // sites of the reference variable ref. If return false then nothing
  // is known about its allocation sites.
  virtual bool get_allocation_sites(const variable_t &ref,
				    std::vector<allocation_site> &out) = 0;

  // If return true then out contains all the *possible* tags
  // associated with the reference variable ref within region rgn. If
  // return false then nothing is known about its tags.
  virtual bool get_tags(const variable_t &rgn, const variable_t &ref,
			std::vector<uint64_t> &out) = 0;
  
};


#include "./abstract_domain_macros.def"
  
} // end namespace domains
} // end namespace crab

