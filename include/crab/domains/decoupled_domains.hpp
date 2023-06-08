#pragma once

/**************************************************************************
 * Decoupled domains based on the paper "Decoupling the Ascending and
 * Descending Phases in Abstract Interpretation" by Arceri, Mastroeni,
 * and Zaffanella published in APLAS'22.
 **************************************************************************/

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/support/stats.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <boost/optional.hpp>
#include <algorithm>
#include <memory>
#include <vector>

namespace crab {
namespace domains {

#define DECOUPLING_DOMAIN_SCOPED_STATS(NAME) \
  CRAB_DOMAIN_SCOPED_STATS(this, NAME, 0)

// A dummy domain pair for the decoupled abstract domain, to be used
// for debugging purposes: we use the very same domain for
// the ascending and descending phases, so that alpha and gamma
// are the identity function.
template <typename Domain>
class dummy_asc_dsc_pair {
public:
  using asc_domain_t = Domain;
  using dsc_domain_t = Domain;
  static asc_domain_t alpha(const dsc_domain_t& dsc) {
    return dsc;
  }
  static dsc_domain_t gamma(const asc_domain_t& asc) {
    return asc;
  }
};

// A domain pair for the decoupled abstract domain where
// both domains are instances of apron_domain: we use method
// SrcDomain::convert_to<DstDomain>() to implement alpha and gamma.
template <typename AscDomain, typename DscDomain>
class apron_asc_dsc_pair {
public:
  using asc_domain_t = AscDomain;
  using dsc_domain_t = DscDomain;
  static asc_domain_t alpha(const dsc_domain_t& dsc) {
    return dsc.template convert_to<asc_domain_t>();
  }
  static dsc_domain_t gamma(const asc_domain_t& asc) {
    return asc.template convert_to<dsc_domain_t>();
  }
};

// The state of a decoupled abstract domain
template <typename AscDscPair>
class asc_dsc_state {
public:
  using asc_dsc_state_t = asc_dsc_state<AscDscPair>;
  using asc_domain_t = typename AscDscPair::asc_domain_t;
  using dsc_domain_t = typename AscDscPair::dsc_domain_t;

  // Note: meant to have public access
  std::unique_ptr<asc_domain_t> m_asc_ptr;
  std::unique_ptr<dsc_domain_t> m_dsc_ptr;
  bool m_asc_phase;

  // By default, start in descending phase
  asc_dsc_state() : m_asc_phase(false) {}

  // Copy ctor needs ad hoc implementation
  asc_dsc_state(const asc_dsc_state_t& other)
    : m_asc_phase(other.m_asc_phase) {
    if (other.m_asc_ptr != nullptr) {
      m_asc_ptr.reset(new asc_domain_t(other.asc()));
    }
    if (other.m_dsc_ptr != nullptr) {
      m_dsc_ptr.reset(new dsc_domain_t(other.dsc()));
    }
  }

  // Copy assignment needs ad hoc implementation
  asc_dsc_state& operator=(const asc_dsc_state_t& other) {
    auto res = other;
    std::swap(*this, res);
    return *this;
  }

  asc_dsc_state(asc_dsc_state_t&&) = default;
  asc_dsc_state& operator=(asc_dsc_state_t&&) = default;
  ~asc_dsc_state() = default;

  bool is_asc_phase() const {
    return m_asc_phase;
  }

  void set_phase(bool asc) {
    m_asc_phase = asc;
  }

  // Logically const: maybe trigger conversion
  void ensure_asc() const {
    if (m_asc_ptr == nullptr) {
      assert(m_dsc_ptr != nullptr);
      auto& ptr = const_cast<asc_dsc_state_t*>(this)->m_asc_ptr;
      DECOUPLING_DOMAIN_SCOPED_STATS(".alpha");
      ptr.reset(new asc_domain_t);
      *ptr = AscDscPair::alpha(*m_dsc_ptr);
    }
  }

  // Logically const: maybe trigger conversion
  void ensure_dsc() const {
    if (m_dsc_ptr == nullptr) {
      assert(m_asc_ptr != nullptr);
      auto& ptr = const_cast<asc_dsc_state_t*>(this)->m_dsc_ptr;
      DECOUPLING_DOMAIN_SCOPED_STATS(".gamma");
      ptr.reset(new dsc_domain_t);
      *ptr = AscDscPair::gamma(*m_asc_ptr);
    }
  }

  // Read-only access
  const asc_domain_t& asc() const {
    ensure_asc();
    return *m_asc_ptr;
  }
  // Read-write access
  asc_domain_t& asc() {
    ensure_asc();
    return *m_asc_ptr;
  }

  // Read-only access
  const dsc_domain_t& dsc() const {
    ensure_dsc();
    return *m_dsc_ptr;
  }
  // Read-write access
  dsc_domain_t& dsc() {
    ensure_dsc();
    return *m_dsc_ptr;
  }

  // Logically const: used just after modifying m_dsc_ptr
  void reset_asc() const {
    const_cast<asc_dsc_state_t*>(this)->m_asc_ptr.reset(nullptr);
  }

  // Logically const: used just after modifying m_asc_ptr
  void reset_dsc() const {
    const_cast<asc_dsc_state_t*>(this)->m_dsc_ptr.reset(nullptr);
  }

  std::string domain_name() const {
    const char* prefix = "Decoupled";
    std::string name1 = asc_domain_t().domain_name();
    std::string name2 = dsc_domain_t().domain_name();
    std::string name;
    name.reserve(name1.size() + name2.size() + 13);
    name.append(prefix);
    name.append("(");
    name.append(name1);
    name.append(",");
    name.append(name2);
    name.append(")");
    return name;
  }
};

class DecoupledDefaultParams {
public:
  enum { implement_inter_transformers = 0 };
};
  
// The decoupled abstract domain, where in the descending (narrowing)
// phase of the analysis we use an abstract domain which is more precise
// that the domain used in the ascending (widening) phase.
template <typename AscDscPair, typename Params = DecoupledDefaultParams>
class decoupled_domain final
  : public abstract_domain_api<decoupled_domain<AscDscPair, Params>> {
public:
  using asc_domain_t = typename AscDscPair::asc_domain_t;
  using dsc_domain_t = typename AscDscPair::dsc_domain_t;
  using decoupled_domain_t = decoupled_domain<AscDscPair, Params>;

  using abstract_domain_t = abstract_domain_api<decoupled_domain_t>;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;
  using number_t = typename asc_domain_t::number_t;
  using varname_t = typename asc_domain_t::varname_t;

  static_assert(std::is_same<number_t, typename dsc_domain_t::number_t>::value,
                "AscDomain and DscDomain must have same type for number_t");
  static_assert(std::is_same<varname_t, typename dsc_domain_t::varname_t>::value,
                "AscDomain and DscDomain must have same type for varname_t");

private:
  using asc_dsc_state_t = asc_dsc_state<AscDscPair>;
  asc_dsc_state_t m_state;

  static decoupled_domain_t
  make_decoupled(bool asc_phase,
                 asc_domain_t* asc_ptr, dsc_domain_t* dsc_ptr) {
    decoupled_domain_t res;
    res.m_state.m_asc_phase = asc_phase;
    res.m_state.m_asc_ptr.reset(asc_ptr);
    res.m_state.m_dsc_ptr.reset(dsc_ptr);
    return res;
  }

  static decoupled_domain_t
  make_from_asc(asc_domain_t&& asc) {
    return make_decoupled(true, new asc_domain_t(std::move(asc)), nullptr);
  }
  static decoupled_domain_t
  make_from_dsc(dsc_domain_t&& dsc) {
    return make_decoupled(false, nullptr, new dsc_domain_t(std::move(dsc)));
  }

  const asc_domain_t& get_asc() const {
    return m_state.asc();
  }
  asc_domain_t& get_asc() {
    return m_state.asc();
  }

  const dsc_domain_t& get_dsc() const {
    return m_state.dsc();
  }
  dsc_domain_t& get_dsc() {
    return m_state.dsc();
  }

  void reset_asc() const {
    assert(m_state.m_dsc_ptr != nullptr);
    m_state.reset_asc();
  }

  void reset_dsc() const {
    assert(m_state.m_asc_ptr != nullptr);
    m_state.reset_dsc();
  }

public:
  decoupled_domain(bool is_asc = false) {
    m_state.set_phase(is_asc);
    set_to_top();
  }
  decoupled_domain(const decoupled_domain_t&) = default;
  decoupled_domain& operator=(const decoupled_domain_t&) = default;
  decoupled_domain(decoupled_domain_t&&) = default;
  decoupled_domain& operator=(decoupled_domain_t&&) = default;
  ~decoupled_domain() = default;

  bool is_asc_phase() const override {
    return m_state.is_asc_phase();
  }
  void set_phase(bool is_ascending) override {
    m_state.set_phase(is_ascending);
  }

  /**************************** Lattice operations ****************************/

  // Return a bottom abstract value
  decoupled_domain_t make_bottom() const override {
    decoupled_domain_t res(is_asc_phase());
    res.set_to_bottom();
    return res;
  }

  // Return a top abstract value
  decoupled_domain_t make_top() const override {
    decoupled_domain_t res(is_asc_phase());
    assert(res.is_top());
    return res;
  }

  // set *this to top
  void set_to_top() override {
    // Do not call get_asc/dsc, as it may trigger a useless conversion
    if (is_asc_phase()) {
      auto asc_ptr = new asc_domain_t;
      asc_ptr->set_to_top();
      m_state.m_asc_ptr.reset(asc_ptr);
      reset_dsc();
    } else {
      auto dsc_ptr = new dsc_domain_t;
      dsc_ptr->set_to_top();
      m_state.m_dsc_ptr.reset(dsc_ptr);
      reset_asc();
    }
  }

  // set *this to bottom
  void set_to_bottom() override {
    // Do not call get_asc/dsc, as it may trigger a useless conversion
    if (is_asc_phase()) {
      auto asc_ptr = new asc_domain_t;
      asc_ptr->set_to_bottom();
      m_state.m_asc_ptr.reset(asc_ptr);
      reset_dsc();
    } else {
      auto dsc_ptr = new dsc_domain_t;
      dsc_ptr->set_to_bottom();
      m_state.m_dsc_ptr.reset(dsc_ptr);
      reset_asc();
    }
  }

  // return true if the abstract state is bottom
  bool is_bottom() const override {
    return is_asc_phase() ? get_asc().is_bottom() : get_dsc().is_bottom();
  }

  // return true if the abstract state is top
  bool is_top() const override {
    return is_asc_phase() ? get_asc().is_top() : get_dsc().is_top();
  }

  // Inclusion operator: return true if *this is equal or more precise than abs
  bool operator<=(const decoupled_domain_t &abs) const override {
    if (is_asc_phase()) {
      return get_asc() <= abs.get_asc();
    } else {
      return get_dsc() <= abs.get_dsc();
    }
  }

  // Join operator: return join(*this, abs)
  decoupled_domain_t operator|(const decoupled_domain_t &abs) const override {
    if (is_asc_phase()) {
      return make_from_asc(get_asc() | abs.get_asc());
    } else {
      return make_from_dsc(get_dsc() | abs.get_dsc());
    }
  }

  // *this = join(*this, abs)
  void operator|=(const decoupled_domain_t &abs) override {
    if (is_asc_phase()) {
      get_asc() |= abs.get_asc();
      reset_dsc();
    } else {
      get_dsc() |= abs.get_dsc();
      reset_asc();
    }
  }

  // *this = meet(*this, abs)
  void operator&=(const decoupled_domain_t &abs) override {
    if (is_asc_phase()) {
      get_asc() &= abs.get_asc();
      reset_dsc();
    } else {
      get_dsc() &= abs.get_dsc();
      reset_asc();
    }
  }

  // Meet operator: return meet(*this, abs)
  decoupled_domain_t operator&(const decoupled_domain_t &abs) const override {
    if (is_asc_phase()) {
      return make_from_asc(get_asc() & abs.get_asc());
    } else {
      return make_from_dsc(get_dsc() & abs.get_dsc());
    }
  }

  // Widening operator: return widening(*this, abs)
  decoupled_domain_t operator||(const decoupled_domain_t &abs) const override {
    assert(is_asc_phase());
    return make_from_asc(get_asc() || abs.get_asc());
  }

  // Narrowing operator: return narrowing(*this, abs)
  decoupled_domain_t operator&&(const decoupled_domain_t &abs) const override {
    assert(!is_asc_phase());
    return make_from_dsc(get_dsc() && abs.get_dsc());
  }

  // Widening with thresholds: return widening(*this, abs) using thresholds ts
  decoupled_domain_t widening_thresholds(
                          const decoupled_domain_t &abs,
                          const thresholds<number_t> &ts) const override {
    assert(is_asc_phase());
    return make_from_asc(get_asc().widening_thresholds(abs.get_asc(), ts));
  }

  /*********************** Results api operations ***********************/
  boolean_value is_null_ref(const variable_t &ref) override {
    return is_asc_phase()
      ? get_asc().is_null_ref(ref)
      : get_dsc().is_null_ref(ref);
  }

  bool get_allocation_sites(const variable_t &ref,
			    std::vector<allocation_site> &out) override {
    return is_asc_phase()
      ? get_asc().get_allocation_sites(ref, out)
      : get_dsc().get_allocation_sites(ref, out);
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
		std::vector<uint64_t> &out) override {
    return is_asc_phase()
      ? get_asc().get_tags(rgn, ref, out)
      : get_dsc().get_tags(rgn, ref, out);
  }

  /**************************** Numerical operations *************************/
  // x := y op z
  void apply(arith_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    if (is_asc_phase()) {
      get_asc().apply(op, x, y, z);
      reset_dsc();
    } else {
      get_dsc().apply(op, x, y, z);
      reset_asc();
    }
  }
  // x := y op k
  void apply(arith_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    if (is_asc_phase()) {
      get_asc().apply(op, x, y, k);
      reset_dsc();
    } else {
      get_dsc().apply(op, x, y, k);
      reset_asc();
    }
  }
  // x := e
  void assign(const variable_t &x, const linear_expression_t &e) override {
    if (is_asc_phase()) {
      get_asc().assign(x, e);
      reset_dsc();
    } else {
      get_dsc().assign(x, e);
      reset_asc();
    }
  }
  // join(*this, copy_of_this(x := e))
  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    if (is_asc_phase()) {
      get_asc().weak_assign(x, e);
      reset_dsc();
    } else {
      get_dsc().weak_assign(x, e);
      reset_asc();
    }
  }
  // add all constraints \in csts
  void operator+=(const linear_constraint_system_t &csts) override {
    if (is_asc_phase()) {
      get_asc() += csts;
      reset_dsc();
    } else {
      get_dsc() += csts;
      reset_asc();
    }
  }
  // return true only if *this entails rhs
  bool entails(const linear_constraint_t &rhs) const override {
    if (is_asc_phase()) {
      return get_asc().entails(rhs);
    } else {
      return get_dsc().entails(rhs);
    }
  }
  // x := y op z
  void apply(bitwise_operation_t op, const variable_t &x,
             const variable_t &y, const variable_t &z) override {
    if (is_asc_phase()) {
      get_asc().apply(op, x, y, z);
      reset_dsc();
    } else {
      get_dsc().apply(op, x, y, z);
      reset_asc();
    }
  }
  // x := y op k
  void apply(bitwise_operation_t op, const variable_t &x,
             const variable_t &y, number_t k) override {
    if (is_asc_phase()) {
      get_asc().apply(op, x, y, k);
      reset_dsc();
    } else {
      get_dsc().apply(op, x, y, k);
      reset_asc();
    }
  }
  // dst := src
  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (is_asc_phase()) {
      get_asc().apply(op, dst, src);
      reset_dsc();
    } else {
      get_dsc().apply(op, dst, src);
      reset_asc();
    }
  }

  // if(cond) lhs := e1 else lhs := e2
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    if (is_asc_phase()) {
      get_asc().select(lhs, cond, e1, e2);
      reset_dsc();
    } else {
      get_dsc().select(lhs, cond, e1, e2);
      reset_asc();
    }
  }

  /**************************** Boolean operations ****************************/
  // lhs := rhs
  void assign_bool_cst(const variable_t &lhs,
                       const linear_constraint_t &rhs) override {
    if (is_asc_phase()) {
      get_asc().assign_bool_cst(lhs, rhs);
      reset_dsc();
    } else {
      get_dsc().assign_bool_cst(lhs, rhs);
      reset_asc();
    }
  }

  void assign_bool_ref_cst(const variable_t &lhs,
                           const reference_constraint_t &rhs) override {
    if (is_asc_phase()) {
      get_asc().assign_bool_ref_cst(lhs, rhs);
      reset_dsc();
    } else {
      get_dsc().assign_bool_ref_cst(lhs, rhs);
      reset_asc();
    }
  }

  // lhs := not(rhs) if is_not_rhs
  // lhs := rhs      otherwise
  void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                       bool is_not_rhs) override {
    if (is_asc_phase()) {
      get_asc().assign_bool_var(lhs, rhs, is_not_rhs);
      reset_dsc();
    } else {
      get_dsc().assign_bool_var(lhs, rhs, is_not_rhs);
      reset_asc();
    }
  }


  // join(*this, copy_of_this(lhs:=rhs))
  void weak_assign_bool_cst(const variable_t &lhs,
                            const linear_constraint_t &rhs) override {
    if (is_asc_phase()) {
      get_asc().weak_assign_bool_cst(lhs, rhs);
      reset_dsc();
    } else {
      get_dsc().weak_assign_bool_cst(lhs, rhs);
      reset_asc();
    }
  }

  // join(*this, copy_of_this(lhs:=not(rhs))) if is_not_rhs
  // join(*this, copy_of_this(lhs:=rhs))      otherwise
  void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                            bool is_not_rhs) override {
    if (is_asc_phase()) {
      get_asc().weak_assign_bool_var(lhs, rhs, is_not_rhs);
      reset_dsc();
    } else {
      get_dsc().weak_assign_bool_var(lhs, rhs, is_not_rhs);
      reset_asc();
    }
  }

  // x := y op z
  void apply_binary_bool(bool_operation_t op, const variable_t &x,
                         const variable_t &y, const variable_t &z) override {
    if (is_asc_phase()) {
      get_asc().apply_binary_bool(op, x, y, z);
      reset_dsc();
    } else {
      get_dsc().apply_binary_bool(op, x, y, z);
      reset_asc();
    }
  }

  // assume(not(v)) if is_negated
  // assume(v)      otherwise
  void assume_bool(const variable_t &v, bool is_negated) override {
    if (is_asc_phase()) {
      get_asc().assume_bool(v, is_negated);
      reset_dsc();
    } else {
      get_dsc().assume_bool(v, is_negated);
      reset_asc();
    }
  }


  // if(cond) lhs := b1 else lhs := b2
  // lhs, cond, b1, and b2 are boolean variables
  void select_bool(const variable_t &lhs, const variable_t &cond,
                   const variable_t &b1, const variable_t &b2) override {
    if (is_asc_phase()) {
      get_asc().select_bool(lhs, cond, b1, b2);
      reset_dsc();
    } else {
      get_dsc().select_bool(lhs, cond, b1, b2);
      reset_asc();
    }
  }

  /**************************** Array operations *****************************/
  // make a fresh array with contents a[j] initialized to val such that
  // j \in [lb_idx,ub_idx] and j % elem_size == val.
  // elem_size is in bytes.
  void array_init(const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &lb_idx,
                  const linear_expression_t &ub_idx,
                  const linear_expression_t &val) override {
    if (is_asc_phase()) {
      get_asc().array_init(a, elem_size, lb_idx, ub_idx, val);
      reset_dsc();
    } else {
      get_dsc().array_init(a, elem_size, lb_idx, ub_idx, val);
      reset_asc();
    }
  }

  // lhs := a[i] where elem_size is in bytes
  void array_load(const variable_t &lhs, const variable_t &a,
                  const linear_expression_t &elem_size,
                  const linear_expression_t &i) override {
    if (is_asc_phase()) {
      get_asc().array_load(lhs, a, elem_size, i);
      reset_dsc();
    } else {
      get_dsc().array_load(lhs, a, elem_size, i);
      reset_asc();
    }
  }

  // a[i] := val where elem_size is in bytes
  void array_store(const variable_t &a,
                   const linear_expression_t &elem_size,
                   const linear_expression_t &i,
                   const linear_expression_t &val,
                   bool is_strong_update) override {
    if (is_asc_phase()) {
      get_asc().array_store(a, elem_size, i, val, is_strong_update);
      reset_dsc();
    } else {
      get_dsc().array_store(a, elem_size, i, val, is_strong_update);
      reset_asc();
    }
  }

  // forall i<=k<j and k % elem_size == 0 :: a[k] := val.
  // elem_size is in bytes
  void array_store_range(const variable_t &a,
                         const linear_expression_t &elem_size,
                         const linear_expression_t &i,
                         const linear_expression_t &j,
                         const linear_expression_t &val) override {
    if (is_asc_phase()) {
      get_asc().array_store_range(a, elem_size, i, j, val);
      reset_dsc();
    } else {
      get_dsc().array_store_range(a, elem_size, i, j, val);
      reset_asc();
    }
  }

  // forall i :: a[i] := b[i]
  void array_assign(const variable_t &a, const variable_t &b) override {
    if (is_asc_phase()) {
      get_asc().array_assign(a, b);
      reset_dsc();
    } else {
      get_dsc().array_assign(a, b);
      reset_asc();
    }
  }

  /***************** Regions and reference operations *****************/
  // Initialize a region
  void region_init(const variable_t &reg) override {
    if (is_asc_phase()) {
      get_asc().region_init(reg);
      reset_dsc();
    } else {
      get_dsc().region_init(reg);
      reset_asc();
    }
  }

  // Make a copy of a region
  void region_copy(const variable_t &lhs_reg,
                   const variable_t &rhs_reg) override {
    if (is_asc_phase()) {
      get_asc().region_copy(lhs_reg, rhs_reg);
      reset_dsc();
    } else {
      get_dsc().region_copy(lhs_reg, rhs_reg);
      reset_asc();
    }
  }

  // Cast between regions of different types
  void region_cast(const variable_t &src_reg,
                   const variable_t &dst_reg) override {
    if (is_asc_phase()) {
      get_asc().region_cast(src_reg, dst_reg);
      reset_dsc();
    } else {
      get_dsc().region_cast(src_reg, dst_reg);
      reset_asc();
    }
  }

  // Create a new reference ref associated with as within region reg
  void ref_make(const variable_t &ref, const variable_t &reg,
                /* size of the allocation in bytes */
                const variable_or_constant_t &size,
                /* identifier for the allocation site */
                const allocation_site &as) override {
    if (is_asc_phase()) {
      get_asc().ref_make(ref, reg, size, as);
      reset_dsc();
    } else {
      get_dsc().ref_make(ref, reg, size, as);
      reset_asc();
    }
  }

  // Remove a reference ref within region reg
  void ref_free(const variable_t &reg, const variable_t &ref) override {
    if (is_asc_phase()) {
      get_asc().ref_free(reg, ref);
      reset_dsc();
    } else {
      get_dsc().ref_free(reg, ref);
      reset_asc();
    }
  }

  // Read the content of reference ref within reg. The content is
  // stored in res.
  void ref_load(const variable_t &ref, const variable_t &reg,
                const variable_t &res) override {
    if (is_asc_phase()) {
      get_asc().ref_load(ref, reg, res);
      reset_dsc();
    } else {
      get_dsc().ref_load(ref, reg, res);
      reset_asc();
    }
  }

  // Write the content of val to the address pointed by ref in region
  // reg.
  void ref_store(const variable_t &ref, const variable_t &reg,
                 const variable_or_constant_t &val) override {
    if (is_asc_phase()) {
      get_asc().ref_store(ref, reg, val);
      reset_dsc();
    } else {
      get_dsc().ref_store(ref, reg, val);
      reset_asc();
    }
  }

  // Create a new reference ref2 to region reg2.
  // The reference ref2 is created by adding offset to ref1.
  void ref_gep(const variable_t &ref1, const variable_t &reg1,
               const variable_t &ref2, const variable_t &reg2,
               const linear_expression_t &offset) override {
    if (is_asc_phase()) {
      get_asc().ref_gep(ref1, reg1, ref2, reg2, offset);
      reset_dsc();
    } else {
      get_dsc().ref_gep(ref1, reg1, ref2, reg2, offset);
      reset_asc();
    }
  }

  // Add constraints between references
  void ref_assume(const reference_constraint_t &cst) override {
    if (is_asc_phase()) {
      get_asc().ref_assume(cst);
      reset_dsc();
    } else {
      get_dsc().ref_assume(cst);
      reset_asc();
    }
  }

  // Convert a reference to an integer variable
  void ref_to_int(const variable_t &reg, const variable_t &ref,
                  const variable_t &int_var) override {
    if (is_asc_phase()) {
      get_asc().ref_to_int(reg, ref, int_var);
      reset_dsc();
    } else {
      get_dsc().ref_to_int(reg, ref, int_var);
      reset_asc();
    }
  }

  // Convert an integer variable to a reference
  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref) override {
    if (is_asc_phase()) {
      get_asc().int_to_ref(int_var, reg, ref);
      reset_dsc();
    } else {
      get_dsc().int_to_ref(int_var, reg, ref);
      reset_asc();
    }
  }

  // if (cond) ref_gep(ref1, rgn1, lhs_ref, lhs_rgn, 0) else
  //           ref_gep(ref2, rgn2, lhs_ref, lhs_rgn, 0)
  // cond is a boolean variable
  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
                  const variable_t &cond,
                  const variable_or_constant_t &ref1,
                  const boost::optional<variable_t> &rgn1,
                  const variable_or_constant_t &ref2,
                  const boost::optional<variable_t> &rgn2) override {
    if (is_asc_phase()) {
      get_asc().select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2, rgn2);
      reset_dsc();
    } else {
      get_dsc().select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2, rgn2);
      reset_asc();
    }
  }

  /**************************** Backward numerical operations ***************/
  // x = y op z
  // Substitute x with y op z in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_apply(op, x, y, z, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_apply(op, x, y, z, invariant.get_dsc());
      reset_asc();
    }
  }

  // x = y op k
  // Substitute x with y op k in the abstract value
  // The result is meet with invariant.
  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t k,
                      const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_apply(op, x, y, k, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_apply(op, x, y, k, invariant.get_dsc());
      reset_asc();
    }
  }

  // x = e
  // Substitute x with e in the abstract value
  // The result is meet with invariant.
  void backward_assign(const variable_t &x,
                       const linear_expression_t &e,
                       const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_assign(x, e, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_assign(x, e, invariant.get_dsc());
      reset_asc();
    }
  }


  /**************************** Backward boolean operations ******************/
  void backward_assign_bool_cst(const variable_t &lhs,
                                const linear_constraint_t &rhs,
                                const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_assign_bool_cst(lhs, rhs, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_assign_bool_cst(lhs, rhs, invariant.get_dsc());
      reset_asc();
    }
  }

  void backward_assign_bool_ref_cst(const variable_t &lhs,
                                    const reference_constraint_t &rhs,
                                    const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_assign_bool_ref_cst(lhs, rhs, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_assign_bool_ref_cst(lhs, rhs, invariant.get_dsc());
      reset_asc();
    }
  }

  void backward_assign_bool_var(const variable_t &lhs,
                                const variable_t &rhs, bool is_not_rhs,
                                const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_assign_bool_var(lhs, rhs, is_not_rhs, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_assign_bool_var(lhs, rhs, is_not_rhs, invariant.get_dsc());
      reset_asc();
    }
	}

  void backward_apply_binary_bool(bool_operation_t op,
                                  const variable_t &x,
                                  const variable_t &y,
                                  const variable_t &z,
                                  const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_apply_binary_bool(op, x, y, z, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_apply_binary_bool(op, x, y, z, invariant.get_dsc());
      reset_asc();
    }
  }


  /**************************** Backward array operations ******************/
  void backward_array_init(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &lb_idx,
                           const linear_expression_t &ub_idx,
                           const linear_expression_t &val,
                           const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_array_init(a, elem_size, lb_idx, ub_idx, val, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_array_init(a, elem_size, lb_idx, ub_idx, val, invariant.get_dsc());
      reset_asc();
    }
  }

  void backward_array_load(const variable_t &lhs, const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &i,
                           const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_array_load(lhs, a, elem_size, i, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_array_load(lhs, a, elem_size, i, invariant.get_dsc());
      reset_asc();
    }
  }

  void backward_array_store(const variable_t &a,
                            const linear_expression_t &elem_size,
                            const linear_expression_t &i,
                            const linear_expression_t &v,
                            bool is_strong_update,
                            const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_array_store(a, elem_size, i, v, is_strong_update,
                                     invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_array_store(a, elem_size, i, v, is_strong_update,
                                     invariant.get_dsc());
      reset_asc();
    }
  }

  void backward_array_store_range(const variable_t &a,
                                  const linear_expression_t &elem_size,
                                  const linear_expression_t &i,
                                  const linear_expression_t &j,
                                  const linear_expression_t &v,
                                  const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_array_store_range(a, elem_size, i, j, v,
                                           invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_array_store_range(a, elem_size, i, j, v,
                                           invariant.get_dsc());
      reset_asc();
    }
  }

  void backward_array_assign(const variable_t &a, const variable_t &b,
                             const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_array_assign(a, b, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_array_assign(a, b, invariant.get_dsc());
      reset_asc();
    }
  }
  /******************** Inter-procedural operations *******************/

  void callee_entry(const crab::domains::callsite_info<variable_t> &callsite,
		    const decoupled_domain_t &caller) override {
    if (is_asc_phase()) {
      inter_abstract_operations<asc_domain_t, Params::implement_inter_transformers>::
	callee_entry(callsite, caller.get_asc(), get_asc());
      reset_dsc();
    } else {
      inter_abstract_operations<dsc_domain_t, Params::implement_inter_transformers>::
	callee_entry(callsite, caller.get_dsc(), get_dsc());
      reset_asc();
    }
  }
  
  void caller_continuation(const crab::domains::callsite_info<variable_t> &callsite,
			   const decoupled_domain_t &callee) override {

    if (is_asc_phase()) {
      inter_abstract_operations<asc_domain_t, Params::implement_inter_transformers>::
	caller_continuation(callsite, callee.get_asc(), get_asc());
      reset_dsc();
    } else {
      inter_abstract_operations<dsc_domain_t, Params::implement_inter_transformers>::
	caller_continuation(callsite, callee.get_dsc(), get_dsc());
      reset_asc();
    }
  }
  
  /**************************** Miscellaneous operations ****************/
  // Forget v
  void operator-=(const variable_t &v) override {
    if (is_asc_phase()) {
      get_asc() -= v;
      reset_dsc();
    } else {
      get_dsc() -= v;
      reset_asc();
    }
  }


  // Return an interval with the possible values of v if such notion
  // exists in the abstract domain. Calling this method might trigger
  // normalization if the underlying domain requires so.
  interval_t operator[](const variable_t &v) override {
    if (is_asc_phase()) {
      return get_asc()[v];
    } else {
      return get_dsc()[v];
    }
  }


  // Similar to operator[] but it doesn't modify the internal state.
  interval_t at(const variable_t &v) const override{
    if (is_asc_phase()) {
      return get_asc().at(v);
    } else {
      return get_dsc().at(v);
    }
  }


  // Convert the abstract state into a conjunction of linear constraints.
  linear_constraint_system_t to_linear_constraint_system() const override {
    if (is_asc_phase()) {
      return get_asc().to_linear_constraint_system();
    } else {
      return get_dsc().to_linear_constraint_system();
    }
  }


  // Convert the abstract state into a disjunction of conjunction
  // of linear constraints.
  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    if (is_asc_phase()) {
      return get_asc().to_disjunctive_linear_constraint_system();
    } else {
      return get_dsc().to_disjunctive_linear_constraint_system();
    }
  }


  // Rename in the abstract state the variables "from" with those from to.
  //
  // If any variable from "to" exists already in the abstract state
  // then an error will be raised. This might be a bit restrictive and
  // it can be relaxed if needed in the future.
  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (is_asc_phase()) {
      get_asc().rename(from, to);
      reset_dsc();
    } else {
      get_dsc().rename(from, to);
      reset_asc();
    }
  }


  // Normalize the abstract domain if such notion exists.
  void normalize() override {
    if  (is_asc_phase()) {
      get_asc().normalize();
      reset_dsc();
    } else {
      get_dsc().normalize();
      reset_asc();
    }
  }


  // Reduce the size of the abstract domain representation.
  void minimize() override {
    if (is_asc_phase()) {
      get_asc().minimize();
      reset_dsc();
    } else {
      get_dsc().minimize();
      reset_asc();
    }
  }


  // Forget variables form the abstract domain
  void forget(const variable_vector_t &variables) override {
    if (is_asc_phase()) {
      get_asc().forget(variables);
      reset_dsc();
    } else {
      get_dsc().forget(variables);
      reset_asc();
    }
  }


  // Project the abstract domain onto variables (dual to forget)
  void project(const variable_vector_t &variables) override {
    if (is_asc_phase()) {
      get_asc().project(variables);
      reset_dsc();
    } else {
      get_dsc().project(variables);
      reset_asc();
    }
  }


  // Make a new copy of var without relating var with new_var
  void expand(const variable_t &var, const variable_t &new_var) override {
    if (is_asc_phase()) {
      get_asc().expand(var, new_var);
      reset_dsc();
    } else {
      get_dsc().expand(var, new_var);
      reset_asc();
    }
  }


  // Function whose semantics is defined by the particular abstract
  // domain
  void intrinsic(std::string name,
                 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    if (is_asc_phase()) {
      get_asc().intrinsic(name, inputs, outputs);
      reset_dsc();
    } else {
      get_dsc().intrinsic(name, inputs, outputs);
      reset_asc();
    }
  }


  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const decoupled_domain_t &invariant) override {
    if (is_asc_phase()) {
      get_asc().backward_intrinsic(name, inputs, outputs, invariant.get_asc());
      reset_dsc();
    } else {
      get_dsc().backward_intrinsic(name, inputs, outputs, invariant.get_dsc());
      reset_asc();
    }
  }


  // Print the internal state of the abstract domain
  void write(crab::crab_os &o) const override {
    if (is_asc_phase()) {
      CRAB_LOG("decoupled", crab::outs() << "phase = ASC : ";);
      get_asc().write(o);
    } else {
      CRAB_LOG("decoupled", crab::outs() << "phase = DSC : ";);
      get_dsc().write(o);
    }
  }


  // Return a string the abstract domain name
  std::string domain_name(void) const override {
    return m_state.domain_name();
  }


}; // class decoupled_domain

template <typename AscDscPair, typename Params>
struct abstract_domain_traits<decoupled_domain<AscDscPair, Params>> {
  using number_t = typename AscDscPair::dsc_domain_t::number_t;
  using varname_t = typename AscDscPair::dsc_domain_t::varname_t;
};

} // end namespace domains
} // namespace crab
