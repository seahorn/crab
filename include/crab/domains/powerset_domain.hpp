#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/interval.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace domains {

/*
 * The powerset domain consists of all possible subsets made from a
 * base domain.
 *
 * There is no a generic way of implementing the widening operation
 * for the powerset domain. The current widening smashes all disjuncts
 * before calling the widening of the base domain.
 */
template <typename Domain>
class powerset_domain final
    : public abstract_domain_api<powerset_domain<Domain>> {

public:
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
  using powerset_domain_t = powerset_domain<Domain>;
  using abstract_domain_t = abstract_domain_api<powerset_domain_t>;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;
  using typename abstract_domain_t::interval_t;
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::reference_constraint_t;
  using typename abstract_domain_t::variable_or_constant_t;
  using typename abstract_domain_t::variable_or_constant_vector_t;  
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;

private:
  using base_dom_vector = std::vector<Domain>;
  /*
    Bottom is represented by a vector of one element whose value is bottom
    Top is represented by a vector of one element whose value is top
    We don't represent the empty powerset.
   */
  base_dom_vector m_disjuncts;

  inline Domain make_base_bottom() const {
    Domain dom;
    return dom.make_bottom();
  }

  inline Domain make_base_top() const {
    Domain dom;
    return dom.make_top();
  }

  powerset_domain(bool is_bottom) {
    if (is_bottom) {
      m_disjuncts.push_back(make_base_bottom());
    } else {
      m_disjuncts.push_back(make_base_top());
    }
  }

  powerset_domain(Domain &&dom) {
    m_disjuncts.emplace_back(std::move(dom));
    normalize_if_top();
  }

  powerset_domain(base_dom_vector &&powerset)
      : m_disjuncts(std::move(powerset)) {
    normalize_if_top();
    if (m_disjuncts.size() > crab_domain_params_man::get().powerset_max_disjuncts()) {
      smash_disjuncts();
    }
  }

  void normalize_if_top() {
    for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
      if (m_disjuncts[i].is_top()) {
        set_to_top();
        return;
      }
    }
  }

  // Remove redundant disjuncts.
  // Expensive operation (quadratic in the number of disjunctions)
  void simplify(base_dom_vector &disjuncts) {
    std::set<unsigned> redundant_set;
    base_dom_vector res;
    res.reserve(disjuncts.size());
    for (unsigned i = 0, sz = disjuncts.size(); i < sz; ++i) {
      Domain &dom = disjuncts[i];
      bool redundant = false;
      for (unsigned j = 0; j < sz; ++j) {
        if (i != j && (redundant_set.count(j) <= 0) && (dom <= disjuncts[j])) {
          redundant_set.insert(i);
          redundant = true;
          break;
        }
      }
      if (!redundant) {
        res.push_back(dom);
      }
    }
    std::swap(disjuncts, res);
  }

  Domain smash_disjuncts(const powerset_domain_t &pw) const {
    if (pw.is_bottom()) {
      return make_base_bottom();
    } else if (pw.is_top()) {
      return make_base_top();
    }

    assert(!pw.m_disjuncts.empty());
    Domain res = pw.m_disjuncts[0];
    for (unsigned i = 1, sz = pw.m_disjuncts.size(); i < sz; ++i) {
      res |= pw.m_disjuncts[i];
    }
    return res;
  }

  Domain smash_disjuncts() {
    if (is_bottom()) {
      set_to_bottom();
    } else if (is_top()) {
      set_to_top();
    }

    CRAB_LOG("powerset", crab::outs() << "Smashing the powerset\n"
                                      << *this << " into \n";);
    assert(!m_disjuncts.empty());
    Domain res = m_disjuncts[0];
    for (unsigned i = 1, sz = m_disjuncts.size(); i < sz; ++i) {
      res |= m_disjuncts[i];
    }
    m_disjuncts.clear();
    m_disjuncts.push_back(res);
    CRAB_LOG("powerset", crab::outs() << *this << "\n";);
    return res;
  }

  void insert(base_dom_vector &vec, Domain dom) const {
    for (unsigned i = 0, sz = vec.size(); i < sz; ++i) {
      if (dom <= vec[i]) {
        return;
      }
    }
    vec.push_back(dom);
  }

  void append(base_dom_vector &vec1, const base_dom_vector &vec2) const {
    for (unsigned i = 0, sz = vec2.size(); i < sz; ++i) {
      insert(vec1, vec2[i]);
    }
  }

  powerset_domain_t powerset_join_with(const powerset_domain_t &other) const {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    }
    base_dom_vector res(m_disjuncts.begin(), m_disjuncts.end());
    append(res, other.m_disjuncts);
    return powerset_domain_t(std::move(res));
  }

  powerset_domain_t powerset_meet_with(const powerset_domain_t &other) const {
    if (is_bottom() || other.is_bottom()) {
      powerset_domain_t bot(true);
      return bot;
    } else if (is_top()) {
      return other;
    } else if (other.is_top()) {
      return *this;
    }

    base_dom_vector res;
    res.reserve(m_disjuncts.size() + other.m_disjuncts.size());
    for (unsigned i = 0, sz_i = m_disjuncts.size(); i < sz_i; ++i) {
      for (unsigned j = 0, sz_j = other.m_disjuncts.size(); j < sz_j; ++j) {
        Domain meet = m_disjuncts[i] & other.m_disjuncts[j];
        if (!meet.is_bottom()) {
          res.emplace_back(std::move(meet));
        }
      }
    }
    return powerset_domain_t(std::move(res));
  }

public:
  powerset_domain() {
    Domain disjunct; // top by default
    m_disjuncts.push_back(disjunct);
  }
  powerset_domain(const powerset_domain_t &other) = default;
  powerset_domain(powerset_domain_t &&other) = default;
  powerset_domain_t &operator=(const powerset_domain_t &other) = default;
  powerset_domain_t &operator=(powerset_domain_t &&other) = default;

  virtual bool is_bottom() const override {
    for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
      if (!m_disjuncts[i].is_bottom()) {
        return false;
      }
    }
    return true;
  }

  virtual bool is_top() const override {
    for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
      if (m_disjuncts[i].is_top()) {
        // set_to_top();
        return true;
      }
    }
    return false;
  }

  virtual powerset_domain_t make_top() const override {
    base_dom_vector disjuncts;
    Domain disjunct; // top by default
    disjuncts.push_back(disjunct);
    return disjuncts;
  }

  virtual powerset_domain_t make_bottom() const override {
    base_dom_vector disjuncts;
    Domain disjunct;
    disjunct.set_to_bottom();
    disjuncts.push_back(disjunct);
    return disjuncts;
  }

  virtual void set_to_top() override {
    m_disjuncts.clear();
    Domain disjunct; // top by default
    m_disjuncts.push_back(disjunct);
  }

  virtual void set_to_bottom() override {
    m_disjuncts.clear();
    Domain disjunct;
    disjunct.set_to_bottom();
    m_disjuncts.push_back(disjunct);
  }

  virtual bool operator<=(const powerset_domain_t &other) const override {
    powerset_domain_t pow_left(*this);
    powerset_domain_t pow_right(other);
    Domain left = smash_disjuncts(pow_left);
    Domain right = smash_disjuncts(pow_right);
    return left <= right;
  }

  virtual void operator|=(const powerset_domain_t &other) override {
    CRAB_LOG("powerset", crab::outs() << "JOIN \n"
                                      << *this << " and\n"
                                      << other << "=\n";);
    if (is_top() || other.is_bottom()) {
      CRAB_LOG("powerset", crab::outs() << *this << "\n";);
      return;
    } else if (is_bottom()) {
      *this = other;
      CRAB_LOG("powerset", crab::outs() << *this << "\n";);
      return;
    } else if (other.is_top()) {
      set_to_top();
      CRAB_LOG("powerset", crab::outs() << *this << "\n";);
      return;
    }

    append(m_disjuncts, other.m_disjuncts);
    if (m_disjuncts.size() > crab_domain_params_man::get().powerset_max_disjuncts()) {
      smash_disjuncts();
    }
    CRAB_LOG("powerset", crab::outs() << *this << "\n";);
  }

  virtual powerset_domain_t
  operator|(const powerset_domain_t &other) const override {
    return powerset_join_with(other);
  }

  virtual powerset_domain_t
  operator&(const powerset_domain_t &other) const override {
    if (crab_domain_params_man::get().powerset_exact_meet()) {
      return powerset_meet_with(other);
    } else {
      Domain left = smash_disjuncts(*this);
      Domain right = smash_disjuncts(other);
      return powerset_domain_t(left & right);
    }
  }

  virtual void operator&=(const powerset_domain_t &other) override {
    // TODO: improve this by avoiding the copy of the left operand.
    *this = *this & other;
  }
  
  virtual powerset_domain_t
  operator||(const powerset_domain_t &other) const override {
    Domain left = smash_disjuncts(*this);
    Domain right = smash_disjuncts(other);
    return powerset_domain_t(left || right);
  }

  virtual powerset_domain_t widening_thresholds(
      const powerset_domain_t &other,
      const thresholds<number_t> &ts) const override {
    Domain left = smash_disjuncts(*this);
    Domain right = smash_disjuncts(other);
    return powerset_domain_t(left.widening_thresholds(right, ts));
  }

  virtual powerset_domain_t
  operator&&(const powerset_domain_t &other) const override {
    Domain left = smash_disjuncts(*this);
    Domain right = smash_disjuncts(other);
    return powerset_domain_t(left && right);
  }

  virtual void assign(const variable_t &x,
                      const linear_expression_t &e) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].assign(x, e);
      }
    }
  }

  virtual void weak_assign(const variable_t &x,
                      const linear_expression_t &e) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].weak_assign(x, e);
      }
    }
  }
  
  virtual void apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, const variable_t &z) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].apply(op, x, y, z);
      }
    }
  }

  virtual void apply(arith_operation_t op, const variable_t &x,
                     const variable_t &y, number_t k) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].apply(op, x, y, k);
      }
    }
  }

  virtual void select(const variable_t &lhs, const linear_constraint_t &cond,
		      const linear_expression_t &e1,  const linear_expression_t &e2) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].select(lhs, cond, e1, e2);
      }
    }
  }
  
  virtual void backward_assign(const variable_t &x,
                               const linear_expression_t &e,
                               const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_apply(arith_operation_t op, const variable_t &x,
                              const variable_t &y, number_t k,
                              const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_apply(arith_operation_t op, const variable_t &x,
                              const variable_t &y, const variable_t &z,
                              const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void operator+=(const linear_constraint_system_t &csts) override {
    if (is_bottom() || csts.is_true()) {
      return;
    }
    if (csts.is_false()) {
      set_to_bottom();
    }
    CRAB_LOG("powerset", crab::outs() << "Adding " << csts << "\n";);
    base_dom_vector vec;
    vec.reserve(m_disjuncts.size());
    for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
      m_disjuncts[i].operator+=(csts);
      if (m_disjuncts[i].is_bottom()) {
        continue;
      }
      vec.emplace_back(m_disjuncts[i]);
    }
    std::swap(m_disjuncts, vec);
    CRAB_LOG("powerset", crab::outs() << "Res=" << *this << "\n";);
  }

  virtual bool entails(const linear_constraint_t &cst) const override {
    if (is_bottom()) {
      return true;
    }
    
    for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
      if (!m_disjuncts[i].entails(cst)) {
	return false;
      } 
    }
    return true;
  }
  
  virtual void operator-=(const variable_t &v) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i] -= v;
        if (m_disjuncts[i].is_top()) {
          set_to_top();
        }
      }
    }
  }

  // cast_operators_api

  virtual void apply(int_conv_operation_t op, const variable_t &dst,
                     const variable_t &src) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].apply(op, dst, src);
      }
    }
  }

  // bitwise_operators_api

  virtual void apply(bitwise_operation_t op, const variable_t &x,
                     const variable_t &y, const variable_t &z) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].apply(op, x, y, z);
      }
    }
  }

  virtual void apply(bitwise_operation_t op, const variable_t &x,
                     const variable_t &y, number_t k) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].apply(op, x, y, k);
      }
    }
  }

  // array_operators_api

  virtual void array_init(const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &lb_idx,
                          const linear_expression_t &ub_idx,
                          const linear_expression_t &val) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].array_init(a, elem_size, lb_idx, ub_idx, val);
      }
    }
  }

  virtual void array_load(const variable_t &lhs, const variable_t &a,
                          const linear_expression_t &elem_size,
                          const linear_expression_t &idx) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].array_load(lhs, a, elem_size, idx);
      }
    }
  }

  virtual void array_store(const variable_t &a,
                           const linear_expression_t &elem_size,
                           const linear_expression_t &idx,
                           const linear_expression_t &val,
                           bool is_strong_update) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].array_store(a, elem_size, idx, val, is_strong_update);
      }
    }
  }

  virtual void array_store_range(const variable_t &a,
                                 const linear_expression_t &elem_size,
                                 const linear_expression_t &lb_idx,
                                 const linear_expression_t &ub_idx,
                                 const linear_expression_t &val) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].array_store_range(a, elem_size, lb_idx, ub_idx, val);
      }
    }
  }

  virtual void array_assign(const variable_t &lhs,
                            const variable_t &rhs) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].array_assign(lhs, rhs);
      }
    }
  }

  // backward array operations

  virtual void
  backward_array_init(const variable_t &a, const linear_expression_t &elem_size,
                      const linear_expression_t &lb_idx,
                      const linear_expression_t &ub_idx,
                      const linear_expression_t &val,
                      const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void
  backward_array_load(const variable_t &lhs, const variable_t &a,
                      const linear_expression_t &elem_size,
                      const linear_expression_t &idx,
                      const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_array_store(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &idx, const linear_expression_t &v,
      bool is_strong_update, const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_array_store_range(
      const variable_t &a, const linear_expression_t &elem_size,
      const linear_expression_t &lb_idx, const linear_expression_t &ub_idx,
      const linear_expression_t &v,
      const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void
  backward_array_assign(const variable_t &a, const variable_t &b,
                        const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  // region/references
  virtual void region_init(const variable_t &reg) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].region_init(reg);
      }
    }
  }

  virtual void region_copy(const variable_t &lhs_reg,
                           const variable_t &rhs_reg) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].region_copy(lhs_reg, rhs_reg);
      }
    }
  }

  virtual void region_cast(const variable_t &src_reg,
                           const variable_t &dst_reg) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].region_cast(src_reg, dst_reg);
      }
    }
  }
  
  virtual void ref_make(const variable_t &ref, const variable_t &reg,
			const variable_or_constant_t &size,
			const allocation_site &as) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].ref_make(ref, reg, size, as);
      }
    }
  }

  virtual void ref_free(const variable_t &reg, const variable_t &ref) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].ref_free(reg, ref);
      }
    }
  }
  
  virtual void ref_load(const variable_t &ref, const variable_t &reg,
                        const variable_t &res) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].ref_load(ref, reg, res);
      }
    }
  }

  virtual void ref_store(const variable_t &ref, const variable_t &reg,
                         const variable_or_constant_t &val) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].ref_store(ref, reg, val);
      }
    }
  }

  virtual void ref_gep(const variable_t &ref1, const variable_t &reg1,
                       const variable_t &ref2, const variable_t &reg2,
                       const linear_expression_t &offset) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].ref_gep(ref1, reg1, ref2, reg2, offset);
      }
    }
  }

  virtual void ref_assume(const reference_constraint_t &cst) override {
    if (!is_bottom()) {
      base_dom_vector vec;
      vec.reserve(m_disjuncts.size());
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].ref_assume(cst);
        if (m_disjuncts[i].is_bottom()) {
          continue;
        }
        vec.emplace_back(m_disjuncts[i]);
      }
      std::swap(m_disjuncts, vec);
    }
  }

  void ref_to_int(const variable_t &reg, const variable_t &ref_var,
                  const variable_t &int_var) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].ref_to_int(reg, ref_var, int_var);
      }
    }
  }

  void int_to_ref(const variable_t &int_var, const variable_t &reg,
                  const variable_t &ref_var) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].int_to_ref(int_var, reg, ref_var);
      }
    }
  }

  void select_ref(const variable_t &lhs_ref, const variable_t &lhs_rgn,
		  const variable_t &cond,
		  const variable_or_constant_t &ref1,
		  const boost::optional<variable_t> &rgn1,
		  const variable_or_constant_t &ref2,
		  const boost::optional<variable_t> &rgn2) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].select_ref(lhs_ref, lhs_rgn, cond, ref1, rgn1, ref2, rgn2);
      }
    }
  }

  boolean_value is_null_ref(const variable_t &ref) override {
    boolean_value r = boolean_value::bottom();
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
	r = r |  m_disjuncts[i].is_null_ref(ref);
      }
    }
    return r;
  }
  bool get_allocation_sites(const variable_t &ref,
			    std::vector<allocation_site> &alloc_sites) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        if (!m_disjuncts[i].get_allocation_sites(ref, alloc_sites)) {
	  return false;
	}
      }
      return true;
    }
    return false;
  }

  bool get_tags(const variable_t &rgn, const variable_t &ref,
		std::vector<uint64_t> &tags) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        if (!m_disjuncts[i].get_tags(rgn, ref, tags)) {
	  return false;
	}
      }
      return true;
    }
    return false;
  }
  
  // boolean operators
  virtual void assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].assign_bool_cst(lhs, rhs);
      }
    }
  }

  virtual void assign_bool_ref_cst(const variable_t &lhs,
                                   const reference_constraint_t &rhs) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].assign_bool_ref_cst(lhs, rhs);
      }
    }
  }

  virtual void assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].assign_bool_var(lhs, rhs, is_not_rhs);
      }
    }
  }

  virtual void weak_assign_bool_cst(const variable_t &lhs,
                               const linear_constraint_t &rhs) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].weak_assign_bool_cst(lhs, rhs);
      }
    }
  }

  virtual void weak_assign_bool_var(const variable_t &lhs, const variable_t &rhs,
                               bool is_not_rhs) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].weak_assign_bool_var(lhs, rhs, is_not_rhs);
      }
    }
  }
  
  virtual void apply_binary_bool(bool_operation_t op, const variable_t &x,
                                 const variable_t &y,
                                 const variable_t &z) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].apply_binary_bool(op, x, y, z);
      }
    }
  }

  virtual void assume_bool(const variable_t &v, bool is_negated) override {
    if (!is_bottom()) {
      base_dom_vector vec;
      vec.reserve(m_disjuncts.size());
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].assume_bool(v, is_negated);
        if (m_disjuncts[i].is_bottom()) {
          continue;
        }
        vec.emplace_back(m_disjuncts[i]);
      }
      std::swap(m_disjuncts, vec);
    }
  }

  void select_bool(const variable_t &lhs, const variable_t &cond,
		   const variable_t &b1, const variable_t &b2) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].select_bool(lhs, cond, b1, b2);
      }
    }
  }
  
  // backward boolean operators
  virtual void backward_assign_bool_cst(const variable_t &lhs,
                                        const linear_constraint_t &rhs,
                                        const powerset_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void
  backward_assign_bool_ref_cst(const variable_t &lhs,
                               const reference_constraint_t &rhs,
                               const powerset_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void backward_assign_bool_var(const variable_t &lhs,
                                        const variable_t &rhs, bool is_not_rhs,
                                        const powerset_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  virtual void
  backward_apply_binary_bool(bool_operation_t op, const variable_t &x,
                             const variable_t &y, const variable_t &z,
                             const powerset_domain_t &inv) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  // Intrinsics

  virtual void intrinsic(std::string name,
			 const variable_or_constant_vector_t &inputs,
                         const variable_vector_t &outputs) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].intrinsic(name, inputs, outputs);
      }
    }
  }

  virtual void backward_intrinsic(std::string name,
                                  const variable_or_constant_vector_t &inputs,
                                  const variable_vector_t &outputs,
                                  const powerset_domain_t &invariant) override {
    CRAB_WARN(domain_name(), " does not implement backward operations");
  }

  // Miscellaneous operations

  interval_t operator[](const variable_t &v) override {
    return at(v);
  }

  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    }
    if (m_disjuncts.empty()) {
      return interval_t::top();
    }
    
    interval_t res = m_disjuncts[0].at(v);
    for (unsigned i = 1, sz = m_disjuncts.size(); i < sz; ++i) {
      res = res | m_disjuncts[i].at(v);
    }
    return res;
  }
  
  void set(const variable_t &v, interval_t intv) {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].set(v, intv);
      }
    }
  }

  virtual void normalize() override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].normalize();
      }
    }
  }

  virtual void minimize() override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].minimize();
      }
    }
  }

  virtual void rename(const variable_vector_t &from,
                      const variable_vector_t &to) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].rename(from, to);
      }
    }
  }

  virtual void expand(const variable_t &x, const variable_t &new_x) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].expand(x, new_x);
      }
    }
  }

  virtual void forget(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].forget(variables);
        if (m_disjuncts[i].is_top()) {
          set_to_top();
        }
      }
    }
  }

  virtual void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        m_disjuncts[i].project(variables);
      }
    }
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    Domain smashed = smash_disjuncts(*this);
    return smashed.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    if (is_bottom()) {
      disjunctive_linear_constraint_system_t res(true);
      return res;
    } else if (is_top()) {
      disjunctive_linear_constraint_system_t res(false);
      return res;
    } else {
      disjunctive_linear_constraint_system_t res(true);
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
        res += m_disjuncts[i].to_linear_constraint_system();
      }
      return res;
    }
  }

  void write(crab::crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      for (unsigned i = 0, sz = m_disjuncts.size(); i < sz;) {
        o << m_disjuncts[i];
        ++i;
        if (i < sz) {
          o << " or \n";
        }
      }
    }
  }

  void dump(void) {
    crab::outs() << "== Begin powerset internal representation === \n";
    if (m_disjuncts.empty()) {
      crab::outs() << "empty\n";
    }
    for (unsigned i = 0, sz = m_disjuncts.size(); i < sz; ++i) {
      crab::outs() << m_disjuncts[i] << " || ";
    }
    crab::outs() << "== End powerset internal representation === \n";
  }

  std::string domain_name() const override {
    Domain dom;
    return std::string("Powerset(") + dom.domain_name() + ")";
  }
};

template <typename Domain>
struct abstract_domain_traits<powerset_domain<Domain>> {
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
};

} // namespace domains
} // namespace crab
