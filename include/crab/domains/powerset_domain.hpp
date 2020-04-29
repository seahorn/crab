#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/backward_assign_operations.hpp>
#include <crab/domains/linear_constraints.hpp>
#include <crab/domains/interval.hpp>

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
namespace powerset_impl {
class DefaultParams {
public:
  /// Exact meet 
  enum { exact_meet = 0};
  /// Smash if the number of disjunctions exceeds this threshold.
  enum { max_disjuncts = 99999};
};
}
  
template<typename Domain, class Params = powerset_impl::DefaultParams>
class powerset_domain final
  : public abstract_domain<powerset_domain<Domain, Params>> {
  
public:
  
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
  using powerset_domain_t = powerset_domain<Domain, Params>;
  using abstract_domain_t = abstract_domain<powerset_domain_t>;
  using typename abstract_domain_t::disjunctive_linear_constraint_system_t;  
  using typename abstract_domain_t::linear_constraint_system_t;
  using typename abstract_domain_t::linear_constraint_t;
  using typename abstract_domain_t::linear_expression_t;
  using typename abstract_domain_t::variable_t;
  using typename abstract_domain_t::variable_vector_t;
  using ptr_cst_t = crab::pointer_constraint<variable_t>;
  using interval_t = interval<number_t>;
  
private:

  using base_dom_vector =  std::vector<Domain>;
  /*
    Bottom is represented by a vector of one element whose value is bottom
    Top is represented by a vector of one element whose value is top
    We don't represent the empty powerset.
   */
  base_dom_vector m_disjuncts;
  
  powerset_domain(bool is_bottom) {
    if (is_bottom) {
      m_disjuncts.push_back(Domain::bottom());
    } else {
      m_disjuncts.push_back(Domain::top());
    }
  }

  powerset_domain(Domain &&dom) {
    m_disjuncts.emplace_back(std::move(dom));
  }
  
  powerset_domain(base_dom_vector &&powerset)
    : m_disjuncts(std::move(powerset)) {
    if (m_disjuncts.size() > Params::max_disjuncts) {
      smash_disjuncts();
    }
  }

  // Remove redundant disjuncts.
  // Expensive operation (quadratic in the number of disjunctions)
  void simplify(base_dom_vector &disjuncts) {
    std::set<unsigned> redundant_set;    
    base_dom_vector res;
    res.reserve(disjuncts.size());
    for (unsigned i=0, sz=disjuncts.size(); i<sz; ++i) {
      Domain &dom = disjuncts[i];
      bool redundant = false;
      for (unsigned j=0; j<sz; ++j) {
	if (i != j && (redundant_set.count(j) <=0) && (dom <= disjuncts[j])) {
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
  
  // powerset is not modified but it cannot passed as const.
  Domain smash_disjuncts(powerset_domain_t &pw) const {
    if (pw.is_bottom()) {
      return Domain::bottom();
    } else if (pw.is_top()) {
      return Domain::top();
    }

    assert(!pw.m_disjuncts.empty());
    Domain res = pw.m_disjuncts[0];
    for (unsigned i=1, sz= pw.m_disjuncts.size();i<sz;++i) {
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

    CRAB_LOG("powerset", crab::outs() << "Smashing the powerset\n" << *this << " into \n";);          
    assert(!m_disjuncts.empty());
    Domain res = m_disjuncts[0];
    for (unsigned i=1, sz= m_disjuncts.size();i<sz;++i) {
      res |= m_disjuncts[i];
    }
    m_disjuncts.clear();
    m_disjuncts.push_back(res);
    CRAB_LOG("powerset", crab::outs() << *this << "\n";);
    return res;
  }

  void insert(base_dom_vector &vec, Domain dom) {
    for (unsigned i=0,sz=vec.size();i<sz;++i) {
      if (dom <= vec[i]) {
	return;
      }
    }
    vec.push_back(dom);
  }

  void append(base_dom_vector &vec1, const base_dom_vector &vec2) {
    for (unsigned i=0, sz=vec2.size();i<sz;++i) {
      insert(vec1, vec2[i]);
    }
  }

  powerset_domain_t powerset_join_with(powerset_domain_t other) {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    }
    base_dom_vector res(m_disjuncts.begin(), m_disjuncts.end());
    append(res, other.m_disjuncts);
    return powerset_domain_t(std::move(res));
  }

  powerset_domain_t powerset_meet_with(powerset_domain_t other) {
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
    for(unsigned i=0, sz_i=m_disjuncts.size(); i<sz_i; ++i) {
      for(unsigned j=0, sz_j=other.m_disjuncts.size(); j<sz_j; ++j) {
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
    m_disjuncts.push_back(Domain::top());
  }
  powerset_domain(const powerset_domain_t &other) = default;
  powerset_domain(powerset_domain_t &&other) = default;
  powerset_domain_t &operator=(const powerset_domain_t &other) = default;
  powerset_domain_t &operator=(powerset_domain_t &&other) = default;  

  virtual bool is_bottom() override {
    for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
      if (!m_disjuncts[i].is_bottom()) {
	return false;
      }
    }
    return true;
  }

  virtual bool is_top() override {
    for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
      if (m_disjuncts[i].is_top()) {
	set_to_top();
	return true;
      }
    }
    return false;
  }

  virtual void set_to_top() override {
    m_disjuncts.clear();
    m_disjuncts.push_back(Domain::top());
  }

  virtual void set_to_bottom() override {
    m_disjuncts.clear();
    m_disjuncts.push_back(Domain::bottom());
  }
  
  
  virtual bool operator<=(powerset_domain_t other) override  {
    Domain left  = smash_disjuncts(*this);
    Domain right = smash_disjuncts(other);
    return left <= right;
  }

  virtual void operator|=(powerset_domain_t other) override {
    CRAB_LOG("powerset", crab::outs() << "JOIN \n" << *this << " and\n" << other << "=\n";);    
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
    if (m_disjuncts.size() > Params::max_disjuncts) {
      smash_disjuncts();
    } 
    CRAB_LOG("powerset", crab::outs() << *this << "\n";);      
  }
  
  virtual powerset_domain_t operator|(powerset_domain_t other) override {
    return powerset_join_with(other);
  }

  virtual powerset_domain_t operator&(powerset_domain_t other) override {
    if (Params::exact_meet) {
      return powerset_meet_with(other);
    } else {
      Domain left  = smash_disjuncts(*this);
      Domain right = smash_disjuncts(other);
      return powerset_domain_t(left & right);
    }
  }

  virtual powerset_domain_t operator||(powerset_domain_t other) override {
    Domain left  = smash_disjuncts(*this);
    Domain right = smash_disjuncts(other);
    return powerset_domain_t(left || right);
  }

  virtual powerset_domain_t widening_thresholds(powerset_domain_t other,
			    const crab::iterators::thresholds<number_t> &ts) override {
    Domain left  = smash_disjuncts(*this);
    Domain right = smash_disjuncts(other);
    return powerset_domain_t(left.widening_thresholds(right, ts));
  }

  virtual powerset_domain_t operator&&(powerset_domain_t other) override {
    Domain left  = smash_disjuncts(*this);
    Domain right = smash_disjuncts(other);
    return powerset_domain_t(left && right);
  }

  virtual void assign(variable_t x, linear_expression_t e) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].assign(x, e);
      }
    }
  }

  virtual void apply(operation_t op, variable_t x, variable_t y, variable_t z) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].apply(op, x, y, z);
      }
    }
  }
  

  virtual void apply(operation_t op, variable_t x, variable_t y, number_t k) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].apply(op, x, y, k);
      }
    }
  }

  virtual void backward_assign(variable_t x, linear_expression_t e,
			       powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");
  }

  virtual void backward_apply(operation_t op, variable_t x, variable_t y, number_t k,
			      powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");    
  }

  virtual void backward_apply(operation_t op, variable_t x, variable_t y, variable_t z,
			      powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");        
  }

  virtual void operator+=(linear_constraint_system_t csts) override {
    if (is_bottom() || csts.is_true()) {
      return;
    }
    if (csts.is_false()) {
      set_to_bottom();
    }
    CRAB_LOG("powerset", crab::outs() << "Adding " << csts << "\n";);
    base_dom_vector vec;
    vec.reserve(m_disjuncts.size());
    for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
      m_disjuncts[i].operator+=(csts);
      if (m_disjuncts[i].is_bottom()) {
	continue;
      }
      vec.emplace_back(m_disjuncts[i]);
    }
    std::swap(m_disjuncts, vec);
    CRAB_LOG("powerset", crab::outs() << "Res=" << *this << "\n";);
  }

  virtual void operator-=(variable_t v) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i] -= v;
	if (m_disjuncts[i].is_top()) {
	  set_to_top();
	}
      }
    }            
  }

  // cast_operators_api
  
  virtual void apply(int_conv_operation_t op, variable_t dst, variable_t src) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].apply(op, dst, src);
      }
    }
  }

  // bitwise_operators_api

  virtual void apply(bitwise_operation_t op, variable_t x, variable_t y, variable_t z) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].apply(op, x, y, z);
      }
    }
  }

  virtual void apply(bitwise_operation_t op, variable_t x, variable_t y, number_t k) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].apply(op, x, y, k);
      }
    }
  }

  // array_operators_api

  virtual void array_init(variable_t a, linear_expression_t elem_size,
                          linear_expression_t lb_idx,
                          linear_expression_t ub_idx,			  
                          linear_expression_t val) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].array_init(a, elem_size, lb_idx, ub_idx, val);
      }                
    }
  }

  virtual void array_load(variable_t lhs, variable_t a,
                          linear_expression_t elem_size,
                          linear_expression_t idx) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].array_load(lhs, a, elem_size, idx);
      }
    }
  }

  virtual void array_store(variable_t a, linear_expression_t elem_size,
                           linear_expression_t idx, linear_expression_t val,
                           bool is_strong_update) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].array_store(a, elem_size, idx, val, is_strong_update);
      }
    }
  }

  virtual void array_store(variable_t a_new, variable_t a_old,
			   linear_expression_t elem_size,
                           linear_expression_t idx, linear_expression_t val,
                           bool is_strong_update) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].array_store(a_new, a_old, elem_size, idx, val, is_strong_update);
      }
    }
  }
  
  virtual void array_store_range(variable_t a, linear_expression_t elem_size,
                                 linear_expression_t lb_idx, linear_expression_t ub_idx,
                                 linear_expression_t val) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].array_store_range(a, elem_size, lb_idx, ub_idx, val);
      }
    }
  }
  
  virtual void array_store_range(variable_t a_new, variable_t a_old,
                                 linear_expression_t elem_size,
                                 linear_expression_t lb_idx, linear_expression_t ub_idx,
                                 linear_expression_t val) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].array_store_range(a_new, a_old, elem_size, lb_idx, ub_idx, val);
      }
    }
  }

  virtual void array_assign(variable_t lhs, variable_t rhs) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].array_assign(lhs, rhs);
      }
    }
  }

  // backward array operations

  virtual void backward_array_init(variable_t a, linear_expression_t elem_size,
                                   linear_expression_t lb_idx,
                                   linear_expression_t ub_idx,
                                   linear_expression_t val, powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                
  }
  
  virtual void backward_array_load(variable_t lhs, variable_t a,
                                   linear_expression_t elem_size,
                                   linear_expression_t idx, powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                
  }
  
  virtual void backward_array_store(variable_t a, linear_expression_t elem_size,
                                    linear_expression_t idx,
                                    linear_expression_t v,
                                    bool is_strong_update, powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                
  }
  
  virtual void backward_array_store(variable_t a_new, variable_t a_old,
                                    linear_expression_t elem_size,
                                    linear_expression_t idx,
                                    linear_expression_t v,
                                    bool is_strong_update, powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                
  }
  
  virtual void
  backward_array_store_range(variable_t a, linear_expression_t elem_size,
                             linear_expression_t lb_idx, linear_expression_t ub_idx,
                             linear_expression_t v, powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                    
  }
  
  virtual void backward_array_store_range(variable_t a_new, variable_t a_old,
                                          linear_expression_t elem_size,
                                          linear_expression_t lb_idx,
                                          linear_expression_t ub_idx,
                                          linear_expression_t v,
                                          powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                
  }
  
  virtual void backward_array_assign(variable_t a, variable_t b,
                                     powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");     
  }

  
  // pointer_operators_api
  virtual void pointer_load(variable_t lhs, variable_t rhs, linear_expression_t elem_size) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_load(lhs, rhs, elem_size);
      }                        
    }
  }

  virtual void pointer_store(variable_t lhs, variable_t rhs, linear_expression_t elem_size) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_store(lhs, rhs, elem_size);
      }                            
    }
  }

  virtual void pointer_assign(variable_t lhs, variable_t rhs, linear_expression_t offset) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_assign(lhs, rhs, offset);
      }
    }
  }

  virtual void pointer_mk_obj(variable_t lhs, ikos::index_t address) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_mk_obj(lhs, address);
      }
    }
  }

  virtual void pointer_function(variable_t lhs, varname_t func) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_function(lhs, func);
      }
    }
  }

  virtual void pointer_mk_null(variable_t lhs) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_mk_null(lhs);
      }
    }
  }

  virtual void pointer_assume(ptr_cst_t cst) override {
    if (!is_bottom()) {
      base_dom_vector vec;
      vec.reserve(m_disjuncts.size());
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_assume(cst);
	if (m_disjuncts[i].is_bottom()) {
	  continue;
	}
	vec.emplace_back(m_disjuncts[i]);      
      }
      std::swap(m_disjuncts, vec);
    }
  }

  virtual void pointer_assert(ptr_cst_t cst) override {
    if (!is_bottom()) {
      base_dom_vector vec;
      vec.reserve(m_disjuncts.size());    
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].pointer_assert(cst);
	if (m_disjuncts[i].is_bottom()) {
	continue;
	}
	vec.emplace_back(m_disjuncts[i]);            
      }
      std::swap(m_disjuncts, vec);
    }
  }

  // boolean operators
  virtual void assign_bool_cst(variable_t lhs,
                               linear_constraint_t rhs) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].assign_bool_cst(lhs, rhs);
      }
    }
  }

  virtual void assign_bool_var(variable_t lhs, variable_t rhs,
                               bool is_not_rhs) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].assign_bool_var(lhs, rhs, is_not_rhs);
      }
    }
  }

  virtual void apply_binary_bool(bool_operation_t op, variable_t x,
                                 variable_t y, variable_t z) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].apply_binary_bool(op, x, y, z);
      }
    }
  }

  virtual void assume_bool(variable_t v, bool is_negated) override {
    if (!is_bottom()) {
      base_dom_vector vec;
      vec.reserve(m_disjuncts.size());        
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].assume_bool(v, is_negated);
	if (m_disjuncts[i].is_bottom()) {
	  continue;
	}
	vec.emplace_back(m_disjuncts[i]);
      }
      std::swap(m_disjuncts, vec);    
    }
  }

  // backward boolean operators
  virtual void backward_assign_bool_cst(variable_t lhs, linear_constraint_t rhs,
                                        powerset_domain_t inv) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");            
  }

  virtual void backward_assign_bool_var(variable_t lhs, variable_t rhs,
                                        bool is_not_rhs,
                                        powerset_domain_t inv) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                
  }

  virtual void backward_apply_binary_bool(bool_operation_t op, variable_t x,
                                          variable_t y, variable_t z,
                                          powerset_domain_t inv) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");                    
  }

  // Intrinsics

  virtual void intrinsic(std::string name,
			 const variable_vector_t &inputs,
			 const variable_vector_t &outputs) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].intrinsic(name, inputs, outputs);
      }
    }
  }

  virtual void backward_intrinsic(std::string name,
				  const variable_vector_t &inputs,
				  const variable_vector_t &outputs,
				  powerset_domain_t invariant) override {
    CRAB_WARN(getDomainName(), " does not implement backward operations");            
  }
  
  // Miscellaneous operations 

  interval_t operator[](variable_t v) {
    Domain smashed = smash_disjuncts(*this);
    return smashed[v];
  }

  void set(variable_t v, interval_t intv) {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].set(v, intv);
      }
    }
  }
  
  virtual void normalize() override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].normalize();
      }
    }
    
  }
  
  virtual void minimize() override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].minimize();
      }
    }
  }
  
  virtual void rename(const variable_vector_t &from,
                      const variable_vector_t &to) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].rename(from, to);
      }
    }
  }


  virtual void expand(variable_t x, variable_t new_x) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].expand(x, new_x);
      }
    }
  }

  virtual void forget(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].forget(variables);
	if (m_disjuncts[i].is_top()) {
	  set_to_top();
	}
      }
    }
  }

  virtual void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	m_disjuncts[i].project(variables);
      }
    }
  }
  
  linear_constraint_system_t to_linear_constraint_system() {
    Domain smashed = smash_disjuncts(*this);
    return smashed.to_linear_constraint_system();
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() {
    if (is_bottom()) {
      disjunctive_linear_constraint_system_t res(true);
      return res;
    } else if (is_top()) {
      disjunctive_linear_constraint_system_t res(false);
      return res;
    } else {
      disjunctive_linear_constraint_system_t res(true);
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz; ++i) {
	res += m_disjuncts[i].to_linear_constraint_system();
      }
      return res;
    } 
  }
  
  void write(crab::crab_os &o) {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      for (unsigned i=0, sz=m_disjuncts.size(); i<sz;) {
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
    for (unsigned i=0, sz=m_disjuncts.size(); i<sz;++i) {
      crab::outs() << m_disjuncts[i] << " || ";
    }
    crab::outs() << "== End powerset internal representation === \n";    
  }
  
  static std::string getDomainName() {
    return std::string("Powerset(") + Domain::getDomainName() + ")";
  }

}; 
  
template <typename Domain, class Params>
struct abstract_domain_traits<powerset_domain<Domain, Params>> {
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
};

} // namespace domains
} // namespace crab
