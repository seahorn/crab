#pragma once

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/union_find_domain.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>
#include <string>

namespace crab {
namespace domains {

/*
 * Variable packing domain for numerical domains.
 *
 * Inspired by the paper "Precise and Efficient Static Array Bound
 * Checking for Large Embedding C Programs" by A. Venet and G. Brat.
 *
 * The idea is to partition the set of numerical variables into
 * multiple packs and model each pack separately using an abstract
 * state of type NumDom. The partitioning is done dynamically using an
 * union-find datastructure which has been adapted to have lattice
 * operations. Note that we use intersection semantics in the
 * union-find because the concretization of the variable packing
 * domain is the intersection of the concretization of each pack.
 */
template <class NumDom>
class numerical_packing_domain
    : public abstract_domain_api<numerical_packing_domain<NumDom>> {
public:
  using this_type = numerical_packing_domain<NumDom>;
  using abstract_domain_api_t = abstract_domain_api<this_type>;
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

private:
  /* begin type definitions */
  using base_domain_t = NumDom;

  struct uf_forget_element_in_domain {
    void operator()(base_domain_t &dom, const variable_t &v) const { dom -= v; }
  };
  struct uf_rename_element_in_domain {
    void operator()(base_domain_t &dom, const variable_t &x,
                    const variable_t &y) const {
      std::vector<variable_t> old_vars{x};
      std::vector<variable_t> new_vars{y};
      dom.rename(old_vars, new_vars);
    }
  };
  using union_find_domain_t = union_find_domain<
      variable_t, base_domain_t, uf_intersection_semantics<base_domain_t>,
      uf_forget_element_in_domain, uf_rename_element_in_domain>;
  using pack_t = typename union_find_domain_t::equivalence_class_t;
  using pack_vars_t = typename union_find_domain_t::equivalence_class_elems_t;
  /* end type definitions */

  union_find_domain_t m_packs;

  // return null if the merge causes bottom.
  std::shared_ptr<base_domain_t> merge(const variable_vector_t &vars) {
    assert(!m_packs.is_bottom());
    // shouldn't be top even if the whole abstract value is top.
    assert(!m_packs.is_top());
    assert(!vars.empty());

    std::shared_ptr<base_domain_t> top_val = std::make_shared<base_domain_t>();

    variable_t x = vars[0];
    if (!m_packs.contains(x)) {
      m_packs.make(x, top_val);
    }

    for (unsigned i = 1, sz = vars.size(); i < sz; ++i) {
      if (!m_packs.contains(vars[i])) {
        m_packs.make(vars[i], top_val);
      }
      // note that this join has intersection semantics
      if (!m_packs.join(x, vars[i])) {
        // join returns none if the result of join is bottom
        return nullptr;
      }
    }

    pack_t &pack = m_packs.get_equiv_class(x);
    return pack.detach_and_get_absval();
  }

  base_domain_t merge(const variable_vector_t &vars) const {
    assert(!m_packs.is_bottom());
    assert(!m_packs.is_top());
    assert(!vars.empty());
    base_domain_t base;
    for (auto const &v : vars) {
      if (m_packs.contains(v)) {
        base &= *(m_packs.get_equiv_class(v).get_absval());
      }
    }
    return std::move(base);
  }

  // x := f(y)
  std::shared_ptr<base_domain_t> apply_packs(const variable_t &x, const variable_t &y) {
    std::shared_ptr<base_domain_t> absval = nullptr;
    if (x != y) {
      m_packs.forget(x);
      variable_vector_t vars{x, y};
      absval = merge(vars);
    } else {
      variable_vector_t vars{x};
      absval = merge(vars);
    }
    assert(absval);
    return absval;
  }

  // x := f(y,z)
  std::shared_ptr<base_domain_t> apply_packs(const variable_t &x, const variable_t &y, const variable_t &z) {
    std::shared_ptr<base_domain_t> absval = nullptr;
    if (x != y && x != z) {
      m_packs.forget(x);
      variable_vector_t vars{x, y, z};
      absval = merge(vars);
    } else {
      if (x == y) {
	variable_vector_t vars{x, z};
	absval = merge(vars);
      } else {
	variable_vector_t vars{x, y};
	absval = merge(vars);
      } 
    }
    assert(absval);
    return absval;
  }
  
  numerical_packing_domain(union_find_domain_t &&packs)
      : m_packs(std::move(packs)) {}

public:
  numerical_packing_domain(bool is_bottom = false)
      : m_packs(is_bottom ? union_find_domain_t::bottom()
                          : union_find_domain_t()) {}

  numerical_packing_domain(const this_type &o) = default;
  numerical_packing_domain(this_type &&o) = default;
  this_type &operator=(const this_type &o) = default;
  this_type &operator=(this_type &&o) = default;

  void set_to_top() override {
    m_packs = union_find_domain_t(); // empty union-find
  }

  void set_to_bottom() override { m_packs.set_to_bottom(); }

  this_type make_bottom() const override {
    this_type res(true);
    return res;
  }

  this_type make_top() const override {
    this_type res(false);
    return res;
  }

  bool is_bottom() const override { return m_packs.is_bottom(); }

  bool is_top() const override {
    if (is_bottom()) {
      return false;
    }
    for (std::shared_ptr<const base_domain_t> dom : m_packs.domains()) {
      if (!dom->is_top()) {
        return false;
      }
    }
    return true;
  }

  bool operator<=(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return true;
    } else if (is_top() || other.is_bottom()) {
      return false;
    } else {
      return m_packs <= other.m_packs;
    }
  }

  void operator|=(const this_type &other) override {
    if (is_bottom() || other.is_top()) {
      *this = other;
    } else if (other.is_bottom() || is_top()) {
      // do nothing
    } else {
      m_packs = m_packs | other.m_packs;
    }
  }

  this_type operator|(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      return this_type(m_packs | other.m_packs);
    }
  }

  void operator&=(const this_type &other) override {
    if (is_bottom() || other.is_top()) {
      // do nothing
    } else if (other.is_bottom() || is_top()) {
      *this = other;
    } else {
      m_packs = m_packs & other.m_packs;
    }
  }

  this_type operator&(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      return this_type(m_packs & other.m_packs);
    }
  }

  this_type operator||(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      return this_type(m_packs || other.m_packs);
    }
  }

  this_type
  widening_thresholds(const this_type &other,
                      const thresholds<number_t> & /*ts*/) const override {
    if (is_bottom() || other.is_top()) {
      return other;
    } else if (other.is_bottom() || is_top()) {
      return *this;
    } else {
      // TODO: widening w/ thresholds
      return this_type(m_packs || other.m_packs);
    }
  }

  this_type operator&&(const this_type &other) const override {
    if (is_bottom() || other.is_top()) {
      return *this;
    } else if (other.is_bottom() || is_top()) {
      return other;
    } else {
      return this_type(m_packs && other.m_packs);
    }
  }

  void operator-=(const variable_t &var) override {
    if (!(is_bottom() || is_top())) {
      if (m_packs.contains(var)) {
        pack_t &pack = m_packs.get_equiv_class(var);
        pack.detach_and_get_absval()->operator-=(var);
      }
      m_packs.forget(var);
    }
  }

  interval_t operator[](const variable_t &v) override {
    if (is_bottom()) {
      return interval_t::bottom();
    }

    if (is_top() || !m_packs.contains(v)) {
      return interval_t::top();
    }

    pack_t &pack = m_packs.get_equiv_class(v);
    // it might trigger normalization in the base domain
    return pack.detach_and_get_absval()->operator[](v);
  }

  interval_t at(const variable_t &v) const override {
    if (is_bottom()) {
      return interval_t::bottom();
    }

    if (is_top() || !m_packs.contains(v)) {
      return interval_t::top();
    }

    const pack_t &pack = m_packs.get_equiv_class(v);
    return pack.get_absval()->at(v);
  }

  void operator+=(const linear_constraint_system_t &csts) override {
    if (!is_bottom()) {
      for (auto const &cst : csts) {
        if (cst.is_contradiction()) {
          set_to_bottom();
          return;
        }
        if (cst.is_tautology()) {
          continue;
        }
        variable_vector_t vars(cst.variables().begin(), cst.variables().end());
        if (!vars.empty()) {
          if (std::shared_ptr<base_domain_t> absval = merge(vars)) {
            *absval += cst;
            if (absval->is_bottom()) {
              set_to_bottom();
              return;
            }
          }
        }
      }
    }
  }

  bool entails(const linear_constraint_t &cst) const override {
    if (is_bottom()) {
      return true;
    } else if (cst.is_tautology()) {
      return true;
    } else if (cst.is_contradiction()) {
      return false;
    }
    variable_vector_t vars(cst.variables().begin(), cst.variables().end());
    if (!vars.empty()) {
      base_domain_t absval = merge(vars);
      return absval.entails(cst);
    }
    return false;
  }

  void assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      std::shared_ptr<base_domain_t> absval = nullptr;
      variable_vector_t vars(e.variables().begin(), e.variables().end());      
      if (std::find(vars.begin(), vars.end(), x) == vars.end()) {
	m_packs.forget(x);
	vars.push_back(x);
      }       
      absval = merge(vars);
      if (absval) {
        absval->assign(x, e);
      } else {
	CRAB_ERROR(domain_name(), "::assign produced bottom!");
      }
      
    }
  }

  void weak_assign(const variable_t &x, const linear_expression_t &e) override {
    if (!is_bottom()) {
      std::shared_ptr<base_domain_t> absval = nullptr;
      variable_vector_t vars(e.variables().begin(), e.variables().end());      
      if (std::find(vars.begin(), vars.end(), x) == vars.end()) {
	m_packs.forget(x);
	vars.push_back(x);
      }       
      absval = merge(vars);
      if (absval) {
        absval->weak_assign(x, e);
      } else {
	CRAB_ERROR(domain_name(), "::weak_assign produced bottom!");
      }
    }
  }
  
  
  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    if (!is_bottom()) {
      if (std::shared_ptr<base_domain_t> absval = apply_packs(x, y)) {
	absval->apply(op, x, y, z);
      } else {
	CRAB_ERROR(domain_name(), "::apply 1 produced bottom!");
      }
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      if (std::shared_ptr<base_domain_t> absval = apply_packs(x, y, z)) {
	absval->apply(op, x, y, z);
      } else {
	CRAB_ERROR(domain_name(), "::apply 2 produced bottom!");
      }
    }
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    if (!is_bottom() && (src != dst)) {
      m_packs.forget(dst);
      variable_vector_t vars{src, dst};
      if (std::shared_ptr<base_domain_t> absval = merge(vars)) {
        absval->apply(op, dst, src);
      } else {
	CRAB_ERROR(domain_name(), "::apply 3 produced bottom!");
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    if (!is_bottom()) {
      if (std::shared_ptr<base_domain_t> absval = apply_packs(x, y)) {
        absval->apply(op, x, y, k);
      } else {
	CRAB_ERROR(domain_name(), "::apply 4 produced bottom!");
      }
    }
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    if (!is_bottom()) {
      if (std::shared_ptr<base_domain_t> absval = apply_packs(x, y, z)) {
        absval->apply(op, x, y, z);
      } else {
	CRAB_ERROR(domain_name(), "::apply 5 produced bottom!");
      }
    }
  }
  
  void select(const variable_t &lhs, const linear_constraint_t &cond,
              const linear_expression_t &e1,
              const linear_expression_t &e2) override {
    if (!is_bottom()) {
      m_packs.forget(lhs);

      variable_vector_t vars;
      size_t num_vars =
          std::distance(cond.variables().begin(), cond.variables().end()) +
          std::distance(e1.variables_begin(), e1.variables_end()) +
          std::distance(e2.variables_begin(), e2.variables_end());
      vars.reserve(num_vars + 1);

      vars.push_back(lhs);
      vars.insert(vars.end(), cond.variables().begin(), cond.variables().end());
      vars.insert(vars.end(), e1.variables_begin(), e1.variables_end());
      vars.insert(vars.end(), e2.variables_begin(), e2.variables_end());

      if (std::shared_ptr<base_domain_t> absval = merge(vars)) {
        absval->select(lhs, cond, e1, e2);
      }
    }
  }

  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_assign not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_apply not implemented");
  }

  BOOL_OPERATIONS_NOT_IMPLEMENTED(this_type)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(this_type)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(this_type)

  linear_constraint_system_t to_linear_constraint_system() const override {

    if (is_bottom()) {
      return linear_constraint_system_t(linear_constraint_t::get_false());
    }

    if (is_top()) {
      return linear_constraint_system_t(linear_constraint_t::get_true());
    }

    linear_constraint_system_t res;
    for (std::shared_ptr<const base_domain_t> absval : m_packs.domains()) {
      res += absval->to_linear_constraint_system();
    }
    return res;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    CRAB_WARN(domain_name(),
              "::to_disjunctive_linear_constraint_system not implemented");
    disjunctive_linear_constraint_system_t res;
    return res;
  }

  void forget(const variable_vector_t &variables) override {
    if (!(is_bottom() || is_top())) {
      for (const variable_t &v : variables) {
        if (m_packs.contains(v)) {
          pack_t &pack = m_packs.get_equiv_class(v);
          pack.detach_and_get_absval()->operator-=(v);
        }
        m_packs.forget(v);
      }
    }
  }

  void project(const variable_vector_t &variables) override {
    if (!is_bottom()) {
      for (std::shared_ptr<base_domain_t> absval : m_packs.domains()) {
        absval->project(variables);
      }
      m_packs.project(variables);
    }
  }

  void expand(const variable_t &var, const variable_t &new_var) override {
    if (!is_bottom()) {
      m_packs.forget(new_var);
      variable_vector_t vars{var, new_var};
      if (std::shared_ptr<base_domain_t> absval = merge(vars)) {
        absval->expand(var, new_var);
      }
    }
  }

  void normalize() override {}
  void minimize() override {}

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    if (!is_bottom()) {
      for (std::shared_ptr<base_domain_t> absval : m_packs.domains()) {
        absval->rename(from, to);
      }
      m_packs.rename(from, to);
    }
  }

  void intrinsic(std::string name, const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {

    auto error_if_not_variable = [&name](const variable_or_constant_t &vc) {
      if (!vc.is_variable()) {
        CRAB_ERROR("Intrinsics ", name, " expected a variable input");
      }
    };

    if (is_bottom()) {
      return;
    }

    if (name == "var_packing_merge") {
      // dom.intrinsic("var_packing_merge", {v1, v2, ...}, {})

      if (!outputs.empty()) {
        CRAB_ERROR("Intrinsics ", name,
                   " unexpected number of output parameters");
      }
      variable_vector_t vars;
      vars.reserve(inputs.size());
      for (auto &v_or_c : inputs) {
        error_if_not_variable(v_or_c);
        variable_t v = v_or_c.get_variable();
        vars.push_back(v);
      }
      if (!vars.empty()) {
        merge(vars);
      }
    } else {
      for (std::shared_ptr<base_domain_t> absval : m_packs.domains()) {
        absval->intrinsic(name, inputs, outputs);
      }
    }
  }

  void backward_intrinsic(std::string name,
                          const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const this_type &invariant) override {
    CRAB_WARN(domain_name(), "::backward_intrinsic not implemented");
  }

  void write(crab_os &o) const override {
    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "top";
    } else {
      o << "{";
      pack_vars_t packs_vars = m_packs.equiv_classes_elems();
      bool first_pack = true;

      /// Sort to print deterministically
      std::vector<typename union_find_domain_t::element_t> sorted_pack_reps;
      sorted_pack_reps.reserve(packs_vars.size());
      for (auto kv : packs_vars) {
	sorted_pack_reps.push_back(kv.first);
      }
      std::sort(sorted_pack_reps.begin(), sorted_pack_reps.end());

      for (auto rep: sorted_pack_reps) {
	auto it = packs_vars.find(rep);
	if (it == packs_vars.end()) {
	  continue;
	}
	
        if (!first_pack) {
          o << ",";
        } else {
          first_pack = false;
        }
        const pack_t &pack = m_packs.get_equiv_class(rep);
        typename union_find_domain_t::element_set_t pack_vars = it->second;
        o << "Pack({";
        bool first_e = true;
        for (auto const &v : pack_vars) {
          if (!first_e) {
            o << ",";
          } else {
            first_e = false;
          }
          o << v;
        }
        o << "}," << *(pack.get_absval()) << ")";
      }
    }
  }

  friend crab_os &operator<<(crab_os &o, const this_type &packs) {
    packs.write(o);
    return o;
  }

  std::string domain_name() const override {
    base_domain_t absval;
    return "NumPackDomain(" + absval.domain_name() + ")";
  }
};

template <typename Domain>
struct abstract_domain_traits<numerical_packing_domain<Domain>> {
  using number_t = typename Domain::number_t;
  using varname_t = typename Domain::varname_t;
};

} // end namespace domains
} // end namespace crab
