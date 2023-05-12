#pragma once

#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/region/ghost_variables.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/os.hpp>

#include <algorithm>
#include <string>
#include <unordered_map>

namespace crab {
namespace domains {
namespace region_domain_impl {

/* API for a ghost variable manager
   typedef ghost_variables_t;

   ghost_variables_t get_or_insert(const variable_t &v);
   boost::optional<ghost_variables_t> get(const variable_t &v) const;

   static bool need_renaming();
   ghost_var_manager_t join(const ghost_var_manager_t &right_man,
                           GhostDomain &left_val, GhostDomain &right_val) const;
   ghost_var_manager_t meet(const ghost_var_manager_t &right_man,
                                 GhostDomain &left_val, GhostDomain &right_val)
   const; void project_and_meet(const variable_vector_t &variables, const
   ghost_var_man_t &right_man, GhostDomain &left_val, const GhostDomain
   &right_val); void forget(const variable_t &var, GhostDomain &val); void
   project(const variable_vector_t &vars, GhostDomain &val); void rename(const
   variable_vector_t &old_vars, const variable_vector_t &new_vars, GhostDomain
   &val); void write(crab_os &os, const GhostDomain &val) const;

   ghost_variables_t dup(const ghost_variables_t &gvars);

   ghost_linear_constraint_t
   ghosting_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst,
                                  ghost_variable_kind kind);

   ghost_variable_or_constant_t
   rename_variable_or_constant(const variable_or_constant_t &v);
   ghost_linear_expression_t rename_linear_expr(const linear_expression_t &e);
   ghost_linear_constraint_t rename_linear_cst(const linear_constraint_t &cst);
   ghost_linear_constraint_system_t
   rename_linear_cst_sys(const linear_constraint_system_t &csts);

   boost::optional<variable_t> rev_rename_var(const ghost_variable_t &v,
                                              bool ignore_references) const;
   boost::optional<linear_expression_t>
   rev_rename_linear_expr(const ghost_linear_expression_t &e,
                          bool ignore_references) const;
   boost::optional<linear_constraint_t>
   rev_rename_linear_cst(const ghost_linear_constraint_t &cst,
                         bool ignore_references) const;
 */

template <typename Domain, typename GhostDomain>
class ghost_variable_manager_with_fixed_naming {
  using number_t = typename Domain::number_t;
  using ghost_var_manager_t =
      ghost_variable_manager_with_fixed_naming<Domain, GhostDomain>;
  // typedefs for Domain
  using linear_constraint_system_t =
      typename Domain::linear_constraint_system_t;
  using linear_constraint_t = typename Domain::linear_constraint_t;
  using linear_expression_t = typename Domain::linear_expression_t;
  using reference_constraint_t = typename Domain::reference_constraint_t;
  using variable_or_constant_t = typename Domain::variable_or_constant_t;
  using variable_t = typename Domain::variable_t;
  using variable_vector_t = typename Domain::variable_vector_t;
  using varname_t = typename Domain::varname_t;
  // typedefs for GhostDomain
  using ghost_linear_constraint_system_t =
      typename GhostDomain::linear_constraint_system_t;
  using ghost_linear_constraint_t = typename GhostDomain::linear_constraint_t;
  using ghost_linear_expression_t = typename GhostDomain::linear_expression_t;
  using ghost_reference_constraint_t =
      typename GhostDomain::reference_constraint_t;
  using ghost_variable_or_constant_t =
      typename GhostDomain::variable_or_constant_t;
  using ghost_variable_t = typename GhostDomain::variable_t;
  using ghost_variable_vector_t = typename GhostDomain::variable_vector_t;
  using ghost_varname_t = typename GhostDomain::varname_t;
  using GhostDomainVarAlloc = typename ghost_varname_t::variable_factory_t;
  static_assert(std::is_same<variable_t, ghost_variable_t>::value,
                "Ghost variables are created by a fixed scheme so varname must "
                "have same type");

public:
  using ghost_variables_t =
      ghost_variables<Domain, GhostDomain, GhostDomainVarAlloc>;
  using ghost_variable_kind = typename ghost_variables_t::ghost_variable_kind;  
  using get_type_fn = std::function<variable_type(const variable_t &)>;

private:
  // function to get the type of each variable
  get_type_fn m_get_type;

public:
  ghost_variable_manager_with_fixed_naming(get_type_fn get_type)
      : m_get_type(get_type) {}
  ghost_variable_manager_with_fixed_naming(const ghost_var_manager_t &o) =
      default;
  ghost_variable_manager_with_fixed_naming(ghost_var_manager_t &&o) = default;
  ghost_var_manager_t &operator=(const ghost_var_manager_t &o) = default;
  ghost_var_manager_t &operator=(ghost_var_manager_t &&o) = default;

  ghost_variables_t get_or_insert(const variable_t &v) {
    auto gvars_opt = get(v);
    assert(gvars_opt);
    return *gvars_opt;
  }

  boost::optional<ghost_variables_t> get(const variable_t &v) const {
    auto &vfac = const_cast<varname_t *>(&(v.name()))->get_var_factory();
    auto alloc_fn = [](GhostDomainVarAlloc &alloc, const varname_t &name,
                       const std::string &role) {
      if (role == "") {
        return name;
      } else {
        return alloc.get(name, role);
      }
    };
    return ghost_variables_t::create(vfac, alloc_fn, v.name(), m_get_type(v));
  }

  static bool need_renaming() { return false; }

  ghost_var_manager_t join(const ghost_var_manager_t &, GhostDomain &,
                           GhostDomain &) const {
    CRAB_ERROR("This manager does not provide join operation");
  }

  ghost_var_manager_t meet(const ghost_var_manager_t &, GhostDomain &,
                           GhostDomain &) const {
    CRAB_ERROR("This manager does not provide meet operation");
  }

  // left_val = meet(left_val, project(right_val, variables))
  void merge(const variable_vector_t &variables,
	     const ghost_var_manager_t &, GhostDomain &left_val,
	     const GhostDomain &right_val) {

    if (!variables.empty()) {
      GhostDomain copy_right(right_val);
      copy_right.project(variables);
      left_val.forget(variables);
      left_val = left_val & copy_right;
    }
  }

  void forget(const variable_t &var, GhostDomain &val) {
    if (boost::optional<ghost_variables_t> gvars = get(var)) {
      (*gvars).forget(val);
    }
  }

  void project(const variable_vector_t &vars, GhostDomain &val) {
    ghost_variable_vector_t gvars_vec;
    for (auto const &v : vars) {
      if (boost::optional<ghost_variables_t> gvars = get(v)) {
        (*gvars).add(gvars_vec);
      }
    }

    val.project(gvars_vec);
  }

  // Note that dup does not return a fresh set of ghost variables each
  // time it is called. Instead, it returns the same set of ghost
  // variables.
  ghost_variables_t dup(const ghost_variables_t &gvars) {
    auto &vfac = gvars.get_vfac();
    auto alloc_fn = [](GhostDomainVarAlloc &alloc, const ghost_varname_t &name,
                       const std::string &role) {
      return alloc.get(name, role + ".dup");
    };
    return ghost_variables_t::create(vfac, alloc_fn, gvars);
  }

  void rename(const variable_vector_t &old_vars,
              const variable_vector_t &new_vars, GhostDomain &val) {

    ghost_variable_vector_t old_ghost_vars, new_ghost_vars;
    for (unsigned i = 0, sz = old_vars.size(); i < sz; ++i) {
      variable_t old_var = old_vars[i];
      if (boost::optional<ghost_variables_t> old_gvars = get(old_var)) {
        (*old_gvars).add(old_ghost_vars);
        variable_t new_var = new_vars[i];
        ghost_variables_t new_gvars = get_or_insert(new_var);
        new_gvars.add(new_ghost_vars);
      } else {
        continue;
      }
    }
    val.rename(old_ghost_vars, new_ghost_vars);
  }

  // The ghost domain doesn't know about references so we convert a
  // reference constraint into a linear constraint.
  // A reference variable can be ghosted by multiple ghost variables.
  // The parameter @kind tells which one.  
  ghost_linear_constraint_t
  ghosting_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst,
				 ghost_variable_kind kind) {

    auto get_ghost_var = [this](const variable_t &var, ghost_variable_kind kind) ->
      boost::optional<ghost_variable_t> {
      auto ghost_vars  = get_or_insert(var);
      if (kind == ghost_variable_kind::ADDRESS) {
	return ghost_vars.get_var();
      } else if (kind == ghost_variable_kind::OFFSET) {
	if (ghost_vars.has_offset_and_size()) {
	  return ghost_vars.get_offset_and_size().get_offset();
	}
      } else if (kind == ghost_variable_kind::SIZE) {
	if (ghost_vars.has_offset_and_size()) {
	  return ghost_vars.get_offset_and_size().get_size();
	}
      }
      return boost::none;
    };
      
    if (ref_cst.is_tautology()) {
      return ghost_linear_constraint_t::get_true();
    } else if (ref_cst.is_contradiction()) {
      return ghost_linear_constraint_t::get_false();
    } else {
      if (ref_cst.is_unary()) {
        assert(ref_cst.lhs().get_type().is_reference());
        ghost_variable_t x = get_or_insert(ref_cst.lhs()).get_var();
        if (ref_cst.is_equality()) {
          return ghost_linear_constraint_t(x == number_t(0));
        } else if (ref_cst.is_disequality()) {
          return ghost_linear_constraint_t(x != number_t(0));
        } else if (ref_cst.is_less_or_equal_than()) {
          return ghost_linear_constraint_t(x <= number_t(0));
        } else if (ref_cst.is_less_than()) {
          return ghost_linear_constraint_t(x < number_t(0));
        } else if (ref_cst.is_greater_or_equal_than()) {
          return ghost_linear_constraint_t(x >= number_t(0));
        } else if (ref_cst.is_greater_than()) {
          return ghost_linear_constraint_t(x > number_t(0));
        }
      } else {
        assert(ref_cst.lhs().get_type().is_reference());
        assert(ref_cst.rhs().get_type().is_reference());
	boost::optional<ghost_variable_t> x = get_ghost_var(ref_cst.lhs(), kind);
	boost::optional<ghost_variable_t> y = get_ghost_var(ref_cst.rhs(), kind);
	if (!x || !y) {
	  return ghost_linear_constraint_t::get_true();
	}
        number_t offset = ref_cst.offset();
        if (ref_cst.is_equality()) {
          return ghost_linear_constraint_t(*x == *y + offset);
        } else if (ref_cst.is_disequality()) {
          return ghost_linear_constraint_t(*x != *y + offset);
        } else if (ref_cst.is_less_or_equal_than()) {
          return ghost_linear_constraint_t(*x <= *y + offset);
        } else if (ref_cst.is_less_than()) {
          return ghost_linear_constraint_t(*x < *y + offset);
        } else if (ref_cst.is_greater_or_equal_than()) {
          return ghost_linear_constraint_t(*x >= *y + offset);
        } else if (ref_cst.is_greater_than()) {
          return ghost_linear_constraint_t(*x > *y + offset);
        }
      }
    }
    CRAB_ERROR("unexpected reference constraint");
  }

  void write(crab_os &os, const GhostDomain &val) const { os << val; }

  /*
   * These functions must be provided because Domain and GhostDomain
   * can have different variable factories with different varname
   * types
   */

  ghost_variable_or_constant_t
  rename_variable_or_constant(const variable_or_constant_t &v) {
    return std::move(v);
  }

  ghost_linear_expression_t rename_linear_expr(const linear_expression_t &e) {
    return std::move(e);
  }

  ghost_linear_constraint_t rename_linear_cst(const linear_constraint_t &cst) {
    return std::move(cst);
  }

  ghost_linear_constraint_system_t
  rename_linear_cst_sys(const linear_constraint_system_t &csts) {
    return std::move(csts);
  }

  boost::optional<variable_t> rev_rename_var(const ghost_variable_t &v,
                                             bool ignore_references) const {

    if (!ignore_references || !v.get_type().is_reference()) {
      return v;
    } else {
      return boost::none;
    }
  }

  boost::optional<linear_expression_t>
  rev_rename_linear_expr(const ghost_linear_expression_t &e,
                         bool ignore_references) const {
    linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const ghost_variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      if (boost::optional<variable_t> v_opt =
              rev_rename_var(v, ignore_references)) {
        out = out + (coef * (*v_opt));
      } else {
        return boost::optional<linear_expression_t>();
      }
    }
    return out;
  }

  boost::optional<linear_constraint_t>
  rev_rename_linear_cst(const ghost_linear_constraint_t &cst,
                        bool ignore_references) const {
    if (boost::optional<linear_expression_t> e =
            rev_rename_linear_expr(cst.expression(), ignore_references)) {

      return linear_constraint_t(
            *e, (typename linear_constraint_t::kind_t)cst.kind());
    } else {
      return boost::optional<linear_constraint_t>();
    }
  }
};

/* Ghost variable manager that allows ghost names to change and allow
 * Domain and GhostDomain to use different variable factories.
 *
 * It assumes that abstract operations are defined over Domain and
 * there is another abstract domain of type GhostDomain to which
 * Domain's abstract transformers are redirected to. This manager also
 * assumes that Domain and GhostDomain use different varname factories
 * which causes a lot of conversions back and
 * forth. GhostDomainVarAlloc is the varname factory for GhostDomain.
 *
 * BREAK OF MODULARITY:
 *
 * However, this ghosting scheme can break the modularity of the
 * subdomains because ghost variables can change after binary
 * operations. A common domain configuration is:
 *
 *                                     kind of ghosting (none, variable, fixed)
 *  ------------------------------
 * | region domain               |     variable
 *  ------------------------------
 * | array_adaptive domain       |     variable on non-smashed arrays and none for smashed 
 *  ------------------------------
 * | array_smashing              |     fixed
 *  ------------------------------
 * | product(boolean, numerical) |     none
 *  ------------------------------
 * | product of numerical domains|     none
 *  ------------------------------
 *
 * This domain configuration breaks the array_smashing domain because
 * it assumes that array variables cannot be renamed. However,
 * ghost_variable_manager_with_variable_naming renames array
 * variables.
 *
 * On the other hand, this configuration is safe:
 *
 *  ------------------------------
 * | array_adaptive domain       |
 *  ------------------------------
 * | array_smashing              |
 *  ------------------------------
 * | product(boolean, numerical) |
 *  ------------------------------
 * | product of numerical domains|
 *  ------------------------------
 */
template <typename Domain, typename GhostDomain, typename GhostDomainVarAlloc>
class ghost_variable_manager_with_variable_naming {

  using number_t = typename Domain::number_t;
  using ghost_var_manager_t =
      ghost_variable_manager_with_variable_naming<Domain, GhostDomain,
                                                  GhostDomainVarAlloc>;

  // typedefs for Domain
  using linear_constraint_system_t =
      typename Domain::linear_constraint_system_t;
  using linear_constraint_t = typename Domain::linear_constraint_t;
  using linear_expression_t = typename Domain::linear_expression_t;
  using reference_constraint_t = typename Domain::reference_constraint_t;
  using variable_or_constant_t = typename Domain::variable_or_constant_t;
  using variable_t = typename Domain::variable_t;
  using variable_vector_t = typename Domain::variable_vector_t;
  using varname_t = typename Domain::varname_t;

  // typedefs for GhostDomain
  using ghost_linear_constraint_system_t =
      typename GhostDomain::linear_constraint_system_t;
  using ghost_linear_constraint_t = typename GhostDomain::linear_constraint_t;
  using ghost_linear_expression_t = typename GhostDomain::linear_expression_t;
  using ghost_reference_constraint_t =
      typename GhostDomain::reference_constraint_t;
  using ghost_variable_or_constant_t =
      typename GhostDomain::variable_or_constant_t;
  using ghost_variable_t = typename GhostDomain::variable_t;
  using ghost_variable_vector_t = typename GhostDomain::variable_vector_t;
  using ghost_varname_t = typename GhostDomain::varname_t;
  
public:
  using get_type_fn = std::function<variable_type(const variable_t &)>;
  // Map from variable to its ghost variables
  using ghost_variables_t =
      ghost_variables<Domain, GhostDomain, GhostDomainVarAlloc>;
  using ghost_variable_kind = typename ghost_variables_t::ghost_variable_kind;

private:
  using ghost_map_t = std::unordered_map<variable_t, ghost_variables_t>;
  // Reverse map used only for pretty printing: a CrabIR variable has
  // associated a set of ghost variables. This reverse map maps back
  // each ghost variable to its CrabIR variable.
  using ghost_var_id = unsigned;
  using rev_ghost_map_t =
      std::unordered_map<ghost_variable_t, std::pair<variable_t, ghost_var_id>>;

  // map Domain's variable to GhostDomain's variables
  ghost_map_t m_map;
  // reverse map
  rev_ghost_map_t m_rev_map;
  // create ghost variables in the ghost domain
  GhostDomainVarAlloc m_alloc;
  // function to get the type of each variable
  get_type_fn m_get_type;

  ghost_variables_t create_gvars(GhostDomainVarAlloc &alloc,
                                 const variable_t &v) const {
    auto alloc_fn = [](GhostDomainVarAlloc &alloc_, const varname_t,
                       const std::string &) { return alloc_.next(); };
    return ghost_variables_t::create(alloc, alloc_fn, v.name(), m_get_type(v));
  }

  ghost_variables_t create_gvars(GhostDomainVarAlloc &alloc,
                                 const ghost_variables_t &gvars) const {
    auto alloc_fn = [](GhostDomainVarAlloc &alloc_, const ghost_varname_t,
                       const std::string &) { return alloc_.next(); };
    return ghost_variables_t::create(alloc, alloc_fn, gvars);
  }

  ghost_variables_t insert_ghost_vars(const variable_t &v) {
    auto gvars = create_gvars(m_alloc, v);
    auto res = m_map.insert({v, gvars});
    if (res.second) {
      gvars.update_rev_varmap(m_rev_map, v);
    }
    return res.first->second;
  }

  ghost_variable_manager_with_variable_naming(ghost_map_t &&map,
                                              rev_ghost_map_t &&rev_map,
                                              GhostDomainVarAlloc &&alloc,
                                              get_type_fn get_type)
      : m_map(std::move(map)), m_rev_map(std::move(rev_map)),
        m_alloc(std::move(alloc)), m_get_type(get_type) {}

public:
  ghost_variable_manager_with_variable_naming(get_type_fn get_type)
      : m_get_type(get_type) {}
  ghost_variable_manager_with_variable_naming(const ghost_var_manager_t &o) =
      default;
  ghost_variable_manager_with_variable_naming(ghost_var_manager_t &&o) =
      default;
  ghost_var_manager_t &operator=(const ghost_var_manager_t &o) = default;
  ghost_var_manager_t &operator=(ghost_var_manager_t &&o) = default;

  ghost_variables_t get_or_insert(const variable_t &v) {
    auto it = m_map.find(v);
    if (it != m_map.end()) {
      return it->second;
    } else {
      return insert_ghost_vars(v);
    }
  }

  boost::optional<ghost_variables_t> get(const variable_t &v) const {
    auto it = m_map.find(v);
    if (it != m_map.end()) {
      return it->second;
    } else {
      return boost::none;
    }
  }

  static bool need_renaming() { return true; }

  ghost_var_manager_t join(const ghost_var_manager_t &right_man,
                           GhostDomain &left_val,
                           GhostDomain &right_val) const {

    assert(need_renaming());

    GhostDomainVarAlloc out_alloc(m_alloc, right_man.m_alloc);
    ghost_map_t out_map;
    rev_ghost_map_t out_rev_map;
    ghost_variable_vector_t left_vars, right_vars, out_vars;

    // upper bound to avoid reallocations
    size_t num_renamings = m_map.size();
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    // perform the common renamings
    for (auto const &kv : m_map) {
      const variable_t &v = kv.first;
      auto it = right_man.m_map.find(v);
      if (it != right_man.m_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          // crab::CrabStats::count(
          // domain_name() + ".count.join.skipped.inconsistent_dynamic_type");
          continue;
        }
        ghost_variables_t out_gvars(create_gvars(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_map, v);
      }
    }

    left_val.project(left_vars);
    left_val.rename(left_vars, out_vars);

    right_val.project(right_vars);
    right_val.rename(right_vars, out_vars);

    return ghost_var_manager_t(std::move(out_map), std::move(out_rev_map),
                               std::move(out_alloc), m_get_type);
  }

  ghost_var_manager_t meet(const ghost_var_manager_t &right_man,
                           GhostDomain &left_val,
                           GhostDomain &right_val) const {

    assert(need_renaming());

    GhostDomainVarAlloc out_alloc(m_alloc, right_man.m_alloc);
    ghost_map_t out_map;
    rev_ghost_map_t out_rev_map;

    ghost_variable_vector_t left_vars, right_vars, out_vars;
    ghost_variable_vector_t only_left_vars, only_left_out_vars;
    ghost_variable_vector_t only_right_vars, only_right_out_vars;
    // upper bound to avoid reallocations
    size_t left_renamings = m_map.size();
    size_t right_renamings = right_man.m_map.size();
    size_t num_renamings = left_renamings + right_renamings;
    left_vars.reserve(num_renamings);
    right_vars.reserve(num_renamings);
    out_vars.reserve(num_renamings);
    only_left_vars.reserve(left_renamings);
    only_left_out_vars.reserve(left_renamings);
    only_right_vars.reserve(right_renamings);
    only_right_out_vars.reserve(right_renamings);

    for (auto const &kv : m_map) {
      const variable_t &v = kv.first;
      auto it = right_man.m_map.find(v);
      if (it != right_man.m_map.end()) {
        // if we ignore unknown regions, types are the same by construction
        if (!crab_domain_params_man::get().region_skip_unknown_regions() &&
            !kv.second.same_type(it->second)) {
          // crab::CrabStats::count(
          //    domain_name() +
          //    ".count.meet_or_narrowing.skipped.inconsistent_dynamic_type");
          continue;
        }
        // common renaming
        ghost_variables_t out_gvars(create_gvars(out_alloc, kv.second));
        kv.second.add(left_vars);
        it->second.add(right_vars);
        out_gvars.add(out_vars);
        out_map.insert({v, out_gvars});
        out_gvars.update_rev_varmap(out_rev_map, v);
      } else {
        // only on left
        ghost_variables_t out_gvars(create_gvars(out_alloc, kv.second));
        kv.second.add(only_left_vars);
        out_gvars.add(only_left_out_vars);
        out_map.insert({kv.first, out_gvars});
        out_gvars.update_rev_varmap(out_rev_map, kv.first);
      }
    }

    // add missing maps from right operand
    for (auto const &kv : right_man.m_map) {
      const variable_t &v = kv.first;
      auto it = m_map.find(v);
      if (it == m_map.end()) {
        // only on right
        ghost_variables_t out_gvars(create_gvars(out_alloc, kv.second));
        kv.second.add(only_right_vars);
        out_gvars.add(only_right_out_vars);
        out_map.insert({kv.first, out_gvars});
        out_gvars.update_rev_varmap(out_rev_map, kv.first);
      }
    }

    // append only_left_vars and left_vars
    only_left_vars.insert(only_left_vars.end(), left_vars.begin(),
                          left_vars.end());
    // append only_left_out_vars end out_vars
    only_left_out_vars.insert(only_left_out_vars.end(), out_vars.begin(),
                              out_vars.end());

    left_val.project(only_left_vars);
    left_val.rename(only_left_vars, only_left_out_vars);

    // append only_right_vars and right_vars
    only_right_vars.insert(only_right_vars.end(), right_vars.begin(),
                           right_vars.end());
    // append only_right_out_vars end out_vars
    only_right_out_vars.insert(only_right_out_vars.end(), out_vars.begin(),
                               out_vars.end());

    right_val.project(only_right_vars);
    right_val.rename(only_right_vars, only_right_out_vars);

    return ghost_var_manager_t(std::move(out_map), std::move(out_rev_map),
                               std::move(out_alloc), m_get_type);
  }

  // left_val = meet(left_val, project(right_val, variables))
  void merge(const variable_vector_t &variables,
	     const ghost_var_manager_t &right_man,
	     GhostDomain &left_val, const GhostDomain &right_val) {

    if (variables.empty()) {
      return;
    }

    /// All of this is needed to avoid variable clashing with the left
    /// operand
    GhostDomainVarAlloc right_alloc(m_alloc, right_man.m_alloc);

    ghost_variable_vector_t new_left_vars, new_right_vars, old_right_vars;
    new_left_vars.reserve(variables.size());
    new_right_vars.reserve(variables.size());
    for (const variable_t &v : variables) {
      if (boost::optional<ghost_variables_t> old_right_gvars =
              right_man.get(v)) {
        ghost_variables_t new_left_gvars = get_or_insert(v);
        ghost_variables_t new_right_gvars = create_gvars(right_alloc, v);
        new_left_gvars.add(new_left_vars);
        new_right_gvars.add(new_right_vars);
        (*old_right_gvars).add(old_right_vars);
      }
    }

    assert(old_right_vars.size() == new_right_vars.size());
    assert(new_left_vars.size() == new_right_vars.size());

    /// Propagate invariants from the right to the left */
    GhostDomain copy(right_val);
    copy.project(old_right_vars);

    // Renaming in two steps to avoid variable clashing with the left
    // operand. The assumption is that old_right_vars and
    // new_left_vars might have common names. (it might be not
    // necessary anymore).
    copy.rename(old_right_vars, new_right_vars);
    copy.rename(new_right_vars, new_left_vars);
    left_val = left_val & copy;
  }

  void forget(const variable_t &var, GhostDomain &val) {
    auto it = m_map.find(var);
    if (it != m_map.end()) {
      it->second.forget(val);
      it->second.remove_rev_varmap(m_rev_map);
      m_map.erase(it);
    }
  }

  void project(const variable_vector_t &vars, GhostDomain &val) {
    ghost_variable_vector_t gvars_vec;
    for (auto const &v : vars) {
      if (boost::optional<ghost_variables_t> gvars = get(v)) {
        (*gvars).add(gvars_vec);
      }
    }

    val.project(gvars_vec);

    variable_vector_t sorted_vars(vars.begin(), vars.end());
    std::sort(sorted_vars.begin(), sorted_vars.end());
    for (auto it = m_map.begin(); it != m_map.end();) {
      if (!std::binary_search(sorted_vars.begin(), sorted_vars.end(),
                              it->first)) {
        it->second.remove_rev_varmap(m_rev_map);
        it = m_map.erase(it);
      } else {
        ++it;
      }
    }
  }

  ghost_variables_t dup(const ghost_variables_t &gvars) {
    return create_gvars(m_alloc, gvars);
  }

  void rename(const variable_vector_t &old_vars,
              const variable_vector_t &new_vars, GhostDomain &val) {

    ghost_variable_vector_t old_ghost_vars, new_ghost_vars;
    for (unsigned i = 0, sz = old_vars.size(); i < sz; ++i) {
      variable_t old_var = old_vars[i];
      auto it = m_map.find(old_var);
      if (it == m_map.end()) {
        continue;
      }
      // update m_var_map and m_rev_var_map to remove old_var and its
      // ghost variables
      const ghost_variables_t &old_gvars = it->second;
      old_gvars.add(old_ghost_vars);
      old_gvars.remove_rev_varmap(m_rev_map);
      m_map.erase(old_var); // add this point don't use old_gvars

      // update m_var_map and m_rev_var_map to accomodate new_var and its
      // ghost variables.
      variable_t new_var = new_vars[i];

      ghost_variables_t new_gvars = get_or_insert(new_var);
      new_gvars.add(new_ghost_vars);
    }
    val.rename(old_ghost_vars, new_ghost_vars);
  }

  // The ghost domain doesn't know about references so we convert a
  // reference constraint into a linear constraint.
  // A reference variable can be ghosted by multiple ghost variables.
  // The parameter @kind tells which one.
  ghost_linear_constraint_t
  ghosting_ref_cst_to_linear_cst(const reference_constraint_t &ref_cst,
				 ghost_variable_kind kind) {
    auto get_ghost_var = [this](const variable_t &var, ghost_variable_kind kind) ->
      boost::optional<ghost_variable_t> {
      auto ghost_vars  = get_or_insert(var);
      if (kind == ghost_variable_kind::ADDRESS) {
	return ghost_vars.get_var();
      } else if (kind == ghost_variable_kind::OFFSET) {
	if (ghost_vars.has_offset_and_size()) {
	  return ghost_vars.get_offset_and_size().get_offset();
	}
      } else if (kind == ghost_variable_kind::SIZE) {
	if (ghost_vars.has_offset_and_size()) {
	  return ghost_vars.get_offset_and_size().get_size();
	}
      }
      return boost::none;
    };
      
    if (ref_cst.is_tautology()) {
      return ghost_linear_constraint_t::get_true();
    } else if (ref_cst.is_contradiction()) {
      return ghost_linear_constraint_t::get_false();
    } else {
      if (ref_cst.is_unary()) {
        assert(ref_cst.lhs().get_type().is_reference());
        ghost_variable_t x = get_or_insert(ref_cst.lhs()).get_var();
        if (ref_cst.is_equality()) {
          return ghost_linear_constraint_t(x == number_t(0));
        } else if (ref_cst.is_disequality()) {
          return ghost_linear_constraint_t(x != number_t(0));
        } else if (ref_cst.is_less_or_equal_than()) {
          return ghost_linear_constraint_t(x <= number_t(0));
        } else if (ref_cst.is_less_than()) {
          return ghost_linear_constraint_t(x < number_t(0));
        } else if (ref_cst.is_greater_or_equal_than()) {
          return ghost_linear_constraint_t(x >= number_t(0));
        } else if (ref_cst.is_greater_than()) {
          return ghost_linear_constraint_t(x > number_t(0));
        }
      } else {
        assert(ref_cst.lhs().get_type().is_reference());
        assert(ref_cst.rhs().get_type().is_reference());
	boost::optional<ghost_variable_t> x = get_ghost_var(ref_cst.lhs(), kind);
	boost::optional<ghost_variable_t> y = get_ghost_var(ref_cst.rhs(), kind);
	if (!x || !y) {
	  return ghost_linear_constraint_t::get_true();
	}	
        number_t offset = ref_cst.offset();
        if (ref_cst.is_equality()) {
          return ghost_linear_constraint_t(*x == *y + offset);
        } else if (ref_cst.is_disequality()) {
          return ghost_linear_constraint_t(*x != *y + offset);
        } else if (ref_cst.is_less_or_equal_than()) {
          return ghost_linear_constraint_t(*x <= *y + offset);
        } else if (ref_cst.is_less_than()) {
          return ghost_linear_constraint_t(*x < *y + offset);
        } else if (ref_cst.is_greater_or_equal_than()) {
          return ghost_linear_constraint_t(*x >= *y + offset);
        } else if (ref_cst.is_greater_than()) {
          return ghost_linear_constraint_t(*x > *y + offset);
        }
      }
    }
    CRAB_ERROR("unexpected reference constraint");
  }

  void write(crab_os &os, const GhostDomain &val) const {
    std::unordered_map<std::string, std::string> renaming_map;
    // Assigns a string name to each ghost variable
    ghost_variables_t::mk_renaming_map(m_rev_map, m_get_type, renaming_map);
    m_alloc.add_renaming_map(renaming_map);
    os << val;
    m_alloc.clear_renaming_map();
  }

  /*
   * These functions must be provided because Domain and GhostDomain
   * can have different variable factories with different varname
   * types
   */

  ghost_variable_or_constant_t
  rename_variable_or_constant(const variable_or_constant_t &v) {
    if (v.is_constant()) {
      return ghost_variable_or_constant_t(v.get_constant(), v.get_type());
    } else {
      ghost_variable_t bv = get_or_insert(v.get_variable()).get_var();
      return ghost_variable_or_constant_t(bv);
    }
  }

  ghost_linear_expression_t rename_linear_expr(const linear_expression_t &e) {
    ghost_linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      ghost_variable_t nv = get_or_insert(v).get_var();
      out = out + (coef * nv);
    }
    return out;
  }

  ghost_linear_constraint_t rename_linear_cst(const linear_constraint_t &cst) {
     return ghost_linear_constraint_t(
          rename_linear_expr(cst.expression()),
          (typename ghost_linear_constraint_t::kind_t)cst.kind());
  }

  ghost_linear_constraint_system_t
  rename_linear_cst_sys(const linear_constraint_system_t &csts) {
    ghost_linear_constraint_system_t out;
    for (auto const &cst : csts) {
      out += rename_linear_cst(cst);
    }
    return out;
  }

  boost::optional<variable_t> rev_rename_var(const ghost_variable_t &v,
                                             bool ignore_references) const {
    auto it = m_rev_map.find(v);
    if (it != m_rev_map.end()) {
      variable_t v = it->second.first;
      ghost_var_id ghost_id = it->second.second;
      if (!ignore_references || !v.get_type().is_reference()) {
        if (ghost_id == 1) { // we only get the first ghost variable
                             // which represents the address of the
                             // reference. Skipping other ghost
                             // variables such as offsets and sizes.
          return v;
        }
      }
    }
    return boost::none;
  }

  boost::optional<linear_expression_t>
  rev_rename_linear_expr(const ghost_linear_expression_t &e,
                         bool ignore_references) const {
    linear_expression_t out(e.constant());
    for (auto it = e.begin(), et = e.end(); it != et; ++it) {
      const ghost_variable_t &v = (*it).second;
      const number_t &coef = (*it).first;
      if (boost::optional<variable_t> v_opt =
              rev_rename_var(v, ignore_references)) {
        out = out + (coef * (*v_opt));
      } else {
        return boost::optional<linear_expression_t>();
      }
    }
    return out;
  }

  boost::optional<linear_constraint_t>
  rev_rename_linear_cst(const ghost_linear_constraint_t &cst,
                        bool ignore_references) const {
    if (boost::optional<linear_expression_t> e =
            rev_rename_linear_expr(cst.expression(), ignore_references)) {

      return linear_constraint_t(
            *e, (typename linear_constraint_t::kind_t)cst.kind());
    } else {
      return boost::optional<linear_constraint_t>();
    }
  }
};

} // end namespace region_domain_impl
} // end namespace domains
} // end namespace crab
