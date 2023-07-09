#pragma once

#include <crab/domains/abstract_domain_params.hpp>
#include <crab/support/debug.hpp>

namespace crab {
namespace domains {
namespace region_domain_impl {

// The region domain models three things about references: (1) the
// address of the reference, (2) its offset (distance to the base
// address of the memory object to which the reference points to), and
// (3) the size of the memory object. The ghost variables (2)-(3) are
// optional.
//
// To make the code of the region domain simpler and extensible, we
// group together all the ghost variables in a class and provide
// abstract operations so that the rest of the code is mostly
// unaware of how many ghost variables the domain is keeping track
// and how they should be updated.
template <class Domain, class GhostDomain, typename GhostVarnameAlloc>
class ghost_variables {
public:
  using domain_t = Domain;
  using ghost_domain_t = GhostDomain;
  using ghost_varname_allocator_t = GhostVarnameAlloc;

  using linear_constraint_t = typename domain_t::linear_constraint_t;
  using linear_expression_t = typename domain_t::linear_expression_t;
  using reference_constraint_t = typename domain_t::reference_constraint_t;
  using variable_or_constant_t = typename domain_t::variable_or_constant_t;
  using variable_or_constant_vector_t =
      typename domain_t::variable_or_constant_vector_t;
  using variable_t = typename domain_t::variable_t;
  using variable_vector_t = typename domain_t::variable_vector_t;
  using number_t = typename domain_t::number_t;
  using varname_t = typename domain_t::varname_t;

  using ghost_variable_vector_t = typename ghost_domain_t::variable_vector_t;
  using ghost_variable_t = typename ghost_domain_t::variable_t;
  using ghost_varname_t = typename ghost_variable_t::varname_t;
  using ghost_variable_or_constant_t =
      typename ghost_domain_t::variable_or_constant_t;
  using ghost_linear_expression_t =
      typename ghost_domain_t::linear_expression_t;
  using ghost_linear_constraint_t =
      typename ghost_domain_t::linear_constraint_t;
  using ghost_linear_constraint_system_t =
      typename ghost_domain_t::linear_constraint_system_t;

  using var_allocator_fn_t = std::function<ghost_varname_t(
      ghost_varname_allocator_t &alloc, const varname_t &v,
      const std::string &role)>;

  using ghost_var_allocator_fn_t = std::function<ghost_varname_t(
      ghost_varname_allocator_t &alloc, const ghost_varname_t &v,
      const std::string &role)>;

  // class for shadowing offset and sizes
  class ghost_offset_and_size {
    friend class ghost_variables;

    ghost_variable_t m_offset;
    ghost_variable_t m_size;

    ghost_offset_and_size(ghost_variable_t &&offset, ghost_variable_t &&size)
        : m_offset(std::move(offset)), m_size(std::move(size)) {
      if (!m_offset.get_type().is_integer() ||
          !m_size.get_type().is_integer()) {
        CRAB_ERROR("ghost offset and size variables must be integers");
      }
    }
    ghost_offset_and_size(ghost_variable_t offset, ghost_variable_t size)
        : m_offset(offset), m_size(size) {
      if (!m_offset.get_type().is_integer() ||
          !m_size.get_type().is_integer()) {
        CRAB_ERROR("ghost offset and size variables must be integers");
      }
    }

  public:
    
    ghost_variable_t get_offset() const { return m_offset; }

    ghost_variable_t get_size() const { return m_size; }

    void init(ghost_domain_t &dom,
              const ghost_variable_or_constant_t &size) const {
      dom.assign(m_offset, number_t(0));
      if (size.is_constant()) {
        dom.assign(m_size, size.get_constant());
        dom.intrinsic("var_packing_merge", {m_offset, m_size}, {});
      } else {
        dom.assign(m_size, size.get_variable());
        dom.intrinsic("var_packing_merge",
                      {m_offset, m_size, size.get_variable()}, {});
      }
    }

    void assign(ghost_domain_t &dom, const ghost_offset_and_size &rhs) const {
      if (m_offset.get_type() != rhs.m_offset.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }
      if (m_size.get_type() != rhs.m_size.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }

      dom.assign(m_offset, rhs.m_offset);
      dom.assign(m_size, rhs.m_size);
    }

    void weak_assign(ghost_domain_t &dom, const ghost_offset_and_size &rhs) const {
      if (m_offset.get_type() != rhs.m_offset.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }
      if (m_size.get_type() != rhs.m_size.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }

      dom.weak_assign(m_offset, rhs.m_offset);
      dom.weak_assign(m_size, rhs.m_size);
    }
    
    void assign(ghost_domain_t &dom, const ghost_offset_and_size &rhs,
                const ghost_linear_expression_t &offset) const {
      if (m_offset.get_type() != rhs.m_offset.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }
      if (m_size.get_type() != rhs.m_size.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }

      dom.assign(m_offset, rhs.m_offset + offset);
      dom.assign(m_size, rhs.m_size);
    }

    void weak_assign(ghost_domain_t &dom, const ghost_offset_and_size &rhs,
                const ghost_linear_expression_t &offset) const {
      if (m_offset.get_type() != rhs.m_offset.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }
      if (m_size.get_type() != rhs.m_size.get_type()) {
        CRAB_ERROR("Type inconsistency while assign ghost object_offset_size "
                   "variables");
      }

      dom.weak_assign(m_offset, rhs.m_offset + offset);
      dom.weak_assign(m_size, rhs.m_size);
    }

    // return true iff dom |= offset + byte_sz <= size
    bool is_deref(ghost_domain_t &dom,
		  // number of dereferenceable bytes
                  const ghost_variable_or_constant_t &byte_sz) const {
      // TODO: we don't check that m_offset >= 0
      ghost_domain_t tmp(dom);
      if (byte_sz.is_constant()) {
        tmp += ghost_linear_constraint_t(m_offset + byte_sz.get_constant() >
                                         m_size);
      } else {
        tmp += ghost_linear_constraint_t(m_offset + byte_sz.get_variable() >
                                         m_size);
      }
      return tmp.is_bottom();
    }

    void expand(ghost_domain_t &dom, const ghost_offset_and_size &rhs) const {
      if (m_offset.get_type() != rhs.m_offset.get_type()) {
        CRAB_ERROR("Type inconsistency while expand ghost object_offset_size "
                   "variables");
      }
      if (m_size.get_type() != rhs.m_size.get_type()) {
        CRAB_ERROR("Type inconsistency while expand ghost object_offset_size "
                   "variables");
      }

      dom.expand(m_offset, rhs.m_offset);
      dom.expand(m_size, rhs.m_size);
    }

    void forget(ghost_domain_t &dom) const {
      dom -= m_offset;
      dom -= m_size;
    }
  }; /* end class ghost_offset_and_size */

  // - if the variable is a reference then `m_var` is its symbolic
  //   address and its type is an integer.  note that the base domain
  //   is completely unaware of regions and references.
  // - if the variable is T=integer|real|bool or arrays of T then
  //   `m_var` is its value and the type is preserved.
  // - if the variable is a region of type T=integer|real|bool then
  //   `m_var` represents its contents and its type is T.
  // - if the variable is a region of type reference then `m_var`
  //   represents its content and its type is an integer.
  ghost_variable_t m_var;
  // if the variable is a reference or a region of references then
  // m_object_offset_size.first is the `offset`, the distance between its
  // address and the base address of the memory object to which
  // address belongs to, and m_object_offset_size.second is the `size` of the
  // memory object.
  boost::optional<ghost_offset_and_size> m_object_offset_size;
  // The CrabIR type of the variable which is being shadowed
  variable_type m_vty;

  void check_types() {
    if (m_var.get_type().is_region() || m_var.get_type().is_reference()) {
      CRAB_ERROR("Ghost variables cannot be regions or references");
    }

    if (m_object_offset_size) {
      if (!(*m_object_offset_size).get_offset().get_type().is_integer()) {
        CRAB_ERROR("ghost offset can only be integer");
      }
      if (!(*m_object_offset_size).get_size().get_type().is_integer()) {
        CRAB_ERROR("ghost size can only be integer");
      }
    }
  }

  ghost_variables(ghost_variable_t bv, variable_type vty)
      : m_var(bv), m_object_offset_size(boost::none), m_vty(vty) {
    check_types();
  }

  ghost_variables(ghost_variable_t bv, ghost_variable_t offset,
                  ghost_variable_t size, variable_type vty)
      : m_var(bv), m_object_offset_size(ghost_offset_and_size(offset, size)),
        m_vty(vty) {
    check_types();
  }

  // Create a ghost variable in the ghost domain for ghosting v
  static ghost_variable_t make_ghost_variable(ghost_varname_t name,
                                              const variable_type &vty) {
    if (vty.is_reference()) {
      return ghost_variable_t(name, crab::INT_TYPE,
                              32 /*should be defined in Params*/);
    } else if (vty.is_region()) {
      variable_type_kind ty;
      unsigned bitwidth = 0;
      if (vty.is_bool_region()) {
        ty = BOOL_TYPE;
        bitwidth = 1;
      } else if (vty.is_integer_region()) {
        ty = INT_TYPE;
        bitwidth = vty.get_integer_region_bitwidth();
      } else if (vty.is_reference_region()) {
        ty = INT_TYPE;
        bitwidth = 32;
      } else if (vty.is_real_region()) {
        ty = REAL_TYPE;
      } else if (vty.is_bool_array_region()) {
        ty = ARR_BOOL_TYPE;
      } else if (vty.is_int_array_region()) {
        ty = ARR_INT_TYPE;
      } else if (vty.is_real_array_region()) {
        ty = ARR_REAL_TYPE;
      } else {
        assert(false);
        CRAB_ERROR("make_ghost_variable: unreachable");
      }
      return ghost_variable_t(name, ty, bitwidth);
    } else {
      return ghost_variable_t(name, vty);
    }
  }

  // VarAllocFn can be either var_allocator_fn_t or ghost_var_allocator_fn_t
  // VariableName can be either varname_t or ghost_varname_t
  template <typename VarAllocFn, typename VariableName>
  static ghost_variable_t
  make_ghost_variable(ghost_varname_allocator_t &alloc, VarAllocFn alloc_fn,
                      const VariableName &name, const std::string &role,
                      const variable_type &vty) {
    ghost_varname_t gname(alloc_fn(alloc, name, role));
    return make_ghost_variable(gname, vty);
  }

public:

  enum class ghost_variable_kind {
    ADDRESS = 0, OFFSET = 1, SIZE = 2
  };
  
  static ghost_variables create(ghost_varname_allocator_t &alloc,
                                var_allocator_fn_t alloc_fn,
                                const varname_t &name, const variable_type &vty,
                                unsigned line = 0) {
    if (vty.is_unknown_region()) {
      CRAB_ERROR("create should not be called at line ", line, " with a type ",
                 vty);
    }

    if (crab_domain_params_man::get().region_is_dereferenceable() &&
        (vty.is_reference() || vty.is_reference_region())) {
      return ghost_variables(
          make_ghost_variable(alloc, alloc_fn, name, ".address", vty),
          make_ghost_variable(alloc, alloc_fn, name, ".offset", vty),
          make_ghost_variable(alloc, alloc_fn, name, ".size", vty), vty);
    } else {
      return ghost_variables(
          make_ghost_variable(alloc, alloc_fn, name, "", vty), vty);
    }
  }

  // Create a fresh set of ghost variables from x. We only care about
  // the types of gvars.
  static ghost_variables create(ghost_varname_allocator_t &alloc,
                                ghost_var_allocator_fn_t alloc_fn,
                                const ghost_variables &gvars) {
    if (gvars.has_offset_and_size()) {
      ghost_offset_and_size offset_size = gvars.get_offset_and_size();
      auto addr_gvar = gvars.get_var();
      auto off_gvar = offset_size.get_offset();
      auto sz_gvar = offset_size.get_size();
      return ghost_variables(
          make_ghost_variable(alloc, alloc_fn, addr_gvar.name(), ".address",
                              addr_gvar.get_type()),
          make_ghost_variable(alloc, alloc_fn, off_gvar.name(), ".offset",
                              off_gvar.get_type()),
          make_ghost_variable(alloc, alloc_fn, sz_gvar.name(), ".size",
                              sz_gvar.get_type()),
          gvars.m_vty);
    } else {
      auto gvar = gvars.get_var();
      return ghost_variables(make_ghost_variable(alloc, alloc_fn, gvar.name(),
                                                 "", gvar.get_type()),
                             gvars.m_vty);
    }
  }

  // Return the variable factory used to create the ghost variables
  // (i.e., passed to create)
  ghost_varname_allocator_t &get_vfac() const {
    auto &vfac =
        const_cast<ghost_varname_t *>(&(m_var.name()))->get_var_factory();
    return vfac;
  }

  // Return true if the types of the ghost variables are
  // syntactically equal.
  bool same_type(const ghost_variables &o) const {
    // Check that the types of CrabIR variables are consistent
    if (!((m_vty.is_region() && o.m_vty.is_unknown_region()) ||
          (m_vty.is_unknown_region() && o.m_vty.is_region()) ||
          (m_vty == o.m_vty))) {
      return false;
    }

    // Check that the types of the ghost variables are syntactically
    // equal.
    if (get_var().get_type() != o.get_var().get_type()) {
      return false;
    }

    bool x_has_offset_and_size = has_offset_and_size();
    bool y_has_offset_and_size = o.has_offset_and_size();

    if (x_has_offset_and_size && y_has_offset_and_size) {
      ghost_offset_and_size x_offset_size = get_offset_and_size();
      ghost_offset_and_size y_offset_size = o.get_offset_and_size();
      return (x_offset_size.get_offset().get_type() ==
                  y_offset_size.get_offset().get_type() &&
              x_offset_size.get_size().get_type() ==
                  y_offset_size.get_size().get_type());
    }

    return (!x_has_offset_and_size && !y_has_offset_and_size);
  }

  // Store all ghost variables into a vector
  void add(ghost_variable_vector_t &vec) const {
    vec.push_back(m_var);
    if (m_object_offset_size) {
      vec.push_back((*m_object_offset_size).get_offset());
      vec.push_back((*m_object_offset_size).get_size());
    }
  }

  // Abstract domain operations
  void assign(ghost_domain_t &dom, const ghost_variables &rhs) const {
    if (m_var.get_type() != rhs.m_var.get_type()) {
      CRAB_ERROR("Type inconsistency while assign in ghost variables ",
                 m_var.get_type(), " and ", rhs.m_var.get_type());
    }

    if (m_var.get_type().is_bool()) {
      dom.assign_bool_var(m_var, rhs.m_var, false);
    } else {
      dom.assign(m_var, rhs.m_var);
    }

    if (m_object_offset_size && rhs.m_object_offset_size) {
      (*m_object_offset_size).assign(dom, *(rhs.m_object_offset_size));
    }
  }

  void expand(ghost_domain_t &dom, const ghost_variables &new_gvars) const {
    if (m_var.get_type() != new_gvars.m_var.get_type()) {
      CRAB_ERROR("Type inconsistency while expand in ghost variables ",
                 m_var.get_type(), " and ", new_gvars.m_var.get_type());
    }
    dom.expand(m_var, new_gvars.m_var);
    if (m_object_offset_size && new_gvars.m_object_offset_size) {
      (*m_object_offset_size).expand(dom, *(new_gvars.m_object_offset_size));
    }
  }

  void forget(ghost_domain_t &dom) const {
    dom -= m_var;
    if (m_object_offset_size) {
      (*m_object_offset_size).forget(dom);
    }
  }

  ghost_variable_t get_var() const { return m_var; }

  bool has_offset_and_size() const { return (bool)m_object_offset_size; }

  // Add ghost variables in the reverse map. A reverse map (rev_map)
  // maps ghost variables to CrabIR variables.
  void update_rev_varmap(
      std::unordered_map<ghost_variable_t, std::pair<variable_t, ghost_variable_kind>>
          &rev_map,
      variable_t v) const {
    rev_map.insert({m_var, {v, 1}});
    if (m_object_offset_size) {
      rev_map.insert({(*m_object_offset_size).get_offset(), {v, 2}});
      rev_map.insert({(*m_object_offset_size).get_size(), {v, 3}});
    }
  }

  // Remove ghost variables from the reverse map. A reverse map
  // (rev_map) maps ghost variables to CrabIR variables.
  void remove_rev_varmap(
      std::unordered_map<ghost_variable_t, std::pair<variable_t, ghost_variable_kind>>
          &rev_map) const {
    rev_map.erase(m_var);
    if (m_object_offset_size) {
      rev_map.erase((*m_object_offset_size).get_offset());
      rev_map.erase((*m_object_offset_size).get_size());
    }
  }

  // Assign a string name to each ghost variable
  static void mk_renaming_map(
      const std::unordered_map<ghost_variable_t,
                               std::pair<variable_t, ghost_variable_kind>> &rev_map,
      std::function<variable_type(const variable_t &)> get_type_fn,
      std::unordered_map<std::string, std::string> &out_str_map) {
    for (auto &kv : rev_map) {
      ghost_variable_t bv = kv.first;
      variable_t v = kv.second.first;
      ghost_variable_kind ghost_id = kv.second.second;

      variable_type vty = get_type_fn(v);
      if (crab_domain_params_man::get().region_is_dereferenceable() &&
          (vty.is_reference() || vty.is_reference_region())) {
        if (ghost_id == ghost_variable_kind::ADDRESS) {
          out_str_map[bv.name().str()] = v.name().str() + ".address";
        } else if (ghost_id == ghost_variable_kind::OFFSET) {
          out_str_map[bv.name().str()] = v.name().str() + ".offset";
        } else if (ghost_id == ghost_variable_kind::SIZE) {
          out_str_map[bv.name().str()] = v.name().str() + ".size";
        }
      } else {
        out_str_map[bv.name().str()] = v.name().str();
      }
    }
  }

  ghost_offset_and_size get_offset_and_size() const {
    if (!m_object_offset_size) {
      CRAB_ERROR("ghost variables for offset and size have not been allocated");
    }
    return *m_object_offset_size;
  }

  void write(crab_os &o) const {
    if (!m_object_offset_size) {
      o << m_var;
    } else {
      o << "{reference:" << m_var
        << ",offset:" << (*m_object_offset_size).get_offset()
        << ",size:" << (*m_object_offset_size).get_size() << "}";
    }
  }

  friend crab_os &operator<<(crab_os &o, const ghost_variables &gvars) {
    gvars.write(o);
    return o;
  }
}; /* end class ghost_variables */

} // end namespace region_domain_impl
} // end namespace domains
} // end namespace crab
