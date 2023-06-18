#pragma once

/**
 * Implement generic backward assignments.
 *
 * Be aware that unless the assignment is invertible the result is an
 * over-approximation so we need to adapt these operations in case we
 * need under-approximations.
 **/

#include <crab/config.h>

#include <crab/domains/abstract_domain_operators.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace domains {

template <class AbsDom> class BackwardAssignOps {
public:
  using number_t = typename AbsDom::number_t;
  using varname_t = typename AbsDom::varname_t;
  using variable_t = typename AbsDom::variable_t;
  using linear_constraint_t = typename AbsDom::linear_constraint_t;
  using linear_expression_t = typename AbsDom::linear_expression_t;

  /*
   * Backward x := e
   *
   *  General case:
   *   if x does not appear in e
   *      1) add constraint x = e
   *      2) forget x
   *   else
   *      1) add new variable x'
   *      2) add constraint x = e[x'/x]
   *      3) forget x
   *      4) rename x' as x
   *
   *  Invertible operation (y can be equal to x):
   *    x = y + k <--> y = x - k
   *    x = y - k <--> y = x + k
   *    x = y * k <--> y = x / k  if (k != 0)
   *    x = y / k <--> y = x * k  if (k != 0)
   *
   *  Fallback case:
   *   forget(x)
   **/

  // x := e
  static void assign(AbsDom &dom, const variable_t &x,
                     const linear_expression_t &e, const AbsDom &inv) {
    crab::CrabStats::count(dom.domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(dom.domain_name() + ".backward_assign");

    if (dom.is_bottom())
      return;

    if (std::find(e.variables_begin(), e.variables_end(), x) !=
        e.variables_end()) {
      auto &vfac = const_cast<varname_t *>(&(x.name()))->get_var_factory();
      variable_t old_x(vfac.get(), x.get_type());
      std::map<variable_t, variable_t> renaming_map;
      renaming_map.insert({x, old_x});
      linear_expression_t renamed_e = e.rename(renaming_map);
      dom += linear_constraint_t(renamed_e - x, linear_constraint_t::EQUALITY);
      dom -= x;
      dom.rename({old_x}, {x});
    } else {
      dom += linear_constraint_t(e - x, linear_constraint_t::EQUALITY);
      dom -= x;
    }
    dom = dom & inv;
  }

  // x := y op k
  static void apply(AbsDom &dom, arith_operation_t op, const variable_t &x,
                    const variable_t &y, number_t k, const AbsDom &inv) {
    crab::CrabStats::count(dom.domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(dom.domain_name() + ".backward_apply");

    if (dom.is_bottom()) {
      return;
    }

    CRAB_LOG("backward", crab::outs()
                             << x << ":=" << y << " " << op << " " << k << "\n"
                             << "BEFORE " << dom << "\n";);

    switch (op) {
    case OP_ADDITION:
      dom.apply(OP_SUBTRACTION, y, x, k);
      if (!(x == y)) {
        dom -= x;
      }
      break;
    case OP_SUBTRACTION:
      dom.apply(OP_ADDITION, y, x, k);
      if (!(x == y)) {
        dom -= x;
      }
      break;
    case OP_MULTIPLICATION:
      if (k != 0) {
        dom.apply(OP_SDIV, y, x, k);
        if (!(x == y)) {
          dom -= x;
        }
      } else {
        dom -= x;
      }
      break;
    case OP_SDIV:
      if (k != 0) {
        dom.apply(OP_MULTIPLICATION, y, x, k);
        if (!(x == y)) {
          dom -= x;
        }
      } else {
        dom -= x;
      }
      break;    
    default:
      //case OP_UDIV:
      //case OP_SREM:
      //case OP_UREM:      
      CRAB_WARN("backwards x:= y ", op, " k is not implemented");
      dom -= x;
    }
    dom = dom & inv;

    CRAB_LOG("backward", crab::outs() << "AFTER " << dom << "\n");
    return;
  }

  // x = y op z
  static void apply(AbsDom &dom, arith_operation_t op, const variable_t &x,
                    const variable_t &y, const variable_t &z,
                    const AbsDom &inv) {
    crab::CrabStats::count(dom.domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(dom.domain_name() + ".backward_apply");

    if (dom.is_bottom()) {
      return;
    }

    CRAB_LOG("backward", crab::outs()
                             << x << ":=" << y << " " << op << " " << z << "\n"
                             << "BEFORE " << dom << "\n";);

    switch (op) {
    case OP_ADDITION:
      assign(dom, x, linear_expression_t(y + z), inv);
      break;
    case OP_SUBTRACTION:
      assign(dom, x, linear_expression_t(y - z), inv);
      break;
    case OP_MULTIPLICATION:
    case OP_SDIV:
    case OP_UDIV:
    case OP_SREM:
    case OP_UREM:
      CRAB_WARN("backwards x = y ", op, " z not implemented");
      dom -= x;
      dom = dom & inv;
      break;
    }
    CRAB_LOG("backward", crab::outs() << "AFTER " << dom << "\n");
  }
};

} // end namespace domains
} // end namespace crab
