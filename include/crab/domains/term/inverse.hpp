#pragma once

/* To propagate down, our domain needs to provide
 * operations for computing predecessor states of operations. */

namespace crab {
namespace domains {
namespace term {
// Default implementations - assign becomes a lincons,
// and everything else is a nop.
template <class Num, class Var, class Dom> class InverseOps {
public:
  typedef Num num_t;
  typedef Var varname_t;

  typedef Dom dom_t;
  typedef typename dom_t::variable_t var_t;

  typedef typename dom_t::linear_constraint_t lincst_t;
  typedef typename dom_t::linear_expression_t linexp_t;

  // x = sum c_i y_i -> -x + sum c_i y_i = 0.
  static void assign(dom_t &dom, varname_t x, linexp_t expr) {
    linexp_t x_expr = expr - var_t(x);
    dom += lincst_t(x_expr, lincst_t::EQUALITY);
  }

  // x = y op z
  static void apply(dom_t &dom, operation_t op, varname_t x, varname_t y,
                    varname_t z) {
    switch (op) {
    case OP_ADDITION:
      dom += lincst_t(var_t(y) + var_t(z) - var_t(x), lincst_t::EQUALITY);
      break;
    case OP_SUBTRACTION:
      dom += lincst_t(var_t(y) - var_t(z) - var_t(x), lincst_t::EQUALITY);
      break;
    case OP_MULTIPLICATION:
    case OP_SDIV:
    case OP_UDIV:
    case OP_SREM:
    case OP_UREM:
      break;
    }
    return;
  }

  // x = y op k
  static void apply(dom_t &dom, operation_t op, varname_t x, varname_t y,
                    num_t k) {
    switch (op) {
    case OP_ADDITION:
      dom += lincst_t(var_t(y) + k - var_t(x), lincst_t::EQUALITY);
      break;
    case OP_SUBTRACTION:
      dom += lincst_t(var_t(y) - k - var_t(x), lincst_t::EQUALITY);
      break;
    case OP_MULTIPLICATION:
      dom += lincst_t(var_t(y) * k - var_t(x), lincst_t::EQUALITY);
      break;
    case OP_SDIV:
      // dom += lincst_t(var_t(y) - var_t(x)*k, INEQUALITY);
      // What are the precise semantics of division?
      break;
    case OP_UDIV:
    case OP_SREM:
    case OP_UREM:
      break;
    }
    return;
  }
};

template <class Num, class Var>
class InverseOps<Num, Var, interval_domain<Num, Var>> {
public:
  typedef Num num_t;
  typedef Var varname_t;

  typedef interval_domain<num_t, varname_t> dom_t;
  typedef typename dom_t::variable_t var_t;

  typedef typename dom_t::linear_constraint_t lincst_t;
  typedef typename dom_t::linear_expression_t linexp_t;

  typedef ikos::bound<num_t> bound_t;
  typedef ikos::interval<num_t> interval_t;

  // x = sum c_i y_i -> -x + sum c_i y_i = 0.
  static void assign(dom_t &dom, varname_t x, linexp_t expr) {
    linexp_t x_expr = expr - var_t(x);
    dom += lincst_t(x_expr, lincst_t::EQUALITY);
  }

  static void apply(dom_t &dom, operation_t op, varname_t x, varname_t y,
                    varname_t z) {
    switch (op) {
    case OP_ADDITION:
      dom += lincst_t(var_t(y) + var_t(z) - var_t(x), lincst_t::EQUALITY);
      break;
    case OP_SUBTRACTION:
      dom += lincst_t(var_t(y) - var_t(z) - var_t(x), lincst_t::EQUALITY);
      break;
      // FIXME: check this reflects the semantics correctly.
    case OP_MULTIPLICATION: {
      interval_t x0 = dom[x];
      interval_t y0 = dom[y];
      interval_t z0 = dom[z];

      // For x = yz, if 0 in x and 0 in y, we can't infer anything about z.
      if (!x0[0] || !z0[0])
        dom.set(y, y0 & (x0 / z0));
      if (!x0[0] || !y0[0])
        dom.set(z, z0 & (x0 / y0));
      break;
    }
    // This assumes we're always rounding towards -inf.
    // Will have to fix if semantics round towards 0.
    case OP_SDIV: {
      // Adjust for rounding.
      interval_t x0_relax = dom[x] + interval_t(num_t(0), num_t(1));
      interval_t y0 = dom[y];
      interval_t z0 = dom[z];

      dom.set(y, y0 & (x0_relax * z0));
      dom.set(z, z0 & (y0 / x0_relax));
      break;
    }
    case OP_UDIV:
    case OP_SREM:
    case OP_UREM:
      break;
    }
    return;
  }

  // x = y op k
  static void apply(dom_t &dom, operation_t op, varname_t x, varname_t y,
                    num_t k) {
    switch (op) {
    case OP_ADDITION:
      dom.set(y, dom[y] & (dom[x] - k));
      break;
    case OP_SUBTRACTION:
      dom.set(y, dom[y] & (dom[x] + k));
      break;
    case OP_MULTIPLICATION:
      if (k != 0)
        dom.set(y, dom[y] & (dom[x] / k));
      break;
    case OP_SDIV:
      // GKG: Again, potentially correct for round-towards-0.
      dom.set(y, dom[y] & ((dom[x] + interval_t(num_t(0, 1))) * k));
      break;
    case OP_UDIV:
    case OP_SREM:
    case OP_UREM:
      break;
    }
    return;
  }
};

} // end namespace term
} // end namespace domains
} // end namespace crab
