#pragma once

#include <crab/config.h>

#include <crab/domains/abstract_domain.hpp>
#include <crab/domains/abstract_domain_specialized_traits.hpp>
#include <crab/domains/intervals.hpp>
#include <crab/support/debug.hpp>
#include <crab/support/stats.hpp>

namespace crab {
namespace domains {
using elina_domain_id_t = enum { ELINA_ZONES, ELINA_OCT, ELINA_PK };

template <typename Number> class ElinaDefaultParams {
public:
  // use integers with truncation rounding
  enum { use_integers = 1 };
};

template <> class ElinaDefaultParams<ikos::q_number> {
public:
  // use reals
  enum { use_integers = 0 };
};
} // namespace domains
} // namespace crab

#ifndef HAVE_ELINA
/*
 * Dummy implementation if Elina not found
 */
#include <crab/domains/dummy_abstract_domain.hpp>
namespace crab {
namespace domains {
template <typename N, typename V, elina_domain_id_t Dom,
          class Params = ElinaDefaultParams<N>>
class elina_domain final
  : public dummy_abstract_domain<elina_domain<N, V, Dom, Params>> {
public:
  std::string not_implemented_msg() const override {
    return "No Elina. Run cmake with -DCRAB_USE_ELINA=ON";
  }
};
} // namespace domains
} // namespace crab
#else

/*
 * Real implementation starts here.
 */

#include <crab/domains/abstract_domain_params.hpp>
#include <crab/domains/elina/elina.hpp>
#include <crab/domains/linear_interval_solver.hpp>

#include <algorithm>
#include <boost/bimap.hpp>
#include <boost/optional.hpp>

/**
 * If template parameter Number is ikos::q_number then the Elina
 * domain will use real variables. Otherwise, all variables are
 * integers.
 **/
namespace crab {
namespace domains {

using namespace elina;

template <typename Number, typename VariableName, elina_domain_id_t ElinaDom,
          class Params = ElinaDefaultParams<Number>>
class elina_domain final
    : public abstract_domain_api<
          elina_domain<Number, VariableName, ElinaDom, Params>> {
  using elina_domain_t = elina_domain<Number, VariableName, ElinaDom, Params>;
  using abstract_domain_t = abstract_domain_api<elina_domain_t>;

public:
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
  using number_t = Number;
  using varname_t = VariableName;

private:
  using interval_domain_t = ikos::interval_domain<number_t, varname_t>;
  using bound_t = ikos::bound<number_t>;
  using var_map_t = boost::bimap<variable_t, elina_dim_t>;
  using binding_t = typename var_map_t::value_type;

  static elina_manager_t *s_apman;

  elina_state_ptr m_apstate;
  var_map_t m_var_map;

  bool is_real() const { return std::is_same<Number, ikos::q_number>::value; }

  bool is_integer() const { return !is_real(); }

  static elina_manager_t *get_man() {
    if (!s_apman) {
      switch (ElinaDom) {
      case ELINA_ZONES:
        s_apman = opt_zones_manager_alloc();
        break;
      case ELINA_OCT:
        s_apman = opt_oct_manager_alloc();
        break;
      case ELINA_PK:
        s_apman = opt_pk_manager_alloc(false);
        break;
      }
    }
    return s_apman;
  }

  // Elina assume dimensions [0..intdim-1] correspond to integer
  // variables, and dimensions [intdim..intdim+realdim-1] to real
  // variables. In our case, we don't allow mix.
  size_t get_dims(elina_state_ptr s) const {
    elina_dimension_t dims = _elina_abstract0_dimension(&*s);
    if (is_real()) {
      assert(dims.intdim == 0);
      return dims.realdim;
    } else {
      assert(dims.realdim == 0);
      return dims.intdim;
    }
  }

  size_t get_dims() const { return get_dims(m_apstate); }

  // If v is in the map then it maps v to a dimension, otherwise null
  boost::optional<elina_dim_t> get_var_dim(const var_map_t &m,
                                           variable_t v) const {
    auto it = m.left.find(v);
    if (it != m.left.end())
      return it->second;
    else
      return boost::optional<elina_dim_t>();
  }

  boost::optional<elina_dim_t> get_var_dim(variable_t v) const {
    return get_var_dim(m_var_map, v);
  }

  elina_dim_t get_var_dim_insert(variable_t v) {
    assert(m_var_map.size() == get_dims());
    if (auto dim = get_var_dim(v))
      return *dim;
    else {
      elina_dim_t i = m_var_map.size();
      m_var_map.insert(binding_t(v, i));
      add_dimensions(m_apstate, 1);
      assert(m_var_map.size() == get_dims());
      return i;
    }
  }

  bool has_variable(const var_map_t &m, elina_dim_t i) const {
    return m.right.find(i) != m.right.end();
  }

  bool has_variable(elina_dim_t i) const { return has_variable(m_var_map, i); }

  variable_t get_variable(const var_map_t &m, elina_dim_t i) const {
    auto it = m.right.find(i);
    if (it != m.right.end())
      return it->second;
    CRAB_ERROR("Elina dimension ", i, " is not used!");
  }

  variable_t get_variable(elina_dim_t i) const {
    return get_variable(m_var_map, i);
  }

  void add_dimensions(elina_state_ptr &s, size_t dims) const {
    if (dims <= 0)
      return;
    elina_dimchange_t *dimchange = nullptr;
    if (is_integer()) {
      dimchange = elina_dimchange_alloc(dims, 0);
    } else {
      assert(is_real());
      if (ElinaDom == ELINA_PK) {
        CRAB_ERROR("Elina polyhedra does not support rational coefficients");
      }
      dimchange = elina_dimchange_alloc(0, dims);
    }
    for (unsigned i = 0; i < dims; i++)
      dimchange->dim[i] = get_dims(s); // add dimension at the end

    s = elinaPtr(get_man(), elina_abstract0_add_dimensions(
                                get_man(), false, &*s, dimchange, false));
    elina_dimchange_free(dimchange);
  }

  void remove_dimensions(elina_state_ptr &s,
                         std::vector<elina_dim_t> dims) const {
    if (dims.empty())
      return;

    // Elina assumption: make sure that the removing dimensions
    //                   are in ascending order.
    std::sort(dims.begin(), dims.end());

    elina_dimchange_t *dimchange = nullptr;
    if (is_integer()) {
      dimchange = elina_dimchange_alloc(dims.size(), 0);
    } else {
      assert(is_real());
      dimchange = elina_dimchange_alloc(0, dims.size());
    }
    for (unsigned i = 0; i < dims.size(); i++) {
      // remove dimension dims[i] and shift to the left all the
      // dimensions greater than dims[i]
      dimchange->dim[i] = dims[i];
    }

    s = elinaPtr(get_man(), elina_abstract0_remove_dimensions(get_man(), false,
                                                              &*s, dimchange));
    elina_dimchange_free(dimchange);

#if 0
      crab::outs() << "Removed " << dims.size() << " dimensions\n";
      crab::outs() << "Size = " << get_dims(s) << "\n";
#endif
  }

  bool check_perm(const elina_dimperm_t *perm, size_t size) const {
    // it does not check injectivity
    if (perm->size != size)
      return false;
    for (unsigned i = 0; i < perm->size; i++) {
      if (perm->dim[i] >= size) {
        return false;
      }
    }
    return true;
  }

  var_map_t merge_var_map(const var_map_t &m_x, elina_state_ptr &s_x,
                          const var_map_t &m_y, elina_state_ptr &s_y) const {

    assert(m_x.size() == get_dims(s_x));
    assert(m_y.size() == get_dims(s_y));

    // -- collect all vars from the two maps
    std::set<variable_t> vars;
    for (auto const &px : m_x.left)
      vars.insert(px.first);
    for (auto const &py : m_y.left)
      vars.insert(py.first);

    assert(vars.size() >= get_dims(s_x));
    assert(vars.size() >= get_dims(s_y));

    add_dimensions(s_x, vars.size() - get_dims(s_x));
    add_dimensions(s_y, vars.size() - get_dims(s_y));

    assert(get_dims(s_x) == get_dims(s_y));

    // -- create a fresh map
    var_map_t res;
    for (auto v : vars) {
      elina_dim_t i = res.size();
      assert(i < get_dims(s_x));
      res.insert(binding_t(v, i));
    }

    // build the permutations maps
    elina_dimperm_t *perm_x = elina_dimperm_alloc(get_dims(s_x));
    elina_dimperm_t *perm_y = elina_dimperm_alloc(get_dims(s_x));
    char *xmap1 = (char *)calloc(get_dims(s_x), sizeof(char));
    if (!xmap1)
      CRAB_ERROR("calloc does not have more available memory");
    char *xmap2 = (char *)calloc(get_dims(s_x), sizeof(char));
    if (!xmap2)
      CRAB_ERROR("calloc does not have more available memory");
    for (auto const &px : m_x.left) {
      elina_dim_t ind = res.left.at(px.first);
      perm_x->dim[px.second] = ind;
      // This sets 1 if the index that has been assigned
      assert(px.second < get_dims(s_x));
      xmap1[px.second] = 1;
      // This sets 1 if the value has been assigned
      assert(ind < get_dims(s_x));
      xmap2[ind] = 1;
    }
    elina_dim_t i, counter = 0;
    for (i = 0; i < get_dims(s_x); i++) {
      // If the index has been assigned, skip
      if (xmap1[i])
        continue;
      // Find the next available element that has not been assigned
      while (xmap2[counter])
        counter++;
      perm_x->dim[i] = counter;
      counter++;
    }
    free(xmap1);
    free(xmap2);

    char *ymap1 = (char *)calloc(get_dims(s_x), sizeof(char));
    if (!ymap1)
      CRAB_ERROR("calloc does not have more available memory");
    char *ymap2 = (char *)calloc(get_dims(s_x), sizeof(char));
    if (!ymap2)
      CRAB_ERROR("calloc does not have more available memory");
    for (auto const &py : m_y.left) {
      elina_dim_t ind = res.left.at(py.first);
      perm_y->dim[py.second] = ind;
      assert(py.second < get_dims(s_x));
      ymap1[py.second] = 1;
      assert(ind < get_dims(s_x));
      ymap2[ind] = 1;
    }

    counter = 0;
    for (i = 0; i < get_dims(s_x); i++) {
      if (ymap1[i])
        continue;
      while (ymap2[counter])
        counter++;
      perm_y->dim[i] = counter;
      counter++;
    }

    free(ymap1);
    free(ymap2);

#if 0          
      crab::outs() << "Permutations \n";
      elina_dimperm_fprint(stdout, perm_x);          
      crab::outs() << "Permutations \n";
      elina_dimperm_fprint(stdout, perm_y);
#endif

    assert(check_perm(perm_x, get_dims(s_x)));
    assert(check_perm(perm_y, get_dims(s_x)));

    // apply the permutations
    s_x = elinaPtr(get_man(), elina_abstract0_permute_dimensions(
                                  get_man(), false, &*s_x, perm_x));
    s_y = elinaPtr(get_man(), elina_abstract0_permute_dimensions(
                                  get_man(), false, &*s_y, perm_y));

    elina_dimperm_free(perm_x);
    elina_dimperm_free(perm_y);

    assert(res.size() == get_dims(s_x));
    assert(res.size() == get_dims(s_y));

    return res;
  }

  // --- from crab to elina

  //// using tree expressions
  inline elina_texpr0_t *var2texpr(variable_t v) {
    return elina_texpr0_dim(get_var_dim_insert(v));
  }

  inline elina_texpr0_t *num2texpr(number_t i) const {
    ikos::q_number qi(i);
    return elina_texpr0_cst_scalar_mpq(qi.get_mpq_t());
  }

  inline elina_texpr0_t *intv2texpr(number_t a, number_t b) const {
    ikos::q_number qa(a);
    ikos::q_number qb(b);
    return elina_texpr0_cst_interval_mpq(qa.get_mpq_t(), qb.get_mpq_t());
  }

  static elina_texpr0_t *texpr_add(elina_texpr0_t *a, elina_texpr0_t *b) {
    if (Params::use_integers) {
      return elina_texpr0_binop(ELINA_TEXPR_ADD, a, b, ELINA_RTYPE_INT,
                                ELINA_RDIR_ZERO);
    } else {
      // ELINA_RTYPE_REAL does not have rounding so we choose an
      // arbitrary one
      return elina_texpr0_binop(ELINA_TEXPR_ADD, a, b, ELINA_RTYPE_REAL,
                                ELINA_RDIR_ZERO);
    }
  }
  static elina_texpr0_t *texpr_sub(elina_texpr0_t *a, elina_texpr0_t *b) {
    if (Params::use_integers) {
      return elina_texpr0_binop(ELINA_TEXPR_SUB, a, b, ELINA_RTYPE_INT,
                                ELINA_RDIR_ZERO);
    } else {
      // ELINA_RTYPE_REAL does not have rounding so we choose an
      // arbitrary one
      return elina_texpr0_binop(ELINA_TEXPR_SUB, a, b, ELINA_RTYPE_REAL,
                                ELINA_RDIR_ZERO);
    }
  }
  static elina_texpr0_t *texpr_mul(elina_texpr0_t *a, elina_texpr0_t *b) {
    if (Params::use_integers) {
      return elina_texpr0_binop(ELINA_TEXPR_MUL, a, b, ELINA_RTYPE_INT,
                                ELINA_RDIR_ZERO);
    } else {
      // ELINA_RTYPE_REAL does not have rounding so we choose an
      // arbitrary one
      return elina_texpr0_binop(ELINA_TEXPR_MUL, a, b, ELINA_RTYPE_REAL,
                                ELINA_RDIR_ZERO);
    }
  }
  static elina_texpr0_t *texpr_div(elina_texpr0_t *a, elina_texpr0_t *b) {
    if (Params::use_integers) {
      return elina_texpr0_binop(ELINA_TEXPR_DIV, a, b, ELINA_RTYPE_INT,
                                ELINA_RDIR_ZERO);
    } else {
      // ELINA_RTYPE_REAL does not have rounding so we choose an
      // arbitrary one
      return elina_texpr0_binop(ELINA_TEXPR_DIV, a, b, ELINA_RTYPE_REAL,
                                ELINA_RDIR_ZERO);
    }
  }

  inline elina_texpr0_t *expr2texpr(const linear_expression_t &e) {
    number_t cst = e.constant();
    elina_texpr0_t *res = num2texpr(cst);
    for (auto p : e) {
      elina_texpr0_t *term = texpr_mul(num2texpr(p.first), var2texpr(p.second));
      res = texpr_add(res, term);
    }
    return res;
  }

  inline elina_tcons0_t const2tconst(const linear_constraint_t &cst) {
    linear_expression_t exp = cst.expression();
    if (cst.is_equality()) {
      return elina_tcons0_make(ELINA_CONS_EQ, expr2texpr(exp), NULL);
    } else if (cst.is_inequality()) {
      return elina_tcons0_make(ELINA_CONS_SUPEQ, expr2texpr(-exp), NULL);
    } else if (cst.is_strict_inequality()) {
      return elina_tcons0_make(ELINA_CONS_SUP, expr2texpr(-exp), NULL);
    } else {
      assert(cst.is_disequation());
      return elina_tcons0_make(ELINA_CONS_DISEQ, expr2texpr(exp), NULL);
    }
  }

  //// Using linear expressions and constraints
  elina_linexpr0_t *crab_linexpr_to_elina(const linear_expression_t &e) {
    elina_linexpr0_t *linexpr =
        elina_linexpr0_alloc(ELINA_LINEXPR_SPARSE, e.size());
    if (!linexpr) {
      CRAB_ERROR("elina cannot allocate linear expression");
    }

    ikos::q_number q(e.constant());
    elina_scalar_set_mpq(linexpr->cst.val.scalar, q.get_mpq_t());

    unsigned i = 0;
    for (auto it = e.begin(), et = e.end(); it != et; ++it, ++i) {
      number_t crab_coeff = it->first;
      const variable_t &crab_var = it->second;
      elina_linterm_t *linterm = &linexpr->p.linterm[i];
      linterm->dim = get_var_dim_insert(crab_var);
      ikos::q_number qcoef(crab_coeff);
      elina_scalar_set_mpq(linterm->coeff.val.scalar, qcoef.get_mpq_t());
    }

    return linexpr;
  }

  // for debugging
  void dump(elina_linexpr0_t *e, const var_map_t &m,
            elina_state_ptr apstate) const {
    std::vector<char *> names;
    for (unsigned i = 0; i < get_dims(apstate); i++) {
      std::string varname;
      if (has_variable(m, i))
        varname = get_variable(m, i).name().str();
      else // unused dimension
        varname = std::string("_x") + std::to_string(i);
      char *name = new char[varname.length() + 1];
      strcpy(name, varname.c_str());
      names.push_back(name);
    }
    elina_linexpr0_print(e, &names[0]);
    for (auto n : names) {
      delete[] n;
    }
  }

  elina_lincons0_t crab_linconst_to_elina(const linear_constraint_t &c) {
    const linear_expression_t &e = c.expression();
    if (c.is_equality()) {
      return elina_lincons0_make(ELINA_CONS_EQ, crab_linexpr_to_elina(e),
                                 nullptr);
    } else if (c.is_inequality()) {
      return elina_lincons0_make(ELINA_CONS_SUPEQ, crab_linexpr_to_elina(-e),
                                 nullptr);
    } else if (c.is_strict_inequality()) {
      return elina_lincons0_make(ELINA_CONS_SUP, crab_linexpr_to_elina(-e),
                                 nullptr);
    } else if (c.is_disequation()) {
      return elina_lincons0_make(ELINA_CONS_DISEQ, crab_linexpr_to_elina(e),
                                 nullptr);
    } else {
      CRAB_ERROR("unsupported linear constraint kind");
    }
  }

  // --- from elina to crab

  inline void convert_elina_number(double n, ikos::z_number &res) const {
    res = ikos::z_number((long)n);
  }

  inline void convert_elina_number(double n, ikos::q_number &res) const {
    res = ikos::q_number(n);
  }

  inline void convert_elina_number(mpq_ptr mp, ikos::z_number &res) const {
    ikos::q_number q = ikos::q_number::from_mpq_srcptr(mp);
    res = q.round_to_lower();
  }
  inline void convert_elina_number(mpq_ptr mp, ikos::q_number &res) const {
    res = ikos::q_number::from_mpq_srcptr(mp);
  }

  number_t coeff_to_num(elina_coeff_t *coeff) const {
    assert(coeff->discr == ELINA_COEFF_SCALAR);

    elina_scalar_t *scalar = coeff->val.scalar;
    if (scalar->discr == ELINA_SCALAR_DOUBLE) { // elina uses double
      number_t res;
      convert_elina_number(scalar->val.dbl, res);
      return res;
    } else if (scalar->discr == ELINA_SCALAR_MPQ) {
      number_t res;
      convert_elina_number(scalar->val.mpq, res);
      return res;
    } else {
      CRAB_ERROR("elina translation only covers double or mpq scalars");
    }
  }

  linear_expression_t term2expr(elina_coeff_t *coeff, elina_dim_t i) const {
    return variable_t(get_variable(i)) * coeff_to_num(coeff);
  }

  linear_constraint_t tconst2const(elina_lincons0_t cons) const {
    assert(cons.scalar == NULL); // Not modulo form
    elina_linexpr0_t *linexp = cons.linexpr0;
    assert(elina_linexpr0_is_linear(linexp));

    unsigned i;
    elina_dim_t dim;
    elina_coeff_t *coef;
    linear_expression_t e(0);
    elina_linexpr0_ForeachLinterm(linexp, i, dim, coef) {
      if (elina_coeff_zero(coef))
        continue;
      e = e + term2expr(coef, dim);
    }

    // add constant
    elina_coeff_t *cst = elina_linexpr0_cstref(linexp);
    if (!elina_coeff_zero(cst))
      e = e + coeff_to_num(cst);

    linear_constraint_t res;
    switch (cons.constyp) {
    case ELINA_CONS_EQ:
      // e == k
      res = linear_constraint_t(e, linear_constraint_t::kind_t::EQUALITY);
      break;
    case ELINA_CONS_SUPEQ:
      // e >= k
      e = -e;
      res = linear_constraint_t(e, linear_constraint_t::kind_t::INEQUALITY);
      break;
    case ELINA_CONS_SUP:
      // e > k
      e = -e;
      res = linear_constraint_t(e,
                                linear_constraint_t::kind_t::STRICT_INEQUALITY);
      break;
    case ELINA_CONS_EQMOD:
      res = linear_constraint_t::get_true();
      break;
    case ELINA_CONS_DISEQ:
      // e != k
      res = linear_constraint_t(e, linear_constraint_t::kind_t::DISEQUATION);
      break;
    }
    return res;
  }

  bound_t scalar2bound(elina_scalar_t *scalar) const {
    switch(elina_scalar_infty(scalar)) {
    case -1: /* -oo */
      return bound_t::minus_infinity();
    case 1:  /* +oo */
      return bound_t::plus_infinity();
    default: /* finite */
      assert(elina_scalar_infty(scalar) == 0);
      
      switch(scalar->discr) {
      case ELINA_SCALAR_DOUBLE: {
	number_t val;
	convert_elina_number(scalar->val.dbl, val);
	return bound_t(val);
      }
      case ELINA_SCALAR_MPQ: {
	number_t val;
	convert_elina_number(scalar->val.mpq, val);
	return bound_t(val);
      }
      default:
	CRAB_ERROR("elina translation only covers double or mpq scalars");
      }
    }
  }    

  
  void dump(const var_map_t &m, elina_state_ptr apstate) const {
    crab::outs() << "\nNumber of dimensions=" << get_dims(apstate) << "\n";
    crab::outs() << "variable map [";
    std::vector<char *> names;
    for (unsigned i = 0; i < get_dims(apstate); i++) {
      std::string varname;
      if (has_variable(m, i))
        varname = get_variable(m, i).name().str();
      else // unused dimension
        varname = std::string("_x") + std::to_string(i);
      crab::outs() << i << " -> " << varname << ";";
      char *name = new char[varname.length() + 1];
      strcpy(name, varname.c_str());
      names.push_back(name);
    }
    crab::outs() << "]\n";
    elina_abstract0_fprint(stdout, get_man(), &*apstate, &names[0]);
    for (auto n : names) {
      delete[] n;
    }
  }

  // x != n
  void inequalities_from_disequation(const variable_t &x, number_t n,
                                     linear_constraint_system_t &out) {
    interval_t i = this->operator[](x);
    interval_t ni(n);
    interval_t new_i =
        ikos::linear_interval_solver_impl::trim_interval<interval_t>(
            i, ni);
    if (new_i.is_bottom()) {
      out += linear_constraint_t::get_false();
    } else if (!new_i.is_top() && (new_i <= i)) {
      if (new_i.lb().is_finite()) {
        // strenghten lb
        out += linear_constraint_t(x >= *(new_i.lb().number()));
      }
      if (new_i.ub().is_finite()) {
        // strenghten ub
        out += linear_constraint_t(x <= *(new_i.ub().number()));
      }
    }
  }

  interval_t compute_residual(const linear_expression_t &e,
                              const variable_t &pivot) {
    interval_t residual(-e.constant());
    for (auto kv : e) {
      const variable_t &v = kv.second;
      if (v.index() != pivot.index()) {
        residual = residual - (interval_t(kv.first) * this->operator[](v));
      }
    }
    return residual;
  }

  void inequalities_from_disequation(const linear_expression_t &e,
                                     linear_constraint_system_t &o) {
    for (auto kv : e) {
      variable_t pivot = kv.second;
      interval_t i = compute_residual(e, pivot) / interval_t(kv.first);
      if (auto k = i.singleton()) {
        inequalities_from_disequation(pivot, *k, o);
      }
    }
  }

  // Using elina tree expressions
  void assign_texpr(const variable_t &x, const linear_expression_t &e) {
    elina_texpr0_t *t = expr2texpr(e);
    assert(t);
    auto dim_x = get_var_dim_insert(x);
    m_apstate =
        elinaPtr(get_man(), elina_abstract0_assign_texpr(
                                get_man(), false, &*m_apstate, dim_x, t, NULL));
    elina_texpr0_free(t);
  }

  // Using elina linear expresions
  void assign_linexpr(const variable_t &x, const linear_expression_t &e) {
    elina_linexpr0_t *t = crab_linexpr_to_elina(e);
    assert(t);
    auto dim_x = get_var_dim_insert(x);
    m_apstate =
        elinaPtr(get_man(), elina_abstract0_assign_linexpr(
                                get_man(), false, &*m_apstate, dim_x, t, NULL));
    elina_linexpr0_free(t);
  }

  // Using elina tree expressions
  void apply_texpr(arith_operation_t op, const variable_t &x,
                   const variable_t &y, number_t z) {
    elina_texpr0_t *a = var2texpr(y);
    elina_texpr0_t *b = num2texpr(z);
    elina_texpr0_t *res = nullptr;

    switch (op) {
    case OP_ADDITION:
      res = texpr_add(a, b);
      break;
    case OP_SUBTRACTION:
      res = texpr_sub(a, b);
      break;
    case OP_MULTIPLICATION:
      res = texpr_mul(a, b);
      break;
    case OP_SDIV:
      res = texpr_div(a, b);
      break;
    default:
      CRAB_ERROR("elina operation ", op, " not supported");
    }
    assert(res);
    auto dim_x = get_var_dim_insert(x);
    m_apstate = elinaPtr(
        get_man(), elina_abstract0_assign_texpr(get_man(), false, &*m_apstate,
                                                dim_x, res, NULL));
    elina_texpr0_free(res);
  }

  // Using elina linear expresions
  void apply_linexpr(arith_operation_t op, const variable_t &x,
                     const variable_t &y, number_t z) {
    elina_linexpr0_t *rhs = nullptr;
    switch (op) {
    case OP_ADDITION:
      rhs = crab_linexpr_to_elina(y + z);
      break;
    case OP_SUBTRACTION:
      rhs = crab_linexpr_to_elina(y - z);
      break;
    case OP_MULTIPLICATION:
      rhs = crab_linexpr_to_elina(y * z);
      break;
    case OP_SDIV:
      apply_texpr(op, x, y, z);
      CRAB_LOG("elina", crab::outs() << "--- " << x << ":=" << y << op << z
                                     << " --> " << *this << "\n";);
      return;
    default:
      CRAB_ERROR("elina operation ", op, " not supported");
    }
    assert(rhs);
    auto dim_x = get_var_dim_insert(x);
    m_apstate = elinaPtr(
        get_man(), elina_abstract0_assign_linexpr(get_man(), false, &*m_apstate,
                                                  dim_x, rhs, NULL));
    elina_linexpr0_free(rhs);
  }

  // Using elina tree expressions
  void apply_texpr(arith_operation_t op, const variable_t &x,
                   const variable_t &y, const variable_t &z) {
    elina_texpr0_t *a = var2texpr(y);
    elina_texpr0_t *b = var2texpr(z);
    elina_texpr0_t *res = nullptr;

    switch (op) {
    case OP_ADDITION:
      res = texpr_add(a, b);
      break;
    case OP_SUBTRACTION:
      res = texpr_sub(a, b);
      break;
    case OP_MULTIPLICATION:
      res = texpr_mul(a, b);
      break;
    case OP_SDIV:
      res = texpr_div(a, b);
      break;
    default:
      CRAB_ERROR("elina operation ", op, " not supported");
    }
    assert(res);
    auto dim_x = get_var_dim_insert(x);
    m_apstate = elinaPtr(
        get_man(), elina_abstract0_assign_texpr(get_man(), false, &*m_apstate,
                                                dim_x, res, NULL));
    elina_texpr0_free(res);
  }

  // Using elina linear expresions
  void apply_linexpr(arith_operation_t op, const variable_t &x,
                     const variable_t &y, const variable_t &z) {
    elina_linexpr0_t *rhs = nullptr;
    switch (op) {
    case OP_ADDITION:
      rhs = crab_linexpr_to_elina(y + z);
      break;
    case OP_SUBTRACTION:
      rhs = crab_linexpr_to_elina(y - z);
      break;
    case OP_MULTIPLICATION:
    case OP_SDIV:
      apply_texpr(op, x, y, z);
      CRAB_LOG("elina", crab::outs() << "--- " << x << ":=" << y << op << z
                                     << " --> " << *this << "\n";);
      return;
    default:
      CRAB_ERROR("elina operation ", op, " not supported");
    }
    assert(rhs);
    auto dim_x = get_var_dim_insert(x);
    m_apstate = elinaPtr(
        get_man(), elina_abstract0_assign_linexpr(get_man(), false, &*m_apstate,
                                                  dim_x, rhs, NULL));
    elina_linexpr0_free(rhs);
  }

  // using tree expressions
  void add_tcons(const linear_constraint_system_t &csts) {
    elina_tcons0_array_t array = elina_tcons0_array_make(csts.size());
    unsigned i = 0;
    for (auto cst : csts) {
      elina_tcons0_t tcons = const2tconst(cst);
      array.p[i] = tcons;
      ++i;
    }

    // ///// debugging
    // std::vector<char*> names;
    // for (unsigned i=0; i < get_dims(m_apstate) ; i++){
    // 	std::string varname;
    // 	if (has_variable(m_var_map, i))
    // 	  varname = get_variable(m_var_map, i).name().str();
    // 	else // unused dimension
    // 	  varname = std::string("_x") + std::to_string(i);
    // 	char* name = new char [varname.length() + 1];
    // 	strcpy(name, varname.c_str());
    // 	names.push_back(name);
    // }
    // elina_tcons0_array_fprint(stdout, &array, &names[0]);
    // for (auto n : names) { delete[] n; }

    m_apstate = elinaPtr(get_man(), elina_abstract0_meet_tcons_array(
                                        get_man(), false, &*m_apstate, &array));

    elina_tcons0_array_clear(&array);
  }

  // using linear expressions
  void add_lincons(const linear_constraint_system_t &csts) {
    elina_lincons0_array_t array = elina_lincons0_array_make(csts.size());
    unsigned i = 0;
    for (auto cst : csts) {
      elina_lincons0_t tcons = crab_linconst_to_elina(cst);
      array.p[i] = tcons;
      ++i;
    }
    m_apstate = elinaPtr(get_man(), elina_abstract0_meet_lincons_array(
                                        get_man(), false, &*m_apstate, &array));
    elina_lincons0_array_clear(&array);
  }

  // public:
  //   void print_stats() { elina_abstract0_fprint(stdout, get_man(),
  //   &*m_apstate, NULL); }
  // private:

  elina_domain(elina_state_ptr &&apState, var_map_t &&varMap,
               bool &&compact = true)
      : m_apstate(std::move(apState)), m_var_map(std::move(varMap)) {

    if (compact) {
      std::vector<elina_dim_t> dims;
      var_map_t res;
      /// XXX: we must iterate on the dimension id's to preserve
      /// order between them
      for (auto const &p : m_var_map.right) {
        if (elina_abstract0_is_dimension_unconstrained(get_man(), &*m_apstate,
                                                       p.first)) {
          dims.push_back(p.first);
        } else {
          elina_dim_t i = res.size();
          res.insert(binding_t(p.second, i));
        }
      }
      remove_dimensions(m_apstate, dims);
      std::swap(m_var_map, res);
      assert(m_var_map.size() == get_dims());
    }
  }

public:
  elina_domain(bool isBot = false)
      : m_apstate(elinaPtr(get_man(),
                           (isBot ? elina_abstract0_bottom(get_man(), 0, 0)
                                  : elina_abstract0_top(get_man(), 0, 0)))) {}

  ~elina_domain() = default;

  elina_domain(const elina_domain_t &o)
      : m_apstate(elinaPtr(get_man(),
                           elina_abstract0_copy(get_man(), &*(o.m_apstate)))),
        m_var_map(o.m_var_map) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
  }

  elina_domain(elina_domain_t &&o)
      : m_apstate(std::move(o.m_apstate)), m_var_map(std::move(o.m_var_map)) {}

  elina_domain_t &operator=(const elina_domain_t &o) {
    crab::CrabStats::count(domain_name() + ".count.copy");
    crab::ScopedCrabStats __st__(domain_name() + ".copy");
    if (this != &o) {
      m_apstate =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*(o.m_apstate)));
      m_var_map = o.m_var_map;
    }
    return *this;
  }

  elina_domain_t &operator=(elina_domain_t &&o) {
    if (this != &o) {
      m_apstate = std::move(o.m_apstate);
      m_var_map = std::move(o.m_var_map);
    }
    return *this;
  }

  elina_domain_t make_top() const override {
    elina_domain_t out(false);
    return out;
  }

  elina_domain_t make_bottom() const override {
    elina_domain_t out(true);
    return out;
  }

  void set_to_top() override {
    elina_domain_t abs(false);
    std::swap(*this, abs);
  }

  void set_to_bottom() override {
    elina_domain_t abs(true);
    std::swap(*this, abs);
  }

  bool is_bottom() const override {
    return elina_abstract0_is_bottom(get_man(), &*m_apstate);
  }

  bool is_top() const override {
    return elina_abstract0_is_top(get_man(), &*m_apstate);
  }

  bool operator<=(const elina_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.leq");
    crab::ScopedCrabStats __st__(domain_name() + ".leq");

    if (is_bottom())
      return true;
    else if (o.is_bottom())
      return false;
    else if (o.is_top())
      return true;
    else if (is_top() && !o.is_top())
      return false;
    else if (is_top() && o.is_top())
      return true;
    else {
      elina_state_ptr x =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*m_apstate));
      elina_state_ptr y =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));
      merge_var_map(m_var_map, x, o.m_var_map, y);
      return elina_abstract0_is_leq(get_man(), &*x, &*y);
    }
  }

  void operator|=(const elina_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom()) {
      *this = o;
    } else if (o.is_top()) {
      set_to_top();
    } else if (is_top() || o.is_bottom()) {
      // do nothing
    } else {
      CRAB_LOG("elina", crab::outs()
                            << "JOIN \n\t" << *this << "\n\t" << o << "\n";);
      elina_state_ptr x =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));
      m_var_map =
          std::move(merge_var_map(m_var_map, m_apstate, o.m_var_map, x));
      m_apstate = elinaPtr(
          get_man(), elina_abstract0_join(get_man(), false, &*m_apstate, &*x));
      CRAB_LOG("elina", crab::outs() << "RESULT " << *this << "\n");
    }
  }

  elina_domain_t operator|(const elina_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.join");
    crab::ScopedCrabStats __st__(domain_name() + ".join");

    if (is_bottom()) {
      return o;
    } else if (o.is_top() || is_top()) {
      elina_domain_t res;
      return res;
    } else if (o.is_bottom()) {
      return *this;
    } else {
      elina_state_ptr x =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*m_apstate));
      elina_state_ptr y =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));
      var_map_t m = merge_var_map(m_var_map, x, o.m_var_map, y);
      return elina_domain_t(
          elinaPtr(get_man(), elina_abstract0_join(get_man(), false, &*x, &*y)),
          std::move(m));
    }
  }
  
  void operator&=(const elina_domain_t &o) override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      // do nothing
    } else if (is_top() || o.is_bottom()) {
      *this = o;
    } else {
      elina_state_ptr x =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));
      m_var_map =
	std::move(merge_var_map(m_var_map, m_apstate, o.m_var_map, x));
      m_apstate = elinaPtr(
	get_man(), elina_abstract0_meet(get_man(), false, &*m_apstate, &*x));
    }
  }
  
  elina_domain_t operator&(const elina_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.meet");
    crab::ScopedCrabStats __st__(domain_name() + ".meet");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      elina_state_ptr x =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*m_apstate));
      elina_state_ptr y =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));

      var_map_t m = merge_var_map(m_var_map, x, o.m_var_map, y);
      return elina_domain_t(
          elinaPtr(get_man(), elina_abstract0_meet(get_man(), false, &*x, &*y)),
          std::move(m));
    }
  }

  elina_domain_t operator||(const elina_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");
    // XXX: is_bottom will close the left operand
    // if (is_bottom())
    //	return o;
    // else
    // if (o.is_bottom())
    // return *this;
    // else {
    elina_state_ptr x =
        elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*m_apstate));
    elina_state_ptr y =
        elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));

    var_map_t m = merge_var_map(m_var_map, x, o.m_var_map, y);

    // widening precondition: the old value is included in the new value.
    // widen(old, new) = widen(old,(join(old,new)))    
    elina_state_ptr x_join_y = elinaPtr(get_man(),
					elina_abstract0_join(get_man(), false, &*x, &*y));
    
    return elina_domain_t(
        elinaPtr(get_man(), elina_abstract0_widening(get_man(), &*x, &*x_join_y)),
        std::move(m), false /* do not compact */);
    //}
  }

  elina_lincons0_array_t
  make_thresholds(const elina_domain_t &o,
                  const thresholds<number_t> &ts) const {
    // TODO: make some constraints using the constants from ts
    elina_lincons0_array_t csts = elina_lincons0_array_make(0);
    return csts;
  }

  elina_domain_t widening_thresholds(
      const elina_domain_t &o,
      const thresholds<number_t> &ts) const override {
    crab::CrabStats::count(domain_name() + ".count.widening");
    crab::ScopedCrabStats __st__(domain_name() + ".widening");

    // XXX: is_bottom will close the left operand
    // if (is_bottom())
    //	return o;
    // else
    // if (o.is_bottom())
    // return *this;
    // else {
    elina_state_ptr x =
        elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*m_apstate));
    elina_state_ptr y =
        elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));

    var_map_t m = merge_var_map(m_var_map, x, o.m_var_map, y);
//////
// We cannot refine the result of widening with
// widening w/ thresholds over intervals because it might
// cause non-termination.
/////
// This causes a loss of precision in a couple of tests:
// - tests/domains/test2-rat.cc
// - tests/domains/test3-rat.cc
/////
#if 0
      // widening w/o thresholds in the elina domain
      elina_domain_t res(elinaPtr(get_man(), 
				  elina_abstract0_widening(get_man(), 
							   &*x, &*y)),
			 std::move(m));
      // widening w/ thresholds in the interval domain
      auto intv_this  = this->to_interval_domain();
      auto intv_o     = o.to_interval_domain();
      auto intv_widen = intv_this.widening_thresholds(intv_o, ts);	    
      // refine the elina domain using the widen intervals
      elina_domain_t elina_intv_widen;
      elina_intv_widen += intv_widen.to_linear_constraint_system();
      return res & elina_intv_widen;
#else
    // widening precondition: the old value is included in the new value.
    // widen(old, new) = widen(old,(join(old,new)))    
    elina_state_ptr x_join_y = elinaPtr(get_man(),
					elina_abstract0_join(get_man(), false, &*x, &*y));
      
    elina_lincons0_array_t csts = make_thresholds(o, ts);
    elina_domain_t res(elinaPtr(get_man(), elina_abstract0_widening_threshold(
                                               get_man(), &*x, &*x_join_y, &csts)),
                       std::move(m), false /* do not compact */);
    elina_lincons0_array_clear(&csts);
    return res;
#endif
    //}
  }

  elina_domain_t operator&&(const elina_domain_t &o) const override {
    crab::CrabStats::count(domain_name() + ".count.narrowing");
    crab::ScopedCrabStats __st__(domain_name() + ".narrowing");

    if (is_bottom() || o.is_top()) {
      return *this;
    } else if (is_top() || o.is_bottom()) {
      return o;
    } else {
      elina_state_ptr x =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*m_apstate));
      elina_state_ptr y =
          elinaPtr(get_man(), elina_abstract0_copy(get_man(), &*o.m_apstate));

      var_map_t m = merge_var_map(m_var_map, x, o.m_var_map, y);
      switch (ElinaDom) {
      case ELINA_OCT:
        return elina_domain_t(
            elinaPtr(get_man(),
                     elina_abstract0_opt_oct_narrowing(get_man(), &*x, &*y)),
            std::move(m));
      case ELINA_ZONES:
      case ELINA_PK:
        // CRAB_WARN("used meet instead of narrowing: \n",
        //           "make sure only a finite number of descending iterations
        //           are run.");
        return elina_domain_t(
            elinaPtr(get_man(),
                     elina_abstract0_meet(get_man(), false, &*x, &*y)),
            std::move(m));
      }
    }
  }

  void project(const variable_vector_t &vars) override {
    crab::CrabStats::count(domain_name() + ".count.project");
    crab::ScopedCrabStats __st__(domain_name() + ".project");

    if (is_bottom() || is_top())
      return;

    if (vars.empty()) {
      set_to_top();
      return;
    }

    std::set<variable_t> s1, s2;
    variable_vector_t s3;
    for (auto p : m_var_map.left)
      s1.insert(p.first);
    s2.insert(vars.begin(), vars.end());
    std::set_difference(s1.begin(), s1.end(), s2.begin(), s2.end(),
                        std::back_inserter(s3));
    forget(s3);
  }

  void forget(const variable_vector_t &vars) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top())
      return;

    std::vector<elina_dim_t> vector_dims;
    std::set<elina_dim_t> set_dims;

    for (variable_t v : vars) {
      if (auto dim = get_var_dim(v)) {
        vector_dims.push_back(*dim);
        set_dims.insert(*dim);
      }
    }

    if (vector_dims.empty())
      return;

    m_apstate =
        elinaPtr(get_man(), elina_abstract0_forget_array(
                                get_man(), false, &*m_apstate, &vector_dims[0],
                                vector_dims.size(), false));

    // -- Remove forgotten dimensions while compacting
    var_map_t res;
    /// XXX: we must iterate on the dimension id's to preserve
    /// order between them
    for (auto const &p : m_var_map.right) {
      if (set_dims.count(p.first) <= 0) {
        elina_dim_t i = res.size();
        res.insert(binding_t(p.second, i));
      }
    }

    remove_dimensions(m_apstate, vector_dims);
    std::swap(m_var_map, res);
    CRAB_LOG("elina", crab::outs() << "--- "
                                   << "Forget {";
             for (variable_t v
                  : vars) { crab::outs() << v << ";"; }
             // crab::outs()<< "} Elina dimensions ={";
             // for(unsigned i=0,e=vector_dims.size();i<e;++i) {
             // 	 crab::outs() << vector_dims[i] << ";";
             // }
             crab::outs()
             << "}\n";
             crab::outs() << *this << "\n";);
  }

  void operator-=(const variable_t &var) override {
    crab::CrabStats::count(domain_name() + ".count.forget");
    crab::ScopedCrabStats __st__(domain_name() + ".forget");

    if (is_bottom() || is_top())
      return;

    std::vector<elina_dim_t> vector_dims;
    if (auto dim = get_var_dim(var)) {
      vector_dims.push_back(*dim);
      m_apstate =
          elinaPtr(get_man(), elina_abstract0_forget_array(
                                  get_man(), false, &*m_apstate,
                                  &vector_dims[0], vector_dims.size(), false));
      // -- Remove forgotten dimensions while compacting
      var_map_t res;
      /// XXX: we must iterate on the dimension id's to preserve
      /// order between them
      for (auto const &p : m_var_map.right) {
        if (p.first != *dim) {
          elina_dim_t i = res.size();
          res.insert(binding_t(p.second, i));
        }
      }
      remove_dimensions(m_apstate, vector_dims);
      std::swap(m_var_map, res);

      CRAB_LOG("elina", crab::outs() << "--- "
                                     << "Forget "
                                     << var
                                     //<< " Elina dimension=" << *dim
                                     << "\n"
                                     << *this << "\n";);
    }
  }

  interval_t operator[](const variable_t &v) override {
    return at(v);
  }
  
  interval_t at(const variable_t &v) const override {
    crab::CrabStats::count(domain_name() + ".count.to_intervals");
    crab::ScopedCrabStats __st__(domain_name() + ".to_intervals");

    if (is_bottom()) {
      return interval_t::bottom();
    } else if (is_top()) {
      return interval_t::top();
    }

    if (auto dim = get_var_dim(v)) {
      elina_interval_t *intv =
          elina_abstract0_bound_dimension(get_man(), &*m_apstate, *dim);

      if (elina_interval_is_top(intv)) {
        elina_interval_free(intv);
        return interval_t::top();
      } else {
	bound_t lb = scalar2bound(intv->inf);
	bound_t ub = scalar2bound(intv->sup);
	elina_interval_free(intv);
	return interval_t(lb,ub);
      }
    } else {
      return interval_t::top();
    }
  }

  void set(const variable_t &v, interval_t ival) {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (is_bottom()) {
      return;
    }

    // -- forget v
    *this -= v;

    // -- add constraints v >= lb and v <= ub
    linear_constraint_system_t csts;
    auto lb = ival.lb();
    if (lb.is_finite()) {
      // v >= lb <--> -v + lb <= 0
      assert(lb.number());
      linear_expression_t e = (number_t(-1) * v) + *(lb.number());
      csts += (linear_constraint_t(e, linear_constraint_t::kind_t::INEQUALITY));
    }
    auto ub = ival.ub();
    if (ub.is_finite()) {
      // v <= ub <--> v - ub <= 0
      assert(ub.number());
      linear_expression_t e = (v - *(ub.number()));
      csts += (linear_constraint_t(e, linear_constraint_t::kind_t::INEQUALITY));
    }

    if (csts.size() > 0)
      *this += csts;

    CRAB_LOG("elina", crab::outs() << "--- " << v << " := " << ival << "\n"
                                   << *this << "\n";);
  }

  
  void operator+=(const linear_constraint_system_t &_csts) override {
    crab::CrabStats::count(domain_name() + ".count.add_constraints");
    crab::ScopedCrabStats __st__(domain_name() + ".add_constraints");

    if (is_bottom())
      return;

    if (_csts.is_false()) {
      set_to_bottom();
      return;
    }

    if (_csts.is_true()) {
      return;
    }

    // XXX: filter out unsigned linear inequalities, and analyze
    //      separately disequalities because elina does not support
    //      them.
    linear_constraint_system_t csts;
    for (auto const &c : _csts) {
      if (c.is_strict_inequality()) {
        // We try to convert a strict to non-strict.
        csts +=
            ikos::linear_constraint_impl::strict_to_non_strict_inequality(c);
      } else if (c.is_disequation()) {
	// We try to convert a disequation into a strict inequality
	constraint_simp_domain_traits<elina_domain_t>::lower_disequality(*this, c, csts);
        // We try to convert a disequation into conjunctive inequalities
        inequalities_from_disequation(c.expression(), csts);
      } else {
        csts += c;
      }
    }

    if (csts.is_false()) {
      // csts can be false after breaking disequalities into
      // inequalities
      set_to_bottom();
      return;
    }

    if (csts.is_true()) {
      return;
    }

    CRAB_LOG("elina", crab::outs() << "--- "
                                   << "Assume " << csts << " --> ";);

    if (crab_domain_params_man::get().elina_use_tree_expressions()) {
      add_tcons(csts);
    } else {
      add_lincons(csts);
    }

    CRAB_LOG("elina", crab::outs() << *this << "\n";);
  }

  DEFAULT_ENTAILS(elina_domain_t)
  
  void assign(const variable_t &x, const linear_expression_t &e) override {
    crab::CrabStats::count(domain_name() + ".count.assign");
    crab::ScopedCrabStats __st__(domain_name() + ".assign");

    if (is_bottom())
      return;

    CRAB_LOG("elina", crab::outs() << "--- " << x << ":=" << e << " --> ";);

    if (ElinaDom == ELINA_ZONES) {
      // Zones domain does not implement assignment of linear
      // expression
      assign_texpr(x, e);
    } else {
      if (crab_domain_params_man::get().elina_use_tree_expressions()) {
        assign_texpr(x, e);
      } else {
        assign_linexpr(x, e);
      }
    }

    CRAB_LOG("elina", crab::outs() << *this << "\n";);
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             number_t z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;

    if (op >= OP_ADDITION && op <= OP_SDIV) {
      if (ElinaDom == ELINA_ZONES) {
        // Zones domain does not implement assignment of linear
        // expression
        apply_texpr(op, x, y, z);
      } else {
        if (crab_domain_params_man::get().elina_use_tree_expressions()) {
          apply_texpr(op, x, y, z);
        } else {
          apply_linexpr(op, x, y, z);
        }
      }

      CRAB_LOG("elina", crab::outs() << "--- " << x << ":=" << y << op << z
                                     << " --> " << *this << "\n";);
    } else {
      // Convert to intervals and perform the operation
      interval_t yi = operator[](y);
      interval_t zi(z);
      interval_t xi = interval_t::top();
      switch (op) {
      case OP_UDIV:
        xi = yi.UDiv(zi);
        break;
      case OP_SREM:
        xi = yi.SRem(zi);
        break;
      case OP_UREM:
        xi = yi.URem(zi);
        break;
      default:
        CRAB_ERROR("elina operation ", op, " not supported");
      }
      set(x, xi);
    }
  }

  void apply(arith_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom())
      return;

    if (op >= OP_ADDITION && op <= OP_SDIV) {
      if (ElinaDom == ELINA_ZONES) {
        // Zones domain does not implement assignment of linear
        // expression
        apply_texpr(op, x, y, z);
      } else {
        if (crab_domain_params_man::get().elina_use_tree_expressions()) {
          apply_texpr(op, x, y, z);
        } else {
          apply_linexpr(op, x, y, z);
        }
      }

      CRAB_LOG("elina", crab::outs() << "--- " << x << ":=" << y << op << z
                                     << " --> " << *this << "\n";);
    } else {
      // Convert to intervals and perform the operation
      interval_t yi = operator[](y);
      interval_t zi = operator[](z);
      interval_t xi = interval_t::top();

      switch (op) {
      case OP_UDIV:
        xi = yi.UDiv(zi);
        break;
      case OP_SREM:
        xi = yi.SRem(zi);
        break;
      case OP_UREM:
        xi = yi.URem(zi);
        break;
      default:
        CRAB_ERROR("elina operation ", op, " not supported");
      }
      set(x, xi);
    }
  }

  void apply(int_conv_operation_t op, const variable_t &dst,
             const variable_t &src) override {
    int_cast_domain_traits<elina_domain_t>::apply(*this, op, dst, src);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             const variable_t &z) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom()) {
      return;
    }

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi = operator[](z);
    interval_t xi = interval_t::top();
    switch (op) {
    case OP_AND:
      xi = yi.And(zi);
      break;
    case OP_OR:
      xi = yi.Or(zi);
      break;
    case OP_XOR:
      xi = yi.Xor(zi);
      break;
    case OP_SHL:
      xi = yi.Shl(zi);
      break;
    case OP_LSHR:
      xi = yi.LShr(zi);
      break;
    case OP_ASHR:
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }

  void apply(bitwise_operation_t op, const variable_t &x, const variable_t &y,
             number_t k) override {
    crab::CrabStats::count(domain_name() + ".count.apply");
    crab::ScopedCrabStats __st__(domain_name() + ".apply");

    if (is_bottom()) {
      return;
    }

    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi(k);
    interval_t xi = interval_t::top();
    switch (op) {
    case OP_AND:
      xi = yi.And(zi);
      break;
    case OP_OR:
      xi = yi.Or(zi);
      break;
    case OP_XOR:
      xi = yi.Xor(zi);
      break;
    case OP_SHL:
      xi = yi.Shl(zi);
      break;
    case OP_LSHR:
      xi = yi.LShr(zi);
      break;
    case OP_ASHR:
      xi = yi.AShr(zi);
      break;
    }
    set(x, xi);
  }

  DEFAULT_SELECT(elina_domain_t)
  DEFAULT_WEAK_ASSIGN(elina_domain_t)      
  
  void backward_assign(const variable_t &x, const linear_expression_t &e,
                       const elina_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_assign");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_assign");

    if (is_bottom()) {
      return;
    }

    if (invariant.is_bottom()) {
      set_to_bottom();
      return;
    }

    /// Elina substitution only implemented for linear expressions.
    elina_linexpr0_t *rhs = crab_linexpr_to_elina(e);
    elina_dim_t dim_x = get_var_dim_insert(x);

    // ensure that m_apstate and invariant.m_apstate have the same
    // dimensions.
    elina_state_ptr fwd_inv = elinaPtr(
        get_man(), elina_abstract0_copy(get_man(), &*invariant.m_apstate));
    m_var_map =
        merge_var_map(m_var_map, m_apstate, invariant.m_var_map, fwd_inv);
    m_apstate = elinaPtr(
        get_man(), elina_abstract0_substitute_linexpr(
                       get_man(), false, &*m_apstate, dim_x, rhs, &*fwd_inv));
    elina_linexpr0_free(rhs);
    CRAB_LOG("elina", crab::outs() << "--- " << x << " :=_bwd " << e << " --> "
                                   << *this << "\n";);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, number_t z,
                      const elina_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    CRAB_LOG("elina", crab::outs() << "--- " << x << " :=_bwd " << y << op << z
                                   << "\n";);

    if (is_bottom()) {
      return;
    }

    if (invariant.is_bottom()) {
      set_to_bottom();
      return;
    }

    /// Elina substitution only implemented for linear expressions.
    elina_linexpr0_t *rhs = nullptr;
    switch (op) {
    case OP_ADDITION:
      rhs = crab_linexpr_to_elina(y + z);
      break;
    case OP_SUBTRACTION:
      rhs = crab_linexpr_to_elina(y - z);
      break;
    case OP_MULTIPLICATION:
      rhs = crab_linexpr_to_elina(y * z);
      break;
    case OP_SDIV:
      CRAB_WARN("backward signed division not implemented in elina");
      operator-=(x);
      return;
    default:
      CRAB_ERROR("elina operation ", op, " not supported");
    }
    assert(rhs);
    elina_dim_t dim_x = get_var_dim_insert(x);

    // ensure that m_apstate and invariant.m_apstate have the same
    // dimensions.
    elina_state_ptr fwd_inv = elinaPtr(
        get_man(), elina_abstract0_copy(get_man(), &*invariant.m_apstate));

    m_var_map =
        merge_var_map(m_var_map, m_apstate, invariant.m_var_map, fwd_inv);

    m_apstate = elinaPtr(
        get_man(), elina_abstract0_substitute_linexpr(
                       get_man(), false, &*m_apstate, dim_x, rhs, &*fwd_inv));

    elina_linexpr0_free(rhs);
    CRAB_LOG("elina", crab::outs() << " --> " << *this << "\n";);
  }

  void backward_apply(arith_operation_t op, const variable_t &x,
                      const variable_t &y, const variable_t &z,
                      const elina_domain_t &invariant) override {
    crab::CrabStats::count(domain_name() + ".count.backward_apply");
    crab::ScopedCrabStats __st__(domain_name() + ".backward_apply");

    if (is_bottom()) {
      return;
    }

    if (invariant.is_bottom()) {
      set_to_bottom();
      return;
    }

    /// Elina substitution only implemented for linear expressions.
    elina_linexpr0_t *rhs = nullptr;
    switch (op) {
    case OP_ADDITION:
      rhs = crab_linexpr_to_elina(y + z);
      break;
    case OP_SUBTRACTION:
      rhs = crab_linexpr_to_elina(y - z);
      break;
    case OP_MULTIPLICATION:
    case OP_SDIV:
      CRAB_WARN("backward non-linear ", op, " not implemented in elina");
      operator-=(x);
      return;
    default:
      CRAB_ERROR("elina operation ", op, " not supported");
    }
    assert(rhs);
    elina_dim_t dim_x = get_var_dim_insert(x);

    // ensure that m_apstate and invariant.m_apstate have the same
    // dimensions.
    elina_state_ptr fwd_inv = elinaPtr(
        get_man(), elina_abstract0_copy(get_man(), &*invariant.m_apstate));

    m_var_map =
        merge_var_map(m_var_map, m_apstate, invariant.m_var_map, fwd_inv);
    m_apstate = elinaPtr(
        get_man(), elina_abstract0_substitute_linexpr(
                       get_man(), false, &*m_apstate, dim_x, rhs, &*fwd_inv));
    elina_linexpr0_free(rhs);
    CRAB_LOG("elina", crab::outs() << "--- " << x << ":=_bwd " << y << op << z
                                   << " --> " << *this << "\n";);
  }

  /// Elina domains implement only standard abstract operations of a
  /// numerical domain so it is intended to be used as a leaf domain
  /// in the hierarchy of domains.
  BOOL_OPERATIONS_NOT_IMPLEMENTED(elina_domain_t)
  ARRAY_OPERATIONS_NOT_IMPLEMENTED(elina_domain_t)
  REGION_AND_REFERENCE_OPERATIONS_NOT_IMPLEMENTED(elina_domain_t)
  
  interval_domain_t to_interval_domain() {
    crab::CrabStats::count(domain_name() + ".count.to_interval_domain");
    crab::ScopedCrabStats __st__(domain_name() + ".to_interval_domain");

    if (is_bottom())
      return interval_domain_t::bottom();
    if (is_top())
      return interval_domain_t::top();

    interval_domain_t res;
    for (auto &px : m_var_map.left)
      res.set(px.first, this->operator[](px.first));
    return res;
  }

  linear_constraint_system_t to_linear_constraint_system() const override {
    crab::CrabStats::count(domain_name() +
                           ".count.to_linear_constraint_system");
    crab::ScopedCrabStats __st__(domain_name() +
                                 ".to_linear_constraint_system");

    linear_constraint_system_t csts;
    if (is_bottom()) {
      csts += linear_constraint_t::get_false();
    } else if (is_top()) {
      csts += linear_constraint_t::get_true();
    } else {
      // to_lincons_array calls closure
      elina_lincons0_array_t lcons_arr =
          elina_abstract0_to_lincons_array(get_man(), &*m_apstate);
      for (unsigned i = 0; i < lcons_arr.size; i++) {
        csts += tconst2const(lcons_arr.p[i]);
      }

      elina_lincons0_array_clear(&lcons_arr);
    }
    return csts;
  }

  disjunctive_linear_constraint_system_t
  to_disjunctive_linear_constraint_system() const override {
    auto lin_csts = to_linear_constraint_system();
    if (lin_csts.is_false()) {
      return disjunctive_linear_constraint_system_t(true /*is_false*/);
    } else if (lin_csts.is_true()) {
      return disjunctive_linear_constraint_system_t(false /*is_false*/);
    } else {
      return disjunctive_linear_constraint_system_t(lin_csts);
    }
  }

  void rename(const variable_vector_t &from,
              const variable_vector_t &to) override {
    crab::CrabStats::count(domain_name() + ".count.rename");
    crab::ScopedCrabStats __st__(domain_name() + ".rename");

    if (is_top() || is_bottom())
      return;

    CRAB_LOG("elina", crab::outs() << "Renaming {"; for (auto v
                                                         : from) crab::outs()
                                                    << v << ";";
             crab::outs() << "} with "; for (auto v
                                             : to) crab::outs()
                                        << v << ";";
             crab::outs() << "}:\n"; crab::outs() << *this << "\n";);

    for (unsigned i = 0, sz = from.size(); i < sz; ++i) {
      variable_t v = from[i];
      variable_t new_v = to[i];
      if (v == new_v) { // nothing to rename
        continue;
      }

      { // We do garbage collection of unconstrained variables only
        // after joins so it's possible to find new_v but we are ok as
        // long as it's unconstrained.
        auto it = m_var_map.left.find(new_v);
        if (it != m_var_map.left.end()) {
          if (auto dim = get_var_dim(new_v)) {
            if (!elina_abstract0_is_dimension_unconstrained(
                    get_man(), &*m_apstate, *dim)) {
              CRAB_ERROR(domain_name() + "::rename assumes that ", new_v,
                         " does not exist");
            }
          }
        }
      }

      auto it = m_var_map.left.find(v);
      if (it != m_var_map.left.end()) {
        elina_dim_t dim = it->second;
        m_var_map.left.erase(it);
        m_var_map.insert(binding_t(new_v, dim));
      }
    }

    CRAB_LOG("elina", crab::outs() << "RESULT=" << *this << "\n");
  }

  void expand(const variable_t &x, const variable_t &dup) override {
    crab::CrabStats::count(domain_name() + ".count.expand");
    crab::ScopedCrabStats __st__(domain_name() + ".expand");

    if (is_bottom() || is_top())
      return;

    if (get_var_dim(dup)) {
      CRAB_ERROR("expand second parameter ", dup,
                 " cannot be already a variable in the elina domain ", *this);
    }

    // --- increases number of dimensions by one
    auto dim_x = get_var_dim_insert(x);
    m_apstate =
        elinaPtr(get_man(), elina_abstract0_expand(get_man(), false,
                                                   &*m_apstate, dim_x, 1));

    // --- the additional dimension is put at the end of integer
    //     dimensions.
    m_var_map.insert(binding_t(dup, get_dims() - 1));
  }

  void normalize() override {
    crab::CrabStats::count(domain_name() + ".count.normalize");
    crab::ScopedCrabStats __st__(domain_name() + ".normalize");

    elina_abstract0_canonicalize(get_man(), &*m_apstate);
  }

  // reduce the size of the internal representation
  void minimize() override {
    crab::CrabStats::count(domain_name() + ".count.minimize");
    crab::ScopedCrabStats __st__(domain_name() + ".minimize");

    std::vector<elina_dim_t> dims;
    var_map_t res;
    for (auto const &p : m_var_map.right) {
      if (elina_abstract0_is_dimension_unconstrained(get_man(), &*m_apstate,
                                                     p.first)) {
        dims.push_back(p.first);
      } else {
        elina_dim_t i = res.size();
        res.insert(binding_t(p.second, i));
      }
    }
    remove_dimensions(m_apstate, dims);
    std::swap(m_var_map, res);
    assert(m_var_map.size() == get_dims());
  }

  /* begin intrinsics operations */
  void intrinsic(std::string name,
		 const variable_or_constant_vector_t &inputs,
                 const variable_vector_t &outputs) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }

  void backward_intrinsic(std::string name,
			  const variable_or_constant_vector_t &inputs,
                          const variable_vector_t &outputs,
                          const elina_domain_t &invariant) override {
    CRAB_WARN("Intrinsics ", name, " not implemented by ", domain_name());
  }
  /* end intrinsics operations */

  void dump() const { dump(m_var_map, m_apstate); }

  void write(crab_os &o) const override {
    crab::CrabStats::count(domain_name() + ".count.write");
    crab::ScopedCrabStats __st__(domain_name() + ".write");

#if 0
      o << "\nElina internal representation:";
      dump();
#endif

    if (is_bottom()) {
      o << "_|_";
    } else if (is_top()) {
      o << "{}";
    } else {
      linear_constraint_system_t inv = to_linear_constraint_system();
      o << inv;
    }
  }

  std::string domain_name() const override {
    switch (ElinaDom) {
    case ELINA_ZONES:
      return "ElinaZones";
    case ELINA_OCT:
      return "ElinaOctagon";
    case ELINA_PK:
      return "ElinaPolyhedra";
    }
  }
};

// --- global datastructures
template <typename N, typename V, elina_domain_id_t D, class P>
elina_manager_t *elina_domain<N, V, D, P>::s_apman = nullptr;

} // namespace domains
} // namespace crab
#endif

namespace crab {
namespace domains {
template <typename Number, typename VariableName, elina_domain_id_t ElinaDom,
          class Params>
struct abstract_domain_traits<
    elina_domain<Number, VariableName, ElinaDom, Params>> {
  using number_t = Number;
  using varname_t = VariableName;
};
} // namespace domains
} // namespace crab
