#include <crab/config.h>

#include "../common.hpp"
#include "../program_options.hpp"
#ifdef HAVE_LDD
#include <crab/domains/ldd/ldd.hpp>
#endif
#include <boost/optional.hpp>
#include <crab/domains/boxes.hpp>
using namespace std;
using namespace ikos;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
#ifdef HAVE_LDD
using namespace crab::domains::ldd;
#endif

#define RATIONALS

#ifdef RATIONALS
using number_t = q_number;
#else
using number_t = z_number;
#endif

using linear_constraint_t = linear_constraint<number_t, varname_t>;
using linear_expression_t = linear_expression<number_t, varname_t>;
using interval_t = interval<number_t>;

#ifdef HAVE_LDD

LddManager *create_ldd_man(size_t num_vars) {
  DdManager *cudd = Cudd_Init(0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
#ifndef RATIONALS
  theory_t *theory = tvpi_create_boxz_theory(num_vars);
#else
  theory_t *theory = tvpi_create_box_theory(num_vars);
#endif

  LddManager *ldd = Ldd_Init(cudd, theory);

  const bool dvo = true;
  if (dvo)
    Cudd_AutodynEnable(cudd, CUDD_REORDER_GROUP_SIFT);
  return ldd;
}

void destroy_ldd_man(LddManager *ldd) {
  DdManager *cudd = NULL;
  theory_t *theory = NULL;
  if (ldd) {
    cudd = Ldd_GetCudd(ldd);
    theory = Ldd_GetTheory(ldd);
    Ldd_Quit(ldd);
  }

  if (theory)
    tvpi_destroy_theory(theory);
  if (cudd)
    Cudd_Quit(cudd);
}

void dump(LddNodePtr val, LddManager *man) {
  DdManager *cudd = Ldd_GetCudd(man);
  FILE *fp = Cudd_ReadStdout(cudd);
  Cudd_SetStdout(cudd, stderr);
  if (val.get() == Ldd_GetTrue(man))
    crab::outs() << "true\n";
  else if (val.get() == Ldd_GetFalse(man))
    crab::outs() << "false\n";
  else
    Ldd_PrintMinterm(man, val.get());
  Cudd_SetStdout(cudd, fp);
}

LddNodePtr top(LddManager *man) { return lddPtr(man, Ldd_GetTrue(man)); }

LddNodePtr bot(LddManager *man) { return lddPtr(man, Ldd_GetFalse(man)); }

bool isBot(LddManager *man, LddNodePtr v) { return &*v == Ldd_GetFalse(man); }

bool isTop(LddManager *man, LddNodePtr v) { return &*v == Ldd_GetTrue(man); }

LddNodePtr approx(LddManager *man, LddNodePtr v) {
  return lddPtr(man, Ldd_TermMinmaxApprox(man, &*v));
}

LddNodePtr meet(LddManager *man, LddNodePtr n1, LddNodePtr n2) {
  return lddPtr(man, Ldd_And(man, &*n1, &*n2));
}

LddNodePtr convex_join(LddManager *man, LddNodePtr n1, LddNodePtr n2) {
  return approx(man, lddPtr(man, Ldd_Or(man, &*n1, &*n2)));
}

LddNodePtr join(LddManager *man, LddNodePtr n1, LddNodePtr n2) {
  return lddPtr(man, Ldd_Or(man, &*n1, &*n2));
}

bool isLeq(LddManager *man, LddNodePtr n1, LddNodePtr n2) {
  return Ldd_TermLeq(man, &(*n1), &(*n2));
}

LddNodePtr widen(LddManager *man, LddNodePtr n1, LddNodePtr n2) {
  return lddPtr(man, Ldd_IntervalWiden(man, &*n1, &*n2));
}

int getVarId(varname_t v) { return v.index(); }

/** return term for variable v, neg for negation of variable */
linterm_t termForVal(LddManager *man, varname_t v, bool neg = false) {
  int varId = getVarId(v);
  int sgn = neg ? -1 : 1;
  linterm_t term =
      Ldd_GetTheory(man)->create_linterm_sparse_si(&varId, &sgn, 1);
  return term;
}

/** v := k, where k is a constant */
LddNodePtr assign(LddManager *man, LddNodePtr n, varname_t v, number_t k) {
  if (isBot(man, n))
    return n;

  q_number q(k);
  linterm_t t = termForVal(man, v);
  constant_t kkk = (constant_t)tvpi_create_cst(q.get_mpq_t());
  LddNodePtr newVal =
      lddPtr(man, Ldd_TermReplace(man, &(*n), t, NULL, NULL, kkk, kkk));
  Ldd_GetTheory(man)->destroy_cst(kkk);

  assert(!isBot(man, newVal));
  return newVal;
}

/** v := k, where k is an interval */
LddNodePtr assign(LddManager *man, LddNodePtr n, varname_t v, interval_t ival) {

  if (isBot(man, n))
    return n;

  constant_t kmin = NULL, kmax = NULL;

  if (boost::optional<number_t> l = ival.lb().number()) {
    q_number q(*l);
    kmin = (constant_t)tvpi_create_cst(q.get_mpq_t());
  }

  if (boost::optional<number_t> u = ival.ub().number()) {
    q_number q(*u);
    kmax = (constant_t)tvpi_create_cst(q.get_mpq_t());
  }

  linterm_t t = termForVal(man, v);
  LddNodePtr newVal =
      lddPtr(man, Ldd_TermReplace(man, &(*n), t, NULL, NULL, kmin, kmax));

  if (kmin)
    Ldd_GetTheory(man)->destroy_cst(kmin);
  if (kmax)
    Ldd_GetTheory(man)->destroy_cst(kmax);

  assert(!isBot(man, newVal));
  return newVal;
}

/** v := a * u + k, where a, k are constants and u variable */
LddNodePtr apply(LddManager *man, const LddNodePtr n, varname_t v, varname_t u,
                 number_t a, number_t k) {
  if (isTop(man, n) || isBot(man, n))
    return n;

  linterm_t t = termForVal(man, v);
  linterm_t r = termForVal(man, u);

  q_number qa(a);
  q_number qk(k);

  constant_t aaa = (constant_t)tvpi_create_cst(qa.get_mpq_t());
  constant_t kkk = (constant_t)tvpi_create_cst(qk.get_mpq_t());
  LddNodePtr newVal =
      lddPtr(man, Ldd_TermReplace(man, &(*n), t, r, aaa, kkk, kkk));
  Ldd_GetTheory(man)->destroy_cst(aaa);
  Ldd_GetTheory(man)->destroy_cst(kkk);

  assert(!isBot(man, newVal));
  return newVal;
}

/** v := u */
LddNodePtr assign(LddManager *man, const LddNodePtr n, varname_t v,
                  varname_t u) {
  return apply(man, n, v, u, 1, 0);
}
#endif

int main(int argc, char **argv) {
#ifdef HAVE_LDD

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  LddManager *man = create_ldd_man(3000);
  Ldd_SanityCheck(man);

  variable_factory_t vfac;
  varname_t x = vfac["x"];
  varname_t y = vfac["y"];

  {
    LddNodePtr s1 = top(man);
    LddNodePtr s11 = assign(man, s1, x, 5.5);
    LddNodePtr s111 = assign(man, s11, y, 8.3);

    LddNodePtr s2 = top(man);
    LddNodePtr s22 = assign(man, s2, x, 5.5);
    LddNodePtr s222 = assign(man, s22, y, 12.2);

    dump(s111, man);
    dump(s222, man);
    crab::outs() << "Join \n";
    LddNodePtr s3 = join(man, s111, s222);
    dump(s3, man);
    crab::outs() << "Widening \n";
    LddNodePtr s4 = widen(man, s111, s222);
    dump(s4, man);
  }

#if 1
  {
    // To reproduce an old bug with boxes.
    // Assertion failed: (level < (unsigned)
    // cuddI(unique,Cudd_Regular(E)->index)), function cuddUniqueInter, file
    // /Users/E30338/Repos/crab/build_crab2/ldd/src/ldd/cudd-2.4.2/cudd/cuddTable.c,
    // line 1143.

    using boxes_domain_t = boxes_domain<number_t, varname_t>;
    using var_t = typename boxes_domain_t::variable_t;

    boxes_domain_t s1;
    s1 += linear_constraint_t(var_t(x) >= number_t(0));
    s1 += linear_constraint_t(var_t(x) <= number_t(1, 4));
    s1 += linear_constraint_t(var_t(y) >= number_t(1, 4));
    s1 += linear_constraint_t(var_t(y) <= number_t(1, 2));
    boxes_domain_t s2;
    s2 += linear_constraint_t(var_t(x) == number_t(0));
    s2 += linear_constraint_t(var_t(y) == number_t(1, 2));

    boxes_domain_t s3 = s2 | s1;
    ///////////////////
    boxes_domain_t s4;
    s4 += linear_constraint_t(var_t(x) >= number_t(1, 4));
    s4 += linear_constraint_t(var_t(x) <= number_t(1));
    s4 += linear_constraint_t(var_t(y) >= number_t(1, 8));
    s4 += linear_constraint_t(var_t(y) <= number_t(1, 4));

    boxes_domain_t s5 = s4 | s2;
    //////
    crab::outs() << "Widening of \n" << s3;
    crab::outs() << "and\n" << s5 << " = \n";
    boxes_domain_t s6 = s3 || s5;
    crab::outs() << s6 << "\n";
  }
#endif

  // Ldd_NodeSanityCheck(man, &(*s4));

  // FIXME: seg fault here
  // destroy_ldd_man(man);
#endif
  return 0;
}
