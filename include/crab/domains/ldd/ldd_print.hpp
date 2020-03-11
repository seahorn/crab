#pragma once

#include "lddInt.h"
#include "util.h"

// ugly but it disables cudd macro that might conflict with boost::lexical_cast
#ifdef fail
#undef fail
#endif

namespace crab {
namespace domains {
namespace ldd {

static void Ldd_print_minterm_smtlib_aux(LddManager *ldd, LddNode *node,
                                         int *list, char **vnames);

/**
 * Prints disjoint sum of product cover for the function in smtlibv1
 * format. Based on Ldd_PrintMinterm.
 */
static inline int Ldd_PrintMintermSmtLibv1(LddManager *ldd, LddNode *node,
                                           char **vnames) {
  int i, *list;
  // DdNode *zero;
  // zero = Cudd_Not(DD_ONE(CUDD));

  list = ALLOC(int, CUDD->size);
  if (list == NULL) {
    CUDD->errorCode = CUDD_MEMORY_OUT;
    return 0;
  }

  for (i = 0; i < CUDD->size; i++)
    list[i] = 2;
  Ldd_print_minterm_smtlib_aux(ldd, node, list, vnames);
  FREE(list);
  return (1);
}

static void Ldd_print_minterm_smtlib_aux(LddManager *ldd, LddNode *n, int *list,
                                         char **vnames) {
  DdNode *N, *Nv, *Nnv;
  int i, v, index, p;

  /**
   * the latest negative constraint to be printed.
   * is NULL if there isn't one
   */
  lincons_t negc = NULL;

  N = Cudd_Regular(n);

  if (cuddIsConstant(N)) {
    /* n == N here implies that n is one */
    if (n == N) {
      /* for each level */
      for (i = 0; i < CUDD->size; i++) {
        /* let p be the index of level i */
        p = CUDD->invperm[i];
        /* let v be the value of p */
        v = list[p];

        /* skip don't care */
        if (v == 2)
          continue;

        if (v == 0 && ldd->ddVars[p] != NULL) {
          lincons_t c;
          c = THEORY->negate_cons(ldd->ddVars[p]);

          if (negc != NULL) {
            /* print negative constraint if it is not
               implied by c
            */
            if (!THEORY->is_stronger_cons(c, negc)) {
              THEORY->print_lincons_smtlibv1(CUDD->out, negc, vnames);
              fprintf(CUDD->out, " ");
            }
            THEORY->destroy_lincons(negc);
          }

          /* store the current constraint to be printed later */
          negc = c;
          continue;
        }

        /* if there is a negative constraint waiting to be printed,
           print it now
        */
        if (negc != NULL) {
          THEORY->print_lincons_smtlibv1(CUDD->out, negc, vnames);
          fprintf(CUDD->out, " ");
          THEORY->destroy_lincons(negc);
          negc = NULL;
        }

        /* if v is not a don't care but p does not correspond to
         * a constraint, print it as a Boolean variable */
        if (v != 2 && ldd->ddVars[p] == NULL)
          fprintf(stderr, "%sb%d", (v == 0 ? "!" : " "), p);
        /* v is true */
        else if (v == 1) {
          THEORY->print_lincons_smtlibv1(CUDD->out, ldd->ddVars[p], vnames);
          fprintf(CUDD->out, " ");
        }
      }

      /* if there is a constraint waiting to be printed, do it now */
      if (negc != NULL) {
        THEORY->print_lincons_smtlibv1(CUDD->out, negc, vnames);
        THEORY->destroy_lincons(negc);
        negc = NULL;
      }
      fprintf(CUDD->out, "\n");
    }
  } else {
    Nv = Cudd_NotCond(cuddT(N), N != n);
    Nnv = Cudd_NotCond(cuddE(N), N != n);
    index = N->index;
    list[index] = 0;
    Ldd_print_minterm_smtlib_aux(ldd, Nnv, list, vnames);
    list[index] = 1;
    Ldd_print_minterm_smtlib_aux(ldd, Nv, list, vnames);
    list[index] = 2;
  }
  return;
}
} // namespace ldd
} // namespace domains
} // namespace crab
