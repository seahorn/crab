/*******************************************************************************
 * Types and functions to represent program expressions.
 ******************************************************************************/

#ifndef EXPR_H__
#define EXPR_H__

#include <iostream>

// Data structures for representing expression trees.
typedef enum { E_VAR, E_CONST, E_ADD, E_SUB, E_MUL, E_DIV } ExpT; 
typedef enum { B_LB, B_UB } BoundT;
typedef enum { L_LEQ, L_LT, L_GEQ, L_GT, L_EQ, L_DISEQ } LinConT;

typedef enum { U_VAR, U_CONST } UTermT;
typedef enum { U_LEQ, U_LT, U_EQ, U_DIS } UConT;

typedef enum { D_LB, D_UB, D_DIFF } DExprT;

typedef struct {
  ExpT kind;
  int refs;
  void* args[2];
} exp_data;
typedef exp_data* exp_t;

// Either a constant or
// variable with unit-coefficient.
typedef struct {
  int kind;
  int val;
} uterm;

// Binary constraint with
// unit coefficients.
typedef struct {
  UConT kind;
  uterm args[2]; 
} ucon;

typedef struct {
  DExprT kind;
  int args[2];
  int konst;
} dexpr;

typedef struct {
  int r_from;
  int r_to;
} rmap;

exp_t exp_var(int var);
exp_t exp_const(int val);
exp_t exp_add(exp_t x, exp_t y);
exp_t exp_sub(exp_t x, exp_t y);
exp_t exp_mul(exp_t x, exp_t y);
exp_t exp_div(exp_t x, exp_t y);

void exp_dealloc(exp_t x);
void exp_collect_ptr(exp_t* x);

// void linexp_dealloc(linexp lexp);

void exp_print_to(std::ostream &o, exp_t e);

uterm uvar(int v);
uterm uconst(int c);
ucon mk_ucon(uterm x, UConT kind, uterm y);

#endif
