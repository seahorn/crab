/*******************************************************************************
 * Types and functions to represent program expressions.
 ******************************************************************************/

#ifndef EXPR_H__
#define EXPR_H__

#include <iostream>
#include <crab/common/bignums.hpp>

/////
// Type to represent values
/////
// Ideally, we would like to use z_number_t to represent values with
// unlimited precision. However, val_t needs to be casted to (void *)
// so it cannot be bigger than the size of a pointer. 
typedef signed long val_t;
// VAL_MAX and VAL_MIN must be the maximum and minimum possible values
// of val_t
#define VAL_MAX LONG_MAX
#define VAL_MIN (-LONG_MAX) // We need -VAL_MAX = VAL_MIN.
/////


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
  val_t val;
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
  val_t konst;
} dexpr;

typedef struct {
  int r_from;
  int r_to;
} rmap;

exp_t exp_var(int var);
exp_t exp_const(val_t val);
exp_t exp_add(exp_t x, exp_t y);
exp_t exp_sub(exp_t x, exp_t y);
exp_t exp_mul(exp_t x, exp_t y);
exp_t exp_div(exp_t x, exp_t y);

void exp_dealloc(exp_t x);
void exp_collect_ptr(exp_t* x);

void exp_print_to(std::ostream &o, exp_t e);

uterm uvar(int v);
uterm uconst(val_t c);
ucon mk_ucon(uterm x, UConT kind, uterm y);

#endif
