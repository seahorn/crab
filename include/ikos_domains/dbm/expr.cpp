#include <ikos_domains/dbm/expr.h>

exp_t exp_var(int var)
{
  exp_t r = new exp_data;
  r->kind = E_VAR;
  r->args[0] = (void*) var;
  r->refs = 0;
  return r;
}

exp_t exp_const(int val)
{
  exp_t r = new exp_data;
  r->kind = E_CONST;
  r->args[0] = (void*) val;
  r->refs = 0;
  return r;
}

exp_t exp_binop(ExpT kind, exp_t x, exp_t y)
{
  exp_t r = new exp_data;
  r->kind = kind;
  r->args[0] = (void*) x;
  r->args[1] = (void*) y;
  r->refs = 0;
  x->refs++;
  y->refs++;
  return r;
}

void exp_print_to(std::ostream &o, exp_t e){
  if (e->kind == E_CONST || e->kind == E_VAR)
    std::cout << (size_t) e->args[0];
  else{
    std::cout << "(";
    exp_print_to(o, (exp_t) e->args[0]);
    switch (e->kind){
      case E_ADD:
        std::cout << " + ";
      break;
      case E_SUB:
        std::cout << " - ";
      break;
      case E_MUL:
        std::cout << " * ";
      break;
      case E_DIV:
        std::cout << " / ";
      break;
      default: ;;
    }
    exp_print_to(o, (exp_t) e->args[1]);        
    std::cout << ")";
  }
}

exp_t exp_add(exp_t x, exp_t y) { return exp_binop(E_ADD, x, y); }
exp_t exp_sub(exp_t x, exp_t y) { return exp_binop(E_SUB, x, y); }
exp_t exp_mul(exp_t x, exp_t y) { return exp_binop(E_MUL, x, y); }
exp_t exp_div(exp_t x, exp_t y) { return exp_binop(E_DIV, x, y); }

void exp_dealloc(exp_t x)
{
  if(!x)
    return;

  switch(x->kind)
  {
    // Terminals
    case E_VAR:
    case E_CONST:
      break;
    default:
      exp_dealloc((exp_t) (x->args[0]));
      exp_dealloc((exp_t) (x->args[1]));
      break;
  }
  delete x;
}

void exp_collect(exp_t x)
{
  if(!x)
    return;
  if(x->refs)
    return;
  switch(x->kind)
  {
    case E_VAR:
    case E_CONST:
      break;
    default:
      {
        exp_t l = (exp_t) (x->args[0]);
        exp_t r = (exp_t) (x->args[1]);
        l->refs--;
        r->refs--;
        exp_collect(l);
        exp_collect(r);
      }
      break;
  }
  delete x;
}
void exp_collect_ptr(exp_t* xptr) {
  exp_collect(*xptr);
}
  

/*
void linexp_dealloc(linexp lexp)
{
  delete lexp;
}
*/

uterm uvar(int v) { uterm t = { U_VAR, v }; return t; }
uterm uconst(int c) { uterm t = { U_CONST, c}; return t; }
ucon mk_ucon(uterm x, UConT kind, uterm y) {
  ucon c = { kind };
  c.args[0] = x;
  c.args[1] = y;
  return c;
}
