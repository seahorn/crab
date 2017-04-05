#include <crab/config.h>

#include "../common.hpp"
#ifdef HAVE_LDD
#include <crab/domains/ldd/ldd.hpp>
#endif 
#include <boost/optional.hpp>

using namespace std;
using namespace ikos;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
#ifdef HAVE_LDD
using namespace crab::domains::ldd;
#endif 

#define RATIONALS

#ifdef RATIONALS
typedef q_number number_t;
#else
typedef z_number number_t;
#endif 

typedef linear_constraint<number_t, varname_t> linear_constraint_t;
typedef linear_expression<number_t, varname_t> linear_expression_t;
typedef interval <number_t> interval_t;

#ifdef HAVE_LDD

LddManager* create_ldd_man (size_t num_vars) {
  DdManager* cudd = Cudd_Init (0, 0, CUDD_UNIQUE_SLOTS, 127, 0);
  #ifndef RATIONALS
  theory_t * theory = tvpi_create_boxz_theory (num_vars);
  #else
  theory_t * theory = tvpi_create_box_theory (num_vars);
  #endif 
  
  LddManager* ldd = Ldd_Init (cudd, theory);

  const bool dvo = true;
  if(dvo) Cudd_AutodynEnable (cudd, CUDD_REORDER_GROUP_SIFT);
  return ldd;
}

void destroy_ldd_man (LddManager* ldd) {
  DdManager *cudd = NULL;
  theory_t *theory = NULL;
  if (ldd) {
    cudd = Ldd_GetCudd (ldd);
    theory = Ldd_GetTheory (ldd);
    Ldd_Quit(ldd);
  }
  
  if (theory) tvpi_destroy_theory(theory);
  if (cudd) Cudd_Quit(cudd);
}

void dump (LddNodePtr val) {
  LddManager *ldd = getLddManager (val);
  DdManager *cudd = Ldd_GetCudd (ldd);
  FILE *fp = Cudd_ReadStdout(cudd);
  Cudd_SetStdout(cudd,stderr);
  if (val.get () == Ldd_GetTrue (ldd)) crab::outs () << "true\n";
  else if (val.get () == Ldd_GetFalse (ldd)) crab::outs () << "false\n";	
  else Ldd_PrintMinterm(ldd,val.get ());
  Cudd_SetStdout(cudd,fp);      
}

LddNodePtr top (LddManager *ldd) 
{
  return lddPtr (ldd, Ldd_GetTrue (ldd));
}

LddNodePtr bot (LddManager *ldd)
{
  return lddPtr (ldd, Ldd_GetFalse (ldd));
}

bool isBot (LddManager * ldd, LddNodePtr v)
{
  return &*v == Ldd_GetFalse (ldd);
}

bool isTop (LddManager * ldd, LddNodePtr v)
{
  return &*v == Ldd_GetTrue (ldd);
}

LddNodePtr approx (LddManager *ldd, LddNodePtr v) 
{
  return lddPtr (ldd, Ldd_TermMinmaxApprox(ldd, &*v));
}

LddNodePtr meet (LddManager *ldd, LddNodePtr n1, LddNodePtr n2)
{
  return lddPtr (ldd, Ldd_And (ldd, &*n1, &*n2));
}

LddNodePtr join (LddManager *ldd, LddNodePtr n1, LddNodePtr n2)
{
  return approx (ldd, lddPtr (ldd, Ldd_Or (ldd, &*n1, &*n2)));
}

bool isLeq (LddManager *ldd, LddNodePtr n1, LddNodePtr n2)
{
  return Ldd_TermLeq(ldd, &(*n1), &(*n2));
}

LddNodePtr widen (LddManager *ldd, LddNodePtr n1, LddNodePtr n2)
{
  return lddPtr (ldd, Ldd_IntervalWiden (ldd, &*n1, &*n2));
}
 
int getVarId (varname_t v) {
  return v.index ();
}

/** return term for variable v, neg for negation of variable */
linterm_t termForVal(LddManager* ldd, varname_t v, bool neg = false) 
{
  int varId = getVarId (v);
  int sgn = neg ? -1 : 1;
  linterm_t term = Ldd_GetTheory (ldd)->create_linterm_sparse_si (&varId, &sgn, 1);
  return term; 
}


/** v := k, where k is a constant */
LddNodePtr assign (LddManager *ldd, LddNodePtr n, varname_t v, number_t k) {
  if (isBot (ldd, n)) return n;

  #ifndef RATIONALS
  mpq_class kk ((mpz_class) k);
  #else
  mpq_class kk = ((mpq_class) k);
  #endif 
  linterm_t t = termForVal (ldd, v);
  constant_t kkk = (constant_t) tvpi_create_cst (kk.get_mpq_t ());
  LddNodePtr newVal = 
      lddPtr(ldd, Ldd_TermReplace(ldd, &(*n), t, NULL, NULL, kkk, kkk));
  Ldd_GetTheory (ldd)->destroy_cst(kkk);
  
  assert (!isBot (ldd, newVal));
  return newVal;
}

/** v := k, where k is an interval */
LddNodePtr assign (LddManager* ldd, LddNodePtr n, varname_t v, interval_t ival) {
                   
  if (isBot (ldd, n)) return n;
  
  constant_t kmin = NULL, kmax = NULL;
  
  if (boost::optional <number_t> l = ival.lb ().number ()) {
    #ifndef RATIONALS    
    mpq_class val ((mpz_class) (*l));
    #else
    mpq_class val ((mpq_class) (*l));    
    #endif 
    kmin = (constant_t)tvpi_create_cst (val.get_mpq_t ());
  }

  if (boost::optional <number_t> u = ival.ub ().number ()) {
    #ifndef RATIONALS        
    mpq_class val ((mpz_class) (*u));
    #else
    mpq_class val ((mpq_class) (*u));    
    #endif 
    kmax = (constant_t)tvpi_create_cst (val.get_mpq_t ());
  }
       
  linterm_t t = termForVal(ldd, v);
  LddNodePtr newVal = 
      lddPtr(ldd, Ldd_TermReplace(ldd, &(*n), t, NULL, NULL, kmin, kmax));
  
  if (kmin) Ldd_GetTheory (ldd)->destroy_cst(kmin);
  if (kmax) Ldd_GetTheory (ldd)->destroy_cst(kmax);
  
  assert (!isBot (ldd, newVal));
  return newVal;
}

/** v := a * u + k, where a, k are constants and u variable */
LddNodePtr apply (LddManager* ldd,
                  const LddNodePtr n, 
                  varname_t v, varname_t u,
                  number_t a, number_t k)
{
  if (isTop (ldd, n) || isBot (ldd, n)) return n;
  
  linterm_t t = termForVal (ldd, v);
  linterm_t r = termForVal (ldd, u);

  #ifndef RATIONALS          
  mpq_class aa ((mpz_class) a); 
  mpq_class kk ((mpz_class) k);
  #else
  mpq_class aa ((mpq_class) a); 
  mpq_class kk ((mpq_class) k); 
  #endif   
  
  constant_t aaa = (constant_t) tvpi_create_cst (aa.get_mpq_t ());
  constant_t kkk = (constant_t) tvpi_create_cst (kk.get_mpq_t ());        
  LddNodePtr newVal = 
      lddPtr(ldd, Ldd_TermReplace(ldd, &(*n), t, r, aaa, kkk, kkk));
  Ldd_GetTheory (ldd)->destroy_cst(aaa);
  Ldd_GetTheory (ldd)->destroy_cst(kkk);
  
  assert (!isBot (ldd, newVal));
  return newVal;
}

/** v := u */
LddNodePtr assign (LddManager* ldd, const LddNodePtr n, 
                   varname_t v, varname_t u)
{
  return apply (ldd, n, v, u, 1, 0);      
}
#endif     
   
int main (int argc, char** argv )
{
#ifdef HAVE_LDD
  SET_TEST_OPTIONS(argc,argv)

  LddManager* man = create_ldd_man (3000);
  Ldd_SanityCheck (man);
	  
  variable_factory_t vfac;
  varname_t x = vfac["x"];
  varname_t y = vfac["y"];

  {
    LddNodePtr s1  = top (man);
    LddNodePtr s11 = assign (man, s1, x, 5.5);
    LddNodePtr s111 = assign (man, s11, y, 8.3);
    
    LddNodePtr s2  = top (man);
    LddNodePtr s22 = assign (man, s2, x, 5.5);
    LddNodePtr s222 = assign (man, s22, y, 12.2);
    
    dump (s111);
    dump (s222);
    crab::outs () << "Join \n";
    LddNodePtr s3 = join (man, s111, s222);
    dump (s3);
    crab::outs () << "Widening \n";  
    LddNodePtr s4 = widen (man, s111, s222);
    dump (s4);
  }

  // Ldd_NodeSanityCheck (man, &(*s4));
  
  // FIXME: seg fault here
  // destroy_ldd_man (man);
#endif 
  return 0;
}
