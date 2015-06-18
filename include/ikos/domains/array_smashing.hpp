/*******************************************************************************
 * Simple array smashing domain
 ******************************************************************************/

#ifndef IKOS_ARRAY_SMASHING_HPP
#define IKOS_ARRAY_SMASHING_HPP

#include <ikos/common/types.hpp>
#include <ikos/domains/uninitialized_domain.hpp>
#include <ikos/domains/numerical_domains_api.hpp>
#include <ikos/domains/domain_traits_impl.hpp>

namespace ikos {

using namespace std;
using namespace boost;
 
template<typename NumDomain, typename Number, typename VariableName>
class array_smashing: 
      public ikos::writeable, 
      public numerical_domain< Number, VariableName>
{

 public:

  typedef typename NumDomain::linear_constraint_t linear_constraint_t;
  typedef typename NumDomain::linear_constraint_system_t linear_constraint_system_t;
  typedef typename NumDomain::linear_expression_t linear_expression_t;
  typedef typename NumDomain::variable_t variable_t;
  
  typedef array_smashing <NumDomain, Number, VariableName> array_smashing_t;

 private:

  typedef uninitialized_domain <VariableName> UninitDomain;

  NumDomain _inv;        
  UninitDomain _uninit; //! only for array variables

  array_smashing (NumDomain inv, UninitDomain uninit): 
      ikos::writeable (), _inv (inv), _uninit (uninit) { }

  void strong_update (VariableName lhs, linear_expression_t rhs )
  {
    _inv.assign (lhs, rhs);
  }

  void weak_update (VariableName lhs, linear_expression_t rhs )
  {
    NumDomain other (_inv);
    other.assign (lhs, rhs);
    _inv = _inv | other;
  }

public:
  
  array_smashing(): 
      ikos::writeable(), 
      _inv (NumDomain::top()), _uninit (UninitDomain::top ())  { }    

  static array_smashing_t top() 
  {
    return array_smashing (NumDomain::top (), UninitDomain::top ());
  }
  
  static array_smashing_t bottom() 
  {
    return array_smashing (NumDomain::bottom (), UninitDomain::bottom ());
  }
      
  array_smashing (const array_smashing_t& other): ikos::writeable(), 
                                                  _inv (other._inv), 
                                                  _uninit (other._uninit)
  { }
  
  array_smashing_t& operator=(array_smashing_t other) 
  {
    if (this != &other)
    {
      _inv = other._inv;
      _uninit = other._uninit;
    }
    return *this;
  }
  
  bool is_bottom() { return _inv.is_bottom () || _uninit.is_bottom (); }

  bool is_top() { return _inv.is_top() && _uninit.is_top (); }
  
  bool operator<=(array_smashing_t other) 
  {
    return (_inv <= other._inv && _uninit <= other._uninit);
  }
  
  array_smashing_t operator|(array_smashing_t other) 
  {
    // Reduction between the numerical and uninitialized domain.
    // 
    // Before we perform the join if an array variable is
    // uninitialized in one of the environments we take the value from
    // the other one.
    if (!_uninit.is_top () && !_uninit.is_bottom ())
    {
      for (auto p : _uninit)
      {
        if (p.second.is_uninitialized () && 
            !other._uninit [p.first].is_uninitialized ())
        {// FIXME: this does not keep relational invariants
          _inv.set (p.first, other._inv [p.first]);
        }
      }
      for (auto p : other._uninit)
      {
        if (!p.second.is_uninitialized () && 
            _uninit [p.first].is_uninitialized ())
        {// FIXME: this does not keep relational invariants
          other._inv.set (p.first, _inv[p.first]);
        }
      }
    }

    return array_smashing_t (_inv | other._inv, _uninit | other._uninit);
  }
  
  array_smashing_t operator&(array_smashing_t other) 
  {
    return array_smashing_t (_inv & other._inv, _uninit & other._uninit);
  }
  
  array_smashing_t operator||(array_smashing_t other) 
  {
    return array_smashing_t (_inv || other._inv, _uninit || other._uninit);
  }
  
  array_smashing_t operator&& (array_smashing_t other) 
  {
    return array_smashing_t (_inv && other._inv, _uninit && other._uninit);
  }
  
  void operator-=(VariableName var)
  {
    _inv -= var;
    _uninit -= var;
  }
  
  void operator += (linear_constraint_system_t csts) 
  {
    _inv += csts;
  }

  void assign (VariableName x, linear_expression_t e) 
  {
    _inv.assign (x, e);
    auto v = e.get_variable ();
    if (v) _uninit.assign (x, (*v).name ());      
  }

  void apply (operation_t op, VariableName x, VariableName y, Number z) 
  {
    _inv.apply (op, x, y, z);
  }

  void apply(operation_t op, VariableName x, VariableName y, VariableName z) 
  {
    _inv.apply (op, x, y, z);
  }

  void apply(operation_t op, VariableName x, Number k) 
  {
    _inv.apply (op, x, k);
  }

  void array_init (VariableName arr)
  {
    // We need to mark an array variable as uninitialized first so
    // when we do the first store we can perform a strong
    // update. Otherwise, if we perform a weak update we will get
    // always top since variables are top by default.
    _uninit.set (arr, uninitializedValue::uninitialized ());
  }

  void load (VariableName lhs, VariableName arr, VariableName idx)
  {
    // We need to be careful when assigning a summarized variable arr
    // into a non-summarized variable lhs. 
    // Simply do _inv.assign (lhs, arr) is wrong.

    VariableName arr_prime = arr.getVarFactory().get(); /* ask for a temp var */
    domain_traits::expand (_inv, arr, arr_prime);
    _inv.assign (lhs, linear_expression_t (arr_prime));
    _inv -= arr_prime; 
  }


  void store (VariableName arr, VariableName /*idx*/,
              linear_expression_t val, bool is_singleton)
  {

    if (is_singleton || _uninit [arr].is_uninitialized ())
    {
      strong_update (arr, val);
      _uninit.set (arr, uninitializedValue::initialized ());
    }
    else
    {
      weak_update (arr, val);
    }
  }

  linear_constraint_system_t to_linear_constraint_system (){
    return _inv.to_linear_constraint_system ();
  }
    
  ostream& write(ostream& o) 
  {
    o << _inv;
    return o;
  }

  const char* getDomainName () const {return "Array smashing";}  

}; // end array_smashing

namespace domain_traits
{

template <typename BaseDomain, typename VariableName, typename Number>
void array_init (array_smashing<BaseDomain,Number,VariableName>& inv, 
                 VariableName arr) {
  inv.array_init (arr);
}

template <typename BaseDomain, typename VariableName, typename Number>
void array_load (array_smashing<BaseDomain, Number, VariableName>& inv, 
                 VariableName lhs, VariableName arr, VariableName idx) {
   inv.load (lhs, arr, idx);
}

template <typename BaseDomain, typename VariableName, typename Number>
void array_store (array_smashing<BaseDomain, Number, VariableName>& inv, 
                  VariableName arr, VariableName idx,
                  typename BaseDomain::linear_expression_t val,
                  bool is_singleton) {
   inv.store (arr, idx, val, is_singleton);
}
} // namespace domain_traits
} // namespace ikos

#endif 
