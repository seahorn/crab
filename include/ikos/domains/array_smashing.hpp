/*******************************************************************************
 * Simple array smashing domain
 ******************************************************************************/

#ifndef IKOS_ARRAY_SMASHING_HPP
#define IKOS_ARRAY_SMASHING_HPP

#include <ikos/common/types.hpp>
#include <ikos/domains/uninitialized_domain.hpp>
#include <ikos/domains/numerical_domains_api.hpp>

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

  void undefined (VariableName arr)
  {
    // We need to mark an array variable as uninitialized first so
    // when we do the first store we can at least perform a strong
    // update. Otherwise, if we perform a weak update we will get
    // always top since variables are top by default.
    _uninit.set (arr, uninitializedValue::uninitialized ());
  }

  void load (VariableName lhs, VariableName arr, VariableName idx)
  {
    // FIXME: if a relational domain we can be more precise by making
    // a copy of arr and store it into lhs but *without* establishing
    // any relationship between lhs and arr (i.e., without using apply).
    _inv.set (lhs, _inv[arr]);
  }


  void store (VariableName arr_out, VariableName arr_in, 
              VariableName /*idx*/, linear_expression_t val,
              bool is_singleton)
  {

    if (is_singleton || _uninit [arr_in].is_uninitialized ())
    {
      strong_update (arr_in, val);
      _uninit.set (arr_in, uninitializedValue::initialized ());
    }
    else
    {
      weak_update (arr_in, val);
    }

    if (!(arr_out == arr_in))
    {
      _inv.assign (arr_out, linear_expression_t (arr_in));
      _uninit.set (arr_out, _uninit [arr_in]);
    }    
  }
    
  ostream& write(ostream& o) 
  {
    o << _inv;
    return o;
  }
  
}; // end array_smashing

} // namespace ikos

#endif 
