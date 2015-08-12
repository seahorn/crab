/*******************************************************************************
 * Array smashing domain
 * 
 * FIXME: assume all array accesses are aligned wrt to the size of the
 * array element (e.g., if the size of the array element is 4 bytes
 * then all array accesses must be multiple of 4). Note that this
 * assumption does not hold in real programs!
 ******************************************************************************/

#ifndef IKOS_ARRAY_SMASHING_HPP
#define IKOS_ARRAY_SMASHING_HPP

#include <ikos/common/types.hpp>
#include <ikos/domains/uninitialized_domain.hpp>
#include <ikos/algorithms/patricia_trees.hpp>
#include <ikos/domains/numerical_domains_api.hpp>
#include <ikos/domains/domain_traits_impl.hpp>
#include "boost/range/algorithm/set_algorithm.hpp"
//#define DEBUG
#include <ikos/common/dbg.hpp>

namespace ikos {

using namespace std;
using namespace boost;

//! Abstract domain to reason about summarized variables. All array
//  elements are `smashed` into a single cell.
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
  typedef interval <Number> interval_t;

 private:
  typedef mergeable_map<VariableName, VariableName> sz_map_t;
  typedef set<VariableName> varname_set_t;
  typedef bound <Number> bound_t; 
  //! scalar and array variables        
  NumDomain _inv; 

  array_smashing (NumDomain inv): 
      ikos::writeable (), 
      _inv (inv) { }

  void strong_update (VariableName lhs, linear_expression_t rhs ) {
    _inv.assign (lhs, rhs);
  }

  void weak_update (VariableName lhs, linear_expression_t rhs) {
    NumDomain other (_inv);
    other.assign (lhs, rhs);
    _inv = _inv | other;
  }

  template <typename Set>
  Set set_union (Set &s1, Set &s2)
  {
    Set s3;
    boost::set_union (s1, s2, std::inserter (s3, s3.end ()));
    return s3;
  }

  template <typename Set, typename Elem>
  void set_insert (Set &s1, Elem e)
  {
    Set s2, s3;
    s2.insert (e);
    boost::set_union (s1, s2, std::inserter (s3, s3.end ()));
    std::swap (s3, s1);
  }

  template <typename Set>
  void set_difference (Set &s1, Set &s2)
  {
    Set s3;
    boost::set_difference (s1, s2, std::inserter (s3, s3.end ()));
    std::swap (s3, s1);
  }
  
  template<typename Set>
  NumDomain project (NumDomain inv, Set all, Set vs) 
  {
    set_difference (all, vs);
    domain_traits::forget (inv, all.begin (), all.end ()); 
    return inv;
  }

public:
  
  array_smashing(): ikos::writeable(), _inv (NumDomain::top()) { }    
      
  static array_smashing_t top() { 
    return array_smashing (NumDomain::top ()); 
  }
  
  static array_smashing_t bottom() {
    return array_smashing (NumDomain::bottom ());
  }
      
  array_smashing (const array_smashing_t& other): 
    ikos::writeable(), 
    _inv (other._inv) { }
  
  array_smashing_t& operator=(array_smashing_t other) {
    _inv = other._inv;
    return *this;
  }
  
  bool is_bottom() { 
    return (_inv.is_bottom ());
  }

  bool is_top() { 
    return (_inv.is_top());
  }
  
  bool operator<=(array_smashing_t other) {
    return (_inv <= other._inv);
  }

  array_smashing_t operator|(array_smashing_t other) 
  {
    return array_smashing_t (_inv | other._inv);
  }
  
  array_smashing_t operator&(array_smashing_t other) {
    return array_smashing_t (_inv & other._inv);
  }
  
  array_smashing_t operator||(array_smashing_t other) {
    return array_smashing_t (_inv || other._inv);
  }
  
  array_smashing_t operator&& (array_smashing_t other) {
    return array_smashing_t (_inv && other._inv);
  }
  
  void operator-=(VariableName var) {
    _inv -= var;
  }
  
  void operator += (linear_constraint_system_t csts) {
    _inv += csts;
  }

  void assign (VariableName x, linear_expression_t e) {
    _inv.assign (x, e);
  }

  void apply (operation_t op, VariableName x, VariableName y, Number z) {
    _inv.apply (op, x, y, z);
  }

  void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
    _inv.apply (op, x, y, z);
  }

  void apply(operation_t op, VariableName x, Number k) {
    _inv.apply (op, x, k);
  }

  void array_init (VariableName a, 
                   const vector<ikos::z_number>& values) {
    if (values.empty ()) return;
    
    interval_t init = interval_t::bottom ();
    for (auto const &v: values) {
      // assume automatic conversion from z_number to bound_t
      init = init | interval_t (bound_t (v)); 
    }
    _inv.set (a, init);
  }

  // All the array elements are initialized to val
  void assume_array (VariableName a, interval_t val) {
    _inv.set (a, val);
  }

  void load (VariableName lhs, VariableName a, 
             VariableName /*i*/, z_number /*n_bytes*/) {

    // We need to be careful when assigning a summarized variable a
    // into a non-summarized variable lhs. Simply _inv.assign (lhs,
    // a) is not sound.
    /* ask for a temp var */
    VariableName a_prime = a.getVarFactory().get(); 
    domain_traits::expand (_inv, a, a_prime);
    _inv.assign (lhs, linear_expression_t (a_prime));
    _inv -= a_prime; 
  }


  void store (VariableName a, VariableName /*i*/,
              linear_expression_t val, z_number /*n_bytes*/,
              bool is_singleton) {

    if (is_singleton)
      strong_update (a, val);
    else 
      weak_update (a, val);

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
                 VariableName a, 
                 const vector<ikos::z_number> &values) {
  inv.array_init (a, values);
}

template <typename BaseDomain, typename VariableName, typename Number>
void assume_array (array_smashing<BaseDomain,Number,VariableName>& inv, 
                   VariableName a, Number val) {
  inv.assume_array (a, interval<Number> (bound <Number> (val)));
}

template <typename BaseDomain, typename VariableName, typename Number>
void assume_array (array_smashing<BaseDomain,Number,VariableName>& inv, 
                   VariableName a, interval<Number> val) {
  inv.assume_array (a, val);
}

template <typename BaseDomain, typename VariableName, typename Number>
void array_load (array_smashing<BaseDomain, Number, VariableName>& inv, 
                 VariableName lhs, VariableName a, 
                 VariableName i, z_number n_bytes) {
  inv.load (lhs, a, i, n_bytes);
}

template <typename BaseDomain, typename VariableName, typename Number>
void array_store (array_smashing<BaseDomain, Number, VariableName>& inv, 
                  VariableName a, VariableName i,
                  typename BaseDomain::linear_expression_t val,
                  z_number n_bytes, bool is_singleton) {
  inv.store (a, i, val, n_bytes, is_singleton);
}
} // namespace domain_traits
} // namespace ikos

#endif 
