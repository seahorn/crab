#pragma once

#include <crab/domains/numerical_domains_api.hpp>
#include <crab/domains/bitwise_operators_api.hpp>
#include <crab/domains/cast_operators_api.hpp>
#include <crab/domains/array_operators_api.hpp>
#include <crab/domains/pointer_operators_api.hpp>
#include <crab/domains/boolean_operators_api.hpp>

#include <crab/common/types.hpp>

namespace crab {

  namespace domains {
    
    /** 
     * All abstract domains must inherit this class and expose
     * publicly all the public typedef's.
     **/
    template<typename N, typename V, typename Dom>
    class abstract_domain:
      public ikos::writeable,
      public ikos::numerical_domain<N,V>,
      public backward_numerical_domain<N,V,Dom>,
      public ikos::bitwise_operators<N,V>,
      public int_cast_operators<N,V>,
      public array_operators<N,V>,
      public pointer_operators<N,V>,
      public boolean_operators<N,V>,
      public backward_boolean_operators<N,V,Dom> {
      
    public:
      
      typedef ikos::linear_expression<N,V> linear_expression_t;
      typedef ikos::linear_constraint<N,V> linear_constraint_t;
      typedef ikos::linear_constraint_system<N,V> linear_constraint_system_t;
      typedef ikos::disjunctive_linear_constraint_system<N,V>
      disjunctive_linear_constraint_system_t;      
      typedef ikos::variable<N,V> variable_t;
      typedef N number_t;
      typedef V varname_t;
      typedef std::vector<variable_t> variable_vector_t;
            
      abstract_domain(): ikos::writeable() {}
      
      virtual ~abstract_domain() {};

      virtual linear_constraint_system_t to_linear_constraint_system() = 0;
      
      virtual disjunctive_linear_constraint_system_t
      to_disjunctive_linear_constraint_system() = 0;
      
      virtual void rename(const variable_vector_t &from, const variable_vector_t &to) {
	// If you see this error message then implement this operation
	// in the corresponding abstract domain.	
	CRAB_ERROR("rename operation not implemented");
      }

      // Normalize the abstract domain if such notion exists.
      virtual void normalize() { }
      
      // Forget variables form the abstract domain
      virtual void forget(const variable_vector_t& variables) = 0;

      // Project the abstract domain onto variables (dual to forget)
      virtual void project(const variable_vector_t& variables) = 0;

      // Make a new copy of var without relating var with new_var
      virtual void expand(variable_t var, variable_t new_var) = 0;            
      
    };
    
  }
}
