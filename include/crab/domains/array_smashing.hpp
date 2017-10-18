/*******************************************************************************
 * Array smashing domain
 * 
 * FIXME: assume all array accesses are aligned wrt to the size of the
 * array element (e.g., if the size of the array element is 4 bytes
 * then all array accesses must be multiple of 4). Note that this
 * assumption does not hold in real programs.
 ******************************************************************************/

#ifndef ARRAY_SMASHING_HPP
#define ARRAY_SMASHING_HPP

#include <crab/common/types.hpp>
#include <crab/common/debug.hpp>
#include <crab/common/stats.hpp>
#include <crab/domains/operators_api.hpp>
#include <crab/domains/domain_traits.hpp>

namespace crab {

   namespace domains {

      //! Abstract domain to reason about summarized variables. All array
      //  elements are `smashed` into a single cell.
      template<typename NumDomain>
      class array_smashing:
       public abstract_domain<typename NumDomain::number_t,
			      typename NumDomain::varname_t,
			      array_smashing<NumDomain> > {
	
      public:
	
	typedef typename NumDomain::number_t Number;
        typedef typename NumDomain::varname_t VariableName;

      private:
	
        typedef array_smashing<NumDomain> array_smashing_t;
	typedef abstract_domain<Number,VariableName,array_smashing_t> abstract_domain_t;
	  
      public:

        using typename abstract_domain_t::linear_expression_t;
        using typename abstract_domain_t::linear_constraint_t;
        using typename abstract_domain_t::linear_constraint_system_t;
        using typename abstract_domain_t::disjunctive_linear_constraint_system_t;	
        using typename abstract_domain_t::variable_t;
        using typename abstract_domain_t::number_t;
        using typename abstract_domain_t::varname_t;
	using typename abstract_domain_t::varname_vector_t;	
        typedef crab::pointer_constraint<VariableName> ptr_cst_t;
        typedef NumDomain content_domain_t;
        typedef interval <Number> interval_t;
        
       private:
        
        typedef bound <Number> bound_t; 
        
        //! scalar and summarized array variables        
        NumDomain _inv; 
        
        array_smashing (NumDomain inv): _inv (inv) { }
	  
        void strong_update (VariableName a, crab::variable_type a_ty,
			    linear_expression_t rhs ) {
	  if (a_ty == ARR_BOOL_TYPE) {
	    if (rhs.is_constant()) {
	      if (rhs.constant() >= Number(1))
		_inv.assign_bool_cst(a, linear_constraint_t::get_true());
	      else
		_inv.assign_bool_cst(a, linear_constraint_t::get_false());
	    } else if (auto rhs_v = rhs.get_variable ()) {
	      _inv.assign_bool_var (a, (*rhs_v).name(), false);
	    }
	  } else if (a_ty == ARR_INT_TYPE || a_ty == ARR_REAL_TYPE) {
            _inv.assign (a, rhs);
	  } else if (a_ty == ARR_PTR_TYPE) {
	    if (rhs.is_constant() && rhs.constant() == Number(0))
	      _inv.pointer_mk_null(a);
	    else if (auto rhs_v = rhs.get_variable())	     
	      _inv.pointer_assign (a, (*rhs_v).name(), Number(0));
	  }
        }
        
        void weak_update (VariableName a, crab::variable_type a_ty,
			  linear_expression_t rhs) {
          NumDomain other (_inv);

	  if (a_ty == ARR_BOOL_TYPE) {
	    if (rhs.is_constant()) {
	      if (rhs.constant() >= Number(1))
		other.assign_bool_cst(a, linear_constraint_t::get_true());
	      else
		other.assign_bool_cst(a, linear_constraint_t::get_false());
	    } else if (auto rhs_v = rhs.get_variable ()) {
	      other.assign_bool_var (a, (*rhs_v).name(), false);
	    }
	  } else if (a_ty == ARR_INT_TYPE || a_ty == ARR_REAL_TYPE) {
            other.assign (a, rhs);
	  } else if (a_ty == ARR_PTR_TYPE) {
	    if (rhs.is_constant() && rhs.constant() == Number(0))
	      other.pointer_mk_null(a);
	    else if (auto rhs_v = rhs.get_variable())	     
	      other.pointer_assign (a, (*rhs_v).name(), Number(0));
	  }
	  
          _inv = _inv | other;
        }
        
       public:
        
        array_smashing(): _inv (NumDomain::top()) { }    
        
        static array_smashing_t top() { 
          return array_smashing (NumDomain::top ()); 
        }
        
        static array_smashing_t bottom() {
          return array_smashing (NumDomain::bottom ());
        }
        
        array_smashing (const array_smashing_t& other): 
	  _inv (other._inv) { 
          crab::CrabStats::count (getDomainName() + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
        }
        
        array_smashing_t& operator=(const array_smashing_t& other) {
          crab::CrabStats::count (getDomainName() + ".count.copy");
          crab::ScopedCrabStats __st__(getDomainName() + ".copy");
          if (this != &other)
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

        void operator|=(array_smashing_t other) {
          _inv |= other._inv;
        }
        
        array_smashing_t operator|(array_smashing_t other) {
          return array_smashing_t (_inv | other._inv);
        }
        
        array_smashing_t operator&(array_smashing_t other) {
          return array_smashing_t (_inv & other._inv);
        }
        
        array_smashing_t operator||(array_smashing_t other) {
          return array_smashing_t (_inv || other._inv);
        }

        template<typename Thresholds>
        array_smashing_t widening_thresholds (array_smashing_t other, 
                                              const Thresholds &ts) {
          return array_smashing_t (_inv.widening_thresholds (other._inv, ts));
        }
        
        array_smashing_t operator&& (array_smashing_t other) {
          return array_smashing_t (_inv && other._inv);
        }
        

        // remove all variables [begin,...end)
        template<typename Iterator>
        void forget (Iterator begin, Iterator end) {
          domain_traits<NumDomain>::forget (_inv, begin, end);
        }

        // dual of forget: remove all variables except [begin,...end)
        template<typename Iterator>
        void project (Iterator begin, Iterator end) {
          domain_traits<NumDomain>::project (_inv, begin, end);
        }

        void operator += (linear_constraint_system_t csts) {
          _inv += csts;
        }

        void operator-=(VariableName var) {
          _inv -= var;
        }
        
        void assign (VariableName x, linear_expression_t e) {
          _inv.assign (x, e);
          
          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< e<< *this <<"\n";);
        }
        
        void apply (operation_t op, VariableName x, VariableName y, Number z) {
          _inv.apply (op, x, y, z);
          
          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< y<< " "<< op<< " "<< z<< *this <<"\n";);
        }
        
        void apply(operation_t op, VariableName x, VariableName y, VariableName z) {
          _inv.apply (op, x, y, z);
          
          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< y<< " "<< op<< " "<< z<< *this <<"\n";);
        }
        
        void apply(operation_t op, VariableName x, Number k) {
          _inv.apply (op, x, k);
          
          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< x<< " "<< op<< " "<< k<< *this <<"\n";);
        }
	   	
	void backward_assign (VariableName x, linear_expression_t e,
			      array_smashing_t inv) {
	  _inv.backward_assign(x, e, inv.get_content_domain());
	}
	
	void backward_apply (operation_t op,
			     VariableName x, VariableName y, Number z,
			     array_smashing_t inv) {
	  _inv.backward_apply(op, x, y, z, inv.get_content_domain());
	}
	
	void backward_apply(operation_t op,
			    VariableName x, VariableName y, VariableName z,
			    array_smashing_t inv) {
	  _inv.backward_apply(op, x, y, z, inv.get_content_domain());
	}
	
        void apply(int_conv_operation_t op,
		   VariableName dst, unsigned dst_width, VariableName src, unsigned src_width) {
          _inv.apply (op, dst, dst_width, src, src_width);
        }
                
        void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
          _inv.apply (op, x, y, z);

          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< y<< " "
		                << op<< " "<< z<< *this <<"\n";);
        }
        
        void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
          _inv.apply (op, x, y, k);

          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< y<< " "
		                << op<< " "<< k<< *this <<"\n";);
        }
        
        // division_operators_api
        void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
          _inv.apply (op, x, y, z);

          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< y<< " "<< op<< " "<< z<< *this <<"\n";);
        }
        
        void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
          _inv.apply (op, x, y, k);

          CRAB_LOG("smashing",
                   crab::outs() << "apply "<< x<< " := "<< y<< " "<< op<< " "<< k<< *this <<"\n";);
        }

	// boolean operators
	virtual void assign_bool_cst (VariableName lhs, linear_constraint_t rhs) override {
	  _inv.assign_bool_cst (lhs, rhs);
	}    
	
	virtual void assign_bool_var(VariableName lhs, VariableName rhs, bool is_not_rhs) override {
	  _inv.assign_bool_var (lhs, rhs, is_not_rhs);
	}    
	
	virtual void apply_binary_bool (bool_operation_t op,VariableName x,
					VariableName y,VariableName z) override {
	  _inv.apply_binary_bool (op, x, y, z);
	}    
	
	virtual void assume_bool (VariableName v, bool is_negated) override {
	  _inv.assume_bool (v, is_negated);
	}    

	// backward boolean operators
	virtual void backward_assign_bool_cst(VariableName lhs, linear_constraint_t rhs,
					      array_smashing_t inv){
	  _inv.backward_assign_bool_cst(lhs, rhs, inv.get_content_domain());	  
	}
	
	virtual void backward_assign_bool_var(VariableName lhs, VariableName rhs, bool is_not_rhs,
					      array_smashing_t inv) {
	  _inv.backward_assign_bool_var(lhs, rhs, is_not_rhs, inv.get_content_domain());	  	  
	}
	
	virtual void backward_apply_binary_bool(bool_operation_t op,
						VariableName x,VariableName y,VariableName z,
						array_smashing_t inv) {
	  _inv.backward_apply_binary_bool(op, x, y, z, inv.get_content_domain());
	}
	
        // pointer_operators_api
        virtual void pointer_load (VariableName lhs, VariableName rhs) override {
          _inv.pointer_load(lhs,rhs);
        }
        
        virtual void pointer_store (VariableName lhs, VariableName rhs) override {
          _inv.pointer_store(lhs,rhs);
        } 
        
        virtual void pointer_assign (VariableName lhs, VariableName rhs, linear_expression_t offset) override {
          _inv.pointer_assign (lhs,rhs,offset);
        }
        
        virtual void pointer_mk_obj (VariableName lhs, ikos::index_t address) override {
          _inv.pointer_mk_obj (lhs, address);
        }
        
        virtual void pointer_function (VariableName lhs, VariableName func) override {
          _inv.pointer_function (lhs, func);
        }
        
        virtual void pointer_mk_null (VariableName lhs) override {
          _inv.pointer_mk_null (lhs);
        }
        
        virtual void pointer_assume (ptr_cst_t cst) override {
          _inv.pointer_assume (cst);
        }    
        
        virtual void pointer_assert (ptr_cst_t cst) override {
          _inv.pointer_assert (cst);
        }    
        
        // array_operators_api 

        // All the array elements are assumed to be equal to val
        virtual void array_assume (VariableName a, variable_type a_ty, 
                                   linear_expression_t /*lb_idx*/,
				   linear_expression_t /*ub_idx*/, 
                                   linear_expression_t val) override {
	  if (a_ty == ARR_BOOL_TYPE)  {
	    if (val.is_constant()) {
	      if (val.constant() >= Number(1))
		_inv.assign_bool_cst(a, linear_constraint_t::get_true());
	      else
		_inv.assign_bool_cst(a, linear_constraint_t::get_false());
	    } else if (auto var = val.get_variable ()) {
	      _inv.assign_bool_var (a, (*var).name(), false);
	    }
	  } else if (a_ty == ARR_INT_TYPE || a_ty == ARR_REAL_TYPE) {
            _inv.assign (a, val);
	  } else if (a_ty == ARR_PTR_TYPE) {
	    if (val.is_constant() && val.constant() == Number(0))
	      _inv.pointer_mk_null(a);
	    else if (auto var = val.get_variable ()) {
	      _inv.pointer_assign (a, (*var).name(), Number(0));
	    }
	  }
          
          CRAB_LOG("smashing",
		   crab::outs() << "forall i:: " << a << "[i]==" << val
		                << " -- " << *this <<"\n";);
        }
        
        virtual void array_load (VariableName lhs, VariableName a,
				 crab::variable_type a_ty,
                                 linear_expression_t i,
				 z_number /*bytes*/) override {

          crab::CrabStats::count (getDomainName() + ".count.load");
          crab::ScopedCrabStats __st__(getDomainName() + ".load");
          
          // We need to be careful when assigning a summarized variable a
          // into a non-summarized variable lhs. Simply _inv.assign (lhs,
          // a) is not sound.
          /* ask for a temp var */
          VariableName a_prime = a.get_var_factory().get(); 
          domain_traits<NumDomain>::expand (_inv, a, a_prime);
	  if (a_ty == ARR_BOOL_TYPE) {
	    _inv.assign_bool_var(lhs, a_prime, false);
	  } else if (a_ty == ARR_INT_TYPE || a_ty == ARR_REAL_TYPE) {
            _inv.assign (lhs, linear_expression_t (a_prime));
	  } else if (a_ty == ARR_PTR_TYPE) {
            _inv.pointer_assign (lhs, a_prime, Number(0));
	  }

          _inv -= a_prime; 
          
          CRAB_LOG("smashing",
		   crab::outs() << lhs << ":=" << a <<"[" << i << "]  -- "
		                << *this <<"\n";);
        }
        
        
        virtual void array_store (VariableName a, crab::variable_type a_ty,
                                  linear_expression_t i,
				  linear_expression_t val, 
                                  z_number /*bytes*/,
				  bool is_singleton) override {
                                  

          crab::CrabStats::count (getDomainName() + ".count.store");
          crab::ScopedCrabStats __st__(getDomainName() + ".store");
          
          if (is_singleton)
            strong_update (a, a_ty, val);
          else 
            weak_update (a, a_ty, val);
          
          CRAB_LOG("smashing",
		   crab::outs() << a << "[" << i << "]:="
		                << val << " -- " << *this <<"\n";);
        }

        virtual void array_assign (VariableName lhs, VariableName rhs, 
                                   crab::variable_type ty) override {
	  if (ty == ARR_BOOL_TYPE) {
	    _inv.assign_bool_var(lhs, rhs, false);
	  } else if (ty == ARR_INT_TYPE || ty == ARR_REAL_TYPE) {
            _inv.assign (lhs, linear_expression_t(rhs));
	  } else  if (ty == ARR_PTR_TYPE) {
            _inv.pointer_assign (lhs, rhs, Number(0));
	  }
        }
        
        linear_constraint_system_t to_linear_constraint_system (){
          return _inv.to_linear_constraint_system ();
        }

        disjunctive_linear_constraint_system_t
	to_disjunctive_linear_constraint_system (){
          return _inv.to_disjunctive_linear_constraint_system ();
        }
	
        NumDomain  get_content_domain () const {      
          return _inv;
        }

        void write(crab_os& o) {
          o << _inv;
        }
        
        static std::string getDomainName () {
          std::string name ("ArraySmashing(" + 
                            NumDomain::getDomainName () + 
                            ")");
          return name;
        }  

	void rename(const varname_vector_t& from, const varname_vector_t &to){
	  _inv.rename(from, to);
	}
	
      }; // end array_smashing
   
     template<typename BaseDomain>
     class domain_traits <array_smashing<BaseDomain> > {
      public:
       
       typedef array_smashing<BaseDomain> array_smashing_t;
       typedef typename BaseDomain::varname_t VariableName;

       template<class CFG>
       static void do_initialization (CFG cfg) { }

       static void normalize (array_smashing_t& inv) { 
         CRAB_WARN ("array smashing normalize not implemented");
       }
       
       template <typename Iter>
       static void forget (array_smashing_t& inv, Iter it, Iter end) {
         inv.forget (it, end);
       }
       
       template <typename Iter >
       static void project (array_smashing_t& inv, Iter it, Iter end) {
         inv.project (it, end);
       }

       static void expand (array_smashing_t& inv, VariableName x, VariableName new_x) {
         // -- lose precision if relational or disjunctive domain
         CRAB_WARN ("array smashing expand not implemented");
       }

     };

   } // namespace domains
}// namespace crab
#endif 
