#ifndef FLAT_BOOLEAN_DOMAIN_HPP
#define FLAT_BOOLEAN_DOMAIN_HPP

/* 
   A simple flat 3-valued boolean domain and a reduced product of this
   flat bool domain with an arbitrary numerical domain 
*/
   
#include <crab/common/types.hpp>
#include <crab/domains/operators_api.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/discrete_domains.hpp>
#include <crab/domains/combined_domains.hpp>
#include <crab/domains/intervals.hpp>

namespace crab {

  namespace domains {  

  class boolean_value : public ikos::writeable {
    /*
              Top
              / \
             /   \
          True  False
             \   /
              \ /
             Bottom
    */
    typedef enum {
      False  = 0x0,      
      True   = 0x1,
      Bottom = 0x2,      
      Top    = 0x3 
    } kind_t;
    
    kind_t _value;
    
    boolean_value(kind_t v) : _value(v){};
    
   public:
    
    boolean_value() : _value(Top) {}
    
    static boolean_value bottom() { return boolean_value(Bottom); }

    static boolean_value top() { return boolean_value(Top); }

    static boolean_value get_true() { return boolean_value(True); }

    static boolean_value get_false() { return boolean_value(False); }

    boolean_value(const boolean_value& other) : 
        writeable(), _value(other._value) {}

    boolean_value& operator=(const boolean_value& other) {
      if (this != &other) _value = other._value; 
      return *this;
    }

    bool is_bottom() { return (_value == Bottom); }

    bool is_top() { return (_value == Top); }
    
    bool is_true() const { return (_value == True); }

    bool is_false() const { return (_value == False); }
    
    bool operator<=(boolean_value other) {

      if (_value == Bottom || other._value == Top)
        return true;
      else if (_value == Top)
        return (other._value == Top);
      else if (_value == True)
        return ((other._value == True) || (other._value == Top));      
      else if (_value == False)
        return ((other._value == False) || (other._value == Top));      	
      // this should be unreachable
      return false;
    }

    bool operator==(boolean_value other) {
      return (_value == other._value);
    }

    boolean_value operator|(boolean_value other) {
      if (is_bottom ()) return other;
      if (other.is_bottom ()) return *this;
      if (is_top () || other.is_top ()) return top ();
      if (_value == other._value) return *this;

      // othewise true | false or false | true ==> top
      return top ();
    }

    boolean_value operator&(boolean_value other) {
      if (is_bottom ()) return *this;
      if (other.is_bottom ()) return other;
      if (is_top ()) return other;
      if (other.is_top ()) return *this;
      if (_value == other._value) return *this;

      // othewise true & false or false & true ==> bottom
      return bottom ();
    }

    // the lattice satisfy ACC so join is the widening
    boolean_value operator||(boolean_value other) { 
      return this->operator|(other); 
    }

    // the lattice satisfy DCC so meet is the narrowing
    boolean_value operator&&(boolean_value other) { 
      return this->operator&(other); 
    }

    // Boolean operations
    /*  
              And  Or  X0r
         0 0   0   0    0
         0 1   0   1    1
         1 0   0   1    1
         1 1   1   1    0
         0 *   0   *    *
         * 0   0   *    *
         1 *   *   1    *
         * 1   *   1    *
         * *   *   *    *
    */
    
    boolean_value And(boolean_value other) {
      if (is_bottom () || other.is_bottom ())
	return bottom ();
	
      if (!is_top () && !other.is_top ())
	return boolean_value(static_cast< kind_t >(
               static_cast< int >(this->_value) & static_cast< int >(other._value)));

      int x = static_cast< int >(this->_value);
      int y = static_cast< int >(other._value);
      if (x == 0 || y == 0)
	return get_false ();
      else
	return top ();
    }

    boolean_value Or(boolean_value other) {
      if (is_bottom () || other.is_bottom ())
	return bottom ();

      if (!is_top () && !other.is_top ())
	return boolean_value(static_cast< kind_t >(
               static_cast< int >(this->_value) | static_cast< int >(other._value)));

      int x = static_cast< int >(this->_value);
      int y = static_cast< int >(other._value);
      if (x == 1 || y == 1)
	return get_true ();
      else
	return top ();
      
    }

    boolean_value Xor(boolean_value other) {
      if (is_bottom () || other.is_bottom ())
	return bottom ();

      if (!is_top () && !other.is_top ())
	return boolean_value(static_cast< kind_t >(
               static_cast< int >(this->_value) ^ static_cast< int >(other._value)));
      else
	return top ();
    }

    
    void write(crab_os& o) {
      switch (_value) {
        case Bottom:      o << "_|_"; break;
        case Top:         o << "*"; break;
        case True:        o << "true"; break;
        default:/*False*/ o << "false";  
      }
    }

  }; // end class boolean_value


  // A simple abstract domain for booleans
  template <typename Number, typename VariableName>
  class flat_boolean_domain : 
        public ikos::writeable,
	public boolean_operators<Number, VariableName>,
        public numerical_domain<Number, VariableName>,
        public bitwise_operators<Number,VariableName>,
        public division_operators<Number, VariableName>,
        public array_operators< Number, VariableName>,
        public pointer_operators< Number, VariableName> {
    
    typedef separate_domain< VariableName, boolean_value> separate_domain_t;
    typedef flat_boolean_domain<Number, VariableName> flat_boolean_domain_t;

   public:

    typedef VariableName varname_t;
    typedef Number number_t;
    typedef boolean_value bool_t;
    typedef linear_expression<Number, VariableName> linear_expression_t;
    typedef linear_constraint<Number, VariableName> linear_constraint_t;
    typedef linear_constraint_system<Number, VariableName> linear_constraint_system_t;
    typedef interval<Number>  interval_t;    
    typedef typename separate_domain_t::iterator iterator;

   private:

    separate_domain_t _env;
    
    flat_boolean_domain(separate_domain_t env) : _env(env) {}
        
   public:
    
    static flat_boolean_domain_t top() {
      return flat_boolean_domain(separate_domain_t::top());
    }
    
    static flat_boolean_domain_t bottom() {
      return flat_boolean_domain(separate_domain_t::bottom());
    }
    
    flat_boolean_domain() : _env(separate_domain_t::top()) {}

    flat_boolean_domain(const flat_boolean_domain_t& e) : 
      writeable(), _env(e._env) {}
    
    flat_boolean_domain_t& operator=(const flat_boolean_domain_t& o) {
      if (this != &o)
        _env = o._env;
      return *this;
    }
    
    iterator begin() { 
      if (is_bottom ()) CRAB_ERROR ("Cannot return iterator from bottom"); 
      return _env.begin(); 
    }
    
    iterator end() { 
      if (is_bottom ()) CRAB_ERROR ("Cannot return iterator from bottom");
      return _env.end(); 
    }
    
    bool is_bottom() { return _env.is_bottom(); }
    
    bool is_top() { return _env.is_top(); }
    
    bool operator<=(flat_boolean_domain_t o) { return (_env <= o._env); }
    
    flat_boolean_domain_t operator|(flat_boolean_domain_t o) {
      flat_boolean_domain_t res (_env | o._env);
      CRAB_LOG("flat-boolean",
	       crab::outs () << "After join " << *this << " and " << o << "="
	                     << res << "\n";);
      return res;
    }

    void operator|=(flat_boolean_domain_t o) {
      CRAB_LOG("flat-boolean",
	       crab::outs () << "After join " << *this << " and "
	                     << o << "=");
      _env = _env  | o._env;
      CRAB_LOG("flat-boolean", crab::outs () << *this << "\n");
    }

    flat_boolean_domain_t operator&(flat_boolean_domain_t o) {
      flat_boolean_domain_t res (_env & o._env);
      CRAB_LOG("flat-boolean",
	       crab::outs () << "After meet " << *this << " and " << o
	                     << "=" << res << "\n");
      return res;
    }
    
    flat_boolean_domain_t operator||(flat_boolean_domain_t o) {
      flat_boolean_domain_t res (_env || o._env);
      CRAB_LOG("flat-boolean",
	       crab::outs () << "After widening " << *this << " and "
	                     << o << "=" << res << "\n");
      return res;
    }

    template<typename Thresholds>
    flat_boolean_domain_t widening_thresholds(flat_boolean_domain_t o,
					      const Thresholds &) {
    
      flat_boolean_domain_t res (_env || o._env);
      CRAB_LOG("flat-boolean",
	       crab::outs () << "After widening " << *this << " and "
	                     << o << "=" << res << "\n");
      return res;
    }
    
    flat_boolean_domain_t operator&&(flat_boolean_domain_t o)
    { return (_env && o._env); }
      
    void operator-=(VariableName v)
    { 
      if (!is_bottom ())
        _env -= v; 
    }
        
    // flat_boolean_domains_api

    // XXX: the flat boolean domain cannot reason about linear
    // constraints so we assign top to x.
    void assign_bool_cst (VariableName x, linear_constraint_t cst)
    {
      _env -= x;
      CRAB_LOG("flat-boolean",
	       auto bx = _env[x];
	       crab::outs () << x << ":=" << bx << "\n");
    }    

    void assign_bool_var (VariableName x, VariableName y)
    {
      _env.set (x, _env [y]);
      CRAB_LOG("flat-boolean",
	       auto bx = _env[x];
	       crab::outs () << "After " << x << ":=" << y << " --->"
 	                     << x << "=" << bx << "\n");
    }
    
    void apply_binary_bool(bool_operation_t op,
			   VariableName x, VariableName y, VariableName z) {
      switch (op) {
      case OP_BAND:
	_env.set (x, _env [y].And (_env[z]));
	break;
      case OP_BOR:
	_env.set (x, _env [y].Or (_env[z]));
	break;
      case OP_BXOR:
	_env.set (x, _env [y].Xor (_env[z]));
	break;
      default:
	CRAB_ERROR ("Unknown boolean operator");
      }
      
      CRAB_LOG("flat-boolean",
	       auto bx = _env[x];
	       crab::outs () << "After " << x << ":=" << y << " " << op << " " << z
 	                                 << " --->" << x << "=" << bx << "\n");
    }

    void assume_bool (VariableName x, bool is_negated) {
      if (!is_negated)
	_env.set (x,  _env[x] & boolean_value::get_true ());
      else 
	_env.set (x,  _env[x] & boolean_value::get_false ());

      CRAB_LOG("flat-boolean",
	       auto bx = _env[x];
	       if (!is_negated)
		 crab::outs () << "After assume(" << x << ") --> "
			       << x << "=" << bx << "\n";
	       else
		 crab::outs () << "After assume(not(" << x << ")) --> "
			       << x << "=" << bx << "\n";);
    }

    // XXX: these methods are not actually part of boolean_operators
    // api but they are used by flat_boolean_numerical_domain and
    // domain_traits.
    
    void set_bool (VariableName x, boolean_value v)
    { _env.set(x, v); }
    
    boolean_value get_bool (VariableName x)
    { return _env[x];}
        
    // numerical_domains_api
    // XXX: needed for making a reduced product with a numerical domain
    void apply(operation_t op, VariableName x, VariableName y, VariableName z) {}
    void apply(operation_t op, VariableName x, VariableName y, Number k) {}
    void assign(VariableName x, linear_expression_t e) {}
    void operator+=(linear_constraint_system_t csts) {}
    void operator+=(linear_constraint_t cst) {}
    // not part of the numerical_domains api but it should be
    void set (VariableName x, interval_t intv) {}
    interval_t operator[](VariableName x) { return interval_t::top ();}

    // division_operators_api
    // XXX: needed for making a reduced product with a numerical domain
    void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {}
    void apply(div_operation_t op, VariableName x, VariableName y, Number z) {}
    
    // bitwise_operators_api
    // XXX: needed for making a reduced product with a numerical domain
    void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {}
    void apply(conv_operation_t op, VariableName x, Number y, unsigned width) {}
    void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {}
    void apply(bitwise_operation_t op, VariableName x, VariableName y, Number z) {}

    static std::string getDomainName () {return "Boolean"; }
    void write(crab_os& o) { _env.write(o); }
    
   }; // class flat_boolean_domain


   template <typename Number, typename VariableName>
   class domain_traits <flat_boolean_domain<Number,VariableName> > {

     typedef flat_boolean_domain<Number,VariableName> flat_boolean_domain_t;

    public:

     template<class CFG>
     static void do_initialization (CFG cfg) { }

     // Normalize the abstract domain if such notion exists.
     static void normalize (flat_boolean_domain_t& inv) { }

     // Remove all variables [begin, end)
     template<typename Iter>
     static void forget (flat_boolean_domain_t& inv, Iter begin, Iter end) {
       for (auto v : boost::make_iterator_range (begin,end)){
         inv -= v; 
       }
     }

     // Forget all variables except [begin, end)
     template <typename Iter>
     static void project(flat_boolean_domain_t& inv, Iter begin, Iter end){
       flat_boolean_domain_t res = flat_boolean_domain_t::top ();
       for (auto v : boost::make_iterator_range (begin, end))
         res.set_bool (v, inv.get_bool(v)); 
       std::swap (inv, res);
     }
         
     // Make a new copy of x without relating x with new_x
     static void expand (flat_boolean_domain_t& inv, VariableName x, VariableName new_x) {
       inv.set_bool (new_x , inv.get_bool (x));
     }
     
   };


    // Reduced product of the flat boolean domain with an arbitrary
    // numerical domain.
    // The reduction happens in two situations:
    //    (1) when bvar := linear_constraint
    //    (2) when assume (bvar) or assume (!bar)
    // The step (2) is quite weak.
    template <typename NumDom>
    class flat_boolean_numerical_domain:
        public writeable,
        public numerical_domain<typename NumDom::number_t, 
                                typename NumDom::varname_t>,
        public bitwise_operators<typename NumDom::number_t, 
                                 typename NumDom::varname_t>,
        public division_operators<typename NumDom::number_t,
                                  typename NumDom::varname_t>,
        public array_operators<typename NumDom::number_t,
                               typename NumDom::varname_t>,
        public pointer_operators<typename NumDom::number_t,
                                 typename NumDom::varname_t>,
	public boolean_operators <typename NumDom::number_t,
				  typename NumDom::varname_t>
    {
      
      typedef typename NumDom::number_t N;
      typedef typename NumDom::varname_t V;
      
     public:
      
      typedef flat_boolean_domain <N,V> bool_domain_t;            
      typedef flat_boolean_numerical_domain<NumDom> bool_num_domain_t;
      using typename numerical_domain<N, V>::linear_expression_t;
      using typename numerical_domain<N, V>::linear_constraint_t;
      using typename numerical_domain<N, V>::linear_constraint_system_t;
      using typename numerical_domain<N, V>::variable_t;
      typedef typename NumDom::number_t number_t;
      typedef typename NumDom::varname_t varname_t;      
      typedef interval<N> interval_t;
      typedef crab::pointer_constraint<V> ptr_cst_t;
      
     private:

      typedef domain_product2<N,V,bool_domain_t,NumDom> domain_product2_t;

      // XXX: for performing reduction from the boolean domain to the
      // numerical one.
      typedef discrete_domain<linear_constraint_t> lin_cst_set_t;
      // Map bool variables to sets of constraints such that if the
      // bool variable is true then the conjunction of the constraints
      // is satisfiable.
      typedef separate_domain<V, lin_cst_set_t> var_to_csts_map_t;
      
      domain_product2_t _product;
      var_to_csts_map_t _var_to_csts;
      
      flat_boolean_numerical_domain(const domain_product2_t& product,
				    const var_to_csts_map_t& var_to_csts):
	_product(product), _var_to_csts (var_to_csts) {}
      
     public:
      
      static bool_num_domain_t top() {
        return bool_num_domain_t (domain_product2_t::top(), var_to_csts_map_t());
      }
      
      static bool_num_domain_t bottom() {
        return bool_num_domain_t(domain_product2_t::bottom(), var_to_csts_map_t());
      }
      
     public:
       
      flat_boolean_numerical_domain() : _product(), _var_to_csts() {}
      
      flat_boolean_numerical_domain(const bool_num_domain_t& other) :
	_product(other._product), _var_to_csts (other._var_to_csts) { }
      
      bool_num_domain_t& operator=(const bool_num_domain_t& other) {
        if (this != &other) {
	  this->_product = other._product;
	  this->_var_to_csts = other._var_to_csts;
	}
        return *this;
      }
      
      bool is_bottom() { return this->_product.is_bottom(); }
      
      bool is_top() { return this->_product.is_top(); }
      
      bool_domain_t& first() { return this->_product.first(); }
      
      NumDom& second() { return this->_product.second(); }
      
      bool operator<=(bool_num_domain_t other)
      { return this->_product <= other._product; }
      
      bool operator==(bool_num_domain_t other)
      { return this->_product == other._product; }
      
      void operator|=(bool_num_domain_t other)
      {
	this->_product |= other._product;
	this->_var_to_csts = this->_var_to_csts | other._var_to_csts;
      }

      bool_num_domain_t operator|(bool_num_domain_t other)
      { return bool_num_domain_t(this->_product | other._product,
				 this->_var_to_csts | other._var_to_csts);}
      
      bool_num_domain_t operator&(bool_num_domain_t other)
      { return bool_num_domain_t(this->_product & other._product,
				 this->_var_to_csts & other._var_to_csts);}
      
      bool_num_domain_t operator||(bool_num_domain_t other)
      { return bool_num_domain_t(this->_product || other._product,
				 this->_var_to_csts || other._var_to_csts);}
       
      template<typename Thresholds>
      bool_num_domain_t widening_thresholds (bool_num_domain_t other, const Thresholds& ts)
      { return bool_num_domain_t(this->_product.widening_thresholds (other._product, ts),
				 this->_var_to_csts || other._var_to_csts);}
      
      bool_num_domain_t operator&&(bool_num_domain_t other)
      { return bool_num_domain_t(this->_product && other._product,
				 this->_var_to_csts && other._var_to_csts);}
      
      // numerical_domains_api

      void apply(operation_t op, V x, V y, V z)
      { this->_product.apply(op, x, y, z); }
      
      void apply(operation_t op, V x, V y, N k)
      { this->_product.apply(op, x, y, k); }

      void assign(V x, linear_expression_t e)
      { this->_product.assign(x, e); }
      
      void operator+=(linear_constraint_system_t csts)
      { this->_product += csts; }
      
      void set (V v, interval_t intv)
      { this->_product.set(v, intv); }
            
      interval_t operator[](V v)
      { return this->_product[v]; }

      void operator-=(V v)
      {
	this->_product -= v;
	this->_var_to_csts -= v;
      }
      
      // boolean_operators
      
      void assign_bool_cst (V x, linear_constraint_t cst)
      { /// Reduction from the numerical domain to the flat boolean
	/// domain
	
	if (this->_product.is_bottom ())
	  return;

	NumDom inv1 (this->_product.second ());
	inv1 += cst;
	if (inv1.is_bottom ()) {
	  // XXX: definitely false
	  this->_product.first().set_bool (x, boolean_value::get_false ());
	} else {
	  // XXX: definitely true	  
	  NumDom inv2 (this->_product.second ());
	  inv2  += cst.negate();
	  if (inv2.is_bottom ())
	    this->_product.first().set_bool (x, boolean_value::get_true ());
	  else
	    this->_product.first().set_bool (x, boolean_value::top ());
	}
	this->_var_to_csts.set (x, cst);
	CRAB_LOG ("flat-boolean",
		  auto bx = this->_product.first ().get_bool (x);
		  crab::outs () << "Reduction numerical --> boolean: "
		                << x << " := " << bx << "\n";);
      }

      void assign_bool_var (V x, V y)
      {
	this->_product.assign_bool_var (x, y);
	this->_var_to_csts.set (x, this->_var_to_csts [y]);	
      }
      
      void apply_binary_bool(bool_operation_t op, V x, V y, V z)
      {
	this->_product.apply_binary_bool (op, x, y, z);

	// XXX: Collect information for the reduction from the boolean
	// to the numerical domain.
	if (op == OP_BAND)
	{
	  this->_var_to_csts.set (x, this->_var_to_csts [y] | this->_var_to_csts [z]);
	  return;
	}

	if (op == OP_BOR || op == OP_BXOR)
	{
	  if (this->_product.first().get_bool (y).is_false())
	  {
	    this->_var_to_csts.set (x, this->_var_to_csts [z]);
	    return;
	  }
	  if (this->_product.first().get_bool (z).is_false())
	  {
	    this->_var_to_csts.set (x, this->_var_to_csts [y]);
	    return;
	  }
	}
	
	/// otherwise we give up
	this->_var_to_csts -= x;
      }

      void assume_bool(V x, bool is_negated)
      {
	this->_product.assume_bool (x, is_negated);

	if (this->_var_to_csts [x].is_top () ||
	    this->_var_to_csts [x].is_bottom ())
	  return;
	
	CRAB_LOG ("flat-boolean",
		  crab::outs () << "Before reduction boolean --> numerical: "
		                << this->_product.second () << "\n";);
	  
	// Reduction from the flat boolean domain to the numerical
	// domain
	for (auto cst: this->_var_to_csts [x])
	{
	  if (!is_negated)
	    this->_product.second() += cst;
	  else
	    this->_product.second() += cst.negate();
	}

	CRAB_LOG ("flat-boolean",
		  crab::outs () << "After reduction boolean --> numerical: "
				<< this->_product.second () << "\n";);
      }

      // bitwise_operators_api
      
      void apply(conv_operation_t op, V x, V y, unsigned width)
      { this->_product.apply(op, x, y, width); }
      
      void apply(conv_operation_t op, V x, N k, unsigned width)
      { this->_product.apply(op, x, k, width); }
      
      void apply(bitwise_operation_t op, V x, V y, V z)
      { this->_product.apply(op, x, y, z); }
      
      void apply(bitwise_operation_t op, V x, V y, N k)
      { this->_product.apply(op, x, y, k); }
      
      // division_operators_api
      
      void apply(div_operation_t op, V x, V y, V z)
      { this->_product.apply(op, x, y, z); }
      
      void apply(div_operation_t op, V x, V y, N k)
      { this->_product.apply(op, x, y, k); }
      
      // array_operators_api
      
      virtual void array_assume (V a, crab::variable_type a_ty, 
                                 linear_expression_t lb_idx, linear_expression_t ub_idx,
                                 V var) override
      { this->_product.array_assume (a, a_ty, lb_idx, ub_idx, var); }
      
      virtual void array_load (V lhs, V a, crab::variable_type a_ty, 
                               linear_expression_t i, ikos::z_number bytes) override
      { this->_product.array_load (lhs, a, a_ty, i, bytes); }
      
      virtual void array_store (V a, crab::variable_type a_ty, 
                                linear_expression_t i, V val, 
                                ikos::z_number bytes, bool is_singleton) override
      { this->_product.array_store (a, a_ty, i, val, bytes, is_singleton); }
      
      virtual void array_assign (V lhs, V rhs, crab::variable_type ty) override
      { this->_product.array_assign (lhs, rhs, ty); }
      
      // pointer_operators_api
      virtual void pointer_load (V lhs, V rhs) override
      {  this->_product.pointer_load (lhs, rhs); }
      
      virtual void pointer_store (V lhs, V rhs) override
      { this->_product.pointer_store (lhs, rhs); }
      
      virtual void pointer_assign (V lhs, V rhs, linear_expression_t offset) override
      { this->_product.pointer_assign (lhs, rhs, offset); }
      
      virtual void pointer_mk_obj (V lhs, ikos::index_t address) override
      { this->_product.pointer_mk_obj (lhs, address); }
      
      virtual void pointer_function (V lhs, V func) override
      { this->_product.pointer_function (lhs, func); }
      
      virtual void pointer_mk_null (V lhs) override
      { this->_product.pointer_mk_null (lhs); }
      
      virtual void pointer_assume (ptr_cst_t cst) override
      { this->_product.pointer_assume (cst); }
      
      virtual void pointer_assert (ptr_cst_t cst) override
      { this->_product.pointer_assert (cst); }
      
      void write(crab_os& o)
      { this->_product.write(o); }
      
      linear_constraint_system_t to_linear_constraint_system()
      { return this->_product.second().to_linear_constraint_system(); }
      
      static std::string getDomainName()
      { return domain_product2_t::getDomainName (); }

      // domain_traits_api
      
      void expand(V x, V new_x) {
        crab::domains::domain_traits<bool_domain_t>::
	  expand (this->_product.first(), x, new_x);
        crab::domains::domain_traits<NumDom>::
	  expand (this->_product.second(), x, new_x);
	// TODO: var_to_csts_map
      }
      
      void normalize() {
        crab::domains::domain_traits<bool_domain_t>::
	  normalize(this->_product.first());
        crab::domains::domain_traits<NumDom>::
	  normalize(this->_product.second());
      }
      
      template <typename Range>
      void forget(Range vars){
        crab::domains::domain_traits<bool_domain_t>::
	  forget(this->_product.first(), vars.begin (), vars.end());
        crab::domains::domain_traits<NumDom>::
	  forget(this->_product.second(), vars.begin (), vars.end());
	// TODO: var_to_csts_map	
      }
      
      template <typename Range>
      void project(Range vars) {
        crab::domains::domain_traits<bool_domain_t>::
	  project(this->_product.first(), vars.begin(), vars.end());
        crab::domains::domain_traits<NumDom>::
	  project(this->_product.second(), vars.begin(), vars.end());
	// TODO: var_to_csts_map	
      }
      
    }; // class flat_boolean_numerical_domain

    template<typename Num>
    class domain_traits <flat_boolean_numerical_domain<Num> > {
     public:

      typedef flat_boolean_numerical_domain<Num> product_t;
      typedef typename product_t::varname_t V;

      template<class CFG>
      static void do_initialization (CFG cfg) { }

      static void normalize (product_t& inv) {
        inv.normalize();
      }

      static void expand (product_t& inv, V x, V new_x) {
        inv.expand (x, new_x);
      }
      
      template <typename Iter>
      static void forget (product_t& inv, Iter it, Iter end){
        inv.forget (boost::make_iterator_range (it, end));
      }
      
      template <typename Iter>
      static void project (product_t& inv, Iter it, Iter end) {
        inv.project (boost::make_iterator_range (it, end));
      }
    };
    
    
  } // end namespace domains 
} // end namespace crab

#endif 
