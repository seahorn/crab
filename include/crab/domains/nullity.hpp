#ifndef NULLITY_HPP
#define NULLITY_HPP

/* A flat lattice for nullity */

#include <crab/common/types.hpp>
#include <crab/domains/separate_domains.hpp>
#include <crab/domains/operators_api.hpp>

namespace crab {

  namespace domains {  

  class nullity_value : public ikos::writeable {
    /*
              Top
              / \
             /   \
          Null  NonNull
             \   /
              \ /
             Bottom
    */
    typedef enum {
      Bottom = 0x0, /*unused*/
      Null = 0x1,
      NonNull = 0x2,
      Top = 0x3 
    } kind_t;
    
    kind_t _value;
    
    nullity_value(kind_t v) : _value(v){};
    
   public:
    
    nullity_value() : _value(Top) {}
    
    static nullity_value bottom() { return nullity_value(Bottom); }

    static nullity_value top() { return nullity_value(Top); }

    static nullity_value non_null() { return nullity_value(NonNull); }

    static nullity_value null() { return nullity_value(Null); }

    nullity_value(const nullity_value& other) : 
        writeable(), _value(other._value) {}

    nullity_value& operator=(const nullity_value& other) {
      if (this != &other) {
        _value = other._value;
      }
      return *this;
    }

    bool is_bottom() { return (_value == Bottom); }

    bool is_top() { return (_value == Top); }
    
    bool is_non_null() const { return (_value == NonNull); }

    bool is_null() const { return (_value == Null); }
    
    bool operator<=(nullity_value other) {

      if (_value == Bottom || other._value == Top)
        return true;
      else if (_value == Top)
        return (other._value == Top);
      else if (_value == NonNull)
        return (other._value >= NonNull);
      else if (_value == Null)
        return ((other._value == Null) || (other._value == Top));
      
      // this should be unreachable
      return false;
    }

    bool operator==(nullity_value other) {
      return (_value == other._value);
    }

    nullity_value operator|(nullity_value other) {
      return nullity_value(static_cast< kind_t >(
          static_cast< int >(this->_value) | static_cast< int >(other._value)));
    }

    // the lattice satisfy ACC so join is the widening
    nullity_value operator||(nullity_value other) { 
      return this->operator|(other); 
    }

    nullity_value operator&(nullity_value other) {
      return nullity_value(static_cast< kind_t >(
          static_cast< int >(this->_value) & static_cast< int >(other._value)));
    }
    
    // the lattice satisfy DCC so meet is the narrowing
    nullity_value operator&&(nullity_value other) { 
      return this->operator&(other); 
    }

    void write(crab_os& o) {
      switch (_value) {
        case Bottom:     o << "_|_"; break;
        case Top:        o << "T"; break;
        case NonNull:    o << "NN"; break;
        default:/*Null*/ o << "N";  
      }
    }

  }; // end class nullity_value


  // Abstract domain for nullity
  template <typename Number, typename VariableName>
  class nullity_domain : 
        public ikos::writeable,
        public numerical_domain<Number, VariableName >,
        public bitwise_operators<Number,VariableName >,
        public division_operators<Number, VariableName >,
        public array_operators< Number, VariableName >,
        public pointer_operators< Number, VariableName > {
    
    typedef separate_domain< VariableName, nullity_value > separate_domain_t;
    typedef nullity_domain<Number, VariableName> nullity_domain_t;
    // lin_exp_t == linear_expression_t
    typedef typename pointer_operators <Number, VariableName>::lin_exp_t lin_exp_t;
    typedef typename pointer_operators <Number, VariableName>::ptr_cst_t ptr_cst_t;

   public:

    typedef VariableName varname_t;
    typedef Number number_t;
    typedef linear_expression<Number, VariableName> linear_expression_t;
    typedef linear_constraint<Number, VariableName> linear_constraint_t;
    typedef linear_constraint_system<Number, VariableName> linear_constraint_system_t;
    typedef typename separate_domain_t::iterator iterator;
    
   private:

    separate_domain_t _env;
    
    nullity_domain(separate_domain_t env) : _env(env) {}
        
   public:
    
    static nullity_domain_t top() {
      return nullity_domain(separate_domain_t::top());
    }
    
    static nullity_domain_t bottom() {
      return nullity_domain(separate_domain_t::bottom());
    }
    
    nullity_domain() : _env(separate_domain_t::top()) {}

    nullity_domain(const nullity_domain_t& e) : 
        writeable(), _env(e._env) {}
    
    nullity_domain_t& operator=(const nullity_domain_t& o) {
      if (this != &o)
        _env = o._env;
      return *this;
    }
    
    iterator begin() { 
      if (is_bottom ()) 
        CRAB_ERROR ("Cannot return iterator from bottom");
      
      return _env.begin(); 
    }
    
    iterator end() { 
      if (is_bottom ()) 
        CRAB_ERROR ("Cannot return iterator from bottom");
      return _env.end(); 
    }
    
    bool is_bottom() { return _env.is_bottom(); }
    
    bool is_top() { return _env.is_top(); }
    
    bool operator<=(nullity_domain_t o) { return (_env <= o._env); }
    
    nullity_domain_t operator|(nullity_domain_t o) {
      nullity_domain_t res (_env | o._env);

      CRAB_LOG("nullity", crab::outs () << "After join " << *this << " and " << o << "=" << res << "\n";);
      return res;
    }

    void operator|=(nullity_domain_t o) {
      CRAB_LOG("nullity", crab::outs () << "After join " << *this << " and " << o << "=");

      _env = _env  | o._env;

      CRAB_LOG("nullity", crab::outs () << *this << "\n");
    }

    nullity_domain_t operator&(nullity_domain_t o) {
      nullity_domain_t res (_env & o._env);
      
      CRAB_LOG("nullity", crab::outs () << "After meet " << *this << " and " << o << "=" << res << "\n");
      return res;
    }
    
    nullity_domain_t operator||(nullity_domain_t o) {
      nullity_domain_t res (_env || o._env);

      CRAB_LOG("nullity", crab::outs () << "After widening " << *this << " and " << o << "=" << res << "\n");
      return res;
    }

    template<typename Thresholds>
    nullity_domain_t widening_thresholds (nullity_domain_t o, const Thresholds &) {
      nullity_domain_t res (_env || o._env);

      CRAB_LOG("nullity", crab::outs () << "After widening " << *this << " and " << o << "=" << res << "\n");
      return res;
    }
    
    nullity_domain_t operator&&(nullity_domain_t o) {
      return (_env && o._env);
    }
    
    void set(VariableName v, nullity_value n) {
      if (!is_bottom())
        _env.set(v, n);
    }

    void assign(VariableName x, VariableName y) {
      if (!is_bottom())
        _env.set (x, _env[y]);

      CRAB_LOG("nullity", crab::outs () << "After " << x << ":="  << y << "=" << *this <<"\n");
    }
    
    nullity_value operator[](VariableName v) { 
      return _env[v]; 
    }
    
    void operator-=(VariableName v) { 
      if (!is_bottom ())
        _env -= v; 
    }
        
    void equality (VariableName p, VariableName q) {
      if (is_bottom ()) return;

      // if (p == q) ...
      nullity_value p_meet_q = _env[p] & _env[q];
      _env.set(p, p_meet_q);
      _env.set(q, p_meet_q);

      CRAB_LOG("nullity", crab::outs () << "After " << p << "=="  << q << "=" << *this <<"\n");
    }

    void equality (VariableName p, nullity_value v) {
      if (!is_bottom ()) 
        _env.set(p, _env[p] & v);

      CRAB_LOG("nullity", crab::outs () << "After " << p << "=="  << v << "=" << *this <<"\n");
    }
    
    void disequality (VariableName p, VariableName q) {
      if (is_bottom ()) return;

      // if (p != q) ...
      if (_env[p].is_null() && _env[q].is_null()) {
        *this = bottom();
      } else if (_env[p].is_top() && _env[q].is_null()) {
        _env.set(p, nullity_value::non_null()); // refine p
      } else if (_env[q].is_top() && _env[p].is_null()) {
        _env.set(q, nullity_value::non_null()); // refine q
      }

      CRAB_LOG("nullity", crab::outs () << "After " << p << "!="  << q << "=" << *this <<"\n");
    }
    
    void disequality (VariableName p, nullity_value v) {
      if (is_bottom ()) return;

      if (_env[p].is_null() && v.is_null()) {
        *this = bottom();
      } else if (_env[p].is_top() && v.is_null()) { // refine p
        _env.set(p, nullity_value::non_null());
      }

      CRAB_LOG("nullity", crab::outs () << "After " << p << "!="  << v << "=" << *this <<"\n");
    }

    // numerical_domains_api
    // XXX: needed for making a reduced product with a numerical domain
    void apply(operation_t op, VariableName x, VariableName y, VariableName z) {}
    void apply(operation_t op, VariableName x, VariableName y, Number k) {}
    void assign(VariableName x, linear_expression_t e) {}
    void operator+=(linear_constraint_system_t csts) {}
    void operator+=(linear_constraint_t cst) {}

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

    // pointer_operators_api 
    virtual void pointer_load (VariableName /*lhs*/, VariableName rhs) override {
      // XXX: assume after the load the rhs must be non-null otherwise
      // the program failed.
      equality (rhs, nullity_value::non_null ());
    }

    virtual void pointer_store (VariableName lhs, VariableName /*rhs*/) override {
      // XXX: assume after the store the lhs must be non-null
      // otherwise the program failed.
      equality (lhs, nullity_value::non_null ());
    } 

    virtual void pointer_assign (VariableName lhs, VariableName rhs, lin_exp_t /*offset*/) override {
      assign (lhs, rhs);
    }

    virtual void pointer_mk_obj (VariableName lhs, ikos::index_t /*address*/) override {
      set (lhs, nullity_value::non_null ());
    }

    virtual void pointer_function (VariableName lhs, VariableName /*func*/) override {
      set (lhs, nullity_value::non_null ());
    }
    
    virtual void pointer_mk_null (VariableName lhs) override {
      set (lhs, nullity_value::null ());
    }
    
    virtual void pointer_assume (ptr_cst_t cst) override {
      if (cst.is_tautology ()) return;
      
      if (cst.is_contradiction ()) {
        *this = bottom();
        return;
      }

      if (cst.is_unary ()) {
        if (cst.is_equality ()) 
          equality (cst.lhs (), nullity_value::null ());
        else // cst.is_disequality ();
          disequality (cst.lhs (), nullity_value::null ());
      } else { 
        assert (cst.is_binary ());
        if (cst.is_equality ()) 
          equality (cst.lhs (), cst.rhs ());
        else  // cst.is_disequality ();
          disequality (cst.lhs (), cst.rhs ());
      }
    }
    
    virtual void pointer_assert (ptr_cst_t cst) override {
      CRAB_WARN ("nullity pointer_assert not implemented");
    }
        
    static std::string getDomainName () {
      return "Nullity";
    }    

    void write(crab_os& o) {
      _env.write(o); 
    }
    
    }; // class nullity_domain
 
  } // end namespace domains 
} // end namespace crab

#endif 
