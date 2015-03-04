/*******************************************************************************
 *
 * Anti-unification domain -- lifting a value domain 
 * using term equivalences.
 *
 * Author: Graeme Gange (gkgange@unimelb.edu.au)
 ******************************************************************************/

#ifndef IKOS_ANTI_UNIF_HPP
#define IKOS_ANTI_UNIF_HPP

/**

 **/
#include <utility>
#include <algorithm>
#include <iostream>
#include <vector>
#include <sstream>
#include <ikos_domains/common.hpp>
#include <ikos_domains/bignums.hpp>
#include <ikos_domains/linear_constraints.hpp>
#include <ikos_domains/numerical_domains_api.hpp>
#include <ikos_domains/bitwise_operators_api.hpp>
#include <ikos_domains/division_operators_api.hpp>
#include <ikos_domains/intervals.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/optional.hpp>

#define BAIL(x) assert(0 && (x))

using namespace boost;
using namespace std;

// #define VERBOSE 

namespace ikos {
  
  template< typename Info >
  class anti_unif: public writeable,
                 public numerical_domain<typename Info::Number, typename Info::VariableName >, 
                 public bitwise_operators< typename Info::Number, typename Info::VariableName >,
                 public division_operators< typename Info::Number, typename Info::VariableName > {
  private:
    // Underlying (value?) domain.
    typedef typename Info::Number Number;
    typedef typename Info::VariableName VariableName;
    typedef typename Info::domain_t  dom_t;

    typedef bound<Number> bound_t;
    // typedef interval_domain< Number, VariableName > intervals_t;
    
    class term_table: public writeable{
    private:
      // Uses a single state of the underlying domain.
      dom_t impl;
      
    public:
      term_table(): writeable() { }	// Assumption: dom_t() yields top.
      
      term_table(const term_table &other): writeable(), impl(other.impl) { }
      
      term_table operator=(term_table other){
        this->impl = other.impl;
        return *this;
      }
      
      void resize(size_t new_size){	// Expected size value to represent number of contained variables.
        impl.resize(new_size);
      }

      // downsize the matrix
      void operator-=(unsigned int k){
        impl -= k;
      }

      size_t size(){	// Returned size indicates number of variables.
        return impl->size();
      }
      
      /*
      bound_t& operator()(unsigned int i, unsigned int j){	
        if(i < 1 || j < 1 || i > 2*this->num_var || j > 2*this->num_var){
          stringstream err;
          err << "OCTAGON: out-of-bounds error accessing internal octagon matrix at [" 
              << i << ", " << j << "] where 2n = " << 2*num_var << ".";
          throw error(err.str());
        }
        return this->matrix[2*this->num_var*(j- 1) + (i- 1)];
      }
      */

      ostream& write(ostream& o) { 
        /*
        for (unsigned int i=1; i <= 2*this->num_var; i++){
          for (unsigned int j=1; j <= 2*this->num_var; j++){
            bound_t val = this->operator()(i,j);
            o << "M[" << i << ","  <<  j << "]=" << val << ";";
          }
          o << endl;
        }
        */
        return o;
      }
      
    }; // end class term_table
    
   public:
    typedef variable< Number, VariableName >                 variable_t;
    typedef linear_constraint< Number, VariableName >        linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    typedef linear_expression< Number, VariableName >        linear_expression_t;
    typedef anti_unif<Info>        anti_unif_t;
    typedef interval< Number >                               interval_t;
    typedef interval_domain< Number, VariableName >          interval_domain_t;

    typedef patricia_tree_set< VariableName >  varname_set_t;
     
   private:
    typedef term_table ttbl_t;
    typedef container::flat_map< variable_t, size_t > map_t;
    typedef typename map_t::value_type value_type;
    typedef unsigned char BOOL;

   private:
    bool     _is_bottom;
    ttbl_t   _term_table;
    map_t    _map;
    bool     _is_normalized;
    vector< BOOL > _norm_vector; // IMPORTANT: Treat this as a vector of booleans.

    void set_to_bottom (){
      this->_is_bottom = true;
    }

   private:
    anti_unif(bool is_top): _is_bottom(!is_top), _is_normalized(true) { }

#if 0
    void add_var(VariableName x, unsigned int i, unsigned int j, bound_t lb, bound_t ub, bool op_eq){
      BAIL("ANTI_UNIF: add_var not yet implemented.");
      /*
      if (lb.is_minus_infinity() && ub.is_plus_infinity()){
        abstract(x);
        return;
      }

      if(op_eq){
        for(size_t j_idx = 1; j_idx <= 2*_dbm.size(); ++j_idx){
          if(j_idx!= 2*j && j_idx!= 2*j- 1){
            _dbm(2*j- 1, j_idx)-= lb;
            _dbm(2*j, j_idx)+= ub;
          }
        }
        for(size_t i_idx = 1; i_idx <= 2*_dbm.size(); ++i_idx){
          if(i_idx!= 2*j && i_idx!= 2*j- 1){
            _dbm(i_idx, 2*j)-= lb;
            _dbm(i_idx, 2*j- 1)+= ub;
          }
        }
        _dbm(2*j- 1, 2*j)-= bound_t(2)*lb;
        _dbm(2*j, 2*j- 1)+= bound_t(2)*ub;
      }
      else {
        abstract(x);
        _dbm(2*j- 1, 2*i- 1)= -lb;
        _dbm(2*i, 2*j)= -lb;
        _dbm(2*i- 1, 2*j- 1)= ub;
        _dbm(2*j, 2*i)= ub;
      }
      */
    }
    
    void subtract_var(VariableName x, unsigned int i, unsigned int j, bound_t lb, bound_t ub, bool op_eq){
      BAIL("ANTI-UNIF: subtract_var not yet implemented.");
      // add_var(x, i, j, -ub, -lb, op_eq);
    }
    
    void multiply_var(VariableName x, unsigned int i, unsigned int j, bound_t lb, bound_t ub, bool op_eq){
      BAIL("ANTI-UNIF: multiply_var not yet implemented.");
      /*
      if(op_eq){
        bound_t t1(_dbm(2*j- 1, 2*j)), t2(_dbm(2*j, 2*j- 1));
        abstract(x);
        bound_t ll = (t1/-2) * lb;
        bound_t lu = (t1/-2) * ub;
        bound_t ul = (t2/2)  * lb;
        bound_t uu = (t2/2)  * ub;
        _dbm(2*j- 1, 2*j) = bound_t::min(ll, lu, ul, uu) * -2;
        _dbm(2*j, 2*j- 1) = bound_t::max(ll, lu, ul, uu) *  2;        
      }
      else {
        bound_t t1(_dbm(2*i- 1, 2*i).operator-()), t2(_dbm(2*i, 2*i- 1));
        abstract(x);
        bound_t ll = (t1/-2) * lb;
        bound_t lu = (t1/-2) * ub;
        bound_t ul = (t2/2)  * lb;
        bound_t uu = (t2/2)  * ub;
        _dbm(2*j- 1, 2*j) = bound_t::min(ll, lu, ul, uu) * -2;
        _dbm(2*j, 2*j- 1) = bound_t::max(ll, lu, ul, uu) *  2;        
      }
      */
    }
    
    void divide_var(VariableName x, unsigned int i, unsigned int j, bound_t _lb, bound_t _ub, bool op_eq){
      BAIL("ANTI-UNIF: divide_var not yet implemented.");
      /*
      interval_t trim_intv = ikos::intervals_impl::trim_bound(interval_t(_lb,_ub), Number(0));

      if (trim_intv.is_bottom ()){
        // definite division by zero
        set_to_bottom ();
        return;
      }
        
      interval_t zero(Number(0));
      if (zero <= trim_intv){
        abstract(x);
        return;
      }
      bound_t lb = trim_intv.lb();
      bound_t ub = trim_intv.ub();
      if(op_eq){
        bound_t t1(_dbm(2*j- 1, 2*j)), t2(_dbm(2*j, 2*j- 1));
        abstract(x);
        bound_t ll = (t1/-2) / lb;
        bound_t lu = (t1/-2) / ub;
        bound_t ul = (t2/2)  / lb;
        bound_t uu = (t2/2)  / ub;
        _dbm(2*j- 1, 2*j) = bound_t::min(ll, lu, ul, uu) * -2;
        _dbm(2*j, 2*j- 1) = bound_t::max(ll, lu, ul, uu) *  2;        
      }
      else {
        bound_t t1(_dbm(2*i- 1, 2*i).operator-()), t2(_dbm(2*i, 2*i- 1));
        abstract(x);
        bound_t ll = (t1/-2) / lb;
        bound_t lu = (t1/-2) / ub;
        bound_t ul = (t2/2)  / lb;
        bound_t uu = (t2/2)  / ub;
        _dbm(2*j- 1, 2*j) = bound_t::min(ll, lu, ul, uu) * -2;
        _dbm(2*j, 2*j- 1) = bound_t::max(ll, lu, ul, uu) *  2;        
      }
      */
    }
#endif

    // x = y op [lb,ub]
    void apply(operation_t op, VariableName x, VariableName y, bound_t lb, bound_t ub){	
      BAIL("ANTI-UNIF: apply not yet implemented.");

      /*
      // Requires normalization.

      // add x in the DBM if not found
      if(this->_map.find(x) == this->_map.end()){
        this->_map.insert(value_type(x,this->_map.size()+ 1));
        this->resize();
      }
      
      if (this->_map.find(y) == this->_map.end())
      {
        this->abstract(x);
        return;
      }

      unsigned int i(_map.find(x)->second), j(_map.find(y)->second);
      normalize();
      
      switch (op) {
        case OP_ADDITION:
          add_var(x, i, j, lb, ub, i==j);
          break;
        case OP_SUBTRACTION:
          subtract_var(x, i, j, lb, ub, i==j);
          break;
        case OP_MULTIPLICATION:
          multiply_var(x, i, j, lb, ub, i==j);
          break;
        case OP_DIVISION:
          divide_var(x, i, j, lb, ub, i==j);
          break;
        default:
          throw error("OCTAGON: unsupported arithmetic operation.");
      }
      _norm_vector.at(i- 1)= 0;
      // Result is not normalized.
      */
    }
    
    void apply_constraint(unsigned int var, bool is_positive, bound_t constraint){	
      BAIL("ANTI-UNIF: apply_constraint not yet implemented.");
    }

    void apply_constraint(unsigned int i, unsigned int j, bool is1_positive, 
                          bool is2_positive, bound_t constraint){	
      BAIL("ANTI-UNIF: apply_constraint not yet implemented.");
    }

    // check satisfiability of cst using intervals
    // Only to be used if cst is too hard for octagons
    bool check_sat(linear_constraint_t cst)  {
      BAIL("ANTI-UNIF: check_sat not yet implemented.");
      /*
      typename linear_constraint_t::variable_set_t vars = cst.variables();
      intervals_t inv;
      normalize();
      for(typename linear_constraint_t::variable_set_t::iterator it = vars.begin(); it!= vars.end(); ++it){
        inv.set (it->name(), to_interval(it->name(), false));
      }
      inv += cst;
      */
      return (!_is_bottom);
    }

    interval_t to_interval(VariableName x, bool requires_normalization) { 
      // projection requires normalization.
      BAIL("ANTI-UNIF: to_interval not yet implemented.");
      typename map_t::iterator it(_map.find(x));
      if(it == _map.end()){
        return interval_t::top();
      }
      else{
        return interval_t::top();
        /*
        unsigned int idx(it->second);
        if (requires_normalization)
          normalize();
        return interval_t(_dbm(2*idx- 1, 2*idx).operator/(-2), _dbm(2*idx, 2*idx- 1).operator/(2));
        */
      }
    } //Maintains normalization.
    
    void resize(){
      // _dbm.resize(_map.size());
      _norm_vector.resize(_map.size(), 0);
    }
    
    void is_normalized(bool b){
      _is_normalized= b;
      for(vector<unsigned char>::iterator it = _norm_vector.begin(); it!= _norm_vector.end() ; ++it){
        (*it)=(b ? 1 : 0);
      }
    }

    
  public:
    static anti_unif_t top() {
      return anti_unif(true);
    }
    
    static anti_unif_t bottom() {
      return anti_unif(false);
    }
    
  public:
    // Constructs top octagon, represented by a size of 0.
    anti_unif(): _is_bottom(false), _is_normalized(true) { }
    
    anti_unif(const anti_unif_t& o): 
       writeable(), 
       numerical_domain<Number, VariableName >(),
       bitwise_operators< Number, VariableName >(),
       division_operators< Number, VariableName >(),   
       _is_bottom(o._is_bottom), 
       _map(o._map), _is_normalized(o._is_normalized), 
       _norm_vector(o._norm_vector) { }
    
    anti_unif_t operator=(anti_unif_t o) {
      _is_bottom= o.is_bottom();
      _map= o._map;
      _is_normalized= o._is_normalized;
      _norm_vector= o._norm_vector;
      return *this;
    }
    
    bool is_bottom() {
      return _is_bottom;
    }
    
    bool is_top() {
      return !_map.size() && !is_bottom();
    }
    
    bool is_normalized(){
      return _is_normalized;
    }

    varname_set_t get_variables() const {
      varname_set_t vars;
      for (typename map_t::const_iterator it = this->_map.begin(), et=this->_map.end(); it!=et; ++it){
        variable_t v(it->first);
        vars += v.name();
      }
      return vars;
    }

    /*
    inline bound_t C(const bound_t &a, const bound_t &b, const bound_t &c, const bound_t &d, const bound_t &e){
      // Separating the bound_t::min()s yields 20-25% increases in performance.
      return bound_t::min(a, bound_t::min(b, bound_t::min(c, bound_t::min(d,e))));
    }
    */

    // Compute the strong closure algorithm
    void normalize(){
      if(_is_normalized){
        return;
      }

      fprintf(stderr, "WARNING: ANTI_UNIF::normalize not yet implemented.\n");

      /*
      is_normalized(false);

      size_t num_var(_dbm.size());
      bound_t zero(0);
      for(size_t k = 1; k <= num_var; ++k){
        if(_norm_vector.at(k- 1)== 0){
          for(size_t i = 1; i <= 2*num_var; ++i){
            for(size_t j = 1; j <= 2*num_var; ++j){  
              // to ensure the "closed" property
              _dbm(i, j)= C(_dbm(i, j), 
                            _dbm(i, 2*k- 1) + _dbm(2*k- 1, j),
                            _dbm(i, 2*k)    + _dbm(2*k, j),
                            _dbm(i, 2*k- 1) + _dbm(2*k- 1, 2*k) + _dbm(2*k, j),
                            _dbm(i, 2*k)    + _dbm(2*k, 2*k- 1) + _dbm(2*k- 1, j));                             
            }
          }
          // to ensure for all i,j: m_ij <= (m_i+i- + m_j-j+)/2
          for(size_t i = 1; i <= 2*num_var; ++i){
            for(size_t j = 1; j <= 2*num_var; ++j)
              _dbm(i, j)= bound_t::min(_dbm(i, j), 
                                       (_dbm(i, i+ 2*(i%2)- 1) + _dbm(j + 2*(j%2)- 1, j))/Number(2) );
          }
          _norm_vector.at(k- 1)= 1;
        }
      }
      // negative cycle?
      for(size_t i = 1; i <= 2*num_var; ++i){
        if(_dbm(i, i) < 0){
          this->operator=(bottom());
          return;
        }
        _dbm(i, i)= zero;
      }
      _is_normalized= true;
      */
    }

    // Lattice operations
    bool operator<=(anti_unif_t o)  {	
      // Require normalization of the first argument
      this->normalize();

      if (is_bottom()) {
        return true;
      } else if(o.is_bottom()) {
        return false;
      } else {
        BAIL("ANTI-UNIF: <= not yet implemented.");
        return false;
        /*
        map_t _temp; // Used to construct a list of variables that appear in both octagons.
        unsigned int i1, j1, i2, j2;
        for(typename map_t::iterator ito = o._map.begin(); ito!= o._map.end(); ++ito) {
          i2= ito->second;
          if(_map.find(ito->first)== _map.end()) {
            if(!o._dbm(2*i2- 1, 2*i2).is_infinite() || !o._dbm(2*i2, 2*i2- 1).is_infinite()) {	
              // Case: Variable exists and is finite in _other but
              // does not exist in _this.
              return false;
            }
          } else {
            _temp.insert(value_type(ito->first, 0));
          }
        }
        
        for(typename map_t::iterator it=_temp.begin(); it!= _temp.end(); ++it) {	
          // Case: Variable exists in both _this and _other.
          i1= _map.find(it->first)->second;
          i2= o._map.find(it->first)->second;
          
          if(!(_dbm(2*i1- 1, 2*i1- 1) <= o._dbm(2*i2- 1, 2*i2- 1)) || 
             !(_dbm(2*i1- 1, 2*i1)    <= o._dbm(2*i2- 1, 2*i2)) || 
             !(_dbm(2*i1, 2*i1- 1)    <= o._dbm(2*i2, 2*i2- 1)) || 
             !(_dbm(2*i1, 2*i1)       <= o._dbm(2*i2, 2*i2))){
            return false;
          }
          for(typename map_t::iterator it2=it+1; it2!= _temp.end(); ++it2) {
            j1= _map.find(it2->first)->second;
            j2= o._map.find(it2->first)->second;
            
            if(!(_dbm(2*i1- 1, 2*j1- 1) <= o._dbm(2*i2- 1, 2*j2- 1)) || 
               !(_dbm(2*i1- 1, 2*j1)    <= o._dbm(2*i2- 1, 2*j2)) || 
               !(_dbm(2*i1, 2*j1- 1)    <= o._dbm(2*i2, 2*j2- 1)) || 
               !(_dbm(2*i1, 2*j1)       <= o._dbm(2*i2, 2*j2))){
              return false;
            }
            if(!(_dbm(2*j1- 1, 2*i1- 1) <= o._dbm(2*j2- 1, 2*i2- 1)) || 
               !(_dbm(2*j1- 1, 2*i1)    <= o._dbm(2*j2- 1, 2*i2)) || 
               !(_dbm(2*j1, 2*i1- 1)    <= o._dbm(2*j2, 2*i2- 1)) || 
               !(_dbm(2*j1, 2*i1)       <= o._dbm(2*j2, 2*i2))){
              return false;
            }
          }
        }
        return true; 
        */
      }
    }  // Maintains normalization.

    struct join_op{
      bound_t operator()(bound_t v1, bound_t v2) {
        return bound_t::max(v1,v2);
      }
    };

    struct widening_op{
      bound_t operator()(bound_t v1, bound_t v2) {
        bound_t infinite("+oo");
        return (v1 >= v2) ? v1 : infinite;
      }
    };
    
    template<typename op_t>
    anti_unif_t pointwise_binary_op(anti_unif_t o1, anti_unif_t o2){
      BAIL("ANTI-UNIF: pointwise_binary_op not yet implemented.");
      return top();
      /*
      octagon_t n;
      // Set intersection of the two maps
      for(typename map_t::iterator it=o1._map.begin(); it!= o1._map.end(); ++it){
        if(o2._map.find(it->first)!= o2._map.end()){
          n._map.insert(value_type(it->first, n._map.size()+ 1));
        }
      }
      if(n._map.size()== 0){
        return top();
      }
      n.resize();
      n.is_normalized(true);
      
      unsigned int i1, i2, i3, j1, j2, j3;
      i1=i2=i3=j1=j2=j3=0;
      for(typename map_t::iterator it=n._map.begin(); it!= n._map.end(); ++it){	
        // Finds the union of each 2x2 identity matrix.
        i1= o1._map.find(it->first)->second;
        i2= o2._map.find(it->first)->second;
        i3= it->second;
        
        op_t op;
        n._dbm(2*i3- 1, 2*i3- 1)= op(o1._dbm(2*i1- 1, 2*i1- 1), o2._dbm(2*i2- 1, 2*i2- 1));
        n._dbm(2*i3- 1, 2*i3)   = op(o1._dbm(2*i1- 1, 2*i1)   , o2._dbm(2*i2- 1, 2*i2));
        n._dbm(2*i3, 2*i3- 1)   = op(o1._dbm(2*i1, 2*i1- 1)   , o2._dbm(2*i2, 2*i2- 1));
        n._dbm(2*i3, 2*i3)      = op(o1._dbm(2*i1, 2*i1)      , o2._dbm(2*i2, 2*i2));
        
        for(typename map_t::iterator it2=it+1; it2!= n._map.end(); ++it2){	
          // Finds the union of each pair of 2x2 relational matrices.
          j1= o1._map.find(it2->first)->second;
          j2= o2._map.find(it2->first)->second;
          j3= it2->second;
          
          n._dbm(2*i3- 1, 2*j3- 1)= op(o1._dbm(2*i1- 1, 2*j1- 1),o2._dbm(2*i2- 1, 2*j2- 1));
          n._dbm(2*i3- 1, 2*j3)   = op(o1._dbm(2*i1- 1, 2*j1)   ,o2._dbm(2*i2- 1, 2*j2));
          n._dbm(2*i3, 2*j3- 1)   = op(o1._dbm(2*i1, 2*j1- 1)   ,o2._dbm(2*i2, 2*j2- 1));
          n._dbm(2*i3, 2*j3)      = op(o1._dbm(2*i1, 2*j1)      ,o2._dbm(2*i2, 2*j2));
          
          n._dbm(2*j3- 1, 2*i3- 1)= op(o1._dbm(2*j1- 1, 2*i1- 1),o2._dbm(2*j2- 1, 2*i2- 1));
          n._dbm(2*j3- 1, 2*i3)   = op(o1._dbm(2*j1- 1, 2*i1)   ,o2._dbm(2*j2- 1, 2*i2));
          n._dbm(2*j3, 2*i3- 1)   = op(o1._dbm(2*j1, 2*i1- 1)   ,o2._dbm(2*j2, 2*i2- 1));
          n._dbm(2*j3, 2*i3)      = op(o1._dbm(2*j1, 2*i1)      ,o2._dbm(2*j2, 2*i2));
        }
      }
      return n;
      */
    }

    // Join
    anti_unif_t operator|(anti_unif_t o) { 
      // Requires normalization of both operands
      normalize();
      o.normalize();
      if (is_bottom()) {
        return o;
      } 
      else if(o.is_bottom()) {
        return *this;
      } 
      else {
        return pointwise_binary_op<join_op >(*this, o);
      }
    } // Returned matrix is normalized.

    // Widening
    anti_unif_t operator||(anti_unif_t o) {	
      // The left operand of the widenning cannot be closed, otherwise
      // termination is not ensured. However, if the right operand is
      // close precision may be improved.
      o.normalize();
      if (is_bottom()) {
        return o;
      } 
      else if(o.is_bottom()) {
        return *this;
      } 
      else {
        return pointwise_binary_op<widening_op >(*this, o);
      }
    } // Returned matrix is not normalized.


    // Meet
    anti_unif_t operator&(anti_unif_t o) { 
      // Does not require normalization of any of the two operands
      if (is_bottom() || o.is_bottom()) {
        return bottom();
      } else {
        BAIL("ANTI-UNIF: meet not yet implemented.");
        return top();
        /*
        octagon_t n;
        // Set union of the two maps
        for(typename map_t::iterator it=_map.begin(); it!= _map.end(); ++it){
          n._map.insert(value_type(it->first, n._map.size()+ 1));
        }
        for(typename map_t::iterator it=o._map.begin(); it!= o._map.end(); ++it){
          n._map.insert(value_type(it->first, n._map.size()+ 1));
        }
        if(n._map.size()== 0){
          return top();
        }
        n.resize();
        n.is_normalized(false);
        
        typename map_t::iterator testi1, testi2, testj1, testj2;
        unsigned int i1, i2, i3, j1, j2, j3;
        i1=i2=i3=j1=j2=j3=0;
        for(typename map_t::iterator it=n._map.begin(); it!= n._map.end(); ++it){	
          // Finds the intersection on each 2x2 identity matrix.
          testi1= _map.find(it->first);
          testi2= o._map.find(it->first);
          i3= it->second;
          
          if(testi1== _map.end()){
            i1= 0;
          } else {
            i1= testi1->second;
          }
          if(testi2== o._map.end()){
            i2= 0;
          } else {
            i2= testi2->second;
          }
          
          n._dbm(2*i3- 1, 2*i3- 1)= (!i1) ? o._dbm(2*i2- 1, 2*i2- 1) : ((!i2) ? _dbm(2*i1- 1, 2*i1- 1) : bound_t::min(_dbm(2*i1- 1, 2*i1- 1), o._dbm(2*i2- 1, 2*i2- 1)));
          n._dbm(2*i3- 1, 2*i3)= (!i1) ? o._dbm(2*i2- 1, 2*i2) : ((!i2) ? _dbm(2*i1- 1, 2*i1) : bound_t::min(_dbm(2*i1- 1, 2*i1),o._dbm(2*i2- 1, 2*i2)));
          n._dbm(2*i3, 2*i3- 1)= (!i1) ? o._dbm(2*i2, 2*i2- 1) : ((!i2) ? _dbm(2*i1, 2*i1- 1) : bound_t::min(_dbm(2*i1, 2*i1- 1),o._dbm(2*i2, 2*i2- 1)));
          n._dbm(2*i3, 2*i3)= (!i1) ? o._dbm(2*i2, 2*i2) : ((!i2) ? _dbm(2*i1, 2*i1) : bound_t::min(_dbm(2*i1, 2*i1),o._dbm(2*i2, 2*i2)));
          

          for(typename map_t::iterator it2=it+1; it2!= n._map.end(); ++it2){	
            // Finds the intersection of each pair of 2x2 relational matrices.
            testj1= _map.find(it2->first);
            testj2= o._map.find(it2->first);
            j3= it2->second;
            
            if(testj1== _map.end()){
              j1= 0;
            } else {
              j1= testj1->second;
            }
            if(testj2== o._map.end()){
              j2= 0;
            } else {
              j2= testj2->second;
            }
            
            if(((i1 && !j1) || (!i1 && j1)) && ((i2 && !j2) || (!i2 && j2))){
              continue;
            }
            if(i1 > j1){
              swap(i1, j1);
            }
            if(i2 > j2){
              swap(i2, j2);
            }
            
            n._dbm(2*i3- 1, 2*j3- 1)= (!i1 || !j1) ? o._dbm(2*i2- 1, 2*j2- 1) : ((!i2 || !j2) ? _dbm(2*i1- 1, 2*j1- 1) : bound_t::min(_dbm(2*i1- 1, 2*j1- 1),o._dbm(2*i2- 1, 2*j2- 1)));
            n._dbm(2*i3- 1, 2*j3)= (!i1 || !j1) ? o._dbm(2*i2- 1, 2*j2) : ((!i2 || !j2) ? _dbm(2*i1- 1, 2*j1) : bound_t::min(_dbm(2*i1- 1, 2*j1),o._dbm(2*i2- 1, 2*j2)));
            n._dbm(2*i3, 2*j3- 1)= (!i1 || !j1) ? o._dbm(2*i2, 2*j2- 1) : ((!i2 || !j2) ? _dbm(2*i1, 2*j1- 1) : bound_t::min(_dbm(2*i1, 2*j1- 1),o._dbm(2*i2, 2*j2- 1)));
            n._dbm(2*i3, 2*j3)= (!i1 || !j1) ? o._dbm(2*i2, 2*j2) : ((!i2 || !j2) ? _dbm(2*i1, 2*j1) : bound_t::min(_dbm(2*i1, 2*j1),o._dbm(2*i2, 2*j2)));
            
            n._dbm(2*j3- 1, 2*i3- 1)= (!i1 || !j1) ? o._dbm(2*j2- 1, 2*i2- 1) : ((!i2 || !j2) ? _dbm(2*j1- 1, 2*i1- 1) : bound_t::min(_dbm(2*j1- 1, 2*i1- 1),o._dbm(2*j2- 1, 2*i2- 1)));
            n._dbm(2*j3- 1, 2*i3)= (!i1 || !j1) ? o._dbm(2*j2- 1, 2*i2) : ((!i2 || !j2) ? _dbm(2*j1- 1, 2*i1) : bound_t::min(_dbm(2*j1- 1, 2*i1),o._dbm(2*j2- 1, 2*i2)));
            n._dbm(2*j3, 2*i3- 1)= (!i1 || !j1) ? o._dbm(2*j2, 2*i2- 1) : ((!i2 || !j2) ? _dbm(2*j1, 2*i1- 1) : bound_t::min(_dbm(2*j1, 2*i1- 1),o._dbm(2*j2, 2*i2- 1)));
            n._dbm(2*j3, 2*i3)= (!i1 || !j1) ? o._dbm(2*j2, 2*i2) : ((!i2 || !j2) ? _dbm(2*j1, 2*i1) : bound_t::min(_dbm(2*j1, 2*i1),o._dbm(2*j2, 2*i2)));
          }
        }
        
        return n;
        */
      }
    }	// Returned matrix is not normalized.
    
    // Narrowing
    anti_unif_t operator&&(anti_unif_t o) {	
      // Does not require normalization of any of the two operands
      if (is_bottom() || o.is_bottom()) {
        return bottom();
      } 
      else {
        BAIL("ANTI-UNIF: narrowing not yet implemented.");
        anti_unif_t n(*this);
        return n;
        /*
        octagon_t n;
        // Set union of the two maps
        for(typename map_t::iterator it=_map.begin(); it!= _map.end(); ++it){
          n._map.insert(value_type(it->first, n._map.size()+ 1));
        }
        for(typename map_t::iterator it=o._map.begin(); it!= o._map.end(); ++it){
          n._map.insert(value_type(it->first, n._map.size()+ 1));
        }
        if(n._map.size()== 0){
          return top();
        }
        n.resize();
        n.is_normalized(false);
        
        typename map_t::iterator testi1, testi2, testj1, testj2;
        unsigned int i1, i2, i3, j1, j2, j3;
        i1=i2=i3=j1=j2=j3=0;
        for(typename map_t::iterator it=n._map.begin(); it!= n._map.end(); ++it){	
          // Finds the narrowing on each 2x2 identity matrix.
          testi1= _map.find(it->first);
          testi2= o._map.find(it->first);
          i3= it->second;
          
          if(testi1== _map.end()){
            i1= 0;
          } else {
            i1= testi1->second;
          }
          if(testi2== o._map.end()){
            i2= 0;
          } else {
            i2= testi2->second;
          }
          
          n._dbm(2*i3- 1, 2*i3- 1)= (!i2) ? _dbm(2*i1- 1, 2*i1- 1) : ((!i1 || _dbm(2*i1- 1, 2*i1- 1).is_infinite()) ? o._dbm(2*i2- 1, 2*i2- 1) : _dbm(2*i1- 1, 2*i1- 1));
          n._dbm(2*i3- 1, 2*i3)= (!i2) ? _dbm(2*i1- 1, 2*i1) : ((!i1 || _dbm(2*i1- 1, 2*i1).is_infinite()) ? o._dbm(2*i2- 1, 2*i2) : _dbm(2*i1- 1, 2*i1));
          n._dbm(2*i3, 2*i3- 1)= (!i2) ? _dbm(2*i1, 2*i1- 1) : ((!i1 || _dbm(2*i1, 2*i1- 1).is_infinite()) ? o._dbm(2*i2, 2*i2- 1) : _dbm(2*i1, 2*i1- 1));
          n._dbm(2*i3, 2*i3)= (!i2) ? _dbm(2*i1, 2*i1) : ((!i1 || _dbm(2*i1, 2*i1).is_infinite()) ? o._dbm(2*i2, 2*i2) : _dbm(2*i1, 2*i1));
          
          for(typename map_t::iterator it2=it+1; it2!= n._map.end(); ++it2){	
            // Finds the narrowing of each pair of 2x2 relational matrices.
            testj1= _map.find(it2->first);
            testj2= o._map.find(it2->first);
            j3= it2->second;
            
            if(testj1== _map.end()){
              j1= 0;
            } else {
              j1= testj1->second;
            }
            if(testj2== o._map.end()){
              j2= 0;
            } else {
              j2= testj2->second;
            }
            
            if(((i1 && !j1) || (!i1 && j1)) && ((i2 && !j2) || (!i2 && j2))){
              continue;
            }
            if(i1 > j1){
              swap(i1, j1);
            }
            if(i2 > j2){
              swap(i2, j2);
            }
            
            n._dbm(2*i3- 1, 2*j3- 1)= (!i2 || !j2) ? _dbm(2*i1- 1, 2*j1- 1) : ((!i1 || !j1 || _dbm(2*i1- 1, 2*j1- 1).is_infinite()) ? o._dbm(2*i2- 1, 2*j2- 1) : _dbm(2*i1- 1, 2*j1- 1));
            n._dbm(2*i3- 1, 2*j3)= (!i2 || !j2) ? _dbm(2*i1- 1, 2*j1) : ((!i1 || !j1 || _dbm(2*i1- 1, 2*j1).is_infinite()) ? o._dbm(2*i2- 1, 2*j2) : _dbm(2*i1- 1, 2*j1));
            n._dbm(2*i3, 2*j3- 1)= (!i2 || !j2) ? _dbm(2*i1, 2*j1- 1) : ((!i1 || !j1 || _dbm(2*i1, 2*j1- 1).is_infinite()) ? o._dbm(2*i2, 2*j2- 1) : _dbm(2*i1, 2*j1- 1));
            n._dbm(2*i3, 2*j3)= (!i2 || !j2) ? _dbm(2*i1, 2*j1) : ((!i1 || !j1 || _dbm(2*i1, 2*j1).is_infinite()) ? o._dbm(2*i2, 2*j2) : _dbm(2*i1, 2*j1));
            
            n._dbm(2*j3- 1, 2*i3- 1)= (!i2 || !j2) ? _dbm(2*j1- 1, 2*i1- 1) : ((!i1 || !j1 || _dbm(2*j1- 1, 2*i1- 1).is_infinite()) ? o._dbm(2*j2- 1, 2*i2- 1) : _dbm(2*j1- 1, 2*i1- 1));
            n._dbm(2*j3- 1, 2*i3)= (!i2 || !j2) ? _dbm(2*j1- 1, 2*i1) : ((!i1 || !j1 || _dbm(2*j1- 1, 2*i1).is_infinite()) ? o._dbm(2*j2- 1, 2*i2) : _dbm(2*j1- 1, 2*i1));
            n._dbm(2*j3, 2*i3- 1)= (!i2 || !j2) ? _dbm(2*j1, 2*i1- 1) : ((!i1 || !j1 || _dbm(2*j1, 2*i1- 1).is_infinite()) ? o._dbm(2*j2, 2*i2- 1) : _dbm(2*j1, 2*i1- 1));
            n._dbm(2*j3, 2*i3)= (!i2 || !j2) ? _dbm(2*j1, 2*i1) : ((!i1 || !j1 || _dbm(2*j1, 2*i1).is_infinite()) ? o._dbm(2*j2, 2*i2) : _dbm(2*j1, 2*i1));
          }
        }
        
        return n;
        */
      }
    }	// Returned matrix is not normalized.

    void operator-=(VariableName v) {
      BAIL("Anti_UNIF::operator -= not yet implemented.");
      /*
      if (boost::optional<typename map_t::iterator> it = this->abstract(v)){
        size_t n = (*it)->second;
        this->_term_table -= n;
        this->_map.erase(*it);
        // update the values in _map
        for (typename map_t::iterator itz = this->_map.begin(); itz != this->_map.end(); ++itz){
          if (itz->second > n)
            this->_map[itz->first]--;
        }
        this->_norm_vector.resize(this->_map.size(), 0);
        this->_is_normalized = false;
      }
      */
    }

    void assign(VariableName x, linear_expression_t e) {
      if (this->is_bottom())
      { 
        return; 
      }
      BAIL("ANTI_UNIF::assign not yet implemented.");
      return;
      /*
      optional<variable_t> v = e.get_variable();
      if (v && ((*v).name() == x))
      {
        return ; 
      }

      // add x in the matrix if not found
      typename map_t::iterator it = this->_map.find(x);
      unsigned int i;
      if (it == this->_map.end()){
        i = this->_map.insert(value_type(x,this->_map.size()+ 1)).first->second;
        this->resize();
      }
      else
        i = it->second;

      this->abstract(x); // call normalize()

      if (e.is_constant()){
        this->apply_constraint(i, true , bound_t(e.constant()));    // adding  x <=  c
        this->apply_constraint(i, false, bound_t(-(e.constant()))); // adding -x <= -c
      }
      else if (v){
        VariableName y = (*v).name();
        typename map_t::iterator itz = this->_map.find(y);        
        if (itz == this->_map.end()){
          return; // x has been already abstracted
        }
        unsigned int j = itz->second;
        this->apply_constraint(i, j, true , false , bound_t(0));
        this->apply_constraint(i, j, false, true  , bound_t(0));       
      }
      else
        throw error("OCTAGON: only supports constant or variable on the rhs of assignment");

      this->is_normalized(false);
      */
    }

    // Apply operations to variables.

    // x = y op z
    void apply(operation_t op, VariableName x, VariableName y, VariableName z){	
      BAIL("ANTI_UNIF::apply not yet implemented");
      return;
      /*
      // Requires normalization.

      typename map_t::iterator itz(_map.find(z));
      if (itz == this->_map.end())
      {
        this->abstract(x);
        return;
      }
      unsigned int n(itz->second);

      if (!(x == y)){
        assign(x, linear_expression_t(y));
        apply(op, x, x, _dbm(2*n- 1, 2*n).operator/(-2), _dbm(2*n, 2*n- 1).operator/(2));
      }
      else{
        apply(op, x, y, _dbm(2*n- 1, 2*n).operator/(-2), _dbm(2*n, 2*n- 1).operator/(2));
      }
      // Sets state to not normalized.
      */
    }	
    
    // x = y op k
    void apply(operation_t op, VariableName x, VariableName y, Number k){	
      BAIL("ANTI_UNIF::apply not yet implemented");
      return;
      /*
      // Requires normalization.

      if (!(x == y)){
        assign(x, linear_expression_t(y));
        apply(op, x, x, bound_t(k), bound_t(k));
      }
      else{
        apply(op, x, y, bound_t(k), bound_t(k));
      }
      // Sets state to not normalized.
      */
    }	

    // Add a linear constraint
    void operator+=(linear_constraint_t cst) {  
      BAIL("ANTI_UNIF::operator+= not yet implemented");
      return;
#if 0
      // Does not require normalization.
      if(this->is_bottom()){
        return;
      }
      bool v1,v2, is1_positive, is2_positive;
      v1= v2 = is1_positive = is2_positive = false;
      unsigned int i,j;
      i=j=0;
      for(typename linear_expression_t::iterator it=cst.expression().begin(); it!= cst.expression().end(); ++it){
        if(!v1) {	
          // Calculates and loads information for the first variable.
          if(it->first== -1){
            is1_positive= false;
          } 
          else if(it->first== 1){
            is1_positive= true;
          } 
          else 
            throw error("OCTAGON: expr contains unexpected coefficient (accepted values are -1, 0, and 1).");
          i= _map.insert(value_type(it->second,_map.size()+ 1)).first->second;
          v1= true;
        }
        else if(!v2) {	
          // Calculates and loads information for the second variable,
          // if it exists.
          if(it->first== -1){
            is2_positive= false;
          } 
          else if(it->first== 1){
            is2_positive= true;
          } 
          else 
            throw error("OCTAGON: expr contains unexpected coefficient (accepted values are -1, 0, and 1).");
          j= _map.insert(value_type(it->second,_map.size()+ 1)).first->second;
          v2= true;
        }
        else {
          ostringstream buf;
          buf << "OCTAGON: " << cst << " is not an octagon constraint (> 2 variables).";
          throw error(buf.str());
        }
      }
      if(!v1){
        if((cst.is_inequality() && cst.constant() >= 0) || (cst.is_equality() && cst.constant()== 0)){
          return;
        }
        ostringstream buf;
        buf << "OCTAGON: " << cst << " contains no variables.";
        throw error(buf.str());
      }
      resize();
      bound_t constant(cst.constant()), neg_constant(-(cst.constant()));
      
      if(cst.is_inequality()){	// Applies inequality constraints in the form of octagonal constraints.
        if(v1 && !v2){
          apply_constraint(i, is1_positive, constant);
        } 
        else /*if(v1 && v2)*/{
          apply_constraint(i, j, is1_positive, is2_positive, constant);
        }
      } 
      else if(cst.is_equality()){	// Applies equality constraints as two octagonal constraints.
        if(v1 && !v2){
          apply_constraint(i, is1_positive, constant);
          apply_constraint(i, !is1_positive, neg_constant);
        } 
        else /*if(v1 && v2)*/{
          apply_constraint(i, j, is1_positive, is2_positive, constant);
          apply_constraint(i, j, !is1_positive, !is2_positive, neg_constant);
        }
      }
      else if (cst.is_disequation()){
        // we use intervals to reason about disequations
        if (!check_sat(cst))
          this->operator=(bottom());
      }
      _is_normalized= false;
#endif
    } // Sets state to not normalized.
    
    // Add a system of linear constraints
    void operator+=(linear_constraint_system_t cst) {  // Does not require normalization.
      for(typename linear_constraint_system_t::iterator it=cst.begin(); it!= cst.end(); ++it){
        this->operator+=(*it);
      }
    } // Sets state to not normalized.
    
    // abstract the variable
    // GKG: Looks like this returns the index of variable v.
    boost::optional<typename map_t::iterator> abstract(variable_t v) {	
      // Requires normalization.
      BAIL("ANTI-UNIF: abstract not yet fully implemented.");
      typename map_t::iterator it(_map.find(v));
      if(it!= _map.end()) {
        normalize();
        /*
        size_t n = it->second;
        size_t odd  = 2*n- 1;
        size_t even = 2*n;
        bound_t infinite("+oo"), zero(0);
        size_t size = 2*_dbm.size();
        for(unsigned int idx=1; idx <= size; ++idx){
          _dbm(idx, odd)  = infinite;
          _dbm(idx, even) = infinite;
          _dbm(odd, idx)  = infinite;
          _dbm(even, idx) = infinite;
        }
        _dbm(odd, odd)= zero;
        _dbm(even, even)= zero;
        */
        return boost::optional<typename map_t::iterator>(it);
      }
      return boost::optional<typename map_t::iterator>();
    }  //Maintains normalization.
    
    interval_t operator[](VariableName x) { 
      return to_interval(x, true);
    } 

    void set(VariableName x, interval_t intv){
      // add x in the matrix if not found
      BAIL("ANTI-UNIF: set not yet implemented.");
      /*
      typename map_t::iterator it = this->_map.find(x);
      unsigned int idx;
      if(it == this->_map.end()){
        idx = this->_map.insert(value_type(x,this->_map.size()+ 1)).first->second;
        this->resize();
      }
      else{
        idx = it->second;
      }
      this->abstract(x);  // normalize 
      this->apply_constraint(idx, true , intv.ub());  // x <= ub
      this->apply_constraint(idx, false, -intv.lb()); // -x <= -lb
      */
    }


    interval_domain_t to_intervals(){
      // Requires normalization.
      if(this->_is_bottom){
        return interval_domain_t::bottom();
      }
      else{
        interval_domain_t itv = interval_domain_t::top();
        // we normalize just once
        normalize();
        for(typename map_t::iterator it=_map.begin(); it!= _map.end(); ++it){
          itv.set(it->first.name(), to_interval(it->first.name(), false));
        }
        return itv;
      }
    }

    // bitwise_operators_api
    
    void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width){
      // since reasoning about infinite precision we simply assign and
      // ignore the width.
      assign(x, linear_expression_t(y));
    }

    void apply(conv_operation_t op, VariableName x, Number k, unsigned width){
      // since reasoning about infinite precision we simply assign
      // and ignore the width.
      assign(x, k);
    }

    void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z){
      // Convert to intervals and perform the operation
      interval_t yi = this->operator[](y);
      interval_t zi = this->operator[](z);
      interval_t xi = interval_t::bottom();
      switch (op) {
        case OP_AND: {
	xi = yi.And(zi);
	break;
        }
        case OP_OR: {
	xi = yi.Or(zi);
	break;
        }
        case OP_XOR: {
	xi = yi.Xor(zi);
	break;
        }
        case OP_SHL: {
	xi = yi.Shl(zi);
	break;
        }
        case OP_LSHR: {
          xi = yi.LShr(zi);
          break;
        }
        case OP_ASHR: {
	xi = yi.AShr(zi);
	break;
        }
        default: 
          throw error("OCTAGON: unreachable");
      }
      this->set(x, xi);
    }
    
    void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k){
      /*
      // Convert to intervals and perform the operation
      interval_t yi = this->operator[](y);
      interval_t zi(k);
      interval_t xi = interval_t::bottom();

      switch (op) {
        case OP_AND: {
	xi = yi.And(zi);
	break;
        }
        case OP_OR: {
	xi = yi.Or(zi);
	break;
        }
        case OP_XOR: {
	xi = yi.Xor(zi);
	break;
        }
        case OP_SHL: {
	xi = yi.Shl(zi);
	break;
        }
        case OP_LSHR: {
          xi = yi.LShr(zi);
          break;
        }
        case OP_ASHR: {
	xi = yi.AShr(zi);
	break;
        }
        default: 
          throw error("OCTAGON: unreachable");
      }
      this->set(x, xi);
      */
      BAIL("ANTI-UNIF: bitwise operations not yet implemented.");
    }
    
    // division_operators_api
    
    void apply(div_operation_t op, VariableName x, VariableName y, VariableName z){
      /*
      if (op == OP_SDIV){
        apply(OP_DIVISION, x, y, z);
      }
      else{
        // Convert to intervals and perform the operation
        interval_t yi = this->operator[](y);
        interval_t zi = this->operator[](z);
        interval_t xi = interval_t::bottom();
      
        switch (op) {
          case OP_UDIV: {
            xi = yi.UDiv(zi);
            break;
          }
          case OP_SREM: {
            xi = yi.SRem(zi);
            break;
          }
          case OP_UREM: {
            xi = yi.URem(zi);
            break;
          }
          default: 
            throw error("ANTI-UNIF: unreachable");
        }
        this->set(x, xi);
      }
      */
      BAIL("ANTI-UNIF: div not yet implemented.");
    }

    void apply(div_operation_t op, VariableName x, VariableName y, Number k){
      /*
      if (op == OP_SDIV)
      {
        apply(OP_DIVISION, x, y, k);
      }
      else{
        // Convert to intervals and perform the operation
        interval_t yi = this->operator[](y);
        interval_t zi(k);
        interval_t xi = interval_t::bottom();
      
        switch (op) {
          case OP_UDIV: {
            xi = yi.UDiv(zi);
            break;
          }
          case OP_SREM: {
            xi = yi.SRem(zi);
            break;
          }
          case OP_UREM: {
            xi = yi.URem(zi);
            break;
          }
          default: 
            throw error("ANTI-UNIF: unreachable");
        }
        this->set(x, xi);
      }
      */
      BAIL("ANTI-UNIF: div not yet implemented.");
    }
    
    // Output function
    ostream& write(ostream& o) { 

      // Normalization is not enforced in order to maintain accuracy
      // but we force it to display all the relationships.
      normalize();

      if(_is_bottom){
        return o << "_|_";
      }
      if(_map.size()== 0) {
        return o << "{}";
      }      

      o << "{" ;
      o << "}";
     
#ifdef VERBOSE
      /// For debugging purposes     
      { // print intervals
        interval_domain_t intervals = to_intervals();
        cout << intervals;
      }
      { // print internal datastructures
        cout << endl << "term-table: " << endl;
        cout << "{";
        cout << "}" << endl;
      }
#endif 
      return o;
    } // Maintains normalization. 

    const char* getDomainName () const {return "Anti-Unification(T)";}

  }; // class anti_unif
} // namespace ikos

#endif // IKOS_ANTI_UNIF_HPP
