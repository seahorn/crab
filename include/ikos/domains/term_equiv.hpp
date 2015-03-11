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
#include <ikos/common/types.hpp>
#include <ikos/common/bignums.hpp>
#include <ikos/algorithms/linear_constraints.hpp>
#include <ikos/domains/numerical_domains_api.hpp>
#include <ikos/domains/bitwise_operators_api.hpp>
#include <ikos/domains/division_operators_api.hpp>
#include <ikos/domains/intervals.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/optional.hpp>
#include <ikos/domains/term_expr.hpp>

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

    typedef typename dom_t::variable_t dom_var_t;
    typedef typename Info::Alloc dom_var_alloc_t;

    typedef bound<Number> bound_t;
    // typedef interval_domain< Number, VariableName > intervals_t;
    
   public:
    typedef variable< Number, VariableName >                 variable_t;
    typedef linear_constraint< Number, VariableName >        linear_constraint_t;
    typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
    typedef linear_expression< Number, VariableName >        linear_expression_t;
    typedef anti_unif<Info>        anti_unif_t;
    typedef interval< Number >                               interval_t;
    typedef interval_domain< Number, VariableName >          interval_domain_t;

    typedef term::term_table< Number, operation_t > ttbl_t;
    typedef typename ttbl_t::term_id_t term_id_t;
    typedef patricia_tree_set< VariableName >  varname_set_t;
     
   private:
    typedef typename dom_t::linear_constraint_t        dom_lincst_t;
    typedef typename dom_t::linear_expression_t        dom_linexp_t;

    typedef container::flat_map< variable_t, term_id_t > var_map_t;
    typedef container::flat_map< term_id_t, dom_var_t > term_map_t;
    // typedef std::map< term_id_t, dom_var_t > term_map_t;
//    typedef typename map_t::value_type value_type;
    typedef unsigned char BOOL;

   private:
    bool     _is_bottom;
    // Uses a single state of the underlying domain.
    ttbl_t          _ttbl;
    dom_t           _impl;
    dom_var_alloc_t _alloc;
    var_map_t       _var_map;
    term_map_t      _term_map;
    bool            _is_normalized;

    void set_to_bottom (){
      this->_is_bottom = true;
    }

   private:
    anti_unif(bool is_top): _is_bottom(!is_top), _is_normalized(true) { }

    // x = y op [lb,ub]
    term_id_t term_of_itv(bound_t lb, bound_t ub)
    {
      // FIXME: Add special handling for [x, x].
      term_id_t t_itv = _ttbl.fresh_var();
      dom_var_t dom_itv = domvar_of_term(t_itv);
      _impl.apply_constraint(dom_itv, true , ub);    // adding  x <=  c
      _impl.apply_constraint(dom_itv, false, -lb); // adding -x <= -c
      return t_itv;
    }

    term_id_t term_of_expr(operation_t op, term_id_t ty, term_id_t tz)
    {
      optional<term_id_t> opt_tx = _ttbl.find_ftor(op, ty, tz);
      if(opt_tx)
      {
        // If the term already exists, we can learn nothing.
        return *opt_tx;
      } else {
        // Otherwise, assign the term, and evaluate.
        term_id_t tx = _ttbl.apply_ftor(op, ty, tz);
        _impl.assign(op,
            domvar_of_term(tx),
            domvar_of_term(ty), domvar_of_term(tz));
        return tx;
      }
    }

    void apply(operation_t op, VariableName x, VariableName y, bound_t lb, bound_t ub){	
      term_id_t t_x = term_of_expr(op, term_of_var(y), term_of_itv(lb, ub));
      _var_map.insert(std::make_pair(x, t_x));
    }
    
    // check satisfiability of cst using intervals
    // Only to be used if cst is too hard for octagons
    bool check_sat(linear_constraint_t cst)  {
      dom_lincst_t dom_cst(rename_linear_cst(cst));
      return _impl->check_sat(dom_cst);
      return (!_is_bottom);
    }

    interval_t to_interval(VariableName x, bool requires_normalization) { 
      // projection requires normalization.
      optional<dom_var_t> dom_var(domvar_of_var(x));
      if(dom_var)
      {
        return _impl->check_sat(dom_var);
      } else {
        return interval_t::top();
      }
    } //Maintains normalization.
    
    void resize(){
      // _dbm.resize(_map.size());
    }
    
    void is_normalized(bool b){
      _is_normalized= b;
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
       _ttbl(o._ttbl), _impl(o._impl),
       _var_map(o._var_map), _term_map(o._term_map), _is_normalized(o._is_normalized)
    { } 

    anti_unif_t operator=(anti_unif_t o) {
      _is_bottom= o.is_bottom();
      _var_map= o._var_map;
      _term_map= o._term_map;
      _is_normalized= o._is_normalized;
      return *this;
    }
    
    bool is_bottom() {
      return _is_bottom;
    }
    
    bool is_top() {
      return !_var_map.size() && !is_bottom();
    }
    
    bool is_normalized(){
      return _is_normalized;
    }

    varname_set_t get_variables() const {
      varname_set_t vars;
      for(auto& p : _var_map)
      {
        variable_t v(p.first);
        vars += v.name();
      }
      return vars;
    }

    // Compute the strong closure algorithm
    void normalize(){
      if(_is_normalized){
        return;
      }

      fprintf(stderr, "WARNING: ANTI_UNIF::normalize not yet implemented.\n");
      // _impl->normalize();
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
      }
    }  // Maintains normalization.

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
        // FIXME
        return top();
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
        // FIXME
        return top();
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
      }
    }	// Returned matrix is not normalized.

    void operator-=(VariableName v) {
      // BAIL("Anti_UNIF::operator -= not yet implemented.");
      // Remove a variable from the scope
      auto it(_var_map.find(v));
      if(it != _var_map.end())
      {
        term_id_t t = (*it).second;
        _var_map.erase(it); 

        std::vector<term_id_t> forgotten;
        _ttbl.deref(t, forgotten);

        for(term_id_t ft : forgotten)
        {
          dom_var_t dom_v(domvar_of_term(ft));
          _impl -= dom_v.name();
        }
      }
    }

    // Build the tree for a linexpr, and ensure that
    // values for the subterms are sustained.
    term_id_t build_linexpr(dom_linexp_t& e)
    {
      // FIXME - implement.
      term_id_t t = _ttbl.fresh_var();
      dom_var_t v(domvar_of_term(t));
      _impl.assign(v.name(), e);
      return t; 
    }

    void assign(VariableName x_name, linear_expression_t e) {
      if (this->is_bottom())
      { 
        return;
      } else {
        variable_t x(x_name);

        dom_linexp_t dom_e(rename_linear_expr(e));
        term_id_t tx(build_linexpr(e));
        _var_map.insert(std::make_pair(x, tx));

        return;
      }
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
    }	

    term_id_t term_of_var(variable_t v)
    {
      auto it(_var_map.find(v)); 
      if(it != _var_map.end())
      {
        return (*it).second;
      } else {
        // Allocate a fresh term
        term_id_t id(_ttbl.fresh_var());
        _var_map[v] = id;
        return id;
      }
    }

    dom_var_t domvar_of_term(term_id_t id)
    {
      typename term_map_t::iterator it(_term_map.find(id));
      if(it != _term_map.end())
      {
        return (*it).second;
      } else {
        // Allocate a fresh variable
        dom_var_t dvar(_alloc.next());
        _term_map.insert(std::make_pair(id, dvar));
        return dvar;
      }
    }

    dom_var_t domvar_of_var(variable_t v)
    {
      return domvar_of_term(term_of_var(v));
    }

    // Remap a linear constraint to the domain.
    dom_linexp_t rename_linear_expr(linear_expression_t exp)
    {
      dom_linexp_t dom_exp;
      for(auto v : exp.variables())
      {
        dom_exp = dom_exp + exp[v]*domvar_of_var(v);
      }
      return dom_exp;
    }

    dom_lincst_t rename_linear_cst(linear_constraint_t cst)
    {
      return dom_lincst_t(rename_linear_expr(cst.expression()), cst.kind());
    }

    void operator+=(linear_constraint_t cst) {  
      dom_lincst_t cst_rn(rename_linear_cst(cst));
      _impl += cst_rn;
      is_normalized(false); 
      return;
    } // Sets state to not normalized.
    
    // Add a system of linear constraints
    void operator+=(linear_constraint_system_t cst) {  // Does not require normalization.
      for(typename linear_constraint_system_t::iterator it=cst.begin(); it!= cst.end(); ++it){
        this->operator+=(*it);
      }
    } // Sets state to not normalized.
    
    // abstract the variable
    // GKG: Looks like this returns the index of variable v.
    /*
    boost::optional<typename map_t::iterator> abstract(variable_t v) {	
      // Requires normalization.
      BAIL("ANTI-UNIF: abstract not yet fully implemented.");
      typename map_t::iterator it(_map.find(v));
      if(it!= _map.end()) {
        normalize();
        return boost::optional<typename map_t::iterator>(it);
      }
      return boost::optional<typename map_t::iterator>();
    }  //Maintains normalization.
    */
    
    interval_t operator[](VariableName x) { 
      return to_interval(x, true);
    } 

    void set(VariableName x, interval_t intv){
      typename var_map_t::iterator it = this->_var_map.find(x);
      dom_var_t dx(domvar_of_var(x));
      _impl->set(x, intv);
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
        for(typename var_map_t::iterator it=_var_map.begin(); it!= _var_map.end(); ++it){
          variable_t x(it->first.name());
          dom_var_t dx(domvar_of_var(x));
          itv.set(x, _impl->to_interval(dx));
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
      
      /*
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
      */
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
      if(_var_map.size()== 0) {
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
