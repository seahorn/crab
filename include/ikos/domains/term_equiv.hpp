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

#define VERBOSE 

namespace ikos {
  // Function to call D->normalize only if it exists.
  // FIXME: Figure out how to abuse SFINAE to do this.
  /*
  template<class D>
  void _normalize(D& elt)
  { elt->normalize(); }
  */

  template<class D>
  void _normalize(D& elt)
  { }

  // If is_normalized exists, call it. Otherwise, return true.
  /*
  template<class D>
  bool _is_normalized(D& elt)
  { return elt->is_normalized(); }
  */

  template<class D>
  bool _is_normalized(D& elt)
  { return true; }

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
    typedef patricia_tree_set< dom_var_t > domvar_set_t;

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
    // typedef typename dom_t::Number                     dom_number;
    // WARNING: assumes the underlying domain uses the same number type.
    typedef typename Info::Number                     dom_number;
    typedef typename dom_t::linear_constraint_t        dom_lincst_t;
    typedef typename dom_t::linear_expression_t        dom_linexp_t;

    typedef typename linear_expression_t::component_t linterm_t;

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

    void set_to_bottom (){
      this->_is_bottom = true;
    }

   private:
    anti_unif(bool is_top): _is_bottom(!is_top) { }

    anti_unif(dom_var_alloc_t alloc, var_map_t vm, ttbl_t tbl, term_map_t tmap, dom_t impl)
      : _is_bottom(false), _ttbl(tbl), _impl(impl), _alloc(alloc), _var_map(vm), _term_map(tmap)
    { }

    // x = y op [lb,ub]
    term_id_t term_of_itv(bound_t lb, bound_t ub)
    {
//      if(lb == ub)
//        return build_const(lb);

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
      check_terms();
    }
    
    /*
    bool check_sat(linear_constraint_t cst)  {
      dom_lincst_t dom_cst(rename_linear_cst(cst));
      return _impl.check_sat(dom_cst);
      return (!_is_bottom);
    }
    */

    /*
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
    */
    
  public:
    static anti_unif_t top() {
      return anti_unif(true);
    }
    
    static anti_unif_t bottom() {
      return anti_unif(false);
    }
    
  public:
    // Constructs top octagon, represented by a size of 0.
    anti_unif(): _is_bottom(false) { }
    
    anti_unif(const anti_unif_t& o): 
       writeable(), 
       numerical_domain<Number, VariableName >(),
       bitwise_operators< Number, VariableName >(),
       division_operators< Number, VariableName >(),   
       _is_bottom(o._is_bottom), 
       _ttbl(o._ttbl), _impl(o._impl),
       _var_map(o._var_map), _term_map(o._term_map)
    { check_terms(); } 

    anti_unif_t operator=(anti_unif_t o) {
      _is_bottom= o.is_bottom();
      _var_map= o._var_map;
      _term_map= o._term_map;
      _ttbl = o._ttbl;
//      check_terms();
      return *this;
    }
    
    bool is_bottom() {
      return _is_bottom;
    }
    
    bool is_top() {
      return !_var_map.size() && !is_bottom();
    }
    
    bool is_normalized(){
      return _is_normalized(_impl);
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
      _normalize(_impl);
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
        typename ttbl_t::term_map_t gen_map;

        // Build up the mapping of o onto this, variable by variable.
        // Assumption: the set of variables in x & o are common.
        for(auto p : _var_map)
        {
          if(!_ttbl.map_leq(o._ttbl, term_of_var(p.first), o.term_of_var(p.first), gen_map))
            return false;
        }
        // We now have a mapping of reachable y-terms to x-terms.
        // Create copies of _impl and o._impl with a common
        // variable set.
        dom_t x_impl(_impl);
        dom_t y_impl(o._impl);

        // Perform the mapping
        domvar_set_t xvars;
        domvar_set_t yvars;
        for(auto p : gen_map)
        {
          dom_var_t vt = _alloc.next();
          dom_var_t vx = domvar_of_term(p.second); 
          dom_var_t vy = o.domvar_of_term(p.first);

          xvars += vx;
          yvars += vy;

          x_impl.assign(vt.name(), dom_linexp_t(vx));
          y_impl.assign(vt.name(), dom_linexp_t(vy));
        }
        for(auto vx : xvars)
          x_impl -= vx.name();
        for(auto vy : yvars)
          y_impl -= vy.name();

        return x_impl <= y_impl;
      }
    } 

    anti_unif_t operator|(anti_unif_t o) {
      // Requires normalization of both operands
//      std::cout << "SIZES: " << _var_map.size() << ", " <<
//        o._var_map.size() << std::endl;
      normalize();
      o.normalize();
      if (is_bottom()) {
        return o;
      } 
      else if(o.is_bottom()) {
        return *this;
      } 
      else {
        // First, we need to compute the new term table.
        ttbl_t out_tbl;
        // Mapping of (term, term) pairs to terms in the join state
        typename ttbl_t::gener_map_t gener_map;
        
        var_map_t out_vmap;

        // For each program variable in state, compute a generalization
        for(auto p : _var_map)
        {
          variable_t v(p.first);
          term_id_t tx(term_of_var(v));
          term_id_t ty(o.term_of_var(v));

          term_id_t tz = _ttbl.generalize(o._ttbl, tx, ty, out_tbl, gener_map);
          out_vmap[v] = tz;
        }

        // Rename the common terms together
        dom_t x_impl(_impl);
        dom_t y_impl(o._impl);

        // Perform the mapping
        term_map_t out_map;
        domvar_set_t xvars;
        domvar_set_t yvars;
        for(auto p : gener_map)
        {
          auto txy = p.first;
          term_id_t tz = p.second;
          dom_var_t vt = _alloc.next();
          out_map.insert(std::make_pair(tz, vt));

          dom_var_t vx = domvar_of_term(txy.first);
          dom_var_t vy = o.domvar_of_term(txy.second);

          xvars += vx;
          yvars += vy;

          x_impl.assign(vt.name(), dom_linexp_t(vx));
          y_impl.assign(vt.name(), dom_linexp_t(vy));
        }
//        std::cout << "ren_0(x) = " << x_impl << std::endl;
//        std::cout << "ren_0(y) = " << y_impl << std::endl;
        for(auto vx : xvars)
          x_impl -= vx.name();
        for(auto vy : yvars)
          y_impl -= vy.name();

//        std::cout << "ren(x) = " << x_impl << std::endl;
//        std::cout << "ren(y) = " << y_impl << std::endl;
        dom_t x_join_y = x_impl|y_impl;
        return anti_unif(_alloc, out_vmap, out_tbl, out_map, x_join_y);
      }
    }

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
        // First, we need to compute the new term table.
        ttbl_t out_tbl;
        // Mapping of (term, term) pairs to terms in the join state
        typename ttbl_t::gener_map_t gener_map;
        
        var_map_t out_vmap;
        // For each program variable in state, compute a generalization
        for(auto p : _var_map)
        {
          variable_t v(p.first);
          term_id_t tx(term_of_var(v));
          term_id_t ty(o.term_of_var(v));

          term_id_t tz = _ttbl.generalize(o._ttbl, tx, ty, out_tbl, gener_map);
          out_vmap[v] = tz;
        }

        // Rename the common terms together
        dom_t x_impl(_impl);
        dom_t y_impl(o._impl);

        // Perform the mapping
        term_map_t out_map;
        domvar_set_t xvars;
        domvar_set_t yvars;
        for(auto p : gener_map)
        {
          auto txy = p.first;
          term_id_t tz = p.second;
          dom_var_t vt = _alloc.next();
          out_map.insert(std::make_pair(tz, vt));

          dom_var_t vx = domvar_of_term(txy.first);
          dom_var_t vy = o.domvar_of_term(txy.second);

          xvars += vx;
          yvars += vy;

          x_impl.assign(vt.name(), dom_linexp_t(vx));
          y_impl.assign(vt.name(), dom_linexp_t(vy));
        }
        for(auto vx : xvars)
          x_impl -= vx.name();
        for(auto vy : yvars)
          y_impl -= vy.name();

        dom_t x_widen_y = x_impl||y_impl;
        return anti_unif(_alloc, out_vmap, out_tbl, out_map, x_widen_y);
      }
    }

    // Meet
    anti_unif_t operator&(anti_unif_t o) {
      // Does not require normalization of any of the two operands
      if (is_bottom() || o.is_bottom()) {
        return bottom();
      } else {
        BAIL("ANTI-UNIF: meet not yet implemented.");
        return top();
      }
    }
    
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

    void check_terms(void)
    {
      for(auto p : _var_map)
      {
        assert(p.second < _ttbl.size());
      }
    }
    template<class T>
    T check_terms(T& t)
    {
      check_terms();
      return t;
    }

    void rebind_var(variable_t& x, term_id_t tx)
    {
      auto it(_var_map.find(x));
      if(it != _var_map.end())
        _var_map.erase(it);
      _var_map.insert(std::make_pair(x, tx));
    }

    // Build the tree for a linexpr, and ensure that
    // values for the subterms are sustained.
    term_id_t build_const(const Number& n)
    {
      dom_number dom_n(n);
      optional<term_id_t> opt_n(_ttbl.find_const(dom_n));
      if(opt_n)
      {
        return *opt_n;
      } else {
        term_id_t term_n(_ttbl.make_const(dom_n));
        dom_var_t v = domvar_of_term(term_n);
        dom_linexp_t exp(n);
        _impl.assign(v.name(), exp);
        return term_n;
      }
    }

    term_id_t build_linterm(linterm_t term)
    {
      return build_term(OP_MULTIPLICATION,
          build_const(term.first),
          term_of_var(term.second));
    }

    term_id_t build_linexpr(linear_expression_t& e)
    {
      typename linear_expression_t::iterator it = e.begin();
      if(it == e.end())
      {
        Number cst = e.constant();
        return build_const(cst);
      }
      term_id_t t(build_linterm((linterm_t) *it));
      it++;
      for(; it != e.end(); it++)
      {
        t = build_term(OP_ADDITION, t, build_linterm(*it));
      }
      return t; 
    }

    term_id_t build_term(operation_t op, term_id_t ty, term_id_t tz)
    {
      // Check if the term already exists
      optional<term_id_t> eopt(_ttbl.find_ftor(op, ty, tz));
      if(eopt)
      {
        return *eopt;
      } else {
        // Create the term
        term_id_t tx = _ttbl.apply_ftor(op, ty, tz);
        dom_var_t v(domvar_of_term(tx));
        // Set up the evaluation.
        _impl.apply(op, v.name(), domvar_of_term(ty).name(), domvar_of_term(tz).name());
        return tx;
      }
    }

    void assign(VariableName x_name, linear_expression_t e) {
      if (this->is_bottom())
      { 
        return;
      } else {
        variable_t x(x_name);

//        dom_linexp_t dom_e(rename_linear_expr(e));
        term_id_t tx(build_linexpr(e));

        rebind_var(x, tx);

        check_terms();
        return;
      }
    }

    // Apply operations to variables.

    // x = y op z
    void apply(operation_t op, VariableName x, VariableName y, VariableName z){	
      if (this->is_bottom())
      {
        return;   
      } else {
        variable_t vx(x);
          
        term_id_t tx(build_term(op, term_of_var(y), term_of_var(z)));
        rebind_var(vx, tx);
      }
      check_terms();
    }
    
    // x = y op k
    void apply(operation_t op, VariableName x, VariableName y, Number k){	
      if (this->is_bottom())
      {
        return;   
      } else {
        variable_t vx(x);
          
        term_id_t tx(build_term(op, term_of_var(y), build_const(k)));
        rebind_var(vx, tx);
      }

      check_terms();
      return;
    }

    term_id_t term_of_var(variable_t v)
    {
      auto it(_var_map.find(v)); 
      if(it != _var_map.end())
      {
        // assert ((*it).first == v);
        assert(_ttbl.size() > (*it).second);
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
      Number cst(exp.constant());
      dom_linexp_t dom_exp(cst);
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
    
    /*
    interval_t operator[](VariableName x) { 
      return to_interval(x, true);
    } 
    */

    void set(VariableName x, interval_t intv){
      typename var_map_t::iterator it = this->_var_map.find(x);
      dom_var_t dx(domvar_of_var(x));
      _impl->set(x, intv);
    }

    /*
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
          itv.set(x, _impl.to_interval(dx));
        }
        return itv;
      }
    }
    */

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

      bool first = true;
      o << "{" ;
      for(auto p : _var_map)
      {
        if(first)
          first = false;
        else
          o << ", ";
        o << p.first << " -> t" << p.second;
      }
      o << "}";
     
#ifdef VERBOSE
      /// For debugging purposes     
      { // print intervals
        //interval_domain_t intervals = to_intervals();
        // o << intervals;
        o << _impl;
      }
      { // print internal datastructures
        o << endl << "term-table: " << endl;
        o << "{";
        o << _ttbl;
        o << "}" << endl;
      }
#endif 
      return o;
    }

    const char* getDomainName () const {return "Anti-Unification(T)";}

  }; // class anti_unif
} // namespace ikos

#endif // IKOS_ANTI_UNIF_HPP
