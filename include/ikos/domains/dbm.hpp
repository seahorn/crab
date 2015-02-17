/*******************************************************************************
 *
 * Difference Bounds Matrix domain
 *
 * Based on the paper "Fast and Flexible Difference Constraint
 * Propagation for DPLL(T) by Cotton and Maler.
 *
 * Author: Graeme Gange (gkgange@unimelb.edu.au)
 ******************************************************************************/

#ifndef IKOS_DBM_HPP
#define IKOS_DBM_HPP

#include <ikos/common/mergeable_map.hpp>
#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp>
#include <ikos/domains/dbm/dbm.h>
#include <ikos/domains/dbm/expr.h>
#include <ikos/domains/dbm/dbm_iter.h>
#include <ikos/domains/intervals.hpp>
#include <ikos/domains/numerical_domains_api.hpp>
#include <ikos/domains/bitwise_operators_api.hpp>
#include <ikos/domains/division_operators_api.hpp>

#include <boost/optional.hpp>
#include <boost/unordered_set.hpp>
#include <boost/lexical_cast.hpp>

namespace ikos {

using namespace boost;
using namespace std;

namespace DBM_impl
{
   template<typename Number>
   int ntoi(Number n);

   template<>
   int ntoi(z_number k){
     ostringstream buf;
     buf << k;
     return boost::lexical_cast<int>(buf.str());
   }

   template<typename Number>
   Number iton(int n);

   template<>
   z_number iton(int n){
     return z_number(n);
   }
}; // end namespace DBM_impl

template< typename Number, typename VariableName, std::size_t size = 3000 >
class DBM: public writeable,
           public numerical_domain<Number, VariableName >,
           public bitwise_operators< Number, VariableName >,
           public division_operators< Number, VariableName >{
  
 private:

  struct Index: public writeable{
    Index(index_t idx): writeable(), _idx(idx){ }
    Index(int idx): writeable(), _idx(idx){ }
    Index(const Index &o): writeable(), _idx(o._idx) { }
    Index& operator=(Index o){
      this->_idx = o._idx;
      return *this;
    }
    index_t index(){
      return _idx;
      }
    bool operator<(const Index&o) const {
      return (_idx < o._idx);
    }    
    ostream& write(ostream& o) {
      o << _idx;
      return o;
      }
    index_t _idx;
  };
  
 public:
  typedef variable< Number, VariableName >                 variable_t;
  typedef linear_constraint< Number, VariableName >        linear_constraint_t;
  typedef typename linear_constraint_t::kind_t             constraint_kind_t;
  typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
  typedef linear_expression< Number, VariableName >        linear_expression_t;
  typedef DBM< Number, VariableName >                      DBM_t;
  typedef patricia_tree_set< VariableName >                varname_set_t;


 public:
  typedef interval< Number >  interval_t;

 private:
  typedef bound<Number >  bound_t;
  typedef interval_domain< Number, VariableName > intervals_t;

  typedef mergeable_map< Index, VariableName > map_t;    


 private:
  unsigned _sz;
  dbm      _dbm;
  map_t    _map;
  varname_set_t _all_vars;

 private:

  // return true if edge (i,j) exists in the DBM
  bool has_edge(int i, int j)
  {
    return (src_is_live(this->_dbm, i) && 
            in_graph(this->_dbm, i, j));
  }
                         

  // special variables
  inline int ZERO(){ return (int) this->_sz-1; }
  inline int TMP() { return (int) this->_sz-2; }
    
  void add_variable(VariableName x){
    this->_all_vars += x;
    int k = x.index();
    if (!this->_map[k]){
      if ((unsigned) k >= (this->_sz - 1))
        throw ikos::error("DBM: need to enlarge the matrix");
      else
        _map.set(k,x);
    }
  }

 private:

  DBM(bool is_bottom): 
    writeable(), _sz(size), _dbm(dbm_bottom()), _map(), _all_vars() { 
    if (!is_bottom)
      _dbm = dbm_top(_sz);
  }

  DBM(dbm dbm, map_t map, varname_set_t all_vars): 
      writeable(), _sz(dbm->sz), _dbm(dbm), _map(map), _all_vars(all_vars) { }

  void swap(dbm& x, dbm& y){
    dbm tmp = x;
    x = y;
    y = tmp;
  }

  void _assign(VariableName x, exp_t exp){
    dbm ret = NULL;      
    ret = dbm_assign(x.index(), exp, this->_dbm);
    dbm_dealloc(this->_dbm);
    swap(this->_dbm, ret);
    exp_dealloc(exp);
  }

  void _assign_tmp(VariableName x){

    // pre: _forget(TMP())
    dexpr e1;     
    e1.args[0] = x.index();
    e1.args[1] = TMP();
    e1.konst   = 0;
    e1.kind    = D_DIFF;
    _apply_dexpr(e1);
    dexpr e2;     
    e2.args[0] = TMP();
    e2.args[1] = x.index();
    e2.konst   = 0;
    e2.kind    = D_DIFF;
    _apply_dexpr(e2);
  }

  void _apply_cond(ucon con){
    dbm ret = NULL;
    ret = dbm_cond(con, this->_dbm);
    dbm_dealloc(this->_dbm);
    swap(this->_dbm, ret);
  }

  void _apply_dexpr(dexpr d){
    dbm ret = NULL;
    ret = dbm_apply_dexpr(d, this->_dbm);
    dbm_dealloc(this->_dbm);
    swap(this->_dbm, ret);
  }

  void _forget(index_t idx){
    dbm ret = NULL;
    ret = dbm_forget(idx, this->_dbm);
    dbm_dealloc(this->_dbm);
    swap(this->_dbm, ret);
  }
    
  void _forget(vector<int> idxs){
    dbm ret = NULL;
    ret = dbm_forget_array(&idxs[0], idxs.size(), this->_dbm);
    dbm_dealloc(this->_dbm);
    swap(this->_dbm, ret);
  }

  void set_to_bottom(){
    dbm_dealloc(this->_dbm);
    this->_dbm = dbm_bottom();
  }

  VariableName map_index(index_t k){
    boost::optional<VariableName> v = this->_map[k];
    if (!v)
      throw ikos::error("DBM: index cannot be mapped to a variable name") ;
    else 
      return *v;
  }

  void apply(operation_t op, VariableName x, index_t y, index_t z){	
    switch(op){
      case OP_ADDITION:
        {
          exp_t exp = exp_add(exp_var(y), exp_var(z));
          _assign(x, exp);
        }
        break;
      case OP_SUBTRACTION:
        {
          exp_t exp = exp_sub(exp_var(y), exp_var(z));
          _assign(x, exp);
        }
        break;
      case OP_MULTIPLICATION:
        {
          exp_t exp = exp_mul(exp_var(y), exp_var(z));
          _assign(x, exp);
        }
        break;
      case OP_DIVISION:
        {
          exp_t exp = exp_div(exp_var(y), exp_var(z));
          _assign(x, exp);
        }
        break;
      default:;;
    }
  }	

  void apply(operation_t op, VariableName x, index_t y, Number k){	
    int z = DBM_impl::ntoi<Number>(k);
    switch(op){
      case OP_ADDITION:
        {
          exp_t exp = exp_add(exp_var(y), exp_const(z));
          _assign(x, exp);
        }
        break;
      case OP_SUBTRACTION:
        {
          exp_t exp = exp_sub(exp_var(y), exp_const(z));
          _assign(x, exp);
        }
        break;
      case OP_MULTIPLICATION:
        {
          exp_t exp = exp_mul(exp_var(y), exp_const(z));
          _assign(x, exp);
        }
        break;
      case OP_DIVISION:
        {
          exp_t exp = exp_div(exp_var(y), exp_const(z));
          _assign(x, exp);          
        }
        break;
      default:;;
    }
  }	

  // check satisfiability of cst using intervals
  // Only to be used if cst is too hard for dbm
  bool intervals_check_sat(linear_constraint_t cst)  {
    if (is_top())    return true;
    if (is_bottom()) return false;

    typename linear_constraint_t::variable_set_t vars = cst.variables();
    intervals_t inv;
    for(typename linear_constraint_t::variable_set_t::iterator it = vars.begin(); it!= vars.end(); ++it)
      inv.set (it->name(), this->operator[](it->name()));
    inv += cst;
    return !inv.is_bottom();
  }

  // x is unit coefficient
  void add_constraint(int coef_x, variable_t x, Number c, linear_constraint_t cst){

    if (!(coef_x == 1 || coef_x == -1))
      throw ikos::error("DBM: coefficients can be only 1 or -1");

    add_variable(x.name());
    if (coef_x == 1){
      uterm tx = uvar(x.name().index());
      uterm ty = uconst(DBM_impl::ntoi<Number>(c));
      ucon con;
      if (cst.is_inequality()){       // x <= c
        con = mk_ucon(tx, U_LEQ, ty);
        _apply_cond(con);
      }
      else if (cst.is_equality()){    // x == c
        con = mk_ucon(tx, U_EQ, ty);
        _apply_cond(con);
      }
      else if (cst.is_disequation()){ // x != c
        // DBM ignores disequations so we better use intervals
        if (!intervals_check_sat(cst)){ this->set_to_bottom(); }
      }
    }
    else{ 
      if (cst.is_inequality()){       // 0-x <= c
        dexpr e;     
        e.args[0] = ZERO();
        e.args[1] = x.name().index();
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        _apply_dexpr(e);
      }
      else if (cst.is_equality()){    // -x == c iff x == -c
        uterm tx = uvar(x.name().index());
        uterm tc = uconst(-DBM_impl::ntoi<Number>(c));
        ucon con = mk_ucon(tx, U_EQ, tc);
        _apply_cond(con);
      }
      else if (cst.is_disequation()){ // -x != c iff x != -c
        // DBM ignores disequations so we better use intervals
        if (!intervals_check_sat(cst)){ this->set_to_bottom(); }
      }
    }
  }

  // x, y are unit coefficient
  void add_constraint(int coef_x, variable_t x, int coef_y, variable_t y, Number c, linear_constraint_t cst){

    if (!(coef_x == 1 || coef_x == -1))
      throw ikos::error("DBM: coefficients can be only 1 or -1");
    if (!(coef_y == 1 || coef_y == -1))
      throw ikos::error("DBM: coefficients can be only 1 or -1");
    if (coef_x == coef_y)
      throw ikos::error("DBM: same coefficients");        

    add_variable(x.name());
    add_variable(y.name());

    if (coef_x == 1 && coef_y == -1){  
      if (cst.is_inequality()){       // x - y <= c
        dexpr e;     
        e.args[0] = x.name().index();
        e.args[1] = y.name().index();
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        _apply_dexpr(e);
      }
      else if (cst.is_equality()){    // x - y == c iff x - y <= c and y -x <= -c
        dexpr e1;     
        e1.args[0] = x.name().index();
        e1.args[1] = y.name().index();
        e1.konst   = DBM_impl::ntoi<Number>(c);
        e1.kind    = D_DIFF;
        _apply_dexpr(e1);
        dexpr e2;     
        e2.args[0] = y.name().index();
        e2.args[1] = x.name().index();
        e2.konst   = -DBM_impl::ntoi<Number>(c);
        e2.kind    = D_DIFF;
        _apply_dexpr(e2);
      }
      else if (cst.is_disequation())
      {  // x - y != c
        if (has_edge(x.name ().index (), y.name ().index ()) && 
            has_edge(y.name ().index (), x.name ().index ()))
        {
          int k1 = edge_val(this->_dbm, x.name ().index (), y.name ().index ());
          int k2 = edge_val(this->_dbm, y.name ().index (), x.name ().index ());
          if ( (k1 + k2 == 0) && (k1 == DBM_impl::ntoi<Number>(c)))
          {
            set_to_bottom ();
          }
        }
        else if (!intervals_check_sat(cst))
        { 
          this->set_to_bottom(); 
        }
      }
    }      
    else{ // coef_x == -1 && coef_y == 1
      if (cst.is_inequality()){       // y - x <= c
        dexpr e;     
        e.args[0] = y.name().index();
        e.args[1] = x.name().index();
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        _apply_dexpr(e);
      }
      else if (cst.is_equality()){    // y - x == c iff y - x <= c and x -y <= -c
        dexpr e1;     
        e1.args[0] = y.name().index();
        e1.args[1] = x.name().index();
        e1.konst   = DBM_impl::ntoi<Number>(c);
        e1.kind    = D_DIFF;
        _apply_dexpr(e1);
        dexpr e2;     
        e2.args[0] = x.name().index();
        e2.args[1] = y.name().index();
        e2.konst   = -DBM_impl::ntoi<Number>(c);
        e2.kind    = D_DIFF;
        _apply_dexpr(e2);
      }
      else if (cst.is_disequation())
      { // y - x != c
        if (has_edge(y.name ().index (), x.name ().index ()) && 
            has_edge(x.name ().index (), y.name ().index ()))
        {
          int k1 = edge_val(this->_dbm, y.name ().index (), x.name ().index ());
          int k2 = edge_val(this->_dbm, x.name ().index (), y.name ().index ());
          if ( (k1 + k2 == 0) && (k1 == DBM_impl::ntoi<Number>(c)))
          {
            set_to_bottom ();
          }
        }
        else if (!intervals_check_sat(cst))
        { 
          this->set_to_bottom(); 
        }
      }
    }
  }

  linear_expression_t index_to_expr(index_t k)
  {
    if (k == static_cast<index_t>(ZERO()))
      return linear_expression_t(Number (0));
    else
      return linear_expression_t(map_index (k));
  }
  

 public:

  static DBM_t top() { return DBM(false); }
    
  static DBM_t bottom() { return DBM(true); }

  ~DBM(){
    dbm_dealloc(this->_dbm);
  }
    
 public:

  DBM(): writeable(), _sz(size), _dbm(dbm_top(size)), _map(), _all_vars() {} 
           

  DBM(const DBM_t& o): 
    writeable(), 
    numerical_domain<Number, VariableName >(),
    bitwise_operators< Number, VariableName >(),
    division_operators< Number, VariableName >(),
    _sz(o._sz), _dbm(dbm_copy(o._dbm)), _map(o._map), _all_vars(o._all_vars) { }
   
  DBM_t operator=(DBM_t o) {
    this->_sz  = o._sz;
    dbm_dealloc(this->_dbm);
    this->_dbm = dbm_copy(o._dbm); 
    this->_map = o._map;
    this->_all_vars = o._all_vars;
    return *this;
  }

  varname_set_t get_variables() const{
    return this->_all_vars;
  }

  bool is_bottom() {
    return dbm_is_bottom(_dbm);
  }
    
  bool is_top() {
    return dbm_is_top(_dbm);
  }
    
  bool operator<=(DBM_t o)  {       
    if (is_bottom()) 
      return true;
    else if(o.is_bottom())
      return false;
    else 
      return dbm_is_leq(this->_dbm, o._dbm);
  }  

  bool operator==(DBM_t o){
    return (*this <= o && o <= *this);
  }

  DBM_t operator|(DBM_t o) {
    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else{
      return DBM_t(dbm_join(this->_dbm, o._dbm), 
                   this->_map | o._map, 
                   this->_all_vars | o._all_vars);
    }
  } 

  DBM_t operator||(DBM_t o) {	
    return DBM_t(dbm_widen(this->_dbm, o._dbm), 
                 this->_map | o._map,
                 this->_all_vars | o._all_vars);
  } 

  DBM_t operator&(DBM_t o) { 
    if (is_top())
      return o;
    else if (o.is_top())
      return *this;
    else{
      return DBM_t(dbm_meet(this->_dbm, o._dbm), 
                   this->_map | o._map,
                   this->_all_vars | o._all_vars);
                     
    }
  }	
    
  DBM_t operator&&(DBM_t o) {	
    return this->operator&(o);
  }	

  void normalize() {
    dbm_canonical(this->_dbm);
  }

  void operator-=(VariableName v) {
    _forget(v.index());
  }

  void operator-=(vector<VariableName> vs) {
    vector<int> indexes;
    for(unsigned int i=0;i<vs.size();i++)
      indexes.push_back((int) vs[i].index());
    _forget(indexes);
  }

  void assign(VariableName x, linear_expression_t e) {

    if(this->is_bottom())
      return;

    add_variable(x);
    if (e.is_constant()){
      exp_t exp = exp_const(DBM_impl::ntoi<Number>(e.constant()));
      _assign(x,exp);
    }
    else if  (boost::optional<variable_t> v = e.get_variable()){
      VariableName y = (*v).name();
      if (!(x==y)){
        add_variable(y);
        exp_t exp = exp_var(y.index());        
        _assign(x,exp);
      }
    }
    else
      throw ikos::error("DBM: only supports constant or variable on the rhs of assignment");
  }

  void apply(operation_t op, VariableName x, VariableName y, VariableName z){	

    add_variable(x);
    add_variable(y);
    add_variable(z);

    if (x == y){
      // to make sure that lhs does not appear on the rhs
      _assign_tmp(y); 
      apply(op, x, TMP(), z.index());
      _forget(TMP());
    }
    else if (x == z){
      // to make sure that lhs does not appear on the rhs
      _assign_tmp(z); 
      apply(op, x, y.index(), TMP());
      _forget(TMP());
    }
    else{
      if (x == y && y == z)
        throw ikos::error("DBM: does not support x := x + x ");
      else
        apply(op, x, y.index(), z.index());
    }
  }
    
  void apply(operation_t op, VariableName x, VariableName y, Number k){	
    add_variable(x);
    add_variable(y);

    if (x == y){
      // to make sure that lhs does not appear on the rhs
      _assign_tmp(y);
      apply(op, x, TMP(), k);
      _forget(TMP());
    }
    else{
      apply(op, x, y.index(), k);
    }
  }
   
   
  void operator+=(linear_constraint_t cst) {  
    if(this->is_bottom())
      return;

    if (cst.is_tautology())
      return;
      
    if (cst.is_contradiction()){
      set_to_bottom();
      return ;
    }

    linear_expression_t exp = cst.expression();
    if (exp.size() > 2){
      throw ikos::error("DBM supports constraints with at most two variables");
    }

    if (exp.size() == 0){
      cout << cst << endl;
      throw ikos::error("DBM: bad-formed constraint");
    }
      
    Number k = -exp.constant(); 

    typename linear_expression_t::iterator it=exp.begin();

    Number coef_x = it->first;
    variable_t x  = it->second;      
    ++it;
    if (it == exp.end()){
      add_constraint(DBM_impl::ntoi<Number>(coef_x), x, k, cst);
    }
    else {
      Number coef_y = it->first;
      variable_t y  = it->second;      
      add_constraint(DBM_impl::ntoi<Number>(coef_x), 
                     x, 
                     DBM_impl::ntoi<Number>(coef_y), y, k, cst);
    }
  } 
    
  void operator+=(linear_constraint_system_t cst) {  
    for(typename linear_constraint_system_t::iterator it=cst.begin(); it!= cst.end(); ++it){
      this->operator+=(*it);
    }
  }

  interval_t operator[](VariableName x) { 

    if (is_top())    return interval_t::top();
    if (is_bottom()) return interval_t::bottom();

    int k = x.index();
    if (!this->_map[k])
      return interval_t::top();
    else{
      bound_t lb("-oo");
      if (has_edge(ZERO(), k))
        lb = bound_t(DBM_impl::iton<Number>(-edge_val(this->_dbm, ZERO(), k)));

      bound_t ub("+oo");
      if (has_edge(k, ZERO()))
        ub = bound_t(DBM_impl::iton<Number>(edge_val(this->_dbm, k, ZERO())));

      interval_t i(lb,ub);
      return i;
    }
  }

  void set(VariableName x, interval_t intv){
    int k = x.index();
    if (this->_map[k])
      this->operator-=(x);
    else
      add_variable(x);
      
    if (!intv.is_top()){
      boost::optional<Number> lb = intv.lb().number();
      boost::optional<Number> ub = intv.ub().number();
      if (ub){
        // x - 0 <= ub
        dexpr e;     
        e.args[0] = k;
        e.args[1] = ZERO();
        e.konst   = DBM_impl::ntoi<Number>(*ub);
        e.kind    = D_DIFF;
        _apply_dexpr(e);
      }
      if (lb){
        // x >= lb iff 0-x <= -lb 
        dexpr e;     
        e.args[0] = ZERO();
        e.args[1] = k;
        e.konst   = -DBM_impl::ntoi<Number>(*lb);
        e.kind    = D_DIFF;
        _apply_dexpr(e);
      }
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
        throw ikos::error("DBM: unreachable");
    }
    this->set(x, xi);
  }
    
  void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k){
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
        throw ikos::error("DBM: unreachable");
    }
    this->set(x, xi);
  }
    
  // division_operators_api
    
  void apply(div_operation_t op, VariableName x, VariableName y, VariableName z){
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
          throw ikos::error("DBM: unreachable");
      }
      this->set(x, xi);
    }
  }

  void apply(div_operation_t op, VariableName x, VariableName y, Number k){
    if (op == OP_SDIV){
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
          throw ikos::error("DBM: unreachable");
      }
      this->set(x, xi);
    }
  }
    
  // Output function
  ostream& write(ostream& o) { 
    dbm_canonical(this->_dbm);

#if 0
    // print internal representation
    cout << endl;
    cout << "map: " <<  this->_map << endl;
    cout << "csts: ";
    dbm_print_to(cout, this->_dbm);
    cout << endl;
#endif 
    if(is_bottom()){
      o << "_|_";
      return o;
    }
    else if (is_top()){
      o << "{}";
      return o;
    }
    else
    {
      boost::optional<linear_constraint_system_t> inv = to_linear_constraint_system ();
      assert(inv);
      o << *inv;
      return o;
    }
  }

  const char* getDomainName () const {return "DBM";}

  boost::optional<linear_constraint_system_t> 
  to_linear_constraint_system ()
  {

    dbm_canonical(this->_dbm);
    
    if(is_bottom ())
    {
      return boost::optional<linear_constraint_system_t>();
    }

    linear_constraint_system_t csts;
    boost::unordered_set< pair<int,int> > visited; // to only consider half matrix
    for(edge_iter iter=edge_iterator(this->_dbm); !srcs_end(iter); next_src(iter))
    {
      int ii = src(iter);
      if (ii == TMP()) continue;
        
      for(; !dests_end(iter); next_dest(iter))
      {
        int jj = dest(iter);

        if (jj == TMP()) continue;
        if (visited.find(make_pair(jj,ii)) != visited.end()) continue;
          
        int k1 = edge_val(this->_dbm, ii, jj);
        if (in_graph(this->_dbm, jj, ii))
        {
          int k2 = edge_val(this->_dbm, jj, ii);
          if ((k1 + k2) == 0)
          {
            // EQUALITY 
            if (k1 == 0)
            { // ii = jj
              csts += linear_constraint_t(index_to_expr(ii) == index_to_expr(jj));
            }
            else
            {
              if (ii == ZERO())
              { // jj = -k1
                linear_expression_t e(index_to_expr(jj) + Number(k1));
                csts += linear_constraint_t(e == 0);
              }
              else if (jj == ZERO())
              { // ii = k1
                linear_expression_t e(index_to_expr(ii) - Number(k1));
                  csts += linear_constraint_t(e == 0);
                }
                else
                {
                  // ii - jj = k1
                  linear_expression_t e(index_to_expr(ii) - index_to_expr(jj));
                  csts += linear_constraint_t(e == Number(k1));
                }
              }
            visited.insert(make_pair(ii,jj));
          }
          else
          {
            // INEQUALITY
            if (ii == ZERO())
            { // jj >= -k1
              linear_expression_t e(index_to_expr(jj) + Number(k1));
              csts += linear_constraint_t(e >= 0);
            }
            else if (jj == ZERO())
            {
              // ii <= k1
              linear_expression_t e(index_to_expr(ii) - Number(k1));
              csts += linear_constraint_t(e <= 0);
            }
            else
            {
              // ii - jj <= k1
              linear_expression_t e(index_to_expr(ii) - index_to_expr(jj));
              csts += linear_constraint_t(e <= Number(k1));
            }
          }
        }
        else
        {
          if (ii == ZERO())
          {// jj >= -k1
            linear_expression_t e(index_to_expr(jj) + Number(k1));
            csts += linear_constraint_t(e >= 0);
          }
          else if (jj == ZERO())
          {
            // ii <= k1
            linear_expression_t e(index_to_expr(ii) - Number(k1));
            csts += linear_constraint_t(e <= 0);
          }
          else
          {
            // ii - jj <= k1
            linear_expression_t e(index_to_expr(ii) - index_to_expr(jj));
              csts += linear_constraint_t(e <= Number(k1));
          }
        }
      }
    }
    return boost::optional<linear_constraint_system_t>(csts);
  }

}; // class DBM
} // namespace ikos

#endif // IKOS_DBM_HPP
