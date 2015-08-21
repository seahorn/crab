/*******************************************************************************
 *
 * Difference Bounds Matrix domain
 *
 * Based on the paper "Fast and Flexible Difference Constraint
 * Propagation for DPLL(T) by Cotton and Maler.
 *
 * C++ wrapper for a C implementation written by Graeme Gange
 * (gkgange@unimelb.edu.au)
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
   inline int ntoi(z_number k){
     ostringstream buf;
     buf << k;
     return boost::lexical_cast<int>(buf.str());
   }

   template<typename Number>
   Number iton(int n);

   template<>
   inline z_number iton(int n){
     return z_number(n);
   }
}; // end namespace DBM_impl

template< typename Number, typename VariableName, std::size_t Size = 3000 >
class DBM: public writeable,
           public numerical_domain<Number, VariableName >,
           public bitwise_operators< Number, VariableName >,
           public division_operators< Number, VariableName >{
  
 private:

  struct Index: public writeable {
    Index(index_t idx): writeable(), _idx(idx){ }
    Index(int idx): writeable(), _idx(idx){ }
    Index(const Index &o): writeable(), _idx(o._idx) { }
    Index& operator=(Index o){
      _idx = o._idx;
      return *this;
    }
    index_t index(){ return _idx; } 
    bool operator<(const Index&o) const { return (_idx < o._idx); }
    ostream& write(ostream& o) {
      o << _idx;
      return o;
    }
    index_t _idx;
  };
  
 public:
  typedef variable< Number, VariableName > variable_t;
  typedef linear_constraint< Number, VariableName > linear_constraint_t;
  typedef typename linear_constraint_t::kind_t constraint_kind_t;
  typedef linear_constraint_system< Number, VariableName > linear_constraint_system_t;
  typedef linear_expression< Number, VariableName > linear_expression_t;
  typedef interval< Number >  interval_t;
  typedef DBM< Number, VariableName > DBM_t;

 private:
  typedef bound<Number >  bound_t;
  typedef interval_domain< Number, VariableName > intervals_t;
  typedef mergeable_map< Index, VariableName > rev_map_t;    

 private:
  dbm _dbm;
  rev_map_t _rev_map;

  // return true if edge (i,j) exists in the DBM
  bool has_edge(int i, int j) {
    return (src_is_live(_dbm, i) && 
            in_graph(_dbm, i, j));
  }
                         
  // special variables
  inline int getZero(){ return (int) Size-1; }
  inline int getTmp() { return (int) Size-2; }
    
  void add_variable(VariableName x){
    int k = x.index();
    if (!_rev_map[k]){
      if ((unsigned) k > (Size - 2))
        IKOS_ERROR("DBM: need to enlarge the matrix");
      else
        _rev_map.set(k,x);
    }
  }

 private:

  DBM(bool is_bottom): 
    writeable(), _dbm(dbm_bottom()), _rev_map() { 
    if (!is_bottom)
      _dbm = dbm_top(Size);
  }

  DBM(dbm dbm, rev_map_t map): 
      writeable(), _dbm(dbm), _rev_map(map) { }

  void swap(dbm& x, dbm& y) {
    dbm tmp = x;
    x = y;
    y = tmp;
  }

  void assign(VariableName x, exp_t exp) {
    dbm ret = NULL;      
    ret = dbm_assign(x.index(), exp, _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
    exp_dealloc(exp);
  }

  void assign_tmp(VariableName x) {
    // pre: forget(getTmp())
    dexpr e1;     
    e1.args[0] = x.index();
    e1.args[1] = getTmp();
    e1.konst   = 0;
    e1.kind    = D_DIFF;
    apply_dexpr(e1);
    dexpr e2;     
    e2.args[0] = getTmp();
    e2.args[1] = x.index();
    e2.konst   = 0;
    e2.kind    = D_DIFF;
    apply_dexpr(e2);
  }

  void apply_cond(ucon con) {
    dbm ret = NULL;
    ret = dbm_cond(con, _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
  }

  void apply_dexpr(dexpr d) {
    dbm ret = NULL;
    ret = dbm_apply_dexpr(d, _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
  }

  void forget(index_t idx) {
    dbm ret = NULL;
    ret = dbm_forget(idx, _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
  }
    
  void forget(vector<int> idxs) {
    dbm ret = NULL;
    ret = dbm_forget_array(&idxs[0], idxs.size(), _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
  }

  void set_to_bottom() {
    dbm_dealloc(_dbm);
    _dbm = dbm_bottom();
  }

  VariableName map_index(index_t k) {
    optional<VariableName> v = _rev_map[k];
    if (!v)
      IKOS_ERROR("DBM: index cannot be mapped to a variable name") ;
    else 
      return *v;
  }

  void apply(operation_t op, VariableName x, index_t y, index_t z) { 	
    switch(op){
      case OP_ADDITION:
        {
          exp_t exp = exp_add(exp_var(y), exp_var(z));
          assign(x, exp);
        }
        break;
      case OP_SUBTRACTION:
        {
          exp_t exp = exp_sub(exp_var(y), exp_var(z));
          assign(x, exp);
        }
        break;
      case OP_MULTIPLICATION:
        {
          exp_t exp = exp_mul(exp_var(y), exp_var(z));
          assign(x, exp);
        }
        break;
      case OP_DIVISION:
        {
          exp_t exp = exp_div(exp_var(y), exp_var(z));
          assign(x, exp);
        }
        break;
      default:;;
    }
  }	

  void apply(operation_t op, VariableName x, index_t y, Number k) {	
    int z = DBM_impl::ntoi<Number>(k);
    switch(op){
      case OP_ADDITION:
        {
          exp_t exp = exp_add(exp_var(y), exp_const(z));
          assign(x, exp);
        }
        break;
      case OP_SUBTRACTION:
        {
          exp_t exp = exp_sub(exp_var(y), exp_const(z));
          assign(x, exp);
        }
        break;
      case OP_MULTIPLICATION:
        {
          exp_t exp = exp_mul(exp_var(y), exp_const(z));
          assign(x, exp);
        }
        break;
      case OP_DIVISION:
        {
          exp_t exp = exp_div(exp_var(y), exp_const(z));
          assign(x, exp);          
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

    auto vars = cst.variables();
    intervals_t inv;
    for(auto v: vars)
      inv.set (v.name(), operator[](v.name()));
    inv += cst;
    return !inv.is_bottom();
  }

  // x is unit coefficient
  void add_constraint(int coef_x, variable_t x, Number c, 
                      linear_constraint_t cst) {

    if (!(coef_x == 1 || coef_x == -1))
      IKOS_ERROR("DBM: coefficients can be only 1 or -1");

    add_variable(x.name());
    if (coef_x == 1){
      uterm tx = uvar(x.name().index());
      uterm ty = uconst(DBM_impl::ntoi<Number>(c));
      ucon con;
      if (cst.is_inequality()){       // x <= c
        con = mk_ucon(tx, U_LEQ, ty);
        apply_cond(con);
      }
      else if (cst.is_equality()){    // x == c
        con = mk_ucon(tx, U_EQ, ty);
        apply_cond(con);
      }
      else if (cst.is_disequation()){ // x != c
        // DBM ignores disequations so we better use intervals
        if (!intervals_check_sat(cst)){ 
          set_to_bottom(); 
        }
      }
    }
    else{ 
      if (cst.is_inequality()){       // 0-x <= c
        dexpr e;     
        e.args[0] = getZero();
        e.args[1] = x.name().index();
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
      else if (cst.is_equality()){    // -x == c iff x == -c
        uterm tx = uvar(x.name().index());
        uterm tc = uconst(-DBM_impl::ntoi<Number>(c));
        ucon con = mk_ucon(tx, U_EQ, tc);
        apply_cond(con);
      }
      else if (cst.is_disequation()){ // -x != c iff x != -c
        // DBM ignores disequations so we better use intervals
        if (!intervals_check_sat(cst)){ 
          set_to_bottom(); 
        }
      }
    }
  }

  // x, y are unit coefficient
  void add_constraint(int coef_x, variable_t x, int coef_y, variable_t y, 
                      Number c, linear_constraint_t cst) {

    if (!(coef_x == 1 || coef_x == -1))
      IKOS_ERROR("DBM: coefficients can be only 1 or -1");
    if (!(coef_y == 1 || coef_y == -1))
      IKOS_ERROR("DBM: coefficients can be only 1 or -1");
    if (coef_x == coef_y)
      IKOS_ERROR("DBM: same coefficients");        

    add_variable(x.name());
    add_variable(y.name());

    if (coef_x == 1 && coef_y == -1){  
      if (cst.is_inequality()){       // x - y <= c
        dexpr e;     
        e.args[0] = x.name().index();
        e.args[1] = y.name().index();
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
      else if (cst.is_equality()){    // x - y == c iff x - y <= c and y -x <= -c
        dexpr e1;     
        e1.args[0] = x.name().index();
        e1.args[1] = y.name().index();
        e1.konst   = DBM_impl::ntoi<Number>(c);
        e1.kind    = D_DIFF;
        apply_dexpr(e1);
        dexpr e2;     
        e2.args[0] = y.name().index();
        e2.args[1] = x.name().index();
        e2.konst   = -DBM_impl::ntoi<Number>(c);
        e2.kind    = D_DIFF;
        apply_dexpr(e2);
      }
      else if (cst.is_disequation())
      {  // x - y != c
        if (has_edge(x.name ().index (), y.name ().index ()) && 
            has_edge(y.name ().index (), x.name ().index ()))
        {
          int k1 = edge_val(_dbm, x.name ().index (), y.name ().index ());
          int k2 = edge_val(_dbm, y.name ().index (), x.name ().index ());
          if ( (k1 + k2 == 0) && (k1 == DBM_impl::ntoi<Number>(c)))
          {
            set_to_bottom ();
          }
        }
        else if (!intervals_check_sat(cst)) {
          set_to_bottom(); 
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
        apply_dexpr(e);
      }
      else if (cst.is_equality()){    // y - x == c iff y - x <= c and x -y <= -c
        dexpr e1;     
        e1.args[0] = y.name().index();
        e1.args[1] = x.name().index();
        e1.konst   = DBM_impl::ntoi<Number>(c);
        e1.kind    = D_DIFF;
        apply_dexpr(e1);
        dexpr e2;     
        e2.args[0] = x.name().index();
        e2.args[1] = y.name().index();
        e2.konst   = -DBM_impl::ntoi<Number>(c);
        e2.kind    = D_DIFF;
        apply_dexpr(e2);
      }
      else if (cst.is_disequation())
      { // y - x != c
        if (has_edge(y.name ().index (), x.name ().index ()) && 
            has_edge(x.name ().index (), y.name ().index ()))
        {
          int k1 = edge_val(_dbm, y.name ().index (), x.name ().index ());
          int k2 = edge_val(_dbm, x.name ().index (), y.name ().index ());
          if ( (k1 + k2 == 0) && (k1 == DBM_impl::ntoi<Number>(c)))
          {
            set_to_bottom ();
          }
        }
        else if (!intervals_check_sat(cst)) {
          set_to_bottom(); 
        }
      }
    }
  }

  linear_expression_t index_to_expr(index_t k) {
    if (k == static_cast<index_t>(getZero()))
      return linear_expression_t(Number (0));
    else
      return linear_expression_t(map_index (k));
  }
  

 public:

  static DBM_t top() { return DBM(false); }
    
  static DBM_t bottom() { return DBM(true); }

  ~DBM(){ dbm_dealloc(_dbm); }
    
 public:

  DBM(): writeable(), _dbm(dbm_top(Size)), _rev_map() {}  
           
  DBM(const DBM_t& o): 
    writeable(), 
    numerical_domain<Number, VariableName >(),
    bitwise_operators< Number, VariableName >(),
    division_operators< Number, VariableName >(),
    _dbm(dbm_copy(o._dbm)), _rev_map(o._rev_map) { }
   
  DBM_t operator=(DBM_t o) {
    dbm_dealloc(_dbm);
    _dbm = dbm_copy(o._dbm); 
    _rev_map = o._rev_map;
    return *this;
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
      return dbm_is_leq(_dbm, o._dbm);
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
      return DBM_t(dbm_join(_dbm, o._dbm), 
                   _rev_map | o._rev_map);
    }
  } 

  DBM_t operator||(DBM_t o) {	
    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else
      return DBM_t(dbm_widen(_dbm, o._dbm), 
                   _rev_map | o._rev_map);
  } 

  DBM_t operator&(DBM_t o) { 
    if (is_bottom() || o.is_bottom())
      return bottom();
    else if (is_top())
      return o;
    else if (o.is_top())
      return *this;
    else{
      return DBM_t(dbm_meet(_dbm, o._dbm), 
                   _rev_map | o._rev_map);
    }
  }	
    
  DBM_t operator&&(DBM_t o) {	
    return *this & o ;
  }	

  void normalize() {
    dbm_canonical(_dbm);
  }

  void operator-=(VariableName v) {
    forget(v.index());
  }

  template<typename Iterator>
  void forget (Iterator it, Iterator et) {
    // transform
    vector<int> idxs (std::distance (it,et));
    for(int i=0; it!=et; ++it, ++i)
      idxs[i] = (int) (*it).index();
    // --- this is more efficient than calling operator-=
    //     multiple times.
    forget (idxs);
  }

  void assign(VariableName x, linear_expression_t e) {

    if(is_bottom())
      return;

    add_variable(x);
    if (e.is_constant()){
      exp_t exp = exp_const(DBM_impl::ntoi<Number>(e.constant()));
      assign(x,exp);
    }
    else if  (boost::optional<variable_t> v = e.get_variable()){
      VariableName y = (*v).name();
      if (!(x==y)){
        add_variable(y);
        exp_t exp = exp_var(y.index());        
        assign(x,exp);
      }
    }
    else
      IKOS_ERROR("DBM: only supports constant/variable on the rhs of assignment");
  }

  void apply(operation_t op, VariableName x, VariableName y, VariableName z){	

    add_variable(x);
    add_variable(y);
    add_variable(z);

    if (x == y){
      // to make sure that lhs does not appear on the rhs
      assign_tmp(y); 
      apply(op, x, getTmp(), z.index());
      forget(getTmp());
    }
    else if (x == z){
      // to make sure that lhs does not appear on the rhs
      assign_tmp(z); 
      apply(op, x, y.index(), getTmp());
      forget(getTmp());
    }
    else{
      if (x == y && y == z)
        IKOS_ERROR("DBM: does not support x := x + x ");
      else
        apply(op, x, y.index(), z.index());
    }
  }
    
  void apply(operation_t op, VariableName x, VariableName y, Number k) {	
    add_variable(x);
    add_variable(y);

    if (x == y){
      // to make sure that lhs does not appear on the rhs
      assign_tmp(y);
      apply(op, x, getTmp(), k);
      forget(getTmp());
    }
    else{
      apply(op, x, y.index(), k);
    }
  }
   
   
  void operator+=(linear_constraint_t cst) {  
    if(is_bottom())
      return;

    if (cst.is_tautology())
      return;
      
    if (cst.is_contradiction()){
      set_to_bottom();
      return ;
    }

    linear_expression_t exp = cst.expression();
    if (exp.size() > 2)
      IKOS_ERROR("DBM supports constraints with at most two variables");

    if (exp.size() == 0){
      cout << cst << endl;
      IKOS_ERROR("DBM: bad-formed constraint");
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
    
  void operator+=(linear_constraint_system_t csts) {  
    for(auto cst: csts) {
      operator+=(cst);
    }
  }

  interval_t operator[](VariableName x) { 

    if (is_top())    return interval_t::top();
    if (is_bottom()) return interval_t::bottom();

    int k = x.index();
    if (!_rev_map[k])
      return interval_t::top();
    else{
      bound_t lb("-oo");
      if (has_edge(getZero(), k))
        lb = bound_t(DBM_impl::iton<Number>(-edge_val(_dbm, getZero(), k)));

      bound_t ub("+oo");
      if (has_edge(k, getZero()))
        ub = bound_t(DBM_impl::iton<Number>(edge_val(_dbm, k, getZero())));

      interval_t i(lb,ub);
      return i;
    }
  }

  void set(VariableName x, interval_t intv) {
    int k = x.index();
    if (_rev_map[k])
      operator-=(x);
    else
      add_variable(x);
      
    if (!intv.is_top()){
      boost::optional<Number> lb = intv.lb().number();
      boost::optional<Number> ub = intv.ub().number();
      if (ub){
        // x - 0 <= ub
        dexpr e;     
        e.args[0] = k;
        e.args[1] = getZero();
        e.konst   = DBM_impl::ntoi<Number>(*ub);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
      if (lb){
        // x >= lb iff 0-x <= -lb 
        dexpr e;     
        e.args[0] = getZero();
        e.args[1] = k;
        e.konst   = -DBM_impl::ntoi<Number>(*lb);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
    }
  }

  // bitwise_operators_api
    
  void apply(conv_operation_t op, VariableName x, VariableName y, unsigned width) {
    // since reasoning about infinite precision we simply assign and
    // ignore the width.
    assign(x, linear_expression_t(y));
  }

  void apply(conv_operation_t op, VariableName x, Number k, unsigned width) {
    // since reasoning about infinite precision we simply assign
    // and ignore the width.
    assign(x, k);
  }

  void apply(bitwise_operation_t op, VariableName x, VariableName y, VariableName z) {
    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
    interval_t zi = operator[](z);
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
          IKOS_ERROR("DBM: unreachable");
    }
    set(x, xi);
  }
    
  void apply(bitwise_operation_t op, VariableName x, VariableName y, Number k) {
    // Convert to intervals and perform the operation
    interval_t yi = operator[](y);
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
          IKOS_ERROR("DBM: unreachable");
    }
    set(x, xi);
  }
    
  // division_operators_api
    
  void apply(div_operation_t op, VariableName x, VariableName y, VariableName z) {
    if (op == OP_SDIV){
      apply(OP_DIVISION, x, y, z);
    }
    else{
      // Convert to intervals and perform the operation
      interval_t yi = operator[](y);
      interval_t zi = operator[](z);
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
            IKOS_ERROR("DBM: unreachable");
      }
      set(x, xi);
    }
  }

  void apply(div_operation_t op, VariableName x, VariableName y, Number k) {
    if (op == OP_SDIV){
      apply(OP_DIVISION, x, y, k);
    }
    else{
      // Convert to intervals and perform the operation
      interval_t yi = operator[](y);
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
            IKOS_ERROR("DBM: unreachable");
      }
      set(x, xi);
    }
  }

  //! copy of x into a new fresh variable y
  void expand (VariableName x, VariableName y) {
    dbm ret = NULL;      
    ret = dbm_expand(x.index(), y.index (), _dbm);
    add_variable(y);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
  }

  // Output function
  ostream& write(ostream& o) { 
    dbm_canonical(_dbm);

#if 0
    // print internal representation
    cout << endl;
    cout << "map: " <<  _rev_map << endl;
    cout << "csts: ";
    dbm_print_to(cout, _dbm);
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
      linear_constraint_system_t inv = to_linear_constraint_system ();
      o << inv;
      return o;
    }
  }

  const char* getDomainName () const {return "DBM";}

  linear_constraint_system_t to_linear_constraint_system () {
    dbm_canonical(_dbm);
    linear_constraint_system_t csts;
    
    if(is_bottom ())
    {
      csts += linear_constraint_t (linear_expression_t (Number(1)) == 
                                     linear_expression_t (Number(0)));
      return csts;
    }

    boost::unordered_set< pair<int,int> > visited; // to only consider half matrix
    for(edge_iter iter=edge_iterator(_dbm); !srcs_end(iter); next_src(iter))
    {
      int ii = src(iter);
      if (ii == getTmp()) continue;
        
      for(; !dests_end(iter); next_dest(iter))
      {
        int jj = dest(iter);

        if (jj == getTmp()) continue;
        if (visited.find(make_pair(jj,ii)) != visited.end()) continue;
          
        int k1 = edge_val(_dbm, ii, jj);
        if (in_graph(_dbm, jj, ii))
        {
          int k2 = edge_val(_dbm, jj, ii);
          if ((k1 + k2) == 0)
          {
            // EQUALITY 
            if (k1 == 0)
            { // ii = jj
              csts += linear_constraint_t(index_to_expr(ii) == index_to_expr(jj));
            }
            else
            {
              if (ii == getZero())
              { // jj = -k1
                linear_expression_t e(index_to_expr(jj) + Number(k1));
                csts += linear_constraint_t(e == 0);
              }
              else if (jj == getZero())
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
            if (ii == getZero())
            { // jj >= -k1
              linear_expression_t e(index_to_expr(jj) + Number(k1));
              csts += linear_constraint_t(e >= 0);
            }
            else if (jj == getZero())
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
          if (ii == getZero())
          {// jj >= -k1
            linear_expression_t e(index_to_expr(jj) + Number(k1));
            csts += linear_constraint_t(e >= 0);
          }
          else if (jj == getZero())
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
    return csts;
  }

}; // class DBM

namespace domain_traits {

template <typename Number, typename VariableName>
void expand (DBM<Number, VariableName>& inv, 
             VariableName x, VariableName new_x) {
  inv.expand (x, new_x);
}

template <typename Number, typename VariableName>
void normalize (DBM<Number, VariableName>& inv) {
   inv.normalize();
}

template <typename Number, typename VariableName, typename Iterator >
void forget (DBM<Number, VariableName>& inv, Iterator it, Iterator end){
  inv.forget (it, end);
}

} // namespace domain_traits
} // namespace ikos

#endif // IKOS_DBM_HPP
