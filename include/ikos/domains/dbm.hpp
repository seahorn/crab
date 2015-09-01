/*******************************************************************************
 *
 * Difference Bounds Matrix domain
 *
 * Based on the paper "Fast and Flexible Difference Constraint
 * Propagation for DPLL(T)" by Cotton and Maler.
 *
 * C++ wrapper for a C implementation written by 
 * Graeme Gange (gkgange@unimelb.edu.au)
 ******************************************************************************/

#ifndef IKOS_DBM_HPP
#define IKOS_DBM_HPP

// Uncomment for enabling debug information
// #include <ikos/common/dbg.hpp>

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
#include <boost/unordered_map.hpp>
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

template< class Number, class VariableName, 
          unsigned Sz = 100, int Delta = 50> 
class DBM: public writeable,
           public numerical_domain<Number, VariableName >,
           public bitwise_operators<Number,VariableName >,
           public division_operators<Number, VariableName >{
  
 public:

  typedef variable< Number, VariableName> variable_t;
  typedef linear_constraint< Number, VariableName> linear_constraint_t;
  typedef typename linear_constraint_t::kind_t constraint_kind_t;
  typedef linear_constraint_system< Number, VariableName> linear_constraint_system_t;
  typedef linear_expression< Number, VariableName> linear_expression_t;
  typedef interval< Number>  interval_t;
  typedef DBM<Number,VariableName,Sz> DBM_t;

 private:
  typedef bound<Number>  bound_t;
  typedef interval_domain< Number, VariableName > intervals_t;

  typedef int id_t;
  //! map variable name to id 
  typedef boost::unordered_map<VariableName, id_t> var_map_t;
  //! map back id to variable name
  typedef boost::unordered_map <id_t, VariableName> rev_map_t;
  dbm _dbm;
  id_t _id;
  // class invariant: for all variable v. rev_map (var_map(v)) = v
  var_map_t _var_map;
  rev_map_t _rev_map;
  int _max_sz; //maximum allocate size of the dbm matrix

  // return true if edge (i,j) exists in the DBM
  bool has_edge(int i, int j) {
    return (src_is_live(_dbm, i) && in_graph(_dbm, i, j));
  }

  bool is_live (dbm x, int i){
    return (src_is_live(x, i) || dest_is_live (x, i));
  }

  // special variables
  inline int get_zero(){ return (int) _max_sz-1; }
  inline int get_tmp() { return (int) _max_sz-2; }

  //! Perform some some sort of defragmentation by removing matrix
  //  indexes which are not alive. A similar process is done by
  //  merge_var_map. Return true if the defragmentation took place.
  bool dbm_compact () {
     
    if (_id < (int) _max_sz/2 /*roughly half of the matrix size */)  
      return true;
    
    // assign a new id for each variable starting from 0
    id_t id = 0;
    std::set<VariableName> all_live_keys;
    for (auto px: _var_map) {
      if (is_live(_dbm, px.second))
        all_live_keys.insert (px.first);
    }

    if (all_live_keys.size () >= _max_sz  -2) {
      return false;
    }
      
    var_map_t res;
    for (auto k: all_live_keys)
      res.insert (make_pair (k, id++));

    // build the substitution map from the new var map
    vector<rmap> subs_x;
    for (auto const &px: res) {
      if (is_live (_dbm, _var_map[px.first]))
        subs_x.push_back (rmap {_var_map[px.first], px.second});      
    }

    // Update the this' id counter 
    _id = id;
    // Update the this' var map
    _var_map = res;
    // build a new this' reverse map
    _rev_map.clear ();
    for(auto p: _var_map)
      _rev_map.insert (make_pair (p.second, p.first));
    // rename the this' dbm
    dbm ret = NULL;
    ret= dbm_rename (&subs_x[0], subs_x.size(), _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);

    IKOS_DEBUG ("Compacted dbm to size ", _var_map.size (), "\n", *this);
    return true;
  }

  //! Unify the mappings from two different dbms. It performs also
  //! some defragmentation.
  var_map_t merge_var_map (var_map_t &x, 
                           dbm& dbm_x,
                           var_map_t &y, 
                           dbm& dbm_y,
                           id_t &id,
                           vector<rmap>& subs_x, 
                           vector<rmap>& subs_y) { 
    // assign a new id for each variable starting from 0
    id = 0;
    std::set<VariableName> all_keys;
    for (auto px: x) {
      if (is_live(dbm_x, px.second))
        all_keys.insert (px.first);
    }
    for (auto py: y) {
      if (is_live(dbm_y, py.second))
        all_keys.insert (py.first);
    }

    var_map_t res;
    for (auto k: all_keys) {
      // --- if id >= _max_sz then the matrix will be resized by the
      // --- caller.
      res.insert (make_pair (k, id));
      ++id;
    }

    // build the substitution maps
    for (auto const &px: x) {
      if (is_live(dbm_x, px.second))
        subs_x.push_back (rmap {px.second, res[px.first]});      
    }
    for (auto const &py: y) {
      if (is_live(dbm_y, py.second))
        subs_y.push_back (rmap {py.second, res[py.first]});      
    }
    
    return res;
  }

  //! Convert a variable name to a matrix index
  id_t get_dbm_index (VariableName x) {
    auto it = _var_map.find (x);
    if (it != _var_map.end ())
      return it->second;

    id_t k = _id++;
    if (k >= (int) (_max_sz - 2)) {
      // --- This might be possible if Delta is too small
      // IKOS_ERROR("DBM: need to enlarge the matrix. This should not happen!");
      resize (Delta);
    }
    
    _var_map.insert (make_pair (x, k));
    _rev_map.insert (make_pair (k, x));
    return k;
  }

  // Inverse of get_dbm_index
  VariableName get_rev_map (id_t k) {
    auto it = _rev_map.find (k);
    if (it!= _rev_map.end ())
      return it->second;

    IKOS_ERROR("DBM: index ", k," cannot be mapped to a variable name") ;
  }

  // k is the index of a variable produced by get_dbm_index
  linear_expression_t to_expr (id_t k) {
    if (k == get_zero())
      return linear_expression_t (Number (0));
    else
      return linear_expression_t (get_rev_map (k));
  }

 private:

  DBM(bool is_bottom): writeable(), _dbm(dbm_bottom()), _id(0), _max_sz (Sz) { 
    if (!is_bottom)
      _dbm = dbm_top(_max_sz);
  }

  DBM(dbm dbm, id_t id, var_map_t var_map, rev_map_t rev_map, unsigned max_sz): 
      writeable(), 
      _dbm(dbm), 
      _id (id),
      _var_map(var_map), 
      _rev_map(rev_map),
      _max_sz (max_sz) { 
    if (!_dbm)
      set_to_bottom ();
  }

 private:

  void swap(dbm& x, dbm& y) {
    dbm tmp = x;
    x = y;
    y = tmp;
  }

  void assign(VariableName x, exp_t exp) {
    dbm ret = NULL;      
    ret = dbm_assign(get_dbm_index(x), exp, _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
    exp_dealloc(exp);
  }

  // tmp := x
  void assign_tmp(VariableName x) {
    dexpr e1;     
    e1.args[0] = get_dbm_index(x);
    e1.args[1] = get_tmp();
    e1.konst   = 0;
    e1.kind    = D_DIFF;
    apply_dexpr(e1);
    dexpr e2;     
    e2.args[0] = get_tmp();
    e2.args[1] = get_dbm_index(x);
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

  void forget(id_t idx) {
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
    if (_dbm)
      dbm_dealloc (_dbm);
    _dbm = dbm_bottom ();
    _var_map.clear ();
    _rev_map.clear ();
  }

  void apply(operation_t op, VariableName x, id_t y, id_t z) { 	
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

  void apply(operation_t op, VariableName x, id_t y, Number k) {	
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
      return;

    if (coef_x == 1){
      uterm tx = uvar(get_dbm_index (x.name()));
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
        e.args[0] = get_zero();
        e.args[1] = get_dbm_index (x.name());
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
      else if (cst.is_equality()){    // -x == c iff x == -c
        uterm tx = uvar(get_dbm_index (x.name()));
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

    if (coef_x == 1 && coef_y == -1){  
      if (cst.is_inequality()){       // x - y <= c
        dexpr e;     
        e.args[0] = get_dbm_index(x.name());
        e.args[1] = get_dbm_index(y.name());
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
      else if (cst.is_equality()){    // x - y == c iff x - y <= c and y -x <= -c
        dexpr e1;     
        e1.args[0] = get_dbm_index(x.name());
        e1.args[1] = get_dbm_index(y.name());
        e1.konst   = DBM_impl::ntoi<Number>(c);
        e1.kind    = D_DIFF;
        apply_dexpr(e1);
        dexpr e2;     
        e2.args[0] = get_dbm_index(y.name());
        e2.args[1] = get_dbm_index(x.name());
        e2.konst   = -DBM_impl::ntoi<Number>(c);
        e2.kind    = D_DIFF;
        apply_dexpr(e2);
      }
      else if (cst.is_disequation())
      {  // x - y != c
        if (has_edge(get_dbm_index(x.name ()), get_dbm_index(y.name ())) && 
            has_edge(get_dbm_index(y.name ()), get_dbm_index(x.name ())))
        {
          int k1 = edge_val(_dbm, get_dbm_index(x.name ()), get_dbm_index(y.name ()));
          int k2 = edge_val(_dbm, get_dbm_index(y.name ()), get_dbm_index(x.name ()));
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
        e.args[0] = get_dbm_index(y.name());
        e.args[1] = get_dbm_index(x.name());
        e.konst   = DBM_impl::ntoi<Number>(c);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
      else if (cst.is_equality()){    // y - x == c iff y - x <= c and x -y <= -c
        dexpr e1;     
        e1.args[0] = get_dbm_index(y.name());
        e1.args[1] = get_dbm_index(x.name());
        e1.konst   = DBM_impl::ntoi<Number>(c);
        e1.kind    = D_DIFF;
        apply_dexpr(e1);
        dexpr e2;     
        e2.args[0] = get_dbm_index(x.name());
        e2.args[1] = get_dbm_index(y.name());
        e2.konst   = -DBM_impl::ntoi<Number>(c);
        e2.kind    = D_DIFF;
        apply_dexpr(e2);
      }
      else if (cst.is_disequation())
      { // y - x != c
        if (has_edge(get_dbm_index(y.name ()), get_dbm_index(x.name ())) && 
            has_edge(get_dbm_index(x.name ()), get_dbm_index(y.name ())))
        {
          int k1 = edge_val(_dbm, get_dbm_index(y.name ()), get_dbm_index(x.name ()));
          int k2 = edge_val(_dbm, get_dbm_index(x.name ()), get_dbm_index(y.name ()));
          if ( (k1 + k2 == 0) && (k1 == DBM_impl::ntoi<Number>(c))) {
            set_to_bottom ();
          }
        }
        else if (!intervals_check_sat(cst)) {
          set_to_bottom(); 
        }
      }
    }
  }

  // Return a new db with all edges of d but new size sz
  dbm resize (dbm d, int sz) {
    if (dbm_is_bottom (d))
      return d;

    if (sz < 3 || !(sz > d->sz))
      return d;

    IKOS_DEBUG ("Resizing the matrix to ", sz, " ... ");

    dbm x = NULL;      
    x= dbm_resize (d, sz);

    // rename special variable "0"
    vector<rmap>  subs;
    subs.push_back (rmap {d->sz-1,sz-1});
    dbm ret = NULL;
    ret= dbm_rename (&subs[0], subs.size(), x);

    dbm_dealloc (x);
    return ret;
  }

  // resize this->_dbm to this->_max_size + delta
  void resize (int delta) {
    int max_sz = _max_sz + delta;
    dbm tmp = resize (_dbm, max_sz);
    dbm_dealloc(_dbm);    
    swap(_dbm, tmp);
    _max_sz = max_sz;
  }

 public:

  static DBM_t top() { return DBM(false); }
    
  static DBM_t bottom() { return DBM(true); }
    
  DBM(): writeable(), _dbm(dbm_top(Sz)), _id (0), _max_sz (Sz) {}  
           
  DBM(const DBM_t& o): 
    writeable(), 
    numerical_domain<Number, VariableName >(),
    bitwise_operators< Number, VariableName >(),
    division_operators< Number, VariableName >(),
    _dbm(dbm_copy(o._dbm)), 
    _id (o._id),
    _var_map (o._var_map),
    _rev_map(o._rev_map),
    _max_sz (o._max_sz)
  { }
   
  DBM_t operator=(const DBM_t& o) {
    if (this != &o) {
      dbm_dealloc(_dbm);
      _dbm = dbm_copy(o._dbm); 
      _id = o._id;
      _var_map = o._var_map;
      _rev_map = o._rev_map;
      _max_sz = o._max_sz;
    }
    return *this;
  }

  ~DBM() { 
    dbm_dealloc(_dbm); 
    _var_map.clear();
    _rev_map.clear();
  }

 public:

  bool is_bottom() {
    return dbm_is_bottom(_dbm);
  }
    
  bool is_top() {
    return dbm_is_top(_dbm);
  }
    
  bool operator<=(DBM_t o)  {       
    // cover all trivial cases to avoid allocating a dbm matrix
    if (is_bottom()) 
      return true;
    else if(o.is_bottom())
      return false;
    else if (o.is_top ())
      return true;
    else if (is_top () && !o.is_top ())
      return false;
    else if (is_top () && o.is_top ())
      return true;
    else { 

      vector<rmap>  subs_x;
      vector<rmap>  subs_y;
      id_t id = 0;
      var_map_t var_map = merge_var_map (_var_map, _dbm, 
                                         o._var_map, o._dbm,
                                         id, subs_x, subs_y);

      // Build the reverse map
      rev_map_t rev_map;
      for(auto p: var_map)
        rev_map.insert (make_pair (p.second, p.first));

      dbm dbm_x = NULL;
      dbm dbm_y = NULL;

      dbm_x = dbm_rename (&subs_x[0], subs_x.size(), _dbm);
      dbm_y = dbm_rename (&subs_y[0], subs_y.size(), o._dbm);

      int max_sz = std::max (_max_sz, o._max_sz);
      if (max_sz != _max_sz) {
        dbm tmp = resize (dbm_x, max_sz);
        dbm_dealloc(dbm_x);    
        swap(dbm_x, tmp);
      }
      if (max_sz != o._max_sz) {
        dbm tmp = resize (dbm_y, max_sz);
        dbm_dealloc(dbm_y);    
        swap(dbm_y, tmp);
      }
        
      assert (dbm_x->sz == dbm_y->sz);            
      bool res = dbm_is_leq (dbm_x, dbm_y);
      dbm_dealloc(dbm_x);
      dbm_dealloc(dbm_y);

      return res;
    }
  }  

  DBM_t operator|(DBM_t o) {
    if (is_bottom() || o.is_top ())
      return o;
    else if (is_top () || o.is_bottom())
      return *this;
    else {
      IKOS_DEBUG ("Before join:\n","DBM 1\n",*this,"\n","DBM 2\n",o);

      vector<rmap>  subs_x;
      vector<rmap>  subs_y;
      id_t id = 0;
      var_map_t var_map = merge_var_map (_var_map, _dbm, 
                                         o._var_map, o._dbm,
                                         id, subs_x, subs_y);

      // Build the reverse map
      rev_map_t rev_map;
      for(auto p: var_map)
        rev_map.insert (make_pair (p.second, p.first));

      dbm dbm_x = NULL;
      dbm dbm_y = NULL;

      dbm_x = dbm_rename (&subs_x[0], subs_x.size(), _dbm);
      dbm_y = dbm_rename (&subs_y[0], subs_y.size(), o._dbm);

#if 0     
      cout << "After resizing DBMs: \n";
      cout << "DBM 1:\n";
      dbm_print_to(cout, dbm_x);
      cout << "subs map {";
      for (auto p: subs_x)
        cout << "(" << p.r_from << "-" << p.r_to << ");";
      cout << "}\n";
      cout << "DBM 2:\n";
      cout << "subs map {";
      for (auto p: subs_y)
        cout << "(" << p.r_from << "-" << p.r_to << ");";
      cout << "}\n";
      dbm_print_to(cout, dbm_y);
#endif 
      
      int max_sz = std::max (_max_sz, o._max_sz);
      if (max_sz != _max_sz) {
        dbm tmp = resize (dbm_x, max_sz);
        dbm_dealloc(dbm_x);    
        swap(dbm_x, tmp);
      }
      if (max_sz != o._max_sz) {
        dbm tmp = resize (dbm_y, max_sz);
        dbm_dealloc(dbm_y);    
        swap(dbm_y, tmp);
      }

      assert (dbm_x->sz == dbm_y->sz);
      DBM_t res (dbm_join(dbm_x, dbm_y), 
                 id, var_map, rev_map, 
                 max_sz);
      
      dbm_dealloc(dbm_x);
      dbm_dealloc(dbm_y);

      IKOS_DEBUG ("Result join:\n",res);

      return res;
    }
  } 

  DBM_t operator||(DBM_t o) {	
    if (is_bottom())
      return o;
    else if (o.is_bottom())
      return *this;
    else {

      IKOS_DEBUG ("Before widening:\n","DBM 1\n",*this,"\n","DBM 2\n",o);

      vector<rmap>  subs_x;
      vector<rmap>  subs_y;
      id_t id = 0;
      var_map_t var_map = merge_var_map (_var_map, _dbm,
                                         o._var_map, o._dbm,
                                         id, subs_x, subs_y);

      // Build the reverse map
      rev_map_t rev_map;
      for(auto p: var_map)
        rev_map.insert (make_pair (p.second, p.first));

      dbm dbm_x = NULL;
      dbm dbm_y = NULL;

      dbm_x = dbm_rename (&subs_x[0], subs_x.size(), _dbm);
      dbm_y = dbm_rename (&subs_y[0], subs_y.size(), o._dbm);

      int max_sz = std::max (_max_sz, o._max_sz);
      if (max_sz != _max_sz) {
        dbm tmp = resize (dbm_x, max_sz);
        dbm_dealloc(dbm_x);    
        swap(dbm_x, tmp);
      }

      if (max_sz != o._max_sz) {
        dbm tmp = resize (dbm_y, max_sz);
        dbm_dealloc(dbm_y);    
        swap(dbm_y, tmp);
      }
                    
      assert (dbm_x->sz == dbm_y->sz);
      DBM_t res (dbm_widen(dbm_x, dbm_y), 
                 id, var_map, rev_map, 
                 max_sz);
      
      dbm_dealloc(dbm_x);
      dbm_dealloc(dbm_y);

      IKOS_DEBUG ("Result widening:\n",res);
      return res;
    }
  } 

  DBM_t operator&(DBM_t o) { 
    if (is_bottom() || o.is_bottom())
      return bottom();
    else if (is_top())
      return o;
    else if (o.is_top())
      return *this;
    else{

      IKOS_DEBUG ("Before meet:\n","DBM 1\n",*this,"\n","DBM 2\n",o);

      vector<rmap>  subs_x;
      vector<rmap>  subs_y;
      id_t id = 0;
      var_map_t var_map = merge_var_map (_var_map, _dbm,
                                         o._var_map, o._dbm,
                                         id, subs_x, subs_y);

      // Build the reverse map
      rev_map_t rev_map;
      for(auto p: var_map)
        rev_map.insert (make_pair (p.second, p.first));

      dbm dbm_x = NULL;
      dbm dbm_y = NULL;

      dbm_x = dbm_rename (&subs_x[0], subs_x.size(), _dbm);
      dbm_y = dbm_rename (&subs_y[0], subs_y.size(), o._dbm);

      int max_sz = std::max (_max_sz, o._max_sz);
      if (max_sz != _max_sz) {
        dbm tmp = resize (dbm_x, max_sz);
        dbm_dealloc(dbm_x);    
        swap(dbm_x, tmp);
      }

      if (max_sz != o._max_sz) {
        dbm tmp = resize (dbm_y, max_sz);
        dbm_dealloc(dbm_y);    
        swap(dbm_y, tmp);
      }
                    
      assert (dbm_x->sz == dbm_y->sz);
      DBM_t res (dbm_meet(dbm_x, dbm_y), 
                 id, var_map, rev_map, 
                 max_sz); 

      dbm_dealloc(dbm_x);
      dbm_dealloc(dbm_y);

      IKOS_DEBUG ("Result meet:\n",res);
      return res;
    }
  }	
    
  DBM_t operator&&(DBM_t o) {	
    if (is_bottom() || o.is_bottom())
      return bottom();
    else if (is_top ())
      return o;
    else{
      IKOS_DEBUG ("Before narrowing:\n","DBM 1\n",*this,"\n","DBM 2\n",o);

      vector<rmap>  subs_x;
      vector<rmap>  subs_y;
      id_t id = 0;
      var_map_t var_map = merge_var_map (_var_map, _dbm,
                                         o._var_map, o._dbm,
                                         id, subs_x, subs_y);

      // Build the reverse map
      rev_map_t rev_map;
      for(auto p: var_map)
        rev_map.insert (make_pair (p.second, p.first));

      dbm dbm_x = NULL;
      dbm dbm_y = NULL;

      dbm_x = dbm_rename (&subs_x[0], subs_x.size(), _dbm);
      dbm_y = dbm_rename (&subs_y[0], subs_y.size(), o._dbm);

      int max_sz = std::max (_max_sz, o._max_sz);
      if (max_sz != _max_sz) {
        dbm tmp = resize (dbm_x, max_sz);
        dbm_dealloc(dbm_x);    
        swap(dbm_x, tmp);
      }
      if (max_sz != o._max_sz) {
        dbm tmp = resize (dbm_y, max_sz);
        dbm_dealloc(dbm_y);    
        swap(dbm_y, tmp);
      }
            
      assert (dbm_x->sz == dbm_y->sz);
      DBM_t res (dbm_narrowing(dbm_x, dbm_y), 
                 id, var_map, rev_map, 
                 max_sz); 

      dbm_dealloc(dbm_x);
      dbm_dealloc(dbm_y);

      IKOS_DEBUG ("Result narrowing:\n",res);
      return res;
    }
  }	

  void normalize() {
    dbm_canonical(_dbm);
  }

  void operator-=(VariableName v) {
    if (is_bottom ())
      return;

    auto it = _var_map.find (v);
    if (it != _var_map.end ()) {
      forget (it->second);
      _var_map.erase (v);
      _rev_map.erase (it->second);
    }
  }

  template<typename Iterator>
  void forget (Iterator vIt, Iterator vEt) {
    if (is_bottom ())
      return;
    vector<int> idxs;
    for (auto v: boost::make_iterator_range (vIt,vEt)) {
      auto it = _var_map.find (v);
      if (it != _var_map.end ()) {
        idxs.push_back  ((int) it->second);
        _var_map.erase (v);
        _rev_map.erase (it->second);
      }
    }
    if (!idxs.empty()) {
      // --- forget all at once is more efficient than calling
      //     operator-= multiple times.
      forget (idxs);
      if (!dbm_compact ())
        resize (Delta);
    }
  }

  void assign(VariableName x, linear_expression_t e) {

    if(is_bottom())
      return;

    if (!dbm_compact ())
      resize (Delta);
    
    if (e.is_constant()){
      exp_t exp = exp_const(DBM_impl::ntoi<Number>(e.constant()));
      assign (x,exp);
    }
    else if  (optional<variable_t> v = e.get_variable()){
      VariableName y = (*v).name();
      if (!(x==y)){
        exp_t exp = exp_var (get_dbm_index(y));        
        assign (x,exp);
      }
    }
    else
      IKOS_ERROR("DBM: only supports constant/variable on the rhs of assignment");

    IKOS_DEBUG("---", x, ":=", e,"\n",*this);
  }

  void apply(operation_t op, VariableName x, VariableName y, VariableName z){	

    if(is_bottom())
      return;

    if (!dbm_compact ())
      resize (Delta);

    if (x == y){
      // --- ensure lhs does not appear on the rhs
      assign_tmp(y); 
      apply(op, x, get_tmp(), get_dbm_index(z));
      forget(get_tmp());
    }
    else if (x == z){
      // --- ensure lhs does not appear on the rhs
      assign_tmp(z); 
      apply(op, x, get_dbm_index (y), get_tmp());
      forget(get_tmp());
    }
    else{
      if (x == y && y == z)
        IKOS_ERROR("DBM: does not support x := x + x ");
      else
        apply(op, x, get_dbm_index(y), get_dbm_index(z));
    }

    IKOS_DEBUG("---", x, ":=", y, op, z,"\n", *this);
  }

    
  void apply(operation_t op, VariableName x, VariableName y, Number k) {	

    if(is_bottom())
      return;

    if (!dbm_compact ())
      resize (Delta);

    if (x == y){
      // to make sure that lhs does not appear on the rhs
      assign_tmp(y);
      apply(op, x, get_tmp(), k);
      forget(get_tmp());
    }
    else{
      apply(op, x, get_dbm_index(y), k);
    }

    IKOS_DEBUG("---", x, ":=", y, op, k,"\n", *this);
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

    if (exp.size() == 0)
      IKOS_ERROR("DBM: bad-formed constraint: ", cst);

    if (!dbm_compact ())
      resize (Delta);

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

    IKOS_DEBUG("---", cst, "\n", *this);
  } 
    
  void operator+=(linear_constraint_system_t csts) {  
    if(is_bottom()) return;

    for(auto cst: csts) {
      operator+=(cst);
    }
  }

  interval_t operator[](VariableName x) { 

    if (is_top())    return interval_t::top();
    if (is_bottom()) return interval_t::bottom();

    int k = get_dbm_index(x);
    bound_t lb("-oo");
    if (has_edge(get_zero(), k))
      lb = bound_t(DBM_impl::iton<Number>(-edge_val(_dbm, get_zero(), k)));
    
    bound_t ub("+oo");
    if (has_edge(k, get_zero()))
      ub = bound_t(DBM_impl::iton<Number>(edge_val(_dbm, k, get_zero())));
    
    interval_t i(lb,ub);
    return i;
  }

  void set(VariableName x, interval_t intv) {
    if(is_bottom())
      return;

    int k = get_dbm_index(x);
    if (!intv.is_top()){
      optional<Number> lb = intv.lb().number();
      optional<Number> ub = intv.ub().number();
      if (ub){
        // x - 0 <= ub
        dexpr e;     
        e.args[0] = k;
        e.args[1] = get_zero();
        e.konst   = DBM_impl::ntoi<Number>(*ub);
        e.kind    = D_DIFF;
        apply_dexpr(e);
      }
      if (lb){
        // x >= lb iff 0-x <= -lb 
        dexpr e;     
        e.args[0] = get_zero();
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
    if(is_bottom())
      return;

    dbm ret = NULL;      
    ret = dbm_expand(get_dbm_index(x), get_dbm_index(y), _dbm);
    dbm_dealloc(_dbm);
    swap(_dbm, ret);
  }

  // Output function
  ostream& write(ostream& o) { 

    normalize ();

#ifdef _IKOS_DEBUG_
    cout << "var_map={";
    for (auto &p: _var_map) 
      cout << p.first << "(" << p.first.index () << ") " << "->" << p.second << ";";
    cout << "}\n";
    cout << "rev_map={";
    for (auto &p: _rev_map) 
      cout << p.first << "->" << p.second << ";";
    cout << "}\n";
    cout << "matrix:\n";
    dbm_print_to(cout, _dbm);
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

  linear_constraint_system_t to_linear_constraint_system () {

    normalize ();

    linear_constraint_system_t csts;
    
    if(is_bottom ()) {
      csts += linear_constraint_t (linear_expression_t (Number(1)) == 
                                   linear_expression_t (Number(0)));
      return csts;
    }

    // --- to only consider half matrix
    boost::unordered_set< pair<int,int> > visited; 
    for(edge_iter iter=edge_iterator(_dbm); !srcs_end(iter); next_src(iter))
    {
      int ii = src(iter);
      if (ii == get_tmp()) continue;
        
      for(; !dests_end(iter); next_dest(iter))
      {
        int jj = dest(iter);

        if (jj == get_tmp()) continue;
        if (visited.find(make_pair(jj,ii)) != visited.end()) continue;
          
        int k1 = edge_val(_dbm, ii, jj);
        if (has_edge (jj,ii))
        {
          int k2 = edge_val(_dbm, jj, ii);
          if ((k1 + k2) == 0)
          {
            // EQUALITY 
            if (k1 == 0)
            { // ii = jj
              csts += linear_constraint_t(to_expr(ii) == to_expr(jj));
            }
            else
            {
              if (ii == get_zero())
              { // jj = -k1
                linear_expression_t e(to_expr(jj) + Number(k1));
                csts += linear_constraint_t(e == 0);
              }
              else if (jj == get_zero())
              { // ii = k1
                linear_expression_t e(to_expr(ii) - Number(k1));
                csts += linear_constraint_t(e == 0);
              }
              else
              {
                // ii - jj = k1
                linear_expression_t e(to_expr(ii) - to_expr(jj));
                csts += linear_constraint_t(e == Number(k1));
              }
            }
            visited.insert(make_pair(ii,jj));
          }
          else
          {
            // INEQUALITY
            if (ii == get_zero())
            { // jj >= -k1
              linear_expression_t e(to_expr(jj) + Number(k1));
              csts += linear_constraint_t(e >= 0);
            }
            else if (jj == get_zero())
            {
              // ii <= k1
              linear_expression_t e(to_expr(ii) - Number(k1));
              csts += linear_constraint_t(e <= 0);
            }
            else
            {
              // ii - jj <= k1
              linear_expression_t e(to_expr(ii) - to_expr(jj));
              csts += linear_constraint_t(e <= Number(k1));
            }
          }
        }
        else
        {
          if (ii == get_zero())
          {// jj >= -k1
            linear_expression_t e(to_expr(jj) + Number(k1));
            csts += linear_constraint_t(e >= 0);
          }
          else if (jj == get_zero())
          {
            // ii <= k1
            linear_expression_t e(to_expr(ii) - Number(k1));
            csts += linear_constraint_t(e <= 0);
          }
          else
          {
            // ii - jj <= k1
            linear_expression_t e(to_expr(ii) - to_expr(jj));
              csts += linear_constraint_t(e <= Number(k1));
          }
        }
      }
    }
    return csts;
  }

  const char* getDomainName () const {return "DBM";}

}; // class DBM

namespace domain_traits {

template <typename Number, typename VariableName>
void expand (DBM<Number,VariableName>& inv, 
             VariableName x, VariableName new_x) {
  inv.expand (x, new_x);
}

template <typename Number, typename VariableName>
void normalize (DBM<Number,VariableName>& inv) {
   inv.normalize();
}

template <typename Number, typename VariableName, typename Iterator >
void forget (DBM<Number,VariableName>& inv, Iterator it, Iterator end){
  inv.forget (it, end);
}

} // namespace domain_traits
} // namespace ikos

#endif // IKOS_DBM_HPP
