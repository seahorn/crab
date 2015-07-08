#include <iostream>
#include <boost/tuple/tuple.hpp>
#include <ikos/common/types.hpp>
// Array graph
#include <ikos/domains/array_graph.hpp>
// Scalar/Weight domains
#include <ikos/algorithms/linear_constraints.hpp>
#include <ikos/common/bignums.hpp>
#include <ikos/domains/intervals.hpp>
#include <ikos/domains/octagons.hpp>
#include <ikos/domains/dbm.hpp>
// CFG and fixpoint
#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/FwdAnalyzer.hpp>

namespace cfg_impl
{
  using namespace cfg;

  template<> inline std::string get_label_str(std::string e) 
  { return e; }

  class StrVariableFactory : public boost::noncopyable  
  {
    typedef var_factory_impl::VariableFactory< std::string > StrVariableFactory_t;
    std::unique_ptr< StrVariableFactory_t > m_factory; 
    
   public: 

    typedef StrVariableFactory_t::variable_t varname_t;

    StrVariableFactory(): m_factory (new StrVariableFactory_t()){ }

    varname_t operator[](std::string v)

    { return (*m_factory)[v];}
  }; 

  // A variable factory based on strings
  typedef StrVariableFactory VariableFactory;
  typedef typename VariableFactory::varname_t varname_t;

  // CFG
  typedef variable< z_number, varname_t >      z_var;
  typedef std::string                          basic_block_label_t;
  typedef Cfg< basic_block_label_t, varname_t> cfg_t;
  typedef cfg_t::basic_block_t                 basic_block_t;
} // end namespace

namespace domain_impl
{
  using namespace cfg_impl;
  // Numerical domains  
  typedef interval_domain< z_number, varname_t > interval_domain_t;
  typedef DBM< z_number, varname_t > dbm_domain_t;
} // end namespace

using namespace cfg_impl;
using namespace domain_impl;
using namespace analyzer;


// using namespace std;
// using namespace ikos;

typedef interval< z_number> interval_t;
typedef linear_constraint< z_number, varname_t>  linear_constraint_t;
typedef linear_constraint_system< z_number, varname_t>  linear_constraint_system_t;
typedef linear_expression< z_number, varname_t > linear_expression_t;

// HERE TO CHOOSE THE SCALAR DOMAIN
/*
typedef interval_domain_t scalar_domain_t;
typedef octagon_domain_t scalar_domain_t;
*/
typedef dbm_domain_t scalar_domain_t;

typedef array_graph<varname_t, interval_domain_t, false, scalar_domain_t> array_graph_t;

typedef array_graph_domain<scalar_domain_t,
                           z_number,varname_t,
                           interval_domain_t, false> array_graph_domain_t;

/////////////////////////////////////////////////////////////////////////
/// TESTS FOR ARRAY_GRAPH ///////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////

// void test_0()
// {

//     cout << "TEST #0: array graph " << endl;
//     VariableFactory vfac;
//     varname_t v_A  = vfac["A"];
//     varname_t v_0  = vfac["0"];
//     varname_t v_i  = vfac["i"];
//     varname_t v_n  = vfac["n"];
      
//     array_graph_t g = array_graph_t::top();
//     g.insert_vertex(v_0);
//     g.insert_vertex(v_i);
//     g.insert_vertex(v_n);
      
//     interval_domain_t inv1 = interval_domain_t::top();
//     inv1.set (v_A, z_number (5));
      
//     interval_domain_t inv2 = interval_domain_t::bottom();
//     g.update_weight(v_0,v_i, inv1);
//     g.update_weight(v_i,v_n, inv2);
//     g.update_weight(v_n,v_i, inv2);
      
//     cout << g << endl;

// }

// void test_1(){
//   try{
//     cout << "TEST #1: array_graph reduce " << endl;
//     VariableFactory vfac;
//     varname_t v_A  = vfac["A"];
//     varname_t v_0  = vfac["0"];
//     varname_t v_i  = vfac["i"];
//     varname_t v_i1 = vfac["i+"];
//     varname_t v_n  = vfac["n"];
      
      
//     array_graph_t g = array_graph_t::top();
//     g.insert_vertex(v_0);
//     g.insert_vertex(v_i);
//     g.insert_vertex(v_i1);
//     g.insert_vertex(v_n);
      
//     uninitialized_domain_t inv1 = uninitialized_domain_t::top();
//     inv1.set(v_A, uninitialized_value_t::initialized());
      
//     uninitialized_domain_t inv2 = uninitialized_domain_t::bottom();
//     g.update_weight(v_0,v_i , inv1);
//     g.update_weight(v_i,v_i1, inv1);
//     g.update_weight(v_i1,v_i, inv2);
      
//     cout << g << endl;
      
//     scalar_domain_t inv;
//     inv += linear_constraint_t( var_t(v_0) == 0);
//     inv += linear_constraint_t( var_t(v_i) >= var_t(v_0));
//     inv += linear_constraint_t( var_t(v_n) == 11);
//     inv += linear_constraint_t( var_t(v_n) >= var_t(v_0) + 1);
//     inv += linear_constraint_t( var_t(v_i) <= var_t(v_n) - 1);
//     inv += linear_constraint_t( var_t(v_i1) == var_t(v_i) + 1);
      
//     cout << "After reduction with " << inv << endl;
//     g.reduce(inv);
//     cout << g << endl;
      
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_2(){
//   try{
//     cout << "TEST #2: array_graph forget " << endl;
//     VariableFactory vfac;
//     varname_t v_A  = vfac["A"];
//     varname_t v_0 = vfac["0"];
//     varname_t v_i = vfac["i"];
//     varname_t v_n = vfac["n"];

//     uninitialized_domain_t inv_top = uninitialized_domain_t::top();
//     uninitialized_domain_t inv_bot = uninitialized_domain_t::bottom();
//     uninitialized_domain_t inv1 = uninitialized_domain_t::top();
//     inv1.set(v_A, uninitialized_value_t::initialized());
//     uninitialized_domain_t inv2 = uninitialized_domain_t::top();
//     inv2.set(v_A, uninitialized_value_t::uninitialized());

//     array_graph_t g1(array_graph_t::top());
//     g1.insert_vertex(v_0);
//     g1.insert_vertex(v_i);
//     g1.insert_vertex(v_n);
//     g1.update_weight(v_0,v_i,inv1);
//     g1.update_weight(v_i,v_n,inv2);
//     g1.update_weight(v_0,v_n,inv_top);
//     g1.update_weight(v_i,v_0,inv_bot);
//     g1.update_weight(v_n,v_i,inv_bot);
//     g1.update_weight(v_n,v_0,inv_bot);

//     array_graph_t g2(g1);
//     cout << "g1:" << g1 << endl;
//     cout << "g2:" << g2 << endl;

//     g2 -= v_i;
//     cout << "After removing " << v_i << " from g2:" << endl;
//     cout << "g1:" << g1 << endl;
//     cout << "g2:" << g2 << endl;

//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_3(){
//   try{
//     cout << "TEST #3: array_graph <=, |, and & " << endl;
//     VariableFactory vfac;
//     varname_t v_A  = vfac["A"];
//     varname_t v_0 = vfac["0"];
//     varname_t v_i = vfac["i"];
//     varname_t v_n = vfac["n"];

//     uninitialized_domain_t inv_top = uninitialized_domain_t::top();
//     uninitialized_domain_t inv_bot = uninitialized_domain_t::bottom();
//     uninitialized_domain_t inv1 = uninitialized_domain_t::top();
//     inv1.set(v_A, uninitialized_value_t::initialized());
//     uninitialized_domain_t inv2 = uninitialized_domain_t::top();
//     inv2.set(v_A, uninitialized_value_t::uninitialized());

//     array_graph_t g1(array_graph_t::top());
//     g1.insert_vertex(v_0);
//     g1.insert_vertex(v_i);
//     g1.insert_vertex(v_n);
//     g1.update_weight(v_0,v_i,inv1);
//     g1.update_weight(v_i,v_n,inv2);
//     g1.update_weight(v_0,v_n,inv_top);
//     g1.update_weight(v_i,v_0,inv_bot);
//     g1.update_weight(v_n,v_i,inv_bot);
//     g1.update_weight(v_n,v_0,inv_bot);

//     array_graph_t g2(array_graph_t::top());
//     g2.insert_vertex(v_0);
//     g2.insert_vertex(v_i);
//     g2.insert_vertex(v_n);

//     cout << "g1:" << g1 << endl;
//     cout << "g2:" << g2 << endl;
//     bool f1  = (g1 <= g2);
//     bool f2  = (g2 <= g1);
//     cout << "(g1 <= g2): " << f1 << endl;
//     cout << "(g2 <= g1): " << f2 << endl;
//     array_graph_t g3 = g1 | g2;
//     cout << "join(g1,g2): " << g3  << endl;
//     array_graph_t g4(g1 & g2);
//     cout << "meet(g1,g2): " << g4  << endl;
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_4a(){
//   try{
//     cout << "TEST #4a: array_graph forget " << endl;
//     VariableFactory vfac;
//     varname_t v_0 = vfac["0"];
//     varname_t v_i = vfac["i"];
//     varname_t v_n = vfac["n"];
//     varname_t v_j = vfac["j"];

//     array_graph_t g1(array_graph_t::top());
//     g1.insert_vertex(v_0);
//     g1.insert_vertex(v_i);
//     g1.insert_vertex(v_n);
//     cout << "g1:" << g1 << endl;
//     g1 -= v_i;
//     cout << "g1 after removing " << v_i << ":" << g1 << endl;

//     array_graph_t g2(array_graph_t::top());
//     g2.insert_vertex(v_0);
//     g2.insert_vertex(v_i);
//     g2.insert_vertex(v_j);
//     g2.insert_vertex(v_n);
//     cout << "g2:" << g2 << endl;
//     g2 -= v_i;
//     cout << "g2 after removing " << v_i << ":" << g2 << endl;
//     g2 -= v_j;
//     cout << "g2 after removing " << v_i << " and " << v_j << ":" << g2 << endl;

//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_4b(){
//   try{
//     cout << "TEST #4b: array_graph forget " << endl;
//     VariableFactory vfac;
//     varname_t v_0 = vfac["0"];
//     varname_t v_i = vfac["i"];
//     varname_t v_n = vfac["n"];
//     varname_t v_A = vfac["A"];

//     array_graph_t g(array_graph_t::top());
//     g.insert_vertex(v_0);
//     g.insert_vertex(v_i);
//     g.insert_vertex(v_n);

//     uninitialized_domain_t w = uninitialized_domain_t::top();
//     w.set(v_A, uninitialized_value_t::initialized());
//     g.update_weight(v_0,v_n,w);

//     cout <<  g << endl;
//     g -= v_i;
//     cout << "after removing " << v_i << ":" << g << endl;

//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
// }


////////////////////////////////////////////////////////////////////////////////
/// TESTS FOR ARRAY_GRAPH_DOMAIN ///////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

void test_5()
{
    cout << "TEST #5: scalar domain x array_graph " << endl;
    VariableFactory vfac;
    array_graph_domain_t inv;
    varname_t vA = vfac["A"];
    varname_t vB = vfac["B"];
    varname_t v0 = vfac["0"];
    varname_t vi = vfac["i"];
    varname_t vj = vfac["j"];
    varname_t vn = vfac["n"];

    interval_domain_t w1, w2;
    w1.set (vA, interval_t (0,5));
    w2.set (vB, interval_t (6,10));

    inv.meet_weight(0,vi,w1);
    inv.meet_weight(0,vi,w2);
    inv.meet_weight(vi,vn,w2);
    cout << "Before adding constraints: "<< endl << inv << endl;

    linear_constraint_system_t csts;
    csts += linear_constraint_t( z_var(vi) >= 0);
    csts += linear_constraint_t( z_var(vj) >= 0);
    csts += linear_constraint_t( z_var(vi) <= z_var(vn) - 1);
    csts += linear_constraint_t( z_var(vj) <= z_var(vn) - 1);

    inv += csts;

    cout << "After adding constraints: "<< endl << inv << endl;
}

// void test_6(){
//   try{
//     cout << "TEST #6: (scalar domain x array_graph) load/store, assign, forget " << endl;
//     VariableFactory vfac;
//     varname_t vA  = vfac["A"];
//     varname_t vi  = vfac["i"];
//     varname_t vn  = vfac["n"];

//     uninitialized_domain_t w = uninitialized_domain_t::top();
//     w.set(vA, uninitialized_value_t::initialized());

//     array_graph_domain_t inv = array_graph_domain_t::top();
//     linear_constraint_system_t csts;
//     csts += linear_constraint_t( var_t(vi) >= 0);
//     csts += linear_constraint_t( var_t(vi) <= var_t(vn) - 1);
//     inv.assertion(csts,vfac);
//     inv.meet_weight(0,1,w,vfac); // ad-hoc
//     inv.meet_weight(1,vi,w,vfac); // ad-hoc
//     cout << inv << endl;

//     // v = A[i]
//     uninitialized_domain_t x = inv[vi];
//     uninitialized_value_t  v = x[vA]; 
//     cout << "array read  " << vA << "["<< vi << "]" <<  " = " << v  << endl;
//     // A[i] = v
//     inv.store(vi, w);
//     cout << "array write " <<  vA << "["<< vi << "]" << ":= \"initialized\"" << endl;
//     // v = A[i]
//     x = inv[vi];
//     v = x[vA]; 
//     cout << "array read  " << vA << "["<< vi << "]" <<  " = " << v  << endl;

//     varname_t vj  = vfac["j"];
//     inv.assign(vj, linear_expression_t(vi), vfac);
//     cout << vj << ":=" << vi << endl << inv << endl;    
//     inv -= vi;
//     cout << "after removing " << vi << endl << inv << endl;    
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_7(){
//   try{
//     cout << "TEST #7: (scalar domain x array_graph) arithmetic operations " << endl;
//     VariableFactory vfac;
//     varname_t vA  = vfac["A"];
//     varname_t vi  = vfac["i"];
//     varname_t vn  = vfac["n"];

//     array_graph_domain_t inv = array_graph_domain_t::top();
//     linear_constraint_system_t csts;

//     csts += linear_constraint_t( var_t(vi) >= 0);
//     csts += linear_constraint_t( var_t(vn) == 10);
//     csts += linear_constraint_t( var_t(vi) <= var_t(vn) - 1);
//     inv.assertion(csts, vfac);
//     uninitialized_domain_t w = uninitialized_domain_t::top();
//     w.set(vA, uninitialized_value_t::initialized());
    
//     inv.meet_weight(z_number(0),vi,w, vfac); // ad-hoc
//     inv.store(vi, w);    

//     cout << "Before operations: "<< inv << endl;
//     array_graph_domain_t tmp3(inv);
//     cout << vi << ":=" << vi << " + 0" << endl;
//     tmp3.apply(OP_ADDITION, vi, z_number(0), vfac);
//     cout << tmp3 << endl;
//     array_graph_domain_t tmp1(inv);
//     cout << vi << ":=" << vi << " + 1" << endl;
//     tmp1.apply(OP_ADDITION, vi, z_number(1), vfac);
//     cout << tmp1 << endl;
//     array_graph_domain_t tmp2(inv);
//     cout << vi << ":=" << vi << " + 2" << endl;
//     tmp2.apply(OP_ADDITION, vi, z_number(2), vfac);
//     cout << tmp2 << endl;
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
//   catch(ikos_error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_8(){
//   try{
//     cout << "TEST #8: (scalar domain x array_graph) more arithmetic operations " << endl;
//     VariableFactory vfac;
//     varname_t vA  = vfac["A"];
//     varname_t vi  = vfac["i"];
//     varname_t vj  = vfac["j"];
//     varname_t vn  = vfac["n"];

//     array_graph_domain_t inv = array_graph_domain_t::top();
//     linear_constraint_system_t csts;
//     csts += linear_constraint_t( var_t(vi) >= 0);
//     csts += linear_constraint_t( var_t(vn) == 15);
//     csts += linear_constraint_t( var_t(vi) <= 9);
//     csts += linear_constraint_t( var_t(vj) == var_t(vi)+2);
//     inv.assertion(csts, vfac);

//     uninitialized_domain_t w = uninitialized_domain_t::top();
//     w.set(vA, uninitialized_value_t::initialized());

//     inv.meet_weight(z_number(0),vi,w, vfac); // ad-hoc
//     inv.meet_weight(vi,vj, w, vfac); // ad-hoc
//     inv.store(vi, w);    

//     cout << "Before operations:" << inv << endl;
//     array_graph_domain_t tmp3(inv);
//     cout << vi << ":=" << vi << " + 0" << endl;
//     tmp3.apply(OP_ADDITION, vi, z_number(0), vfac);
//     // tmp3 -= vj;
//     cout << tmp3 << endl;
//     array_graph_domain_t tmp1(inv);
//     cout << vi << ":=" << vi << " + 1" << endl;
//     tmp1.apply(OP_ADDITION, vi, z_number(1), vfac);
//     // tmp1 -= vj;
//     cout << tmp1 << endl;
//     array_graph_domain_t tmp2(inv);
//     cout << vi << ":=" << vi << " + 2" << endl;
//     tmp2.apply(OP_ADDITION, vi, z_number(2), vfac);
//     // tmp2 -= vj;
//     cout << tmp2 << endl;
//     array_graph_domain_t tmp4(inv);
//     cout << vi << ":=" << vi << " + 4" << endl;
//     tmp4.apply(OP_ADDITION, vi, z_number(3), vfac);
//     // tmp4.forget(vj);
//     cout << tmp4 << endl;
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
//   catch(ikos_error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_9(){
//   try{
//     cout << "TEST #9: (scalar domain x array_graph) substraction " << endl;
//     VariableFactory vfac;
//     varname_t vA  = vfac["A"];
//     varname_t vi  = vfac["i"];
//     varname_t vk  = vfac["k"];
//     varname_t vn  = vfac["n"];

//     array_graph_domain_t inv = array_graph_domain_t::top();
//     linear_constraint_system_t csts;    
//     csts += linear_constraint_t( var_t(vi) >= 5);
//     csts += linear_constraint_t( var_t(vn) == 10);
//     csts += linear_constraint_t( var_t(vi) <= 9);
//     csts += linear_constraint_t( var_t(vk) == var_t(vi)-1);
//     inv.assertion(csts, vfac);
//     uninitialized_domain_t w = uninitialized_domain_t::top();
//     w.set(vA, uninitialized_value_t::initialized());

//     inv.meet_weight(vi,vn, w, vfac); // ad-hoc
//     inv.meet_weight(vk,vi, w, vfac); // ad-hoc
//     cout << "Before operations:" << inv << endl;
//     inv.apply(OP_SUBTRACTION, vi, z_number(1), vfac);
//     inv -= vk;
//     cout << vi << ":=" << vi << " - 1" << endl << inv << endl;
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
//   catch(ikos_error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_10(){
//   try{
//     cout << "TEST #10: (scalar domain x array_graph) substraction " << endl;
//     VariableFactory vfac;
//     varname_t vA  = vfac["A"];
//     varname_t vi  = vfac["i"];
//     varname_t vk  = vfac["k"];

//     array_graph_domain_t inv = array_graph_domain_t::top();
//     linear_constraint_system_t csts;    
//     csts += linear_constraint_t( var_t(vi) >= 1);
//     inv.assertion(csts, vfac);
//     uninitialized_domain_t w = uninitialized_domain_t::top();
//     w.set(vA, uninitialized_value_t::initialized());
//     inv.meet_weight(0, vi, w, vfac); // ad-hoc
//     cout << "Before operations:" << inv << endl;
//     inv.apply(OP_SUBTRACTION, vk, vi, z_number(1), vfac);
//     cout << vk << ":=" << vi << " - 1" << endl << inv << endl;
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
//   catch(ikos_error &e){ 
//     cout << e << endl; 
//   }
// }

// void test_11(){
//   try{
//     cout << "TEST #11 " << endl;
//     VariableFactory vfac;
//     varname_t vi  = vfac["i"]; varname_t vn  = vfac["n"]; varname_t vA  = vfac["A"];
//     ///////////////////////////////////////////////////////
//     // 1st
//     ///////////////////////////////////////////////////////
//     array_graph_domain_t inv1 = array_graph_domain_t::top();
//     /// assume n>0
//     inv1.assertion(linear_constraint_t( var_t(vn) >= 1), vfac);
//     /// i = 0
//     inv1.assign(vi, 0, vfac);
//     /// i < n
//     inv1.assertion(linear_constraint_t( var_t(vi) <= var_t(vn) - 1), vfac);
//     /// A[i] = 5;

//     uninitialized_domain_t w;
//     w.set(vA, uninitialized_value_t::initialized());
//     inv1.store(vi, w);
//     /// i++
//     inv1.apply(OP_ADDITION, vi, 1, vfac);    
//     cout << "after 1st iteration: " << endl << inv1 << endl;
//     ///////////////////////////////////////////////////////
//     // 2nd
//     ///////////////////////////////////////////////////////
//     array_graph_domain_t inv2(inv1);
//     /// i < n
//     inv2.assertion(linear_constraint_t( var_t(vi) <= var_t(vn) - 1), vfac);
//     /// A[i] = 5;
//     inv2.store(vi, w);
//     /// i++
//     inv2.apply(OP_ADDITION, vi, 1, vfac);    
//     cout << "after 2nd iteration: " << endl << inv2 << endl;
//     array_graph_domain_t inv3 = inv1 || inv2;
//     cout << "Widening 1st and 2nd: " << endl << inv3 << endl;
//     ///////////////////////////////////////////////////////
//     // 3rd
//     ///////////////////////////////////////////////////////
//     /// i < n
//     inv3.assertion(linear_constraint_t( var_t(vi) <= var_t(vn) -1), vfac);
//     /// A[i] = 5;
//     inv3.store(vi, w);
//     /// i++
//     inv3.apply(OP_ADDITION, vi, 1, vfac);    
//     cout << "after 3rd iteration: " << endl << inv3 << endl;
//     // ///////////////////////////////////////////////////////
//     // // 4th join with initial state
//     // ///////////////////////////////////////////////////////
//     // array_graph_domain_t inv4 = array_graph_domain_t::top();
//     // /// i = 0
//     // inv4.assign(vi, 0, vfac);
//     // inv4.assertion(linear_constraint_t( var_t(vn) >= 1), vfac);
//     // array_graph_domain_t inv5(inv3 | inv4);
//     // cout << "after join with initial state: " << endl << inv5 << endl;
//     // ///////////////////////////////////////////////////////
//     // // 5th
//     // ///////////////////////////////////////////////////////
//     // array_graph_domain_t inv6(inv5);
//     // /// i >= n
//     // inv6.assertion(linear_constraint_t( var_t(vi) >= var_t(vn)), vfac);
//     // inv6 -= vi;
//     // cout << "End: "<< inv6 << endl;
//   }
//   catch(ikos_cc::error &e){ 
//     cout << e << endl; 
//   }
//   catch(ikos_error &e){ 
//     cout << e << endl; 
//   }
// }


int main()
{
  // test_0();
  // test_1();
  // test_2();
  // test_3();
  // test_4a();
  // test_4b();
  test_5();
  // test_6();
  // test_7();
  // test_8();
  // test_9();
  // test_10();
  // test_11();
  return 42;
}
