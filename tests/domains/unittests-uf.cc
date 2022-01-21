#include "../program_options.hpp"
#include "../common.hpp"
#include <crab/domains/uf_domain.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

namespace {
  using z_uf_domain_t = uf_domain<z_number, varname_t>;
}

int main(int argc, char** argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }
  variable_factory_t vfac;
  
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var w(vfac["w"], crab::INT_TYPE, 32);
  z_var z(vfac["z"], crab::INT_TYPE, 32);

  crab::outs() << "=== Join === \n";  
  {
    z_uf_domain_t left;
    z_uf_domain_t right;
    left.assign(y, z_number(8));
    left.apply(OP_MULTIPLICATION, x, y, z_number(5));
    right.apply(OP_MULTIPLICATION, x, y, z_number(5));
    z_uf_domain_t l_join_r = left | right;
    crab::outs() << left << " | " << right << " = " << l_join_r << "\n";
  }
  
  {
    z_uf_domain_t left;
    left.assign(x, z_number(1));
    z_uf_domain_t right = left;
    right.apply(OP_ADDITION, x, x, z_number(1));
    z_uf_domain_t l_join_r = left | right;
    crab::outs() << left << " | " << right << " = " << l_join_r << "\n";
  }

  crab::outs() << "=== Equalities/Disequalities === \n";
  {
    z_var x1(vfac["x1"], crab::INT_TYPE, 32);
    z_var x2(vfac["x2"], crab::INT_TYPE, 32);
    z_var x3(vfac["x3"], crab::INT_TYPE, 32);

    z_uf_domain_t dom;
    dom.apply(OP_MULTIPLICATION, x1, y, z_number(5));
    dom.apply(OP_MULTIPLICATION, x2, y, z_number(5));
    dom.apply(OP_MULTIPLICATION, x3, y, z_number(8));
    
    crab::outs() << "Extracting equalities from " << dom << "="
		 << dom.to_linear_constraint_system() << "\n";

    z_lin_cst_t c1(z_lin_exp_t(x1) == z_lin_exp_t(x2));    
    dom += c1;
    crab::outs() << "After adding " << c1 << " --> " << dom << "\n";
    z_lin_cst_t c2(z_lin_exp_t(x1) != z_lin_exp_t(x2));    
    dom += c2;
    crab::outs() << "After adding " << c2 << " --> " << dom << "\n";
  }

  crab::outs() << "=== Project === \n";
  {
    z_uf_domain_t inv;
    inv.assign(x, 5);
    inv.assign(y, 5);
    inv.assign(z, 9);
    std::vector<z_var> vs;
    vs.push_back(x);
    vs.push_back(y);
    crab::outs() << "Before project " << inv << "\n";
    inv.project(vs);
    crab::outs() << "After project " << inv << "\n";
  }

  crab::outs() << "==== Meet ====\n";
  {
    z_uf_domain_t left, right;

    // {w = 5, x = 5, y = '+'(5,3), z = 3}
    left.assign(x, 5);
    left.assign(w, x);
    left.assign(z, 3);
    left.apply(OP_ADDITION, y, x, z);
    
    // {w = 8,  x = '+'(8,2), y = 8, z = 2}
    right.assign(y, 8); 
    right.assign(w, y);
    right.assign(z, 2); 
    right.apply(OP_ADDITION, x, w, z);

    // meet = {x=y=w= '+'(VAR_0,VAR_1), z=VAR_1}
    crab::outs() << "Meet" << "\n  " << left << " \n  " << right << "\n"; 
    crab::outs() << "Result=" << (left & right) << "\n";
  }

  {
    z_uf_domain_t left, right;

    // {w = 5, x = 5, y = '+' (5,10), z = 10} --> w=x and y=+(x,z) and z=10
    left.assign(x, 5);
    left.assign(w, x);
    left.assign(z, 10);
    left.apply(OP_ADDITION, y, x, z);
    
    // {w = 7, y = 7, x = '+' (7,5), z = 5} --> w=y and x=+(w,z) and z=5
    right.assign(y, 7); 
    right.assign(w, y);
    right.assign(z, 10); 
    right.apply(OP_ADDITION, x, w, z);

    // meet = ({x=y=w= '+'(VAR_0,VAR_1), z=VAR_1}
    crab::outs() << "Meet" << "\n  " << left << " \n  " << right << "\n"; 
    crab::outs() << "Result=" << (left & right) << "\n";
  }

  // Nothing wrong with this but different printing depending on OS.
  #if 0
  crab::outs() << "=== Arrays ==== \n";
  {
    z_uf_domain_t dom;
    z_var a(vfac["A"], crab::ARR_INT_TYPE, 32);    
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var w(vfac["w"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);        

    // Create the term w = f(a, y)
    dom.array_load(w, a, 4, y);
    // Create the term z = f(a, y)    
    dom.array_load(z, a, 4, y);
    // At this point, we know w=z
    crab::outs() << dom << "\n";
    
    z_lin_cst_t c1(z_lin_exp_t(w) == z_lin_exp_t(z));
    crab::outs() << "Added " << c1 << "\n";
    dom += c1;
    crab::outs() << "Result=" << dom << "\n";
    
    z_lin_cst_t c2(z_lin_exp_t(w) != z_lin_exp_t(z));    
    crab::outs() << "Added " << c2 << "\n";
    dom += c2;
    crab::outs() << "Result=" << dom << "\n";    
  }
  #endif
  
  crab::outs() << "==== Boolean operations ====\n";
  {
    z_var b1(vfac["b1"], crab::BOOL_TYPE, 1);
    z_var b2(vfac["b2"], crab::BOOL_TYPE, 1);
    z_var b3(vfac["b3"], crab::BOOL_TYPE, 1);
    z_var b4(vfac["b4"], crab::BOOL_TYPE, 1);
    z_var b5(vfac["b5"], crab::BOOL_TYPE, 1);
    z_var b6(vfac["b6"], crab::BOOL_TYPE, 1);
    z_var b7(vfac["b7"], crab::BOOL_TYPE, 1);            

    z_uf_domain_t dom;
    dom.assign_bool_var(b4, b5, false);
    dom.apply_binary_bool(OP_BAND, b2, b1, b4);
    dom.apply_binary_bool(OP_BAND, b3, b1, b5);
    dom.apply_binary_bool(OP_BOR, b6, b1, b5);
    dom.apply_binary_bool(OP_BXOR, b7, b1, b5);        
    crab::outs() << "Extracting equalities from " << dom << "="
		 << dom.to_linear_constraint_system() << "\n";
  }
  
  return 0;
}
