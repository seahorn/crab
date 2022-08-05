#include "../program_options.hpp"
#include "../common.hpp"


using namespace crab::analyzer;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;
using namespace ikos;

/* Example of how to build a CFG */
z_cfg_t *prog1(variable_factory_t &vfac) {

  /*
    int i,x;
    int N = nd_int();
    __CRAB_assume(N > 0);
    i = 0;
    x = 0; 
    y = 0;
    while (i < N) {
      i++;
      x' = x + 4;
      x  = x'
      y  = y + 8;
    }
    __CRAB_assert(x == 4*N ); // OK and provable by fixed_tvpi_domain
    __CRAB_assert(y == 8*N ); // OK and provable by fixed_tvpi_domain
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);  
  z_var x_next(vfac["x.next"], crab::INT_TYPE, 32);  
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> loop_header;
  loop_header >> loop_body;
  loop_header >> loop_exit;
  loop_body >> loop_header;
  loop_exit >> exit;
  // adding statements
  entry.havoc(n);
  entry.assume(n >= 1);
  entry.assign(i, 0);
  entry.assign(x, 0);
  entry.assign(y, 0);  
  loop_body.assume(z_lin_exp_t(i) < n);
  loop_exit.assume(z_lin_exp_t(i) >= n);
  loop_exit.assertion(x == 4*n);
  loop_exit.assertion(y == 8*n);  
  loop_body.add(i, i, 1);
  loop_body.add(x_next, x, 4);
  loop_body.add(y, y, 8);  
  loop_body.assign(x, x_next);

  return cfg;
}


z_cfg_t *prog2(variable_factory_t &vfac) {

  /*
    int i,x;
    int N = nd_int();
    __CRAB_assume(N > 0);
    i = 0;
    x = 0; 
    while (i < N) {
      i++;
      if (*) {
       x = x+2;
      } else {
       x = x+3;
      }
    }
    __CRAB_assert(x >= 2*N); // OK but NOT provable by fixed_tvpi_domain
    __CRAB_assert(x <= 3*N); // OK but NOT provable by fixed_tvpi_domain
   */
  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var n(vfac["n"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_header = cfg->insert("loop_header");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &loop_body_then = cfg->insert("loop_body_then");
  z_basic_block_t &loop_body_else = cfg->insert("loop_body_else");
  z_basic_block_t &loop_body_tail = cfg->insert("loop_body_tail");      
  z_basic_block_t &loop_exit = cfg->insert("loop_exit");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> loop_header;
  loop_header >> loop_body;
  loop_header >> loop_exit;
  loop_body >> loop_body_then;
  loop_body >> loop_body_else;
  loop_body_then >> loop_body_tail;
  loop_body_else >> loop_body_tail;    
  loop_body_tail >> loop_header;
  loop_exit >> exit;
  // adding statements
  entry.havoc(n);
  entry.assume(n >= 1);
  entry.assign(i, 0);
  entry.assign(x, 0);
  loop_body.assume(z_lin_exp_t(i) < n);
  loop_exit.assume(z_lin_exp_t(i) >= n);
  loop_exit.assertion(x >= 2*n);
  loop_exit.assertion(x <= 3*n);  
  loop_body.add(i, i, 1);
  loop_body_then.add(x, x, 2);
  loop_body_else.add(x, x, 3);
  return cfg;
}


int main (int argc, char** argv) {
  bool stats_enabled = false;  
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  
  {
    crab_domain_params_man::get().coefficients().push_back(4);
    crab_domain_params_man::get().coefficients().push_back(8); 
    variable_factory_t vfac;
    z_cfg_t *cfg = prog1(vfac);
    z_fixed_tvpi_domain_t init;
    run(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    crab_domain_params_man::get().coefficients().clear();
  }

  {
    crab_domain_params_man::get().coefficients().push_back(2);
    crab_domain_params_man::get().coefficients().push_back(3); 
    variable_factory_t vfac;
    z_cfg_t *cfg = prog2(vfac);
    //crab::outs() << *cfg << "\n";
    z_fixed_tvpi_domain_t init;
    //crab_domain_params_man::get().write(crab::outs());
    run(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    run_and_check(cfg, cfg->entry(), init, false, 2, 1, 20, stats_enabled);
    crab_domain_params_man::get().coefficients().clear();
  }

  {
    crab_domain_params_man::get().coefficients().push_back(2);
    variable_factory_t vfac;
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);  
    
    z_fixed_tvpi_domain_t val1;
    
    val1.assign(y, 2*x);
    val1 += (x >= z_number(0));
    val1 += (x <= z_number(3));
    
    bool check1 = val1.entails(y >= z_number(0));
    if (check1) {
      crab::outs() << "Check1 passed\n";
      val1 += (y >= z_number(0));
    } else {
      crab::outs() << "Check1 failed\n";
    }
    bool check2 = val1.entails(y <= z_number(6));
    if (check2) {
      crab::outs() << "Check2 passed\n";      
      val1 += (y <= z_number(6));
    } else {
      crab::outs() << "Check2 failed\n";            
    }
    bool check3 = val1.entails(2*x == y);
    if (check3) {
      crab::outs() << "Check3 passed\n";      
    } else {
      crab::outs() << "Check3 failed\n";            
    }
    
    crab::outs() <<  val1 << "\n";
    crab_domain_params_man::get().coefficients().clear();
  }
  
  // {
  //   variable_factory_t vfac;
  //   z_var x(vfac["x"], crab::INT_TYPE, 32);
  //   z_var i(vfac["i"], crab::INT_TYPE, 32);
  //   z_var x_next(vfac["x.next"], crab::INT_TYPE, 32);
  //   z_var i_next(vfac["i.next"], crab::INT_TYPE, 32);
  //   z_var n(vfac["n"], crab::INT_TYPE, 32);    
  //   z_fixed_tvpi_domain_t inv;

  //   inv += (n >= z_number(1));
  //   inv.assign(i, 0);
  //   inv.assign(x, 0);
  //   inv.apply(OP_ADDITION, i, i, 1);
  //   inv.apply(OP_ADDITION, x_next, x, 4);
  //   inv.assign(x, x_next);
  //   inv += (x == 4*z_lin_exp_t(n));
  //   crab::outs() << inv << "\n";        
  // }

  return 0;
}
