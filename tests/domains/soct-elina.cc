#include "../common.hpp"
#include "../program_options.hpp"

#include <crab/analysis/fwd_analyzer.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<z_cfg_t> prog1(variable_factory_t &vfac)  {
  // Definining program variables
  z_var i(vfac ["i"], crab::INT_TYPE, 32);
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var x1(vfac ["x1"], crab::INT_TYPE, 32);
  z_var x2(vfac ["x2"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry","ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, bb2);
  BB(cfg, ret);
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  //  entry.assign(x1, 1);
  entry.assign(k, 0);
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(i, i, 1);
  //bb2.add(x2, x1, 1);
  bb2.add(k, k, 1);
  return cfg;
}

std::unique_ptr<z_cfg_t> prog2(variable_factory_t &vfac) {
  auto cfg = std::make_unique<z_cfg_t>("loop1_entry","ret");
  BB(cfg, loop1_entry);
  BB(cfg, loop1_bb1);
  BB(cfg, loop1_bb1_t);
  BB(cfg, loop1_bb1_f);
  BB(cfg, loop1_bb2);
  BB(cfg, loop2_entry);
  BB(cfg, loop2_bb1);
  BB(cfg, loop2_bb1_t);
  BB(cfg, loop2_bb1_f);
  BB(cfg, loop2_bb2);
  BB(cfg, ret);

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t; loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2; loop1_bb2 >> loop1_bb1; loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t; loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2; loop2_bb2 >> loop2_bb1; loop2_bb1_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var j(vfac["j"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);

  loop1_entry.assign(i, 0);
  loop1_entry.assign(k, 30);
  loop1_bb1_t.assume(i <= 9);
  loop1_bb1_f.assume(i >= 10);
  loop1_bb2.add(i, i, 1);

  loop2_entry.assign(j, 0);
  loop2_bb1_t.assume(j <= 9);
  loop2_bb1_f.assume(j >= 10);
  loop2_bb2.add(j, j, 1);
  return cfg;
}

#ifdef HAVE_ELINA
template<class Dom1, class Dom2>
void run_and_compare(crab::cfg_impl::z_cfg_t* cfg,
		     crab::cfg_impl::basic_block_label_t entry,
		     unsigned widening, 
		     unsigned narrowing, 
		     unsigned jump_set_size,
		     bool enable_stats){

  using namespace crab::analyzer;
  typedef intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom1> analyzer_1_t;
  typedef intra_fwd_analyzer<crab::cfg_impl::z_cfg_ref_t, Dom2> analyzer_2_t;  
  typename analyzer_1_t::assumption_map_t assumptions1;
  typename analyzer_2_t::assumption_map_t assumptions2;

  crab::fixpoint_parameters fixpo_params;
  fixpo_params.get_widening_delay() = widening;
  fixpo_params.get_descending_iterations() = narrowing;
  fixpo_params.get_max_thresholds() = jump_set_size;

  
  // Run with Dom1
  Dom1 absval_fac1, inv1;
  analyzer_1_t a1(*cfg, absval_fac1, nullptr, fixpo_params);
  a1.run(entry, inv1, assumptions1);
  
  // Run with Dom2
  Dom2 absval_fac2, inv2;        
  analyzer_2_t a2(*cfg, absval_fac2, nullptr, fixpo_params);
  a2.run(entry, inv2, assumptions2);

  for (auto &bb: *cfg) {
    Dom1 inv1 = a1.get_pre(bb.label());
    Dom2 inv2 = a2.get_pre(bb.label());
    
    // use elina as baseline
    z_oct_elina_domain_t baseline1, baseline2;
    
    baseline1 += inv1.to_linear_constraint_system();
    baseline2 += inv2.to_linear_constraint_system();
    bool error = false;
    if (!(baseline1 <= baseline2)) {
      crab::outs() << "FAIL. Different invariant at block " << bb << ". ";
      crab::outs() << inv1.domain_name() << " is not included in "
		   << inv2.domain_name() << ":\n";
      crab::outs() << inv1.domain_name() << ":" << inv1 << "\n";
      crab::outs() << inv2.domain_name() << ":" << inv2 << "\n";
      error = true;
    }
    if (!(baseline2 <= baseline1)) {
      crab::outs() << "FAIL. Different invariant at block " << bb << ". ";
      crab::outs() << inv2.domain_name() << " is not included in "
		   << inv1.domain_name() << ":\n";
      crab::outs() << inv1.domain_name() << ":" << inv1 << "\n";
      crab::outs() << inv2.domain_name() << ":" << inv2 << "\n";
      error = true;      
    }
    if (error) {
      crab::outs() << "ERROR: Domains " << inv1.domain_name() << " and "
		   << inv2.domain_name() << " disagree agree on block " 
		   << bb.label () << "\n";
    }    
  }
}
#endif 


std::unique_ptr<z_cfg_t> prog3(variable_factory_t &vfac) {
  auto cfg = std::make_unique<z_cfg_t>("entry","ret");
  BB(cfg, entry);
  BB(cfg, loop1_head);
  BB(cfg, loop1_t);
  BB(cfg, loop1_f);
  BB(cfg, loop1_body);

  BB(cfg, loop1_body_t);
  BB(cfg, loop1_body_f);
  BB(cfg, loop1_body_x);

  BB(cfg, cont);
  BB(cfg, loop2_head);
  BB(cfg, loop2_t);
  BB(cfg, loop2_f);
  BB(cfg, loop2_body);
  BB(cfg, ret);

  entry >> loop1_head;
  loop1_head >> loop1_t; 
  loop1_head >> loop1_f; 
  loop1_t >>    loop1_body; 

  loop1_body >> loop1_body_t;
  loop1_body >> loop1_body_f;
  loop1_body_t >> loop1_body_x;
  loop1_body_f >> loop1_body_x;
  loop1_body_x >> loop1_head;

  loop1_f >> cont;
  cont >> loop2_head;
  loop2_head >> loop2_t; 
  loop2_head >> loop2_f; 
  loop2_t >>    loop2_body; 
  loop2_body >> loop2_head;
  loop2_f >> ret;
  
  z_var i(vfac["i"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  loop1_t.assume(i <= 10);
  loop1_f.assume(i >= 11);
  loop1_body.add(i, i, 1);

  loop1_body_t.assume(i >= 9);
  loop1_body_t.assign(i , 0);
  loop1_body_f.assume(i <= 8);

  loop2_t.assume(i <= 100);
  loop2_f.assume(i >= 101);
  loop2_body.sub(i, i, 1);
  return cfg;
}

std::unique_ptr<z_cfg_t> prog4(variable_factory_t &vfac) {
  auto cfg = std::make_unique<z_cfg_t>("entry","ret");
  BB(cfg, entry);
  BB(cfg, loop_head);
  BB(cfg, loop_t);
  BB(cfg, loop_f);
  BB(cfg, loop_body);
  BB(cfg, ret);

  entry >> loop_head;
  loop_head >> loop_t; 
  loop_head >> loop_f; 
  loop_t >> loop_body; 
  loop_body >> loop_head;
  loop_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var p(vfac["p"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  entry.assign(p, 0);

  loop_t.assume(i <= 9);
  loop_f.assume(i >= 10);
  loop_body.add(i, i, 1);
  loop_body.add(p, p, 4);

  return cfg;
}

/* Example of how to build a CFG */
std::unique_ptr<z_cfg_t> prog5(variable_factory_t &vfac)  {
  // Definining program variables
  z_var i(vfac ["i"], crab::INT_TYPE, 32);
  z_var k(vfac ["k"], crab::INT_TYPE, 32);
  z_var nd(vfac ["nd"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry","ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, bb2);
  BB(cfg, ret);
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  entry.assign(k, 0);
  entry.assign(i, 0);
  bb1_t.assume(i != 9);
  bb1_f.assume(i == 9);
  bb2.add(i, i, 1);
  bb2.add(k, k, 1);
  return cfg;
}

int main(int argc, char** argv) {
#ifdef HAVE_ELINA
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc,argv,stats_enabled)) {
      return 0;
  }

  // From zones test
  {
    variable_factory_t vfac;
    auto cfg = prog1(vfac);
    // crab::outs() << *cfg << "\n";
    z_oct_elina_domain_t inv;
    run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg.get(), cfg->entry(), 1, 2, 20, stats_enabled);
  }

  {
    variable_factory_t vfac;
    auto cfg = prog2(vfac);
    // crab::outs() << *cfg << "\n";
    z_oct_elina_domain_t inv;
    run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg.get(), cfg->entry(), 1, 2, 20, stats_enabled);
  }

  {
    variable_factory_t vfac;
    auto cfg = prog3(vfac);
    // crab::outs() << *cfg << "\n";
    z_oct_elina_domain_t inv;
    run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg.get(), cfg->entry(), 1, 2, 20, stats_enabled);
  }

  {
    variable_factory_t vfac;
    auto cfg = prog4(vfac);
    // crab::outs() << *cfg << "\n";
    z_oct_elina_domain_t inv;
    run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg.get(), cfg->entry(), 1, 2, 20, stats_enabled);
  }

  {
    variable_factory_t vfac;
    auto cfg = prog5(vfac);
    // crab::outs() << *cfg << "\n";
    z_oct_elina_domain_t inv;
    run_and_compare<z_oct_elina_domain_t, z_soct_domain_t>(cfg.get(), cfg->entry(), 1, 2, 20, stats_enabled);
  }
  
  #endif
  return 0;
}
