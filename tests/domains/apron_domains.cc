#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog1(variable_factory_t &vfac) {

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var x1(vfac["x1"], crab::INT_TYPE, 32);
  z_var x2(vfac["x2"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  //  entry.assign(x1, 1);
  entry.assign(k, 0);
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(i, i, 1);
  // bb2.add(x2, x1, 1);
  bb2.add(k, k, 1);
  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac) {

  z_cfg_t *cfg = new z_cfg_t("loop1_entry", "ret");
  z_basic_block_t &loop1_entry = cfg->insert("loop1_entry");
  z_basic_block_t &loop1_bb1 = cfg->insert("loop1_bb1");
  z_basic_block_t &loop1_bb1_t = cfg->insert("loop1_bb1_t");
  z_basic_block_t &loop1_bb1_f = cfg->insert("loop1_bb1_f");
  z_basic_block_t &loop1_bb2 = cfg->insert("loop1_bb2");
  z_basic_block_t &loop2_entry = cfg->insert("loop2_entry");
  z_basic_block_t &loop2_bb1 = cfg->insert("loop2_bb1");
  z_basic_block_t &loop2_bb1_t = cfg->insert("loop2_bb1_t");
  z_basic_block_t &loop2_bb1_f = cfg->insert("loop2_bb1_f");
  z_basic_block_t &loop2_bb2 = cfg->insert("loop2_bb2");
  z_basic_block_t &ret = cfg->insert("ret");

  loop1_entry >> loop1_bb1;
  loop1_bb1 >> loop1_bb1_t;
  loop1_bb1 >> loop1_bb1_f;
  loop1_bb1_t >> loop1_bb2;
  loop1_bb2 >> loop1_bb1;
  loop1_bb1_f >> loop2_entry;

  loop2_entry >> loop2_bb1;
  loop2_bb1 >> loop2_bb1_t;
  loop2_bb1 >> loop2_bb1_f;
  loop2_bb1_t >> loop2_bb2;
  loop2_bb2 >> loop2_bb1;
  loop2_bb1_f >> ret;

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

z_cfg_t *prog3(variable_factory_t &vfac) {

  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop1_head = cfg->insert("loop1_head");
  z_basic_block_t &loop1_t = cfg->insert("loop1_t");
  z_basic_block_t &loop1_f = cfg->insert("loop1_f");
  z_basic_block_t &loop1_body = cfg->insert("loop1_body");

  z_basic_block_t &loop1_body_t = cfg->insert("loop1_body_t");
  z_basic_block_t &loop1_body_f = cfg->insert("loop1_body_f");
  z_basic_block_t &loop1_body_x = cfg->insert("loop1_body_x");

  z_basic_block_t &cont = cfg->insert("cont");
  z_basic_block_t &loop2_head = cfg->insert("loop2_head");
  z_basic_block_t &loop2_t = cfg->insert("loop2_t");
  z_basic_block_t &loop2_f = cfg->insert("loop2_f");
  z_basic_block_t &loop2_body = cfg->insert("loop2_body");
  z_basic_block_t &ret = cfg->insert("ret");

  entry >> loop1_head;
  loop1_head >> loop1_t;
  loop1_head >> loop1_f;
  loop1_t >> loop1_body;

  loop1_body >> loop1_body_t;
  loop1_body >> loop1_body_f;
  loop1_body_t >> loop1_body_x;
  loop1_body_f >> loop1_body_x;
  loop1_body_x >> loop1_head;

  loop1_f >> cont;
  cont >> loop2_head;
  loop2_head >> loop2_t;
  loop2_head >> loop2_f;
  loop2_t >> loop2_body;
  loop2_body >> loop2_head;
  loop2_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  loop1_t.assume(i <= 10);
  loop1_f.assume(i >= 11);
  loop1_body.add(i, i, 1);

  loop1_body_t.assume(i >= 9);
  loop1_body_t.assign(i, 0);
  loop1_body_f.assume(i <= 8);

  loop2_t.assume(i <= 100);
  loop2_f.assume(i >= 101);
  loop2_body.sub(i, i, 1);
  return cfg;
}

z_cfg_t *prog4(variable_factory_t &vfac) {

  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &loop_head = cfg->insert("loop_head");
  z_basic_block_t &loop_t = cfg->insert("loop_t");
  z_basic_block_t &loop_f = cfg->insert("loop_f");
  z_basic_block_t &loop_body = cfg->insert("loop_body");
  z_basic_block_t &ret = cfg->insert("ret");

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
z_cfg_t *prog5(variable_factory_t &vfac) {

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "ret");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  entry.assign(k, 0);
  entry.assign(i, 0);
  bb1_t.assume(i != 9);
  bb1_f.assume(i == 9);
  bb2.add(i, i, 1);
  bb2.add(k, k, 1);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
#ifdef HAVE_APRON

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

#if 1

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog1(vfac);
    crab::outs() << *cfg << "\n";
    {
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_box_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_oct_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_pk_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog2(vfac);
    crab::outs() << *cfg << "\n";
    {
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_box_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_oct_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_pk_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog3(vfac);
    crab::outs() << *cfg << "\n";
    {
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_box_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_oct_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_pk_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog4(vfac);
    crab::outs() << *cfg << "\n";
    {
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_box_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_oct_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_pk_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    delete cfg;
  }

  {
    variable_factory_t vfac;
    z_cfg_t *cfg = prog5(vfac);
    crab::outs() << *cfg << "\n";
    {
      z_interval_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_box_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_oct_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    {
      z_pk_apron_domain_t init;
      run(cfg, cfg->entry(), init, false, 1, 2, 20, stats_enabled);
    }
    delete cfg;
  }
#endif

  /////
  // testing operations
  /////

#if 0
  {
    variable_factory_t vfac;
    
    ap_manager_t* man = box_manager_alloc();;
    // x:0 y:1 z:2
    ap_state_ptr ap1 = apPtr(man, ap_abstract0_top(man, 3, 0));
    ap1 = apPtr(man,
                 ap_abstract0_assign_texpr(man, false, 
                                           &*ap1, 
                                           0, ap_texpr0_cst_scalar_int((int) 5), 
                                           NULL));
    ap1 = apPtr(man,
                 ap_abstract0_assign_texpr(man, false, 
                                           &*ap1, 
                                           1, ap_texpr0_cst_scalar_int((int) 2), 
                                           NULL));
    
    ap_abstract0_fprint(stdout, man, &*ap1, NULL);
    
    ap_dimperm_t*  p = ap_dimperm_alloc(3);
    
    p->dim[0] = 2;
    p->dim[1] = 1;
    p->dim[2] = 0;
    
    ap_dimperm_fprint(stdout, p);
    ap_state_ptr ap2 = apPtr(man, ap_abstract0_permute_dimensions(man, false, &*ap1, p));
    ap_abstract0_fprint(stdout, man, &*ap2, NULL);
    ap_dimperm_free(p);
  }
  
  { 
    variable_factory_t vfac;
    z_pk_apron_domain_t inv1;
    inv1.assign(vfac ["x"], 5);
    z_lin_cst_sys_t csts;
    csts +=(z_lin_exp_t(vfac ["x"]) == z_lin_exp_t(vfac ["y"]));
    inv1 += csts;
    z_pk_apron_domain_t inv2(inv1);
    crab::outs() << "Before expand x into z:" << inv1 << "\n";
    inv1.expand(vfac ["x"], vfac["z"]);
    crab::outs() << "After expand x into z: " << inv1 << "\n";
    crab::outs() << "Copy before: " << inv2 << "\n";

    pk_apron_domain_t inv3 = inv1 | inv2;
    crab::outs() << "Join: " << inv3 << "\n";
  }
  
  { 
    variable_factory_t vfac;
    z_pk_apron_domain_t inv1;
    inv1.assign(vfac ["x"], 5);

    z_pk_apron_domain_t inv2(inv1);
    inv2.apply(OP_ADDITION, vfac ["x"], vfac ["x"], 1);

    z_pk_apron_domain_t inv3 = inv1;
    crab::outs() << inv1 << "\n";
    crab::outs() << inv2 << "\n";
    crab::outs() << inv3 << "\n";
  }

  {
    variable_factory_t vfac;
    z_oct_apron_domain_t inv1;
    z_oct_apron_domain_t inv2;
    z_var i(vfac ["i"], crab::INT_TYPE, 32);
    z_var k(vfac ["k"], crab::INT_TYPE, 32);

    {
      z_lin_cst_sys_t csts;
      csts += (z_lin_exp_t(k) <= 100);
      csts += (- z_lin_exp_t(i) + z_lin_exp_t(k) <= 0);
      csts += (z_lin_exp_t(i) + z_lin_exp_t(k)  <= 200);
      csts += (-z_lin_exp_t(k) <= 0);
      csts += (-z_lin_exp_t(i) - z_lin_exp_t(k)  <= 0);
      csts += (z_lin_exp_t(i) - z_lin_exp_t(k)  <= 0);
      csts += (z_lin_exp_t(i) <= 100);
      csts += (-z_lin_exp_t(i) <= 0);
      inv1 += csts;
    }

    {
      z_lin_cst_sys_t csts;
      csts += (-z_lin_exp_t(i) + z_lin_exp_t(k)  <= 0);
      csts += (-z_lin_exp_t(k) <= 0);
      csts += (-z_lin_exp_t(i) - z_lin_exp_t(k)  <= 0);
      csts += (z_lin_exp_t(i) - z_lin_exp_t(k)  <= 0);
      csts += (-z_lin_exp_t(i) <= 0);
      inv2 += csts;
    }

    bool res = inv1 <= inv2;
    crab::outs() << "Checking " << inv1 << " <= " << inv2 << "\nRes=" << res << "\n"; 
  }
#endif
#endif

  return 0;
}
