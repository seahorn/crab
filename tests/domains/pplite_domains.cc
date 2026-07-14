#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<z_cfg_t> prog1(variable_factory_t &vfac) {

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var x1(vfac["x1"], crab::INT_TYPE, 32);
  z_var x2(vfac["x2"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, bb2);
  BB(cfg, ret);
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

std::unique_ptr<z_cfg_t> prog2(variable_factory_t &vfac) {

  auto cfg = std::make_unique<z_cfg_t>("loop1_entry", "ret");
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

std::unique_ptr<z_cfg_t> prog3(variable_factory_t &vfac) {

  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
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

std::unique_ptr<z_cfg_t> prog4(variable_factory_t &vfac) {

  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
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
std::unique_ptr<z_cfg_t> prog5(variable_factory_t &vfac) {

  // Definining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, bb2);
  BB(cfg, ret);
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

#ifdef HAVE_PPLITE

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  {
    variable_factory_t vfac;
    auto cfg = prog1(vfac);
    crab::outs() << *cfg << "\n";
    run_all<z_pk_apron_domain_t, z_poly_pplite_domain_t, z_fpoly_pplite_domain_t, z_pset_pplite_domain_t>(cfg, stats_enabled);
#ifdef HAVE_APRON
    {
      z_decoupled_box_poly_domain_t init;
      run(cfg, init, stats_enabled);
    }
#endif     
  }

  {
    variable_factory_t vfac;
    auto cfg = prog2(vfac);
    crab::outs() << *cfg << "\n";
    run_all<z_pk_apron_domain_t, z_poly_pplite_domain_t, z_fpoly_pplite_domain_t, z_pset_pplite_domain_t>(cfg, stats_enabled);
  }

  {
    variable_factory_t vfac;
    auto cfg = prog3(vfac);
    crab::outs() << *cfg << "\n";
    run_all<z_pk_apron_domain_t, z_poly_pplite_domain_t, z_fpoly_pplite_domain_t, z_pset_pplite_domain_t>(cfg, stats_enabled);
  }

  {
    variable_factory_t vfac;
    auto cfg = prog4(vfac);
    crab::outs() << *cfg << "\n";
    run_all<z_pk_apron_domain_t, z_poly_pplite_domain_t, z_fpoly_pplite_domain_t, z_pset_pplite_domain_t>(cfg, stats_enabled);
  }

  {
    variable_factory_t vfac;
    auto cfg = prog5(vfac);
    crab::outs() << *cfg << "\n";
    run_all<z_poly_pplite_domain_t, z_poly_pplite_domain_t, z_fpoly_pplite_domain_t, z_pset_pplite_domain_t>(cfg, stats_enabled);
  }

  /////
  // testing operations
  /////

  auto make_var
    = [](variable_factory_t& vfac, const char* vname) {
        z_var res(vfac[vname], crab::INT_TYPE, 32);
        return res;
      };

  {
    variable_factory_t vfac;
    z_poly_pplite_domain_t inv1;
    auto x = make_var(vfac, "x");
    auto y = make_var(vfac, "y");
    inv1.assign(x, 5);
    z_lin_cst_sys_t csts;
    csts += (z_lin_exp_t(x) == z_lin_exp_t(y));
    inv1 += csts;
    z_poly_pplite_domain_t inv2(inv1);
    crab::outs() << "Before expand x into z:" << inv1 << "\n";
    auto z = make_var(vfac, "z");
    inv1.expand(x, z);
    crab::outs() << "After expand x into z: " << inv1 << "\n";
    crab::outs() << "Copy before: " << inv2 << "\n";

    z_poly_pplite_domain_t inv3 = inv1 | inv2;
    crab::outs() << "Join: " << inv3 << "\n";
  }

  {
    variable_factory_t vfac;
    z_poly_pplite_domain_t inv1;
    auto x = make_var(vfac, "x");
    inv1.assign(x, 5);

    z_poly_pplite_domain_t inv2(inv1);
    inv2.apply(OP_ADDITION, x, x, 1);

    z_poly_pplite_domain_t inv3 = inv1;
    crab::outs() << inv1 << "\n";
    crab::outs() << inv2 << "\n";
    crab::outs() << inv3 << "\n";
  }

#endif // HAVE_PPLITE

  return 0;
}
