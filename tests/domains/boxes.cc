#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog1(variable_factory_t &vfac) {

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
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.add(i, i, 1);
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

z_cfg_t *prog6(variable_factory_t &vfac) {

  // Definining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  // entry and exit block
  z_cfg_t *cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &header = cfg->insert("header");
  z_basic_block_t &body = cfg->insert("body");
  z_basic_block_t &exit = cfg->insert("exit");

  // adding control flow
  entry >> header;
  header >> body;
  body >> header;
  header >> exit;

  // adding statements
  entry.assign(x, z_number(1));
  entry.assign(y, z_number(0));

  body.add(x, x, y);
  body.add(y, y, z_number(1));

  exit.assertion(x >= y);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
#ifdef HAVE_LDD

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;

  z_cfg_t *cfg1 = prog1(vfac);
  z_cfg_t *cfg2 = prog2(vfac);
  z_cfg_t *cfg3 = prog3(vfac);
  z_cfg_t *cfg4 = prog4(vfac);
  z_cfg_t *cfg5 = prog5(vfac);
  z_cfg_t *cfg6 = prog6(vfac);

  z_boxes_domain_t init;

  cfg1->simplify(); // this is optional
  crab::outs() << *cfg1 << "\n";
  run(cfg1, cfg1->entry(), init, false, 10, 2, 20, stats_enabled);

  cfg2->simplify(); // this is optional
  crab::outs() << *cfg2 << "\n";
  run(cfg2, cfg2->entry(), init, false, 10, 2, 20, stats_enabled);

  cfg3->simplify(); // this is optional
  crab::outs() << *cfg3 << "\n";
  run(cfg3, cfg3->entry(), init, false, 10, 2, 20, stats_enabled);

  cfg4->simplify(); // this is optional
  crab::outs() << *cfg4 << "\n";
  run(cfg4, cfg4->entry(), init, false, 10, 2, 20, stats_enabled);

  crab::outs() << *cfg5 << "\n";
  run(cfg5, cfg5->entry(), init, false, 10, 2, 20, stats_enabled);

  crab::outs() << *cfg6 << "\n";
  run(cfg6, cfg6->entry(), init, false, 1, 2, 20, stats_enabled);

  delete cfg1;
  delete cfg2;
  delete cfg3;
  delete cfg4;
  delete cfg5;
  delete cfg6;

  if (true) {
    crab::outs() << "Testing some boxes operations ...\n";
    z_var x(vfac["x"], crab::INT_TYPE, 32);
    z_var y(vfac["y"], crab::INT_TYPE, 32);
    z_var z(vfac["z"], crab::INT_TYPE, 32);

    z_boxes_domain_t inv1;

    inv1.assign(y, 6);
    inv1.assign(z, 7);

    z_boxes_domain_t inv2;
    inv2.assign(y, 3);
    inv2.assign(z, 4);

    z_boxes_domain_t inv3 = inv1 | inv2;

    crab::outs() << inv3 << "\n";

    crab::outs() << x << ":=" << y << " + " << z << "= \n";
    inv3.apply(OP_ADDITION, x, y, z);
    crab::outs() << inv3 << "\n";

    z_boxes_domain_t inv4(inv3); // for later
    z_boxes_domain_t inv5(inv3); // for later

    crab::outs() << x << ":=" << y << " - " << z << "= \n";
    inv3.apply(OP_SUBTRACTION, x, y, z);
    crab::outs() << inv3 << "\n";

    crab::outs() << x << ":=" << y << " * " << z << "= \n";
    inv3.apply(OP_MULTIPLICATION, x, y, z);
    crab::outs() << inv3 << "\n";

    crab::outs() << x << ":=" << y << " / " << z << "= \n";
    inv3.apply(OP_SDIV, x, y, z);
    crab::outs() << inv3 << "\n";

    z_var cx(vfac["x"], crab::INT_TYPE, 32);
    z_var cy(vfac["y"], crab::INT_TYPE, 32);
    z_var cz(vfac["z"], crab::INT_TYPE, 32);

    crab::outs() << "INV: " << inv4 << "\n";
    inv4 += (cx >= cy);
    crab::outs() << "ADDED x >= y \n" << inv4 << "\n";
    crab::outs() << "INV: " << inv5 << "\n";
    inv5 += (cx <= cy + cz - 1);
    crab::outs() << "ADDED x <= y + z -1\n" << inv5 << "\n";

    z_boxes_domain_t inv6;
    crab::outs() << "INV: " << inv6 << "\n";
    inv6 += (cx != 9);
    crab::outs() << "ADDED x != 9\n" << inv6 << "\n";
    inv6 += (cy >= 9);
    crab::outs() << "ADDED y >= 9\n" << inv6 << "\n";
    inv6 += (cy <= 9);
    crab::outs() << "ADDED y <= 9\n" << inv6 << "\n";
    inv6 += (cz >= 10);
    crab::outs() << "ADDED z > 9\n" << inv6 << "\n";
    inv6 += (cz <= 8);
    crab::outs() << "ADDED z < 9\n" << inv6 << "\n";
  }
#endif

  return 0;
}
