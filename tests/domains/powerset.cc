#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

std::unique_ptr<z_cfg_t> prog1(variable_factory_t &vfac) {

  /*
    x = 0;
    y = 1;
    if (nd()) {
      x += 1;
      y += 2;
    }
    if (nd()) {
      x += 2;
      y += 3;
    }
    if (nd()) {
      y += 1;
    }
    assert(y > x);
   */

  // Defining program variables
  z_var x(vfac["x"], crab::INT_TYPE, 32);
  z_var y(vfac["y"], crab::INT_TYPE, 32);
  z_var nd1(vfac["nd1"], crab::INT_TYPE, 32);
  z_var nd2(vfac["nd2"], crab::INT_TYPE, 32);
  z_var nd3(vfac["nd3"], crab::INT_TYPE, 32);

  // entry and exit block
  auto cfg = std::make_unique<z_cfg_t>("entry", "exit");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, b1);
  BB(cfg, b1_tt);
  BB(cfg, b1_ff);
  BB(cfg, b2);
  BB(cfg, b2_tt);
  BB(cfg, b2_ff);
  BB(cfg, b3);
  BB(cfg, b3_tt);
  BB(cfg, b3_ff);
  BB(cfg, exit);

  // adding control flow
  entry >> b1;
  b1 >> b1_tt;
  b1 >> b1_ff;
  b1_tt >> b2;
  b1_ff >> b2;

  b2 >> b2_tt;
  b2 >> b2_ff;
  b2_tt >> b3;
  b2_ff >> b3;

  b3 >> b3_tt;
  b3 >> b3_ff;
  b3_tt >> exit;
  b3_ff >> exit;

  // adding statements
  entry.assign(x, 0);
  entry.assign(y, 1);

  b1_tt.assume(nd1 >= 1);
  b1_tt.add(x, x, 1);
  b1_tt.add(y, y, 2);
  b1_ff.assume(nd1 <= 0);

  b2_tt.assume(nd2 >= 1);
  b2_tt.add(x, x, 2);
  b2_tt.add(y, y, 3);
  b2_ff.assume(nd2 <= 0);

  b3_tt.assume(nd3 >= 1);
  b3_tt.add(y, y, 1);
  b3_ff.assume(nd3 <= 0);

  exit.assertion(y >= x + 1);
  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {

  array_adaptive_domain_params p(true/*is_smashable*/,
				 false/*smash_at_nonzero_offset*/,
				 64/*max_smashable_cells*/,
				 512/*max_array_size*/);
  crab_domain_params_man::get().update_params(p);

  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }

  variable_factory_t vfac;
  auto cfg = prog1(vfac);
  crab::outs() << *cfg << "\n";

  check_all<z_pow_aa_int_t, z_interval_domain_t>(cfg, stats_enabled);

  // free the CFG

  return 0;
}
