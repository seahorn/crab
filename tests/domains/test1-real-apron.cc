#include "../common.hpp"
#include "../program_options.hpp"

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<q_cfg_t> prog(variable_factory_t &vfac) {

  // Definining program variables
  q_var i(vfac["i"], crab::REAL_TYPE);
  // entry and exit block
  auto cfg = std::make_unique<q_cfg_t>("entry", "ret");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, bb1);
  BB(cfg, bb1_t);
  BB(cfg, bb1_f);
  BB(cfg, ret);
  // adding control flow
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb1;
  bb1_f >> ret;
  // adding statements
  entry.assign(i, q_number(0.0));
  bb1_t.assume(i <= q_number(9.9));
  bb1_f.assume(i >= q_number(10));
  bb1_t.add(i, i, q_number(1));

  return cfg;
}

/* Example of how to infer invariants from the above CFG */
int main(int argc, char **argv) {
#ifdef HAVE_APRON  
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;
  auto cfg = prog(vfac);
  crab::outs() << *cfg << "\n";

  q_pk_apron_domain_t init;
  run(cfg, init, stats_enabled);
#endif
  return 0;
}
