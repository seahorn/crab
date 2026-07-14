#include "../common.hpp"
#include "../program_options.hpp"

// To run abstract domains defined over reals

using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
std::unique_ptr<q_cfg_t> prog(variable_factory_t &vfac) {

  // Definining program variables
  q_var x(vfac["x"], crab::REAL_TYPE);
  q_var y(vfac["y"], crab::REAL_TYPE);
  // entry and exit block
  auto cfg = std::make_unique<q_cfg_t>("entry", "exit");
  // adding blocks
  BB(cfg, entry);
  BB(cfg, header);
  BB(cfg, body);
  BB(cfg, exit);

  // adding control flow
  entry >> header;
  header >> body;
  body >> header;
  header >> exit;

  // adding statements
  entry.assign(x, q_number(1));
  entry.assign(y, q_number(0));

  body.add(x, x, y);
  body.add(y, y, q_number(1));

  exit.assertion(x >= y);
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
