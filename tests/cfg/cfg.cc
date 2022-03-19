#include "../common.hpp"
#include "../program_options.hpp"
#include <crab/analysis/graphs/cdg.hpp>
#include <crab/fixpoint/wto.hpp>

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
z_cfg_t *prog(variable_factory_t &vfac) {

  // Defining program variables
  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var k(vfac["k"], crab::INT_TYPE, 32);
  z_var nd(vfac["nd"], crab::INT_TYPE, 32);
  z_var inc(vfac["inc"], crab::INT_TYPE, 32);
  // entry and exit block
  auto cfg = new z_cfg_t("x0", "ret");
  // adding blocks
  z_basic_block_t &x0 = cfg->insert("x0");
  z_basic_block_t &x1 = cfg->insert("x1");
  z_basic_block_t &x2 = cfg->insert("x2");
  z_basic_block_t &x3 = cfg->insert("x3");
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &bb1 = cfg->insert("bb1");
  z_basic_block_t &bb1_t = cfg->insert("bb1_t");
  z_basic_block_t &bb1_f = cfg->insert("bb1_f");
  z_basic_block_t &bb2 = cfg->insert("bb2");
  z_basic_block_t &ret = cfg->insert("ret");
  // adding control flow
  x0 >> x1;
  x1 >> x2;
  x2 >> x3;
  x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t;
  bb1 >> bb1_f;
  bb1_t >> bb2;
  bb2 >> bb1;
  bb1_f >> ret;
  // adding statements
  x0.assign(k, 2147483648);
  entry.assign(i, 0);
  bb1_t.assume(i <= 99);
  bb1_f.assume(i >= 100);
  bb2.havoc(nd);
  bb2.select(inc, nd, 1, 2);
  bb2.add(i, i, inc);

  return cfg;
}

z_cfg_t *prog2(variable_factory_t &vfac) {
  // entry and exit block
  auto cfg = new z_cfg_t("entry", "exit");
  // adding blocks
  z_basic_block_t &entry = cfg->insert("entry");
  z_basic_block_t &b1 = cfg->insert("x1");
  z_basic_block_t &b2 = cfg->insert("x2");
  z_basic_block_t &b3 = cfg->insert("x3");
  z_basic_block_t &b4 = cfg->insert("x4");
  z_basic_block_t &b5 = cfg->insert("x5");
  z_basic_block_t &b6 = cfg->insert("x6");
  z_basic_block_t &b7 = cfg->insert("x7");
  z_basic_block_t &b8 = cfg->insert("x8");
  z_basic_block_t &b9 = cfg->insert("x9");
  z_basic_block_t &exit = cfg->insert("exit");
  // adding control flow
  entry >> b1;
  b1 >> b2;
  b2 >> b3;
  b2 >> b4;
  b3 >> b5;
  b4 >> b5;
  b5 >> b1;
  b1 >> b6;
  b6 >> b7;
  b6 >> b8;
  b7 >> b9, b8 >> b9;
  b9 >> exit;
  // if added this edge then b6 should control dependent on b1.
  b1 >> exit;
  return cfg;
}

int main(int argc, char **argv) {
  bool stats_enabled = false;
  if (!crab_tests::parse_user_options(argc, argv, stats_enabled)) {
    return 0;
  }
  variable_factory_t vfac;

  {
    z_cfg_t *cfg = prog(vfac);

    crab::outs() << "CFG\n" << *cfg << "\n";
    z_cfg_rev_t rev_cfg(*cfg);
    crab::outs() << "Reversed CFG\n" << rev_cfg << "\n";

    crab::outs() << "Weak reversed topological order of CFG \n";
    bool first = true;
    for (auto &N :
         crab::analyzer::graph_algo::weak_rev_topo_sort(z_cfg_ref_t(*cfg))) {
      if (!first)
        crab::outs() << " -- ";
      crab::outs() << N;
      first = false;
    }
    crab::outs() << "\n";

    crab::outs() << "Weak topological order of the reversed CFG \n";
    first = true;
    for (auto &N : crab::analyzer::graph_algo::weak_topo_sort(rev_cfg)) {
      if (!first)
        crab::outs() << " -- ";
      crab::outs() << N;
      first = false;
    }
    crab::outs() << "\n";

    crab::outs() << "Bourdoncle WTO of the reversed CFG\n";
    ikos::wto<z_cfg_rev_t> wto_g(rev_cfg);
    crab::outs() << wto_g << "\n";

    cfg->simplify();
    crab::outs() << "Simplified CFG\n" << *cfg << "\n";

    rev_cfg = z_cfg_rev_t(*cfg);
    crab::outs() << "Reversed simplified CFG\n" << rev_cfg << "\n";

    delete cfg;
  }

  {

    z_cfg_t *cfg = prog2(vfac);
    crab::outs() << "CFG\n" << *cfg << "\n";

    std::map<typename z_cfg_t::node_t, std::vector<typename z_cfg_t::node_t>>
        cdg;
    crab::analyzer::graph_algo::control_dep_graph(z_cfg_ref_t(*cfg), cdg);
    crab::outs() << "Control-dependence graph \n";
    for (auto &kv : cdg) {
      crab::outs() << "{";
      for (auto v : kv.second) {
        crab::outs() << crab::basic_block_traits<z_basic_block_t>::to_string(v)
		     << ";";
      }
      crab::outs() << "} "
                   << " control-dependent on ";
      crab::outs() << crab::basic_block_traits<z_basic_block_t>::to_string(kv.first)
		   << "\n";
    }

    delete cfg;
  }
  return 0;
}
