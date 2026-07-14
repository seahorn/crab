#include "../common.hpp"
#include "../program_options.hpp"

using namespace std;
using namespace crab::cfg;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

// To test the interface to BGL graphs
void write(z_cfg_ref_t g) {
  crab::outs() << "Num of vertices: " << num_vertices(g) << "\n";
  for (auto v : boost::make_iterator_range(vertices(g))) {
    crab::outs() << "Vertex: " << v << "\n";
    crab::outs() << "Num of predecessors=" << in_degree(v, g) << "\n";
    crab::outs() << "Num of successors  =" << out_degree(v, g) << "\n";
    crab::outs() << "Num of neighbors   =" << degree(v, g) << "\n";
    crab::outs() << "Succs={";
    {
      auto p = out_edges(v, g);
      auto succIt = p.first;
      auto succEnd = p.second;
      for (; succIt != succEnd; ++succIt) {
        crab::outs() << target(*succIt, g) << ";";
      }
    }
    crab::outs() << "}"
                 << "\n";
    crab::outs() << "Preds={ ";
    {
      auto p = in_edges(v, g);
      auto predIt = p.first;
      auto predEnd = p.second;
      for (; predIt != predEnd; ++predIt) {
        crab::outs() << source(*predIt, g) << ";";
      }
    }
    crab::outs() << "}"
                 << "\n";
  }
}

std::unique_ptr<z_cfg_t> prog(variable_factory_t &vfac) {

  /*
     i:=0;
     p:=0;
     while (i <= 9) {
        i := i + 1;
        p := p + 4;
     }
   */
  auto cfg = std::make_unique<z_cfg_t>("entry", "ret");
  BB(cfg, entry);
  BB(cfg, loop_head);
  BB(cfg, loop_t);
  BB(cfg, loop_f);
  BB(cfg, ret);

  entry >> loop_head;
  loop_head >> loop_t;
  loop_head >> loop_f;
  loop_t >> loop_head;
  loop_f >> ret;

  z_var i(vfac["i"], crab::INT_TYPE, 32);
  z_var p(vfac["p"], crab::INT_TYPE, 32);

  entry.assign(i, 0);
  entry.assign(p, 0);

  loop_t.assume(i <= 9);
  loop_f.assume(i >= 10);
  loop_t.add(i, i, 1);
  loop_t.add(p, p, 4);

  return cfg;
}

int main(int argc, char **argv) {
  return crab_tests::test_main(argc, argv, [](bool stats_enabled) -> int {
    variable_factory_t vfac;
    auto cfg = prog(vfac);
    crab::outs() << *cfg << "\n";

    {
      z_interval_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_dbm_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_sdbm_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_ric_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_term_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }
    {
      z_dis_interval_domain_t init;
      run(cfg, init, stats_enabled, run_config().with_liveness());
    }

    return 0;
  });
}
