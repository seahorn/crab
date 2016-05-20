#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

/* Example of how to build a CFG */
cfg_t* prog (VariableFactory &vfac)  {

  // Defining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  z_var nd (vfac ["nd"]);
  z_var inc (vfac ["inc"]);
  // entry and exit block
  auto cfg = new cfg_t("x0","ret");
  // adding blocks
  basic_block_t& x0 = cfg->insert ("x0");
  basic_block_t& x1 = cfg->insert ("x1");
  basic_block_t& x2 = cfg->insert ("x2");
  basic_block_t& x3 = cfg->insert ("x3");
  basic_block_t& entry = cfg->insert ("entry");
  basic_block_t& bb1   = cfg->insert ("bb1");
  basic_block_t& bb1_t = cfg->insert ("bb1_t");
  basic_block_t& bb1_f = cfg->insert ("bb1_f");
  basic_block_t& bb2   = cfg->insert ("bb2");
  basic_block_t& ret   = cfg->insert ("ret");
  // adding control flow
  x0 >> x1; x1 >> x2; x2 >> x3; x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  x0.assign (k, 2147483648);
  entry.assign (i, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.havoc(nd.name());
  bb2.select(inc,nd,1,2);
  bb2.add(i, i, inc);

  return cfg;
}

template<class CFG> 
std::vector<typename CFG::node_t> rev_topo(CFG cfg) {
  std::vector<typename CFG::node_t> rev_sccg_order;
  crab::analyzer::graph_algo::SccGraph <CFG> Scc_g (cfg);
  crab::analyzer::graph_algo::rev_topo_sort (Scc_g, rev_sccg_order);
  
  std::vector<typename CFG::node_t> order;
  for (auto &n : rev_sccg_order) {
    auto &members = Scc_g.getComponentMembers (n);
    order.insert (order.end (), members.begin (), members.end ());
  }
  return order;
}

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  VariableFactory vfac;
  cfg_t* cfg = prog(vfac);

  crab::outs() << "CFG\n" << *cfg << endl;

  std::cout << "Reversed topological order of CFG \n";
  bool first=true;
  for (auto &N: rev_topo(cfg_ref_t(*cfg))) {
    if (!first)
      std::cout << " -- ";
    std::cout << N;
    first = false;
  }
  std::cout << "\n";
  
  cfg_rev_t rev_cfg(*cfg);
  crab::outs() << "Reversed CFG\n" << rev_cfg << endl;

  std::cout << "Reversed topological order of the reversed CFG \n";
  first=true;
  for (auto &N: rev_topo(rev_cfg)) {
    if (!first)
      std::cout << " -- ";
    std::cout << N;
    first = false;
  }
  std::cout << "\n";
  
  cfg->simplify ();
  crab::outs() << "Simplified CFG\n" << *cfg << endl;

  rev_cfg = cfg_rev_t(*cfg);
  crab::outs() << "Reversed simplified CFG\n" << rev_cfg << endl;
  
  // free the CFG
  delete cfg;
  return 0;
}
