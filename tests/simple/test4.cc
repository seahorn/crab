#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

// To test the interface to BGL graphs
void write (cfg_ref_t g)
{
  crab::outs() << "Num of vertices: " << num_vertices (g) << "\n";
  for (auto v: boost::make_iterator_range (vertices (g)))
  {
    crab::outs() << "Vertex: " << v << "\n"; 
    crab::outs() << "Num of predecessors=" << in_degree (v, g) << "\n";
    crab::outs() << "Num of successors  =" << out_degree (v, g) << "\n";
    crab::outs() << "Num of neighbors   =" << degree (v, g) << "\n";
    crab::outs() << "Succs={";
    {
      auto p = out_edges (v, g);
      auto succIt  = p.first;
      auto succEnd = p.second;
      for(; succIt != succEnd; ++succIt){
        crab::outs() << target (*succIt, g) << ";";
      }
    }
    crab::outs() << "}" << "\n";
    crab::outs() << "Preds={ ";
    {
      auto p = in_edges (v, g);
      auto predIt  = p.first;
      auto predEnd = p.second;
      for(; predIt != predEnd; ++predIt){
        crab::outs() << source (*predIt, g) << ";";
      }
    }
    crab::outs() << "}" << "\n";
  }
}

cfg_t* prog (variable_factory_t &vfac) 
{

  cfg_t* cfg = new cfg_t("entry","ret");
  basic_block_t& entry      = cfg->insert ("entry");
  basic_block_t& loop_head  = cfg->insert ("loop_head");
  basic_block_t& loop_t     = cfg->insert ("loop_t");
  basic_block_t& loop_f     = cfg->insert ("loop_f");
  basic_block_t& loop_body  = cfg->insert ("loop_body");
  basic_block_t& ret        = cfg->insert ("ret");

  entry >> loop_head;
  loop_head >> loop_t; 
  loop_head >> loop_f; 
  loop_t >> loop_body; 
  loop_body >> loop_head;
  loop_f >> ret;

  z_var i(vfac["i"]);
  z_var p(vfac["p"]);

  entry.assign (i, 0);
  entry.assign (p, 0);

  loop_t.assume (i <= 9);
  loop_f.assume (i >= 10);
  loop_body.add (i, i, 1);
  loop_body.add (p, p, 4);

  return cfg;
}

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  variable_factory_t vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << "\n";

  run<interval_domain_t>(cfg,vfac,true,1,2,20);
  run<dbm_domain_t>(cfg,vfac,true,1,2,20);
  run<sdbm_domain_t>(cfg,vfac,true,1,2,20);
  run<ric_domain_t>(cfg,vfac,true,1,2,20);
  run<term_domain_t>(cfg,vfac,true,1,2,20);
  run<dis_interval_domain_t>(cfg,vfac,true,1,2,20);
  delete cfg;
  return 0;
}
