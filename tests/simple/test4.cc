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

cfg_t* prog (VariableFactory &vfac) 
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


template <typename Domain, typename Live>
void run(cfg_ref_t cfg, VariableFactory &vfac, Live *live)
{
  typename NumFwdAnalyzer <cfg_ref_t, Domain, VariableFactory>::type 
      It (cfg, vfac, live, 1, 2, 20);
  Domain inv = Domain::top ();
  It.Run (inv);
  crab::outs() << "Invariants using " << Domain::getDomainName () << ":\n";

  for (auto &b : cfg)
  {
    // invariants at the entry of the block
    auto inv = It [b.label ()];
    crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
  }
  crab::outs() << "\n";
  if (stats_enabled) {
    crab::CrabStats::Print(crab::outs());
    crab::CrabStats::reset();
  }  
}

int main (int argc, char** argv )
{
  SET_TEST_OPTIONS(argc,argv)

  VariableFactory vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << "\n";
  //write (*cfg);
  //crab::outs() << "\n";

  Liveness<cfg_ref_t> live (*cfg);
  live.exec ();

  run<interval_domain_t>(*cfg,vfac,&live);
  run<dbm_domain_t>(*cfg,vfac,&live);
  run<sdbm_domain_t>(*cfg,vfac,&live);
  run<ric_domain_t>(*cfg,vfac,&live);
  run<term_domain_t>(*cfg,vfac,&live);
  run<dis_interval_domain_t>(*cfg,vfac,&live);
  delete cfg;
  return 0;
}
