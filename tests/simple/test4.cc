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
    crab::outs() << "Vertex: " << v << endl; 
    crab::outs() << "Num of predecessors=" << in_degree (v, g) << endl;
    crab::outs() << "Num of successors  =" << out_degree (v, g) << endl;
    crab::outs() << "Num of neighbors   =" << degree (v, g) << endl;
    crab::outs() << "Succs={";
    {
      auto p = out_edges (v, g);
      auto succIt  = p.first;
      auto succEnd = p.second;
      for(; succIt != succEnd; ++succIt){
        crab::outs() << target (*succIt, g) << ";";
      }
    }
    crab::outs() << "}" << endl;
    crab::outs() << "Preds={ ";
    {
      auto p = in_edges (v, g);
      auto predIt  = p.first;
      auto predEnd = p.second;
      for(; predIt != predEnd; ++predIt){
        crab::outs() << source (*predIt, g) << ";";
      }
    }
    crab::outs() << "}" << endl;
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


int main (int argc, char** argv )
{
  SET_LOGGER(argc,argv)

  VariableFactory vfac;
  cfg_t* cfg = prog (vfac);
  cfg->simplify ();
  crab::outs() << *cfg << endl;
  write (*cfg);
  crab::outs() << endl;

  Liveness<cfg_ref_t> live (*cfg);
  live.exec ();

  {
    NumFwdAnalyzer <cfg_ref_t, interval_domain_t,VariableFactory>::type a (*cfg,vfac,&live);
    interval_domain_t inv = interval_domain_t::top ();
    a.Run (inv);
    crab::outs() << "Invariants using " << interval_domain_t::getDomainName () << "\n";
    for (auto &b : *cfg)
    {
      auto inv = a [b.label ()];
      crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_ref_t, dbm_domain_t,VariableFactory>::type a (*cfg,vfac,&live);
    dbm_domain_t inv = dbm_domain_t::top ();
    a.Run (inv);
    crab::outs() << "Invariants using " << dbm_domain_t::getDomainName () << "\n";
    for (auto &b : *cfg)
    {
      auto inv = a [b.label ()];
      crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_ref_t, sdbm_domain_t,VariableFactory>::type a (*cfg,vfac,&live);
    sdbm_domain_t inv = sdbm_domain_t::top ();
    a.Run (inv);
    crab::outs() << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : *cfg)
    {
      auto inv = a [b.label ()];
      crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_ref_t, ric_domain_t,VariableFactory>::type a (*cfg,vfac,&live);
    ric_domain_t inv = ric_domain_t::top ();
    a.Run (inv);
    crab::outs() << "Invariants using " << ric_domain_t::getDomainName () << "\n";
    for (auto &b : *cfg)
    {
      auto inv = a [b.label ()];
      crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_ref_t, term_domain_t,VariableFactory>::type a (*cfg,vfac,&live);
    term_domain_t inv = term_domain_t::top ();
    a.Run (inv);
    crab::outs() << "Invariants using " << term_domain_t::getDomainName () << "\n";
    for (auto &b : *cfg)
    {
      auto inv = a [b.label ()];
      crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_ref_t, dis_interval_domain_t,VariableFactory>::type 
        a (*cfg,vfac,&live, 1, 2, 20);
    a.Run (dis_interval_domain_t::top ());
    crab::outs() << "Invariants using " << dis_interval_domain_t::getDomainName () << "\n";
    for (auto &b : *cfg) 
    {
      auto inv = a [b.label ()];
      crab::outs() << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  delete cfg;
  return 0;
}
