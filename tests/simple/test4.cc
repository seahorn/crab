#include "../common.hpp"

using namespace std;
using namespace crab::analyzer;
using namespace crab::cfg_impl;
using namespace crab::domain_impl;

// To test the interface to BGL graphs
void write (cfg_t g)
{
  cout << "Num of vertices: " << boost::num_vertices (g) << "\n";
  for (auto v: boost::make_iterator_range (boost::vertices (g)))
  {
    cout << "Vertex: " << v << endl; 
    cout << "Num of predecessors=" << boost::in_degree (v, g) << endl;
    cout << "Num of successors  =" << boost::out_degree (v, g) << endl;
    cout << "Num of neighbors   =" << boost::degree (v, g) << endl;
    cout << "Succs={";
    {
      auto p = boost::out_edges (v, g);
      auto succIt  = p.first;
      auto succEnd = p.second;
      for(; succIt != succEnd; ++succIt){
        cout << boost::target (*succIt, g) << ";";
      }
    }
    cout << "}" << endl;
    cout << "Preds={ ";
    {
      auto p = boost::in_edges (v, g);
      auto predIt  = p.first;
      auto predEnd = p.second;
      for(; predIt != predEnd; ++predIt){
        cout << boost::source (*predIt, g) << ";";
      }
    }
    cout << "}" << endl;
  }
}

cfg_t prog (VariableFactory &vfac) 
{

  cfg_t cfg ("entry","ret");
  basic_block_t& entry      = cfg.insert ("entry");
  basic_block_t& loop_head  = cfg.insert ("loop_head");
  basic_block_t& loop_t     = cfg.insert ("loop_t");
  basic_block_t& loop_f     = cfg.insert ("loop_f");
  basic_block_t& loop_body  = cfg.insert ("loop_body");
  basic_block_t& ret        = cfg.insert ("ret");

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
  VariableFactory vfac;
  cfg_t cfg = prog (vfac);
  cfg.simplify ();
  cout << cfg << endl;
  write (cfg);
  cout << endl;

  Liveness<cfg_t> live (cfg);
  live.exec ();

  {
    NumFwdAnalyzer <cfg_t, interval_domain_t,VariableFactory>::type a (cfg,vfac,&live);
    interval_domain_t inv = interval_domain_t::top ();
    a.Run (inv);
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg)
    {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_t, dbm_domain_t,VariableFactory>::type a (cfg,vfac,&live);
    dbm_domain_t inv = dbm_domain_t::top ();
    a.Run (inv);
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg)
    {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_t, ric_domain_t,VariableFactory>::type a (cfg,vfac,&live);
    ric_domain_t inv = ric_domain_t::top ();
    a.Run (inv);
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg)
    {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  {
    NumFwdAnalyzer <cfg_t, term_domain_t,VariableFactory>::type a (cfg,vfac,&live);
    term_domain_t inv = term_domain_t::top ();
    a.Run (inv);
    cout << "Invariants using " << inv.getDomainName () << "\n";
    for (auto &b : cfg)
    {
      auto inv = a [b.label ()];
      std::cout << get_label_str (b.label ()) << "=" << inv << "\n";
    }
  }

  return 0;
}
