#include "boost/tuple/tuple.hpp"

#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/FwdAnalyzer.hpp>
#include <ikos/cfg/CfgBgl.hpp>

#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp> 
#include <ikos/domains/intervals.hpp>                      
#include <ikos/domains/intervals_congruences.hpp>                      
#include <ikos/domains/octagons.hpp>                      
#include <ikos/domains/dbm.hpp>                      

using namespace std;

namespace cfg_impl
{
  using namespace cfg;

  template<> inline std::string get_label_str(std::string e) 
  { return e; }

  class StrVariableFactory : public boost::noncopyable  
  {
    typedef var_factory_impl::VariableFactory< std::string > StrVariableFactory_t;
    std::unique_ptr< StrVariableFactory_t > m_factory; 
    
   public: 

    typedef StrVariableFactory_t::variable_t varname_t;

    StrVariableFactory(): m_factory (new StrVariableFactory_t()){ }

    varname_t operator[](std::string v)

    { return (*m_factory)[v];}
  }; 

  // A variable factory based on strings
  typedef StrVariableFactory VariableFactory;
  typedef typename VariableFactory::varname_t varname_t;

  // CFG
  typedef variable< z_number, varname_t >      z_var;
  typedef std::string                          basic_block_label_t;
  typedef Cfg< basic_block_label_t, varname_t> cfg_t;
  typedef cfg_t::basic_block_t                 basic_block_t;
} // end namespace

namespace domain_impl
{
  using namespace ikos;
  // Numerical domains
  typedef interval_domain< z_number, cfg_impl::varname_t > interval_domain_t;
  typedef interval_congruence_domain< z_number, cfg_impl::varname_t > interval_congruences_domain_t;
  typedef DBM< z_number, cfg_impl::varname_t > dbm_domain_t;
  typedef octagon< z_number, cfg_impl::varname_t > octagon_domain_t;
} // end namespace


// To test the interface to BGL graphs
void write (cfg_impl::cfg_t g)
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

using namespace cfg_impl;
using namespace domain_impl;
using namespace analyzer;

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

  const bool run_live = true;

  NumFwdAnalyzer <cfg_t, octagon_domain_t>::type fwd_anal (cfg, run_live);
  fwd_anal.Run (octagon_domain_t::top ());
  cout << "Results:\n";
  for (auto &b : cfg)
  {
    octagon_domain_t inv = fwd_anal [b.label ()];
    std::cout << cfg_impl::get_label_str (b.label ()) << "=" << inv << "\n";
  }

  return 0;
}
