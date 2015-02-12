#include <ikos_cfg/Cfg.hpp>
#include <ikos_cfg/VarFactory.hpp>
#include <ikos_analysis/FwdAnalyzer.hpp>

#include <ikos_domains/common.hpp>
#include <ikos_domains/linear_constraints.hpp> 
#include <ikos_domains/intervals.hpp>                      
#include <ikos_domains/intervals_congruences.hpp>                      
#include <ikos_domains/octagons.hpp>                      
#include <ikos_domains/dbm.hpp>                      

using namespace std;

namespace cfg_impl
{
  using namespace cfg;

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
  typedef cfg_t::BasicBlock_t                  basic_block_t;
} // end namespace

namespace domain_impl
{
  using namespace cfg_impl;
  // Numerical domains
  typedef interval_domain< z_number, varname_t >             interval_domain_t;
  typedef interval_congruence_domain< z_number, varname_t >  interval_congruences_domain_t;
  typedef DBM< z_number, varname_t >                         dbm_domain_t;
  typedef octagon< z_number, varname_t >                     octagon_domain_t;
} // end namespace

using namespace cfg_impl;
using namespace domain_impl;
using namespace analyzer;

cfg_t prog (VariableFactory &vfac) 
{

  ////
  // Building the CFG
  ////

  // Definining program variables
  z_var i (vfac ["i"]);
  z_var k (vfac ["k"]);
  // entry and exit block
  cfg_t cfg ("x0","ret");
  // adding blocks
  basic_block_t& x0 = cfg.insert ("x0");
  basic_block_t& x1 = cfg.insert ("x1");
  basic_block_t& x2 = cfg.insert ("x2");
  basic_block_t& x3 = cfg.insert ("x3");
  basic_block_t& entry = cfg.insert ("entry");
  basic_block_t& bb1   = cfg.insert ("bb1");
  basic_block_t& bb1_t = cfg.insert ("bb1_t");
  basic_block_t& bb1_f = cfg.insert ("bb1_f");
  basic_block_t& bb2   = cfg.insert ("bb2");
  basic_block_t& ret   = cfg.insert ("ret");
  // adding control flow
  x0 >> x1; x1 >> x2; x2 >> x3; x3 >> entry;
  entry >> bb1;
  bb1 >> bb1_t; bb1 >> bb1_f;
  bb1_t >> bb2; bb2 >> bb1; bb1_f >> ret;
  // adding statements
  x0.assign (k, 50);
  entry.assign (i, 0);
  bb1_t.assume (i <= 99);
  bb1_f.assume (i >= 100);
  bb2.add(i, i, 1);

  return cfg;
}


int main (int argc, char** argv )
{
  VariableFactory vfac;
  cfg_t cfg = prog (vfac);
  cfg.simplify ();
  cout << cfg << endl;

  const bool run_live = true;
  interval_domain_t intervals = interval_domain_t::top ();
  FwdAnalyzer <basic_block_label_t, varname_t, cfg_t, VariableFactory, interval_domain_t> 
      itv_a (cfg, vfac, run_live);
  itv_a.Run (intervals);
  cout << "Results with intervals:\n";
  for (auto &b : cfg)
  {
    interval_domain_t inv = itv_a [b.label ()];
    std::cout << cfg_impl::get_str (b.label ()) << "=" << inv << "\n";
  }

  dbm_domain_t dbm = dbm_domain_t::top ();
  FwdAnalyzer <basic_block_label_t, varname_t, cfg_t, VariableFactory, dbm_domain_t> 
      dbm_a (cfg, vfac, run_live);
  dbm_a.Run (dbm);
  cout << "Results with DBMs:\n";
  for (auto &b : cfg)
  {
    dbm_domain_t inv = dbm_a [b.label ()];
    std::cout << cfg_impl::get_str (b.label ()) << "=" << inv << "\n";
  }

  return 0;
}
