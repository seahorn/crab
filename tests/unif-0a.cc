#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/analysis/FwdAnalyzer.hpp>

#include <ikos/common/types.hpp>
#include <ikos/algorithms/linear_constraints.hpp> 
#include <ikos/domains/intervals.hpp>                      
#include <ikos/domains/intervals_congruences.hpp>                      
#include <ikos/domains/octagons.hpp>                      
#include <ikos/domains/dbm.hpp>                      
#include <ikos/domains/term_equiv.hpp>

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
  typedef cfg_t::BasicBlock_t                  basic_block_t;
} // end namespace

/*
class StrVarAlloc {
public:
  StrVarAlloc()
    : vfac(new typename cfg_impl::StrVariableFactory)
  { }

  StrVarAlloc(StrVarAlloc& o)
    : vfac(o.vfac)
  { }
   
  cfg_impl::StrVariableFactory::varname_t next(void) {
    std::stringstream ss;
    ss << "x" << next_id++;
    return (*(vfac.get()))[ss.str()];
  }
protected:
  static int next_id;

  std::shared_ptr<cfg_impl::StrVariableFactory> vfac;
};
int StrVarAlloc::next_id = 0;
*/

class StrVarAlloc_col {
  static const char** col_prefix;
public:
  static cfg_impl::StrVariableFactory vfac;

  StrVarAlloc_col()
    : colour(0), next_id(0)
  { }

  StrVarAlloc_col(const StrVarAlloc_col& o)
    : colour(o.colour), next_id(o.next_id)
  { }

  StrVarAlloc_col(const StrVarAlloc_col& x, const StrVarAlloc_col& y)
    : colour(fresh_colour(x.colour, y.colour)),
      next_id(0)
  {
    assert(colour != x.colour);
    assert(colour != y.colour);
  }

  StrVarAlloc_col& operator=(const StrVarAlloc_col& x)
  {
    colour = x.colour;
    next_id = x.next_id;
    return *this;
  }

  cfg_impl::StrVariableFactory::varname_t next(void) {
    std::stringstream ss;
    ss << col_prefix[colour] << next_id++;
    return vfac[ss.str()];
  }
    
protected:
  int colour;
  int next_id;

  int fresh_colour(int col_x, int col_y)
  {
    switch(col_x)
    {
      case 0:
      {
        return col_y == 1 ? 2 : 1;
      }
      case 1:
      {
        return col_y == 0 ? 2 : 0;
      }
      case 2:
      {
        return col_y == 0 ? 1 : 0;
      }
      default:
        assert(0 && "Not reachable.");
        return 0;
    }
  }
};
static const char* col_prefix_data[] = { "_x", "_y", "_z" };
const char** StrVarAlloc_col::col_prefix = col_prefix_data;

cfg_impl::StrVariableFactory StrVarAlloc_col::vfac;


namespace domain_impl
{
  using namespace cfg_impl;
  // Numerical domains
  typedef interval_domain< z_number, varname_t >             interval_domain_t;
  /*
  typedef interval_congruence_domain< z_number, varname_t >  interval_congruences_domain_t;
  typedef DBM< z_number, varname_t >                         dbm_domain_t;
  typedef octagon< z_number, varname_t >                     octagon_domain_t;
  */

  class TDomInfo {
  public:
    typedef z_number Number;
    typedef varname_t VariableName;
    typedef StrVarAlloc_col Alloc;
    typedef interval_domain_t domain_t;
  };
  typedef anti_unif<TDomInfo>::anti_unif_t term_domain_t;
} // end namespace

using namespace cfg_impl;
using namespace domain_impl;
using namespace analyzer;

int main (int argc, char** argv )
{
  VariableFactory vfac;

  varname_t x = vfac["x"];
  varname_t y = vfac["y"];

  term_domain_t dom_left = term_domain_t::top ();
  dom_left.assign(x, z_number(1));

  term_domain_t dom_right = dom_left;
  dom_right.apply(OP_ADDITION, x, x, z_number(1));

  term_domain_t l_join_r = dom_left | dom_right;

  std::cout << dom_left << " | " << dom_right << " = " << l_join_r << std::endl;
//  std::cout << dom_left;
  return 0;
}
