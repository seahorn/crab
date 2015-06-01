#ifndef __IKOS_TERM_UTIL_H__
#define __IKOS_TERM_UTIL_H__
// Some boilerplate stuff.
#include <ikos/cfg/Cfg.hpp>
#include <ikos/cfg/VarFactory.hpp>
#include <ikos/domains/term_equiv.hpp>

namespace ikos {
  namespace term {

// Three-coloured variable allocation
// So the number of variables is bounded by 3|Tbl|,
// rather than always increasing.
class StrVariableFactory : public boost::noncopyable  
{
  typedef cfg::var_factory_impl::VariableFactory< std::string > StrVariableFactory_t;
  std::unique_ptr< StrVariableFactory_t > m_factory; 
  
 public: 

  typedef StrVariableFactory_t::variable_t varname_t;

  StrVariableFactory(): m_factory (new StrVariableFactory_t()){ }

  varname_t operator[](std::string v)

  { return (*m_factory)[v];}
}; 

class IntVariableFactory : public boost::noncopyable  
{
 public: 
  typedef int varname_t;

  IntVariableFactory() { }

  varname_t operator[](int v)
  { return v; }
}; 

inline int fresh_colour(int col_x, int col_y)
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

class StrVarAlloc_col {
  static const char** col_prefix;
public:
  typedef StrVariableFactory::varname_t varname_t;
  static StrVariableFactory vfac;

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

  StrVariableFactory::varname_t next(void) {
    std::stringstream ss;
    ss << col_prefix[colour] << next_id++;
    return vfac[ss.str()];
  }
    
protected:
  int colour;
  int next_id;
};

// GKG: This is likely to fail horribly if this header
//      is included in multiple source files.
StrVariableFactory StrVarAlloc_col::vfac;
static const char* col_prefix_data[] = { "_x", "_y", "_z" };
const char** StrVarAlloc_col::col_prefix = col_prefix_data;

class IntVarAlloc_col {
public:
  static IntVariableFactory vfac;
  typedef int varname_t;

  IntVarAlloc_col()
    : colour(0), next_id(0)
  { }

  IntVarAlloc_col(const IntVarAlloc_col& o)
    : colour(o.colour), next_id(o.next_id)
  { }

  IntVarAlloc_col(const IntVarAlloc_col& x, const IntVarAlloc_col& y)
    : colour(fresh_colour(x.colour, y.colour)),
      next_id(0)
  {
    assert(colour != x.colour);
    assert(colour != y.colour);
  }

  IntVarAlloc_col& operator=(const IntVarAlloc_col& x)
  {
    colour = x.colour;
    next_id = x.next_id;
    return *this;
  }

  IntVariableFactory::varname_t next(void) {
    int id = next_id++;
    return 3*id + colour;
  }
    
protected:
  int colour;
  int next_id;

};


template<class Num, class VName, class Abs>
class TDomInfo {
  public:
    typedef Num Number;
    typedef VName VariableName;
    typedef StrVarAlloc_col Alloc;
    typedef Abs domain_t;
};

  } // namespace term
} // namespace ikos

#endif
