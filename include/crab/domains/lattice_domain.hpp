#pragma once

#include <crab/support/os.hpp>
#include <string>

namespace crab {
namespace domains {  
template<class LatticeDomain>
class lattice_domain_api {
public:
  virtual ~lattice_domain_api() {}
  virtual LatticeDomain make_top() const = 0;
  virtual LatticeDomain make_bottom() const = 0;
  virtual void set_to_top() = 0;
  virtual void set_to_bottom() = 0;
  virtual bool is_top() const = 0;
  virtual bool is_bottom() const = 0;
  virtual bool operator<=(const LatticeDomain &other) const = 0;
  virtual void operator|=(const LatticeDomain &other) = 0;  
  virtual LatticeDomain operator|(const LatticeDomain &other) const = 0;
  virtual LatticeDomain operator||(const LatticeDomain &other) const = 0;
  virtual LatticeDomain operator&(const LatticeDomain &other) const = 0;
  virtual LatticeDomain operator&&(const LatticeDomain &other) const = 0;
  virtual void write(crab::crab_os &o) const = 0;
  virtual std::string domain_name() const = 0;
};
} // end namespace domains
} // end namespace crab
