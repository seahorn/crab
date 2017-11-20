#include <crab/common/debug.hpp>

#ifndef NCRABLOG
#include <set>

namespace crab {
  bool CrabLogFlag = false;
  std::set<std::string> CrabLog;

  void CrabEnableLog (std::string x)  {
    if (x.empty ()) return;
    CrabLogFlag = true;
    CrabLog.insert (x);   
  }
}
  

#else
namespace crab {
  void CrabEnableLog (std::string x) { }
}
#endif

namespace crab {
  unsigned CrabVerbosity = 0;
  void CrabEnableVerbosity (unsigned v) { CrabVerbosity=v;}
}

