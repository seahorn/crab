#include <crab/common/debug.hpp>

#ifndef NCRABLOG
#include <set>
using namespace crab;
bool crab::CrabLogFlag = false;
std::set<std::string> crab::CrabLog;

void crab::CrabEnableLog (std::string x) 
{
  if (x.empty ()) return;
  CrabLogFlag = true;
  CrabLog.insert (x); 
}

namespace crab
{
  // struct LogOpt
  // { void operator=(const std::string &tag) const 
  //   { 
  //     CrabEnableLog (tag); 
  //   } 
  // };  

  // LogOpt loc;

} // end namespace crab

#else
void crab::CrabEnableLog (std::string x) { }
#endif

using namespace crab;
unsigned crab::CrabVerbosity = 0;
void crab::CrabEnableVerbosity (unsigned v) { CrabVerbosity=v;}

