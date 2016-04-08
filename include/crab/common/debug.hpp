#ifndef CRAB__DEBUG__HPP_
#define CRAB__DEBUG__HPP_

/* Code from Avy */

#include <string>
#include <set>

namespace crab
{
#ifndef NCRABLOG
#define CRAB_LOG(TAG,CODE) do { if (::crab::CrabLogFlag && ::crab::CrabLog.count (TAG) > 0) { CODE; } } while (0)
extern bool CrabLogFlag;
extern std::set<std::string> CrabLog;
void CrabEnableLog (std::string x);
#else
#define CRAB_LOG(TAG,CODE) do { } while (0)
void CrabEnableLog (std::string x);
#endif
}
#endif
