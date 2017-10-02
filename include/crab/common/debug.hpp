#ifndef CRAB__DEBUG__HPP_
#define CRAB__DEBUG__HPP_

/* Logging and debug messages */

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


extern unsigned CrabVerbosity;
void CrabEnableVerbosity(unsigned v);
#define CRAB_VERBOSE_IF(LEVEL,CODE) do { if (::crab::CrabVerbosity >= LEVEL) { CODE; } } while (0) 
}
#endif
