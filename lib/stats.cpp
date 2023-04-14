#include <crab/config.h>
#include <crab/support/stats.hpp>
#include <vector>

namespace crab {
bool CrabStatsFlag = false;
void CrabEnableStats(bool v) { CrabStatsFlag = v; }
} // namespace crab

#ifndef CRAB_STATS
namespace crab {
Stopwatch::Stopwatch() {
  (void) started;
  (void) finished;
  (void) timeElapsed;
}
void Stopwatch::start() {}
void Stopwatch::stop() {}
void Stopwatch::resume() {}
long Stopwatch::systemTime() const { return (long)0; }
long Stopwatch::getTimeElapsed() const { return (long)0; }
double Stopwatch::toSeconds() { return (double)0; }
void Stopwatch::Print(crab_os &out) const {}
void CrabStats::reset() {}
void CrabStats::count(const std::string &name) {}
void CrabStats::count_max(const std::string &name, unsigned v) {}
unsigned CrabStats::uset(const std::string &n, unsigned v) {
  return (unsigned)0;
}
unsigned CrabStats::get(const std::string &n) { return (unsigned)0; }
void CrabStats::start(const std::string &name) {}
void CrabStats::stop(const std::string &name) {}
void CrabStats::resume(const std::string &name) {}

/** Outputs all statistics to std output */
void CrabStats::Print(crab_os &OS) {
  OS << "\n\n************** STATS ***************** \n";
  OS << "Crab compiled with support for gathering stats. "
     << "Compile Crab with -DCRAB_ENABLE_STATS=ON\n";
  OS << "************** STATS END ***************** \n";
}

void CrabStats::PrintBrunch(crab_os &OS) {
  OS << "\n\n************** BRUNCH STATS ***************** \n";
  OS << "Crab compiled with support for gathering stats. "
     << "Compile Crab with -DCRAB_ENABLE_STATS=ON\n";
  OS << "************** BRUNCH STATS END ***************** \n";
}

ScopedCrabStats::ScopedCrabStats(std::string &&name, const char* suffix, bool use_count)
  : m_name("") {}

ScopedCrabStats::ScopedCrabStats(const char* name, bool use_count)
  : m_name("") {}
  
ScopedCrabStats::~ScopedCrabStats() {}
} // namespace crab
#else
/* Real implementation stars here */
#include <sys/resource.h>
#include <sys/time.h>

namespace crab {

long Stopwatch::systemTime() const {
  struct rusage ru;
  getrusage(RUSAGE_SELF, &ru);
  long r = ru.ru_utime.tv_sec * 1000000L + ru.ru_utime.tv_usec;
  return r;
}

Stopwatch::Stopwatch() { start(); }

void Stopwatch::start() {
  started = systemTime();
  finished = -1;
  timeElapsed = 0;
}

void Stopwatch::stop() {
  if (finished < started) {
    finished = systemTime();
  }
}

void Stopwatch::resume() {
  if (finished >= started) {
    timeElapsed += finished - started;
    started = systemTime();
    finished = -1;
  }
}

long Stopwatch::getTimeElapsed() const {
  if (finished < started)
    return timeElapsed + systemTime() - started;
  else
    return timeElapsed + finished - started;
}

double Stopwatch::toSeconds() {
  double time = ((double)getTimeElapsed() / 1000000);
  return time;
}

void Stopwatch::Print(crab_os &out) const {
  long time = getTimeElapsed();
  long h = time / 3600000000L;
  long m = time / 60000000L - h * 60;
  float s = ((float)time / 1000000L) - m * 60 - h * 3600;

  if (h > 0)
    out << h << "h";
  if (m > 0)
    out << m << "m";
  out << s << "s";
}
  

std::unordered_map<std::string, unsigned>& CrabStats::getCounters() {
  static std::unordered_map<std::string, unsigned> counters;
  return counters;
}

std::unordered_map<std::string, Stopwatch>& CrabStats::getTimers() {
  static std::unordered_map<std::string, Stopwatch> timers;
  return timers;
}
    
void CrabStats::reset() {
  if (crab::CrabStatsFlag) {  
    getCounters().clear();
    getTimers().clear();
  }
}

void CrabStats::count(const std::string &name) {
  if (crab::CrabStatsFlag) {
    ++getCounters()[name];
  }
}
void CrabStats::count_max(const std::string &name, unsigned v) {
  if (crab::CrabStatsFlag) {  
    getCounters()[name] = std::max(getCounters()[name], v);
  }
}

unsigned CrabStats::uset(const std::string &n, unsigned v) {
  if (crab::CrabStatsFlag) {  
    return getCounters()[n] = v;
  } else {
    return 0;
  }
}
unsigned CrabStats::get(const std::string &n) {
  if (crab::CrabStatsFlag) {  
    return getCounters()[n];
  } else {
    return 0;
  }
}

void CrabStats::start(const std::string &name) {
  if (crab::CrabStatsFlag) {  
    getTimers()[name].start();
  }
}
void CrabStats::stop(const std::string &name) {
  if (crab::CrabStatsFlag) {  
    getTimers()[name].stop();
  }
}
void CrabStats::resume(const std::string &name) {
  if (crab::CrabStatsFlag) {  
    getTimers()[name].resume();
  }
}

/** Outputs all statistics to std output */
void CrabStats::Print(crab_os &OS) {
  if (crab::CrabStatsFlag) {  
    std::vector<std::pair<std::string, unsigned>> sorted_counters;
    sorted_counters.reserve(getCounters().size());
    std::vector<std::pair<std::string, Stopwatch>> sorted_timers;
    sorted_timers.reserve(getTimers().size());
    for (auto &kv : getCounters())
      sorted_counters.push_back(kv);
    for (auto &kv : getTimers())
      sorted_timers.push_back(kv);
    std::sort(sorted_counters.begin(), sorted_counters.end(),
	      [](const std::pair<std::string, unsigned> &p1,
		 const std::pair<std::string, unsigned> &p2) {
		return p1.first < p2.first;
	      });
    std::sort(sorted_timers.begin(), sorted_timers.end(),
	      [](const std::pair<std::string, Stopwatch> &p1,
		 const std::pair<std::string, Stopwatch> &p2) {
		return p1.first < p2.first;
	      });
    OS << "\n\n************** STATS ***************** \n";
    for (auto &kv : sorted_counters)
      OS << kv.first << ": " << kv.second << "\n";
    for (auto &kv : sorted_timers)
      OS << kv.first << ": " << kv.second << "\n";
    OS << "************** STATS END ***************** \n";
  } else {
    OS << "\n\n************** STATS ***************** \n";
    OS << "Need to call CrabEnableStats()\n";    
    OS << "************** STATS END ***************** \n";    
  } 
}

void CrabStats::PrintBrunch(crab_os &OS) {
  if (crab::CrabStatsFlag) {    
    std::vector<std::pair<std::string, unsigned>> sorted_counters;
    sorted_counters.reserve(getCounters().size());
    std::vector<std::pair<std::string, Stopwatch>> sorted_timers;
    sorted_timers.reserve(getTimers().size());
    for (auto &kv : getCounters())
      sorted_counters.push_back(kv);
    for (auto &kv : getTimers())
      sorted_timers.push_back(kv);
    std::sort(sorted_counters.begin(), sorted_counters.end(),
	      [](const std::pair<std::string, unsigned> &p1,
		 const std::pair<std::string, unsigned> &p2) {
		return p1.first < p2.first;
	      });
    std::sort(sorted_timers.begin(), sorted_timers.end(),
	      [](const std::pair<std::string, Stopwatch> &p1,
		 const std::pair<std::string, Stopwatch> &p2) {
		return p1.first < p2.first;
	      });
    
    OS << "\n\n************** BRUNCH STATS ***************** \n";
    for (auto &kv : sorted_counters)
      OS << "BRUNCH_STAT " << kv.first << " " << kv.second << "\n";
    for (auto &kv : sorted_timers)
      OS << "BRUNCH_STAT " << kv.first << " " << (kv.second).toSeconds()
	 << "sec \n";
    OS << "************** BRUNCH STATS END ***************** \n";
  } else {
    OS << "\n\n************** STATS ***************** \n";
    OS << "Need to call CrabEnableStats()\n";    
    OS << "************** STATS END ***************** \n";        
  }
}

ScopedCrabStats::ScopedCrabStats(std::string &&name, const char* suffix, bool use_count)
    : m_name(std::move(name)) {
  if (crab::CrabStatsFlag) {
    m_name.append(suffix);
    CrabStats::resume(m_name);
    if (use_count) {
      CrabStats::count(m_name);
    }
  }
}

ScopedCrabStats::ScopedCrabStats(const char* name, bool use_count)
    : m_name(name) {
  if (crab::CrabStatsFlag) {
    CrabStats::resume(m_name);
    if (use_count) {
      CrabStats::count(m_name);
    }
  }
}
  
ScopedCrabStats::~ScopedCrabStats() {
  if (crab::CrabStatsFlag) {  
    CrabStats::stop(m_name);
  }
}

} // namespace crab
#endif
