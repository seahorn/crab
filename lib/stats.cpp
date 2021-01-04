#include <crab/config.h>
#include <crab/support/stats.hpp>

namespace crab {
bool CrabStatsFlag = false;
void CrabEnableStats(bool v) { CrabStatsFlag = v; }
} // namespace crab

#ifndef CRAB_STATS
namespace crab {
Stopwatch::Stopwatch() {}
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

ScopedCrabStats::ScopedCrabStats(const std::string &name, bool reset)
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

std::map<std::string, unsigned> &CrabStats::getCounters() {
  static std::map<std::string, unsigned> counters;
  return counters;
}

std::map<std::string, Stopwatch> &CrabStats::getTimers() {
  static std::map<std::string, Stopwatch> timers;
  return timers;
}

void CrabStats::reset() {
  getCounters().clear();
  getTimers().clear();
}

void CrabStats::count(const std::string &name) {
  if (!crab::CrabStatsFlag)
    return;
  ++getCounters()[name];
}
void CrabStats::count_max(const std::string &name, unsigned v) {
  if (!crab::CrabStatsFlag)
    return;
  getCounters()[name] = std::max(getCounters()[name], v);
}

unsigned CrabStats::uset(const std::string &n, unsigned v) {
  if (!crab::CrabStatsFlag)
    return 0;
  return getCounters()[n] = v;
}
unsigned CrabStats::get(const std::string &n) {
  if (!crab::CrabStatsFlag)
    return 0;
  return getCounters()[n];
}

void CrabStats::start(const std::string &name) {
  if (!crab::CrabStatsFlag)
    return;
  getTimers()[name].start();
}
void CrabStats::stop(const std::string &name) {
  if (!crab::CrabStatsFlag)
    return;
  getTimers()[name].stop();
}
void CrabStats::resume(const std::string &name) {
  if (!crab::CrabStatsFlag)
    return;
  getTimers()[name].resume();
}

/** Outputs all statistics to std output */
void CrabStats::Print(crab_os &OS) {
  OS << "\n\n************** STATS ***************** \n";
  if (!crab::CrabStatsFlag) {
    OS << "Need to call CrabEnableStats()\n";
  } else {
    for (auto &kv : getCounters())
      OS << kv.first << ": " << kv.second << "\n";
    for (auto &kv : getTimers())
      OS << kv.first << ": " << kv.second << "\n";
  }
  OS << "************** STATS END ***************** \n";
}

void CrabStats::PrintBrunch(crab_os &OS) {
  OS << "\n\n************** BRUNCH STATS ***************** \n";
  if (!crab::CrabStatsFlag) {
    OS << "Need to call CrabEnableStats()\n";
  } else {
    for (auto &kv : getCounters())
      OS << "BRUNCH_STAT " << kv.first << " " << kv.second << "\n";
    for (auto &kv : getTimers())
      OS << "BRUNCH_STAT " << kv.first << " " << (kv.second).toSeconds()
         << "sec \n";
  }
  OS << "************** BRUNCH STATS END ***************** \n";
}

ScopedCrabStats::ScopedCrabStats(const std::string &name, bool reset)
    : m_name("") {
  if (crab::CrabStatsFlag) {
    m_name = name;
    if (reset) {
      m_name += ".last";
      CrabStats::start(m_name);
    } else {
      CrabStats::resume(m_name);
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
