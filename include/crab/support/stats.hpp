#pragma once

#include <crab/support/os.hpp>

#include <unordered_map>
#include <string>

namespace crab {

extern bool CrabStatsFlag;
void CrabEnableStats(bool v = true);
  
class Stopwatch {
  long started;
  long finished;
  long timeElapsed;

  long systemTime() const;

public:
  Stopwatch();
  void start();
  void stop();
  void resume();
  long getTimeElapsed() const;
  void Print(crab_os &out) const;
  double toSeconds();
};

inline crab_os &operator<<(crab_os &OS, const Stopwatch &sw) {
  sw.Print(OS);
  return OS;
}

class CrabStats {
  
  static std::unordered_map<std::string, unsigned> &getCounters();
  static std::unordered_map<std::string, Stopwatch> &getTimers();

public:
  static void reset();

  /* counters */
  static unsigned get(const std::string &n);
  static unsigned uset(const std::string &n, unsigned v);
  static void count(const std::string &name);
  static void count_max(const std::string &name, unsigned v);

  /* stop watch */
  static void start(const std::string &name);
  static void stop(const std::string &name);
  static void resume(const std::string &name);

  /** Outputs all statistics to std output */
  static void Print(crab_os &OS);
  static void PrintBrunch(crab_os &OS);
};

class ScopedCrabStats {
  std::string m_name;

public:
  // Call count and resume on name+suffix
  ScopedCrabStats(std::string &&name, const char* suffix, bool use_count = true);
  // Call count and resume on name
  ScopedCrabStats(const char* name, bool use_count = true);
  ~ScopedCrabStats();
};
} // namespace crab


/**
 * Some convenient macros
 *
 *   CRAB_SCOPED_STATS(name, active)
 *     increase **both** timer and counter for name if active=1
 *   CRAB_COUNT_STATS(name, active)
 *     increment only counter for name if active=1
 *   CRAB_SCOPED_TIMER_STATS(name, active)
 *     increment only timer for name if active=1
 *
 * The CRAB_DOMAIN_* variants below prefix the name with an abstract
 * domain's domain_name():
 *
 *   CRAB_DOMAIN_SCOPED_STATS(who, suffix, active)
 *     increase **both** timer and counter for who->domain_name() + suffix
 *     if active=1. This macro avoids string creation if active=0.
 *     Moreover, it concats efficiently two strings.
 *   CRAB_DOMAIN_COUNT_STATS(suffix, active)
 *     increase only counter for this->domain_name() + suffix if active=1.
 *     This macro avoids string creation if active=0, but it concatenates
 *     strings in an expensive way.
**/
#include <crab/config.h>
#ifdef CRAB_STATS
#define CRAB_SCOPED_STATS(name, active) \
  CRAB_SCOPED_STATS_(name, active)
#define CRAB_SCOPED_STATS_(name, active) \
  CRAB_SCOPED_STATS_ ## active(name)
#define CRAB_SCOPED_STATS_0(name) 
#define CRAB_SCOPED_STATS_1(name) \
  crab::ScopedCrabStats __st__(name, true);

#define CRAB_COUNT_STATS(name, active)		\
  CRAB_COUNT_STATS_(name, active)
#define CRAB_COUNT_STATS_(name, active) \
  CRAB_COUNT_STATS_ ## active(name)
#define CRAB_COUNT_STATS_0(name) 
#define CRAB_COUNT_STATS_1(name) \
  crab::CrabStats::count(name);

#define CRAB_SCOPED_TIMER_STATS(name, active) \
  CRAB_SCOPED_TIMER_STATS_(name, active)
#define CRAB_SCOPED_TIMER_STATS_(name, active) \
  CRAB_SCOPED_TIMER_STATS_ ## active(name)
#define CRAB_SCOPED_TIMER_STATS_0(name) 
#define CRAB_SCOPED_TIMER_STATS_1(name) \
  crab::ScopedCrabStats __st__(name, false);

#define CRAB_DOMAIN_SCOPED_STATS(who, suffix, active) \
  CRAB_DOMAIN_SCOPED_STATS_(who, suffix, active)
#define CRAB_DOMAIN_SCOPED_STATS_(who, suffix, active) \
  CRAB_DOMAIN_SCOPED_STATS_ ## active(who, suffix)
#define CRAB_DOMAIN_SCOPED_STATS_0(who, suffix) 
#define CRAB_DOMAIN_SCOPED_STATS_1(who, suffix) \
  crab::ScopedCrabStats __st__((who)->domain_name(), suffix, true);

#define CRAB_DOMAIN_COUNT_STATS(suffix, active) \
  CRAB_DOMAIN_COUNT_STATS_(suffix, active)
#define CRAB_DOMAIN_COUNT_STATS_(suffix, active) \
  CRAB_DOMAIN_COUNT_STATS_ ## active(suffix)
#define CRAB_DOMAIN_COUNT_STATS_0(suffix) 
#define CRAB_DOMAIN_COUNT_STATS_1(suffix) \
  crab::CrabStats::count(this->domain_name() + suffix);
#else
#define CRAB_SCOPED_STATS(name, active) 
#define CRAB_COUNT_STATS(name, active)
#define CRAB_SCOPED_TIMER_STATS(name, active)
#define CRAB_DOMAIN_SCOPED_STATS(who, suffix, active)
#define CRAB_DOMAIN_COUNT_STATS(suffix, active)
#endif
