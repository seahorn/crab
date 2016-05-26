#include <crab/config.h>

#ifdef HAVE_STATS
#include "crab/common/stats.hpp"

namespace crab
{
  std::map<std::string,unsigned> CrabStats::counters;
  std::map<std::string,Stopwatch> CrabStats::sw;
  std::map<std::string,Averager> CrabStats::av;
  std::map<std::string,std::string> CrabStats::ss;

  void CrabStats::reset () {
    counters.clear();
    sw.clear();
    av.clear();
    ss.clear();
  }

  void CrabStats::count (const std::string &name) { ++counters[name]; }
  void CrabStats::count_max (const std::string &name, unsigned v) {
      counters[name] = std::max (counters[name], v);
  }

  double CrabStats::avg (const std::string &n, double v) { return av[n].add (v); }
  unsigned CrabStats::uset (const std::string &n, unsigned v)
  { return counters [n] = v; }
  unsigned CrabStats::get (const std::string &n) { return counters [n]; }

  void CrabStats::sset (const std::string &n, std::string v) {ss [n] = v;}
  std::string& CrabStats::sget (const std::string &n) {return ss[n];}
  
  void CrabStats::start (const std::string &name) { sw[name].start (); }
  void CrabStats::stop (const std::string &name) { sw[name].stop (); }
  void CrabStats::resume (const std::string &name) { sw[name].resume (); }

  /** Outputs all statistics to std output */
  void CrabStats::Print (crab_os &OS) {
    OS << "\n\n************** STATS ***************** \n";
    for (auto &kv : ss)
      OS << kv.first << ": " << kv.second << "\n";
    for (auto &kv : counters)
      OS << kv.first << ": " << kv.second << "\n";

    for (auto &kv : sw)
      OS << kv.first << ": " << kv.second << "\n";

    for (auto &kv : av)
      OS << kv.first << ": " << kv.second << "\n";


    OS << "************** STATS END ***************** \n";
  }

  void CrabStats::PrintBrunch (crab_os &OS)
  {
    OS << "\n\n************** BRUNCH STATS ***************** \n";
    for (auto &kv : ss) 
      OS << "BRUNCH_STAT " << kv.first << " " << kv.second << "\n";
    
    for (auto &kv : counters)
      OS << "BRUNCH_STAT " << kv.first << " " << kv.second << "\n";

    for (auto &kv : sw)
      OS << "BRUNCH_STAT " << kv.first << " " 
         << (kv.second).toSeconds() << "\n";

    for (auto &kv : av)
      OS << "BRUNCH_STAT " << kv.first << " " << kv.second << "\n";

    OS << "************** BRUNCH STATS END ***************** \n";
  }


  void Stopwatch::Print (crab_os &out) const
  {
    long time = getTimeElapsed ();
    long h = time/3600000000L;
    long m = time/60000000L - h*60;
    float s = ((float)time/1000000L) - m*60 - h*3600;

    if (h > 0) out << h << "h";
    if (m > 0) out << m << "m";
    out << s << "s";
  }  
    
  void Averager::Print (crab_os &out) const { out << avg; }
}
#endif 
