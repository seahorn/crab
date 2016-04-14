#ifndef CRAB_STATS__HPP_
#define CRAB_STATS__HPP_

/* Code from SeaHorn */

#include <map>
#include <sys/time.h>
#include <sys/resource.h>
#include <iostream>
#include <iosfwd>
#include <crab/config.h>

namespace crab
{

#ifdef HAVE_STATS
  class Stopwatch
  {
  private:
    long started;
    long finished;
    long timeElapsed;

    long systemTime () const
    {
      struct rusage ru;
      getrusage (RUSAGE_SELF, &ru);
      long r = ru.ru_utime.tv_sec * 1000000L + ru.ru_utime.tv_usec;
      return r;

    }

  public:
    Stopwatch () { start (); }

    void start ()
    {
      started = systemTime ();
      finished = -1;
      timeElapsed = 0;
    }

    void stop ()
    {
      if (finished < started)
	{
	  finished = systemTime ();
	}

    }

    void resume ()
    {
      if (finished >= started)
	{
	  timeElapsed += finished - started;
	  started = systemTime ();
	  finished = -1;
	}
    }

    long getTimeElapsed () const
    {
      if (finished < started) return timeElapsed + systemTime () - started;
      return timeElapsed + finished - started;
    }


    void Print (std::ostream &out) const;

    double toSeconds(){
      double time = ((double) getTimeElapsed () / 1000000) ;
      return time;
    }

  };


  /** Computes running average */
  class Averager
  {
  private:
    size_t count;
    double avg;

  public :
    Averager () : count(0), avg (0) {}

    double add (double k)
    {
      avg += (k - avg)/(++count);
      return avg;
    }

    void Print (std::ostream &out) const;

  };
#else
  struct Stopwatch
  {
    Stopwatch () {}
    void start () {}
    void stop () {}
    void resume () {}
    long getTimeElapsed () const { return 0;}
    void Print (std::ostream &out) {}
    double toSeconds(){ return 0.0;}

  };
  struct Averager
  {
  public :
   Averager () {}
   double add (double k) { return 0.0;}
   void Print (std::ostream &out) const {}
  };

#endif 
}



namespace crab
{

#ifdef HAVE_STATS
  inline std::ostream &operator<< (std::ostream &OS, const Stopwatch &sw)
  {
    sw.Print (OS);
    return OS;
  }

  inline std::ostream &operator<< (std::ostream &OS, const Averager &av)
  {
    av.Print (OS);
    return OS;
  }

  class CrabStats
  {
  private:
    static std::map<std::string,unsigned> counters;
    static std::map<std::string,Stopwatch> sw;
    static std::map<std::string,Averager> av;
    static std::map<std::string,std::string> ss;

  public:
    static unsigned  get (const std::string &n);
    static double avg (const std::string &n, double v);
    static unsigned uset (const std::string &n, unsigned v);

    static void sset (const std::string &n, std::string v);
    static std::string& sget (const std::string &n);
    
    static void count (const std::string &name);
    static void count_max (const std::string &name, unsigned v);

    static void start (const std::string &name);
    static void stop (const std::string &name);
    static void resume (const std::string &name);

    /** Outputs all statistics to std output */
    static void Print (std::ostream &OS);
    static void PrintBrunch (std::ostream &OS);
  };



  /**
      Usage: add
         auto_timer X("foo.bar");
      at the beginning of a block (i.e., function body, conditional, etc)
      to measure the time it takes to execute.
   */
  // class auto_timer
  // {
  // private:
  //   std::string n;
  // public:
  //   auto_timer (const std::string &name) : n(name) { CrabStats::resume (n); }
  //   ~auto_timer () { CrabStats::stop (n); }
  // };

  template <typename Output>
  class TimeIt
  {
    const char* m_msg;
    Output &m_out;
    Stopwatch m_sw;
    double m_min;
    
  public:
    TimeIt (const char *msg, Output out, double min = 0.0) :
      m_msg (msg), m_out (out), m_min (min) {}
    ~TimeIt () 
    {
      m_sw.stop ();
      if (m_sw.toSeconds () >= m_min)
        m_out << "TimeIt: " << m_msg << " " << m_sw << "\n";
    }
    
  };
  
  class ScopedCrabStats 
  {
    std::string m_name;
  public:
    ScopedCrabStats (const std::string &name, bool reset = false) : m_name(name) 
    { 
      if (reset) 
        { 
          m_name += ".last";
          CrabStats::start (m_name);
        }
      else
        CrabStats::resume (m_name); 
    }
    ~ScopedCrabStats () { CrabStats::stop (m_name); }
  };  
#else
  inline std::ostream &operator<< (std::ostream &OS, const Stopwatch &sw){}
  inline std::ostream &operator<< (std::ostream &OS, const Averager &av){}
  struct CrabStats
  {
    static unsigned  get (const std::string &n){}
    static double avg (const std::string &n, double v){}
    static unsigned uset (const std::string &n, unsigned v){}
    static void sset (const std::string &n, std::string v){}
    static std::string& sget (const std::string &n){}
    static void count (const std::string &name){}
    static void count_max (const std::string &name, unsigned v){}
    static void start (const std::string &name){}
    static void stop (const std::string &name){}
    static void resume (const std::string &name){}
    static void Print (std::ostream &OS){}
    static void PrintBrunch (std::ostream &OS){}
  };
  template <typename Output>
  struct TimeIt
  {
    TimeIt (const char *msg, Output out, double min = 0.0) {}
  };
  struct ScopedCrabStats 
  {
    ScopedCrabStats (const std::string &name, bool reset = false) {}
  };  
#endif 
}

#define CRAB_MEASURE_FN crab::ScopedCrabStats __stats__(__FUNCTION__)
#define CRAB_MEASURE_FN_LAST crab::ScopedCrabStats __stats_last__(__FUNCTION__, true)


#endif
