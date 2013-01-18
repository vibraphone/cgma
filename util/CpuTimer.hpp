//- Class: CpuTimer

#if !defined(CPU_TIMER)
#define CPU_TIMER
#include <sys/types.h>
#ifndef WIN32
#include <sys/time.h>
#else
#include <time.h>
#endif
#include "CubitUtilConfigure.h"


class CUBIT_UTIL_EXPORT CpuTimer {
public:
  CpuTimer();			//- initialise to current system time
  double cpu_secs();		        //- return CPU time in seconds since last
                                //- call to cpu_secs();
  double clock_secs();          //- return wall clock time in seconds since last
                                //- call to clock_secs();
  double elapsed(bool wall_time = false);    //- return CPU time in seconds since 'birth' if 
                                             //- wall_time is false, else returns the wall clock time.
  
private:
  time_t cpu;			//- cpu time in 1/HZ units since last call
                                //- to cpu_secs()
  time_t cpuInitial;             //- cpu time in 1/HZ units since construction.


    // Added by Cat for NT port
  #ifdef WIN32
  void  nt_times(struct tms *);
  clock_t wallTimeInitial;
  clock_t wallTime;
  #else
  timeval wallTimeInitial;     //- Time at construction.
  timeval wallTime;            //- Time since last call.
  #endif
};

#endif

