#include <cstdio>
#include "CpuTimer.hpp"

#include "CubitMessage.hpp"

// Added by Cat for WIN32 port
#ifndef WIN32
  #include <sys/param.h>
  #include <sys/times.h>
#endif


#include <time.h>
#ifdef SOLARIS
#include <unistd.h>
#endif

#ifndef WIN32
#ifndef HZ
#ifdef CLK_TCK
#define HZ CLK_TCK
#else
#define HZ 60
#endif
#endif
#else
#define HZ CLOCKS_PER_SEC
#endif



// Added by Cat for WIN32 port
#ifdef WIN32
struct tms {
	   clock_t tms_utime;		/* user time */
	   clock_t tms_stime;		/* system time */
	   clock_t tms_cutime;		/* user time, children */
	   clock_t tms_cstime;		/* system time, children */
      };

  #include <windows.h>

#endif


CpuTimer::CpuTimer()
{
  tms current;

  // Added by Cat for WIN32 port
#ifdef WIN32
  nt_times(&current);
  wallTimeInitial = clock();
  wallTime = wallTimeInitial;
#else
  gettimeofday(&wallTimeInitial, NULL);
  wallTime = wallTimeInitial;
  times( &current );
#endif
  
  // Store current value of cpu time
  cpu = current.tms_utime +
        current.tms_stime +
        current.tms_cutime +
        current.tms_cstime;
  cpuInitial = cpu;
}

// Return time values

double
CpuTimer::cpu_secs()
{
  tms current;

  // Added by Cat for WIN32 port
#ifdef WIN32
  nt_times(&current);    
#else
  times( &current );
#endif
  
  time_t cpu_now = current.tms_utime +
    current.tms_stime +
    current.tms_cutime +
    current.tms_cstime;
  time_t delta = cpu_now - cpu;
  cpu   = cpu_now;
#ifdef WIN32
  if( delta == 0 )
     delta = 1;
#endif
  return (double) delta / HZ;
}

double
CpuTimer::clock_secs()
{
  double elapsed_time;
#ifdef WIN32
  clock_t current;
  current = clock();
  elapsed_time = double( current - wallTime );
  elapsed_time = elapsed_time / HZ;

  wallTime = current;
#else
  timeval current;
  gettimeofday( &current, NULL );
  double prev_time = wallTime.tv_sec;
  prev_time += double( wallTime.tv_usec ) / 1000000.0;
  
  elapsed_time = current.tv_sec;
  elapsed_time += double( current.tv_usec ) / 1000000.0;
  elapsed_time = elapsed_time - prev_time;
  
  wallTime = current;
#endif
  
  return elapsed_time;
}

double
CpuTimer::elapsed(bool wall_time)
{
  if( wall_time )
  {
    double elapsed;
#ifdef WIN32
    clock_t current;
    current = clock();
    elapsed = double( current - wallTimeInitial);
    elapsed = elapsed / HZ;
#else    
    double current_time;
    timeval current;
    gettimeofday( &current, NULL );
    double initial_time;
    initial_time = wallTimeInitial.tv_sec;
    initial_time += double(wallTimeInitial.tv_usec) / 1000000.0;
    current_time = current.tv_sec;
    current_time += double( current.tv_usec ) / 1000000.0;
    
    elapsed = current_time - initial_time;
#endif
    return elapsed;
  }
  else
  {
    tms current;

// Added by Cat for WIN32 port
#ifdef WIN32
    nt_times(&current);    
#else
    times( &current );
#endif
// Store totals
  
    time_t cpu_now = current.tms_utime +
        current.tms_stime +
        current.tms_cutime +
        current.tms_cstime;
    time_t elapsed = cpu_now - cpuInitial;
    return (double) elapsed / HZ;
  }
}


// Added by Cat for WIN32 port
#ifdef WIN32

void CpuTimer::nt_times(tms *sys)
{  
  HANDLE current_process = GetCurrentProcess();

  FILETIME ftCreation, ftExit, ftKernel, ftUser;  
  GetProcessTimes(current_process, &ftCreation, &ftExit, &ftKernel, &ftUser);  

  //kernel and user times are amounts of time rather than an actual time period. 
  //The value in the FILETIME structure is expressed in 100-nanosecond units.
  //Divide by 10 million to convert to seconds.

  __int64 i64 = *((__int64 *) &ftKernel);  
  double dw_kernel = (i64 / 10000000.0);

  i64 = *((__int64 *) &ftUser);
  double dw_user = (i64 / 10000000.0);

	sys->tms_utime = dw_user*HZ;
	sys->tms_stime = dw_kernel*HZ;
	sys->tms_cutime = 0;		
	sys->tms_cstime = 0;
        //PRINT_DEBUG( CUBIT_DEBUG_3,"clock returned %d\n", sys->tms_stime );  
}
#endif
