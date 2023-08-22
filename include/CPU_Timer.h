#pragma once

#ifndef _PDS_CPU_TIMER_H
#define _PDS_CPU_TIMER_H

#include <time.h>

// CPU_Timer
//
// This class implements a simple processing timer.  Note that the class uses the clock()
// function as a basis for timing. Clock() uses the process tick counter rather than the
// more conventional timeval/gettimeofday() interface.
//-------------------------------------------------------------------------------------------

class CPU_Timer
{
   public :
      clock_t  t0;   ///< start   time (process-ticks)
      clock_t  t1;   ///< end     time (process-ticks)
      clock_t  dt;   ///< elapsed time (process-ticks)

   public :
      CPU_Timer(void) { t0 = t1 = dt = 0LL; }
     ~CPU_Timer()     { }

      void Reset  (void) { t0 = t1 = dt = 0LL; }                     // full reset

      void Start  (void) { t0 = clock(); t1 = t0; dt  = 0;       }   // start, resets accumulated time
      void Restart(void) { t0 = clock(); t1 = t0;                }   // start, continues accumulated time
      void Stop   (void) { t1 = clock();          dt += (t1-t0); }   // stop , accumulates time

      double Usecs (void) const { return( 1000000.0*((double) dt)/CLOCKS_PER_SEC ); }
      double Msecs (void) const { return(    1000.0*((double) dt)/CLOCKS_PER_SEC ); }
      double Secs  (void) const { return(           ((double) dt)/CLOCKS_PER_SEC ); }
};

#define CPU_TIMER_START(msec,tmr) if(msec) { *msec=0.0; tmr.Start(); }
#define CPU_TIMER_STOP(msec,tmr)  if(msec) {            tmr.Stop(); *msec=tmr.Msecs(); }

#endif
