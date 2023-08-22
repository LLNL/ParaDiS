#pragma once

#ifndef _PDS_GPU_TIMER_H
#define _PDS_GPU_TIMER_H

#include "cuda_portability.h"

class GPU_Timer
{
   public :
#ifdef GPU_ENABLED
      cudaEvent_t tmr_beg       ;  ///< points to start cuda event
      cudaEvent_t tmr_end       ;  ///< points to end   cuda event
#endif
      double      tmr_elapsed_ms;  ///< elapsed msecs

   public :
      GPU_Timer(void);
     ~GPU_Timer();

      void    Start (void);
      double  Stop  (void);

      double  Msecs (void) const { return(tmr_elapsed_ms       ); }
      double  Secs  (void) const { return(tmr_elapsed_ms/1000.0); }
};

#endif
