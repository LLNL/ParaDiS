
#include "cuda_portability.h"
#include "GPU_Timer.h"

#ifdef GPU_ENABLED
GPU_Timer::GPU_Timer(void)
{
   tmr_beg=0;
   tmr_end=0;

   cudaEventCreate(&tmr_beg);
   cudaEventCreate(&tmr_end);

   tmr_elapsed_ms = 0.0;
}

GPU_Timer::~GPU_Timer()
{
   if (tmr_beg) { cudaEventDestroy(tmr_beg); }
   if (tmr_end) { cudaEventDestroy(tmr_end); }
}

void GPU_Timer::Start (void)
{
   if (tmr_beg) { cudaEventRecord(tmr_beg,0); }

   tmr_elapsed_ms = 0.0;
}

double GPU_Timer::Stop (void)
{
   if (tmr_end) { cudaEventRecord     (tmr_end,0); 
                  cudaEventSynchronize(tmr_end  ); }

   if (tmr_beg && tmr_end) 
   { 
      float elapsed=0.0;
      cudaEventElapsedTime(&elapsed,tmr_beg,tmr_end); 

      tmr_elapsed_ms = elapsed;
   }

   return(tmr_elapsed_ms);
}
#endif   // GPU_ENABLED

// (stubs when GPU not enabled...)

#ifndef GPU_ENABLED
       GPU_Timer:: GPU_Timer(void) { tmr_elapsed_ms = 0.0;   } 
       GPU_Timer::~GPU_Timer()     {} 
void   GPU_Timer:: Start    (void) { tmr_elapsed_ms = 0.0;   } 
double GPU_Timer:: Stop     (void) { return(tmr_elapsed_ms); }
#endif   // !GPU_ENABLED
