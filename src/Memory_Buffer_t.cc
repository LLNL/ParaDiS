#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Memory_Buffer_t.h"

//---------------------------------------------------------------------------------------------------------

#define MIN(a,b) ( (a)<(b) ? (a) : (b) )
#define MAX(a,b) ( (a)>(b) ? (a) : (b) )

//---------------------------------------------------------------------------------------------------------

int       Memory_Buffer_t::mbuf_cnt          = 0;  ///< number of active buffers (diagnostic)

long long Memory_Buffer_t::mbuf_bytes_heap   = 0;  ///< total allocated on host (heap)
long long Memory_Buffer_t::mbuf_bytes_pinned = 0;  ///< total allocated on host (page-locked)
long long Memory_Buffer_t::mbuf_bytes_device = 0;  ///< total allocated on device

long long Memory_Buffer_t::mbuf_hwm_heap     = 0;  ///< high-water mark (heap)
long long Memory_Buffer_t::mbuf_hwm_pinned   = 0;  ///< high-water mark (page-locked)
long long Memory_Buffer_t::mbuf_hwm_device   = 0;  ///< high-water mark (device)

//---------------------------------------------------------------------------------------------------------

static void _init(Memory_Buffer_t &m)
{
   m.mbuf_mode  = MBUF_HEAP;
   m.mbuf_bytes = 0;
   m.mbuf       = 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

Memory_Buffer_t::Memory_Buffer_t(void)                                                           { _init(*this);                                              mbuf_cnt++; }
Memory_Buffer_t::Memory_Buffer_t(             const size_t bytes, const Memory_Buffer_Mode mode) { _init(*this); Allocate(    bytes,mode);                    mbuf_cnt++; }
Memory_Buffer_t::Memory_Buffer_t(const int n, const size_t bytes, const Memory_Buffer_Mode mode) { _init(*this); Allocate(n  ,bytes,mode);                    mbuf_cnt++; }
Memory_Buffer_t::Memory_Buffer_t(const Memory_Buffer_t & m)                                      { _init(*this); Allocate(m.mbuf_bytes,m.mbuf_mode); Copy(m); mbuf_cnt++; }

//---------------------------------------------------------------------------------------------------------------------------------

Memory_Buffer_t::~Memory_Buffer_t() { Recycle(); mbuf_cnt--; }

//---------------------------------------------------------------------------------------------------------------------------------

void  Memory_Buffer_t::Recycle  (void)
{
   if (mbuf && (mbuf_bytes>0))
   {
      if (mbuf_mode==MBUF_HEAP  ) {                          delete [] mbuf  ; mbuf_bytes_heap   -= mbuf_bytes; }
      if (mbuf_mode==MBUF_PINNED) { CUDART_CHECK(cudaFreeHost((void *) mbuf)); mbuf_bytes_pinned -= mbuf_bytes; }
      if (mbuf_mode==MBUF_DEVICE) { CUDART_CHECK(cudaFree    ((void *) mbuf)); mbuf_bytes_device -= mbuf_bytes; }
   }

   mbuf_mode  = MBUF_HEAP;
   mbuf_bytes = 0;
   mbuf       = 0;
}

//---------------------------------------------------------------------------------------------------------------------------------

void Memory_Buffer_t::Allocate 
(
   const size_t             bytes,   ///< allocation size (bytes)
   const Memory_Buffer_Mode mode     ///< allocation location
)
{
   if ((bytes>0) && (bytes>mbuf_bytes) )
   {
      Recycle();

      if (mode==MBUF_HEAP  ) { mbuf = new unsigned char[bytes];                                            mbuf_mode=MBUF_HEAP  ; }
#ifdef GPU_ENABLED
      if (mode==MBUF_PINNED) { CUDART_CHECK(cudaHostAlloc((void **) &mbuf, bytes, cudaHostAllocPortable)); mbuf_mode=MBUF_PINNED; }
      if (mode==MBUF_DEVICE) { CUDART_CHECK(cudaMalloc   ((void **) &mbuf, bytes));                        mbuf_mode=MBUF_DEVICE; }
#else
      if (mode==MBUF_PINNED) { mbuf = new unsigned char[bytes];                                            mbuf_mode=MBUF_HEAP  ; }
#endif

      mbuf_bytes = ( mbuf ? bytes : 0 );

      if (mode==MBUF_HEAP  ) { mbuf_bytes_heap   += mbuf_bytes; mbuf_hwm_heap   = MAX( mbuf_bytes_heap  , mbuf_hwm_heap   ); } 
      if (mode==MBUF_PINNED) { mbuf_bytes_pinned += mbuf_bytes; mbuf_hwm_pinned = MAX( mbuf_bytes_pinned, mbuf_hwm_pinned ); } 
      if (mode==MBUF_DEVICE) { mbuf_bytes_device += mbuf_bytes; mbuf_hwm_device = MAX( mbuf_bytes_device, mbuf_hwm_device ); } 
   }
}

void Memory_Buffer_t::Allocate (const int n, const size_t bytes, const Memory_Buffer_Mode mode) { Allocate( (size_t) n * bytes, mode); }

//---------------------------------------------------------------------------------------------------------------------------------

void Memory_Buffer_t::Clear (void)
{
   if (mbuf && (mbuf_bytes>0))
   {
      if (mbuf_mode==MBUF_HEAP  ) { memset(mbuf,0,mbuf_bytes); }
      if (mbuf_mode==MBUF_PINNED) { memset(mbuf,0,mbuf_bytes); }
      if (mbuf_mode==MBUF_DEVICE) { CUDART_CHECK(cudaMemset((void *) mbuf, 0, mbuf_bytes)); }
   }
}

// Copy()
//
// Implements a synchronous copy between two buffers.  
//
// Note - if the buffers are not the same size, only copies the smaller of the two byte counts.
//---------------------------------------------------------------------------------------------------------------------------------

void Memory_Buffer_t::Copy (const Memory_Buffer_t & m)
{
   if (this!=&m)
   {
      size_t bytes = MIN(mbuf_bytes, m.mbuf_bytes);

      if (mbuf && m.mbuf && (bytes>0))
      {
              if (    ((mbuf_mode==MBUF_HEAP  ) && (m.mbuf_mode==MBUF_HEAP  ))
                   || ((mbuf_mode==MBUF_HEAP  ) && (m.mbuf_mode==MBUF_PINNED))
                   || ((mbuf_mode==MBUF_PINNED) && (m.mbuf_mode==MBUF_HEAP  ))
                   || ((mbuf_mode==MBUF_PINNED) && (m.mbuf_mode==MBUF_PINNED)) ) { memcpy(mbuf,m.mbuf,bytes); }

         else if (    ((mbuf_mode==MBUF_HEAP  ) && (m.mbuf_mode==MBUF_DEVICE))
                   || ((mbuf_mode==MBUF_PINNED) && (m.mbuf_mode==MBUF_DEVICE)) ) { CUDART_CHECK(cudaMemcpy((void *) mbuf,(void *) m.mbuf,bytes,cudaMemcpyDeviceToHost  )); }

         else if (    ((mbuf_mode==MBUF_DEVICE) && (m.mbuf_mode==MBUF_HEAP  ))
                   || ((mbuf_mode==MBUF_DEVICE) && (m.mbuf_mode==MBUF_PINNED)) ) { CUDART_CHECK(cudaMemcpy((void *) mbuf,(void *) m.mbuf,bytes,cudaMemcpyHostToDevice  )); }

         else if (    ((mbuf_mode==MBUF_DEVICE) && (m.mbuf_mode==MBUF_DEVICE)) ) { CUDART_CHECK(cudaMemcpy((void *) mbuf,(void *) m.mbuf,bytes,cudaMemcpyDeviceToDevice)); }
      }
   } 
}

// Copy_Async()
//
// Implements an asynchronous copy between cpu,gpu. Host-to-host transfers remain synchronous.
//---------------------------------------------------------------------------------------------------------------------------------

#ifdef GPU_ENABLED

void Memory_Buffer_t::Copy_Async 
(
   const Memory_Buffer_t & m    ,  ///< source buffer
         cudaStream_t      strm    ///< cuda stream identifier
)
{
   if (this!=&m)
   {
      size_t bytes = MIN(mbuf_bytes, m.mbuf_bytes);

      if (mbuf && m.mbuf && (bytes>0))
      {
         if (    (mbuf_mode==MBUF_HEAP  ) && (m.mbuf_mode==MBUF_HEAP  )
              || (mbuf_mode==MBUF_HEAP  ) && (m.mbuf_mode==MBUF_PINNED)
              || (mbuf_mode==MBUF_PINNED) && (m.mbuf_mode==MBUF_HEAP  )
              || (mbuf_mode==MBUF_PINNED) && (m.mbuf_mode==MBUF_PINNED) ) { memcpy(mbuf,m.mbuf,bytes); }

         if (    (mbuf_mode==MBUF_HEAP  ) && (m.mbuf_mode==MBUF_DEVICE)
              || (mbuf_mode==MBUF_PINNED) && (m.mbuf_mode==MBUF_DEVICE) ) { CUDART_CHECK(cudaMemcpyAsync((void *) mbuf,(void *) m.mbuf,bytes,cudaMemcpyDeviceToHost  )); }

         if (    (mbuf_mode==MBUF_DEVICE) && (m.mbuf_mode==MBUF_HEAP  )
              || (mbuf_mode==MBUF_DEVICE) && (m.mbuf_mode==MBUF_PINNED) ) { CUDART_CHECK(cudaMemcpyAsync((void *) mbuf,(void *) m.mbuf,bytes,cudaMemcpyHostToDevice  )); }

         if (    (mbuf_mode==MBUF_DEVICE) && (m.mbuf_mode==MBUF_DEVICE) ) { CUDART_CHECK(cudaMemcpyAsync((void *) mbuf,(void *) m.mbuf,bytes,cudaMemcpyDeviceToDevice)); }
      }
   } 
}

#endif

//---------------------------------------------------------------------------------------------------------------------------------

