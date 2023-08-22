#pragma once

#ifndef __PDS_MEMORY_BUFFER_H
#define __PDS_MEMORY_BUFFER_H

#include "cuda_portability.h"

typedef enum
{
   MBUF_HEAP    ,     ///< buffer allocated on host (heap)
   MBUF_PINNED  ,     ///< buffer allocated on host (page-locked)
   MBUF_DEVICE        ///< buffer allocated on device
} Memory_Buffer_Mode;

class Memory_Buffer_t
{
   public : 
      Memory_Buffer_Mode   mbuf_mode         ;  ///< location of the allocation
      size_t               mbuf_bytes        ;  ///< size of the allocation (bytes)
      unsigned char       *mbuf              ;  ///< pointer to the allocation

      // static diagnostics...

      static int           mbuf_cnt          ;  ///< number of active buffers (diagnostic)
      static long long     mbuf_bytes_heap   ;  ///< total allocated on host (heap)
      static long long     mbuf_bytes_pinned ;  ///< total allocated on host (page-locked)
      static long long     mbuf_bytes_device ;  ///< total allocated on device

      static long long     mbuf_hwm_heap     ;  ///< high-water mark (heap)
      static long long     mbuf_hwm_pinned   ;  ///< high-water mark (page-locked)
      static long long     mbuf_hwm_device   ;  ///< high-water mark (device)

   public : 
       Memory_Buffer_t(void);
       Memory_Buffer_t(             const size_t bytes, const Memory_Buffer_Mode mode);
       Memory_Buffer_t(const int n, const size_t bytes, const Memory_Buffer_Mode mode);
       Memory_Buffer_t(const Memory_Buffer_t & m);

      ~Memory_Buffer_t();

      void  Recycle  (void);

      void  Allocate (              const size_t bytes, const Memory_Buffer_Mode mode);
      void  Allocate (const int n,  const size_t bytes, const Memory_Buffer_Mode mode);

      void  Clear    (void);

             size_t     Bytes        (void) const { return (mbuf_bytes       ); }

      // access to static diagnostics...

      static int        Count        (void)       { return (mbuf_cnt         ); }
      static long long  Bytes_Heap   (void)       { return (mbuf_bytes_heap  ); }
      static long long  Bytes_Pinned (void)       { return (mbuf_bytes_pinned); }
      static long long  Bytes_Device (void)       { return (mbuf_bytes_device); }

      static long long  HWM_Heap     (void)       { return (mbuf_hwm_heap    ); }
      static long long  HWM_Pinned   (void)       { return (mbuf_hwm_pinned  ); }
      static long long  HWM_Device   (void)       { return (mbuf_hwm_device  ); }

      // casting...

      operator          char   *() const { return ( (         char   *) mbuf ); }
      operator unsigned char   *() const { return ( (unsigned char   *) mbuf ); }
      operator          short  *() const { return ( (         short  *) mbuf ); }
      operator unsigned short  *() const { return ( (unsigned short  *) mbuf ); }
      operator          int    *() const { return ( (         int    *) mbuf ); }
      operator unsigned int    *() const { return ( (unsigned int    *) mbuf ); }
      operator          float  *() const { return ( (         float  *) mbuf ); }
      operator          double *() const { return ( (         double *) mbuf ); }

      void  Copy       (const Memory_Buffer_t & m);
#ifdef GPU_ENABLED
      void  Copy_Async (const Memory_Buffer_t & m, cudaStream_t strm);
#endif
};

#endif
