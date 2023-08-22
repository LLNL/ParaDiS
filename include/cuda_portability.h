#pragma once

#ifndef _CUDA_PORTABILITY_H
#define _CUDA_PORTABILITY_H

#include <sys/types.h>

#ifdef GPU_ENABLED
#include <cuda_runtime_api.h>
#endif

//----------------------------------------------------------------------------------------
// __host__, __device__, and __global__ are defined via the runtime API.
// If they haven't been defined at this point, we aren't building with CUDA.
// They are defined below as null macros to enable compilation in non-CUDA environs.
//----------------------------------------------------------------------------------------

#ifdef GPU_ENABLED
#define __cuda_host__     __host__
#define __cuda_device__   __device__
#define __cuda_hdev__     __host__ __device__
#define __cuda_global__   __global__
#define __cuda_constant__ __constant__

#define CUDART_CHECK(err) CUDART_Check(err,__FILE__,__func__,__LINE__);

#else
#define __cuda_host__
#define __cuda_device__
#define __cuda_hdev__
#define __cuda_global__
#define __cuda_constant__

#define CUDART_CHECK(err)

#endif

extern __cuda_host__ int             CUDA_Device_Count      (void);
extern __cuda_host__ void            CUDA_Set_MPI_Device    (void);
extern __cuda_host__ void            CUDA_Set_Device        (const int      indx=0 );

#ifdef GPU_ENABLED
extern __cuda_host__ void            CUDART_Check           (const cudaError_t err , ///< cuda error code
                                                             const char       *file, ///< file string     (e.g. __FILE__)
                                                             const char       *func, ///< function string (e.g. __func__)
                                                             const int         ln    ///< line number     (e.g. __LINE__)
                                                            );

extern __cuda_host__ cudaDeviceProp *CUDA_Device_Properties (const int      indx=0 );
extern __cuda_host__ cudaStream_t   *CUDA_Streams_Malloc    (const int      nstrm  );
extern __cuda_host__ void            CUDA_Streams_Free      (cudaStream_t  *strms, int nstrm);
#endif

extern __cuda_host__ unsigned char  *CUDA_Malloc            (const size_t   bytes);
extern __cuda_host__ unsigned char  *CUDA_Malloc_Host       (const size_t   bytes);
extern __cuda_host__ unsigned char  *CUDA_Malloc_Heap       (const size_t   bytes);

extern __cuda_host__ unsigned char  *CUDA_Realloc           (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes);
extern __cuda_host__ unsigned char  *CUDA_Realloc_Host      (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes);
extern __cuda_host__ unsigned char  *CUDA_Realloc_Heap      (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes);

extern __cuda_host__ void            CUDA_Free              (unsigned char *mbuf);
extern __cuda_host__ void            CUDA_Free_Host         (unsigned char *mbuf);
extern __cuda_host__ void            CUDA_Free_Heap         (unsigned char *mbuf);

// define thread sync support macros for device compiles (only active during device compiles)
//------------------------------------------------------------------------------------------------------------

#ifdef __CUDA_ARCH__
#define cuda_syncthreads() __syncthreads()
#else
#define cuda_syncthreads()
#endif

// atomicAdd
//
// Instruction level support for double precision atomic add's are not available on older hardware.
// This routine implements a software version using a spin loop on the address compare and swap (CAS).
//------------------------------------------------------------------------------------------------------------

#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 600)
static __inline__ __cuda_device__ double atomicAdd(double *addr, double val)
{
   unsigned long long int *addr_as_ull = (unsigned long long int*) addr;
   unsigned long long int  old         = *addr_as_ull, assumed;

   if (val==0.0) { return (__longlong_as_double(old)); }

   do {
      assumed = old;
      old     = atomicCAS(addr_as_ull, assumed, __double_as_longlong(val + __longlong_as_double(assumed)) );
   } while (assumed != old);

   return (__longlong_as_double(old));
}
#endif

#endif  // _CUDA_PORTABILITY_H
