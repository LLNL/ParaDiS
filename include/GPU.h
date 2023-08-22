#pragma once

#ifndef _PDS_GPU_DEVICE_H
#define _PDS_GPU_DEVICE_H

#include "cuda_portability.h"

#if defined(GPU_ENABLED)

class GPU_Device 
{
   public:
      char             gpu_name[256]       ;   ///< GPU device name
      int              gpu_cnt             ;   ///< total number of GPUs available on this machine/node
      int              gpu_indx            ;   ///< GPU device index
      cudaDeviceProp   gpu_props           ;   ///< GPU device properties
   
      size_t           gpu_global_mem      ;   ///< total global memory available on the device
      size_t           gpu_global_mem_alloc;   ///< total global memory currently allocated
      size_t           gpu_global_mem_free ;   ///< total global memory currently available 

      size_t           gpu_heap_bytes      ;   ///< size of the current stack (bytes)
      size_t           gpu_stack_bytes     ;   ///< size of the current heap  (bytes)

      cudaError_t      gpu_err             ;   ///< most recently set error code
      char            *gpu_err_str         ;   ///< most recently set error text

      FILE            *gpu_err_fd          ;   ///< points to a file descriptor to redirect error messages (default=stdout)
      void           (*gpu_err_cback) (const char *file, const char *func, const int ln, const char *msg, const cudaError_t err);
                                               ///< error callback - if set, invoked whenever errors are checked within the class

   public:
      GPU_Device(void);
      GPU_Device(const int indx);
     ~GPU_Device();

      cudaError_t    Check_Error      (const char *file, const char *func, const int ln);

      cudaError_t    Set_Device_Flags (const int    flags);

      int            Get_Device_Count (void);

      size_t         Get_Heap_Size    (void);
      size_t         Get_Stack_Size   (void);
      size_t         Set_Heap_Size    (const size_t bytes);
      size_t         Set_Stack_Size   (const size_t bytes);

      unsigned char *Malloc           (const size_t bytes);
      unsigned char *Malloc_Host      (const size_t bytes, unsigned int flags=cudaHostAllocPortable);

      cudaError_t    Free             (const unsigned char   *p, const size_t bytes);
      cudaError_t    Free_Host        (const unsigned char   *p);

      cudaError_t    Memset           (const unsigned char   *q, const int val, const size_t bytes, float *msec=0);
      cudaError_t    Memset           (const          void   *q, const int val, const size_t bytes, float *msec=0);

      cudaError_t    Memcpy           (const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);
      cudaError_t    Memcpy           (const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0, float *rate=0);

      cudaError_t    Memcpy_Async     (cudaStream_t strm, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir);
      cudaError_t    Memcpy_Async     (cudaStream_t strm, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir);

      void           Print            (FILE *fd=stdout, const char *title=0);
};

extern GPU_Device     *GPU_New           (const int indx=0);
extern          void   GPU_Delete        (GPU_Device *gpu);

extern          void   GPU_Init          (GPU_Device *gpu, const int indx=0);
extern          int    GPU_Set_Device    (GPU_Device *gpu);
extern unsigned char  *GPU_Malloc        (GPU_Device *gpu, const size_t bytes);
extern          int    GPU_Free          (GPU_Device *gpu, const unsigned char *p, const size_t bytes);

extern          int    GPU_Memset        (GPU_Device *gpu, const unsigned char   *q, const int val, const size_t bytes, float *msec=0);
extern          int    GPU_Memset        (GPU_Device *gpu, const          void   *q, const int val, const size_t bytes, float *msec=0);

extern          int    GPU_Memcpy        (GPU_Device *gpu, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy        (GPU_Device *gpu, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);

extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir);

extern          void   GPU_Print         (GPU_Device *gpu, const char *title);

#endif  // GPU_ENABLED

#endif  // _PDS_GPU_DEVICE_H
