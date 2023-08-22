#pragma once

#ifndef _PDS_GPU_MISC_H
#define _PDS_GPU_MISC_H

#include "cuda_portability.h"

//-------------------------------------------------------------------------------------------
// Externals defined for many GPU/CUDA helper routines...
//-------------------------------------------------------------------------------------------

#if defined(GPU_ENABLED)

extern          int    GPU_Check_Error      (void);

extern          void   GPU_Lock_Device      (const int indx=0);
extern          void   GPU_Unlock_Device    (const int indx=0);

extern          int    GPU_Get_Device_Count (void);
extern          int    GPU_Get_Device_ID    (void);
extern          int    GPU_Get_Device_Props (const int gndx, cudaDeviceProp *props);
extern cudaDeviceProp *GPU_Get_Device_Props (const int gndx);

extern          int    GPU_Set_Device_ID    (const int id);
extern          int    GPU_Set_Device_Flags (const int flags);

extern          char  *GPU_Name             (const int gndx);
extern          int    GPU_Warp_Size        (const int gndx);

extern          int    GPU_Max_Threads      (const int gndx);
extern          int    GPU_Max_Threads      (const int gndx, int &x, int &y, int &z);

extern unsigned char  *GPU_Malloc           (const int gndx, const size_t bytes);
extern unsigned char  *GPU_Malloc           (const int gndx, const int w, const int h, const int d, const int ch, const int hlen);

extern unsigned char  *GPU_Malloc_Host      (                const size_t bytes,                                                  unsigned int flags=cudaHostAllocPortable);
extern unsigned char  *GPU_Malloc_Host      (                const int w, const int h, const int d, const int ch, const int hlen, unsigned int flags=cudaHostAllocPortable);

extern          int    GPU_Free             (const int gndx, const unsigned char *p);
extern          int    GPU_Free_Host        (                const unsigned char *p);

extern          int    GPU_Memset           (const int gndx, const unsigned char *p, const int val, const size_t bytes);

extern          int    GPU_Memcpy           (const int gndx, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);
extern          int    GPU_Memcpy           (const int gndx, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec=0);

extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir);
extern          int    GPU_Memcpy_Async     (cudaStream_t strm, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir);

extern size_t          GPU_Aligned_Bytes    (                                          int sz, int align);
extern size_t          GPU_Aligned_Bytes    (                  int nx,                 int sz, int align);
extern size_t          GPU_Aligned_Bytes    (                  int nx, int ny,         int sz, int align);
extern size_t          GPU_Aligned_Bytes    (                  int nx, int ny, int nz, int sz, int align);
extern unsigned char  *GPU_Aligned_Bytes    (unsigned char *p,                         int sz, int align);
extern unsigned char  *GPU_Aligned_Bytes    (unsigned char *p, int nx,                 int sz, int align);
extern unsigned char  *GPU_Aligned_Bytes    (unsigned char *p, int nx, int ny,         int sz, int align);
extern unsigned char  *GPU_Aligned_Bytes    (unsigned char *p, int nx, int ny, int nz, int sz, int align);

#endif

#endif
