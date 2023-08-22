//---------------------------------------------------------------------------------------------------------
// Implements a number of GPU-related helper functions with simple error reporting...
//---------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <stdlib.h>

#if defined(GPU_ENABLED)

#include <cuda_runtime_api.h>

#include "GPU_Misc.h"
#include "GPU_Timer.h"

// GPU statics...
//
// There are a number of GPU parameters that need to be periodically interrogated.  Given that these
// parameters remain unchanged because most are physical properties of the devices, there's only 
// a need to initialize them once at startup.
//---------------------------------------------------------------------------------------------------------

static int             gpu_cnt   = 0;       // will be initialized to the number of available GPUs
static cudaDeviceProp *gpu_props = 0;       // will be an array or device properties

//---------------------------------------------------------------------------------------------------------

int GPU_Check_Error (void) 
{ 
   cudaError_t  err =         cudaGetLastError();
   const char  *p   = ( err ? cudaGetErrorString(err) : 0 );

   if (p)
   {
      printf( "GPU subsystem error - \"%s\"\n", p); 
   }

   return( (int) err );
}

int GPU_Get_Device_Count (void)
{
   static int init=0;

   if (!init)
   {
      int err = cudaGetDeviceCount(&gpu_cnt);

      if (err) { GPU_Check_Error(); printf("GPU_Get_Device_Count() - error obtaining cuda device count\n"); }

      gpu_cnt = ( err ? 0 : gpu_cnt );
      init    = 1;
   }

   return(gpu_cnt);
}

int GPU_Get_Device_Props (const int gndx, cudaDeviceProp *props)
{
   int cnt = 0;
   cudaError_t err = cudaGetDeviceCount(&cnt);
 
   if (err) { printf("GPU_Get_Device_Props() = error obtaining cuda device count\n"); }

   if ( (cnt>0) && (0<=gndx) && (gndx<cnt) && props )
   {
       err = cudaGetDeviceProperties(props,gndx);
       if (err) { printf("GPU_Get_Device_Props() - error retrieving GPU device properties\n"); }
   }

   return( (int) err);
}

cudaDeviceProp *GPU_Get_Device_Props (const int gndx)
{
   if (!gpu_props)
   {
      int n = GPU_Get_Device_Count();

      gpu_props = (cudaDeviceProp *) ( (n>0) ? calloc(n,sizeof(cudaDeviceProp)) : 0 );

      if (gpu_props)
      {
         for (int i=0; (i<n); i++)
         {
            int err = cudaGetDeviceProperties(gpu_props+i,i);

            if (err) { GPU_Check_Error(); printf("GPU_Get_Device_Props() - error retrieving GPU device properties\n"); }
         }
      }
   }

   return( gpu_props && (0<=gndx) && (gndx<gpu_cnt) ? (gpu_props+gndx) : 0 );
}

int GPU_Get_Device_ID (void)
{
   int id  = 0;
   int err = cudaGetDevice(&id);

   if (err) { GPU_Check_Error(); printf("GPU_Get_Device_ID() - error obtaining cuda device id\n"); }

   id = ( err ? 0 : id );

   return(id);
}

int GPU_Set_Device_ID (const int id)
{
   int err=0;
   int cnt=GPU_Get_Device_Count();

   if (cnt<1)             { printf("GPU_Set_Device_ID() - no GPUs are available\n"      ); return(-1); } 
   if ((id<0)||(cnt<=id)) { printf("GPU_Set_Device_ID() - requested GPU not available\n"); return(-1); }

   if (id!=GPU_Get_Device_ID())
   {
      err = cudaSetDevice(id);
   
      if (err) { printf("GPU_Set_Device_ID() - error setting cuda device id\n"); }
   } 

   return(err);
}

int GPU_Set_Device_Flags (const int flags)
{
   int err = cudaSetDeviceFlags(flags);

   if (err) { GPU_Check_Error(); printf("GPU_Set_Device_Flags() - error setting cuda device flags\n"); }

   return(err);
}

//---------------------------------------------------------------------------------------------------------

char *GPU_Name (const int gndx)
{
   cudaDeviceProp *props = GPU_Get_Device_Props(gndx);

   return ( props ? props->name : 0 );
}

int GPU_Warp_Size (const int gndx)
{
   cudaDeviceProp *props = GPU_Get_Device_Props(gndx);

   return ( props ? props->warpSize : 0 );
}

int GPU_Max_Threads (const int gndx)
{
   cudaDeviceProp *props = GPU_Get_Device_Props(gndx);

   return ( props ? props->maxThreadsPerBlock : 0 );
}

int GPU_Max_Threads (const int gndx, int &x, int &y, int &z)
{
   x=y=z=1;

   cudaDeviceProp *props = GPU_Get_Device_Props(gndx);

   if (props)
   {
      x = props->maxThreadsDim[0];
      y = props->maxThreadsDim[1];
      z = props->maxThreadsDim[2];
   }

   return ( props ? props->maxThreadsPerBlock : 0 );
}

// GPU_Get_Thread_Args()
//
// Given a GPU index and requested thread count, will return a block and 
// thread count appropriate for the requested GPU.
//---------------------------------------------------------------------------------------------------------

void GPU_Get_Thread_Args (const int gndx, int &nblks, int &nthrds, const int n)
{
   int nwarp = GPU_Warp_Size  (gndx);
   int nmax  = GPU_Max_Threads(gndx);

   if ( (n>0) && (nwarp>0) && (nmax>0) )
   {
      nthrds = ( (n<nwarp) ? nwarp : ( (n<nmax) ? n : nmax ) );  // nthrds = constrain to valid range for the requested gpu
      nthrds = (nwarp*(nthrds-1))/nwarp;                         // nthrds = multiple of warp size
      nblks  = (n    +(nthrds-1))/nthrds;                        // nblks  = number of blocks to launch
   }
}

//---------------------------------------------------------------------------------------------------------

size_t GPU_Aligned_Bytes(                        int sz, int align) { return(GPU_Aligned_Bytes( 1, 1,1,sz,align)); }
size_t GPU_Aligned_Bytes(int nx,                 int sz, int align) { return(GPU_Aligned_Bytes(nx, 1,1,sz,align)); }
size_t GPU_Aligned_Bytes(int nx, int ny,         int sz, int align) { return(GPU_Aligned_Bytes(nx,ny,1,sz,align)); }
size_t GPU_Aligned_Bytes(int nx, int ny, int nz, int sz, int align)
{
   size_t bytes  = (size_t) nx * (size_t) ny * (size_t) nz * (size_t) ((sz<1) ? 1 : sz);
          bytes += ( (align>0) && (bytes%align) ? (align-(bytes%align)) : 0 );

   return(bytes);
}

//---------------------------------------------------------------------------------------------------------

unsigned char *GPU_Aligned_Bytes(unsigned char *p,                         int sz, int align) { return(GPU_Aligned_Bytes(p, 1, 1,1,sz,align)); }
unsigned char *GPU_Aligned_Bytes(unsigned char *p, int nx,                 int sz, int align) { return(GPU_Aligned_Bytes(p,nx, 1,1,sz,align)); }
unsigned char *GPU_Aligned_Bytes(unsigned char *p, int nx, int ny,         int sz, int align) { return(GPU_Aligned_Bytes(p,nx,ny,1,sz,align)); }
unsigned char *GPU_Aligned_Bytes(unsigned char *p, int nx, int ny, int nz, int sz, int align)
{
   size_t bytes  = (size_t) nx * (size_t) ny * (size_t) nz * (size_t) ((sz<1) ? 1 : sz);
          bytes += ( (align>0) && (bytes%align) ? (align-(bytes%align)) : 0 );

   return( p ? p+bytes : 0 );
}

//---------------------------------------------------------------------------------------------------------

unsigned char *GPU_Malloc (const int gndx, const size_t bytes)
{
   void *p = 0;

   if (bytes>0)
   {
      int err=GPU_Set_Device_ID(gndx);

      if (!err) { err=cudaMalloc(&p,bytes); }
      if ( err) { GPU_Check_Error(); printf("GPU_Malloc() - error allocating memory on gpu\n"); }
   }

   return ( (unsigned char *) p );
}

unsigned char *GPU_Malloc (const int gndx, const int w, const int h, const int d, const int ch, const int hlen)
{
   size_t bytes = (size_t)(h*(w*d*ch/8)+hlen); 
   return( GPU_Malloc(gndx,bytes) );
}

//---------------------------------------------------------------------------------------------------------

unsigned char *GPU_Malloc_Host (const size_t bytes, unsigned int flags)
{
   void *p=0;

   if (bytes>0)
   {
      int err = cudaHostAlloc(&p,bytes,flags); 
   
      if (err) { GPU_Check_Error(); printf("GPU_Malloc_Host() - error allocating memory on host\n"); }
   }

   return ( (unsigned char *) p );
}

unsigned char *GPU_Malloc_Host (const int w, const int h, const int d, const int ch, const int hlen, unsigned int flags)
{
   size_t bytes = (size_t)(h*(w*d*ch/8)+hlen); 
   return( GPU_Malloc_Host(bytes,flags) );
}

//---------------------------------------------------------------------------------------------------------

int GPU_Free (const int gndx, const unsigned char *p)
{
   int err=0;

   if (p)
   {
      err=GPU_Set_Device_ID(gndx);

      if (!err) { err=cudaFree( (void **) p); }
      if ( err) { printf("GPU_Free() - error releasing memory from cuda device\n"); }
   }

   return(err);
}

int GPU_Free_Host (const unsigned char *p)
{
   int err = ( p ? (int) cudaFreeHost((void *) p) : 0 );

   if (err) { printf("GPU_Free_Host() - error releasing host/pinned memory\n"); }

   return(err);
}

//---------------------------------------------------------------------------------------------------------

int GPU_Memset (const int gndx, const unsigned char *p, const int val, const size_t bytes)
{
   int err=0;

   if (p && (bytes>0) )
   {
      err = GPU_Set_Device_ID(gndx);

      if (!err) { err = cudaMemset( (void *) p,val,bytes); }
      if ( err) { printf("GPU_Memset() - error setting memory on cuda device id\n"); }
   }

   return(err);
}

//---------------------------------------------------------------------------------------------------------

int GPU_Memcpy (const int gndx, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }
int GPU_Memcpy (const int gndx, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gndx, (void *) q, (void *) p, bytes, dir, msec ) ); }

int GPU_Memcpy (const int gndx, const void *q, const void *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec)
{
   int err=0;

   if (msec) *msec = 0.0f;

   if (p && q && (bytes>0) )
   {
      err = GPU_Set_Device_ID(gndx);

      GPU_Timer tmr;

      if (msec) { tmr.Start(); } 
      {
         if (!err) { err = cudaMemcpy( (void *) q, (void *) p, bytes, dir); }
         if ( err) { GPU_Check_Error(); printf("GPU_Memcpy() - error copying memory\n"); }
      }
      if (msec) { *msec = tmr.Stop(); }
   }

   return(err);
}

//---------------------------------------------------------------------------------------------------------

int GPU_Memcpy_Async (cudaStream_t strm, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }
int GPU_Memcpy_Async (cudaStream_t strm, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir) { return( GPU_Memcpy_Async(strm,(void *) q, (void *) p, bytes, dir) ); }

int GPU_Memcpy_Async (cudaStream_t strm, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir)
{
   int err=0;

   if (p && q && (bytes>0) )
   {
      err = cudaMemcpyAsync((void *) q, (void *) p, bytes, dir, strm);

      if (err) { printf("GPU_Memcpy_Async() - error copying memory\n"); }
   }

   return(err);
}

#endif
