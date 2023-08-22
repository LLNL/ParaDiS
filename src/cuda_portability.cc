#include <stdio.h>
#include <stdlib.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "MPI_Utils.h"

#ifndef GPU_ENABLED
//------------------------------------------------------------------------------------------------------------
// stubs if CUDA not currently active...
//------------------------------------------------------------------------------------------------------------

__cuda_host__ int             CUDA_Device_Count      (void) { return(0); }
__cuda_host__ void            CUDA_Set_MPI_Device    (void) {}
__cuda_host__ void            CUDA_Set_Device        (const int      indx ) {}

__cuda_host__ unsigned char  *CUDA_Malloc            (const size_t   bytes) { return(0); }
__cuda_host__ unsigned char  *CUDA_Malloc_Host       (const size_t   bytes) { return( (bytes>0) ? new unsigned char[bytes] : 0 ); }
__cuda_host__ unsigned char  *CUDA_Malloc_Heap       (const size_t   bytes) { return( (bytes>0) ? new unsigned char[bytes] : 0 ); }

__cuda_host__ unsigned char  *CUDA_Realloc           (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes) { return(0); }

__cuda_host__ void            CUDA_Free              (unsigned char *mbuf) {}
__cuda_host__ void            CUDA_Free_Host         (unsigned char *mbuf) { if (mbuf) delete [] mbuf; }
__cuda_host__ void            CUDA_Free_Heap         (unsigned char *mbuf) { if (mbuf) delete [] mbuf; }

__cuda_host__
unsigned char *CUDA_Realloc_Host (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes)
{
   if ( !mbuf || (bytes>mbuf_bytes) )
   {
      if (mbuf) delete [] mbuf;

      bytes += (bytes/8);

      mbuf       = ( (bytes>0) ? new unsigned char[bytes] : 0 );
      mbuf_bytes = ( mbuf ? bytes : 0 );
   }

   return(mbuf);
}

__cuda_host__
unsigned char *CUDA_Realloc_Heap (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes)
{
   if ( !mbuf || (bytes>mbuf_bytes) )
   {
      if (mbuf) delete [] mbuf;

      bytes += (bytes/8);

      mbuf       = ( (bytes>0) ? new unsigned char[bytes] : 0 );
      mbuf_bytes = ( mbuf ? bytes : 0 );
   }

   return(mbuf);
}

#endif // !GPU_ENABLED

#ifdef GPU_ENABLED

// CUDA_Device_Count()
//
// Returns the total count of GPU's on the system.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
int CUDA_Device_Count (void)
{
   static int gcnt = -1;

   if (gcnt<0)
   {
      CUDART_CHECK(cudaGetDeviceCount(&gcnt));
   }

   return(gcnt);
}

// CUDART_Check()
//
// Checks the cuda error and prints to the console.  Generally called via the CUDART_CHECK() preprocessor
// macro that adds the file, function, and line info to the error.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void CUDART_Check(const cudaError_t err, const char *file, const char *func, const int ln)
{
   if (err!=cudaSuccess)
   {
      printf("%s::%s(ln=%d) cuda error %d (%s)\n", file, func, ln, (int) err, (char *) cudaGetErrorString(err) );
      MPI_Abort(MPI_COMM_WORLD,-1);
      exit(0);
   }
}

// CUDA_Device_Properties()
//
// Returns a pointer to the device properties of the indexed GPU.
// Note - on first call, this routine initializes and loads ALL the device info into a static array.
// Subsequent calls will simply return pointerd to those initialized structures.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
cudaDeviceProp *CUDA_Device_Properties (const int indx)
{
   static cudaDeviceProp *props=0;

   int gcnt = CUDA_Device_Count();

   if (!props && (gcnt>0))
   {
      props = new cudaDeviceProp[gcnt];

      for (int i=0; (i<gcnt); i++)
      { CUDART_CHECK(cudaGetDeviceProperties(props+i,i)); }
   }

   return( props && (0<=indx) && (indx<gcnt) ? (props+indx) : 0 );
}

// CUDA_Set_Device()
//
// Sets the current device.  Note that on multiple GPU systems, the current device is maintained across
// calls to the CUDA runtime. All interactions with the runtime are presumed to be with the most recent
// device set - or the default GPU (0) if the device was never set.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void CUDA_Set_Device (const int indx)
{
   int gcnt = CUDA_Device_Count();

   if ((gcnt>0) && (0<=indx) && (indx<gcnt))
   { CUDART_CHECK(cudaSetDevice(indx)); }
}

// CUDA_Malloc()
//
// Allocates space on the current device.
//------------------------------------------------------------------------------------------------------------

__cuda_host__ unsigned char *CUDA_Malloc (const size_t bytes)
{
   unsigned char *mbuf=0;

   if (bytes>0)
   { CUDART_CHECK(cudaMalloc((void **) &mbuf, bytes)); }

   return(mbuf);
}

// CUDA_Malloc()
//
// Alternate form. Allocates space on the device.
// Note that this version is intended to maintain a buffer that grows over time.  If the buffer already
// exists and can hold the request, nothing is changed. If the buffer hasn't been allocated or is too
// small, the existing allocation is free'd and reallocated with a small fudge to allow for growth.
//------------------------------------------------------------------------------------------------------------

__cuda_host__ unsigned char *CUDA_Realloc (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes)
{
   if ( !mbuf || (bytes>mbuf_bytes) )
   {
      CUDA_Free(mbuf);

      bytes += (bytes/8);

      mbuf       = CUDA_Malloc(bytes);
      mbuf_bytes = ( mbuf ? bytes : 0 );
   }

   return(mbuf);
}

// CUDA_Malloc_Host()
//
// Allocates page-locked (pinned) memory on the host. Pinned memory is special memory on the host
// that can be DMA'd directly between the host/device.  Pinned memory is required for asynchronous
// memory transfers and cuda streaming.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
unsigned char *CUDA_Malloc_Host (const size_t bytes)
{
   unsigned char *mbuf=0;

   if (bytes>0)
   { CUDART_CHECK(cudaHostAlloc((void **) &mbuf, bytes, cudaHostAllocPortable)); }

   return(mbuf);
}

// CUDA_Malloc_Host()
//
// Alternate form. Allocates page-locked (pinned) memory on the
// Note that this version is intended to maintain a buffer that grows over time.  If the buffer already
// exists and can hold the request, nothing is changed. If the buffer hasn't been allocated or is too
// small, the existing allocation is free'd and reallocated with a small fudge to allow for growth.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
unsigned char *CUDA_Realloc_Host (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes)
{
   if ( !mbuf || (bytes>mbuf_bytes) )
   {
      CUDA_Free_Host(mbuf);

      bytes += (bytes/8);

      mbuf       = CUDA_Malloc_Host(bytes);
      mbuf_bytes = ( mbuf ? bytes : 0 );
   }

   return(mbuf);
}

// CUDA_Malloc_Heap()
//
// Allocates heap memory on the host.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
unsigned char *CUDA_Malloc_Heap (const size_t bytes)
{
   unsigned char *mbuf = ( (bytes>0) ? new unsigned char[bytes] : 0 );

   return(mbuf);
}

// CUDA_Malloc_Heap()
//
// Alternate form. Allocates heap memory on the host.
// Note that this version is intended to maintain a buffer that grows over time.  If the buffer already
// exists and can hold the request, nothing is changed. If the buffer hasn't been allocated or is too
// small, the existing allocation is free'd and reallocated with a small fudge to allow for growth.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
unsigned char *CUDA_Realloc_Heap (unsigned char *mbuf, size_t & mbuf_bytes, size_t bytes)
{
   if ( !mbuf || (bytes>mbuf_bytes) )
   {
      CUDA_Free_Heap(mbuf);

      bytes += (bytes/8);

      mbuf       = CUDA_Malloc_Heap(bytes);
      mbuf_bytes = ( mbuf ? bytes : 0 );
   }

   return(mbuf);
}

// CUDA_Free()      - free's space allocated on the GPU
// CUDA_Free_Host() - free's space allocated in page-locked (pinned) memeory on the host
//------------------------------------------------------------------------------------------------------------

__cuda_host__ void  CUDA_Free      (unsigned char *mbuf) { if (mbuf) { CUDART_CHECK(cudaFree    ((void *) mbuf)); } }
__cuda_host__ void  CUDA_Free_Host (unsigned char *mbuf) { if (mbuf) { CUDART_CHECK(cudaFreeHost((void *) mbuf)); } }
__cuda_host__ void  CUDA_Free_Heap (unsigned char *mbuf) { if (mbuf) { delete [] mbuf;                            } }

// CUDA_Streams_Malloc()
//
// Allocates and initializes an array of cuda streams.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
cudaStream_t *CUDA_Streams_Malloc (const int nstrm)
{
   cudaStream_t *strms=0;

   if (nstrm>0)
   {
      strms = new cudaStream_t[nstrm];

      if (strms)
      {
         for (int i=0; (i<nstrm); i++)
         { CUDART_CHECK(cudaStreamCreate(&strms[i])); }
      }
   }

   return(strms);
}

// CUDA_Streams_Free()
//
// Releases streams allocated via CUDA_Streams_Malloc().
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void CUDA_Streams_Free
(
   cudaStream_t   *strms,   ///< array of initialized cuda streams
   int             nstrm    ///< number of streams
)
{
   if (strms && (nstrm>0))
   {
      for (int i=0; (i<nstrm); i++)
      { CUDART_CHECK(cudaStreamDestroy(strms[i])); }

      delete [] strms; strms=0; nstrm=0;
   }
}

// CUDA_Set_MPI_Device()
//
// When GPUs are used along with MPI, only one MPI process can interact with a GPU.
// This routine will set the current cuda device to the index of the MPI process on the current node.
// If the current MPI domain configuration on the deployed node exceeds the number of GPUs, this
// routine will post an error message to the console and abort.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void CUDA_Set_MPI_Device(void)
{
   static int cu_gcnt =  0;   ///< cu_gcnt = cuda device count
   static int cu_gndx = -1;   ///< cu_gndx = cuda device index

   if (cu_gndx<0)
   {
      cu_gcnt = CUDA_Device_Count();

      if (cu_gcnt<1)
      {
         printf("CUDA_Set_MPI_Device() : there are no GPUs available on this node!\n");
         MPI_Abort(MPI_COMM_WORLD,-1);
         exit(0);
      }

      int mpi_node_rank = MPI_Node_Rank();

      if (mpi_node_rank>=cu_gcnt)
      {
         printf("CUDA_Set_MPI_Device() : MPI node rank (%d) exceeds the available GPUs on the host (%d), resubmit with a different MPI domain configuration\n", mpi_node_rank, cu_gcnt);
         MPI_Abort(MPI_COMM_WORLD,-1);
         exit(0);
      }

      cu_gndx = mpi_node_rank;
   }

   // if the gpu index is valid - focus all CUDA runtime calls to the indexed gpu...

   if ( (cu_gcnt>0) && (cu_gndx>=0) ) { CUDA_Set_Device(cu_gndx); }
}

#endif  // GPU_ENABLED
