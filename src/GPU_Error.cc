#include <string.h>

#include "cuda_portability.h"
#include "GPU_Error.h"

static int  *gpu_errv_cpu=0;

#ifdef GPU_ENABLED
static int  *gpu_errv_gpu=0;
#endif

__cuda_constant__ int *gpu_errv_addr = 0;

__cuda_host__
void GPU_Error_Vector::Allocate(const int n)
{
   size_t bytes = (2+2*n)*sizeof(int);

#ifdef GPU_ENABLED
   if (!gpu_errv_cpu || !gpu_errv_gpu )
   {
      if (gpu_errv_cpu) { CUDART_CHECK(cudaFreeHost((void *) gpu_errv_cpu));  gpu_errv_cpu=0; }
      if (gpu_errv_gpu) { CUDART_CHECK(cudaFree    ((void *) gpu_errv_gpu));  gpu_errv_gpu=0; }
   
      CUDART_CHECK(cudaHostAlloc((void **) &gpu_errv_cpu, bytes, cudaHostAllocPortable));
      CUDART_CHECK(cudaMalloc   ((void **) &gpu_errv_gpu, bytes                       ));
   }
#else
   if (!gpu_errv_cpu) { gpu_errv_cpu = new int[2*n]; }
#endif

   if (gpu_errv_cpu) 
   { 
      memset(gpu_errv_cpu,0,bytes);
      gpu_errv_cpu[0]=n; 
      gpu_errv_cpu[1]=0; 
   }

#ifdef GPU_ENABLED
   if (gpu_errv_cpu && gpu_errv_gpu && (bytes>0)) 
      CUDART_CHECK( cudaMemcpy((void *) gpu_errv_gpu, (void *) gpu_errv_cpu, bytes, cudaMemcpyHostToDevice) );

   if (gpu_errv_gpu)
      CUDART_CHECK( cudaMemcpyToSymbol( gpu_errv_addr, gpu_errv_gpu, sizeof(int *)) );
#endif
}

__cuda_host__
void GPU_Error_Vector::Reset (void)
{
   const int    n     = ( gpu_errv_cpu ? gpu_errv_cpu[0] : 0 );
   const size_t bytes = (2+2*n)*sizeof(int);

   if (gpu_errv_cpu) 
   { 
      memset(gpu_errv_cpu,0,bytes);
      gpu_errv_cpu[0]=n; 
      gpu_errv_cpu[1]=0; 
   }

#ifdef GPU_ENABLED
   if (gpu_errv_cpu && gpu_errv_gpu && (bytes>0)) 
      CUDART_CHECK(cudaMemcpy((void *) gpu_errv_gpu, (void *) gpu_errv_cpu, bytes, cudaMemcpyHostToDevice) );
#endif
}

// Get_Errors()
//
// Will retrieve the error vector from the GPU.
// 
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void GPU_Error_Vector::Get_Errors (void)
{
#ifdef GPU_ENABLED
   const int    n     = ( gpu_errv_cpu ? gpu_errv_cpu[0] : 0 );
   const size_t bytes = (2+2*n)*sizeof(int);

   if (gpu_errv_cpu && gpu_errv_gpu && (bytes>0)) 
      CUDART_CHECK(cudaMemcpy((void *) gpu_errv_cpu, (void *) gpu_errv_gpu, bytes, cudaMemcpyDeviceToHost) );
#endif
}

// Set_Error()
//
// Will append the current thread index and code to the current gpu error vector.
// Note that for code executing on the device, we need to pick up the location 
// of the static error vector from the address supplied in constant memory.
//------------------------------------------------------------------------------------------------------------

__cuda_hdev__
void GPU_Error_Vector::Set_Error (const int tndx, const int code)
{
   int *errv=0;

#ifdef GPU_ENABLED
#ifdef __CUDA_ARCH__
   errv = gpu_errv_addr; // (device only)
#else
   errv = gpu_errv_cpu;  // (host only)
#endif
#endif

   if (errv)
   {
      int  pmax = errv[0];
      int  pcnt = errv[1];

      if (pcnt<pmax)
      {
         int *pvec = errv + 2+2*pcnt;

         pvec[0] = tndx;
         pvec[1] = code;
         errv[1]++;
      }
   }
}
