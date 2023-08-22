#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#if defined(GPU_ENABLED)

#include <cuda_runtime_api.h>

#include "GPU.h"
#include "GPU_Misc.h"
#include "GPU_Timer.h"

//---------------------------------------------------------------------------------------------------------

void _init(GPU_Device & g)
{
   memset( g.gpu_name ,0,sizeof(g.gpu_name ));
   memset(&g.gpu_props,0,sizeof(g.gpu_props));

   g.gpu_cnt              = 0;
   g.gpu_indx             = 0;
   
   g.gpu_global_mem       = 0;
   g.gpu_global_mem_alloc = 0;
   g.gpu_global_mem_free  = 0;

   g.gpu_heap_bytes       = 0;
   g.gpu_stack_bytes      = 0;

   g.gpu_err              = cudaSuccess;
   g.gpu_err_str          = 0;

   g.gpu_err_fd           = stdout;
   g.gpu_err_cback        = 0;

}

GPU_Device::GPU_Device(void)           
{ 
   _init(*this); 

#ifdef GPU_ENABLED
   CUDART_CHECK(cudaGetDeviceCount(&gpu_cnt));

   if (gpu_cnt>0)
   {
      CUDART_CHECK(cudaGetDeviceProperties(&gpu_props,0));

      gpu_global_mem       = gpu_props.totalGlobalMem;
      gpu_global_mem_alloc = 0;
      gpu_global_mem_free  = gpu_props.totalGlobalMem;

      CUDART_CHECK(cudaDeviceGetLimit(&gpu_stack_bytes, cudaLimitStackSize     ));
      CUDART_CHECK(cudaDeviceGetLimit(&gpu_heap_bytes , cudaLimitMallocHeapSize));
   }
#endif
} 

GPU_Device::GPU_Device(const int indx) 
{ 
   _init(*this); 

#ifdef GPU_ENABLED
   CUDART_CHECK(cudaGetDeviceCount(&gpu_cnt));

   if ( (gpu_cnt>0) && (0<=indx) && (indx<gpu_cnt) )
   {
      gpu_indx = indx;

      CUDART_CHECK(cudaGetDeviceProperties(&gpu_props,gpu_indx));

      gpu_global_mem       = gpu_props.totalGlobalMem;
      gpu_global_mem_alloc = 0;
      gpu_global_mem_free  = gpu_props.totalGlobalMem;

      CUDART_CHECK(cudaSetDevice     (gpu_indx));
      CUDART_CHECK(cudaDeviceGetLimit(&gpu_stack_bytes, cudaLimitStackSize     ));
      CUDART_CHECK(cudaDeviceGetLimit(&gpu_heap_bytes , cudaLimitMallocHeapSize));
   }
#endif
}

//-----------------------------------------------------------------------------------------------------------------------------------------

GPU_Device::~GPU_Device() {} // (nothing currently allocated)

//-----------------------------------------------------------------------------------------------------------------------------------------

cudaError_t GPU_Device::Check_Error (const char *file, const char *func, const int ln) 
{ 
#ifdef GPU_ENABLED
   gpu_err     = cudaGetLastError();
   gpu_err_str = (char *) ( (gpu_err!=cudaSuccess) ? cudaGetErrorString(gpu_err) : 0 );
#endif

   if (gpu_err!=cudaSuccess)
   {
      if (gpu_err_fd   ) fprintf(gpu_err_fd, "%s::%s(ln=%d) err=%d \"%s\"\n", file, func, ln, (int) gpu_err, gpu_err_str );  
      if (gpu_err_cback) (*gpu_err_cback)(file,func,ln,gpu_err_str,gpu_err);
   }

   return(gpu_err);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

size_t GPU_Device::Get_Heap_Size (void)
{
#ifdef GPU_ENABLED
   CUDART_CHECK(cudaSetDevice     (gpu_indx));
   CUDART_CHECK(cudaDeviceGetLimit(&gpu_heap_bytes,cudaLimitMallocHeapSize));
#endif

   return(gpu_heap_bytes);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

size_t GPU_Device::Get_Stack_Size (void)
{
#ifdef GPU_ENABLED
   CUDART_CHECK(cudaSetDevice     (gpu_indx));
   CUDART_CHECK(cudaDeviceGetLimit(&gpu_stack_bytes,cudaLimitStackSize));
#endif

   return(gpu_stack_bytes);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

size_t GPU_Device::Set_Heap_Size (const size_t bytes)
{
#ifdef GPU_ENABLED
   CUDART_CHECK(cudaSetDevice(gpu_indx));
   CUDART_CHECK(cudaDeviceSetLimit(cudaLimitMallocHeapSize, gpu_heap_bytes));
#endif

   return(gpu_heap_bytes);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

size_t GPU_Device::Set_Stack_Size (const size_t bytes)
{
#ifdef GPU_ENABLED
   CUDART_CHECK(cudaSetDevice(gpu_indx));
   CUDART_CHECK(cudaDeviceSetLimit(cudaLimitStackSize, gpu_stack_bytes));
#endif

   return(gpu_stack_bytes);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

cudaError_t GPU_Device::Set_Device_Flags (const int flags)
{  
#ifdef GPU_ENABLED
   CUDART_CHECK(cudaSetDevice(gpu_indx));
   CUDART_CHECK(cudaSetDeviceFlags(flags));
#endif
   
   return(gpu_err);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

int GPU_Device::Get_Device_Count(void) { return(gpu_cnt); }

//-----------------------------------------------------------------------------------------------------------------------------------------

unsigned char *GPU_Device::Malloc (const size_t bytes)
{  
   void *p = 0;
   
   if (bytes>0 && (bytes>gpu_global_mem_free))
      printf("%s::%s(ln=%d) : allocation request exceeds available GPU memory\n", __FILE__, __func__, __LINE__);

   if (bytes>0 && (bytes<gpu_global_mem_free))
   {  
      CUDART_CHECK(cudaSetDevice(gpu_indx));
      CUDART_CHECK(cudaMalloc   (&p,bytes));
      
      if (gpu_err==cudaSuccess)
      {  
         gpu_global_mem_alloc += bytes;
         gpu_global_mem_free  -= bytes;
      }
   }
   
   return ( (unsigned char *) p );
}

//-----------------------------------------------------------------------------------------------------------------------------------------

unsigned char *GPU_Device::Malloc_Host (const size_t bytes, unsigned int flags)
{  
   void *p=0;
   
   if (bytes>0)
   {  
      gpu_err = cudaHostAlloc(&p,bytes,flags);
      
      if (gpu_err != cudaSuccess) { printf("%s::%s(ln=%d) : error allocating host/pinned memory\n", __FILE__, __func__, __LINE__); }
   }
   
   return ( (unsigned char *) p );
}

//-----------------------------------------------------------------------------------------------------------------------------------------

cudaError_t GPU_Device::Free (const unsigned char *p, const size_t bytes)
{
   if (p) 
   {   
      CUDART_CHECK(cudaSetDevice(gpu_indx));
      CUDART_CHECK(cudaFree( (void **) p));

      if (gpu_err==cudaSuccess)
      {   
         gpu_global_mem_alloc -= bytes;
         gpu_global_mem_free  += bytes;
      }   
   }   

   return(gpu_err);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

cudaError_t GPU_Device::Free_Host (const unsigned char *p) 
{
   if (p) { CUDART_CHECK(cudaFreeHost((void *) p)); }   

   return(gpu_err);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

cudaError_t GPU_Device::Memset (const unsigned char *q, const int val, const size_t bytes, float *msec)  { return( Memset((void *) q,val,bytes,msec) ); }
cudaError_t GPU_Device::Memset (const void          *q, const int val, const size_t bytes, float *msec)
{
   if (q && (bytes>0))
   {
      GPU_Timer tmr;

      if (msec) { tmr.Start(); }

      CUDART_CHECK(cudaMemset((void *) q, val, bytes));

      if (msec) { tmr.Stop(); *msec = tmr.Msecs(); }
   }

   return(gpu_err);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

cudaError_t GPU_Device::Memcpy (const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }
cudaError_t GPU_Device::Memcpy (const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate) { return( Memcpy( (void *) q, (void *) p, bytes,dir,msec,rate) ); }

cudaError_t GPU_Device::Memcpy (const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec, float *rate)
{
   if (msec) { *msec = 0.0f; }
   if (rate) { *rate = 0.0f; }

   if (p && q && (bytes>0) )
   {
      GPU_Timer tmr;

      if (msec) { tmr.Start(); }

      CUDART_CHECK(cudaMemcpy((void *) q, (void *) p, bytes, dir));

      if (msec) { tmr.Stop(); *msec = tmr.Msecs(); }
      if (msec && rate && (*msec>0.0)) { *rate = (1000.0*bytes)/(*msec); }
   }

   return(gpu_err);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }
cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir) { return( Memcpy_Async(strm, (void *) q, (void *) p, bytes,dir) ); }

cudaError_t GPU_Device::Memcpy_Async  (cudaStream_t strm, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir)
{
   if (p && q && (bytes>0) )
   {
      CUDART_CHECK(cudaMemcpyAsync( (void *) q, (void *) p, bytes, dir, strm ));
   }

   return(gpu_err);
}

//-----------------------------------------------------------------------------------------------------------------------------------------

void GPU_Device::Print (FILE *fd, const char *title)
{
#ifdef GPU_ENABLED
   cudaDeviceProp *p  = &gpu_props;
   int            *pi = 0;
                            
                            if (title) { fprintf(fd,"gpu_info title=\"%s\"\n", title); }
                            else       { fprintf(fd,"gpu_info\n");                     }
                                         
                                         fprintf(fd,"   name                            = \"%s\" \n" ,           p->name                        );
                                         fprintf(fd,"   gpu_cnt                         = \"%d\" \n" ,           gpu_cnt                        );
                                         fprintf(fd,"   gpu_indx                        = \"%d\" \n" ,           gpu_indx                       );
                                         fprintf(fd,"   gpu_global_mem                  = \"%lu (%luMB)\" \n",   gpu_global_mem           , gpu_global_mem       >>20 );
                                         fprintf(fd,"   gpu_global_mem_alloc            = \"%lu (%luMB)\" \n",   gpu_global_mem_alloc     , gpu_global_mem_alloc >>20 );
                                         fprintf(fd,"   gpu_global_mem_free             = \"%lu (%luMB)\" \n",   gpu_global_mem_free      , gpu_global_mem_free  >>20 );
                                         fprintf(fd,"   shared_mem_per_block            = \"%lu (%luK)\" \n",    p->sharedMemPerBlock     , p->sharedMemPerBlock >>10 );
                                         fprintf(fd,"   regs_per_block                  = \"%d (%dK)\" \n" ,     p->regsPerBlock          , p->regsPerBlock      >>10 );
                                         fprintf(fd,"   warp_size                       = \"%d\" \n" ,           p->warpSize                                          );
                                         fprintf(fd,"   mem_pitch                       = \"%lu (%luMB)\" \n",   p->memPitch              , p->memPitch          >>20 );
                                         fprintf(fd,"   max_threads_per_block           = \"%d\" \n" ,           p->maxThreadsPerBlock          );
      pi=p->maxThreadsDim           ;    fprintf(fd,"   max_threads_dim                 = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxGridSize             ;    fprintf(fd,"   max_grid_size                   = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
                                         fprintf(fd,"   clock_rate                      = \"%d (%d KHz)\" \n" ,  p->clockRate             , p->clockRate         >>10 );
                                         fprintf(fd,"   total_const_mem                 = \"%lu (%luK)\" \n",    p->totalConstMem         , p->totalConstMem     >>10 );
                                         fprintf(fd,"   major                           = \"%d\" \n" ,           p->major                       );
                                         fprintf(fd,"   minor                           = \"%d\" \n" ,           p->minor                       );
                                         fprintf(fd,"   texture_alignment               = \"%lu\" \n",           p->textureAlignment            );
                                         fprintf(fd,"   texture_pitch_alignment         = \"%lu\" \n",           p->texturePitchAlignment       );
                                         fprintf(fd,"   device_overlap                  = \"%d\" \n" ,           p->deviceOverlap               );
                                         fprintf(fd,"   multi_processor_count           = \"%d\" \n" ,           p->multiProcessorCount         );
                                         fprintf(fd,"   kernel_exec_timeout_enabled     = \"%d\" \n" ,           p->kernelExecTimeoutEnabled    );
                                         fprintf(fd,"   integrated                      = \"%d\" \n" ,           p->integrated                  );
                                         fprintf(fd,"   can_map_host_memory             = \"%d\" \n" ,           p->canMapHostMemory            );
                                         fprintf(fd,"   compute_mode                    = \"%d\" \n" ,           p->computeMode                 );
                                         fprintf(fd,"   max_texture_1D                  = \"[%d]\" \n" ,         p->maxTexture1D                );
                                         fprintf(fd,"   max_texture_1D_mipmap           = \"[%d]\" \n" ,         p->maxTexture1DMipmap          );
                                         fprintf(fd,"   max_texture_1D_linear           = \"[%d]\" \n" ,         p->maxTexture1DLinear          );
      pi=p->maxTexture2D            ;    fprintf(fd,"   max_texture_2D                  = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture2DMipmap      ;    fprintf(fd,"   max_texture_2D_mipmap           = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture2DLinear      ;    fprintf(fd,"   max_texture_2D_linear           = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxTexture2DGather      ;    fprintf(fd,"   max_texture_2D_gather           = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture3D            ;    fprintf(fd,"   max_texture_3D                  = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxTexture3DAlt         ;    fprintf(fd,"   max_texture_3D_alt              = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
                                         fprintf(fd,"   max_texture_cube_map            = \"[%d]\" \n" ,         p->maxTextureCubemap           );
      pi=p->maxTexture1DLayered     ;    fprintf(fd,"   max_texture_1D_layered          = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture2DLayered     ;    fprintf(fd,"   max_texture_2D_layered          = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxTextureCubemapLayered;    fprintf(fd,"   max_texture_cube_map_layered    = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
                                         fprintf(fd,"   max_surface_1D                  = \"[%d]\" \n" ,         p->maxSurface1D                );
      pi=p->maxSurface2D            ;    fprintf(fd,"   max_surface_2D                  = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
      pi=p->maxSurface3D            ;    fprintf(fd,"   max_surface_3D                  = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxSurface1DLayered     ;    fprintf(fd,"   max_surface_1D_layered          = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
      pi=p->maxSurface2DLayered     ;    fprintf(fd,"   max_surface_2D_layered          = \"[%d %d %d]\" \n" ,   pi[0], pi[1], pi[2]            );
                                         fprintf(fd,"   max_surface_cube_map            = \"[%d]\" \n" ,         p->maxSurfaceCubemap           );
      pi=p->maxSurfaceCubemapLayered;    fprintf(fd,"   max_surface_cube_map_layered    = \"[%d %d]\" \n" ,      pi[0], pi[1]                   );
                                         fprintf(fd,"   surface_alignment               = \"%lu\" \n",           p->surfaceAlignment            );
                                         fprintf(fd,"   concurrent_kernels              = \"%d\" \n" ,           p->concurrentKernels           );
                                         fprintf(fd,"   ecc_enabled                     = \"%d\" \n" ,           p->ECCEnabled                  );
                                         fprintf(fd,"   pci_bus_id                      = \"%d\" \n" ,           p->pciBusID                    );
                                         fprintf(fd,"   pci_device_id                   = \"%d\" \n" ,           p->pciDeviceID                 );
                                         fprintf(fd,"   pci_domain_id                   = \"%d\" \n" ,           p->pciDomainID                 );
                                         fprintf(fd,"   tcc_driver                      = \"%d\" \n" ,           p->tccDriver                   );
                                         fprintf(fd,"   async_engine_count              = \"%d\" \n" ,           p->asyncEngineCount            );
                                         fprintf(fd,"   unified_addressing              = \"%d\" \n" ,           p->unifiedAddressing           );
                                         fprintf(fd,"   memory_clock_rate               = \"%d (%d KHz)\" \n" ,  p->memoryClockRate    , p->memoryClockRate>>10 );
                                         fprintf(fd,"   memory_bus_width                = \"%d\" \n" ,           p->memoryBusWidth              );
                                         fprintf(fd,"   l2_cache_size                   = \"%d\" \n" ,           p->l2CacheSize                 );
                                         fprintf(fd,"   max_threads_per_multi_processor = \"%d\" \n" ,           p->maxThreadsPerMultiProcessor );
                                         fprintf(fd,"   stream_priorities_supported     = \"%d\" \n" ,           p->streamPrioritiesSupported   );
                                         fprintf(fd,"   global_l1_cache_supported       = \"%d\" \n" ,           p->globalL1CacheSupported      );
                                         fprintf(fd,"   local_l1_cache_supported        = \"%d\" \n" ,           p->localL1CacheSupported       );
                                         fprintf(fd,"   shared_mem_per_multiprocessor   = \"%lu\" \n",           p->sharedMemPerMultiprocessor  );
                                         fprintf(fd,"   regs_per_multiprocessor         = \"%d\" \n" ,           p->regsPerMultiprocessor       );
                                         fprintf(fd,"   managed_memory                  = \"%d\" \n" ,           p->managedMemory               );
                                         fprintf(fd,"   is_multi_gpu_board              = \"%d\" \n" ,           p->isMultiGpuBoard             );
                                         fprintf(fd,"   multi_gpu_board_group_id        = \"%d\" \n" ,           p->multiGpuBoardGroupID        );
#endif // (GPU_ENABLED)

}

//---------------------------------------------------------------------------------------------------------
// GPU_New()
//
// Will create and initialize a GPU device structure. Returns null if no devices are present.
//---------------------------------------------------------------------------------------------------------

GPU_Device *GPU_New (const int indx)
{
   GPU_Device *gpu = 0;
   int         cnt = 0;
   int         err = cudaGetDeviceCount(&cnt);

   if ( (err==cudaSuccess) && (cnt>0) && (indx<cnt) )
   {
      gpu = (GPU_Device *) calloc(1,sizeof(GPU_Device));

      GPU_Init(gpu,indx);
   }

   return(gpu);
}

//---------------------------------------------------------------------------------------------------------
// GPU_Delete()
//
// Will release the memory allocated by a GPU_New(). Upon exit, the pointer is invalid.
//---------------------------------------------------------------------------------------------------------

void GPU_Delete (GPU_Device *gpu)
{
   if (gpu) free(gpu);
}

//---------------------------------------------------------------------------------------------------------
// GPU_Init()
//
// Will initialize the elements of a GPU device structure.
//---------------------------------------------------------------------------------------------------------

void GPU_Init (GPU_Device *gpu, const int indx)
{
   int cnt = GPU_Get_Device_Count();

   if (gpu && (cnt>0) && (indx<cnt))
   {
      GPU_Set_Device_ID(indx);

      gpu->gpu_cnt   = cnt;
      gpu->gpu_indx  = GPU_Get_Device_ID();

      GPU_Get_Device_Props(gpu->gpu_indx,&gpu->gpu_props);

      strncpy(gpu->gpu_name, gpu->gpu_props.name, sizeof(gpu->gpu_name));

      gpu->gpu_global_mem       = gpu->gpu_props.totalGlobalMem;
      gpu->gpu_global_mem_alloc = 0;
      gpu->gpu_global_mem_free  = gpu->gpu_props.totalGlobalMem;
   }
}

//---------------------------------------------------------------------------------------------------------
// GPU_Set_Device()
//
// On platforms with more than one GPU, this routine will set the current GPU device to the device
// indicated by the pointer provided in the argument.
//---------------------------------------------------------------------------------------------------------

int GPU_Set_Device (GPU_Device *gpu)
{
   return ( gpu ? GPU_Set_Device_ID(gpu->gpu_indx) : -1 );
}

//---------------------------------------------------------------------------------------------------------
// GPU_Malloc()
//
// Will allocate gpu device memory on the specified device. Upon success, returns a pointer to the
// allocated memory.  Prints diagnostics and returns null on failure.
//---------------------------------------------------------------------------------------------------------

unsigned char *GPU_Malloc (GPU_Device *gpu, const size_t bytes)
{
   void *p=0;

   if (!gpu) { printf("GPU_Malloc() - cannot allocate memory on NULL GPU device\n"); return(0); }

   if (gpu && (bytes>0) )
   {
      if (bytes>gpu->gpu_global_mem_free)
      {
         printf("GPU_Malloc() - request exceeds available memory on GPU (req=%lu > free=%lu)\n", bytes, gpu->gpu_global_mem_free);
         return(0); 
      }

      int err = GPU_Set_Device(gpu);

      if (!err) { err=cudaMalloc(&p,bytes); }
      if ( err) { printf("GPU_Malloc() - error allocating memory on GPU (req=%lu)\n", bytes); }

      if (p && !err)
      {
         gpu->gpu_global_mem_alloc += bytes;
         gpu->gpu_global_mem_free  -= bytes;
      }
   }

   return ( (unsigned char *) p );
}

//---------------------------------------------------------------------------------------------------------
// GPU_Free()
//
// Releases memory allocated via the above GPU_Malloc() call above.  
// Also updates the state of the GPU device.
//---------------------------------------------------------------------------------------------------------

int GPU_Free (GPU_Device *gpu, const unsigned char *p, const size_t bytes)
{
   int err=0;

   if (!gpu) { printf("GPU_Free() - cannot free memory on NULL GPU device\n"); return(-1); }

   if (gpu && p)
   {
      err=GPU_Set_Device(gpu);

      if (!err) { err=cudaFree( (void *) p); }
      if ( err) { printf("GPU_Free() - error releasing memory from cuda device\n"); }

      if (!err)
      {
         gpu->gpu_global_mem_alloc -= bytes;
         gpu->gpu_global_mem_free  += bytes;
      }
   }

   return(err);
}

//---------------------------------------------------------------------------------------------------------

int GPU_Memset (GPU_Device *gpu, const unsigned char *q, const int val, const size_t bytes, float *msec) { return(GPU_Memset(gpu,(void *) q, val, bytes, msec)); }

int GPU_Memset (GPU_Device *gpu, const void *q, const int val, const size_t bytes, float *msec)
{
   int err=0;

   if (msec) *msec = 0.0f;

   if (!gpu) { printf("GPU_Memset() - cannot initialize memory on NULL GPU device\n"); return(-1); }

   if (gpu && q && (bytes>0) )
   {
      err = GPU_Set_Device(gpu);

      GPU_Timer tmr;

      if (msec) { tmr.Start(); }
      {
         if (!err) { err = cudaMemset( (void *) q, val, bytes); }
         if ( err) { GPU_Check_Error(); printf("GPU_Memset() - error initializing memory\n"); }
      }
      if (msec) { *msec = tmr.Stop(); }
   }

   return(err);
}

// GPU_Memcpy()
//
// Various overloaded GPU memory copy operations. 
// Note that the form of all the calls copy from the source pointer (p) to the destination pointer (q).
// The direction is specified via the cudaMemcpyKind enum. Also note that the direction needs to match
// the actual pointers. Use caution since there's no way to test if the pointers are valid.
//
// A safer alternative to using these memory copy calls is to use the GPU_Memory_Buffer class which
// retains state info as to whether and where a buffer is allocated (e.g. host, device, pinned, etc).
//---------------------------------------------------------------------------------------------------------

int GPU_Memcpy (GPU_Device *gpu, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }
int GPU_Memcpy (GPU_Device *gpu, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec) { return( GPU_Memcpy(gpu, (void *) q, (void *) p, bytes, dir, msec) ); }

int GPU_Memcpy (GPU_Device *gpu, const void *q, const void *p, const size_t bytes, enum cudaMemcpyKind dir, float *msec)
{
   int err=0;

   if (msec) *msec = 0.0f;

   if (!gpu) { printf("GPU_Memcpy() - cannot copy memory to/from NULL GPU device\n"); return(-1); }

   if (gpu && p && q && (bytes>0) )
   {
      err = GPU_Set_Device(gpu);

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

// GPU_Memcpy_Async()
//
// Same as above, but asynchronous. No timing is available because it's meaningless.
//---------------------------------------------------------------------------------------------------------

int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned char   *q, const unsigned char   *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned short  *q, const unsigned short  *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          short  *q, const          short  *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned int    *q, const unsigned int    *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          int    *q, const          int    *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const unsigned long   *q, const unsigned long   *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          long   *q, const          long   *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          float  *q, const          float  *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }
int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          double *q, const          double *p, const size_t bytes, enum cudaMemcpyKind dir) { return(GPU_Memcpy_Async(gpu, strm, (void *) q, (void *) p, bytes, dir)); }

int GPU_Memcpy_Async  (GPU_Device *gpu, cudaStream_t strm, const          void   *q, const          void   *p, const size_t bytes, enum cudaMemcpyKind dir)
{
   int err=0;

   if (!gpu) { printf("%s:%s() - cannot copy memory to/from NULL GPU device\n", __FILE__, __func__); return(-1); }

   if (gpu && p && q && (bytes>0) )
   {
      err = GPU_Set_Device(gpu);

      if (!err) { err = cudaMemcpyAsync( (void *) q, (void *) p, bytes, dir, strm ); }
      if ( err) { GPU_Check_Error(); printf("%s:%s() - error copying memory\n", __FILE__, __func__); }
   }

   return(err);
}

// GPU_Print()
//
// Prints a human-readable summary of the capabilities of a specific GPU.
//---------------------------------------------------------------------------------------------------------

void GPU_Print (GPU_Device *gpu, const char *title)
{
   cudaDeviceProp  *p = ( gpu ? &gpu->gpu_props : 0 );
   int             *pi = 0;
   
   if (p)
   {
                            if (title) { printf("<gpu_info title=\"%s\">\n", title); }
                            else       { printf("<gpu_info>\n");                     }
                                         printf("   <name                            = \"%s\" >\n" ,           p->name                        );
                                         printf("   <gpu_cnt                         = \"%d\" >\n" ,           gpu->gpu_cnt                   );
                                         printf("   <gpu_indx                        = \"%d\" >\n" ,           gpu->gpu_indx                  );
                                         printf("   <gpu_global_mem                  = \"%lu (%luMB)\" >\n",   gpu->gpu_global_mem      , gpu->gpu_global_mem       >>20 );
                                         printf("   <gpu_global_mem_alloc            = \"%lu (%luMB)\" >\n",   gpu->gpu_global_mem_alloc, gpu->gpu_global_mem_alloc >>20 );
                                         printf("   <gpu_global_mem_free             = \"%lu (%luMB)\" >\n",   gpu->gpu_global_mem_free , gpu->gpu_global_mem_free  >>20 );
                                         printf("   <shared_mem_per_block            = \"%lu (%luK)\" >\n",    p->sharedMemPerBlock     , p->sharedMemPerBlock      >>10 );
                                         printf("   <regs_per_block                  = \"%d (%dK)\" >\n" ,     p->regsPerBlock          , p->regsPerBlock           >>10 );
                                         printf("   <warp_size                       = \"%d\" >\n" ,           p->warpSize                                               );
                                         printf("   <mem_pitch                       = \"%lu (%luMB)\" >\n",   p->memPitch              , p->memPitch               >>20 );
                                         printf("   <max_threads_per_block           = \"%d\" >\n" ,           p->maxThreadsPerBlock          );
      pi=p->maxThreadsDim           ;    printf("   <max_threads_dim                 = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxGridSize             ;    printf("   <max_grid_size                   = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
                                         printf("   <clock_rate                      = \"%d (%d KHz)\" >\n" ,  p->clockRate             , p->clockRate              >>10 );
                                         printf("   <total_const_mem                 = \"%lu (%luK)\" >\n",    p->totalConstMem         , p->totalConstMem          >>10 );
                                         printf("   <major                           = \"%d\" >\n" ,           p->major                       );
                                         printf("   <minor                           = \"%d\" >\n" ,           p->minor                       );
                                         printf("   <texture_alignment               = \"%lu\" >\n",           p->textureAlignment            );
                                         printf("   <texture_pitch_alignment         = \"%lu\" >\n",           p->texturePitchAlignment       );
                                         printf("   <device_overlap                  = \"%d\" >\n" ,           p->deviceOverlap               );
                                         printf("   <multi_processor_count           = \"%d\" >\n" ,           p->multiProcessorCount         );
                                         printf("   <kernel_exec_timeout_enabled     = \"%d\" >\n" ,           p->kernelExecTimeoutEnabled    );
                                         printf("   <integrated                      = \"%d\" >\n" ,           p->integrated                  );
                                         printf("   <can_map_host_memory             = \"%d\" >\n" ,           p->canMapHostMemory            );
                                         printf("   <compute_mode                    = \"%d\" >\n" ,           p->computeMode                 );
                                         printf("   <max_texture_1D                  = \"[%d]\" >\n" ,         p->maxTexture1D                );
                                         printf("   <max_texture_1D_mipmap           = \"[%d]\" >\n" ,         p->maxTexture1DMipmap          );
                                         printf("   <max_texture_1D_linear           = \"[%d]\" >\n" ,         p->maxTexture1DLinear          );
      pi=p->maxTexture2D            ;    printf("   <max_texture_2D                  = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture2DMipmap      ;    printf("   <max_texture_2D_mipmap           = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture2DLinear      ;    printf("   <max_texture_2D_linear           = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxTexture2DGather      ;    printf("   <max_texture_2D_gather           = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture3D            ;    printf("   <max_texture_3D                  = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxTexture3DAlt         ;    printf("   <max_texture_3D_alt              = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
                                         printf("   <max_texture_cube_map            = \"[%d]\" >\n" ,         p->maxTextureCubemap           );
      pi=p->maxTexture1DLayered     ;    printf("   <max_texture_1D_layered          = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
      pi=p->maxTexture2DLayered     ;    printf("   <max_texture_2D_layered          = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxTextureCubemapLayered;    printf("   <max_texture_cube_map_layered    = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
                                         printf("   <max_surface_1D                  = \"[%d]\" >\n" ,         p->maxSurface1D                );
      pi=p->maxSurface2D            ;    printf("   <max_surface_2D                  = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
      pi=p->maxSurface3D            ;    printf("   <max_surface_3D                  = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
      pi=p->maxSurface1DLayered     ;    printf("   <max_surface_1D_layered          = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
      pi=p->maxSurface2DLayered     ;    printf("   <max_surface_2D_layered          = \"[%d %d %d]\" >\n" ,   pi[0], pi[1], pi[2]            );
                                         printf("   <max_surface_cube_map            = \"[%d]\" >\n" ,         p->maxSurfaceCubemap           );
      pi=p->maxSurfaceCubemapLayered;    printf("   <max_surface_cube_map_layered    = \"[%d %d]\" >\n" ,      pi[0], pi[1]                   );
                                         printf("   <surface_alignment               = \"%lu\" >\n",           p->surfaceAlignment            );
                                         printf("   <concurrent_kernels              = \"%d\" >\n" ,           p->concurrentKernels           );
                                         printf("   <ecc_enabled                     = \"%d\" >\n" ,           p->ECCEnabled                  );
                                         printf("   <pci_bus_id                      = \"%d\" >\n" ,           p->pciBusID                    );
                                         printf("   <pci_device_id                   = \"%d\" >\n" ,           p->pciDeviceID                 );
                                         printf("   <pci_domain_id                   = \"%d\" >\n" ,           p->pciDomainID                 );
                                         printf("   <tcc_driver                      = \"%d\" >\n" ,           p->tccDriver                   );
                                         printf("   <async_engine_count              = \"%d\" >\n" ,           p->asyncEngineCount            );
                                         printf("   <unified_addressing              = \"%d\" >\n" ,           p->unifiedAddressing           );
                                         printf("   <memory_clock_rate               = \"%d (%d KHz)\" >\n" ,  p->memoryClockRate    , p->memoryClockRate>>10 );
                                         printf("   <memory_bus_width                = \"%d\" >\n" ,           p->memoryBusWidth              );
                                         printf("   <l2_cache_size                   = \"%d\" >\n" ,           p->l2CacheSize                 );
                                         printf("   <max_threads_per_multi_processor = \"%d\" >\n" ,           p->maxThreadsPerMultiProcessor );
                                         printf("   <stream_priorities_supported     = \"%d\" >\n" ,           p->streamPrioritiesSupported   );
                                         printf("   <global_l1_cache_supported       = \"%d\" >\n" ,           p->globalL1CacheSupported      );
                                         printf("   <local_l1_cache_supported        = \"%d\" >\n" ,           p->localL1CacheSupported       );
                                         printf("   <shared_mem_per_multiprocessor   = \"%lu\" >\n",           p->sharedMemPerMultiprocessor  );
                                         printf("   <regs_per_multiprocessor         = \"%d\" >\n" ,           p->regsPerMultiprocessor       );
                                         printf("   <managed_memory                  = \"%d\" >\n" ,           p->managedMemory               );
                                         printf("   <is_multi_gpu_board              = \"%d\" >\n" ,           p->isMultiGpuBoard             );
                                         printf("   <multi_gpu_board_group_id        = \"%d\" >\n" ,           p->multiGpuBoardGroupID        );
                                         printf("</gpu_info>\n");
   }
}

#endif

