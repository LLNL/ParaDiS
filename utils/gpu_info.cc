#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <cuda_runtime_api.h>

#include "Args.h"
#include "CPU_Timer.h"

int Gbytes(unsigned long bytes) { unsigned long gb = 1024*1024*1024UL; return ( (int) ( (bytes/gb) + ( (bytes%gb) ? 1 : 0 ) ) ); }
int Mbytes(unsigned long bytes) { unsigned long mb =      1024*1024UL; return ( (int) ( (bytes/mb) + ( (bytes%mb) ? 1 : 0 ) ) ); }
int Kbytes(unsigned long bytes) { unsigned long kb =           1024UL; return ( (int) ( (bytes/kb) + ( (bytes%kb) ? 1 : 0 ) ) ); }

// GPU_Count()
//
// Returns the number of currently available GPUs
//-------------------------------------------------------------------------------------------------

static int GPU_Count (void)
{
   int cnt=0;

   cudaGetDeviceCount(&cnt);

   return (cnt);
}

// GPU_Print()
//
// Will dump some basic info on an indexed GPU.
//-------------------------------------------------------------------------------------------------

static void GPU_Print (const int gndx)
{
   cudaError_t    err = cudaSetDevice(gndx);
   cudaDeviceProp props;

   if (err) { printf("%s::%s(ln=%d) : error setting cuda device id=%d\n", __FILE__, __func__, __LINE__, gndx); return; }

   err = cudaGetDeviceProperties(&props,gndx);

   if (err) { printf("%s::%s(ln=%d) : error retrieving GPU device properties\n", __FILE__, __func__, __LINE__); return; }

   else
   {
      cudaDeviceProp *p  = &props;
      int            *pi = 0;

                                         printf("GPU #%d %s :\n", gndx, p->name );
                                         printf("   name                            = %s \n"               , p->name                                                               );
                                         printf("   global memory                   = %dMB (%dGB)\n"       , Mbytes(p->totalGlobalMem), Gbytes(p->totalGlobalMem)                  );
                                         printf("   shared_mem_per_block            = %lu (%dK) \n"        , p->sharedMemPerBlock     , Kbytes(p->sharedMemPerBlock)               );
                                         printf("   regs_per_block                  = %d (%dK) \n"         , p->regsPerBlock          , Kbytes(p->regsPerBlock     )               );
                                         printf("   warp_size                       = %d \n"               , p->warpSize                                                           );
                                         printf("   mem_pitch                       = %lu (%dGB) \n"       , p->memPitch              , Gbytes(p->memPitch         )               );
                                         printf("   max_threads_per_block           = %d \n"               , p->maxThreadsPerBlock                                                 );
      pi=p->maxThreadsDim           ;    printf("   max_threads_dim                 = [%d x %d x %d] \n"   , pi[0], pi[1], pi[2]                                                   );
      pi=p->maxGridSize             ;    printf("   max_grid_size                   = [%dM x %dK x %dK] \n", Mbytes(pi[0]), Kbytes(pi[1]), Kbytes(pi[2])                           );
                                         printf("   clock_rate                      = %d (%d KHz) \n"      , p->clockRate             , Kbytes(p->clockRate        )               );
                                         printf("   total_const_mem                 = %lu (%dK) \n"        , p->totalConstMem         , Kbytes(p->totalConstMem    )               );
                                         printf("   major/minor                     = %d.%d \n"            , p->major, p->minor                                                    );
                                         printf("   texture_alignment               = %lu \n"              , p->textureAlignment                                                   );
                                         printf("   texture_pitch_alignment         = %lu \n"              , p->texturePitchAlignment                                              );
                                         printf("   device_overlap                  = %s \n"               , p->deviceOverlap             ? "yes" : "no"                           );
                                         printf("   multi_processor_count           = %d \n"               , p->multiProcessorCount                                                );
                                         printf("   kernel_exec_timeout_enabled     = %s \n"               , p->kernelExecTimeoutEnabled  ? "yes" : "no"                           );
                                         printf("   integrated                      = %s \n"               , p->integrated                ? "yes" : "no"                           );
                                         printf("   can_map_host_memory             = %s \n"               , p->canMapHostMemory          ? "yes" : "no"                           );
                                         printf("   compute_mode                    = %d \n"               , p->computeMode                                                        );
                                         printf("   max_texture_1D                  = [%dK] \n"            , Kbytes(p->maxTexture1D)                                               );
                                         printf("   max_texture_1D_mipmap           = [%dK] \n"            , Kbytes(p->maxTexture1DMipmap)                                         );
                                         printf("   max_texture_1D_linear           = [%dM] \n"            , Mbytes(p->maxTexture1DLinear)                                         );
      pi=p->maxTexture2D            ;    printf("   max_texture_2D                  = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
      pi=p->maxTexture2DMipmap      ;    printf("   max_texture_2D_mipmap           = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
      pi=p->maxTexture2DLinear      ;    printf("   max_texture_2D_linear           = [%dK x %dK x %dK] \n", Kbytes(pi[0]), Kbytes(pi[1]), Kbytes(pi[2])                           );
      pi=p->maxTexture2DGather      ;    printf("   max_texture_2D_gather           = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
      pi=p->maxTexture3D            ;    printf("   max_texture_3D                  = [%dK x %dK x %dK] \n", Kbytes(pi[0]), Kbytes(pi[1]), Kbytes(pi[2])                           );
      pi=p->maxTexture3DAlt         ;    printf("   max_texture_3D_alt              = [%dK x %dK x %dK] \n", Kbytes(pi[0]), Kbytes(pi[1]), Kbytes(pi[2])                           );
                                         printf("   max_texture_cube_map            = [%dK] \n"            , Kbytes(p->maxTextureCubemap)                                          );
      pi=p->maxTexture1DLayered     ;    printf("   max_texture_1D_layered          = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
      pi=p->maxTexture2DLayered     ;    printf("   max_texture_2D_layered          = [%dK x %dK x %dK] \n", Kbytes(pi[0]), Kbytes(pi[1]), Kbytes(pi[2])                           );
      pi=p->maxTextureCubemapLayered;    printf("   max_texture_cube_map_layered    = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
                                         printf("   max_surface_1D                  = [%dK] \n"            , Kbytes(p->maxSurface1D)                                               );
      pi=p->maxSurface2D            ;    printf("   max_surface_2D                  = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
      pi=p->maxSurface3D            ;    printf("   max_surface_3D                  = [%dK x %dK x %dK] \n", Kbytes(pi[0]), Kbytes(pi[1]), Kbytes(pi[2])                           );
      pi=p->maxSurface1DLayered     ;    printf("   max_surface_1D_layered          = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
      pi=p->maxSurface2DLayered     ;    printf("   max_surface_2D_layered          = [%dK x %dK x %dK] \n", Kbytes(pi[0]), Kbytes(pi[1]), Kbytes(pi[2])                           );
                                         printf("   max_surface_cube_map            = [%dK] \n"            , Kbytes(p->maxSurfaceCubemap)                                          );
      pi=p->maxSurfaceCubemapLayered;    printf("   max_surface_cube_map_layered    = [%dK x %dK] \n"      , Kbytes(pi[0]), Kbytes(pi[1])                                          );
                                         printf("   surface_alignment               = %lu \n"              , p->surfaceAlignment                                                   );
                                         printf("   concurrent_kernels              = %s \n"               , p->concurrentKernels         ? "yes" : "no"                           );
                                         printf("   ecc_enabled                     = %s \n"               , p->ECCEnabled                ? "yes" : "no"                           );
                                         printf("   pci_bus_id                      = %d \n"               , p->pciBusID                                                           );
                                         printf("   pci_device_id                   = %d \n"               , p->pciDeviceID                                                        );
                                         printf("   pci_domain_id                   = %d \n"               , p->pciDomainID                                                        );
                                         printf("   tcc_driver                      = %d \n"               , p->tccDriver                                                          );
                                         printf("   async_engine_count              = %d \n"               , p->asyncEngineCount                                                   );
                                         printf("   unified_addressing              = %s \n"               , p->unifiedAddressing         ? "yes" : "no"                           );
                                         printf("   memory_clock_rate               = %d (%d KHz) \n"      , p->memoryClockRate           , Kbytes(p->memoryClockRate)             );
                                         printf("   memory_bus_width                = %d \n"               , p->memoryBusWidth                                                     );
                                         printf("   l2_cache_size                   = %d (%dK)\n"          , p->l2CacheSize               , Kbytes(p->l2CacheSize    )             );
                                         printf("   max_threads_per_multi_processor = %d \n"               , p->maxThreadsPerMultiProcessor                                        );
                                         printf("   stream_priorities_supported     = %s \n"               , p->streamPrioritiesSupported ? "yes" : "no"                           );
                                         printf("   global_l1_cache_supported       = %s \n"               , p->globalL1CacheSupported    ? "yes" : "no"                           );
                                         printf("   local_l1_cache_supported        = %s \n"               , p->localL1CacheSupported     ? "yes" : "no"                           );
                                         printf("   shared_mem_per_multiprocessor   = %lu (%dK)\n"         , p->sharedMemPerMultiprocessor, Kbytes(p->sharedMemPerMultiprocessor)  );
                                         printf("   regs_per_multiprocessor         = %d (%dK)\n"          , p->regsPerMultiprocessor     , Kbytes(p->regsPerMultiprocessor     )  );
                                         printf("   managed_memory                  = %s \n"               , p->managedMemory             ? "yes" : "no"                           );
                                         printf("   is_multi_gpu_board              = %s \n"               , p->isMultiGpuBoard           ? "yes" : "no"                           );
                                         printf("   multi_gpu_board_group_id        = %d \n"               , p->multiGpuBoardGroupID                                               );
                                         printf("\n");
   }

}

//-------------------------------------------------------------------------------------------------

static void GPU_Rate_Test
(
   unsigned char    *dst      ,
   unsigned char    *src      ,
   unsigned long     bytes    ,
   unsigned long     bytes_max,
   unsigned long     bsize    ,
   enum cudaMemcpyKind dir    ,
            double & secs     ,
            double & rate )
{
   secs = rate = 0.0;

   CPU_Timer tmr;

   tmr.Start();
   {
      unsigned char *p = src, *pe = src+bytes;
      unsigned char *q = dst, *qe = dst+bytes;

      int n = bytes_max/bsize;

      for (int i=0; (i<n); i++)
      {
         cudaMemcpy( (void *) q, (void *) p, bsize, dir);

         p += bsize; p = ( (p<pe) ? p : src );
         q += bsize; q = ( (q<qe) ? q : dst );
      }
   }
   tmr.Stop();

   secs = tmr.Secs();
   rate = ( (secs>0.0) ? ((double) bytes_max)/(secs*1.0e6) : 0.0 );
}

//-------------------------------------------------------------------------------------------------

static void GPU_Rate_Test (void)
{
   unsigned long mbytes    = 1024*1024;
   unsigned long bytes     =     256*mbytes;  // (256 MB)
   unsigned long bytes_max =  2*1024*mbytes;  // (  2 GB)

   unsigned char *cpu_buf = (unsigned char *) malloc(bytes);
   unsigned char *pin_buf = 0; int err = cudaHostAlloc((void **) &pin_buf, bytes, cudaHostAllocPortable);
   unsigned char *gpu_buf = 0;     err = cudaMalloc   ((void **) &gpu_buf, bytes);

   if (cpu_buf && gpu_buf && pin_buf)
   {
      printf("\n");
      printf("#                      data transfer rates (secs, MB/S)             \n");
      printf("#  bsize (K)    cpu->gpu      gpu->cpu      pin->gpu      gpu->pin  \n");
      printf("# ----------  ------------  ------------  ------------  ------------\n");

      for (unsigned long bsize=1024; (bsize<=bytes); bsize*=2 )
      {
         printf(" %11d", (int)(bsize>>10));

         for (int j=0; (j<4); j++)
         {
            unsigned char *src=0, *dst=0; enum cudaMemcpyKind dir;

            if (j==0) { dst=gpu_buf; src=cpu_buf; dir=cudaMemcpyHostToDevice; }
            if (j==1) { dst=cpu_buf; src=gpu_buf; dir=cudaMemcpyDeviceToHost; }
            if (j==2) { dst=gpu_buf; src=pin_buf; dir=cudaMemcpyHostToDevice; }
            if (j==3) { dst=pin_buf; src=gpu_buf; dir=cudaMemcpyDeviceToHost; }

            double secs=0.0, rate=0.0;
            GPU_Rate_Test(dst,src,bytes,bytes_max,bsize,dir,secs,rate);

            printf(" %7.3lf %5d", secs, (int) rate);
         }

         printf("\n");
      }
   }

   if (cpu_buf) { free        (         cpu_buf);  cpu_buf=0; }
   if (pin_buf) { cudaFreeHost((void *) pin_buf);  pin_buf=0; }
   if (gpu_buf) { cudaFree    ((void *) gpu_buf);  gpu_buf=0; }
}

//-------------------------------------------------------------------------------------------------

int main(int argc, char **argv)
{
   int gmem = Get_Flag ("-memtst", argc, argv);

   int gcnt = GPU_Count();

   if (gcnt>0)
   {
      for (int i=0; (i<gcnt); i++) { GPU_Print(i); }
   }
   else
   { printf("No GPUs are detected on this platform\n"); }

   if ( (gcnt>0) && gmem)
      GPU_Rate_Test();

   return(0);
   exit(0);
}

