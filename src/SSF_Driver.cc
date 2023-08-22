//------------------------------------------------------------------------------------------------------------
//
// This module is the main driver for computing gpu-enabled isotropic and anisotropic segment-to-segment forces.
// This module consists of mostly driver routines and kernel modules that invoke the respective host/device
// routines provided in the companion modules...
//
//    SSF_Iso.cc                : computes isotropic   forces for a segment pair (both parallel and non-parallel)
//    SSF_Aniso_NonParallel.cc  : computes anisotropic forces for a non-parallel segment pair
//    SSF_Aniso_Parallel.cc     : computes anisotropic forces for a parallel     segment pair
//
// Short description of routines within this module...
//
//    SSF_Allocate              : allocates buffers on host/device for computing forces on gpu
//    SSF_Allocate_Nodes        : allocates buffers to retain node structures on host/device
//    SSF_Free                  : free's all host/device buffers allocated within this module
//    SSF_Initialize_Iso        : initializes static constants for isotropic forces (isotropic only)
//    SSF_Initialize_Aniso      : initializes static constants and tables for both isotropic and anisotropic forces
//    SSF_Initialize            : initializes static constants and tables using Home data structure
//    SSF_Update_Constants      : updates the static constants
//    SSF_Nodes_Serialize       : serializes an array of nodes that can be used on both host/device
//    SSF_Aniso                 : computes anisotropic forces between two segments (host/device)
//    SSF_Iso                   : computes isotropic   forces between two segments (host/device)
//    SSF_Aniso_CPU             : computes anisotropic forces for an array of segment pairs (host only)
//    SSF_Iso_CPU               : computes isotropic   forces for an array of segment pairs (host only)
//    SSF_Aniso_K               : kernel callback for computing anisotropic forces between two segments (device only)
//    SSF_Iso_K                 : kernel callback for computing isotropic   forces between two segments (device only)
//    SSF_Aniso_Reset_Forces_K  : kernel callback for resetting forces in an array of nodes (device only)
//    SSF_Aniso_Force_Sum_K     : kernel callback for accumulating forces in an array of nodes (device only)
//    SSF_Aniso_Force_Reset_CPU : resets      forces in an array of nodes (host only)
//    SSF_Aniso_Force_Sum_CPU   : accumulates forces in an array of nodes (host only)
//
//    SSF_GPU                   : main driver to compute segment forces on gpu
//    SSF_GPU_Strm              : main driver to compute segment forces on gpu (streams version)
//    SSF_NV_Send               :    sends a node vector to the device     (for streams version)
//    SSF_PV_Send               :    sends a pair vector to the device     (for streams version)
//    SSF_FV_Recv               :    receives computed segment forces      (for streams version)
//    SSF_Launch                :    launches the force kernel computation (for streams version)
//
//    SSF_Gather_Nodes          : creates and indexes an array of nodes suitable for use on device
//    SSF_Send_Nodes            : sends    an array of nodes to   the device
//    SSF_Recv_Nodes            : receives an array of nodes from the device
//    SSF_Gather_Pairs          : serializes the pair vector suitable for use on gpu
//    SSF_Scatter_Forces        : accumulates forces in an array of segments
//
//    SSF_Compute_Forces        : main entry point for computing gpu-enabled segment-to-segment forces
//                              :    drives both isotropic and anisotropic forces
//------------------------------------------------------------------------------------------------------------

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>

#include "mpi_portability.h"
#include "cuda_portability.h"
#include "omp_portability.h"

#include "Typedefs.h"
#include "Complex_t.h"
#include "V3.h"
#include "PBC.h"
#include "Anisotropic.h"
#include "SSF_PV_t.h"
#include "SSF_FV_t.h"
#include "SSF_Iso.h"
#include "SSF_Aniso.h"
#include "SSF_Driver.h"

#ifdef ANISO_QMAX
#define QMAX ANISO_QMAX
#else
#define QMAX 20
#endif

//------------------------------------------------------------------------------------------------------------

#define VDELETE(a)   if (a)      { delete [] a; a=0; }
#define VZERO(a,n)   if (a)      { for (int i=0; (i<(n)); ++i) { (a)[i]=0.0;    } }
#define VCOPY(a,b,n) if (a && b) { for (int i=0; (i<(n)); ++i) { (a)[i]=(b)[i]; } }

#define CONSTRAIN(a,n,m) ( ((a)<(n)) ? (n) : ( ((a)<=(m)) ? (a) : (m) ) )

#ifndef MAX
#define MAX(a,b)         ( ((a)>(b)) ? (a) : (b) )
#endif

// local statics...
//------------------------------------------------------------------------------------------------------------

static int            ssf_init        = 0;        ///< are the local statics initialized?

static real8          a               = 0.0;      ///< core radius
static real8          mu              = 0.0;      ///< shear modulus
static real8          nu              = 0.0;      ///< poisson ratio
static real8          ecrit           = 0.0;      ///< critical angle for parallelism test

static int            pbx             = 0;        ///< pbc active (x) (1=yes,0=no)
static int            pby             = 0;        ///< pbc active (y) (1=yes,0=no)
static int            pbz             = 0;        ///< pbc active (z) (1=yes,0=no)

static real8          lx              = 0.0;      ///< simulation box size (x)
static real8          ly              = 0.0;      ///< simulation box size (y)
static real8          lz              = 0.0;      ///< simulation box size (z)

static real8          sx              = 0.0;      ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
static real8          sy              = 0.0;      ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
static real8          sz              = 0.0;      ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)

static int            qmax            = 0;        ///< spherical harmonics expansion factor (aniso only)
static int            mtype           = 0;        ///< material type index                  (aniso only)

static real8         *c66c            = 0;        ///< elastic constants matrix (6x6)            (host)
static complex8      *c12c            = 0;        ///< complex elements of aniso rotation matrix (host)
static real8         *e3c             = 0;        ///< real    elements of aniso rotation matrix (host)
static real8         *fqrc            = 0;        ///< spherical harmonic <Fq> table (real)      (host)
static real8         *fqic            = 0;        ///< spherical harmonic <Fq> table (imag)      (host)

static real8         *c66g            = 0;        ///< elastic constants matrix (6x6)            (device)
static complex8      *c12g            = 0;        ///< complex elements of aniso rotation matrix (device)
static real8         *e3g             = 0;        ///< real    elements of aniso rotation matrix (device)
static real8         *fqrg            = 0;        ///< spherical harmonic <Fq> table (real)      (device)
static real8         *fqig            = 0;        ///< spherical harmonic <Fq> table (imag)      (device)

static size_t         mbuf_cpu_bytes  = 0;        ///< current size of the allocated host   buffer
static size_t         mbuf_gpu_bytes  = 0;        ///< current size of the allocated device buffer
static unsigned char *mbuf_cpu        = 0;        ///< allocated host   buffer (pinned memory)
static unsigned char *mbuf_gpu        = 0;        ///< allocated device buffer

static size_t         nodes_cpu_bytes = 0;        ///< size of the serialized nodes (host)
static size_t         nodes_gpu_bytes = 0;        ///< size of the serialized nodes (device)
static unsigned char *nodes_cpu       = 0;        ///< array of serialized nodes    (host)
static unsigned char *nodes_gpu       = 0;        ///< array of serialized nodes    (device)

static int            cu_warp         = 32;       ///< number of cuda threads per warp
static int            cu_tmax         = 1024;     ///< max number of cuda threads per warp
static int            cu_nthrd        = 128;      ///< number of cuda threads per kernel launch
static int            cu_npblk        = 64*1024;  ///< number of segment pairs per cuda stream block

#ifdef __CUDACC__
static int             cu_nstrm       = 0;        ///< number of allocated cuda streams
static cudaStream_t   *cu_strms       = 0;        ///< array of allocated cuda streams
static cudaDeviceProp *cu_props       = 0;        ///< points to cuda device properties structure
#endif

// expose the local statics....
//------------------------------------------------------------------------------------------------------------

__cuda_host__ int            SSF_Initialized     (void) { return(::ssf_init); }

__cuda_host__ size_t         SSF_Mbuf_CPU_Bytes  (void) { return(::mbuf_cpu_bytes); }
__cuda_host__ size_t         SSF_Mbuf_GPU_Bytes  (void) { return(::mbuf_gpu_bytes); }
__cuda_host__ unsigned char *SSF_Mbuf_CPU        (void) { return(::mbuf_cpu      ); }
__cuda_host__ unsigned char *SSF_Mbuf_GPU        (void) { return(::mbuf_gpu      ); }

__cuda_host__ real8          SSF_A               (void) { return(::a    ); }
__cuda_host__ real8          SSF_MU              (void) { return(::mu   ); }
__cuda_host__ real8          SSF_NU              (void) { return(::nu   ); }
__cuda_host__ int            SSF_Qmax            (void) { return(::qmax ); }
__cuda_host__ int            SSF_Mtype           (void) { return(::mtype); }
__cuda_host__ real8          SSF_Ecrit           (void) { return(::ecrit); }

__cuda_host__ int            SSF_PBX             (void) { return(::pbx  ); }
__cuda_host__ int            SSF_PBY             (void) { return(::pby  ); }
__cuda_host__ int            SSF_PBZ             (void) { return(::pbz  ); }

__cuda_host__ real8          SSF_LX              (void) { return(::lx   ); }
__cuda_host__ real8          SSF_LY              (void) { return(::ly   ); }
__cuda_host__ real8          SSF_LZ              (void) { return(::lz   ); }

__cuda_host__ real8          SSF_SX              (void) { return(::sx   ); }
__cuda_host__ real8          SSF_SY              (void) { return(::sy   ); }
__cuda_host__ real8          SSF_SZ              (void) { return(::sz   ); }

__cuda_host__ real8         *SSF_Aniso_C66C      (void) { return(::c66c ); }
__cuda_host__ complex8      *SSF_Aniso_C12C      (void) { return(::c12c ); }
__cuda_host__ real8         *SSF_Aniso_E3C       (void) { return(::e3c  ); }
__cuda_host__ real8         *SSF_Aniso_FQRC      (void) { return(::fqrc ); }
__cuda_host__ real8         *SSF_Aniso_FQIC      (void) { return(::fqic ); }

__cuda_host__ real8         *SSF_Aniso_C66G      (void) { return(::c66g ); }
__cuda_host__ complex8      *SSF_Aniso_C12G      (void) { return(::c12g ); }
__cuda_host__ real8         *SSF_Aniso_E3G       (void) { return(::e3g  ); }
__cuda_host__ real8         *SSF_Aniso_FQRG      (void) { return(::fqrg ); }
__cuda_host__ real8         *SSF_Aniso_FQIG      (void) { return(::fqig ); }

__cuda_host__ size_t         SSF_Nodes_CPU_Bytes (void) { return(::nodes_cpu_bytes); }
__cuda_host__ size_t         SSF_Nodes_GPU_Bytes (void) { return(::nodes_gpu_bytes); }
__cuda_host__ unsigned char *SSF_Nodes_CPU       (void) { return(::nodes_cpu      ); }
__cuda_host__ unsigned char *SSF_Nodes_GPU       (void) { return(::nodes_gpu      ); }

// SSF_Allocate()
//
// Allocates the dynamic elements of the GPU-enabled version.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Allocate
(
   const size_t bytes   ///< host and device buffer sizes
)
{
    mbuf_cpu = CUDA_Realloc_Host(mbuf_cpu, mbuf_cpu_bytes, bytes);
    mbuf_gpu = CUDA_Realloc     (mbuf_gpu, mbuf_gpu_bytes, bytes);
}

__cuda_host__
void SSF_Allocate
(
   const int nn,   ///< number of nodes
   const int np    ///< number of segment pairs
)
{
    size_t bytes  =   np*     sizeof(SSF_FV_t);                      // fv (resulting forces  ) (f1(xyz),f2(xyz),f3(xyz),f4(xyz),...)
           bytes += 3*nn*     sizeof(real8   );                      // nv (node positions    ) (n1(xyz),n2(xyz),n3(xyz),n4(xyz),...)
           bytes +=   np* MAX(sizeof(SSF_PV_t), 6*3*sizeof(real8));  // pv (pair pv structures) (pv1,pv2,pv3,pv4,pv5,pv6,...)

    mbuf_cpu = CUDA_Realloc_Host(mbuf_cpu, mbuf_cpu_bytes, bytes);
    mbuf_gpu = CUDA_Realloc     (mbuf_gpu, mbuf_gpu_bytes, bytes);
}

// SSF_Allocate_Nodes
//
// Allocates the dynamic buffers needed for the serialized nodes.
// The space required is a function of the total node and arm count.
//
// Note that all the routines that use the reallocation API add fudges
// to each request to reduce the delete and allocation cycles.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Allocate_Nodes
(
         Node_t **narr     ,    ///< array  of node pointers (may be sparse)
   const int      narr_cnt      ///< number of node pointers
)
{
   int ncnt = 0;  // ncnt = total node count
   int acnt = 0;  // ncnt = total arm  count

   if (narr)
   {
      for (int i=0; (i<narr_cnt); i++)
      { if (narr[i]) { ncnt++; acnt+=narr[i]->numNbrs; } }
   }

   // reserve space for an array of node pointers and serialized nodes...

   size_t bytes  =   ncnt*sizeof(SSF_Node_t  *);   // reserve space for node pointers
          bytes +=   ncnt*sizeof(SSF_Node_t   );   // reserve space for node structures
          bytes +=   acnt*sizeof(unsigned long);   // reserve space for arm tags
          bytes += 3*acnt*sizeof(real8        );   // reserve space for arm burgers vectors
          bytes += 3*acnt*sizeof(real8        );   // reserve space for arm forces

   nodes_cpu = CUDA_Realloc_Host(nodes_cpu, nodes_cpu_bytes, bytes);   // allocate (host)
   nodes_gpu = CUDA_Realloc     (nodes_gpu, nodes_gpu_bytes, bytes);   // allocate (device)
}

__cuda_host__
void SSF_Free (void)
{
    // free the data exchange buffers...

    CUDA_Free_Host(mbuf_cpu );  mbuf_cpu =0; mbuf_cpu_bytes =0;
    CUDA_Free     (mbuf_gpu );  mbuf_gpu =0; mbuf_gpu_bytes =0;

    CUDA_Free_Host(nodes_cpu);  nodes_cpu=0; nodes_cpu_bytes=0;
    CUDA_Free     (nodes_gpu);  nodes_gpu=0; nodes_gpu_bytes=0;

    // also free static elements...

    VDELETE(c66c);  c66c = 0;
    VDELETE(c12c);  c12c = 0;
    VDELETE(e3c );  e3c  = 0;
    VDELETE(fqrc);  fqrc = 0;
    VDELETE(fqic);  fqic = 0;

    CUDA_Free((unsigned char *) c66g);  c66g = 0;
    CUDA_Free((unsigned char *) c12g);  c12g = 0;
    CUDA_Free((unsigned char *) e3g );  e3g  = 0;
    CUDA_Free((unsigned char *) fqrg);  fqrg = 0;
    CUDA_Free((unsigned char *) fqig);  fqig = 0;

#ifdef __CUDACC__
    CUDA_Streams_Free(::cu_strms,::cu_nstrm);  ::cu_strms=0; ::cu_nstrm=0;
#endif
}

// SSF_Initialize_GPU()
//
// Initializes the CUDA thread, stream, and blocking parameters...
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Initialize_CUDA
(
   const int            nthrd ,  ///< number of cuda threads per stream block
   const int            npblk ,  ///< number of pairs per stream block
   const int            nstrm    ///< number of cuda streams
)
{
#ifdef __CUDACC__
    ::cu_props = CUDA_Device_Properties();

    if (::cu_props)
    {
        ::cu_warp = ::cu_props->warpSize;
        ::cu_tmax = ::cu_props->maxThreadsPerBlock;
    }
#endif

    ::cu_nthrd  = CONSTRAIN(nthrd,::cu_warp,::cu_tmax);  // constrain the thread count
    ::cu_npblk  = CONSTRAIN(npblk,1024,100*1024);        // constrain the segment pair count per stream block [1K:100K]

#ifdef __CUDACC__
    if (!::cu_strms || (nstrm!=::cu_nstrm) )             // initialize the active cuda streams
    {                                                    //
        CUDA_Streams_Free(::cu_strms,::cu_nstrm);        // free previous streams (if any)
                                                         //
        ::cu_nstrm = CONSTRAIN(nstrm,1,8);               // constrain active streams to [1:8]
        ::cu_strms = CUDA_Streams_Malloc(::cu_nstrm);    // allocate and initialize the streams
    }                                                    //
#endif
}

// SSF_Initialize_Iso
//
// Initialization mechanism for constants used to compute isotropic forces
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Initialize_Iso
(
   const real8     a     ,  ///< core radius (b)
   const real8     mu    ,  ///< shear modulus
   const real8     nu    ,  ///< poisson ratio
   const real8     ecrit ,  ///< critical angle for parallelism test
   const int       pbx   ,  ///< pbc active (x) (1=yes,0=no)
   const int       pby   ,  ///< pbc active (y) (1=yes,0=no)
   const int       pbz   ,  ///< pbc active (z) (1=yes,0=no)
   const real8     lx    ,  ///< simulation box size (x)
   const real8     ly    ,  ///< simulation box size (y)
   const real8     lz    ,  ///< simulation box size (z)
   const int       nthrd ,  ///< number of cuda threads per stream block
   const int       npblk ,  ///< number of pairs per stream block
   const int       nstrm    ///< number of cuda streams
)
{
    // set host-side statics...

    ::a     = a;                                       // (set simulation constants)
    ::mu    = mu;                                      //
    ::nu    = nu;                                      //
    ::ecrit = ecrit;                                   //
                                                       //
    ::pbx   = pbx;                                     // (set PBC parameters...)
    ::pby   = pby;                                     //
    ::pbz   = pbz;                                     //
                                                       //
    ::lx    = lx;                                      //
    ::ly    = ly;                                      //
    ::lz    = lz;                                      //
                                                       //
    ::sx    = ( (pbx && (lx>0.0)) ? (1.0/lx) : 0.0 );  //
    ::sy    = ( (pby && (ly>0.0)) ? (1.0/ly) : 0.0 );  //
    ::sz    = ( (pbz && (lz>0.0)) ? (1.0/lz) : 0.0 );  //

    SSF_Initialize_CUDA(nthrd,npblk,nstrm);

    ::ssf_init = 1;
}

// SSF_Initialize_Aniso()
//
// Initialization mechanism for both isotropic and anisotropic constants and tables used to compute forces.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Initialize_Aniso
(
    const  real8     a    ,   ///< core radius (b)
    const  real8     mu   ,   ///< shear modulus
    const  real8     nu   ,   ///< poisson ratio
    const  real8     ecrit,   ///< critical angle for parallelism test
    const  int       pbx  ,   ///< pbc active (x) (1=yes,0=no)
    const  int       pby  ,   ///< pbc active (y) (1=yes,0=no)
    const  int       pbz  ,   ///< pbc active (z) (1=yes,0=no)
    const  real8     lx   ,   ///< simulation box size (x)
    const  real8     ly   ,   ///< simulation box size (y)
    const  real8     lz   ,   ///< simulation box size (z)
    const  int       qmax ,   ///< spherical harmonics expansion factor
    const  int       mtype,   ///< material type index
    const  real8    *c66  ,   ///< elastic constants matrix (6x6)
    const  complex8 *c12  ,   ///< complex elements of aniso rotation matrix (1x3,complex)
    const  real8    *e3   ,   ///< real    elements of aniso rotation matrix (1x3,real   )
    const  real8    *fqr  ,   ///< Fq table (real)      size=18*(qmax+2)*(qmax+1)
    const  real8    *fqi  ,   ///< Fq table (imaginary) size=18*(qmax+2)*(qmax+1)
    const  int       nthrd,   ///< number of cuda threads per stream block
    const  int       npblk,   ///< number of pairs per stream block
    const  int       nstrm    ///< number of cuda streams
)
{
    if ( qmax>QMAX )
    {
       printf("%s::%s() ln=%d error - qmax(%d) larger than QMAX(%d) - recompile with larger QMAX!\n", __FILE__, __func__, __LINE__, qmax, QMAX );
       exit(0);
    }

    int pn = (qmax+2)*(qmax+1);   // pn = actual    size of spherical harmonics tables (qmax)
    int qn = (QMAX+2)*(QMAX+1);   // qn = allocated size of spherical harmonics tables (QMAX)
                                  // note - works for (qmax<=QMAX)

    // set host-side statics...

    ::a     = a;                                       // (set simulation constants)
    ::mu    = mu;                                      //
    ::nu    = nu;                                      //
    ::ecrit = ecrit;                                   //
                                                       //
    ::pbx   = pbx;                                     // (set PBC parameters...)
    ::pby   = pby;                                     //
    ::pbz   = pbz;                                     //
                                                       //
    ::lx    = lx;                                      //
    ::ly    = ly;                                      //
    ::lz    = lz;                                      //
                                                       //
    ::sx    = ( (pbx && (lx>0.0)) ? (1.0/lx) : 0.0 );  //
    ::sy    = ( (pby && (ly>0.0)) ? (1.0/ly) : 0.0 );  //
    ::sz    = ( (pbz && (lz>0.0)) ? (1.0/lz) : 0.0 );  //

    ::qmax  = qmax;                                    //
    ::mtype = mtype;                                   //

    // allocate and initialize the host-side aniso tables...

    if (!c66c) { c66c = new real8   [  6*6]; }  VZERO(c66c,  6*6);
    if (!c12c) { c12c = new complex8[    3]; }  VZERO(c12c,    3);
    if (!e3c ) { e3c  = new real8   [    3]; }  VZERO(e3c ,    3);
    if (!fqrc) { fqrc = new real8   [qn*18]; }  VZERO(fqrc,qn*18);
    if (!fqic) { fqic = new real8   [qn*18]; }  VZERO(fqic,qn*18);

    if (c66c && c66) { memcpy(c66c,c66,   6*6*sizeof(real8   ) ); }
    if (c12c && c12) { memcpy(c12c,c12,     3*sizeof(complex8) ); }
    if ( e3c &&  e3) { memcpy( e3c, e3,     3*sizeof(real8   ) ); }

    if (fqrc && fqr) { real8 *q=fqrc, *p=(real8*) fqr; for (int i=0; (i<18); i++) { memcpy(q,p,pn*sizeof(real8)); p+=pn; q+=qn; } }
    if (fqic && fqi) { real8 *q=fqic, *p=(real8*) fqi; for (int i=0; (i<18); i++) { memcpy(q,p,pn*sizeof(real8)); p+=pn; q+=qn; } }

    // allocate and initialize the device-side tables...

    if (!c66g ) { c66g = (real8    *) CUDA_Malloc (  6*6*sizeof(real8   ) ); }
    if (!c12g ) { c12g = (complex8 *) CUDA_Malloc (    3*sizeof(complex8) ); }
    if (!e3g  ) { e3g  = (real8    *) CUDA_Malloc (    3*sizeof(real8   ) ); }
    if (!fqrg ) { fqrg = (real8    *) CUDA_Malloc (qn*18*sizeof(real8   ) ); }
    if (!fqig ) { fqig = (real8    *) CUDA_Malloc (qn*18*sizeof(real8   ) ); }

    CUDART_CHECK( cudaMemcpy( (void *) c66g, (void *) c66c,   6*6*sizeof(real8   ), cudaMemcpyHostToDevice ));
    CUDART_CHECK( cudaMemcpy( (void *) c12g, (void *) c12c,     3*sizeof(complex8), cudaMemcpyHostToDevice ));
    CUDART_CHECK( cudaMemcpy( (void *) e3g , (void *) e3c ,     3*sizeof(real8   ), cudaMemcpyHostToDevice ));
    CUDART_CHECK( cudaMemcpy( (void *) fqrg, (void *) fqrc, qn*18*sizeof(real8   ), cudaMemcpyHostToDevice ));
    CUDART_CHECK( cudaMemcpy( (void *) fqig, (void *) fqic, qn*18*sizeof(real8   ), cudaMemcpyHostToDevice ));

    SSF_Initialize_CUDA(nthrd,npblk,nstrm);

    ::ssf_init = 1;
}

// SSF_Initialize()
//
// Will harvest the current simulation parameters from Home and initialize the SSF constants and tables.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Initialize
(
   Home_t *home  ///< points to Home data structure
)
{
   const Param_t           *param = ( home ? home->param : 0 );          // param = points to home param structure

#ifdef ANISOTROPIC
   const AnisotropicVars_t *avars = ( home ? &(home->anisoVars) : 0 );   // avars = points to the initialized aniso structure (via AnisotropicInit)
#else
   const AnisotropicVars_t *avars = 0;                                   // avars = NULL (default for isotropic forces)
#endif

   const real8     a      = ( param->rc           );                     // a     = core radius
   const real8     mu     = ( param->shearModulus );                     // mu    = shear modulus
   const real8     nu     = ( param->pois         );                     // nu    = poisson ratio
   const real8     ecrit  = ( param->ecrit        );                     // ecrit = critical angle for parallelism test
   const int       pbx    = ( param->xBoundType==Periodic ? 1 : 0 );     // pbx   = periodic boundary active (x)
   const int       pby    = ( param->yBoundType==Periodic ? 1 : 0 );     // pby   = periodic boundary active (y)
   const int       pbz    = ( param->zBoundType==Periodic ? 1 : 0 );     // pbz   = periodic boundary active (z)
   const real8     lx     = ( param->Lx );                               // lx    = simulation box size (x)
   const real8     ly     = ( param->Ly );                               // ly    = simulation box size (y)
   const real8     lz     = ( param->Lz );                               // lz    = simulation box size (z)
   const int       qmax   = ( avars ? avars->qMax : 0 );                 // qmax  = base factor for aniso spherical harmonics expansion
   const int       mtype  = ( param->materialType     );                 // mtype = material type index

         complex8  c12[3] = { Complex_t(1.0,0.0),                        // c12   = complex elements of aniso rotation matrix
                              Complex_t(0.0,1.0),                        //         (default is identity)
                              Complex_t(0.0,0.0) };                      //
         real8     e3 [3] = { 0.0, 0.0, 1.0 };                           // e3    = real elements of aniso rotation matrix
         real8    *c66    = 0;                                           // c66   = elastic constant matrix (6x6)
         real8    *fqr    = 0;                                           // fqr   = Fq table (real)
         real8    *fqi    = 0;                                           // fqi   = Fq table (imag)
         int       nthrd  = ( param->gpu_nthrd );                        // nthrd = number of active cuda threads (per kernel launch)
         int       npblk  = ( param->gpu_npblk );                        // npblk = number of segment pairs per cuda stream block
         int       nstrm  = ( param->gpu_nstrm );                        // nstrm = number of active cuda streams

   if (avars)
   {
      c66    = (real8 *) avars->elasticConstantMatrix2D;
      fqr    = (real8 *) avars->FqReal_v3;
      fqi    = (real8 *) avars->FqImag_v3;

      c12[0] = Complex_t( avars->anisoRotMatrix[0][0], avars->anisoRotMatrix[1][0] );
      c12[1] = Complex_t( avars->anisoRotMatrix[0][1], avars->anisoRotMatrix[1][1] );
      c12[2] = Complex_t( avars->anisoRotMatrix[0][2], avars->anisoRotMatrix[1][2] );

      e3 [0] = avars->anisoRotMatrix[2][0];
      e3 [1] = avars->anisoRotMatrix[2][1];
      e3 [2] = avars->anisoRotMatrix[2][2];
   }

   SSF_Initialize_Aniso(a,mu,nu,ecrit, pbx,pby,pbz, lx,ly,lz, qmax,mtype, c66,c12,e3,fqr,fqi, nthrd,npblk,nstrm);

   ::ssf_init = 1;
}

// SSF_Update_Constants()
//
// Will update the host-side simulation constants used in this module.
//
// Originally only done once, but the latest version of ParaDiS allows for
// a changing simulation box size - hence these need to be updated at each
// simulation step.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Update_Constants
(
   Home_t *home  ///< points to Home data structure
)
{
    // update the host side simulation constants...

   const Param_t           *param = ( home ? home->param : 0 );
#ifdef ANISOTROPIC
   const AnisotropicVars_t *avars = ( home ? &(home->anisoVars) : 0 );
#else
   const AnisotropicVars_t *avars = 0;
#endif

   ::a      = ( param->rc           );
   ::mu     = ( param->shearModulus );
   ::nu     = ( param->pois         );
   ::ecrit  = ( param->ecrit        );

   ::pbx    = ( param->xBoundType==Periodic ? 1 : 0 );
   ::pby    = ( param->yBoundType==Periodic ? 1 : 0 );
   ::pbz    = ( param->zBoundType==Periodic ? 1 : 0 );

   ::lx     = ( param->Lx );
   ::ly     = ( param->Ly );
   ::lz     = ( param->Lz );

   ::sx     = ( (pbx && (::lx>0.0)) ? (1.0/::lx) : 0.0 );
   ::sy     = ( (pby && (::ly>0.0)) ? (1.0/::ly) : 0.0 );
   ::sz     = ( (pbz && (::lz>0.0)) ? (1.0/::lz) : 0.0 );

   ::qmax   = ( avars ? avars->qMax : 0 );
   ::mtype  = ( param->materialType     );
}

__cuda_host__
void SSF_Update_Constants
(
   const real8   a     ,  ///< core radius (b)
   const real8   mu    ,  ///< shear modulus
   const real8   nu    ,  ///< poisson ratio
   const real8   ecrit ,  ///< critical angle for parallelism test
   const int     pbx   ,  ///< pbc active (x) (1=yes,0=no)
   const int     pby   ,  ///< pbc active (y) (1=yes,0=no)
   const int     pbz   ,  ///< pbc active (z) (1=yes,0=no)
   const real8   lx    ,  ///< simulation box size (x)
   const real8   ly    ,  ///< simulation box size (y)
   const real8   lz    ,  ///< simulation box size (z)
   const int     qmax  ,  ///< spherical harmonics expansion factor
   const int     mtype    ///< material type index
)
{
    // update the host side simulation constants...

   ::a     = a;
   ::mu    = mu;
   ::nu    = nu;
   ::ecrit = ecrit;

   ::pbx   = pbx;
   ::pby   = pby;
   ::pbz   = pbz;

   ::lx    = lx;
   ::ly    = ly;
   ::lz    = lz;

   ::sx    = ( (pbx && (lx>0.0)) ? (1.0/lx) : 0.0 );
   ::sy    = ( (pby && (ly>0.0)) ? (1.0/ly) : 0.0 );
   ::sz    = ( (pbz && (lz>0.0)) ? (1.0/lz) : 0.0 );

   ::qmax  = qmax ;
   ::mtype = mtype;
}

//------------------------------------------------------------------------------------------------------------
// In order to move the node array seamlessly between the host and gpu memory spaces, we need to be able
// to serialze the nodes into a buffer that will work on both memory spaced.
//
// The structure of the serialized buffers is basically an array of node pointers, immediately followed
// by the serialized nodes...
//
// +-----------------------------------------+
// | node pointers   (n*sizeof(SSF_Node_t *) |
// +-----------------------------------------+
// | serialized node 0                       |
// | serialized node 1                       |
// |  :      :    :                          |
// |  :      :    :                          |
// | serialized node (n-1)                   |
// +-----------------------------------------+
//
//------------------------------------------------------------------------------------------------------------

// SSF_Nodes_Serialize()
//
// Somewhat different than above in that this version manages the static serialization buffers within
// this module.  The buffers are allocated and the nodes are serialized into the host-side buffer.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Nodes_Serialize
(
   const SSF_Node_t *nodes,    ///< array of nodes
   const int         ncnt      ///< number of nodes
)
{
   size_t bytes = SSF_Nodes_Serialized_Bytes(nodes,ncnt);              // bytes = total bytes needed to serialize the nodes

   nodes_cpu = CUDA_Realloc_Host(nodes_cpu, nodes_cpu_bytes, bytes);   // allocate (host)
   nodes_gpu = CUDA_Realloc     (nodes_gpu, nodes_gpu_bytes, bytes);   // allocate (device)

   SSF_Nodes_Serialize(nodes_cpu,nodes,ncnt);  // (serialize the nodes)
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Nodes_Serialize
(
   const Node_t **narr     ,    ///< array of node pointers (may be sparse)
   const int      narr_cnt      ///< length of the node pointer array
)
{
   // determine the actual number of nodes

   int ncnt=0;
   for (int i=0; (i<narr_cnt); i++)
   { if (narr[i]) ncnt++; }

   // create and populate a temporary array of SSF_Nodes...

   SSF_Node_t *nodes = ( (ncnt>0) ? new SSF_Node_t[ncnt] : 0 );

   for (int i=0,j=0; (i<narr_cnt); i++)
   { if (narr[i]) { nodes[j]= *(narr[i]); j++; } }

   // allocate space for the serialized nodes...

   size_t bytes = ncnt*sizeof(SSF_Node_t *);     // reserve space for node pointers (based on offsets)

   for (int i=0; (i<ncnt); i++)
   {      bytes += nodes[i].Serialized_Bytes(); } // reserve space for serialized nodes

   nodes_cpu = CUDA_Realloc_Host(nodes_cpu, nodes_cpu_bytes, bytes);   // allocate (host)
   nodes_gpu = CUDA_Realloc     (nodes_gpu, nodes_gpu_bytes, bytes);   // allocate (device)

   // serialize the nodes...

   SSF_Nodes_Serialize(nodes_cpu,nodes,ncnt);

   VDELETE(nodes);  // (delete local temporary)
}

//------------------------------------------------------------------------------------------------------------
// Isotropic support...
//------------------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------------------
// SSF_Aniso_CPU() - vector API for anisotropic forces
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Aniso_CPU
(
          real8    *fv,  ///< returned forces                      (f1,f2,f3,f4,....)
          real8    *pv,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
    const int       np   ///< number of segment pairs
)
{
   if (fv && pv && (np>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<np); i++)
      { SSF_Aniso((fv+i*4*3),(pv+i*6*3), a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66c,c12c,e3c,fqrc,fqic); }
   }
}

__cuda_host__
void SSF_Aniso_CPU
(
          SSF_FV_t *fv,  ///< returned forces        (f1,f2,f3,f4,...)
    const real8    *nv,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
    const SSF_PV_t *pv,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,pv6,...)
    const int       np   ///< number of segment pairs
)
{
   if (fv && nv && pv && (np>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<np); i++)
      { SSF_Aniso((fv+i),nv,(pv+i), a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66c,c12c,e3c,fqrc,fqic); }
   }
}

//------------------------------------------------------------------------------------------------------------
// SSF_Iso_CPU() - vector API for isotropic forces
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Iso_CPU
(
          real8    *fv,  ///< returned forces                      (f1,f2,f3,f4,....)
          real8    *pv,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
    const int       np   ///< number of segment pairs
)
{
   if (fv && pv && (np>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<np); i++)
      { SSF_Iso((fv+i*4*3),(pv+i*6*3), a,mu,nu,ecrit, lx,ly,lz, sx,sy,sz); }
   }
}

__cuda_host__
void SSF_Iso_CPU
(
          SSF_FV_t *fv,  ///< returned forces        (f1,f2,f3,f4,...)
    const real8    *nv,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
    const SSF_PV_t *pv,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,pv6,...)
    const int       np   ///< number of segment pairs
)
{
   if (fv && nv && pv && (np>0))
   {
#ifdef _OPENMP
#pragma omp parallel for
#endif
      for (int i=0; (i<np); i++)
      { SSF_Iso((fv+i),nv,(pv+i), a,mu,nu,ecrit, lx,ly,lz, sx,sy,sz); }
   }
}


#ifdef __CUDACC__
//------------------------------------------------------------------------------------------------------------
// SSF_Aniso_K() - kernel routine to drive GPU-enabled aniso-forces
//------------------------------------------------------------------------------------------------------------

__cuda_global__
void SSF_Aniso_K
(
          real8    *fv   ,  ///< resulting forces                     (f1,f2,f3,f4,...)
    const real8    *nv   ,  ///< source node positions                (n1,n2,n3,n4,n5,n6,n7,...)
          real8    *pv   ,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
    const int       np   ,  ///< number of segment pairs
    const real8     a    ,  ///< core radius
    const int       qmax ,  ///< base factor for spherical harmonics expansion
    const int       mtype,  ///< material type index
    const real8     ecrit,  ///< critical angle for parallelism test
    const real8     lx   ,  ///< simulation box size (x)
    const real8     ly   ,  ///< simulation box size (y)
    const real8     lz   ,  ///< simulation box size (z)
    const real8     sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
    const real8     sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
    const real8     sz   ,  ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
    const real8    *c66  ,  ///< elastic constants matrix (6x6)
    const complex8 *c12  ,  ///< complex elements of aniso rotation matrix (1x3,complex)
    const real8    *e3   ,  ///< real    elements of aniso rotation matrix (1x3,real   )
    const real8    *fqr  ,  ///< Fq table (real)
    const real8    *fqi     ///< Fq table (imaginary)
)
{
   int bw = (blockDim .x);  // bw = cuda block  width
   int bx = (blockIdx .x);  // bx = cuda block  index
   int tx = (threadIdx.x);  // tx = cuda thread index
   int i  = bx*bw+tx;       // i  = pair index for this kernel

   if (fv && pv && (i<np) )
   {
      fv += 4*3*i;
      pv += 6*3*i;

      real8 *f1 = fv   ;
      real8 *f2 = fv+ 3;
      real8 *f3 = fv+ 6;
      real8 *f4 = fv+ 9;

      real8 *p1 = pv   ;
      real8 *p2 = pv+ 3;
      real8 *p3 = pv+ 6;
      real8 *p4 = pv+ 9;
      real8 *b1 = pv+12;
      real8 *b3 = pv+15;

      SSF_Aniso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66,c12,e3,fqr,fqi);
   }
}
#endif

#ifdef __CUDACC__
__cuda_global__
void SSF_Aniso_K
(
          SSF_FV_t *fv   ,  ///< returned forces        (f1,f2,f3,f4,...)
    const real8    *nv   ,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
    const SSF_PV_t *pv   ,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,pv6,...)
    const int       np   ,  ///< number of segment pairs
    const real8     a    ,  ///< core radius
    const int       qmax ,  ///< base factor for spherical harmonics expansion
    const int       mtype,  ///< material type index
    const real8     ecrit,  ///< critical angle for parallelism test
    const real8     lx   ,  ///< simulation box size (x)
    const real8     ly   ,  ///< simulation box size (y)
    const real8     lz   ,  ///< simulation box size (z)
    const real8     sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
    const real8     sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
    const real8     sz   ,  ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
    const real8    *c66  ,  ///< elastic constants matrix (6x6)
    const complex8 *c12  ,  ///< complex elements of aniso rotation matrix (1x3,complex)
    const real8    *e3   ,  ///< real    elements of aniso rotation matrix (1x3,real   )
    const real8    *fqr  ,  ///< Fq table (real)
    const real8    *fqi     ///< Fq table (imaginary)
)
{
   int bw = (blockDim .x);  // bw = cuda block  width
   int bx = (blockIdx .x);  // bx = cuda block  index
   int tx = (threadIdx.x);  // tx = cuda thread index
   int i  = bx*bw+tx;       // i  = pair index for this kernel

   if (fv && pv && (i<np) )
   {
            real8 *f1 =        fv[i].f1;
            real8 *f2 =        fv[i].f2;
            real8 *f3 =        fv[i].f3;
            real8 *f4 =        fv[i].f4;

      const real8 *p1 = nv + 3*pv[i].n1;
      const real8 *p2 = nv + 3*pv[i].n2;
      const real8 *p3 = nv + 3*pv[i].n3;
      const real8 *p4 = nv + 3*pv[i].n4;
      const real8 *b1 =        pv[i].b1;
      const real8 *b3 =        pv[i].b3;

      SSF_Aniso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66,c12,e3,fqr,fqi);
   }
}
#endif

#ifdef __CUDACC__
//------------------------------------------------------------------------------------------------------------
// SSF_Iso_K() - kernel routine to drive GPU-enabled aniso-forces
//------------------------------------------------------------------------------------------------------------

__cuda_global__
void SSF_Iso_K
(
          real8    *fv   ,  ///< resulting forces                     (f1,f2,f3,f4,...)
    const real8    *nv   ,  ///< source node positions                (n1,n2,n3,n4,n5,n6,n7,...)
          real8    *pv   ,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
    const int       np   ,  ///< number of segment pairs
    const real8     a    ,  ///< core radius
    const real8     mu   ,  ///< shear modulus
    const real8     nu   ,  ///< poisson ratio
    const real8     ecrit,  ///< critical angle for parallelism test
    const real8     lx   ,  ///< simulation box size (x)
    const real8     ly   ,  ///< simulation box size (y)
    const real8     lz   ,  ///< simulation box size (z)
    const real8     sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
    const real8     sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
    const real8     sz      ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
)
{
   int bw = (blockDim .x);  // bw = cuda block  width
   int bx = (blockIdx .x);  // bx = cuda block  index
   int tx = (threadIdx.x);  // tx = cuda thread index
   int i  = bx*bw+tx;       // i  = pair index for this kernel

   if (fv && pv && (i<np) )
   {
      fv += 4*3*i;
      pv += 6*3*i;

      real8 *f1 = fv   ;
      real8 *f2 = fv+ 3;
      real8 *f3 = fv+ 6;
      real8 *f4 = fv+ 9;

      real8 *p1 = pv   ;
      real8 *p2 = pv+ 3;
      real8 *p3 = pv+ 6;
      real8 *p4 = pv+ 9;
      real8 *b1 = pv+12;
      real8 *b3 = pv+15;

      SSF_Iso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit, lx,ly,lz, sx,sy,sz);
   }
}
#endif

#ifdef __CUDACC__
__cuda_global__
void SSF_Iso_K
(
          SSF_FV_t *fv   ,  ///< returned forces        (f1,f2,f3,f4,...)
    const real8    *nv   ,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
    const SSF_PV_t *pv   ,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,pv6,...)
    const int       np   ,  ///< number of segment pairs
    const real8     a    ,  ///< core radius
    const real8     mu   ,  ///< shear modulus
    const real8     nu   ,  ///< poisson ratio
    const real8     ecrit,  ///< critical angle for parallelism test
    const real8     lx   ,  ///< simulation box size (x)
    const real8     ly   ,  ///< simulation box size (y)
    const real8     lz   ,  ///< simulation box size (z)
    const real8     sx   ,  ///< reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
    const real8     sy   ,  ///< reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
    const real8     sz      ///< reciprocal of simulation box size (1.0/lz) (zero if PBC not active)
)
{
   int bw = (blockDim .x);  // bw = cuda block  width
   int bx = (blockIdx .x);  // bx = cuda block  index
   int tx = (threadIdx.x);  // tx = cuda thread index
   int i  = bx*bw+tx;       // i  = pair index for this kernel

   if (fv && pv && (i<np) )
   {
            real8 *f1 = fv[i].f1;
            real8 *f2 = fv[i].f2;
            real8 *f3 = fv[i].f3;
            real8 *f4 = fv[i].f4;

      const real8 *p1 = nv + 3*pv[i].n1;
      const real8 *p2 = nv + 3*pv[i].n2;
      const real8 *p3 = nv + 3*pv[i].n3;
      const real8 *p4 = nv + 3*pv[i].n4;
      const real8 *b1 =        pv[i].b1;
      const real8 *b3 =        pv[i].b3;

      SSF_Iso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit, lx,ly,lz, sx,sy,sz);
   }
}
#endif

#ifdef __CUDACC__
__cuda_global__
void SSF_Aniso_Reset_Forces_K
(
          SSF_Node_t *nv     ,  ///< source nodes
    const int         ncnt      ///< number of nodes
)
{
   int bw = (blockDim .x);  // bw = cuda block  width
   int bx = (blockIdx .x);  // bx = cuda block  index
   int tx = (threadIdx.x);  // tx = cuda thread index
   int i  = bx*bw+tx;       // i  = pair index for this kernel

   if (nv && (i<ncnt) )
   { nv[i].Reset_Forces(); }
}
#endif

#ifdef __CUDACC__
__cuda_global__
void SSF_Aniso_Force_Sum_K
(
          SSF_Node_t *nv   ,  ///< node structures
    const SSF_FV_t   *fv   ,  ///< source computed forces     (f1,f2,f3,f4,...)
    const SSF_PV_t   *pv   ,  ///< source pair vector structs
    const int         np      ///< number of segment pairs
)
{
   int bw = (blockDim .x);  // bw = cuda block  width
   int bx = (blockIdx .x);  // bx = cuda block  index
   int tx = (threadIdx.x);  // tx = cuda thread index
   int i  = bx*bw+tx;       // i  = pair index for this kernel

   if (nv && fv && pv && (i<np) )
   {
      const    real8 *f1 = fv[i].f1;     // f1 = force on node 1
      const    real8 *f2 = fv[i].f2;     // f2 = force on node 2
      const    real8 *f3 = fv[i].f3;     // f3 = force on node 3
      const    real8 *f4 = fv[i].f4;     // f4 = force on node 4

      unsigned short  n1 = pv[i].n1;     // n1 = node 1 index
      unsigned short  n2 = pv[i].n2;     // n2 = node 2 index
      unsigned short  n3 = pv[i].n3;     // n3 = node 3 index
      unsigned short  n4 = pv[i].n4;     // n4 = node 4 index

      unsigned long   t1 = nv[n1].tag;   // t1 = node 1 tag
      unsigned long   t2 = nv[n2].tag;   // t2 = node 2 tag
      unsigned long   t3 = nv[n3].tag;   // t3 = node 3 tag
      unsigned long   t4 = nv[n4].tag;   // t4 = node 4 tag

      nv[n1].Add_Arm_Force(t2,f2);       // add force contribution to node 1
      nv[n2].Add_Arm_Force(t1,f1);       // add force contribution to node 2
      nv[n3].Add_Arm_Force(t4,f4);       // add force contribution to node 3
      nv[n4].Add_Arm_Force(t3,f3);       // add force contribution to node 4
   }
}
#endif

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Aniso_Force_Reset_CPU
(
          SSF_Node_t *nv   ,  ///< array  of node structures
    const int         ncnt    ///< number of nodes
)
{
   if (nv && (ncnt>0))
   {
      for (int i=0; (i<ncnt); i++)
      { nv[i].Reset_Forces(); }
   }
}

__cuda_host__
void SSF_Aniso_Force_Sum_CPU
(
          SSF_Node_t *nv   ,  ///< node structures
    const SSF_FV_t   *fv   ,  ///< source computed forces     (f1,f2,f3,f4,...)
    const SSF_PV_t   *pv   ,  ///< source pair vector structs
    const int         ncnt ,  ///< number of nodes in the node vector
    const int         pcnt    ///< number of segment pairs
)
{
   SSF_Aniso_Force_Reset_CPU(nv,ncnt);

   if (nv && fv && pv && (pcnt>0) )
   {
      for (int i=0; (i<pcnt); i++)
      {
         const    real8 *f1 = fv[i].f1;     // f1 = force on node 1
         const    real8 *f2 = fv[i].f2;     // f2 = force on node 2
         const    real8 *f3 = fv[i].f3;     // f3 = force on node 3
         const    real8 *f4 = fv[i].f4;     // f4 = force on node 4

         unsigned short  n1 = pv[i].n1;     // n1 = node 1 index
         unsigned short  n2 = pv[i].n2;     // n2 = node 2 index
         unsigned short  n3 = pv[i].n3;     // n3 = node 3 index
         unsigned short  n4 = pv[i].n4;     // n4 = node 4 index

         unsigned long   t1 = nv[n1].tag;   // t1 = node 1 tag
         unsigned long   t2 = nv[n2].tag;   // t2 = node 2 tag
         unsigned long   t3 = nv[n3].tag;   // t3 = node 3 tag
         unsigned long   t4 = nv[n4].tag;   // t4 = node 4 tag

         nv[n1].Add_Arm_Force(t2,f2);       // add force contribution to node 1
         nv[n2].Add_Arm_Force(t1,f1);       // add force contribution to node 2
         nv[n3].Add_Arm_Force(t4,f4);       // add force contribution to node 3
         nv[n4].Add_Arm_Force(t3,f3);       // add force contribution to node 4
      }
   }
}

//------------------------------------------------------------------------------------------------------------
// Main driver for GPU-enabled anisotropic forces.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_GPU
(
          real8 *fv   ,  ///< resulting forces            (f1,f2,f3,f4,....)
    const real8 *nv   ,  ///< source node positions       (n1,n2,n3,n4,n5,n6,...)
    const real8 *pv   ,  ///< source pair vector structs  (p1,p2,p3,p4,b1,b3,...)
    const int    nn   ,  ///< number of source nodes
    const int    np   ,  ///< number of segment pairs
    const int    mode    ///< isotropy mode (0=isotropic,1=anisotropic)
)
{
#ifndef __CUDACC__
    printf("Error - %s::%s() ln=%d application not built with nvcc\n", __FILE__,__func__,__LINE__);
    exit(0);
#endif

#ifdef __CUDACC__
    // allocate host and device buffers...

    size_t bytes  = 4*3*np*sizeof(real8);
           bytes +=   3*nn*sizeof(real8);
           bytes += 6*3*np*sizeof(real8);

    SSF_Allocate(bytes);

    if (mbuf_cpu && mbuf_gpu)
    {
        unsigned char *pc = mbuf_cpu;
        unsigned char *pg = mbuf_gpu;

        real8 *fvc = (real8 *) pc;  pc+=4*3*np*sizeof(real8);  // fvc = resulting forces (host  )
        real8 *nvc = (real8 *) pc;  pc+=  3*nn*sizeof(real8);  // nvc = node vector      (host  )
        real8 *pvc = (real8 *) pc;  pc+=6*3*np*sizeof(real8);  // pvc = position vector  (host  )

        real8 *fvg = (real8 *) pg;  pg+=4*3*np*sizeof(real8);  // fvg = resulting forces (device)
        real8 *nvg = (real8 *) pg;  pg+=  3*nn*sizeof(real8);  // nvg = node vector      (device)
        real8 *pvg = (real8 *) pg;  pg+=6*3*np*sizeof(real8);  // pvg = position vector  (device)

        // copy the node and p-vector from heap to pinned (if they haven't already been staged)

        if ( nv && nvc && (nv!=nvc) ) { memcpy(nvc,nv,  3*nn*sizeof(real8)); }
        if ( pv && pvc && (pv!=pvc) ) { memcpy(pvc,pv,6*3*np*sizeof(real8)); }

        // send node positions and p-vectors to gpu...

        CUDART_CHECK(cudaMemcpy( (void *) nvg, (void *) nvc,   3*nn*sizeof(real8), cudaMemcpyHostToDevice ));
        CUDART_CHECK(cudaMemcpy( (void *) pvg, (void *) pvc, 6*3*np*sizeof(real8), cudaMemcpyHostToDevice ));

        // launch the compute kernel...
        {
           int  nthrd = ::cu_nthrd;
           dim3 cu_blks (((np+nthrd-1)/nthrd),1,1);
           dim3 cu_thrds(nthrd,1,1);

           if (mode==0) {  SSF_Iso_K  <<<cu_blks,cu_thrds>>>(fvg,nvg,pvg,np, a,mu,nu,     ecrit, lx,ly,lz, sx,sy,sz ); }
           else         {  SSF_Aniso_K<<<cu_blks,cu_thrds>>>(fvg,nvg,pvg,np, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66g,c12g,e3g,fqrg,fqig); }

           CUDART_CHECK(cudaGetLastError());
        }

        // get resulting forces from gpu...

        CUDART_CHECK(cudaMemcpy( (void *) fvc, (void *) fvg, 4*3*np*sizeof(real8), cudaMemcpyDeviceToHost ));

        // copy the resulting forces back to the heap (if needed)...

        if ( fv && fvc && (fv!=fvc) ) { memcpy(fv,fvc,4*3*np*sizeof(real8)); }
    }

#endif  // __CUDACC__
}

__cuda_host__
void SSF_GPU
(
          SSF_FV_t *fv  ,  ///< resulting forces            (f1,f2,f3,f4,....)
    const real8    *nv  ,  ///< source node positions       (n1,n2,n3,n4,n5,n6,...)
    const SSF_PV_t *pv  ,  ///< source pair vector structs  (p1,p2,p3,p4,p5,p6,...)
    const int       nn  ,  ///< number of source nodes
    const int       np  ,  ///< number of segment pairs
    const int       mode   ///< isotropy mode (0=isotropic,1=anisotropic)
)
{
#ifndef __CUDACC__
    printf("Error - %s::%s() ln=%d application not built with nvcc\n", __FILE__,__func__,__LINE__);
    exit(0);
#endif

#ifdef __CUDACC__
    // allocate host and device buffers...

    size_t bytes  =   np*sizeof(SSF_FV_t);
           bytes += 3*nn*sizeof(real8   );
           bytes +=   np*sizeof(SSF_PV_t);

    SSF_Allocate(bytes);

    if (mbuf_cpu && mbuf_gpu)
    {
        unsigned char *pc = mbuf_cpu;
        unsigned char *pg = mbuf_gpu;

        SSF_FV_t *fvc = (SSF_FV_t *) pc;  pc+=  np*sizeof(SSF_FV_t);  // fvc = resulting forces (host)
        real8    *nvc = (real8    *) pc;  pc+=3*nn*sizeof(real8   );  // nvc = node vector      (host)
        SSF_PV_t *pvc = (SSF_PV_t *) pc;  pc+=  np*sizeof(SSF_PV_t);  // pvc = position vector  (host)

        SSF_FV_t *fvg = (SSF_FV_t *) pg;  pg+=  np*sizeof(SSF_FV_t);  // fvg = resulting forces (device)
        real8    *nvg = (real8    *) pg;  pg+=3*nn*sizeof(real8   );  // nvg = node vector      (device)
        SSF_PV_t *pvg = (SSF_PV_t *) pg;  pg+=  np*sizeof(SSF_PV_t);  // pvg = position vector  (device)

        // copy the node and p-vector from heap to pinned (if they haven't already been staged)

        if ( nv && nvc && (nv!=nvc) ) { memcpy(nvc,nv,3*nn*sizeof(real8   )); }
        if ( pv && pvc && (pv!=pvc) ) { memcpy(pvc,pv,  np*sizeof(SSF_PV_t)); }

        // send node positions and p-vectors to gpu...

        CUDART_CHECK(cudaMemcpy( (void *) nvg, (void *) nvc, 3*nn*sizeof(real8   ), cudaMemcpyHostToDevice ));
        CUDART_CHECK(cudaMemcpy( (void *) pvg, (void *) pvc,   np*sizeof(SSF_PV_t), cudaMemcpyHostToDevice ));

        // launch the compute kernel...
        {
           int  nthrd = ::cu_nthrd;
           dim3 cu_blks (((np+nthrd-1)/nthrd),1,1);
           dim3 cu_thrds(nthrd,1,1);

           if (mode==0) { SSF_Iso_K  <<<cu_blks,cu_thrds>>>(fvg,nvg,pvg,np, a,mu,nu,     ecrit, lx,ly,lz, sx,sy,sz ); }
           else         { SSF_Aniso_K<<<cu_blks,cu_thrds>>>(fvg,nvg,pvg,np, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66g,c12g,e3g,fqrg,fqig); }

           CUDART_CHECK(cudaGetLastError());
        }

        // get resulting forces from gpu...

        CUDART_CHECK(cudaMemcpy( (void *) fvc, (void *) fvg, np*sizeof(SSF_FV_t), cudaMemcpyDeviceToHost ));

        // copy the resulting forces back to the heap (if needed)...

        if ( fv && fvc && (fv!=fvc) ) { memcpy(fv,fvc, np*sizeof(SSF_FV_t)); }
    }

#endif  // __CUDACC__
}

//------------------------------------------------------------------------------------------------------------
// Main driver for GPU-enabled anisotropic forces (streams version).
//------------------------------------------------------------------------------------------------------------

#ifdef __CUDACC__

__cuda_host__
static void SSF_NV_Send
(
          real8  *nvg,  ///< source node positions (device) (n1(xyz),n2(xyz),n3(xyz),n4(xyz),n5(xyz),n6(xyz),n7(xyz),...)
    const real8  *nvc,  ///< source node positions (host)   (n1(xyz),n2(xyz),n3(xyz),n4(xyz),n5(xyz),n6(xyz),n7(xyz),...)
    const int     nn    ///< number of nodes
)
{
    if (nvg && nvc && (nn>0))
    { CUDART_CHECK(cudaMemcpy( (void *) nvg, (void *) nvc, 3*nn*sizeof(real8), cudaMemcpyHostToDevice)); }
}

__cuda_host__
__attribute__((unused))
static void SSF_PV_Send
(
          real8        *pvg  ,  ///< source positions (device) (p1,p2,p3,p4,b1,b3,...)
          real8        *pvc  ,  ///< source positions (host)   (p1,p2,p3,p4,b1,b3,...)
    const int           n0   ,  ///< base pair index
    const int           npmax,  ///< max number of segment pairs being processed
    const int           npblk,  ///< number of segment pairs per stream block
          cudaStream_t  strm    ///< active stream
)
{
    if (pvc && pvg && (n0<npmax) )
    {
        pvc += n0*6*3;
        pvg += n0*6*3;

        int np = ( ((n0+::cu_npblk)<npmax) ? ::cu_npblk : (npmax-n0) );

        CUDART_CHECK(cudaMemcpyAsync( (void *) pvg, (void *) pvc, 6*3*np*sizeof(real8), cudaMemcpyHostToDevice, strm ));
    }
}

__cuda_host__
static void SSF_PV_Send
(
          SSF_PV_t     *pvg  ,  ///< source positions (device) (p1,p2,p3,p4,p5,p6,...)
          SSF_PV_t     *pvc  ,  ///< source positions (host)   (p1,p2,p3,p4,p5,p6,...)
    const int           n0   ,  ///< base pair index
    const int           npmax,  ///< max number of segment pairs being processed
    const int           npblk,  ///< number of segment pairs per stream block
          cudaStream_t  strm    ///< active stream
)
{
    if (pvc && pvg && (n0<npmax) )
    {
        pvc += n0;
        pvg += n0;

        int np = ( ((n0+::cu_npblk)<npmax) ? ::cu_npblk : (npmax-n0) );

        CUDART_CHECK(cudaMemcpyAsync( (void *) pvg, (void *) pvc, np*sizeof(SSF_PV_t), cudaMemcpyHostToDevice, strm ));
    }
}

__cuda_host__
static void SSF_FV_Recv
(
          SSF_FV_t     *fvc  ,  ///< resulting forces (host)   (f1,f2,f3,f4,...)
          SSF_FV_t     *fvg  ,  ///< resulting forces (device) (f1,f2,f3,f4,...)
    const int           n0   ,  ///< base pair index
    const int           npmax,  ///< max number of segment pairs being processed
    const int           npblk,  ///< number of segment pairs per stream block
          cudaStream_t  strm    ///< active stream
)
{
    if (fvc && fvg && (n0<npmax) )
    {
        fvc += n0;
        fvg += n0;

        int np = ( ((n0+npblk)<npmax) ? npblk : (npmax-n0) );

        CUDART_CHECK(cudaMemcpyAsync( (void *) fvc, (void *) fvg, np*sizeof(SSF_FV_t), cudaMemcpyDeviceToHost, strm ));
    }
}

__cuda_host__
static void SSF_Launch
(
          SSF_FV_t     *fvg  ,  ///< resulting forces        (f1,f2,f3,f4,....)
          real8        *nvg  ,  ///< source node positions   (n1,n2,n3,n4,n5,n6,...)
          SSF_PV_t     *pvg  ,  ///< source positions        (p1,p2,p3,p4,p5,p6,...)
    const int           n0   ,  ///< base pair index
    const int           npmax,  ///< number of segment pairs to compute
    const int           npblk,  ///< number of segment pairs per stream block
    const int           nthrd,  ///< number of active cuda threads per block
    const int           mode ,  ///< isotropy mode (0=isotropic,1=anisotropic)
          cudaStream_t  strm    ///< active stream
)
{
    if (fvg && pvg && (n0<npmax) )
    {
        fvg += n0;
        pvg += n0;

        int np = ( ((n0+npblk)<npmax) ? npblk : (npmax-n0) );

        dim3 cu_blks (((np+nthrd-1)/nthrd),1,1);
        dim3 cu_thrds(nthrd,1,1);

        if (mode==0) {  SSF_Iso_K  <<<cu_blks,cu_thrds,0,strm>>>(fvg,nvg,pvg,np, a,mu,nu,     ecrit, lx,ly,lz, sx,sy,sz ); }
        else         {  SSF_Aniso_K<<<cu_blks,cu_thrds,0,strm>>>(fvg,nvg,pvg,np, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66g,c12g,e3g,fqrg,fqig); }

        CUDART_CHECK(cudaGetLastError());
    }
}

#endif  //  __CUDACC__

__cuda_host__
void SSF_GPU_Strm
(
          SSF_FV_t *fv  ,   ///< resulting forces        (f1,f2,f3,f4,....)
    const real8    *nv  ,   ///< source node positions   (n1,n2,n3,n4,n5,n6,...)
    const SSF_PV_t *pv  ,   ///< source positions        (p1,p2,p3,p4,b12,b34,...)
    const int       nn  ,   ///< number of source nodes
    const int       np  ,   ///< number of segment pairs
    const int       mode    ///< isotropy mode (0=isotropic,1=anisotropic)
)
{
#ifndef __CUDACC__
    printf("Error - %s::%s() ln=%d application not built with nvcc\n", __FILE__,__func__,__LINE__);
    exit(0);
#endif

#ifdef __CUDACC__
    size_t bytes  =   np*sizeof(SSF_FV_t);
           bytes += 3*nn*sizeof(real8   );
           bytes +=   np*sizeof(SSF_PV_t);

    SSF_Allocate(bytes);

    if (mbuf_cpu && mbuf_gpu)
    {
       unsigned char *pc    = mbuf_cpu;                                       // pc    = points to current cpu buffer (host  )
       unsigned char *pg    = mbuf_gpu;                                       // pg    = points to current cpu buffer (device)

       SSF_FV_t      *fvc   = (SSF_FV_t *) pc;  pc+=  np*sizeof(SSF_FV_t);    // fvc   = resulting forces (host)
       real8         *nvc   = (real8    *) pc;  pc+=3*nn*sizeof(real8   );    // nvc   = node vector      (host)
       SSF_PV_t      *pvc   = (SSF_PV_t *) pc;  pc+=  np*sizeof(SSF_PV_t);    // pvc   = position vector  (host)

       SSF_FV_t      *fvg   = (SSF_FV_t *) pg;  pg+=  np*sizeof(SSF_FV_t);    // fvg   = resulting forces (device)
       real8         *nvg   = (real8    *) pg;  pg+=3*nn*sizeof(real8   );    // nvg   = node vector      (device)
       SSF_PV_t      *pvg   = (SSF_PV_t *) pg;  pg+=  np*sizeof(SSF_PV_t);    // pvg   = position vector  (device)

       int            n0    = 0;                                              // n0    = send base pair index
       int            k0    = 0;                                              // k0    = recv/kernel base pair index
       int            sndx  = 0;                                              // sndx  = send stream index
       int            kndx  = 0;                                              // kndx  = recv stream index

       int            nthrd = ::cu_nthrd;                                     // nthrd = number of cuda threads per kernel launch
       int            npblk = ::cu_npblk;                                     // npblk = number of segment pairs per stream block
       int            nstrm = ::cu_nstrm;                                     // nstrm = number of active cuda streams

       // copy the node and p-vector from heap to pinned (if they haven't already been staged)

       if ( nv && nvc && (nv!=nvc) ) { memcpy(nvc,nv,3*nn*sizeof(real8   )); }
       if ( pv && pvc && (pv!=pvc) ) { memcpy(pvc,pv,  np*sizeof(SSF_PV_t)); }

       // send the node vector (note - synchronous, used by ALL streams)...

       SSF_NV_Send(nvg,nvc,nn);

       // send the initial P-vectors to each stream...

       for (int i=0; (i<nstrm); i++)
       {   SSF_PV_Send (pvg,pvc,    n0,np,npblk,            ::cu_strms[sndx%nstrm]);  n0+=npblk; sndx++; }

       while (k0<np)
       {
           SSF_Launch  (fvg,nvg,pvg,k0,np,npblk,nthrd,mode, ::cu_strms[kndx%nstrm]);
           SSF_FV_Recv (fvc,fvg,    k0,np,npblk,            ::cu_strms[kndx%nstrm]);  k0+=npblk; kndx++;
           SSF_PV_Send (pvg,pvc,    n0,np,npblk,            ::cu_strms[sndx%nstrm]);  n0+=npblk; sndx++;
       }

       // block until all the work on the GPU completes and the results have been received...

       for (int i=0; (i<nstrm); i++)
       { CUDART_CHECK(cudaStreamSynchronize(::cu_strms[i])); }

       // copy the resulting forces back to the heap (if needed)

       if ( fv && fvc && (fv!=fvc) ) { memcpy(fv,fvc,np*sizeof(SSF_FV_t)); }
    }

#endif  // __CUDACC__
}

//------------------------------------------------------------------------------------------------------------
// SSF_Index_Nodes()
//
// Given a sparse array of node pointers, will traverse the array and set the indices of each
// node to a monotonically increasing index.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Index_Nodes
(
         Node_t **nodes ,  ///< sparse array of node pointers (can include null pointers)
   const int      ncnt     ///< number of node pointers
)
{
   if (nodes && (ncnt>0))
   {
      for (int i=0,j=0; (i<ncnt); i++)
         if (nodes[i]) { nodes[i]->myIndex = j++;  }
   }
}

//------------------------------------------------------------------------------------------------------------
// SSF_Gather_Nodes
//
// Will construct an N-vector of ALL the nodes used by this process (native+ghost).
// After constructing the list - indexes the nodes from [0:(n-1)].  The resulting
// list is used to copy the node positions to the GPU.
//
// Note that the source native and ghost pointer arrays can contain null pointers as a
// result of topology operations on the nodes. The returned array and node count will
// is a dense array containing only the non-null pointers.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
Node_t **SSF_Gather_Nodes
(
         int    & nn     ,  ///< returns the actual node count (native+ghost)
         Node_t **native ,  ///< sparse array of native node pointers (can include null pointers)
   const int      ncnt   ,  ///< number of native node pointers
         Node_t **ghosts ,  ///< sparse array of ghost  node pointers (can include null pointers)
   const int      gcnt      ///< number of ghost node pointers
)
{
   nn=0;

   // allocate space to hold the pointers...

   Node_t **nodes = ( ((ncnt+gcnt)>0) ? new Node_t *[ncnt+gcnt] : 0 );

   // if allocated - build a contiguous, dense array of node pointers...

   if (nodes)
   {
      memset(nodes,0,(ncnt+gcnt)*sizeof(Node_t *));

      int n=0;
      if (native && (ncnt>0)) { for (int i=0; (i<ncnt); i++) { if (native[i]) { nodes[n++] = native[i]; } } }
      if (ghosts && (gcnt>0)) { for (int i=0; (i<gcnt); i++) { if (ghosts[i]) { nodes[n++] = ghosts[i]; } } }

      // ...index ALL the nodes from [0:(n-1)]...

      for (int i=0; (i<n); i++) { nodes[i]->myIndex=i; }

      nn=n;  // (return actual node count)
   }

   return(nodes);
}

__cuda_host__
void SSF_Gather_Nodes
(
         real8   *nv   ,  ///< source node positions   (n1,n2,n3,n4,n5,n6,...)
         Node_t **nodes,  ///< dense array of node pointers (assembled above)
   const int      nn      ///< number of nodes
)
{
   if (nv && nodes && (nn>0))
   {
      for (int i=0; (i<nn); i++)
      {
         Node_t *node = (Node_t *) nodes[i];
         real8  *q    = nv + (3*node->myIndex);

         q[0] = node->x;
         q[1] = node->y;
         q[2] = node->z;
      }
   }
}

// SSF_Send_Nodes()
//
// Will allocate, serialize, and send the node vector to the device.
// Note - this routine will leave the serialized node vector and nodes remapped for device memory.
// You will need to remap them back to the host to access them via the host.
//
// This routine is expected to be used in conjuction with the SSF_Recv_Nodes() routine
// which will remap the serialized nodes back to their host locations.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Send_Nodes
(
         Node_t **narr   ,  ///< array of node pointers
   const int      ncnt      ///< node count
)
{
   SSF_Allocate_Nodes(narr,ncnt);

   SSF_Node_t *nodes = ( (ncnt>0) ? new SSF_Node_t[ncnt] : 0 );

   if (nodes)
   {
      for (int i=0; (i<ncnt); i++)
      { if (narr[i]) { nodes[i] = *narr[i]; } }
   }

   SSF_Nodes_Serialize (nodes_cpu,nodes,ncnt);
   SSF_Nodes_Send      (nodes_cpu,nodes_gpu,nodes_cpu_bytes,ncnt);

   VDELETE(nodes)
}

// SSF_Recv_Nodes()
//
// Will copy the serialized nodes from the device back to the host.
// Note - this routine will remap the serialized node vector and nodes back to host memory.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Recv_Nodes
(
         Node_t **narr   ,  ///< array of node pointers
   const int      ncnt      ///< node count
)
{
   // retrieve the serialized nodes from the gpu...

   if (nodes_cpu && nodes_gpu && (nodes_cpu_bytes>0))
   { CUDART_CHECK( cudaMemcpy( (void *) nodes_cpu, (void *) nodes_gpu, nodes_cpu_bytes, cudaMemcpyDeviceToHost )); }

   SSF_Nodes_Remap(nodes_cpu,ncnt);   // remap them for host memory
}

// Arm_Index()
//
// Returns the arm index of node 2 within node 1 (called by Add_Forces() below).
//------------------------------------------------------------------------------------------------------------

__cuda_host__
static int Arm_Index
(
   Node_t *n1,  ///< points to node 1
   Node_t *n2   ///< points to node 2
)
{
   if (n1 && n2 && n1->nbrTag)
   {
      int    narms = n1->numNbrs;      // narms = number of arms on node 1
      Tag_t *tags  = n1->nbrTag ;      // tags  = array of arm tags on node 1
      Tag_t  t2    = n2->myTag  ;      // t2    = tag of node 2

      for(int i=0; (i<narms); i++)     // loop through arm tags...
      { if (tags[i]==t2) return(i); }  // (found)
   }

   return(-1);  // (not found)
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Gather_Pairs
(
         real8         *pv   ,  ///< source positions        (p1,p2,p3,p4,b12,b34,...)
         SegmentPair_t *sp   ,  ///< array  of source segment pairs
   const int            np      ///< number of source segment pairs
)
{
   if (pv && sp && (np>0))
   {
      for (int i=0,k=0; (i<np); i++)
      {
         Segment_t *s1 = sp[i].seg1;                            // s1 = points to segment 1
         Segment_t *s2 = sp[i].seg2;                            // s2 = points to segment 2

         Node_t    *n1 = s1->node1;                             // n1 = points to node 1
         Node_t    *n2 = s1->node2;                             // n2 = points to node 2
         Node_t    *n3 = s2->node1;                             // n3 = points to node 3
         Node_t    *n4 = s2->node2;                             // n4 = points to node 4

         int j1 = Arm_Index(n1,n2);                             // j1 = arm index of node 1 to node 2
         int j3 = Arm_Index(n3,n4);                             // j3 = arm index of node 3 to node 4

         pv[k]=n1->x;  pv[k+1]=n1->y;  pv[k+2]=n1->z;   k+=3;   // set node 1 position (xyz)
         pv[k]=n2->x;  pv[k+1]=n2->y;  pv[k+2]=n2->z;   k+=3;   // set node 2 position (xyz)
         pv[k]=n3->x;  pv[k+1]=n3->y;  pv[k+2]=n3->z;   k+=3;   // set node 3 position (xyz)
         pv[k]=n4->x;  pv[k+1]=n4->y;  pv[k+2]=n4->z;   k+=3;   // set node 4 position (xyz)

         pv[k+0]=n1->burgX[j1];                                 // set burgers vector (p1 to p2) (xyz)
         pv[k+1]=n1->burgY[j1];                                 //
         pv[k+2]=n1->burgZ[j1];  k+=3;                          //

         pv[k+0]=n3->burgX[j3];                                 // set burgers vector (p3 to p4) (xyz)
         pv[k+1]=n3->burgY[j3];                                 //
         pv[k+2]=n3->burgZ[j3];  k+=3;                          //
      }
   }
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Gather_Pairs
(
         SSF_PV_t      *pv   ,  ///< source positions        (p1,p2,p3,p4,b12,b34,...)
         SegmentPair_t *sp   ,  ///< array  of source segment pairs
   const int            np      ///< number of source segment pairs
)
{
   if (pv && sp && (np>0))
   {
      for (int i=0; (i<np); i++)
      {
         Segment_t *s1 = sp[i].seg1;  // s1 = points to segment 1
         Segment_t *s2 = sp[i].seg2;  // s2 = points to segment 2

         Node_t    *n1 = s1->node1;   // n1 = points to node 1
         Node_t    *n2 = s1->node2;   // n2 = points to node 2
         Node_t    *n3 = s2->node1;   // n3 = points to node 3
         Node_t    *n4 = s2->node2;   // n4 = points to node 4

         pv[i].n1 = n1->myIndex;      // set node indices
         pv[i].n2 = n2->myIndex;      //
         pv[i].n3 = n3->myIndex;      //
         pv[i].n4 = n4->myIndex;      //

         int j1 = Arm_Index(n1,n2);   // j1 = arm index of node 1 to node 2
         int j3 = Arm_Index(n3,n4);   // j3 = arm index of node 3 to node 4

         pv[i].b1[0]=n1->burgX[j1];   // set burgers vector (p1 to p2) (xyz)
         pv[i].b1[1]=n1->burgY[j1];   //
         pv[i].b1[2]=n1->burgZ[j1];   //

         pv[i].b3[0]=n3->burgX[j3];   // set burgers vector (p3 to p4) (xyz)
         pv[i].b3[1]=n3->burgY[j3];   //
         pv[i].b3[2]=n3->burgZ[j3];   //
      }
   }
}

// Add_Forces
//
// Given a pointer to a segment and two forces, will add the forces to their respective arms on the
// nodes identified by the segment.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
static void Add_Forces
(
         Segment_t *seg,  ///< points to segment
   const real8     *f1 ,  ///< points to force contribution of node 2 on node 1 <xyz>
   const real8     *f2    ///< points to force contribution of node 1 on node 2 <xyz>
)
{
   if (seg)
   {
      seg->forcesSet = 1;           // indicate that we've set forces on this segment

      Node_t *n1 = seg->node1;      // n1 = points to node 1
      Node_t *n2 = seg->node2;      // n2 = points to node 2

      int j12 = Arm_Index(n1,n2);   // j12 = index of node 2 arm in node 1
      int j21 = Arm_Index(n2,n1);   // j21 = index of node 1 arm in node 2

      n1->armfx[j12] += f1[0];      // add f1 to node arm
      n1->armfy[j12] += f1[1];      //
      n1->armfz[j12] += f1[2];      //

      n2->armfx[j21] += f2[0];      // add f2 to node arm
      n2->armfy[j21] += f2[1];      //
      n2->armfz[j21] += f2[2];      //

      seg->f1[0]     += f1[0];      // add f1 to segment force 1
      seg->f1[1]     += f1[1];      //
      seg->f1[2]     += f1[2];      //

      seg->f2[0]     += f2[0];      // add f2 to segment force 2
      seg->f2[1]     += f2[1];      //
      seg->f2[2]     += f2[2];      //
   }
}

// SSF_Scatter_Forces()
//
// This routine will unpack the computed segment forces and re-insert them back into the segment and
// node data structures. Note that this routine assumes the segment pairs were packed and computed
// using the gather routines within this module.
//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Scatter_Forces
(
         SegmentPair_t *sp  , ///< array  of source segment pairs
         SSF_FV_t      *fv  , ///< computed segment pair forces
   const int            np    ///< number of source segment pairs
)
{
   if (sp && fv && (np>0))
   {
      for (int i=0; (i<np); i++)
      {
         Segment_t *s1 = sp[i].seg1;  // s1 = points to segment 1
         Segment_t *s2 = sp[i].seg2;  // s2 = points to segment 2

         real8     *f1 = fv[i].f1;    // f1 = points to force on node 1 <xyz>
         real8     *f2 = fv[i].f2;    // f2 = points to force on node 2 <xyz>
         real8     *f3 = fv[i].f3;    // f3 = points to force on node 3 <xyz>
         real8     *f4 = fv[i].f4;    // f4 = points to force on node 4 <xyz>

         if (sp[i].setSeg1Forces) { Add_Forces(s1,f1,f2); }
         if (sp[i].setSeg2Forces) { Add_Forces(s2,f3,f4); }
      }
   }
}

//------------------------------------------------------------------------------------------------------------

__cuda_host__
void SSF_Compute_Forces
(
         Home_t        *home,  ///< points to Home data structure
         SegmentPair_t *sp  ,  ///< array  of source segment pairs
   const int            np  ,  ///< number of source segment pairs
   const int            mode   ///< isotropy mode (0=isotropic,1=anisotropic)
)
{
   // Initialize the aniso tables and constants (host+device)...

   if (!::ssf_init)
   { SSF_Initialize(home); ::ssf_init=1; }

   // Update any constants that may have changed from the previous cycle...

   SSF_Update_Constants(home);

   // Gather a dense array of node pointers for all the nodes (native and ghost) in this process

   Node_t **native = home->nodeKeys      ;   // native = sparse array of native node pointers
   int      ncnt   = home->newNodeKeyPtr ;   // ncnt   = number of native node pointers
   Node_t **ghosts = home->ghostNodeList ;   // ghosts = sparse array of ghost node pointers
   int      gcnt   = home->ghostNodeCount;   // gcnt   = number of ghost node pointers

   int      nn     = 0;
   Node_t **nodes  = SSF_Gather_Nodes(nn,native,ncnt,ghosts,gcnt);

   int gpu_enabled = ( ( home && home->param ) ? home->param->gpu_enabled : 0 );

   SSF_Allocate(nn,np);

   if (mbuf_cpu)
   {
      unsigned char *pc = mbuf_cpu;                                   // pc    = points to current cpu buffer (host)

      SSF_FV_t *fv = (SSF_FV_t *) pc;  pc+=  np*sizeof(SSF_FV_t);     // fvc   = resulting forces      (host)
      real8    *nv = (real8    *) pc;  pc+=3*nn*sizeof(real8   );     // nvc   = node vector           (host)
      SSF_PV_t *pv = (SSF_PV_t *) pc;  pc+=  np*sizeof(SSF_PV_t);     // pvc   = pv structure  vector  (host)

      SSF_Gather_Nodes(nv,nodes,nn);                                  // serialize the node positions
      SSF_Gather_Pairs(pv,sp,np);                                     // serialize the segment pairs

      if (gpu_enabled) { SSF_GPU_Strm (fv,nv,pv,nn,np,mode); }        // compute forces using gpu
      else
      {
         if (mode==0)  { SSF_Iso_CPU  (fv,nv,pv,np); }                // compute iso   forces using cpu
         else          { SSF_Aniso_CPU(fv,nv,pv,np); }                // compute aniso forces using cpu
      }

      SSF_Scatter_Forces(sp,fv,np);                                   // distribute the resulting forces back to nodes
   }

   VDELETE(nodes);
}

