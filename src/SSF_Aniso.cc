#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Complex_t.h"
#include "PBC.h"
#include "V3.h"
#include "Anisotropic.h"
#include "AnisotropicVars_t.h"
#include "SSF_Aniso.h"

// SSF_Seg_Zero()
//
// There are some exceptional cases that can result in a zero-length segment.
// That can be problematic when computing forces. This routine will return 1
// if a segment is too short to compute forces. 
//
// Note - this routine avoids the sqrt() call and compares on the squared distance.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
inline int SSF_Seg_Zero
(
   const real8 *p1   ,         ///< position of node 1 <xyz>
   const real8 *p2             ///< position of node 2 <xyz>
)
{
   const real8 dx  = (p2[0]-p1[0]);
   const real8 dy  = (p2[1]-p1[1]);
   const real8 dz  = (p2[2]-p1[2]);
   const real8 dot = (dx*dx)+(dy*dy)+(dz*dz);
   const real8 eps = 1.0e-20;

   return( (dot<eps) ? 1 : 0 );
}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso
(
          real8    *f1   ,  ///< resulting force on node 1 (xyz)
          real8    *f2   ,  ///< resulting force on node 2 (xyz)
          real8    *f3   ,  ///< resulting force on node 3 (xyz)
          real8    *f4   ,  ///< resulting force on node 4 (xyz)
    const real8    *p1   ,  ///< node 1 position (xyz)
    const real8    *p2   ,  ///< node 2 position (xyz)
    const real8    *p3   ,  ///< node 3 position (xyz)
    const real8    *p4   ,  ///< node 4 position (xyz)
    const real8    *b1   ,  ///< burgers vector (p1->p2) (xyz)
    const real8    *b3   ,  ///< burgers vector (p3->p4) (xyz)
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
   real8   x1[3] = { p1[0], p1[1], p1[2] };                 // x1 = (local copy of p1)
   real8   x2[3] = { p2[0], p2[1], p2[2] };                 // x2 = (local copy of p2)
   real8   x3[3] = { p3[0], p3[1], p3[2] };                 // x3 = (local copy of p3)
   real8   x4[3] = { p4[0], p4[1], p4[2] };                 // x4 = (local copy of p4)

   PBC_Position(x1,x2,x3,x4, lx,ly,lz, sx,sy,sz);           // (correct for periodic boundaries)

   real8   t12[3];  V3_SUB(t12,x2,x1);                      // t12 = line direction    (p1->p2)
   real8   t34[3];  V3_SUB(t34,x4,x3);                      // t34 = line direction    (p3->p4)

   real8   s12    = sqrt( V3_DOT(t12,t12) );                // s12 = length of segment (p1->p2)
   real8   s34    = sqrt( V3_DOT(t34,t34) );                // s34 = length of segment (p3->p4)

           s12    = ( (s12>0.0) ? (1.0/s12) : 0.0 );        // s12 = reciprocol(s12) (or zero)
           s34    = ( (s34>0.0) ? (1.0/s34) : 0.0 );        // s34 = reciprocol(s34) (or zero)

   V3_SCALE(t12,t12,s12);                                   // t12 = normalized line direction (p1->p2) (or zero)
   V3_SCALE(t34,t34,s34);                                   // t34 = normalized line direction (p3->p4) (or zero)

   real8 c        = V3_DOT(t12,t34);
   int   parallel = ( (fabs(1.0-(c*c)) < ecrit) ? 1 : 0 );

   // initialize resulting forces to zero...

   V3_ZERO(f1);  // f1 = [0,0,0]
   V3_ZERO(f2);  // f2 = [0,0,0]
   V3_ZERO(f3);  // f3 = [0,0,0]
   V3_ZERO(f4);  // f4 = [0,0,0]

   // compute forces for non-zero-length pairs...

   if ( (s12>0.0) && (s34>0.0) )
   {
      if (parallel) { SSF_Aniso_Parallel    (f1,f2,f3,f4,  x1,x2,x3,x4, b1,b3, a,qmax,mtype,ecrit, c66,c12,e3,fqr,fqi); }
      else          { SSF_Aniso_NonParallel (f1,f2,f3,f4,  x1,x2,x3,x4, b1,b3, a,qmax,mtype,       c66,c12,e3,fqr,fqi); }
   }
}

__cuda_hdev__
void SSF_Aniso
(
          real8    *fv   ,  ///< resulting forces                     (f1,f2,f3,f4,...)
    const real8    *pv   ,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
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
         real8  *f1 = fv   ;
         real8  *f2 = fv+ 3;
         real8  *f3 = fv+ 6;
         real8  *f4 = fv+ 9;

   const real8  *p1 = pv   ;
   const real8  *p2 = pv+ 3;
   const real8  *p3 = pv+ 6;
   const real8  *p4 = pv+ 9;
   const real8  *b1 = pv+12;
   const real8  *b3 = pv+15;

   SSF_Aniso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66,c12,e3,fqr,fqi);
}

__cuda_hdev__
void SSF_Aniso
(
          SSF_FV_t *fv   ,  ///< returned forces        (f1,f2,f3,f4,...)
    const real8    *nv   ,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
    const SSF_PV_t *pv   ,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,...)
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
         real8 *f1 = fv->f1;
         real8 *f2 = fv->f2;
         real8 *f3 = fv->f3;
         real8 *f4 = fv->f4;

   const real8 *p1 = nv + 3*pv->n1;
   const real8 *p2 = nv + 3*pv->n2;
   const real8 *p3 = nv + 3*pv->n3;
   const real8 *p4 = nv + 3*pv->n4;
   const real8 *b1 =        pv->b1;
   const real8 *b3 =        pv->b3;

   SSF_Aniso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66,c12,e3,fqr,fqi);
}

// SSF_Aniso()
//
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Aniso
(
         real8  *f1        ,    ///< resulting force on node 1 <xyz>
         real8  *f2        ,    ///< resulting force on node 2 <xyz>
         real8  *f3        ,    ///< resulting force on node 3 <xyz>
         real8  *f4        ,    ///< resulting force on node 4 <xyz>
   const real8  *p1        ,    ///< position of node 1 <xyz>
   const real8  *p2        ,    ///< position of node 2 <xyz>
   const real8  *p3        ,    ///< position of node 3 <xyz>
   const real8  *p4        ,    ///< position of node 4 <xyz>
   const real8  *b1        ,    ///< burgers vector for segment p1->p2
   const real8  *b3        ,    ///< burgers vector for segment p3->p4
         Home_t *home           ///< points to home data structure
)
{
   // retrieve the aniso constants and tables from home...

   const Param_t           *param = ( home ? home->param : 0 );                         // param = points to home param structure

#ifdef ANISOTROPIC
   const AnisotropicVars_t *avars = ( home ? &(home->anisoVars) : 0 );                  // avars = points to the initialized aniso structure (via AnisotropicInit)
#else
   const AnisotropicVars_t *avars = 0;                                                  // avars = NULL (default for isotropic forces)
#endif

   const real8     a      = ( param->rc           );                                    // a     = core radius
   const int       qmax   = ( avars ? avars->qMax : 0 );                                // qmax  = base factor for aniso spherical harmonics expansion
   const int       mtype  = ( param->materialType     );                                // mtype = material type index
   const real8     ecrit  = ( param->ecrit        );                                    // ecrit = critical angle for parallelism test

   const int       pbx    = ( param->xBoundType==Periodic ? 1 : 0 );                    // pbx   = periodic boundary active (x)
   const int       pby    = ( param->yBoundType==Periodic ? 1 : 0 );                    // pby   = periodic boundary active (y)
   const int       pbz    = ( param->zBoundType==Periodic ? 1 : 0 );                    // pbz   = periodic boundary active (z)

   const real8     lx     = ( param->Lx );                                              // lx    = simulation box size (x)
   const real8     ly     = ( param->Ly );                                              // ly    = simulation box size (y)
   const real8     lz     = ( param->Lz );                                              // lz    = simulation box size (z)

   const real8     sx     = ( (pbx && (lx>0.0)) ? (1.0/lx) : 0.0 );                     // sx    = reciprocal of simulation box size (1.0/lx) (zero if PBC not active)
   const real8     sy     = ( (pby && (ly>0.0)) ? (1.0/ly) : 0.0 );                     // sy    = reciprocal of simulation box size (1.0/ly) (zero if PBC not active)
   const real8     sz     = ( (pbz && (lz>0.0)) ? (1.0/lz) : 0.0 );                     // sz    = reciprocal of simulation box size (1.0/lz) (zero if PBC not active)

         complex8  c12[3] = { Complex_t(1.0,0.0),                                       // c12   = complex elements of aniso rotation matrix
                              Complex_t(0.0,1.0),                                       //         (default is identity)
                              Complex_t(0.0,0.0) };                                     //
         real8     e3 [3] = { 0.0, 0.0, 1.0 };                                          // e3    = real elements of aniso rotation matrix

   const real8    *c66    = (real8 *) ( avars ? avars->elasticConstantMatrix2D : 0 );   // c66   = elastic constant matrix (6x6)
   const real8    *fqr    = (real8 *) ( avars ? avars->FqReal_v3               : 0 );   // fqr   = Fq table (real)
   const real8    *fqi    = (real8 *) ( avars ? avars->FqImag_v3               : 0 );   // fqi   = Fq table (imag)

   if (avars)
   {
      c12[0] = Complex_t( avars->anisoRotMatrix[0][0], avars->anisoRotMatrix[1][0] );
      c12[1] = Complex_t( avars->anisoRotMatrix[0][1], avars->anisoRotMatrix[1][1] );
      c12[2] = Complex_t( avars->anisoRotMatrix[0][2], avars->anisoRotMatrix[1][2] );

      e3 [0] = avars->anisoRotMatrix[2][0];
      e3 [1] = avars->anisoRotMatrix[2][1];
      e3 [2] = avars->anisoRotMatrix[2][2];
   }

   // set default result (zero)...

   if (f1) { V3_ZERO(f1); }
   if (f2) { V3_ZERO(f2); }
   if (f3) { V3_ZERO(f3); }
   if (f4) { V3_ZERO(f4); }

   // cope with zero-length segments...

   if ( SSF_Seg_Zero(p1,p2) || SSF_Seg_Zero(p3,p4) ) { return; }

   SSF_Aniso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,qmax,mtype,ecrit, lx,ly,lz, sx,sy,sz, c66,c12,e3,fqr,fqi);
}

