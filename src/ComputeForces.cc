//------------------------------------------------------------------------------------------------------------

#include "mpi_portability.h"

#include "Home.h"
#include "CPU_Timer.h"
#include "PBC.h"
#include "V3.h"
#include "ComputeForces.h"

//------------------------------------------------------------------------------------------------------------
// ComputeForces
//
// Obtains nodal information (i.e coordinates, burgers vectors) of the endpoints defining
// the segments, adjusts the coordinates (if necessary) for boundary conditions, invokes
// the function to calculate the interaction between the two segments, and returns the
// resulting fore at all four endpoints.
//------------------------------------------------------------------------------------------------------------

void ComputeForces
(
   Home_t *home,   ///< points to the home data structure (used to extract a,MU,NU parameters).
   Node_t *n1  ,   ///< points to node 1 (segment 1, node 1)
   Node_t *n2  ,   ///< points to node 2 (segment 1, node 2)
   Node_t *n3  ,   ///< points to node 3 (segment 2, node 1)
   Node_t *n4  ,   ///< points to node 4 (segment 2, node 2)
   real8  *f1  ,   ///< resulting force on node 1 <xyz> (returned)
   real8  *f2  ,   ///< resulting force on node 2 <xyz> (returned)
   real8  *f3  ,   ///< resulting force on node 3 <xyz> (returned)
   real8  *f4      ///< resulting force on node 4 <xyz> (returned)
)
{
   const Param_t *param  = home->param;      // param = points to home param structure

   const real8    a  = param->rc;            // a  = core radius (b)
   const real8    mu = param->shearModulus;  // mu = shear modulus
   const real8    nu = param->pois;          // nu = poisson ratio
   const real8    lx = param->Lx;            // lx = simulation box size (x)
   const real8    ly = param->Ly;            // ly = simulation box size (y)
   const real8    lz = param->Lz;            // lz = simulation box size (z)

   const real8    sx = ( (param->xBoundType == Periodic) ? 1.0/lx : 0.0 );  // sx = ( PBC active ? 1.0/lx : 0 );
   const real8    sy = ( (param->yBoundType == Periodic) ? 1.0/ly : 0.0 );  // sy = ( PBC active ? 1.0/ly : 0 );
   const real8    sz = ( (param->zBoundType == Periodic) ? 1.0/lz : 0.0 );  // sz = ( PBC active ? 1.0/lz : 0 );

   ComputeForces(home, n1,n2,n3,n4, f1,f2,f3,f4, a,mu,nu, lx,ly,lz, sx,sy,sz);
}

//------------------------------------------------------------------------------------------------------------

void ComputeForces
(
         Home_t *home,   ///< points to the home data structure (used to extract a,MU,NU parameters).
         Node_t *n1  ,   ///< points to node 1 (segment 1, node 1)
         Node_t *n2  ,   ///< points to node 2 (segment 1, node 2)
         Node_t *n3  ,   ///< points to node 3 (segment 2, node 1)
         Node_t *n4  ,   ///< points to node 4 (segment 2, node 2)
         real8  *f1  ,   ///< resulting force on node 1 <xyz> (returned)
         real8  *f2  ,   ///< resulting force on node 2 <xyz> (returned)
         real8  *f3  ,   ///< resulting force on node 3 <xyz> (returned)
         real8  *f4  ,   ///< resulting force on node 4 <xyz> (returned)
   const real8   a   ,   ///< core radius (b)
   const real8   mu  ,   ///< shear modulus
   const real8   nu  ,   ///< poisson ratio
   const real8   lx  ,   ///< simulation box size (x)
   const real8   ly  ,   ///< simulation box size (y)
   const real8   lz  ,   ///< simulation box size (z)
   const real8   sx  ,   ///< reciprocal of simulation box size (1.0/lx) (or zero if PBC not active)
   const real8   sy  ,   ///< reciprocal of simulation box size (1.0/ly) (or zero if PBC not active)
   const real8   sz      ///< reciprocal of simulation box size (1.0/lz) (or zero if PBC not active)
)
{
   // Initialize the returned forces to zero...

   if (f1) { V3_ZERO(f1); }
   if (f2) { V3_ZERO(f2); }
   if (f3) { V3_ZERO(f3); }
   if (f4) { V3_ZERO(f4); }

   // This function is only used during the full force calculations
   // where forces for any given segment pair are calculated by only
   // one domain in the problem.  This means that even if one of the
   // segments is non-local, we still need its forces... so in this
   // routine, we treat all segments as local so SegSegForce() sets
   // forces for all nodes.

   real8    p1[3] = { n1->x, n1->y, n1->z };
   real8    p2[3] = { n2->x, n2->y, n2->z };
   real8    p3[3] = { n3->x, n3->y, n3->z };
   real8    p4[3] = { n4->x, n4->y, n4->z };

   PBC_Position(p1,p2, lx,ly,lz, sx,sy,sz);    // (shifts p2 to nearest PBC image of p1)
   PBC_Position(p3,p4, lx,ly,lz, sx,sy,sz);    // (shifts p4 to nearest PBC image of p3)

   // Collision handling can cause one of the segments to be zero-length.  If either segment
   // is too small, there will be no resulting seg/seg forces. They can also create a numerical
   // instability problem in the parallelism test. If either segment is too small - we just
   // return (we've already set the forces to zero above).

   if ( V3_DOT(p1,p2) < 1.0e-20 ) { return; }
   if ( V3_DOT(p3,p4) < 1.0e-20 ) { return; }

   // We now need to shift the position of the second segment (nodes 3 and 4) to
   // the positions of the nearest image of segment 2 with respect to segment 1.

   real8 pa[3] = { (p1[0]+p2[0])/2.0,            // pa = midpoint of segment 1
                   (p1[1]+p2[1])/2.0,            // pa = (p1+p2)/2
                   (p1[2]+p2[2])/2.0 };          //

   real8 pb[3] = { (p3[0]+p4[0])/2.0,            // pb = midpoint of segment 2
                   (p3[1]+p4[1])/2.0,            // pb = (p3+p4)/2
                   (p3[2]+p4[2])/2.0 };          //

   real8 pc[3] = { -rint((pb[0]-pa[0])*sx)*lx,   // pc = pbc adjustment for segment 2
                   -rint((pb[1]-pa[1])*sy)*ly,   //
                   -rint((pb[2]-pa[2])*sz)*lz }; //

   V3_ADD(p3,p3,pc);                             // p3 = (p3+pc) (adjust p3 to nearest image)
   V3_ADD(p4,p4,pc);                             // p4 = (p4+pc) (adjust p4 to nearest image)

   // get the burgers vectors...

   int    j1 = GetArmID(n1,n2);                  // j1 = index of node 2 arm wrt node 1
   int    j3 = GetArmID(n3,n4);                  // j3 = index of node 4 arm wrt node 3

   real8  b1[3]   = { n1->burgX[j1],             // b1 = burgers vector (p1->p2)
                      n1->burgY[j1],             //
                      n1->burgZ[j1] };           //

   real8  b3[3]   = { n3->burgX[j3],             // b3 = burgers vector (p3->p4)
                      n3->burgY[j3],             //
                      n3->burgZ[j3] };           //

   // compute the forces...

   int    seg12Local = 1;
   int    seg34Local = 1;

   SegSegForce(home,
               p1[0], p1[1], p1[2],
               p2[0], p2[1], p2[2],
               p3[0], p3[1], p3[2],
               p4[0], p4[1], p4[2],
               b1[0], b1[1], b1[2],
               b3[0], b3[1], b3[2],
               a, mu, nu, seg12Local, seg34Local,
               f1, f1+1, f1+2,
               f2, f2+1, f2+2,
               f3, f3+1, f3+2,
               f4, f4+1, f4+2 );
}

// Update_Forces()
//
// Updates the segment nodes and forces for a given segment
//------------------------------------------------------------------------------------------------------------

void Update_Forces
(
   Segment_t *seg,      ///< points to segment to update
   Node_t    *n1 ,      ///< points to node 1
   Node_t    *n2 ,      ///< points to node 2
   real8     *f1 ,      ///< accumulated force on node 1 <xyz>
   real8     *f2        ///< accumulated force on node 2 <xyz>
)
{
   if (seg)
   {
       LOCK(&seg->segLock);             //  lock the segment mutex
                                        //
       seg->forcesSet = 1;              //  indicate that forces have been set for this segment
                                        //
       int j12 = GetArmID(n1,n2);       //  j12 = arm ID of node2 wrt node 1
       int j21 = GetArmID(n2,n1);       //  j21 = arm ID of node1 wrt node 2
                                        //
       AddtoArmForce(n1,j12,f1);        //  add f1 to the arm force
       AddtoArmForce(n2,j21,f2);        //  add f2 to the arm force
                                        //
       V3_ACCUM(seg->f1,f1);            //  add f1 to the node force
       V3_ACCUM(seg->f2,f2);            //  add f2 to the node force
                                        //
       UNLOCK(&seg->segLock);           //  release the segment mutex
    }                                   //
}

// ComputeForces
//
// Compute the forces on an array of segment pairs.
//------------------------------------------------------------------------------------------------------------

void ComputeForces
(
   Home_t        *home,      ///< points to home data structure
   SegmentPair_t *sp  ,      ///< points to an array of segment pairs
   int            np         ///< number of segment pairs
)
{
   if (home && sp && (np>0))
   {
      const Param_t *param  = home->param;      // param = points to home param structure

      const real8    a  = param->rc;            // a  = core radius (b)
      const real8    mu = param->shearModulus;  // mu = shear modulus
      const real8    nu = param->pois;          // nu = poisson ratio
      const real8    lx = param->Lx;            // lx = simulation box size (x)
      const real8    ly = param->Ly;            // ly = simulation box size (y)
      const real8    lz = param->Lz;            // lz = simulation box size (z)

      const real8    sx = ( (param->xBoundType == Periodic) ? 1.0/lx : 0.0 );  // sx = ( PBC active ? 1.0/lx : 0 );
      const real8    sy = ( (param->yBoundType == Periodic) ? 1.0/ly : 0.0 );  // sy = ( PBC active ? 1.0/ly : 0 );
      const real8    sz = ( (param->zBoundType == Periodic) ? 1.0/lz : 0.0 );  // sz = ( PBC active ? 1.0/lz : 0 );

#ifdef SPECTRAL
      Spectral_t *spectral = home->spectral;
      real8       rnei2    = 0.0;

      if (param->FFTenabled)
      {
          rnei2 = spectral->rcnei * spectral->rcnei;
          if (rnei2 == 0.0) return;
      }
#endif

      // loop through all the segment pairs in the array, updating forces on each node...

#ifdef _OPENMP
#pragma omp parallel for shared(sp)
#endif
      for (int i=0; (i<np); i++)
      {
         Segment_t *s1 = sp[i].seg1;      // s1 = points to segment 1
         Segment_t *s2 = sp[i].seg2;      // s2 = points to segment 2

         Node_t    *n1 = s1->node1;       // n1 = points to segment 1 node 1
         Node_t    *n2 = s1->node2;       // n2 = points to segment 1 node 2
         Node_t    *n3 = s2->node1;       // n3 = points to segment 2 node 1
         Node_t    *n4 = s2->node2;       // n4 = points to segment 2 node 2

         real8      f1[3] = { 0.0 };      // f1 = forces on node 1 (init to zero)
         real8      f2[3] = { 0.0 };      // f2 = forces on node 2 (init to zero)
         real8      f3[3] = { 0.0 };      // f3 = forces on node 3 (init to zero)
         real8      f4[3] = { 0.0 };      // f4 = forces on node 4 (init to zero)

#ifdef SPECTRAL
         if (param->FFTenabled)
         {
             real8  dist2 = 0.0;

             MinSegSegDist(home, n1, n2, n3, n4, &dist2);

             if ( (dist2<0) || (dist2>rnei2) ) continue;
         }
#endif

         ComputeForces(home, n1,n2,n3,n4, f1,f2,f3,f4, a,mu,nu, lx,ly,lz, sx,sy,sz);

#ifdef SPECTRAL
         if (param->FFTenabled)
         {
             real8 f1c[3] = { 0.0 };          // f1c = spectral force correction on node 1
             real8 f2c[3] = { 0.0 };          // f2c = spectral force correction on node 2
             real8 f3c[3] = { 0.0 };          // f3c = spectral force correction on node 3
             real8 f4c[3] = { 0.0 };          // f4c = spectral force correction on node 4

             real8 rc     = spectral->rcgrid; // rc  = core radius (for spectral force correction)

             ComputeForces(home, n1,n2,n3,n4, f1c,f2c,f3c,f4c, rc,mu,nu, lx,ly,lz, sx,sy,sz);

             V3_REMOVE(f1,f1c);               // f1 -= f1c
             V3_REMOVE(f2,f2c);               // f2 -= f2c
             V3_REMOVE(f3,f3c);               // f3 -= f3c
             V3_REMOVE(f4,f4c);               // f4 -= f4c
         }
#endif

         if (sp[i].setSeg1Forces) { Update_Forces(s1,n1,n2,f1,f2); }
         if (sp[i].setSeg2Forces) { Update_Forces(s2,n3,n4,f3,f4); }
      }
   }
}
