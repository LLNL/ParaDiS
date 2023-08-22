/*****************************************************************************
 *
 *      Module:         RetroactiveCollision.c
 *      Description:    This module contains various functions used
 *                      for detecting various types of collisions
 *                      (segment/segment, node/node, zipping) and
 *                      dealing with those collisions.  These functions
 *                      are specific to the type 3 collision handling
 *                      which determines whether two segments come within
 *                      a small distance during a time period [t t+deltat].
 *
 *                      Each domain handles local collisions and
 *                      distributes the necessary topological changes
 *                      to remote domains via the same mechanism
 *                      employed by remesh.
 *
 *                      NOTE: Certain strict rules govern what topological
 *                      changes are permitted based on node and segment
 *                      ownership.  See comments at the beginning of
 *                      the module Topology.c for details of the
 *                      rule.  Additional restrictions may be implemented
 *                      for the collision handling; see code below
 *                      for details of any such restrictions.
 *
 *      Included functions:
 *          RetroactiveCollisions()
 *
 *****************************************************************************/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "mpi_portability.h"

#include "Typedefs.h"
#include "Home.h"
#include "Util.h"
#include "Comm.h"
#include "Mobility.h"
#include "V3.h"
#include "M33.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

#ifdef I
#undef I
#endif

static int dbgDom;

/***************************************************************************
 *
 * The functions in this module are a C/C++ implementation of Tom Arsenlis'
 * sphere collision routines.
 * The source of these routines is ../mlab/CollisionCriterion/...
 *
**************************************************************************/

/***************************************************************************
 *
 *     Function :     InSphere
 *
 *     Description :  Given a sphere located at X0 with radius=XR and a point,
 *                    will return
 *                    1 = point is inside  the sphere
 *                    0 = point is outside the sphere
 *      Arguments:
 *                    x0 : sphere center
 *                    xr : sphere radius
 *                    xp : point
 ***************************************************************************/
int InSphere(const real8 *x0, const real8 xr, const real8 *xp)
{
   real8 dr;
   real8 dv[3] = { xp[0]-x0[0],
                   xp[1]-x0[1],
                   xp[2]-x0[2] };

   dr    = DotProduct(dv,dv);

   return ( (dr<(xr*xr)) ? 1 : 0 );
}

/*---------------------------------------------------------------------------
 *
 *      Function :   GrowSphere
 *
 *     Description : Given a sphere located at x0 with radius xr and a point outside
 *                   the sphere, will adjust the sphere location and radius to include
 *                   the new point.
 *
 *-----------------------------------------------------------------------------*/
void GrowSphere (
         real8 *x0,     ///< center of the sphere <xyz>    (returned)
         real8 *xr,     ///< radius of the sphere <scalar> (returned)
   const real8 *xp )    ///< position outside of the sphere <xyz>
{
   if ( !InSphere(x0,(*xr),xp) )
   {
      real8 dx[3] = { xp[0]-x0[0],
                      xp[1]-x0[1],
                      xp[2]-x0[2] };

      real8 ds    = sqrt( dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );
            (*xr)    = 0.5*((*xr)+ds);

            x0[0] = xp[0]-dx[0]/ds*(*xr);
            x0[1] = xp[1]-dx[1]/ds*(*xr);
            x0[2] = xp[2]-dx[2]/ds*(*xr);
   }
}

/**************************************************************************
 *
 *       Function:       FindSphere
 *
 *      Description :    Returns a sphere located at x0 with radius=xr that
 *                       contains all 4 points provided.
 *
 **************************************************************************/
void FindSphere (
         real8 *x0,     ///< center of the sphere <xyz>    (returned)
         real8 *xr,     ///< radius of the sphere <scalar> (returned)
   const real8 *x1,     ///< position 1 <xyz>
   const real8 *x2,     ///< position 2 <xyz>
   const real8 *x3,     ///< position 3 <xyz>
   const real8 *x4 )    ///< position 4 <xyz>
{
   x0[0] = (x1[0]+x2[0])*0.5;
   x0[1] = (x1[1]+x2[1])*0.5;
   x0[2] = (x1[2]+x2[2])*0.5;

   real8 dx[3] = { x1[0]-x2[0],
                   x1[1]-x2[1],
                   x1[2]-x2[2] };

   *xr    = 0.5*sqrt(  dx[0]*dx[0] + dx[1]*dx[1] + dx[2]*dx[2] );

   GrowSphere(x0,xr,x3);
   GrowSphere(x0,xr,x4);
}


/**************************************************************************
 *
 *       Function:       MinPointSegDist()
 *
 * Calculates the minimum distance between a point and a line segment
 *
 ***************************************************************************/

void MinPointSegDist (
         real8 *MinDist2 ,
         real8 *L1       ,
   const real8 *x0       ,
   const real8 *y0       ,
   const real8 *y1       )
{
   // initialize

   real8 diff[3]  = { x0[0]-y0[0],
                      x0[1]-y0[1],
                      x0[2]-y0[2] };


         (*MinDist2) = V3_DOT(diff,diff);
         (*L1)       = 0.0;

   // check other end

   real8 diff2[3] = { x0[0]-y1[0],
                      x0[1]-y1[1],
                      x0[2]-y1[2] };

   real8 Dist2    = V3_DOT(diff2,diff2);

   if (Dist2<(*MinDist2))
   {
      (*MinDist2) = Dist2;
      (*L1)       = 1.0;
   }

   // check middle

   diff2[0] = y0[0]-y1[0];
   diff2[1] = y0[1]-y1[1];
   diff2[2] = y0[2]-y1[2];

   real8 B = V3_DOT(diff2,diff2);

   real8 eps = 1.0e-12;

   if (B>eps)
   {
      real8 L1temp =  (-diff[0]*diff2[0]/B)
                     +(-diff[1]*diff2[1]/B)
                     +(-diff[2]*diff2[2]/B);

      if ( (0.0<L1temp) && (L1temp<1.0) )
      {
         real8 tempv[3] =  { diff[0]+diff2[0]*L1temp,
                             diff[1]+diff2[1]*L1temp,
                             diff[2]+diff2[2]*L1temp };

               Dist2 = V3_DOT(tempv,tempv);

         if ( Dist2<(*MinDist2) )
         {
            (*MinDist2) = Dist2;
            (*L1)       = L1temp;
         }
      }
   }
}


/***************************************************************************
 *
 *       Function:       MinSegSegDist()
 *
 *       Calculates the minimum distance between two line segments.
 *       segment 1 goes from x0 to x1 and segment 2 goes from y0 to y1
 **************************************************************************/
void MinSegSegDist (
         real8 *MinDist2 ,
         real8 *L1       ,
         real8 *L2       ,
   const real8 *x0       ,
   const real8 *x1       ,
   const real8 *y0       ,
   const real8 *y1       )
{

   // intialize

   real8 diff[3];

   V3_SUB (diff,x0,y0);
   (*MinDist2) = V3_DOT(diff,diff); (*L1)=0.0; (*L2)=0.0;

   // Check the ends of the solution space
   // first check for (*L1) = 0

   real8 Dist2  = 0.0;
   real8 L1temp = 0.0;
   real8 L2temp = 0.0;

   MinPointSegDist(&Dist2,&L2temp,x0,y0,y1);

   if ( Dist2<(*MinDist2) ) { (*MinDist2) = Dist2; (*L1) = 0.0; (*L2) = L2temp; }

   // second check for (*L1) = 1

   MinPointSegDist(&Dist2,&L2temp,x1,y0,y1);

   if ( Dist2<(*MinDist2) ) { (*MinDist2) = Dist2; (*L1) = 1.0; (*L2) = L2temp; }

   // third check for (*L2) = 0

   MinPointSegDist(&Dist2,&L1temp,y0,x0,x1);

   if ( Dist2<(*MinDist2) ) { (*MinDist2) = Dist2; (*L2) = 0.0; (*L1) = L1temp; }

   // fourth check for L2 = 1

   MinPointSegDist(&Dist2,&L1temp,y1,x0,x1);

   if ( Dist2<(*MinDist2) ) { (*MinDist2) = Dist2; (*L2) = 1.0; (*L1) = L1temp; }

   // check for an internal solution

   real8 seg1[3]; V3_SUB(seg1,x1,x0);
   real8 seg2[3]; V3_SUB(seg2,y1,y0);

   real8 A =  V3_DOT(seg1,seg1);
   real8 B =  V3_DOT(seg2,seg2);
   real8 C =  V3_DOT(seg1,seg2);

   real8 E =  V3_DOT(seg1,diff);
   real8 D = -V3_DOT(seg2,diff);
   real8 G =  C*C-A*B;

   real8 eps = 1.0e-12; // this is a hard coded small distance number

   if ( fabs(G)>eps )
   {
      L1temp = (A*D+C*E)/G;
      L2temp = (C*D+B*E)/G;

      if ( (0.0<L1temp) && (L1temp<1.0) && (0.0<L2temp) && (L2temp<1.0) )
      {
         real8 diff2[3] = { diff[0]+seg1[0]*L1temp-seg2[0]*L2temp,
                            diff[1]+seg1[1]*L1temp-seg2[1]*L2temp,
                            diff[2]+seg1[2]*L1temp-seg2[2]*L2temp };

         Dist2 = V3_DOT(diff2,diff2);

         if ( Dist2<(*MinDist2) ) { (*MinDist2) = Dist2; (*L2) = L2temp; (*L1) = L1temp; }
      }
   }

}

/*---------------------------------------------------------------------------
 *
 *      Function:     MinSegSegDist
 *
 *-------------------------------------------------------------------------*/
void MinSegSegDist
(
   Home_t *home ,
   Node_t *n1   ,
   Node_t *n2   ,
   Node_t *n3   ,
   Node_t *n4   ,
   real8  *dist2
)
{
   real8 x1[3] = { n1->x, n1->y, n1->z };
   real8 x2[3] = { n2->x, n2->y, n2->z };
   real8 x3[3] = { n3->x, n3->y, n3->z };
   real8 x4[3] = { n4->x, n4->y, n4->z };

   Param_t *param = home->param;

   PBCPOSITION(param, x1[0], x1[1], x1[2], &x2[0], &x2[1], &x2[2]);
   PBCPOSITION(param, x1[0], x1[1], x1[2], &x3[0], &x3[1], &x3[2]);
   PBCPOSITION(param, x3[0], x3[1], x3[2], &x4[0], &x4[1], &x4[2]);

   real8 dv1[3] = { x2[0]-x1[0],
                    x2[1]-x1[1],
                    x2[2]-x1[2] };

   real8 dv2[3] = { x4[0]-x3[0],
                    x4[1]-x3[1],
                    x4[2]-x3[2] };

   if ( V3_DOT(dv1,dv1) < 1.0e-20) { *dist2 = -1.0; return; }
   if ( V3_DOT(dv2,dv2) < 1.0e-20) { *dist2 = -1.0; return; }

   real8 l1=0.0;
   real8 l2=0.0;

   MinSegSegDist(dist2, &l1, &l2, x1, x2, x3, x4);
}

/***************************************************************************
 *
 *       Function:       MindistPtPtInTime()
 *
 * find the minimum distance between two points during a time increment
 * assuming that both points move with a constant velocity during the time
 * increment from t to tau
 ***************************************************************************/
void MindistPtPtInTime (
         real8 *MinDist2 ,
         real8 *tratio   ,
   const real8 *x1t      ,
   const real8 *x1tau    ,
   const real8 *x2t      ,
   const real8 *x2tau    )
{
   real8 diff[3], diff2[3];

   // initialize with points at time t

              V3_SUB (diff,x1t,x2t);
   (*MinDist2) = V3_DOT(diff,diff);
   (*tratio)   = 0.0;

   //  check points at time tau

                 V3_SUB (diff2,x1tau,x2tau);
   real8 Dist2 = V3_DOT(diff2,diff2);

   if (Dist2<(*MinDist2))
   {
      (*MinDist2) = Dist2;
      (*tratio)   = 1.0;
   }

   // check for an internal solution

               V3_SUB(diff2,diff2,diff);
   real8 B   = V3_DOT(diff2,diff2);
   real8 eps = 1.0e-12;

   if (B>eps)
   {
      real8 ttemp = -V3_DOT(diff,diff2)/B;

      if ( (0.0<ttemp) && (ttemp<1.0) )
      {
         real8 diff3[3] = { diff[0]+diff2[0]*ttemp,
                            diff[1]+diff2[1]*ttemp,
                            diff[2]+diff2[2]*ttemp };

         Dist2 = V3_DOT(diff3,diff3);

         if (Dist2<(*MinDist2))
         {
             (*MinDist2) = Dist2;
             (*tratio)   = ttemp;
         }
      }
   }
}


/***************************************************************************
 *
 *      Function:        MinDistPtSegInTime()
 *
 * Find the closest distance between a moving point and a moving line segment
 * during a time increment.
 *
 ****************************************************************************/
void MinDistPtSegInTime (
         real8 *MinDist2 ,
         real8 *L2       ,
         real8 *tratio   ,
   const real8 *x1t      ,
   const real8 *x1tau    ,
   const real8 *x3t      ,
   const real8 *x3tau    ,
   const real8 *x4t      ,
   const real8 *x4tau    )
{
   //  intialize

   real8 diff[3];

              V3_SUB (diff,x1t,x3t);
   (*MinDist2) = V3_DOT(diff,diff);
   (*L2)       = 0.0;
   (*tratio)   = 0.0;

   // assume that the the minimum distance between the point and the segment
   // at tratio = 0 and tratio = 1 have already been included in the colision
   // calculation so don't repeat them here.  If that assumption is not true
   // then extra check are needed to check the boundaries of the box.
   // find the minimum distance between the point and the first endpoint of the
   // segment during time increment

   real8 Dist2=0.0, ttemp=0.0;

   MindistPtPtInTime(&Dist2,&ttemp,x1t,x1tau,x3t,x3tau);

   if ( Dist2<(*MinDist2) )
   {
       (*MinDist2) = Dist2;
       (*L2) = 0.0;
       (*tratio) = ttemp;
   }

   //  find the mininum distance between the point and the second enpoint of the
   //  segment during the time increment

   MindistPtPtInTime(&Dist2,&ttemp,x1t,x1tau,x4t,x4tau);

   if ( Dist2<(*MinDist2) )
   {
      (*MinDist2) = Dist2;
      (*L2) = 1.0;
      (*tratio) = ttemp;
   }

   // check for an internal solution

   real8 L2temp  = (*L2);
          ttemp  = (*tratio);
   real8  L34[3];  V3_SUB(L34,x3t,x4t);

   real8 dL34[3] = { (x3tau[0]-x3t[0])+(x4t[0]-x4tau[0]),
                     (x3tau[1]-x3t[1])+(x4t[1]-x4tau[1]),
                     (x3tau[2]-x3t[2])+(x4t[2]-x4tau[2]) };

   real8 dL13[3] = { (x1tau[0]-x1t[0])+(x3t[0]-x3tau[0]),
                     (x1tau[1]-x1t[1])+(x3t[1]-x3tau[1]),
                     (x1tau[2]-x1t[2])+(x3t[2]-x3tau[2]) };

   real8 A = V3_DOT(diff,diff);
   real8 B = V3_DOT(diff, L34);
   real8 C = V3_DOT(diff,dL13);
   real8 D = V3_DOT(diff,dL34);
   real8 E = V3_DOT( L34, L34);
   real8 F = V3_DOT( L34,dL13);
   real8 G = V3_DOT( L34,dL34);
   real8 H = V3_DOT(dL13,dL13);
   real8 I = V3_DOT(dL13,dL34);
   real8 J = V3_DOT(dL34,dL34);

   // optimize some loop invariants...

         D = (D+F);

   real8 twoG = 2.0*G;
   real8 twoI = 2.0*I;

   // newton iterate...

   for (int i=0, imax=20; (i<imax); i++)
   {
       real8  L2L2 = L2temp*L2temp;

       real8  Err[2];
       {
          Err[0] = C + D*L2temp + G*L2L2 + ttemp*(H + twoI*L2temp + J*L2L2);
          Err[1] = B + E*L2temp + ttemp*(D + twoG*L2temp)  + ttemp*ttemp*(I + J*L2temp);
       }

       real8  Mat[4];
       {
          Mat[0] = H + twoI*L2temp +   J*L2L2;
          Mat[1] = D + twoG*L2temp + 2.0*ttemp*(I + J*L2temp);
          Mat[2] = D + twoG*L2temp + 2.0*ttemp*(I + J*L2temp);
          Mat[3] = E + twoG* ttemp +     ttemp*ttemp*(J);
       }

       real8 detM   = Mat[0]*Mat[3] - Mat[2]*Mat[1];
       real8 Errmag = Err[0]*Err[0] + Err[1]*Err[1];

       // convergence tolerance set to 1 percent in both Lratio and tratio

       real8 eps = 1e-12;

       if ( (fabs(detM)<eps) || (Errmag<1e-6) )
           i = imax;
       else
       {
           ttemp  = ttemp  + (Mat[3]*Err[0] - Mat[1]*Err[1])/detM;
           L2temp = L2temp + (Mat[0]*Err[1] - Mat[2]*Err[0])/detM;
       }
   }

   if ( (0.0<L2temp) && (L2temp<1.0) && (0.0<ttemp) && (ttemp<1.0) )
   {
       real8  L2L2 = L2temp*L2temp;

       Dist2 =               (A + 2.0*B*L2temp + E*L2L2)
               +     2.0*ttemp*(C +   D*L2temp + G*L2L2)
               + ttemp*ttemp*(H + twoI*L2temp + J*L2L2);

       if ( Dist2<(*MinDist2) )
       {
           (*MinDist2) = Dist2;
           (*L2)       = L2temp;
           (*tratio)   = ttemp;
       }
   }
}


/***************************************************************************
 *
 *      Function:        MinDistSegSegInTime()
 *
 *     Description:      Finds the minimum distance between to moving line segments.
 *                       The line segments are one segment going from x1 to x2 and
 *                       the other from x3 to x4. It assumed that the ends of the
 *                       line segments go from their t position to their
 *                       tau position linear in time.
 *
 ***************************************************************************/

void MinDistSegSegInTime (
         real8 *MinDist2 ,
         real8 *L1ratio  ,
         real8 *L2ratio  ,
         real8 *dtratio  ,
   const real8 *x1t      ,
   const real8 *x1tau    ,
   const real8 *x2t      ,
   const real8 *x2tau    ,
   const real8 *x3t      ,
   const real8 *x3tau    ,
   const real8 *x4t      ,
   const real8 *x4tau    )
{
   real8 eps     = 1.0e-12;
   (*dtratio) = 0.0;
   real8 Dist2   = 0.0;
   real8 L1temp  = 0.0;
   real8 L2temp  = 0.0;
   real8 dttemp  = 0.0;

   (*MinDist2) = 0.0;

   // Find the minimum distance between the segments at time t

   MinSegSegDist(MinDist2,L1ratio,L2ratio,x1t,x2t,x3t,x4t);

   // Find the minimum distance between the segments at time tau

   MinSegSegDist(&Dist2,&L1temp,&L2temp,x1tau,x2tau,x3tau,x4tau);

  if ( Dist2 < (*MinDist2) )
   {
      (*MinDist2) = Dist2;
      (*L1ratio)  = L1temp;
      (*L2ratio)  = L2temp;
      (*dtratio)  = 1.0;
   }

   // Check node 1 against line segment 34 during time increment

   MinDistPtSegInTime(&Dist2,&L2temp,&dttemp,x1t,x1tau,x3t,x3tau,x4t,x4tau);

   if ( Dist2<(*MinDist2) )
   {
      (*MinDist2) = Dist2;
      (*L1ratio)  = 0.0;
      (*L2ratio)  = L2temp;
      (*dtratio)  = dttemp;
   }

   //check node 2 against line segment 34 during time increment

   MinDistPtSegInTime(&Dist2,&L2temp,&dttemp,x2t,x2tau,x3t,x3tau,x4t,x4tau);

   if ( Dist2<(*MinDist2) )
   {
      (*MinDist2) = Dist2;
      (*L1ratio)  = 1.0;
      (*L2ratio)  = L2temp;
      (*dtratio)  = dttemp;
   }

   // check node 3 against line segment 12 during time increment

   MinDistPtSegInTime(&Dist2,&L1temp,&dttemp,x3t,x3tau,x1t,x1tau,x2t,x2tau);

   if ( Dist2<(*MinDist2) )
   {
      (*MinDist2) = Dist2;
      (*L1ratio)  = L1temp;
      (*L2ratio)  = 0.0;
      (*dtratio)  = dttemp;
   }

   // check node 4 against line segment 12 during time increment

   MinDistPtSegInTime(&Dist2,&L1temp,&dttemp,x4t,x4tau,x1t,x1tau,x2t,x2tau);

   if ( Dist2<(*MinDist2) )
   {
      (*MinDist2) = Dist2;
      (*L1ratio)  = L1temp;
      (*L2ratio)  = 1.0;
      (*dtratio)  = dttemp;
   }

   // All surface solutions have been investigated now check for an internal solution
   // Use the best surface solution as the initial guess for the iterative solve of the
   // non-linear internal problem.

   L1temp = (*L1ratio);
   L2temp = (*L2ratio);
   dttemp = (*dtratio);

   real8 L13 [3]= { x1t[0]-x3t[0], x1t[1]-x3t[1], x1t[2]-x3t[2] };
   real8 L21 [3]= { x2t[0]-x1t[0], x2t[1]-x1t[1], x2t[2]-x1t[2] };
   real8 L34 [3]= { x3t[0]-x4t[0], x3t[1]-x4t[1], x3t[2]-x4t[2] };

   real8 dL13[3]= { (x1tau[0]-x1t[0])+(x3t[0]-x3tau[0]),
                    (x1tau[1]-x1t[1])+(x3t[1]-x3tau[1]),
                    (x1tau[2]-x1t[2])+(x3t[2]-x3tau[2]) };

   real8 dL21[3]= { (x2tau[0]-x2t[0])+(x1t[0]-x1tau[0]),
                    (x2tau[1]-x2t[1])+(x1t[1]-x1tau[1]),
                    (x2tau[2]-x2t[2])+(x1t[2]-x1tau[2]) };

   real8 dL34[3]= { (x3tau[0]-x3t[0])+(x4t[0]-x4tau[0]),
                    (x3tau[1]-x3t[1])+(x4t[1]-x4tau[1]),
                    (x3tau[2]-x3t[2])+(x4t[2]-x4tau[2]) };

   real8 A = V3_DOT( L13, L13);
   real8 B = V3_DOT( L13, L21);
   real8 C = V3_DOT( L13, L34);
   real8 D = V3_DOT( L13,dL13);
   real8 E = V3_DOT( L13,dL21);
   real8 F = V3_DOT( L13,dL34);
   real8 G = V3_DOT( L21, L21);
   real8 H = V3_DOT( L21, L34);
   real8 I = V3_DOT( L21,dL13);
   real8 J = V3_DOT( L21,dL21);
   real8 K = V3_DOT( L21,dL34);
   real8 L = V3_DOT( L34, L34);
   real8 M = V3_DOT( L34,dL13);
   real8 N = V3_DOT( L34,dL21);
   real8 O = V3_DOT( L34,dL34);
   real8 P = V3_DOT(dL13,dL13);
   real8 Q = V3_DOT(dL13,dL21);
   real8 R = V3_DOT(dL13,dL34);
   real8 S = V3_DOT(dL21,dL21);
   real8 T = V3_DOT(dL21,dL34);
   real8 U = V3_DOT(dL34,dL34);

   // declare and optimize some loop invariants...

         E = (E+I);
         F = (F+M);
         K = (K+N);

   real8 twoQ = 2.0*Q;
   real8 twoR = 2.0*R;
   real8 twoT = 2.0*T;
   real8 twoJ = 2.0*J;
   real8 twoO = 2.0*O;

   // newton iterate...

   for (int i=0, imax=20; (i<imax); i++)
   {
      real8 dttemp2 = dttemp*dttemp;
      real8 L1L1    = L1temp*L1temp;
      real8 L2L2    = L2temp*L2temp;
      real8 L1L2    = L1temp*L2temp;

      real8 Err[3];
      {
         Err[0] =             D+E*L1temp+    F*L2temp+   K*L1L2+J*L1L1+O*L2L2
                  +dttemp*(P+twoQ*L1temp+ twoR*L2temp+twoT*L1L2+S*L1L1+U*L2L2);

         Err[1] = B+H*L2temp+G*L1temp +dttemp*(E+K*L2temp+twoJ*L1temp) +dttemp2*(Q+T*L2temp+S*L1temp);
         Err[2] = C+H*L1temp+L*L2temp +dttemp*(F+K*L1temp+twoO*L2temp) +dttemp2*(R+T*L1temp+U*L2temp);
      }

      real8 Mat[9];
      {
         Mat[0] = P+twoQ*L1temp+twoR*L2temp+twoT*L1L2+S*L1L1+U*L2L2;
         Mat[1] = E+K*L2temp+twoJ*L1temp +2.0*dttemp*(Q+T*L2temp+S*L1temp);
         Mat[2] = F+K*L1temp+twoO*L2temp +2.0*dttemp*(R+T*L1temp+U*L2temp);
         Mat[3] = Mat[1];
         Mat[4] = G+twoJ*dttemp+S*dttemp2;
         Mat[5] = H+K*dttemp+T*dttemp2;
         Mat[6] = Mat[2];
         Mat[7] = Mat[5];
         Mat[8] = L+twoO*dttemp+U*dttemp2;
      }

      real8 Errmag = V3_DOT(Err,Err);

      if ( (fabs(M33_DET(Mat))<eps) || (Errmag<1.0e-6) )
         i=imax;
      else
      {
         real8 MatInv[9];
         real8 correction[3];

         M33_INV(MatInv,Mat);
         V3_M33_V3_MUL(correction,MatInv,Err);

         dttemp=dttemp-correction[0];
         L1temp=L1temp-correction[1];
         L2temp=L2temp-correction[2];
      }
   }

   if ( (0.0<dttemp) && (dttemp<1.0) && (0.0<L1temp) && (L1temp<1.0) && (0.0<L2temp) && (L2temp<1.0) )
   {
      real8 dttemp2 = dttemp*dttemp;
      real8 L1L1    = L1temp*L1temp;
      real8 L2L2    = L2temp*L2temp;
      real8 L1L2    = L1temp*L2temp;

      Dist2 =               A+  2.0*B*L1temp+  2.0*C*L2temp+  2.0*H*L1L2+G*L1L1+L*L2L2
              +2.0*dttemp *(D+      E*L1temp+      F*L2temp+      K*L1L2+J*L1L1+O*L2L2)
              +    dttemp2*(P+   twoQ*L1temp+   twoR*L2temp+   twoT*L1L2+S*L1L1+U*L2L2);

      if (Dist2<(*MinDist2))
      {
          (*MinDist2) = Dist2 ;
          (*dtratio)  = dttemp;
          (*L1ratio)  = L1temp;
          (*L2ratio)  = L2temp;
      }
   }
}


/***************************************************************************
 *
 *      Function:        CollisionCriterion
 *
 *      Description:     This function determines whether two segments come
 *                       within a distance mindist during a time period from
 *                       t to tau.
 *                       If no, then CollisionCriterionIsMet=0, if yes then
 *                       CollisionCriterionIsMet=1, and the L1ratio,
 *                       L2ratio are given for the two segments.
 *                       The interval of time [t, tau] can be chosen to be
 *                       (1) the previous time step interval. In this case, mindist
 *                       is rann
 *                       or
 *                       (2) the next time step interval. In this case, mindist
 *                       is eps.
 *
 *
 *      Arguments:
 *                       L1ratio     : Distance from first segment to collision point
 *                       L2ratio     : Distance from second segment to collision point
 *                       x1t[3]      : Position of node 1 at time t
 *                       x1tau[3]    : Position of node 1 at time tau
 *                       x2t[3]      : Position of node 2 at time t
 *                       x2tau[3]    : Position of node 2 at time tau
 *                       x3t[3]      : Position of node 3 at time t
 *                       x3tau[3]    : Position of node 2 at time tau
 *                       x4t[3]      : Position of node 4 at time t
 *                       x4tau[3]    : Position of node 2 at time tau
 *                       mindist     : small distance within which the two segments
 *                                     meet.
 *
 ***************************************************************************/
int CollisionCriterion(real8 *dist2, real8 *L1ratio, real8 *L2ratio, const real8  mindist,
                              const real8 *x1t, const real8 *x1tau,
                              const real8 *x2t, const real8 *x2tau,
                              const real8 *x3t, const real8 *x3tau,
                              const real8 *x4t, const real8 *x4tau)
{
   int CollisionCriterionIsMet = 0;

   *L1ratio  = 0.0;
   *L2ratio  = 0.0;

   real8 center1[3]={0,0,0},radius1=0.0;
   FindSphere(center1,&radius1,x1t,x2t,x1tau,x2tau);
   real8 center2[3]={0,0,0},radius2=0.0;
   FindSphere(center2,&radius2,x3t,x4t,x3tau,x4tau);

   real8 diff1[3];
   V3_SUB(diff1,center1,center2);
   real8 Dist2   = V3_DOT(diff1,diff1);

   real8 filterdist=radius1+radius2+mindist;
   filterdist=filterdist*filterdist;

   if ( Dist2<filterdist )
   {
       real8 L1=0.0,L2=0.0,dt=0.0;
       MinDistSegSegInTime(&Dist2,&L1,&L2,&dt,x1t,x1tau,x2t,x2tau,x3t,x3tau,x4t,x4tau);
       if (Dist2<(mindist*mindist))
       {
          CollisionCriterionIsMet=1;
          *L1ratio =L1;
          *L2ratio =L2;
       }
   }

   *dist2 = Dist2;
   return(CollisionCriterionIsMet);
}

/*---------------------------------------------------------------------------
 *
 *      Function:       RetroactiveCollisions
 *      Description:    Loop though all native nodes, identify segments
 *                      or nodes that should be collided (or zipped) by
 *                      the current domain and handle all such collisions.
 *
 *                      Note the following restrictions on collisions:
 *
 *                      - A 'pinned' node be deleted during a collision
 *                      - A node which has been the target of a merge
 *                        (this timestep) may not be deleted during any
 *                        other type of collision in this timestep.
 *
 *      All nodes with arms owned by another domain have previously
 *      been marked as non-deletable nodes for this cycle.
 *
 *-------------------------------------------------------------------------*/
void RetroactiveCollisions(Home_t *home)
{
        int     i, j, k, q, arm12, arm34;
        int     thisDomain, splitStatus, mergeStatus;
        int     armAB, retCode=0;
        int     globalOp = 1, didCollision;
        int     close2node1, close2node2, close2node3, close2node4;
        int     splitSeg1, splitSeg2;
        int     cell2Index, nbrCell2Index, nextIndex;
        int     cell2X, cell2Y, cell2Z, cx, cy, cz;
        int     localCollisionCnt;
        int     collisionConditionIsMet;
        real8   mindist, mindist2, dist2 = 0, L1, L2, half,dt;

        real8   x1[3], vx1, vy1, vz1;
        real8   x2[3], vx2, vy2, vz2;
        real8   x3[3], vx3, vy3, vz3;
        real8   x4[3], vx4, vy4, vz4;

        real8   xold1[3], xold2[3], xold3[3], xold4[3];
        real8   xnext1[3], xnext2[3], xnext3[3], xnext4[3];

        real8   vnextx1, vnexty1, vnextz1;
        real8   vnextx2, vnexty2, vnextz2;
        real8   vnextx3, vnexty3, vnextz3;
        real8   vnextx4, vnexty4, vnextz4;

        real8   newPos[3], newVel[3];
        real8   vec1[3], vec2[3];
        real8   p1[3], p2[3], p3[3], p4[3];
        Tag_t   oldTag1, oldTag2;
        Node_t  *node1, *node2, *node3, *node4, *tmpNbr;
        Node_t  *mergenode1, *mergenode2, *targetNode;
        Param_t *param;
        MobArgs_t mobArgs;
#ifdef _ARLFEM
        int     resetSurfaceProperties;
        real8   femSurfaceNorm[3];
#endif


        thisDomain = home->myDomain;
        param      = home->param;

        mindist  = param->rann;
        mindist2 = mindist * mindist;

        dt = param->deltaTT;
        half     = 0.5;

        localCollisionCnt = 0;

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom = -1;
#endif

        TimerStart(home, COLLISION_HANDLING);

/*
 *      Start looping through native nodes looking for segments to collide...
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if (node1->flags & NO_COLLISIONS) continue;

            didCollision = 0;

/*
 *          Loop through all cells neighboring the node.  Only
 *          nodes in these neighboring cells are candidates for
 *          collisions.
 *
 *          NOTE: All nodes are assigned cell's membership prior to
 *          entering collision handling.  However, nodes added during
 *          node separation and collision handling are not given
 *          cell's membership for this cycle.  In most cases this
 *          is not a problem since those nodes have usually been
 *          exempted from subsequent collisions this cycle.  There
 *          are certain circumstances, though, in which the exemptions
 *          have not been set.  If we encounter one such node that
 *          has no cell's membership, skip it for this cycle, or bad
 *          things happen.
 */
            cell2Index = node1->cell2Idx;
            if (cell2Index < 0) {
                continue;
            }

            DecodeCell2Idx(home, cell2Index, &cell2X, &cell2Y, &cell2Z);

            for (cx = cell2X - 1; cx <= cell2X + 1; cx++) {
             for (cy = cell2Y - 1; cy <= cell2Y + 1; cy++) {
              for (cz = cell2Z - 1; cz <= cell2Z + 1; cz++) {

                nbrCell2Index = EncodeCell2Idx(home, cx, cy, cz);

/*
 *              Loop though all nodes in the neighbor cell2
 */
                nextIndex = home->cell2[nbrCell2Index];

                while (nextIndex >= 0) {

                    if (didCollision) break;

                    node3 = home->cell2QentArray[nextIndex].node;
                    nextIndex = home->cell2QentArray[nextIndex].next;

                    if (node3 == (Node_t *)NULL) continue;
                    if (node3->flags & NO_COLLISIONS) continue;

                    if (CollisionNodeOrder(home, &node1->myTag,
                                           &node3->myTag) >= 0) {
                        continue;
                    }

/*
 *                  Loop over all arms of node1.  Skip any arms that
 *                  terminate at node3 (those hinge arms will be dealt
 *                  with later)
 */
                    for (arm12 = 0; arm12 < node1->numNbrs; arm12++) {

                        if (didCollision) break;

                        node2 = GetNodeFromTag(home, node1->nbrTag[arm12]);

                        if (node2 == (Node_t *)NULL) continue;
                        if (node2->flags & NO_COLLISIONS) continue;

                        if (CollisionNodeOrder(home, &node1->myTag,
                                               &node2->myTag) > 0) {
                            continue;
                        }

                        if ((node2->myTag.domainID == node3->myTag.domainID) &&
                            (node2->myTag.index    == node3->myTag.index   )) {
                            continue;
                        }

/*
 *                      Segment node1/node2 may only be used in a collision
 *                      if the segment is owned by the current domain.
 */
                        if (!DomainOwnsSeg(home, OPCLASS_COLLISION,
                                           thisDomain, &node2->myTag)) {
                            continue;
                        }

/*
 *                      Loop over all arms of node3.  Skip any segments that
 *                      terminate at node2 (those hinge arms will be dealt
 *                      with later), and skip any segments not owned by
 *                      node3 -- we'll deal with those segments when we
 *                      hit the node owning the segment
 *
 */
                        for (arm34 = 0; arm34 < node3->numNbrs; arm34++) {

                            if (didCollision) break;

                            node4 = GetNodeFromTag(home, node3->nbrTag[arm34]);

                            if (node4 == (Node_t *)NULL) continue;
                            if (node4->flags & NO_COLLISIONS) continue;

                            if ((node4->myTag.domainID==node2->myTag.domainID)&&
                                (node4->myTag.index   ==node2->myTag.index)) {
                                continue;
                            }

                            if (CollisionNodeOrder(home, &node3->myTag,
                                                   &node4->myTag) > 0) {
                                continue;
                            }

/*
 *                          At this point, segment node3/node4 is owned by
 *                          node3.  If node3 is not native to this domain,
 *                          the segment may not be used in a collision since
 *                          the domain doing to collision must own both
 *                          segments.
 */
                            if (node3->myTag.domainID != thisDomain) {
                                continue;
                            }

/*
 *                          At this point, we have 4 nodes to test for collision.
 *                          Define their coordinates, their previous position and
 *                          velocity as well as estimates their next position and
 *                          velocity.
 */
                            x1[X] = node1->x; x1[Y] = node1->y; x1[Z] = node1->z;
                            x2[X] = node2->x; x2[Y] = node2->y; x2[Z] = node2->z;
                            x3[X] = node3->x; x3[Y] = node3->y; x3[Z] = node3->z;
                            x4[X] = node4->x; x4[Y] = node4->y; x4[Z] = node4->z;

                            vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
                            vx2 = node2->vX; vy2 = node2->vY; vz2 = node2->vZ;
                            vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
                            vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;

                            PBCPOSITION(param, x1[X], x1[Y], x1[Z], &x2[X], &x2[Y], &x2[Z]);
                            PBCPOSITION(param, x1[X], x1[Y], x1[Z], &x3[X], &x3[Y], &x3[Z]);
                            PBCPOSITION(param, x3[X], x3[Y], x3[Z], &x4[X], &x4[Y], &x4[Z]);

/*
 *                          It is possible to have a zero-length segment
 *                          (created by a previous collision).  If we find
 *                          such a segment, do not try to use it in any
 *                          subsequent collisions.
 */
                            vec1[X] = x2[X] - x1[X];
                            vec1[Y] = x2[Y] - x1[Y];
                            vec1[Z] = x2[Z] - x1[Z];

                            if (DotProduct(vec1, vec1) < 1.0e-20) {
                                continue;
                            }

                            vec2[X] = x4[X] - x3[X];
                            vec2[Y] = x4[Y] - x3[Y];
                            vec2[Z] = x4[Z] - x3[Z];

                            if (DotProduct(vec2, vec2) < 1.0e-20) {
                                continue;
                            }


                            p1[X] = x1[X];  p1[Y] = x1[Y];  p1[Z] = x1[Z];
                            p2[X] = x2[X];  p2[Y] = x2[Y];  p2[Z] = x2[Z];
                            p3[X] = x3[X];  p3[Y] = x3[Y];  p3[Z] = x3[Z];
                            p4[X] = x4[X];  p4[Y] = x4[Y];  p4[Z] = x4[Z];


                            // Where the points were at the previous time step
                            xold1[X] = node1->oldx; xold1[Y] = node1->oldy; xold1[Z] = node1->oldz;
                            xold2[X] = node2->oldx; xold2[Y] = node2->oldy; xold2[Z] = node2->oldz;
                            xold3[X] = node3->oldx; xold3[Y] = node3->oldy; xold3[Z] = node3->oldz;
                            xold4[X] = node4->oldx; xold4[Y] = node4->oldy; xold4[Z] = node4->oldz;

                            PBCPOSITION(param,x1[X],x1[Y],x1[Z],&xold1[X],&xold1[Y],&xold1[Z]);
                            PBCPOSITION(param,x2[X],x2[Y],x2[Z],&xold2[X],&xold2[Y],&xold2[Z]);
                            PBCPOSITION(param,x3[X],x3[Y],x3[Z],&xold3[X],&xold3[Y],&xold3[Z]);
                            PBCPOSITION(param,x4[X],x4[Y],x4[Z],&xold4[X],&xold4[Y],&xold4[Z]);

/*
 *                          Check whether to collide segments or not.
 */

                            /* First interval considered : [t-delta t, t] */
                            collisionConditionIsMet = CollisionCriterion(&dist2,&L1,&L2,mindist,
                                                                         xold1,x1,xold2,x2,
                                                                         xold3,x3,xold4,x4);

                           if (!collisionConditionIsMet)
                            {
                               // No collision in the past, look for a possible collision in the future....

                               // Determine where the points will be at the next time step (estimate only).
                               // Approximation : v_{n+1} equal to v_n
                               vnextx1 = vx1; vnexty1 = vy1; vnextz1 = vz1;
                               vnextx2 = vx2; vnexty2 = vy2; vnextz2 = vz2;
                               vnextx3 = vx3; vnexty3 = vy3; vnextz3 = vz3;
                               vnextx4 = vx4; vnexty4 = vy4; vnextz4 = vz4;

                               xnext1[X] = x1[X] + half*dt*(vx1+vnextx1);
                               xnext1[Y] = x1[Y] + half*dt*(vy1+vnexty1);
                               xnext1[Z] = x1[Z] + half*dt*(vz1+vnextz1);

                               xnext2[X] = x2[X] + half*dt*(vx2+vnextx2);
                               xnext2[Y] = x2[Y] + half*dt*(vy2+vnexty2);
                               xnext2[Z] = x2[Z] + half*dt*(vz2+vnextz2);

                               xnext3[X] = x3[X] + half*dt*(vx3+vnextx3);
                               xnext3[Y] = x3[Y] + half*dt*(vy3+vnexty3);
                               xnext3[Z] = x3[Z] + half*dt*(vz3+vnextz3);

                               xnext4[X] = x4[X] + half*dt*(vx4+vnextx4);
                               xnext4[Y] = x4[Y] + half*dt*(vy4+vnexty4);
                               xnext4[Z] = x4[Z] + half*dt*(vz4+vnextz4);

                               PBCPOSITION(param,x1[X],x1[Y],x1[Z],&xnext1[X],&xnext1[Y],&xnext1[Z]);
                               PBCPOSITION(param,x2[X],x2[Y],x2[Z],&xnext2[X],&xnext2[Y],&xnext2[Z]);
                               PBCPOSITION(param,x3[X],x3[Y],x3[Z],&xnext3[X],&xnext3[Y],&xnext3[Z]);
                               PBCPOSITION(param,x4[X],x4[Y],x4[Z],&xnext4[X],&xnext4[Y],&xnext4[Z]);

                               /* Second trial interval considered is : [t, t + delta t] */
                               collisionConditionIsMet = CollisionCriterion(&dist2,&L1,&L2,1.e-6,
                                                                            x1,xnext1,x2,xnext2,
                                                                            x3,xnext3,x4,xnext4);

                            }

                           if (!collisionConditionIsMet) continue;

                            // Find newPos: the position of the collision point.
                            // Need to determine mergenode1 and mergenode2: the nodes to be merged at the newPos.


/*
 *                              Identify the first node to be merged : mergenode1.
 *                              Procedure : If the collision point is close to one of
 *                              the nodal endpoints, use that node, otherwise insert a
 *                              new node in the segment.
 *
 *                              NOTE: The current domain owns node1 but may
 *                              not own node2.  If it does not own node2, we
 *                              cannot allow the collision to use node2
 *                              even if the collision point is close to
 *                              that node.
 */

                                vec1[X] = p1[X] - p2[X];
                                vec1[Y] = p1[Y] - p2[Y];
                                vec1[Z] = p1[Z] - p2[Z];

                                close2node1 = ( ((L1 * L1)         * DotProduct(vec1,vec1)) <mindist2);
                                close2node2 = ( ((1.0-L1)*(1.0-L1) * DotProduct(vec1,vec1)) <mindist2);

                                if ((node2->myTag.domainID != thisDomain) && close2node2) {
                                    continue;
                                }

                                if (close2node1) {
                                    mergenode1 = node1;
                                    splitSeg1 = 0;
                                } else if (close2node2) {
                                    mergenode1 = node2;
                                    splitSeg1 = 0;
                                } else {
                                    splitSeg1 = 1;
                                }

/*
 *                              If we need to add a new node to the first
 *                              segment, do it now.
 */
                                if (splitSeg1) {

                                   splitStatus = FindMergeNodeOnSegment(home,&node1,&node2,
                                                                        x1[X],x1[Y],x1[Z],x2[X],x2[Y],x2[Z],
                                                                        arm12, L1, &mergenode1);

                                   if (splitStatus == SPLIT_FAILED) continue;

/*
 *                                   When we actually do a collision, we flag
 *                                   it so we don't attempt additional changes
 *                                   on either segment.  It's possible, however
 *                                   that we will split the segment but later
 *                                   determine a collision will not be done.
 *                                   Since the original node1/node2 segment
 *                                   has been bisected, we must treat it as
 *                                   if a collision had occurred so we do not
 *                                   attempt to collide the now non-existent
 *                                   node1/node2 segment.
 */
                                     didCollision = 1;

                                }  /* if (splitSeg1) */

/*
 *                              Identify the second node to be merged: mergenode2
 *
 *                              Note: The current domain owns node3 but may not
 *                              own node4.  If it does not own node4, we
 *                              cannot allow the collision to use node4
 *                              even if the collision point is close to
 *                              that node.
 */
                                vec2[X] = p3[X] - p4[X];
                                vec2[Y] = p3[Y] - p4[Y];
                                vec2[Z] = p3[Z] - p4[Z];

                                close2node3 = ((L2 * L2)              * DotProduct(vec2,vec2) <mindist2);
                                close2node4 = (((1.0- L2) * (1.0-L2)) * DotProduct(vec2,vec2) <mindist2);

                                if ((node4->myTag.domainID != thisDomain) && close2node4) {
                                    continue;
                                }


                                if (close2node3) {
                                    mergenode2 = node3;
                                    splitSeg2 = 0;
                                } else if (close2node4) {
                                    mergenode2 = node4;
                                    splitSeg2 = 0;
                                } else {
                                    splitSeg2 = 1;
                                }

/*
 *                              If we need to add a new node to the second
 *                              segment, do it now.
 */
                                if (splitSeg2) {

                                   splitStatus = FindMergeNodeOnSegment(home,&node3,&node4,
                                                                        x3[X],x3[Y],x3[Z],x4[X],x4[Y],x4[Z],
                                                                        arm34, L2, &mergenode2);

                                   if (splitStatus == SPLIT_FAILED) continue;

                                }  /* if (splitSeg2) */

/*
 *                              Two colliding nodes mergenode1 and mergenode2
 *                              have been identified. Do the merge after finding
 *                              the collision point.
 *
 *                              If the distance between the segments was small,
 *                              the segments are already touching. In this case,
 *                              the new position can be determined without using
 *                              the FindCollisionPoint function.
 */

                                if (dist2 < 1.0)
                                {
                                    newPos[X] = (1.0-L1)*x1[X] + L1 * x2[X];
                                    newPos[Y] = (1.0-L1)*x1[Y] + L1 * x2[Y];
                                    newPos[Z] = (1.0-L1)*x1[Z] + L1 * x2[Z];
                                }
                                else
                                {
                                    retCode = FindCollisionPoint(home, mergenode1, mergenode2, &newPos[X], &newPos[Y], &newPos[Z]);

                                    if (!retCode) continue;
                                }

                                FoldBox(param, &newPos[X],&newPos[Y],&newPos[Z]);
/*
 *                             However, if it looks like the
 *                             node will be thrown a significant distance
 *                             during the collision (due to glide con-
 *                             straints), don't do the collision.
 */

                               vec1[X] = newPos[X] - mergenode1->x;
                               vec1[Y] = newPos[Y] - mergenode1->y;
                               vec1[Z] = newPos[Z] - mergenode1->z;

                               vec2[X] = newPos[X] - mergenode2->x;
                               vec2[Y] = newPos[Y] - mergenode2->y;
                               vec2[Z] = newPos[Z] - mergenode2->z;

                               if ((DotProduct(vec1, vec1) > 16 * mindist2) &&
                                   (DotProduct(vec2, vec2) > 16 * mindist2)) {
				 continue;
			       }


#ifdef _ARLFEM
/*
 *                             If colliding 2 surface nodes, we may have to
 *                             adjust the collision point so it too is on the
 *                             surface.
 */
                               resetSurfaceProperties = 0;

                               if (HAS_ANY_OF_CONSTRAINTS(
                                       mergenode1->constraint, SURFACE_NODE) &&
                                   HAS_ANY_OF_CONSTRAINTS(
                                       mergenode2->constraint, SURFACE_NODE)) {

                                   Node_t *seg1Node2, *seg2Node2;

                                   seg1Node2 = (mergenode1 == node1) ?
                                               node2 : node1;
                                   seg2Node2 = (mergenode2 == node3) ?
                                               node4 : node3;

                                   FEM_AdjustCollisionPoint(mergenode1,
                                                            seg1Node2,
                                                            mergenode2,
                                                            seg2Node2,
                                                            newPos,
                                                            femSurfaceNorm);
                                   resetSurfaceProperties = 1;
                               }
#endif  /* ifdef _ARLFEM */

                               newVel[X] = half * (mergenode1->vX +
                                                   mergenode2->vX);
                               newVel[Y] = half * (mergenode1->vY +
                                                   mergenode2->vY);
                               newVel[Z] = half * (mergenode1->vZ +
                                                   mergenode2->vZ);


                               oldTag1 = mergenode1->myTag;
                               oldTag2 = mergenode2->myTag;

                               MergeNode(home, OPCLASS_COLLISION, mergenode1,
                                         mergenode2, newPos, &targetNode,
                                         &mergeStatus, globalOp);

/*
 *                             If the merge did not succeed, go back and
 *                             continue looking for collision candidates.
 */
                               if ((mergeStatus & MERGE_SUCCESS) == 0) {
                                   continue;
                               }
#ifdef _ARLFEM
/*
 *                             Need to explicitly reset surface properties
 *                             after colliding 2 surface nodes.
 */
                               if (resetSurfaceProperties) {
                                   if (targetNode != (Node_t *)NULL) {
                                       VECTOR_COPY(targetNode->surfaceNorm,
                                                   femSurfaceNorm);
                                   }
                               }
#endif  /* ifdef _ARLFEM */

/*
 *                             When debugging, dump some info on topological
 *                             changes taking place and the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                               if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                                   if (targetNode == (Node_t *)NULL) {
                                       printf("  Merge(SegCollision): "
                                              "(%d,%d) and (%d,%d)\n",
                                              oldTag1.domainID, oldTag1.index,
                                              oldTag2.domainID, oldTag2.index);
                                   } else {
                                       printf("  Merge(SegCollision): "
                                              "(%d,%d) and (%d,%d) at "
                                              "(%d,%d)\n",
                                              oldTag1.domainID, oldTag1.index,
                                              oldTag2.domainID, oldTag2.index,
                                              targetNode->myTag.domainID,
                                              targetNode->myTag.index);
                                       PrintNode(targetNode);
                                   }
                               }
#endif

/*
 *                             If the target node exists after the merge,
 *                             reset its velocity, and reset the topological
 *                             change exemptions for the target node; use
 *                             NodeTopologyExemptions() to get the basic
 *                             exemptions, then exempt all the node's arms
 *                             from additional segment collisions this cycle.
 */
                               if (targetNode != (Node_t *)NULL) {
/*
 *                                 If we are enforcing glide planes but
 *                                 allowing some fuzziness in the planes, we
 *                                 also need to recalculate the glide
 *                                 planes for the segments attched to the
 *                                 collision node.
 */
                                   if (param->enforceGlidePlanes &&
                                       param->allowFuzzyGlidePlanes) {
                                       int n;
                                       for (n=0; n<targetNode->numNbrs; n++) {
                                           tmpNbr = GetNodeFromTag(home,
                                                   targetNode->nbrTag[n]);
                                           RecalcSegGlidePlane(home,
                                                               targetNode,
                                                               tmpNbr, 1);
                                       }
                                   }

/*
 *                                 Estimate velocity so mobility function
 *                                 has a reasonable starting point
 */
                                   targetNode->vX = newVel[X];
                                   targetNode->vY = newVel[Y];
                                   targetNode->vZ = newVel[Z];

                                   (void)EvaluateMobility(home, targetNode,
                                                          &mobArgs);

                                   targetNode->flags |= NO_COLLISIONS;
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                                   targetNode->multiNodeLife = 0;
#endif
                               }

                               didCollision = 1;
                               localCollisionCnt++;

                        }  /* Loop over node3 arms */
                    }  /* Loop over node1 arms */
                }  /* while(nextIndex >= 0) */
            }  /* Loop over neighboring cell2s */
           }
          }
        }  /* for (i = 0 ...) */


/*
 *      Now we have to loop for collisions on hinge joints (i.e zipping)
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node1 = home->nodeKeys[i]) == (Node_t *)NULL) continue;
            if (node1->flags & NO_COLLISIONS) continue;

            for (j = 0; j < node1->numNbrs; j++) {

                if (node1->myTag.domainID != node1->nbrTag[j].domainID)
                    continue;

                for (k = j + 1; k < node1->numNbrs; k++) {

                    if (node1->myTag.domainID != node1->nbrTag[k].domainID)
                        continue;

                    node3 = GetNodeFromTag(home, node1->nbrTag[j]);
                    node4 = GetNodeFromTag(home, node1->nbrTag[k]);

                    x1[X] = node1->x; x1[Y] = node1->y; x1[Z] = node1->z;
                    x3[X] = node3->x; x3[Y] = node3->y; x3[Z] = node3->z;
                    x4[X] = node4->x; x4[Y] = node4->y; x4[Z] = node4->z;

                    vx1 = node1->vX; vy1 = node1->vY; vz1 = node1->vZ;
                    vx3 = node3->vX; vy3 = node3->vY; vz3 = node3->vZ;
                    vx4 = node4->vX; vy4 = node4->vY; vz4 = node4->vZ;

                    PBCPOSITION(param, x1[X], x1[Y], x1[Z], &x3[X], &x3[Y], &x3[Z]);
                    PBCPOSITION(param, x3[X], x3[Y], x3[Z], &x4[X], &x4[Y], &x4[Z]);

                    VECTOR_COPY(x2,x1);
                    vx2 = vx1; vy2 = vy1; vz2 = vz1;

/*
 *                  It is possible to have a zero-length segment
 *                  (created by a previous collision).  If we find
 *                  such a segment, do not try to use it in any
 *                  subsequent collisions.
 */
                    vec1[X] = x1[X] - x3[X];
                    vec1[Y] = x1[Y] - x3[Y];
                    vec1[Z] = x1[Z] - x3[Z];

                    if (DotProduct(vec1, vec1) < 1.0e-20) {
                       continue;
                    }

                    vec2[X] = x1[X] - x4[X];
                    vec2[Y] = x1[Y] - x4[Y];
                    vec2[Z] = x1[Z] - x4[Z];

                    if (DotProduct(vec2, vec2) < 1.0e-20) {
                       continue;
                    }

                    p1[X] = x1[X];  p1[Y] = x1[Y];  p1[Z] = x1[Z];
                    p2[X] = x2[X];  p2[Y] = x2[Y];  p2[Z] = x2[Z];
                    p3[X] = x3[X];  p3[Y] = x3[Y];  p3[Z] = x3[Z];
                    p4[X] = x4[X];  p4[Y] = x4[Y];  p4[Z] = x4[Z];

                    // Define where the points were at the previous time step
                    xold1[X] = node1->oldx; xold1[Y] = node1->oldy; xold1[Z] = node1->oldz;
                    xold3[X] = node3->oldx; xold3[Y] = node3->oldy; xold3[Z] = node3->oldz;
                    xold4[X] = node4->oldx; xold4[Y] = node4->oldy; xold4[Z] = node4->oldz;

                    PBCPOSITION(param,x1[X],x1[Y],x1[Z],&xold1[X],&xold1[Y],&xold1[Z]);
                    PBCPOSITION(param,x3[X],x3[Y],x3[Z],&xold3[X],&xold3[Y],&xold3[Z]);
                    PBCPOSITION(param,x4[X],x4[Y],x4[Z],&xold4[X],&xold4[Y],&xold4[Z]);

                    VECTOR_COPY(xold2,xold1);

/*
 *                  Check whether to collide segments or not.
 */

                    /* First interval considered : [t-delta t, t] */
                    collisionConditionIsMet = CollisionCriterion(&dist2,&L1,&L2,mindist,
                                                                 xold1,x1,xold4,x4,
                                                                 xold3,x3,xold3,x3);


                   if (!collisionConditionIsMet)
                    {
                       // No collision in the past, look for a possible collision in the future....

                       // Determine where the points will be at the next time step (estimate only).
                       // Approximation : v_{n+1} equal to v_n
                       vnextx1 = vx1; vnexty1 = vy1; vnextz1 = vz1;
                       vnextx3 = vx3; vnexty3 = vy3; vnextz3 = vz3;
                       vnextx4 = vx4; vnexty4 = vy4; vnextz4 = vz4;

                       xnext1[X] = x1[X] + half*dt*(vx1+vnextx1);
                       xnext1[Y] = x1[Y] + half*dt*(vy1+vnexty1);
                       xnext1[Z] = x1[Z] + half*dt*(vz1+vnextz1);

                       xnext3[X] = x3[X] + half*dt*(vx3+vnextx3);
                       xnext3[Y] = x3[Y] + half*dt*(vy3+vnexty3);
                       xnext3[Z] = x3[Z] + half*dt*(vz3+vnextz3);

                       xnext4[X] = x4[X] + half*dt*(vx4+vnextx4);
                       xnext4[Y] = x4[Y] + half*dt*(vy4+vnexty4);
                       xnext4[Z] = x4[Z] + half*dt*(vz4+vnextz4);

                       PBCPOSITION(param,x1[X],x1[Y],x1[Z],&xnext1[X],&xnext1[Y],&xnext1[Z]);
                       PBCPOSITION(param,x3[X],x3[Y],x3[Z],&xnext3[X],&xnext3[Y],&xnext3[Z]);
                       PBCPOSITION(param,x4[X],x4[Y],x4[Z],&xnext4[X],&xnext4[Y],&xnext4[Z]);

                       VECTOR_COPY(xnext2,x2);
                       vnextx2=vnextx1; vnexty2=vnexty1; vnextz2=vnextz1;

                       /* Second trial interval considered is : [t, t + deltat] */
                       collisionConditionIsMet = CollisionCriterion(&dist2,&L1,&L2,1.e-6,
                                                                    x1,xnext1,x4,xnext4,
                                                                    x3,xnext3,x3,xnext3);

                    }


                    if (!collisionConditionIsMet) continue;

/*
 *                      Node3 is used as one of the collision points. Find
 *                      the collision point on the node1/node4 segment
 */
                        mergenode1 = node3;

                        vec1[X] = x1[X] - x4[X];
                        vec1[Y] = x1[Y] - x4[Y];
                        vec1[Z] = x1[Z] - x4[Z];

                        close2node1 = ( ((L1*L1) * DotProduct(vec1,vec1)) <mindist2);
                        close2node4 = ( ((1.0-L1)*(1.0-L1) * DotProduct(vec1,vec1) )<mindist2);

/*
 *                      If the collision point is close to one of the
 *                      endpoints of the node1/node4 segment, but the
 *                      node in question is not owned by the current
 *                      domain, skip the collision this time around.
 */
                        if ((close2node1&&(node1->myTag.domainID!=thisDomain))||
                            (close2node4&&(node4->myTag.domainID!=thisDomain))){
                            continue;
                        }

                        if (close2node1) {
                             mergenode2 = node1;
                        } else if (close2node4) {
                             mergenode2 = node4;
                        } else {

                           armAB = GetArmID(node1, node4);
                           splitStatus = FindMergeNodeOnSegment(home,&node1,&node4,
                                                                x1[X],x1[Y],x1[Z],x4[X],x4[Y],x4[Z],
                                                                armAB, L1,&mergenode2);

                           if (splitStatus == SPLIT_FAILED) continue;

                        }

/*
 *                      Two colliding nodes mergenode1 and mergenode2
 *                      have been identified. Do the merge after finding
 *                      the collision point.
 */

                        retCode = FindCollisionPoint(home, mergenode1,
                                                     mergenode2, &newPos[X],
                                                     &newPos[Y], &newPos[Z]);

                        if (!retCode) continue;

                        FoldBox(param, &newPos[X], &newPos[Y], &newPos[Z]);

/*
 *                      The merge is going to happen, so there's a few
 *                      more nodes we'll need to mark for force/velocity
 *                      re-evaluation.
 */
                        MarkNodeForceObsolete(home, node1);
                        MarkNodeForceObsolete(home, node3);
                        MarkNodeForceObsolete(home, node4);

                        for (q = 0; q < node1->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node1->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

                        for (q = 0; q < node3->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node3->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

                        for (q = 0; q < node4->numNbrs; q++) {
                            tmpNbr = GetNodeFromTag(home, node4->nbrTag[q]);
                            if (tmpNbr == (Node_t *)NULL) continue;
                            tmpNbr->flags |= NODE_RESET_FORCES;
                        }

/*
 *                      Try the collision
 */
                        oldTag1 = mergenode1->myTag;
                        oldTag2 = mergenode2->myTag;

                        MergeNode(home, OPCLASS_COLLISION, mergenode1,
                                  mergenode2, newPos, &targetNode,
                                  &mergeStatus, globalOp);

                        localCollisionCnt++;

/*
 *                      If the target node exists after the merge, reevaulate
 *                      the node velocity, mark the forces as obsolete,
 *                      and reset the topological change exemptions for
 *                      the target node and exempt all the node's arms
 *                      from subsequent segment collisions this cycle.
 */
                        if (targetNode != (Node_t *)NULL) {

/*
 *                          If we are enforcing glide planes but
 *                          allowing some fuzziness in the planes, we
 *                          also need to recalculate the glide
 *                          planes for the segments attched to the
 *                          collision node.
 */
                            if (param->enforceGlidePlanes &&
                                param->allowFuzzyGlidePlanes) {
                                int n;
                                for (n = 0; n < targetNode->numNbrs; n++) {
                                    tmpNbr = GetNodeFromTag(home,
                                            targetNode->nbrTag[n]);
                                    RecalcSegGlidePlane(home, targetNode,
                                                        tmpNbr, 1);
                                }
                            }

                            (void)EvaluateMobility(home, targetNode,
                                                   &mobArgs);

                            targetNode->flags |= NODE_RESET_FORCES;
                            targetNode->flags |= NO_COLLISIONS;

#ifdef DEBUG_TOPOLOGY_CHANGES
                            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                                printf("  Merge(HingeCollision): "
                                       "(%d,%d) and (%d,%d) at (%d,%d)\n",
                                       oldTag1.domainID, oldTag1.index,
                                       oldTag2.domainID, oldTag2.index,
                                       targetNode->myTag.domainID,
                                       targetNode->myTag.index);
                                PrintNode(targetNode);
                            }
#endif
#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
                            targetNode->multiNodeLife = 0;
#endif
                        }

                        //    }  /* conditions for zip were met */
                }  /* for (k = 0...) */
            }  /* for (j = 0...) */
        }  /* for (i = 0...) */


#ifdef DEBUG_LOG_COLLISIONS
        int  globalCollisionCnt = 0;

#ifdef PARALLEL
        MPI_Reduce(&localCollisionCnt, &globalCollisionCnt, 1, MPI_INT, MPI_SUM,
                   0, MPI_COMM_WORLD);
#else
        globalCollisionCnt = localCollisionCnt;
#endif
        if (home->myDomain == 0) {
            printf("  Collision count = %d\n", globalCollisionCnt);
        }
#endif

        TimerStop(home, COLLISION_HANDLING);

        return;
}
