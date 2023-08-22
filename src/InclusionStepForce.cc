/*************************************************************************************************
 *      Author:      S. Queyreau and S. Aubry
 *
 *      Module:      InclusionStepForce.c
 *
 *      Description: Contains function for calculating the nodal force associated to the change in
 * 		     step energy induced by the shearing of a coherent precipitate by a dislocation.
 *		     It depends upon a unique value of surface step energy.
 *
 *      Includes public functions:
 *          InclusionStepForce()
 *
 *      Includes private functions:
 *          NA
 *
 *************************************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "V3.h"
#include "M33.h"
#include "Matrix.h"

#ifdef ESHELBY
#ifdef ESHELBYSTEP

static void GetdAdx(real8 b[3], real8 D[3], real8 dDdx1[3][3], real8 dDdx2[3][3],
		    real8 Lambda, real8 dLambdadx1[3], real8 dLambdadx2[3],
		    real8 dAdx[2][3])
{
  real8 A[3],A0;
  cross(D,b,A);
  A0 = Normal(A);

  real8 cc1[3][3], cc2[3][3];
  for (int j=0;j<3;j++)
    {
      real8 tmpin[3],tmpout[3];
      tmpin[X] = dDdx1[X][j];
      tmpin[Y] = dDdx1[Y][j];
      tmpin[Z] = dDdx1[Z][j];
      cross(tmpin,b,tmpout);
      cc1[X][j] = tmpout[X];
      cc1[Y][j] = tmpout[Y];
      cc1[Z][j] = tmpout[Z];

      tmpin[X] = dDdx2[X][j];
      tmpin[Y] = dDdx2[Y][j];
      tmpin[Z] = dDdx2[Z][j];
      cross(tmpin,b,tmpout);
      cc2[X][j] = tmpout[X];
      cc2[Y][j] = tmpout[Y];
      cc2[Z][j] = tmpout[Z];
    }

  real8 eps = 1e-10;
  if (A0 > eps)
  {
     real8 cro[3];
     cro[X] = A[X]/A0;
     cro[Y] = A[Y]/A0;
     cro[Z] = A[Z]/A0;

     for (int j=0;j<3;j++)
     {
        dAdx[0][j] = 0.0;
        dAdx[1][j] = 0.0;
        for (int k=0;k<3;k++)
	{
           dAdx[0][j] += Lambda * cro[k] * cc1[k][j];
           dAdx[1][j] += Lambda * cro[k] * cc2[k][j];
	}
     }


     for (int j=0;j<3;j++)
     {
        dAdx[0][j] += dLambdadx1[j] * A0;
        dAdx[1][j] += dLambdadx2[j] * A0;
     }
  }
  else
  {
     for (int j=0;j<3;j++)
     {
        dAdx[0][j] = 0.0;
        dAdx[1][j] = 0.0;
     }
  }
}

/*---------------------------------------------------------------------------
 *
 *      Function:     IntersectionSphereSegmentandDer
 *      Description:  Determine the intersection points between a segment and
 *                    a sphere as well as the derivative of these points wrt
 *                    end points of the segment.
 *
 *      Note :  Before entering this routine, care must be made to ensure
 *              that PBC is taken into account for pos1, pos2, and C before
 *              rotation in the case of an ellipse.
 *      Note:   This routine returns ratio that are not truncated to [0,1]
 *              contrary to the routine IntersectionSphereSegment in
 *              SegPartIntersect.cc
 *
 *      Arguments:
 *          x1 : node starting the segment
 *          x2 : node ending the segment
 *          C  : sphere center
 *          r  : sphere radius
 *          ratio     : fractions of the line segment inside the inclusion
 *          dratio1dx : derivative of the fractions of the line segment where
 *                      the first intersection point is wrt to segments ends
 *                      point 1 and 2.
 *          dratio2dx : derivative of the fractions of the line segment where
 *                      the second intersection point is wrt to segments ends
 *                      point 1 and 2.
 *
 *-------------------------------------------------------------------------*/
static int IntersectionSphereSegmentandDer(real8 x1[3], real8 x2[3],real8 C[3], real8 r,
                                           int *n, real8 ratio[2], real8 dratiodx1[2][3],
                                           real8 dratiodx2[2][3])
{
   int isInter = 0;
   real8 Diff[3], t[3], R[2][3], s[2], sinter[2];
   real8 invL, d2, test, sdiff;

   // Initialization
   (*n) = 0;
   ratio[0] = 0.0;
   ratio[1] = 0.0;

   for (int i = 0; i < 3; i++)
       {
	 dratiodx1[0][i] = 0.0;
	 dratiodx2[1][i] = 0.0;

	 dratiodx1[0][i] = 0.0;
	 dratiodx2[1][i] = 0.0;
       }
	
   Diff[X] = x2[X]-x1[X];
   Diff[Y] = x2[Y]-x1[Y];
   Diff[Z] = x2[Z]-x1[Z];
   invL = 1./Normal(Diff);

   t[X] = Diff[X]*invL;
   t[Y] = Diff[Y]*invL;
   t[Z] = Diff[Z]*invL;

   // Derivative of line direction
   real8 dtdx[3][3];
   for (int i = 0; i < 3; i++)
     for (int j = 0; j < 3; j++)
       dtdx[i][j] = -((i==j) - t[i] * t[j]) * invL;

   R[0][X] = x1[X] - C[X];
   R[0][Y] = x1[Y] - C[Y];
   R[0][Z] = x1[Z] - C[Z];

   R[1][X] = x2[X] - C[X];
   R[1][Y] = x2[Y] - C[Y];
   R[1][Z] = x2[Z] - C[Z];

   s[0] = DotProduct(R[0],t);
   s[1] = DotProduct(R[1],t);
   sdiff = s[1]-s[0];

   real8 v[3];
   v[X] = R[0][X] - s[0]*t[X];
   v[Y] = R[0][Y] - s[0]*t[Y];
   v[Z] = R[0][Z] - s[0]*t[Z];

   d2 = DotProduct(v,v);
   test = r*r-d2;

   if (test <= 0.0) return 0;

   // the segment intersects the particle
   isInter = 1;
   real8  sqrttest = sqrt(test);

   // Derivative of s[0] and s[1]
   real8 tmp0[3], tmp1[3];
   real8 dsdx1[2][3], dsdx2[2][3];
   for (int i = 0; i < 3; i++)
     {
       tmp0[i] = 0.0;
       tmp1[i] = 0.0;
       for (int j = 0; j < 3; j++)
         {
	   tmp0[i] += dtdx[i][j]*R[0][j];
	   tmp1[i] += dtdx[i][j]*R[1][j];
         }

       dsdx1[0][i] = t[i] + tmp0[i];
       dsdx1[1][i] =        tmp1[i];

       dsdx2[0][i] =      - tmp0[i];
       dsdx2[1][i] = t[i] - tmp1[i];
     }

   // Derivative of v
   real8 dvdx[2][3][3];
   for (int i = 0; i < 3; i++)
     for (int j = 0; j < 3; j++)
       {
	 dvdx[0][i][j] = (i==j) - t[i]*dsdx1[0][j] - s[0]*dtdx[i][j];
	 dvdx[1][i][j] =        - t[i]*dsdx2[0][j] + s[0]*dtdx[i][j];
       }

   real8 tmp[3];
   for (int i = 0; i < 3; i++)  tmp[i] = v[i]/sqrttest;


   // Derivative of intersection fractions
   real8 dsinterdx1[2][3], dsinterdx2[2][3];
   for (int i = 0; i < 3; i++)
     {
       dsinterdx1[0][i] = 0.0;
       dsinterdx2[0][i] = 0.0;
       for (int j = 0; j < 3; j++)
	 {
	   dsinterdx1[0][i] += tmp[j]*dvdx[0][j][i];
	   dsinterdx2[0][i] += tmp[j]*dvdx[1][j][i];
	 }
     }

   for (int i = 0; i < 3; i++)
     {
       dsinterdx1[1][i] =-dsinterdx1[0][i];
       dsinterdx2[1][i] =-dsinterdx2[0][i];
     }

   // Define ratio and their derivatives
   sinter[0] = -sqrttest;
   sinter[1] = -sinter[0];

   ratio[0] = (sinter[0] - s[0])/sdiff;
   ratio[1] = (sinter[1] - s[0])/sdiff;

   if (ratio[0] > 0.0 && ratio[0] <= 1.0) (*n)++;
   if (ratio[1] > 0.0 && ratio[1] <= 1.0) (*n)++;

   for (int i = 0; i < 3; i++)
   {
      dratiodx1[0][i] = (dsinterdx1[0][i] - dsdx1[0][i] + ratio[0] * t[i])*invL;
      dratiodx2[0][i] = (dsinterdx2[0][i] - dsdx2[0][i] - ratio[0] * t[i])*invL;

      dratiodx1[1][i] = (dsinterdx1[1][i] - dsdx1[0][i] + ratio[1] * t[i])*invL;
      dratiodx2[1][i] = (dsinterdx2[1][i] - dsdx2[0][i] - ratio[1] * t[i])*invL;
   }

   return isInter;
}





void InclusionStepForce(real8 gamma, real8  C[3], real8 radius[3],
                        real8 rotation[2][3], real8 x1[3], real8 x2[3],
                        real8 b[3], real8 f1[3], real8 f2[3])
{
   // Initialization
   f1[0] = 0.0;  f1[1] = 0.0;  f1[2] = 0.0;
   f2[0] = 0.0;  f2[1] = 0.0;  f2[2] = 0.0;

   if (fabs(gamma) < 1e-10) return;

   real8 v[3],L, invL;
   v[X] = x2[X]-x1[X];
   v[Y] = x2[Y]-x1[Y];
   v[Z] = x2[Z]-x1[Z];
   L = Normal(v);

   if (L < 1e-10) return;

   invL = 1./L;

   real8 bn, bn2;
   bn = Normal(b);
   bn2 = bn * bn;

   real8 bnvec[3];
   bnvec[X] = b[X]/bn;
   bnvec[Y] = b[Y]/bn;
   bnvec[Z] = b[Z]/bn;

   // Determination of Q1 and Q2
   // e is the direction orthogonal to b, in the plane (b,t)
   real8 t[3];
   t[X] = v[X] * invL;
   t[Y] = v[Y] * invL;
   t[Z] = v[Z] * invL;

   // Derivative of line direction
   real8 dtdx1[3][3];
   for (int i = 0; i < 3; i++)
     for (int j = 0; j < 3; j++)
       dtdx1[i][j] = -((i==j) - t[i] * t[j]) *invL;

   real8 e[3];
   real8 tdotb = DotProduct(t,b)/bn2;
   e[X] = t[X] - tdotb * b[X];
   e[Y] = t[Y] - tdotb * b[Y];
   e[Z] = t[Z] - tdotb * b[Z];

   real8 en = Normal(e);

   if (en < 1e-10)
      return;   // when the segment is screw, no step force
   else
   {
      e[X] /= en;
      e[Y] /= en;
      e[Z] /= en;
   }

   // Find the intersection of the two circles
   real8 nvec[3];
   cross(e,b,nvec);
   real8 nvecn = Normal(nvec);
   if (nvecn < 1e-10)
   {
      return;
   }

   real8 bnv = DotProduct(b,nvec);

    // Define the circle (newC, newR) at the intersection of the two spheres
   real8 rho, newR, newC[3];
   real8 X1C[3];
   real8 R = radius[0];

   X1C[X] = C[X]-x1[X];
   X1C[Y] = C[Y]-x1[Y];
   X1C[Z] = C[Z]-x1[Z];
   rho = DotProduct(X1C,nvec);
   newR = sqrt(R*R-rho*rho);

   newC[X] = C[X] - rho*nvec[X];
   newC[Y] = C[Y] - rho*nvec[Y];
   newC[Z] = C[Z] - rho*nvec[Z];

   real8 Dcenter[3];
   Dcenter[X] = newC[X] + 0.5*(b[X] - bnv*nvec[X]);
   Dcenter[Y] = newC[Y] + 0.5*(b[Y] - bnv*nvec[Y]);
   Dcenter[Z] = newC[Z] + 0.5*(b[Z] - bnv*nvec[Z]);

   real8 d;
   d = 2*sqrt(newR*newR - bn2/4.0);
   real8 Q1[3], Q2[3];
   Q1[X] = Dcenter[X] - 0.5*d*e[X];
   Q1[Y] = Dcenter[Y] - 0.5*d*e[Y];
   Q1[Z] = Dcenter[Z] - 0.5*d*e[Z];

   Q2[X] = Dcenter[X] + 0.5*d*e[X];
   Q2[Y] = Dcenter[Y] + 0.5*d*e[Y];
   Q2[Z] = Dcenter[Z] + 0.5*d*e[Z];

   // Define the intersection points P1 and P2. They can be on the line
   // segment but not on the dislocation segment.
   int isInt=0, n=0;
   real8 ratio[2], dratiodx1[2][3], dratiodx2[2][3];
   Matrix23_Zero(dratiodx1);
   Matrix23_Zero(dratiodx2);
   isInt = IntersectionSphereSegmentandDer(x1,x2,C,R,&n,ratio, dratiodx1, dratiodx2);

   if (isInt ==0) return;

   real8 P1[3], P2[3];
   for (int i=0;i<3;i++)
   {
      P1[i] = x1[i] + ratio[0]*v[i];
      P2[i] = x1[i] + ratio[1]*v[i];
   }

   // Derivative of P1 and P2 wrt x1, x2.
   real8 dP1dx1[3][3], dP1dx2[3][3];
   real8 dP2dx1[3][3], dP2dx2[3][3];
   for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
       {
	 dP1dx1[i][j] = (i==j)*(1-ratio[0]) + v[i] * dratiodx1[0][j];
	 dP1dx2[i][j] = (i==j)*(  ratio[0]) + v[i] * dratiodx2[0][j];

	 dP2dx1[i][j] = (i==j)*(1-ratio[1]) + v[i] * dratiodx1[1][j];
	 dP2dx2[i][j] = (i==j)*(  ratio[1]) + v[i] * dratiodx2[1][j];
       }

   // Derivative of Q1 and Q2 wrt x1, x2.
   real8 depdt[3][3];
   for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
        depdt[i][j] = (i==j) - b[i] * b[j]/bn2;

   real8 dotdee[3];
   for (int j=0;j<3;j++)
     {
       real8 tmp[3];
       tmp[X] = depdt[j][X];
       tmp[Y] = depdt[j][Y];
       tmp[Z] = depdt[j][Z];
       dotdee[j] = DotProduct(tmp,e);
     }

   real8 dedt[3][3];
   for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
       {
	 dedt[i][j] = (depdt[i][j] - dotdee[i] * e[j])/en;
       }

   real8 dedx[3][3];
   for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
       {
	 dedx[i][j] = 0.0;
	 for (int k=0;k<3;k++)
	   {
	     dedx[i][j] +=  dedt[i][k]*dtdx1[k][j];
	   }
       }

   // Derivative of n
   real8 tmpin[3], tmpout[3];
   real8 cc1[3][3], cc2[3][3];
   for (int j=0;j<3;j++)
     {
       tmpin[X] = dedx[X][j];
       tmpin[Y] = dedx[Y][j];
       tmpin[Z] = dedx[Z][j];

       cross(tmpin,b,tmpout);
       cc1[X][j] = tmpout[X]/bn;
       cc1[Y][j] = tmpout[Y]/bn;
       cc1[Z][j] = tmpout[Z]/bn;

       cc2[X][j] =-tmpout[X]/bn;
       cc2[Y][j] =-tmpout[Y]/bn;
       cc2[Z][j] =-tmpout[Z]/bn;
     }

   real8 dnormndx1[3], dnormndx2[3];
   for (int j=0;j<3;j++)
    {
      dnormndx1[j] = 0.0;
      dnormndx2[j] = 0.0;
      for (int k=0;k<3;k++)
	{
	  dnormndx1[j] +=  cc1[k][j] * nvec[k];
	  dnormndx2[j] +=  cc2[k][j] * nvec[k];
	}
    }


   real8 cc[3];
   cross(e,bnvec,cc);
   real8 ncc = Normal(cc);
   real8 ncc2 = ncc * ncc;

   real8 dndx1[3][3], dndx2[3][3];
   if (ncc > 1e-10)
   {
      for (int i=0;i<3;i++)
         for (int j=0;j<3;j++)
         {
            dndx1[i][j] = (cc1[i][j] * ncc - cc[i] * dnormndx1[j])/ncc2;
            dndx2[i][j] = (cc2[i][j] * ncc - cc[i] * dnormndx2[j])/ncc2;
         }
   }
   else
   {
      Matrix33_Zero(dndx1);
      Matrix33_Zero(dndx2);
   }

   // Derivative of Dcenter
   real8 dDcenterdx1[3][3], dDcenterdx2[3][3];
   real8 bdotn = DotProduct(b,nvec);
   for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
       {
	 real8 tmp1[3], tmp2[3];
	 tmp1[X] = dndx1[X][j];
	 tmp1[Y] = dndx1[Y][j];
	 tmp1[Z] = dndx1[Z][j];
	
	 tmp2[X] = dndx2[X][j];
	 tmp2[Y] = dndx2[Y][j];
	 tmp2[Z] = dndx2[Z][j];
	 dDcenterdx1[i][j] = -0.5*DotProduct(b,tmp1)*nvec[i] -0.5*bdotn*dndx1[i][j];
	 dDcenterdx2[i][j] = -0.5*DotProduct(b,tmp2)*nvec[i] -0.5*bdotn*dndx2[i][j];
       }


   // Derivative of rho
   real8 drhodx1[3], drhodx2[3];
   for (int j=0;j<3;j++)
       {
	 real8 tmp1[3];
	 tmp1[X] = dndx1[X][j];
	 tmp1[Y] = dndx1[Y][j];
	 tmp1[Z] = dndx1[Z][j];
	 real8 tmp2[3];
	 tmp2[X] = dndx2[X][j];
	 tmp2[Y] = dndx2[Y][j];
	 tmp2[Z] = dndx2[Z][j];
	 drhodx1[j] = -nvec[j] + DotProduct(X1C,tmp1);
	 drhodx2[j] =            DotProduct(X1C,tmp2);
       }

   // Derivative of newR
   real8 dnewRdx1[3], dnewRdx2[3];
   for (int j=0;j<3;j++)
     {
       dnewRdx1[j] = -rho * drhodx1[j] / newR;
       dnewRdx2[j] = -rho * drhodx2[j] / newR;
     }

   // Derivative of newC
   real8 dnewCdx1[3][3], dnewCdx2[3][3];
   for (int i=0;i<3;i++)
      for (int j=0;j<3;j++)
      {
         dnewCdx1[i][j] = -drhodx1[j] * nvec[i] - rho * dndx1[i][j];
         dnewCdx2[i][j] = -drhodx2[j] * nvec[i] - rho * dndx2[i][j];
      }

    // Derivative of Q1Q2
   real8 dddx1[3], dddx2[3];
   for (int j=0;j<3;j++)
     {
       dddx1[j] = 4*newR*dnewRdx1[j]/d;
       dddx2[j] = 4*newR*dnewRdx2[j]/d;
     }

   real8 dQ1Q2dx1[3][3], dQ1Q2dx2[3][3];
   for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
       {
	 dQ1Q2dx1[i][j] = e[i] * dddx1[j] + d * dedx[i][j];
	 dQ1Q2dx2[i][j] = e[i] * dddx2[j] - d * dedx[i][j];
       }

   real8 dQ1dx1[3][3], dQ1dx2[3][3], dQ2dx1[3][3], dQ2dx2[3][3];
   for (int i=0;i<3;i++)
     for (int j=0;j<3;j++)
       {
	 dQ1dx1[i][j] = dDcenterdx1[i][j]-0.5*dQ1Q2dx1[i][j];
	 dQ1dx2[i][j] = dDcenterdx2[i][j]-0.5*dQ1Q2dx2[i][j];
	
	 dQ2dx1[i][j] = dDcenterdx1[i][j]+0.5*dQ1Q2dx1[i][j];
	 dQ2dx2[i][j] = dDcenterdx2[i][j]+0.5*dQ1Q2dx2[i][j];
       }

   // Correction function to remove singularity
   real8 dLambdadx1[3], dLambdadx2[3];
   real8 const0 = 0.75;
   int po = 5;
   real8 nP1P2 = sqrt( (P2[X]-P1[X])*(P2[X]-P1[X]) +
                       (P2[Y]-P1[Y])*(P2[Y]-P1[Y]) +
                       (P2[Z]-P1[Z])*(P2[Z]-P1[Z]) );
   real8 val = pow(const0 *  nP1P2, po) + 1;
   real8 dconst0 = po * pow(const0, po);
   real8 dP = pow(nP1P2,po-2);
   real8 cst = dconst0 * dP /pow(val,2);

   real8 Lambda = 1.0 - 1.0/val;
   for (int j=0;j<3;j++)
     {
       dLambdadx1[j] = 0.0;
       dLambdadx2[j] = 0.0;

       for (int i=0;i<3;i++)
	 {
            dLambdadx1[j] += cst * (P2[i]-P1[i]) * (dP2dx1[i][j] - dP1dx1[i][j]);
            dLambdadx2[j] += cst * (P2[i]-P1[i]) * (dP2dx2[i][j] - dP1dx2[i][j]);
	 }
     }

   // Test the position of P1 and P2 wrt to Q1 and Q2
   real8 DP1[3],DP2[3];
   DP1[X] = P1[X]-Dcenter[X];
   DP1[Y] = P1[Y]-Dcenter[Y];
   DP1[Z] = P1[Z]-Dcenter[Z];

   DP2[X] = P2[X]-Dcenter[X];
   DP2[Y] = P2[Y]-Dcenter[Y];
   DP2[Z] = P2[Z]-Dcenter[Z];

   real8 d1,d2;
   d1 = DotProduct(DP1,b);
   d2 = DotProduct(DP2,b);

   real8 e1,e2;
   e1 = DotProduct(DP1,e);
   e2 = DotProduct(DP2,e);

   // Determine the step force as a function of the position
   // of the points P1 and P2 wrt Q1 and Q2
   real8 dAdx[2][3];
   for(int i=0;i<3;i++)
   {
      dAdx[0][i] = 0.0;
      dAdx[1][i] = 0.0;
   }

   real8 D[3], dDdx1[3][3], dDdx2[3][3];
   if (d1 <= 0.0 && d2 <= 0.0)
   {
      // P2-P1 contribution only
      real8 dA1dx[2][3];
      Matrix23_Zero(dA1dx);
      if (ratio[1] >= 0 && ratio[1] <= 1)
      {
         V3_SUB(D,P2,newC);                     // D     = P2     - newC
         Matrix33_Sub(dDdx1,dP2dx1,dnewCdx1);   // dDdx1 = dP2dx1 - dnewCdx1
         Matrix33_Sub(dDdx2,dP2dx2,dnewCdx2);   // dDdx2 = dP2dx2 - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA1dx);
      }

      real8 dA2dx[2][3];
      Matrix23_Zero(dA2dx);
      if (ratio[0] >= 0 && ratio[0] <= 1)
      {
         V3_SUB(D,newC,P1);                     // D     = newC     - P1
         Matrix33_Sub(dDdx1,dnewCdx1,dP1dx1);   // dDdx1 = dnewCdx1 - dP1dx1
         Matrix33_Sub(dDdx2,dnewCdx2,dP1dx2);   // dDdx2 = dnewCdx2 - dP1dx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA2dx);
      }

      if (e1<=0 && e2 >=0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA1dx[i][j] + dA2dx[i][j]; }
      }
      else if (e1<0 && e2<0)  // gtest checked
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = -dA1dx[i][j] + dA2dx[i][j]; }
      }
      else if (e1>0 && e2>0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA1dx[i][j] - dA2dx[i][j]; }
      }
      else
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA1dx[i][j] + dA2dx[i][j]; }
      }
   }

   if (d1 <= 0.0 && d2 > 0.0)
   {
      // P1Q2 abd P2Q2 contributions

      real8 dA01dx[2][3], dA1dx[2][3];
      Matrix23_Zero(dA01dx);
      Matrix23_Zero(dA1dx);
      if (ratio[1] >= 0 && ratio[1] <= 1)
      {
         V3_SUB(D,Q2,newC);                     // D     = Q2     - newC
         Matrix33_Sub(dDdx1,dQ2dx1,dnewCdx1);   // dDdx1 = dQ2dx1 - dnewCdx1
         Matrix33_Sub(dDdx2,dQ2dx2,dnewCdx2);   // dDdx2 = dQ2dx2 - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA01dx);

         V3_SUB(D,newC,P2);                     // D     = newC     - P2
         Matrix33_Sub(dDdx1,dnewCdx1,dP2dx1);   // dDdx1 = dnewCdx1 - dP2dx1
         Matrix33_Sub(dDdx2,dnewCdx2,dP2dx2);   // dDdx2 = dnewCdx2 - dP2dx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA1dx);
      }

      real8 dA02dx[2][3], dA2dx[2][3];
      Matrix23_Zero(dA02dx);
      Matrix23_Zero(dA2dx);
      if (ratio[0] >= 0 && ratio[0] <= 1)
      {
         V3_SUB(D,Q2,newC);                     // D     = Q2     - newC
         Matrix33_Sub(dDdx1,dQ2dx1,dnewCdx1);   // dDdx1 = dQ2dx1 - dnewCdx1
         Matrix33_Sub(dDdx2,dQ2dx2,dnewCdx2);   // dDdx2 = dQ2dx2 - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA02dx);

         V3_SUB(D,P1,newC);                     // D     = P1     - newC
         Matrix33_Sub(dDdx1,dP1dx1,dnewCdx1);   // dDdx1 = dP1dx1 - dnewCdx1
         Matrix33_Sub(dDdx2,dP1dx2,dnewCdx2);   // dDdx2 = dP1dx2 - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA2dx);
      }

      if (e1 < 0 && e2 > 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] - dA1dx[i][j] + dA2dx[i][j]; }
      }
      else if (e1 < 0 && e2 < 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] + dA1dx[i][j] + dA2dx[i][j]; }
      }
      else if (e1 > 0 && e2 > 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] - dA1dx[i][j] - dA2dx[i][j]; }
      }
      else
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] + dA1dx[i][j] - dA2dx[i][j]; }
      }

   }

   if (d1 > 0.0 && d2 <= 0.0)
   {
      // Q1-P1 and Q1-P2 contribution
      real8 dA01dx[2][3], dA1dx[2][3];
      Matrix23_Zero(dA01dx);
      Matrix23_Zero(dA1dx);
      if (ratio[1] >= 0 && ratio[1] <= 1)
      {
         V3_SUB(D,Q1,newC);                     // D     = Q1     - newC
         Matrix33_Sub(dDdx1,dQ1dx1,dnewCdx1);   // dDdx1 = dQ1dx1 - dnewCdx1
         Matrix33_Sub(dDdx2,dQ1dx2,dnewCdx2);   // dDdx2 = dQ1dx2 - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA01dx);

         V3_SUB(D,newC,P2);                     // D     = newC   - P2
         Matrix33_Sub(dDdx1,dP1dx1,dQ1dx1);     // dDdx1 = dP1dx1 - dQ1dx1
         Matrix33_Sub(dDdx2,dP1dx2,dQ1dx2);     // dDdx2 = dP1dx2 - dQ1dx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA1dx);
      }

      real8 dA02dx[2][3], dA2dx[2][3];
      Matrix23_Zero(dA02dx);
      Matrix23_Zero(dA2dx);
      if (ratio[0] >= 0 && ratio[0] <= 1)
      {
         V3_SUB(D,Q1,newC);                     // D     = Q1       - newC
         Matrix33_Sub(dDdx1,dQ1dx1,dnewCdx1);   // dDdx1 = dQ1dx1   - dnewCdx1
         Matrix33_Sub(dDdx2,dQ1dx2,dnewCdx2);   // dDdx2 = dQ1dx2   - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA02dx);

         V3_SUB(D,P1,newC);                     // D     = P1       - newC
         Matrix33_Sub(dDdx1,dnewCdx1,dP1dx1);   // dDdx1 = dnewCdx1 - dP1dx1
         Matrix33_Sub(dDdx2,dnewCdx2,dP1dx2);   // dDdx2 = dnewCdx2 - dP1dx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA2dx);
      }

      if (e1 < 0 && e2 > 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] + dA1dx[i][j] - dA2dx[i][j]; }

      }
      else if (e1 > 0 && e2 > 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] + dA1dx[i][j] + dA2dx[i][j]; }

      }
      else if (e1 < 0 && e2 < 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] - dA1dx[i][j] - dA2dx[i][j]; }
      }
      else
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = dA01dx[i][j] + dA02dx[i][j] - dA1dx[i][j] + dA2dx[i][j]; }
      }

    }

   if (d1 > 0.0 && d2 > 0.0)
   {
      // Q2-Q1 and P1-P2 contributions
      real8 dA01dx[2][3], dA1dx[2][3];
      Matrix23_Zero(dA01dx);
      Matrix23_Zero(dA1dx);
      if (ratio[0] >= 0 && ratio[0] <= 1)
      {
         V3_SUB(D,newC,Q1);                     // D     = newC     - Q1
         Matrix33_Sub(dDdx1,dnewCdx1,dQ1dx1);   // dDdx1 = dnewCdx1 - dQ1dx1
         Matrix33_Sub(dDdx2,dnewCdx2,dQ1dx2);   // dDdx2 = dnewCdx2 - dQ1dx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA01dx);

         V3_SUB(D,newC,P1);                     // D     = newC     - P1
         Matrix33_Sub(dDdx1,dP1dx1,dnewCdx1);   // dDdx1 = dP1dx1 - dnewCdx1
         Matrix33_Sub(dDdx2,dP1dx2,dnewCdx2);   // dDdx2 = dP1dx2 - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA1dx);
      }

      real8 dA02dx[2][3], dA2dx[2][3];
      Matrix23_Zero(dA02dx);
      Matrix23_Zero(dA2dx);
      if (ratio[1] >= 0 && ratio[1] <= 1)
      {
         V3_SUB(D,Q2,newC);                     // D     = Q2     - newC
         Matrix33_Sub(dDdx1,dQ2dx1,dnewCdx1);   // dDdx1 = dQ2dx1 - dnewCdx1
         Matrix33_Sub(dDdx2,dQ2dx2,dnewCdx2);   // dDdx2 = dQ2dx2 - dnewCdx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA02dx);

         V3_SUB(D,P2,newC);                     // D     = P2     - newC
         Matrix33_Sub(dDdx1,dnewCdx1,dP2dx1);   // dDdx1 = dnewCdx1 - dP2dx1
         Matrix33_Sub(dDdx2,dnewCdx2,dP2dx2);   // dDdx2 = dnewCdx2 - dP2dx2
         GetdAdx(b,D,dDdx1,dDdx2,Lambda,dLambdadx1,dLambdadx2,dA2dx);
      }

      if (e1 < 0 && e2 > 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = 2*dA01dx[i][j] + 2*dA02dx[i][j] -  dA1dx[i][j] -  dA2dx[i][j]; }
      }
      else if (e1 > 0 && e2 > 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = 2*dA01dx[i][j] + 2*dA02dx[i][j] +  dA1dx[i][j] -  dA2dx[i][j]; }
      }
      else if (e1 < 0 && e2 < 0)
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = 2*dA01dx[i][j] + 2*dA02dx[i][j] -  dA1dx[i][j] +  dA2dx[i][j]; }
      }
      else
      {
         for(int i=0;i<2;i++)
         for(int j=0;j<3;j++)
         { dAdx[i][j] = 2*dA01dx[i][j] + 2*dA02dx[i][j] +  dA1dx[i][j] +  dA2dx[i][j]; }
      }

   }

   f1[X] =-gamma * dAdx[0][X];
   f1[Y] =-gamma * dAdx[0][Y];
   f1[Z] =-gamma * dAdx[0][Z];

   f2[X] =-gamma * dAdx[1][X];
   f2[Y] =-gamma * dAdx[1][Y];
   f2[Z] =-gamma * dAdx[1][Z];
}
#endif // ESHELBYSTEP
#endif // ESHELBY
