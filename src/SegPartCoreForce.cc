#ifdef ESHELBYCORE
/***************************************************************************
 *
 *      Author:      Sylvie Aubry  
 *
 *      Date :       Nov, 1 2016
 *
 *      Module:      SegPartCoreForce
 *
 *      Description: This module contains the functions needed for
 *                   correcting the core force associated with different
 *                   elastic constants inside and outside Eshelby particles
 *
 *      Includes public functions:
 *          CoreEnergyAndDerivative
 *          SegPartCoreForceIsotropic
 *
 ***************************************************************************/

#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

/*--------------------------------------------------------------------------
 *
 *      Function:        SegPartCoreForceIsotropic
 *      Description:     Calculate the correction of core force due to a part 
 *                       of the dislocation segment inside the particle
 *
 *      Parameters:
 *          IN:  home                 To get access to Ecore2, NU2, particle 
 *                                    information
 *          IN:  x1,y1,z1, x2,y2,z2   Dislocation segment
 *          IN:  b                    Burgers vector
 *          IN:  Ecore                Core energy
 *          IN:  NU                   Poisson ratio
 *          OUT: f1, f2               Force due to core difference inside/outside particles
 *
 *--------------------------------------------------------------------------*/
static void SegPartCoreForceIsotropic(Home_t *home, real8 C[3], real8 radius[3], 
                                      real8 rotation[2][3], real8 p1[3], real8 p2[3], 
                                      real8 b[3], real8 f1[3], real8 f2[3])
{
   real8 E1, E2, nu1, nu2;
   real8 epsi;
   Param_t *param;

   // Initialization   
   param = home->param;
   
   f1[0] = 0.0;  f1[1] = 0.0;  f1[2] = 0.0;
   f2[0] = 0.0;  f2[1] = 0.0;  f2[2] = 0.0;

   E1 = param->Ecore;
   E2 = param->Ecore2;

   nu1 = param->pois;
   nu2 = param->pois2;

   // If the elastic constants are the same, no need to account for change
   // of core forces. Do nothing.
   if ( (fabs(E1 - E2) < 1e-5) && (fabs(nu2 - nu1) < 1e-5) ) return;

   real8 R0 = radius[0];
   real8 R02 = R0*R0;

   real8 v[3];
   v[X] = p2[X]-p1[X];
   v[Y] = p2[Y]-p1[Y];
   v[Z] = p2[Z]-p1[Z];
   real8 L;
   L = sqrt(v[X]*v[X] + v[Y]*v[Y] + v[Z]*v[Z]);
   if (L <= 0.0) return;

   real8 r[3];
   r[X] = p1[X] - C[X];
   r[Y] = p1[Y] - C[Y];
   r[Z] = p1[Z] - C[Z];
    
   real8 a,e,c;
   a =     -  DotProduct(v,v);
   e =     -2*DotProduct(v,r);
   c = R02 -  DotProduct(r,r);

   // s1 and s2 definitions
   real8 test;
   test = e*e - 4*a*c;

   if (test < 0.0) return;

   real8 s1,s2,sa,sb;
   
   s1 = 0.0;
   s2 = 1.0;
    
   real8 sqrttest, twoa;
   sqrttest = sqrt(test);
   twoa = 2*a;
   sa = (-e + sqrttest )/twoa;
   sb = (-e - sqrttest )/twoa;
   
   if ( (sa>1.0) || (sb < 0.0)) return;

   if (sa > 0.0) s1 = sa;
   if (sb < 1.0) s2 = sb;
   
   real8 Cst,s12,s13,s22,s23,d1,d2,d3;
   Cst = 1.0/3.0;
    
   s12 = s1*s1;
   s13 = s12*s1;
      
   s22 = s2*s2;
   s23 = s22*s2;
      
   d3 = s23-s13;
   d2 = s22-s12;
   d1 = s2 -s1;
    
   real8 Lambda;
   Lambda = 1.0/R02;
    
   real8 t[3], invL = 1.0/L;
   t[X] = v[X]*invL;
   t[Y] = v[Y]*invL;
   t[Z] = v[Z]*invL;

   real8 bs, bs2;
   real8 be[3],be2;
   bs = b[0]*t[0] + b[1]*t[1] + b[2]*t[2];
   bs2 = bs*bs;
   be[0] = b[0]-bs*t[0]; be[1] = b[1]-bs*t[1]; be[2]=b[2]-bs*t[2];
   be2 = (be[0]*be[0]+be[1]*be[1]+be[2]*be[2]);

   real8 C1,C2,D1,D2;
   C1 = E1/(1-nu1);
   C2 = E2/(1-nu2);
   D1 = nu1*C1;
   D2 = nu2*C2;

   real8 bterm;
   bterm = (E1-E2)*bs2+be2*(C1-C2);

   real8 Ldbtermdx1[3];
   real8 tmp = 2*(D1-D2)*bs;
   Ldbtermdx1[X] = tmp*be[X];
   Ldbtermdx1[Y] = tmp*be[Y];
   Ldbtermdx1[Z] = tmp*be[Z];
    
   real8 polyterm;
   polyterm = Cst * a * d3 + 0.5 * e * d2 + c * d1;
         
   real8 TwoCst = 2 * Cst;
   real8 dpolydx1[3],dpolydx2[3];
   dpolydx1[X] = TwoCst * v[X] * d3 - (v[X]-r[X]) * d2 - 2 * r[X] * d1;              
   dpolydx1[Y] = TwoCst * v[Y] * d3 - (v[Y]-r[Y]) * d2 - 2 * r[Y] * d1;              
   dpolydx1[Z] = TwoCst * v[Z] * d3 - (v[Z]-r[Z]) * d2 - 2 * r[Z] * d1;              
   
   dpolydx2[X] =-TwoCst * v[X] * d3 - r[X] * d2;              
   dpolydx2[Y] =-TwoCst * v[Y] * d3 - r[Y] * d2;              
   dpolydx2[Z] =-TwoCst * v[Z] * d3 - r[Z] * d2;              
   
   f1[X] = Lambda*Ldbtermdx1[X] * polyterm + Lambda*bterm * t[X] * polyterm - Lambda*bterm * L * dpolydx1[X];
   f1[Y] = Lambda*Ldbtermdx1[Y] * polyterm + Lambda*bterm * t[Y] * polyterm - Lambda*bterm * L * dpolydx1[Y];
   f1[Z] =-Lambda*Ldbtermdx1[Z] * polyterm + Lambda*bterm * t[Z] * polyterm - Lambda*bterm * L * dpolydx1[Z];
   
   f2[X] =-Lambda*Ldbtermdx1[X] * polyterm - Lambda*bterm * t[X] * polyterm - Lambda*bterm * L * dpolydx2[X];
   f2[Y] =-Lambda*Ldbtermdx1[Y] * polyterm - Lambda*bterm * t[Y] * polyterm - Lambda*bterm * L * dpolydx2[Y];
   f2[Z] = Lambda*Ldbtermdx1[Z] * polyterm - Lambda*bterm * t[Z] * polyterm - Lambda*bterm * L * dpolydx2[Z];
     
#if 0
   printf("\n");
   Print3("x1",p1);
   Print3("x2",p2);
   Print3("C",C);
   Print3("b",b);
   printf("R0=%f;E1=%f; E2=%f; nu1=%f; nu2=%f;\n",R0,E1,E2,nu1,nu2);
   printf("\n");
   Print3("f1",f1);
   Print3("f2",f2);
   exit(0);
#endif   
}



void SegPartCoreForce(Home_t *home, real8 EPos[3], real8 radius[3], 
                      real8 rotation[2][3], real8 p1[3], real8 p2[3], 
                      real8 b[3], real8 f1[3], real8 f2[3])
{
#ifdef ANISOTROPIC
   printf("Core force calculations for Eshelby particles not supported for anisotropic elasticity yet.\n");
#else
   SegPartCoreForceIsotropic(home, EPos, radius, rotation, p1, p2, b, f1, f2);
#endif
}



#endif

