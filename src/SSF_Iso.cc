#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "mpi_portability.h"
#include "cuda_portability.h"

#include "Typedefs.h"
#include "PBC.h"
#include "V3.h"
#include "SSF_Iso.h"
#include "SSF_Iso_Rational.h"

//-----------------------------------------------------------------------------------------------
// This is an implementation of Tom Arsenlis's reduced memory version of the isotropic
// segment-to-segment forces.
//-----------------------------------------------------------------------------------------------

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

// SSF_Segs_Parallel()
//
// Returns whether two segments are parallel (1) or non-parallel (0).
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
int SSF_Segs_Parallel
(
   const real8 *p1   ,  ///< position of node 1 <xyz> (segment 1)
   const real8 *p2   ,  ///< position of node 2 <xyz> (segment 1)
   const real8 *p3   ,  ///< position of node 3 <xyz> (segment 2)
   const real8 *p4   ,  ///< position of node 4 <xyz> (segment 2)
   const real8  ecrit   ///< critical angle for parallelism test
)
{
   real8 v1[3]     = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };              // v1        = (p2-p1)
   real8 v3[3]     = { p4[0]-p3[0], p4[1]-p3[1], p4[2]-p3[2] };              // v3        = (p4-p3)
   real8 v1dot     = (v1[0]*v1[0]) + (v1[1]*v1[1]) + (v1[2]*v1[2]);          // v1dot     = dot(v1,v1)
   real8 v3dot     = (v3[0]*v3[0]) + (v3[1]*v3[1]) + (v3[2]*v3[2]);          // v3dot     = dot(v3,v3)
   real8 oneoverLp = ( (v1dot<1.0e-20) ? 0.0 : 1.0/sqrt(v1dot) );            // oneoverLp = 1.0/||p2-p1||
   real8 oneoverL  = ( (v3dot<1.0e-20) ? 0.0 : 1.0/sqrt(v3dot) );            // oneoverL  = 1.0/||p4-p3||
   real8 tp[3]     = { v1[0]*oneoverLp, v1[1]*oneoverLp, v1[2]*oneoverLp };  // tp        = (p2-p1)/||p2-p1|| (normalized direction vector of nodes p1->p2)
   real8 t [3]     = { v3[0]*oneoverL , v3[1]*oneoverL , v3[2]*oneoverL  };  // t         = (p4-p3)/||p4-p3|| (normalized direction vector of nodes p3->p4)
   real8 c         = t[0]*tp[0] + t[1]*tp[1] + t[2]*tp[2];                   // c         = dot(t,tp)
   real8 c2        = c*c;                                                    // c2        = c*c
         c2        = ( (c2<1.0) ? c2 : 1.0 );                                // c2        = min(c*c,1.0)
   real8 onemc2    = 1.0-c2;                                                 // onemc2    = 1.0-c*c

   return ( (onemc2<ecrit) ? 1 : 0 );
}

//-----------------------------------------------------------------------------------------------
// SSF_Iso_Force_Integrals()
//
// Will compute the force integrals along the Y and Z axis.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Force_Integrals
(
         real8 *Fyintegrals,       ///< points to an array to receive the Y-integrals (y0...y15)
         real8 *Fzintegrals,       ///< points to an array to receive the Z-integrals (z0...z15)
         real8 *y          ,       ///< points to an array to receive the mean yin values
         real8 *z          ,       ///< points to an array to receive the mean zin values
   const real8  a2         ,       ///<     (a*a)
   const real8  d          ,       ///<
   const real8  c          ,       ///<
   const real8  c2         ,       ///<     (c*c)
   const real8  onemc2     ,       ///<   (1-c*c)
   const real8  onemc2inv  ,       ///< 1/(1-c*c)
   const real8 *R          ,       ///<
   const real8 *t          ,       ///<
   const real8 *tp                 ///<
)
{
    real8  a2_d2 = (a2+d*d*onemc2);
    real8  denom = (1.0/sqrt(onemc2*a2_d2));
    real8 cdenom = (1.0+c)*denom;

    real8 Rdt [4] = {  R[0]*t [0] +R[ 1]*t [1] +R[ 2]*t [2],            // Rdt = Dot(R,t)
                       R[3]*t [0] +R[ 4]*t [1] +R[ 5]*t [2],            //
                       R[6]*t [0] +R[ 7]*t [1] +R[ 8]*t [2],            //
                       R[9]*t [0] +R[10]*t [1] +R[11]*t [2] };          //

    real8 Rdtp[4] = {  R[0]*tp[0] +R[ 1]*tp[1] +R[ 2]*tp[2],            // Rdtp = Dot(R,tp)
                       R[3]*tp[0] +R[ 4]*tp[1] +R[ 5]*tp[2],            //
                       R[6]*tp[0] +R[ 7]*tp[1] +R[ 8]*tp[2],            //
                       R[9]*tp[0] +R[10]*tp[1] +R[11]*tp[2] };          //

    real8 yin[4]  = { (Rdt[0]-c*Rdtp[0])*onemc2inv,                     // yin = (Dot(R,t)-c*Dot(R,tp))/(1-c*c)
                      (Rdt[1]-c*Rdtp[1])*onemc2inv,                     //
                      (Rdt[2]-c*Rdtp[2])*onemc2inv,                     //
                      (Rdt[3]-c*Rdtp[3])*onemc2inv };                   //

    real8 zin[4]  = { (Rdtp[0]-c*Rdt[0])*onemc2inv,                     // zin = (Dot(R,tp)-c*Dot(R,t))/(1-c*c)
                      (Rdtp[1]-c*Rdt[1])*onemc2inv,                     //
                      (Rdtp[2]-c*Rdt[2])*onemc2inv,                     //
                      (Rdtp[3]-c*Rdt[3])*onemc2inv };                   //

    real8 Ra[4]   = { sqrt(R[0]*R[0]+R[ 1]*R[ 1]+R[ 2]*R[ 2]+a2),
                      sqrt(R[3]*R[3]+R[ 4]*R[ 4]+R[ 5]*R[ 5]+a2),
                      sqrt(R[6]*R[6]+R[ 7]*R[ 7]+R[ 8]*R[ 8]+a2),
                      sqrt(R[9]*R[9]+R[10]*R[10]+R[11]*R[11]+a2) };

    real8 Ra_Rdot_tp   [4] = { Rdtp[0]-Ra[0], Rdtp[1]-Ra[1], Rdtp[2]-Ra[2], Rdtp[3]-Ra[3] };
    real8 Ra_Rdot_t    [4] = { Rdt [0]-Ra[0], Rdt [1]-Ra[1], Rdt [2]-Ra[2], Rdt [3]-Ra[3] };

    real8 Rainv        [4] = { 1.0/Ra        [0], 1.0/Ra        [1], 1.0/Ra        [2], 1.0/Ra        [3] };
    real8 Ra_Rdot_tinv [4] = { 1.0/Ra_Rdot_t [0], 1.0/Ra_Rdot_t [1], 1.0/Ra_Rdot_t [2], 1.0/Ra_Rdot_t [3] };
    real8 Ra_Rdot_tpinv[4] = { 1.0/Ra_Rdot_tp[0], 1.0/Ra_Rdot_tp[1], 1.0/Ra_Rdot_tp[2], 1.0/Ra_Rdot_tp[3] };

    real8 temp         [4] = { atan((Ra[0]-yin[0]-zin[0])*cdenom),
                               atan((Ra[1]-yin[1]-zin[1])*cdenom),
                               atan((Ra[2]-yin[2]-zin[2])*cdenom),
                               atan((Ra[3]-yin[3]-zin[3])*cdenom) };

                         y[0] = 0.5*(yin[0]+yin[1]);
                         y[1] = 0.5*(yin[2]+yin[3]);
                         z[0] = 0.5*(zin[0]+zin[2]);
                         z[1] = 0.5*(zin[1]+zin[3]);

    real8 log_Ra_Rdot_tpbi[2] = { -log(Ra_Rdot_tp[1]*Ra_Rdot_tpinv[0]),
                                  -log(Ra_Rdot_tp[3]*Ra_Rdot_tpinv[2])  };

    real8 log_Ra_Rdot_tbi [2] = { -log(Ra_Rdot_t [2]*Ra_Rdot_tinv [0]),
                                  -log(Ra_Rdot_t [3]*Ra_Rdot_tinv [1])  };

    real8 Ra2_R_tpinvbi[2]    = { Rainv[1]*Ra_Rdot_tpinv[1]-Rainv[0]*Ra_Rdot_tpinv[0],
                                  Rainv[3]*Ra_Rdot_tpinv[3]-Rainv[2]*Ra_Rdot_tpinv[2] };

    real8 Ra2_R_tinvbi[2]     = { Rainv[2]*Ra_Rdot_tinv[2]-Rainv[0]*Ra_Rdot_tinv[0],
                                  Rainv[3]*Ra_Rdot_tinv[3]-Rainv[1]*Ra_Rdot_tinv[1] };

    real8  yRa2_R_tpinvbi[2]  = { y[0]* Ra2_R_tpinvbi[0], y[1]* Ra2_R_tpinvbi[1] };
    real8 y2Ra2_R_tpinvbi[2]  = { y[0]*yRa2_R_tpinvbi[0], y[1]*yRa2_R_tpinvbi[1] };
    real8  zRa2_R_tinvbi [2]  = { z[0]* Ra2_R_tinvbi [0], z[1]* Ra2_R_tinvbi [1] };
    real8 z2Ra2_R_tinvbi [2]  = { z[0]*zRa2_R_tinvbi [0], z[1]*zRa2_R_tinvbi [1] };

    real8 commonf025bi[2]     = { c*yRa2_R_tpinvbi[0] - Rainv[1]+Rainv[0],
                                  c*yRa2_R_tpinvbi[1] - Rainv[3]+Rainv[2] };

    real8 commonf205bi[2]     = { c*zRa2_R_tinvbi[0] - Rainv[2]+Rainv[0],
                                  c*zRa2_R_tinvbi[1] - Rainv[3]+Rainv[1] };

    real8 commonf305bi[2]     = { log_Ra_Rdot_tbi[0] - c2*z2Ra2_R_tinvbi[0] -(yin[2]-c *zin[2])*Rainv[2]+(yin[0]-c*zin[0])*Rainv[0],
                                  log_Ra_Rdot_tbi[1] - c2*z2Ra2_R_tinvbi[1] -(yin[3]-c *zin[3])*Rainv[3]+(yin[1]-c*zin[1])*Rainv[1] };

    real8 commonf035bi[2]     = { log_Ra_Rdot_tpbi[0] - c2*y2Ra2_R_tpinvbi[0] -(zin[1]-c *yin[1])*Rainv[1]+(zin[0]-c*yin[0])*Rainv[0],
                                  log_Ra_Rdot_tpbi[1] - c2*y2Ra2_R_tpinvbi[1] -(zin[3]-c *yin[3])*Rainv[3]+(zin[2]-c*yin[2])*Rainv[2] };

    real8 Rasum               = Ra   [3]-Ra   [2]-Ra   [1]+Ra   [0];
    real8 Rainvsum            = Rainv[3]-Rainv[2]-Rainv[1]+Rainv[0];

    real8  log_Ra_Rdot_tpsum  = log_Ra_Rdot_tpbi[1]-log_Ra_Rdot_tpbi[0];
    real8 ylog_Ra_Rdot_tpsum  = y[1]*log_Ra_Rdot_tpbi[1]-y[0]*log_Ra_Rdot_tpbi[0];

    real8  log_Ra_Rdot_tsum   = log_Ra_Rdot_tbi[1]-log_Ra_Rdot_tbi[0];
    real8 zlog_Ra_Rdot_tsum   = z[1]*log_Ra_Rdot_tbi[1]-z[0]*log_Ra_Rdot_tbi[0];

    real8   Ra2_R_tpinvsum    =   Ra2_R_tpinvbi[1]-  Ra2_R_tpinvbi[0];
    real8  yRa2_R_tpinvsum    =  yRa2_R_tpinvbi[1]- yRa2_R_tpinvbi[0];
    real8 y2Ra2_R_tpinvsum    = y2Ra2_R_tpinvbi[1]-y2Ra2_R_tpinvbi[0];

    real8   Ra2_R_tinvsum     =   Ra2_R_tinvbi[1]-  Ra2_R_tinvbi[0];
    real8  zRa2_R_tinvsum     =  zRa2_R_tinvbi[1]- zRa2_R_tinvbi[0];
    real8 z2Ra2_R_tinvsum     = z2Ra2_R_tinvbi[1]-z2Ra2_R_tinvbi[0];

    real8  ycommonf025sum     =      y[1]*commonf025bi[1] -      y[0]*commonf025bi[0];
    real8 y2commonf025sum     = y[1]*y[1]*commonf025bi[1] - y[0]*y[0]*commonf025bi[0];

    real8  zcommonf205sum     =      z[1]*commonf205bi[1] -      z[0]*commonf205bi[0];
    real8 z2commonf205sum     = z[1]*z[1]*commonf205bi[1] - z[0]*z[0]*commonf205bi[0];

    real8  commonf305sum      =      commonf305bi[1]-     commonf305bi[0];
    real8 zcommonf305sum      = z[1]*commonf305bi[1]-z[0]*commonf305bi[0];

    real8  commonf035sum      =      commonf035bi[1]-     commonf035bi[0];
    real8 ycommonf035sum      = y[1]*commonf035bi[1]-y[0]*commonf035bi[0];

    real8 f_003sum            = -2.0*denom*(temp[3]-temp[2]-temp[1]+temp[0]);
    real8 adf_003sum          = a2_d2*f_003sum;
    real8 commonf225sum       = f_003sum-c*Rainvsum;
    real8 commonf223sum       = (c*Rasum-adf_003sum)*onemc2inv;

    real8 f_103sum            = ( c*log_Ra_Rdot_tsum  - log_Ra_Rdot_tpsum )*onemc2inv;
    real8 f_013sum            = ( c*log_Ra_Rdot_tpsum - log_Ra_Rdot_tsum  )*onemc2inv;
    real8 f_113sum            = ( c*adf_003sum - Rasum )*onemc2inv;
    real8 tf_113sum           = 2.0*f_113sum;
    real8 f_203sum            = zlog_Ra_Rdot_tsum  + commonf223sum;
    real8 f_023sum            = ylog_Ra_Rdot_tpsum + commonf223sum;

    real8 f_005sum            = ( f_003sum - yRa2_R_tpinvsum - zRa2_R_tinvsum )/a2_d2;
    real8 f_105sum            = ( Ra2_R_tpinvsum - c*Ra2_R_tinvsum )*onemc2inv;
    real8 f_015sum            = ( Ra2_R_tinvsum  - c*Ra2_R_tpinvsum)*onemc2inv;
    real8 f_115sum            = ( Rainvsum - c*( yRa2_R_tpinvsum + zRa2_R_tinvsum + f_003sum ))*onemc2inv;
    real8 f_205sum            = ( yRa2_R_tpinvsum + c2*zRa2_R_tinvsum  + commonf225sum )*onemc2inv;
    real8 f_025sum            = ( zRa2_R_tinvsum  + c2*yRa2_R_tpinvsum + commonf225sum )*onemc2inv;
    real8 f_215sum            = ( f_013sum - ycommonf025sum + c*(zcommonf205sum-f_103sum) )*onemc2inv;
    real8 f_125sum            = ( f_103sum - zcommonf205sum + c*(ycommonf025sum - f_013sum) )*onemc2inv;
    real8 f_225sum            = ( f_203sum - zcommonf305sum + c*( y2commonf025sum - tf_113sum) )*onemc2inv;
    real8 f_305sum            = (y2Ra2_R_tpinvsum + c*commonf305sum + 2.0*f_103sum)*onemc2inv;
    real8 f_035sum            = (z2Ra2_R_tinvsum  + c*commonf035sum + 2.0*f_013sum)* onemc2inv;
    real8 f_315sum            = (tf_113sum - y2commonf025sum + c*(zcommonf305sum - f_203sum))*onemc2inv;
    real8 f_135sum            = (tf_113sum - z2commonf205sum + c*(ycommonf035sum - f_023sum))*onemc2inv;

    Fyintegrals[ 0] = f_003sum;
    Fyintegrals[ 1] = f_103sum;
    Fyintegrals[ 2] = f_013sum;
    Fyintegrals[ 3] = f_113sum;
    Fyintegrals[ 4] = f_203sum;
    Fyintegrals[ 5] = f_005sum;
    Fyintegrals[ 6] = f_105sum;
    Fyintegrals[ 7] = f_015sum;
    Fyintegrals[ 8] = f_115sum;
    Fyintegrals[ 9] = f_205sum;
    Fyintegrals[10] = f_025sum;
    Fyintegrals[11] = f_215sum;
    Fyintegrals[12] = f_125sum;
    Fyintegrals[13] = f_225sum;
    Fyintegrals[14] = f_305sum;
    Fyintegrals[15] = f_315sum;

    Fzintegrals[ 0] = f_003sum;
    Fzintegrals[ 1] = f_013sum;
    Fzintegrals[ 2] = f_103sum;
    Fzintegrals[ 3] = f_113sum;
    Fzintegrals[ 4] = f_023sum;
    Fzintegrals[ 5] = f_005sum;
    Fzintegrals[ 6] = f_015sum;
    Fzintegrals[ 7] = f_105sum;
    Fzintegrals[ 8] = f_115sum;
    Fzintegrals[ 9] = f_025sum;
    Fzintegrals[10] = f_205sum;
    Fzintegrals[11] = f_125sum;
    Fzintegrals[12] = f_215sum;
    Fzintegrals[13] = f_225sum;
    Fzintegrals[14] = f_035sum;
    Fzintegrals[15] = f_135sum;
}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Forces
(
         real8 *f3        ,     ///< resulting force on node 3 <xyz>
         real8 *f4        ,     ///< resulting force on node 4 <xyz>
   const real8 *Fintegrals,     ///< integrals array
   const real8 *t         ,     ///<
   const real8 *tp        ,     ///<
   const real8 *txtp      ,     ///< Cross(t,tp)
   const real8 *b         ,     ///< burgers vector <xyz>
   const real8 *bp        ,     ///< burgers vector <xyz>
   const real8  c         ,
   const real8  d         ,
   const real8  m4p       ,
   const real8  m4pn      ,
   const real8  a2m4pn    ,
   const real8  a2m8p     ,
   const real8 *y         ,
   const real8  oneoverL
)
{
   real8 m4pd       = m4p   *d;
   real8 m4pnd      = m4pn  *d;
   real8 m4pnd2     = m4pnd *d;
   real8 m4pnd3     = m4pnd2*d;
   real8 a2m4pnd    = a2m4pn*d;
   real8 a2m8pd     = a2m8p *d;

   real8 bxt[3]     = { b [1]*t [2]-b [2]*t [1],
                        b [2]*t [0]-b [0]*t [2],
                        b [0]*t [1]-b [1]*t [0] };

   real8 bpxtp[3]   = { bp[1]*tp[2]-bp[2]*tp[1],
                        bp[2]*tp[0]-bp[0]*tp[2],
                        bp[0]*tp[1]-bp[1]*tp[0] };

   real8 tdb        =  t    [0]*b [0]+t    [1]*b [1]+t    [2]*b [2];
   real8 tdbp       =  t    [0]*bp[0]+t    [1]*bp[1]+t    [2]*bp[2];
   real8 tpdb       =  tp   [0]*b [0]+tp   [1]*b [1]+tp   [2]*b [2];
   real8 tpdbp      =  tp   [0]*bp[0]+tp   [1]*bp[1]+tp   [2]*bp[2];
   real8 txtpdb     =  txtp [0]*b [0]+txtp [1]*b [1]+txtp [2]*b [2];
   real8 txbpdtp    = -txtp [0]*bp[0]-txtp [1]*bp[1]-txtp [2]*bp[2];
   real8 bpxtpdb    =  bpxtp[0]*b [0]+bpxtp[1]*b [1]+bpxtp[2]*b [2];
   real8 txbpdb     =  bxt  [0]*bp[0]+bxt  [1]*bp[1]+bxt  [2]*bp[2];

   real8 txtpxt [3] = { (tp[0]-c   *t[0]), (tp[1]-c   *t[1]), (tp[2]-c   *t[2]) };
   real8 txbpxt [3] = { (bp[0]-tdbp*t[0]), (bp[1]-tdbp*t[1]), (bp[2]-tdbp*t[2]) };
   real8 bpxtpxt[3] = { (tdbp*tp[0]-c*bp[0]), (tdbp*tp[1]-c*bp[1]), (tdbp*tp[2]-c*bp[2]) };
   real8 txtpxbpdtp = tdbp-tpdbp*c;
   real8 txtpxbpdb  = tdbp*tpdb-tpdbp*tdb;

   real8 temp1      = tdbp*tpdb+txtpxbpdb;
   real8 temp2      = 0.0;
   real8 I00a[3]    = { -temp1*txtp[0],
                        -temp1*txtp[1],
                        -temp1*txtp[2] };
   real8 I00b[3]    = { (bxt[0]*txtpxbpdtp),
                        (bxt[1]*txtpxbpdtp),
                        (bxt[2]*txtpxbpdtp) };

         temp1      = (m4pnd*txtpdb);
         temp2      = (m4pnd*bpxtpdb);
   real8 I_003[3]   = { m4pd*I00a[0] - m4pnd*I00b[0] +  temp1*bpxtpxt[0] +  temp2*txtpxt[0],
                        m4pd*I00a[1] - m4pnd*I00b[1] +  temp1*bpxtpxt[1] +  temp2*txtpxt[1],
                        m4pd*I00a[2] - m4pnd*I00b[2] +  temp1*bpxtpxt[2] +  temp2*txtpxt[2] };

         temp1      = (m4pnd3*txtpxbpdtp*txtpdb);
   real8 I_005[3]   = { a2m8pd*I00a[0] - a2m4pnd*I00b[0] - temp1*txtpxt[0],
                        a2m8pd*I00a[1] - a2m4pnd*I00b[1] - temp1*txtpxt[1],
                        a2m8pd*I00a[2] - a2m4pnd*I00b[2] - temp1*txtpxt[2] };

   real8 I10a[3]    = { txbpxt[0]*tpdb - txtp[0]*txbpdb,
                        txbpxt[1]*tpdb - txtp[1]*txbpdb,
                        txbpxt[2]*tpdb - txtp[2]*txbpdb };

   real8 I10b[3]    = { bxt[0]*txbpdtp,
                        bxt[1]*txbpdtp,
                        bxt[2]*txbpdtp };

         temp1      = (m4pn*tdb);
   real8 I_103[3]   = { temp1*bpxtpxt[0] + m4p*I10a[0] -  m4pn*I10b[0],
                        temp1*bpxtpxt[1] + m4p*I10a[1] -  m4pn*I10b[1],
                        temp1*bpxtpxt[2] + m4p*I10a[2] -  m4pn*I10b[2] };

         temp1      = m4pnd2*(txbpdtp*txtpdb+txtpxbpdtp*tdb);
   real8 I_105[3]   = { a2m8p*I10a[0] - a2m4pn*I10b[0] - temp1*txtpxt[0],
                        a2m8p*I10a[1] - a2m4pn*I10b[1] - temp1*txtpxt[1],
                        a2m8p*I10a[2] - a2m4pn*I10b[2] - temp1*txtpxt[2] };

   real8 I01a[3]    = { txtp[0]*bpxtpdb - bpxtpxt[0]*tpdb,
                        txtp[1]*bpxtpdb - bpxtpxt[1]*tpdb,
                        txtp[2]*bpxtpdb - bpxtpxt[2]*tpdb };

         temp1      = (m4pn*tpdb);
         temp2      = (m4pn*bpxtpdb);
   real8 I_013[3]   = { m4p*I01a[0] + temp1*bpxtpxt[0] - temp2*txtp[0],
                        m4p*I01a[1] + temp1*bpxtpxt[1] - temp2*txtp[1],
                        m4p*I01a[2] + temp1*bpxtpxt[2] - temp2*txtp[2] };

         temp1      = (m4pnd2*txtpxbpdtp*tpdb);
         temp2      = (m4pnd2*txtpxbpdtp*txtpdb);
   real8 I_015[3]   = { a2m8p*I01a[0] - temp1*txtpxt[0] + temp2*txtp[0],
                        a2m8p*I01a[1] - temp1*txtpxt[1] + temp2*txtp[1],
                        a2m8p*I01a[2] - temp1*txtpxt[2] + temp2*txtp[2] };

         temp1      = (m4pnd*txbpdtp*tdb);
   real8 I_205[3]   = { -temp1*txtpxt[0],
                        -temp1*txtpxt[1],
                        -temp1*txtpxt[2] };

         temp1      =  (m4pnd*txtpxbpdtp*tpdb) ;
   real8 I_025[3]   =  { temp1*txtp[0],
                         temp1*txtp[1],
                         temp1*txtp[2] };

         temp1      = (m4pnd*(txtpxbpdtp*tdb+txbpdtp*txtpdb));
         temp2      = (m4pnd*txbpdtp*tpdb);
   real8 I_115[3]   = { temp1*txtp[0] - temp2*txtpxt[0],
                        temp1*txtp[1] - temp2*txtpxt[1],
                        temp1*txtp[2] - temp2*txtpxt[2] };

         temp1      = (m4pn*txbpdtp*tdb);
   real8 I_215[3]   = { temp1*txtp[0],
                        temp1*txtp[1],
                        temp1*txtp[2] };

         temp1      = (m4pn*txbpdtp*tpdb);
   real8 I_125[3]   = { temp1*txtp[0],
                        temp1*txtp[1],
                        temp1*txtp[2] };

    if (f3)
    {
       real8 Fint_003 = y[1]*Fintegrals[ 0]-Fintegrals[ 1];
       real8 Fint_103 = y[1]*Fintegrals[ 1]-Fintegrals[ 4];
       real8 Fint_013 = y[1]*Fintegrals[ 2]-Fintegrals[ 3];
       real8 Fint_005 = y[1]*Fintegrals[ 5]-Fintegrals[ 6];
       real8 Fint_105 = y[1]*Fintegrals[ 6]-Fintegrals[ 9];
       real8 Fint_015 = y[1]*Fintegrals[ 7]-Fintegrals[ 8];
       real8 Fint_115 = y[1]*Fintegrals[ 8]-Fintegrals[11];
       real8 Fint_205 = y[1]*Fintegrals[ 9]-Fintegrals[14];
       real8 Fint_025 = y[1]*Fintegrals[10]-Fintegrals[12];
       real8 Fint_215 = y[1]*Fintegrals[11]-Fintegrals[15];
       real8 Fint_125 = y[1]*Fintegrals[12]-Fintegrals[13];

       f3[0]  = I_003[0]*Fint_003;
       f3[0] += I_103[0]*Fint_103;
       f3[0] += I_013[0]*Fint_013;
       f3[0] += I_005[0]*Fint_005;
       f3[0] += I_105[0]*Fint_105;
       f3[0] += I_015[0]*Fint_015;
       f3[0] += I_115[0]*Fint_115;
       f3[0] += I_205[0]*Fint_205;
       f3[0] += I_025[0]*Fint_025;
       f3[0] += I_215[0]*Fint_215;
       f3[0] += I_125[0]*Fint_125;
       f3[0] *= oneoverL;

       f3[1]  = I_003[1]*Fint_003;
       f3[1] += I_103[1]*Fint_103;
       f3[1] += I_013[1]*Fint_013;
       f3[1] += I_005[1]*Fint_005;
       f3[1] += I_105[1]*Fint_105;
       f3[1] += I_015[1]*Fint_015;
       f3[1] += I_115[1]*Fint_115;
       f3[1] += I_205[1]*Fint_205;
       f3[1] += I_025[1]*Fint_025;
       f3[1] += I_215[1]*Fint_215;
       f3[1] += I_125[1]*Fint_125;
       f3[1] *= oneoverL;

       f3[2]  = I_003[2]*Fint_003;
       f3[2] += I_103[2]*Fint_103;
       f3[2] += I_013[2]*Fint_013;
       f3[2] += I_005[2]*Fint_005;
       f3[2] += I_105[2]*Fint_105;
       f3[2] += I_015[2]*Fint_015;
       f3[2] += I_115[2]*Fint_115;
       f3[2] += I_205[2]*Fint_205;
       f3[2] += I_025[2]*Fint_025;
       f3[2] += I_215[2]*Fint_215;
       f3[2] += I_125[2]*Fint_125;
       f3[2] *= oneoverL;
    }

    if (f4)
    {
       real8 Fint_003 = Fintegrals[ 1]-y[0]*Fintegrals[ 0];
       real8 Fint_103 = Fintegrals[ 4]-y[0]*Fintegrals[ 1];
       real8 Fint_013 = Fintegrals[ 3]-y[0]*Fintegrals[ 2];
       real8 Fint_005 = Fintegrals[ 6]-y[0]*Fintegrals[ 5];
       real8 Fint_105 = Fintegrals[ 9]-y[0]*Fintegrals[ 6];
       real8 Fint_015 = Fintegrals[ 8]-y[0]*Fintegrals[ 7];
       real8 Fint_115 = Fintegrals[11]-y[0]*Fintegrals[ 8];
       real8 Fint_205 = Fintegrals[14]-y[0]*Fintegrals[ 9];
       real8 Fint_025 = Fintegrals[12]-y[0]*Fintegrals[10];
       real8 Fint_215 = Fintegrals[15]-y[0]*Fintegrals[11];
       real8 Fint_125 = Fintegrals[13]-y[0]*Fintegrals[12];

       f4[0]  = I_003[0]*Fint_003;
       f4[0] += I_103[0]*Fint_103;
       f4[0] += I_013[0]*Fint_013;
       f4[0] += I_005[0]*Fint_005;
       f4[0] += I_105[0]*Fint_105;
       f4[0] += I_015[0]*Fint_015;
       f4[0] += I_115[0]*Fint_115;
       f4[0] += I_205[0]*Fint_205;
       f4[0] += I_025[0]*Fint_025;
       f4[0] += I_215[0]*Fint_215;
       f4[0] += I_125[0]*Fint_125;
       f4[0] *= oneoverL;

       f4[1]  = I_003[1]*Fint_003;
       f4[1] += I_103[1]*Fint_103;
       f4[1] += I_013[1]*Fint_013;
       f4[1] += I_005[1]*Fint_005;
       f4[1] += I_105[1]*Fint_105;
       f4[1] += I_015[1]*Fint_015;
       f4[1] += I_115[1]*Fint_115;
       f4[1] += I_205[1]*Fint_205;
       f4[1] += I_025[1]*Fint_025;
       f4[1] += I_215[1]*Fint_215;
       f4[1] += I_125[1]*Fint_125;
       f4[1] *= oneoverL;

       f4[2]  = I_003[2]*Fint_003;
       f4[2] += I_103[2]*Fint_103;
       f4[2] += I_013[2]*Fint_013;
       f4[2] += I_005[2]*Fint_005;
       f4[2] += I_105[2]*Fint_105;
       f4[2] += I_015[2]*Fint_015;
       f4[2] += I_115[2]*Fint_115;
       f4[2] += I_205[2]*Fint_205;
       f4[2] += I_025[2]*Fint_025;
       f4[2] += I_215[2]*Fint_215;
       f4[2] += I_125[2]*Fint_125;
       f4[2] *= oneoverL;
    }
}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Special_Force_Integrals
(
         real8 *Fintegrals,
   const real8  a2_d2     ,
   const real8 *y         ,
   const real8 *z
)
{
   real8 ypz         [4] = { y[0]+z[0], y[1]+z[1], y[2]+z[2], y[3]+z[3] };
   real8 ymz         [4] = { y[0]-z[0], y[1]-z[1], y[2]-z[2], y[3]-z[3] };
   real8 Ra          [4] = {  sqrt(a2_d2+ypz[0]*ypz[0]),
                              sqrt(a2_d2+ypz[1]*ypz[1]),
                              sqrt(a2_d2+ypz[2]*ypz[2]),
                              sqrt(a2_d2+ypz[3]*ypz[3]) };
   real8 Rainv       [4] = { 1.0/Ra[0], 1.0/Ra[1], 1.0/Ra[2], 1.0/Ra[3] };
   real8 a2d2inv         =   1.0/a2_d2;
   real8 Rapypz      [4] = { Ra[0]+ypz[0], Ra[1]+ypz[1], Ra[2]+ypz[2], Ra[3]+ypz[3] };

   real8 Rainvbi     [2] = { Rainv[0]-Rainv[1],
                             Rainv[3]-Rainv[2] };

   real8 ypzRainvbi  [2] = { ypz[0]*Rainv[0]-ypz[2]*Rainv[2],
                             ypz[3]*Rainv[3]-ypz[1]*Rainv[1] };

   real8 Log_Ra_ypzbi[2] = { log(Rapypz[0]/Rapypz[2]),
                             log(Rapypz[3]/Rapypz[1]) };

   real8 Rasum           = Ra[0]-Ra[1]-Ra[2]+Ra[3];
   real8 ymzRasum        = ymz[0]*Ra[0]-ymz[1]*Ra[1] -ymz[2]*Ra[2]+ymz[3]*Ra[3];
   real8 Rainvsum        = Rainvbi[0]+Rainvbi[1];
   real8 Log_Ra_ypzsum   = Log_Ra_ypzbi[0]+Log_Ra_ypzbi[1];

   real8 f_003sum        = Rasum*a2d2inv;
   real8 f_103sum        = -0.5*(Log_Ra_ypzsum - a2d2inv*ymzRasum);
   real8 f_113sum        = -Log_Ra_ypzsum;
   real8 f_213sum        = z[0]*Log_Ra_ypzbi[0]+z[1]*Log_Ra_ypzbi[1]-Rasum;
   real8 f_005sum        = a2d2inv*(2*a2d2inv*Rasum - Rainvsum);
   real8 f_105sum        = a2d2inv*(a2d2inv*ymzRasum- (y[0]*Rainvbi[0]+y[2]*Rainvbi[1]));
   real8 f_115sum        = -a2d2inv*(ypzRainvbi[0]+ypzRainvbi[1]);
   real8 f_215sum        = Rainvsum +a2d2inv*(z[0]*ypzRainvbi[0]+z[1]*ypzRainvbi[1]);

   Fintegrals[0] = f_003sum;
   Fintegrals[1] = f_103sum;
   Fintegrals[2] = f_113sum;
   Fintegrals[3] = f_213sum;
   Fintegrals[4] = f_005sum;
   Fintegrals[5] = f_105sum;
   Fintegrals[6] = f_115sum;
   Fintegrals[7] = f_215sum;
}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Special_Forces
(
         real8  *f1        ,
         real8  *f2        ,
   const real8  *Fintegrals,
   const real8  *t         ,
   const real8  *b         ,
   const real8  *bp        ,
   const real8  *nd        ,
   const real8   m4p       ,
   const real8   m4pn      ,
   const real8   a2m4pn    ,
   const real8   a2m8p     ,
   const real8  *y         ,
   const real8   oneoverL
)
{
   real8 tdb          = (t[0]*b [0]) + (t[1]*b [1]) + (t[2]*b [2]);
   real8 tdbp         = (t[0]*bp[0]) + (t[1]*bp[1]) + (t[2]*bp[2]);
   real8 tdbv [3]     = { tdb , tdb , tdb  };
   real8 tdbpv[3]     = { tdbp, tdbp, tdbp };
   real8 nddb         = (nd[0]*b[0]) + (nd[1]*b[1]) + (nd[2]*b[2]);
   real8 nddbv[3]     = { nddb, nddb, nddb, };
   real8 bxt  [3]     = { b [1]*t[2]-b [2]*t[1], b [2]*t[0]-b [0]*t[2], b [0]*t[1]-b [1]*t[0] };
   real8 bpxt [3]     = { bp[1]*t[2]-bp[2]*t[1], bp[2]*t[0]-bp[0]*t[2], bp[0]*t[1]-bp[1]*t[0] };
   real8 ndxt [3]     = { nd[1]*t[2]-nd[2]*t[1], nd[2]*t[0]-nd[0]*t[2], nd[0]*t[1]-nd[1]*t[0] };
   real8 bpxtdb       = (bpxt[0]*b [0]) + (bpxt[1]*b [1]) + (bpxt[2]*b [2]);
   real8 bpxtdnd      = (bpxt[0]*nd[0]) + (bpxt[1]*nd[1]) + (bpxt[2]*nd[2]);
   real8 bpxtdndv[3]  = { bpxtdnd, bpxtdnd, bpxtdnd };
   real8 bpxtxt  [3]  = { tdbpv[0]*t[0]-bp[0],
                          tdbpv[1]*t[1]-bp[1],
                          tdbpv[2]*t[2]-bp[2] };

   real8 I_003[3]     = { m4pn*(nddbv[0]*bpxtxt[0] + bpxtdb*ndxt[0] - bpxtdndv[0]*bxt[0]) - m4p*tdbv[0]*tdbpv[0]*nd[0],
                          m4pn*(nddbv[1]*bpxtxt[1] + bpxtdb*ndxt[1] - bpxtdndv[1]*bxt[1]) - m4p*tdbv[1]*tdbpv[1]*nd[1],
                          m4pn*(nddbv[2]*bpxtxt[2] + bpxtdb*ndxt[2] - bpxtdndv[2]*bxt[2]) - m4p*tdbv[2]*tdbpv[2]*nd[2] };

   real8 I_113[3]     = { ((m4pn-m4p)*tdbv[0]*bpxtxt[0]),
                          ((m4pn-m4p)*tdbv[1]*bpxtxt[1]),
                          ((m4pn-m4p)*tdbv[2]*bpxtxt[2]) };

   real8 I_005[3]     = { -a2m8p*tdbv[0]*tdbpv[0]*nd[0] - a2m4pn*bpxtdndv[0]*bxt[0] - m4pn*bpxtdndv[0]*nddbv[0]*ndxt[0],
                          -a2m8p*tdbv[1]*tdbpv[1]*nd[1] - a2m4pn*bpxtdndv[1]*bxt[1] - m4pn*bpxtdndv[1]*nddbv[1]*ndxt[1],
                          -a2m8p*tdbv[2]*tdbpv[2]*nd[2] - a2m4pn*bpxtdndv[2]*bxt[2] - m4pn*bpxtdndv[2]*nddbv[2]*ndxt[2] };

   real8 I_115[3]     = { -a2m8p*tdbv[0]*bpxtxt[0] - m4pn*bpxtdndv[0]*tdbv[0]*ndxt[0],
                          -a2m8p*tdbv[1]*bpxtxt[1] - m4pn*bpxtdndv[1]*tdbv[1]*ndxt[1],
                          -a2m8p*tdbv[2]*bpxtxt[2] - m4pn*bpxtdndv[2]*tdbv[2]*ndxt[2] };

   if (f1)
   {
      real8 Fint_003  = y[2]*Fintegrals[0]-Fintegrals[1];
      real8 Fint_113  = y[2]*Fintegrals[2]-Fintegrals[3];
      real8 Fint_005  = y[2]*Fintegrals[4]-Fintegrals[5];
      real8 Fint_115  = y[2]*Fintegrals[6]-Fintegrals[7];

            f1[0]     = (I_003[0]*Fint_003 + I_113[0]*Fint_113 + I_005[0]*Fint_005 + I_115[0]*Fint_115)*oneoverL;
            f1[1]     = (I_003[1]*Fint_003 + I_113[1]*Fint_113 + I_005[1]*Fint_005 + I_115[1]*Fint_115)*oneoverL;
            f1[2]     = (I_003[2]*Fint_003 + I_113[2]*Fint_113 + I_005[2]*Fint_005 + I_115[2]*Fint_115)*oneoverL;
   }

   if (f2)
   {
      real8 Fint_003  = Fintegrals[1]-y[0]*Fintegrals[0];
      real8 Fint_113  = Fintegrals[3]-y[0]*Fintegrals[2];
      real8 Fint_005  = Fintegrals[5]-y[0]*Fintegrals[4];
      real8 Fint_115  = Fintegrals[7]-y[0]*Fintegrals[6];

            f2[0]     = (I_003[0]*Fint_003 + I_113[0]*Fint_113 + I_005[0]*Fint_005 + I_115[0]*Fint_115)*oneoverL;
            f2[1]     = (I_003[1]*Fint_003 + I_113[1]*Fint_113 + I_005[1]*Fint_005 + I_115[1]*Fint_115)*oneoverL;
            f2[2]     = (I_003[2]*Fint_003 + I_113[2]*Fint_113 + I_005[2]*Fint_005 + I_115[2]*Fint_115)*oneoverL;
   }

}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Special_Remote_Node_Force_Half
(
         real8 *f3   ,
         real8 *f4   ,
   const real8 *x1   ,
   const real8 *x2   ,
   const real8 *x3   ,
   const real8 *x4   ,
   const real8 *b1   ,
   const real8 *b3   ,
   const real8  a    ,
   const real8  mu   ,
   const real8  nu   ,
   const real8  ecrit
)
{
   // compute constants of integration...

   const real8 a2     = a*a         ;
   const real8 m4p    = 0.25*mu/M_PI;
   const real8 m4pn   = m4p/(1.0-nu);
   const real8 a2m4pn = a2*m4pn     ;
   const real8 a2m8p  = a2*m4p/2.0  ;

   real8  eps        = 1.0e-10;
   real8  tanthetac  = sqrt((ecrit*1.01)/(1.0-ecrit*(1.01)));

   real8  vec1[3]    = { x4[0]-x3[0], x4[1]-x3[1], x4[2]-x3[2] };
   real8  vec2[3]    = { x2[0]-x1[0], x2[1]-x1[1], x2[2]-x1[2] };
   real8  x34bar[3]  = { 0.5*(x3[0]+x4[0]), 0.5*(x3[1]+x4[1]), 0.5*(x3[2]+x4[2]) };
   real8  x12bar[3]  = { 0.5*(x1[0]+x2[0]), 0.5*(x1[1]+x2[1]), 0.5*(x1[2]+x2[2]) };

   real8  L          = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1] + vec1[2]*vec1[2]);
   real8  oneoverL   = 1.0/L;
   real8  t[3]       = { vec1[0]*oneoverL, vec1[1]*oneoverL, vec1[2]*oneoverL };
   real8  Lp         = vec2[0]*t[0] + vec2[1]*t[1] + vec2[2]*t[2];
   real8  halfLp     = 0.5*Lp;
   real8  x1mod[3]   = { x12bar[0]-halfLp*t[0], x12bar[1]-halfLp*t[1], x12bar[2]-halfLp*t[2] };
   real8  x2mod[3]   = { x12bar[0]+halfLp*t[0], x12bar[1]+halfLp*t[1], x12bar[2]+halfLp*t[2] };

   real8  R[3]       = { x34bar[0]-x12bar[0], x34bar[1]-x12bar[1], x34bar[2]-x12bar[2] };

   real8  Rdt        = R[0]*t[0] + R[1]*t[1] + R[2]*t[2];
   real8  nd[3]      = { R[0]-Rdt*t[0], R[1]-Rdt*t[1], R[2]-Rdt*t[2] };

   real8  d2         = nd[0]*nd[0] + nd[1]*nd[1] + nd[2]*nd[2];
   real8  a2_d2      = a2+d2;

   real8  r4         = x4[0]*t[0] + x4[1]*t[1] + x4[2]*t[2];
   real8  r3         = x3[0]*t[0] + x3[1]*t[1] + x3[2]*t[2];

   real8  s2         = x2mod[0]*t[0] + x2mod[1]*t[1] + x2mod[2]*t[2];
   real8  s1         = x1mod[0]*t[0] + x1mod[1]*t[1] + x1mod[2]*t[2];

   real8  y[4]       = {  r3,  r3,  r4,  r4 };
   real8  z[4]       = { -s1, -s2, -s1, -s2 };

   real8  Fintegrals[8] = { 0,0,0,0, 0,0,0,0 };
   real8  f3p[3]        = { 0,0,0 };
   real8  f4p[3]        = { 0,0,0 };

   SSF_Iso_Special_Force_Integrals(Fintegrals,a2_d2,y,z);
   SSF_Iso_Special_Forces         (f3p,f4p,Fintegrals,t,b3,b1,nd,m4p,m4pn,a2m4pn,a2m8p,y,oneoverL);

   // Quadratic Approximation

   real8  Diff[3]     = { vec2[0]-Lp*t[0], vec2[1]-Lp*t[1], vec2[2]-Lp*t[2] };
   real8  DiffMag2    = Diff[0]*Diff[0] + Diff[1]*Diff[1] + Diff[2]*Diff[2];

   real8  oneoverLp   = 1.0/Lp;
   real8  DiffMag     = sqrt(DiffMag2);

   real8  temp=0.0, temp1=0.0;

   if ( (DiffMag2*oneoverLp*oneoverLp) > eps )
      temp = 1.0/DiffMag;
   else
   {
      Diff[0] = -t[0]*t[0];
      Diff[1] = -t[0]*t[1];
      Diff[2] = -t[0]*t[2];

      Diff[0] = Diff[0]+1.0;
      temp1   = Diff[0]*Diff[0] + Diff[1]*Diff[1] + Diff[2]*Diff[2];

      if ( (temp1*oneoverLp*oneoverLp) < eps )
      {
         Diff[0] = -t[1]*t[0];
         Diff[1] = -t[1]*t[1];
         Diff[2] = -t[1]*t[2];

         Diff[1] = Diff[1] + 1.0;
         temp1   = Diff[0]*Diff[0] + Diff[1]*Diff[1] + Diff[2]*Diff[2];
      }

      temp = 1.0/sqrt(temp1);
   }

   real8  tp[3]  = { Diff[0]*temp, Diff[1]*temp, Diff[2]*temp };
          temp   = halfLp*tanthetac;

   real8  x1a[3] = { x1mod[0]-temp*tp[0], x1mod[1]-temp*tp[1], x1mod[2]-temp*tp[2] };
   real8  x2a[3] = { x2mod[0]+temp*tp[0], x2mod[1]+temp*tp[1], x2mod[2]+temp*tp[2] };
   real8  x1b[3] = { x1mod[0]+temp*tp[0], x1mod[1]+temp*tp[1], x1mod[2]+temp*tp[2] };
   real8  x2b[3] = { x2mod[0]-temp*tp[0], x2mod[1]-temp*tp[1], x2mod[2]-temp*tp[2] };

   real8  f3a[3] = { 0.0, 0.0, 0.0 };
   real8  f4a[3] = { 0.0, 0.0, 0.0 };
   real8  f3b[3] = { 0.0, 0.0, 0.0 };
   real8  f4b[3] = { 0.0, 0.0, 0.0 };

   SSF_Iso_Non_Parallel(0,0,f3a,f4a, x1a,x2a,x3,x4, b1,b3, a,mu,nu,ecrit);
   SSF_Iso_Non_Parallel(0,0,f3b,f4b, x1b,x2b,x3,x4, b1,b3, a,mu,nu,ecrit);

   real8  shape  = 0.25*DiffMag/temp;
   real8  shapea = (2.0*shape+1.0)*shape;
   real8  shapeb = (2.0*shape-1.0)*shape;
   real8  shapec = (1.0-2.0*shape)*(1.0+2.0*shape);

   if (f3) { f3[0] = shapea*f3a[0] + shapeb*f3b[0] + shapec*f3p[0];	
             f3[1] = shapea*f3a[1] + shapeb*f3b[1] + shapec*f3p[1];	
             f3[2] = shapea*f3a[2] + shapeb*f3b[2] + shapec*f3p[2]; }

   if (f4) { f4[0] = shapea*f4a[0] + shapeb*f4b[0] + shapec*f4p[0];
             f4[1] = shapea*f4a[1] + shapeb*f4b[1] + shapec*f4p[1];
             f4[2] = shapea*f4a[2] + shapeb*f4b[2] + shapec*f4p[2]; }
}

// SSF_Iso_Non_Parallel()
//
// Given two non-parallel segments and their associated burgers vectors, will compute
// the forces at each node.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Non_Parallel
(
         real8 *f1    ,      ///< force on node 1 <xyz>  (returned if not null)
         real8 *f2    ,      ///< force on node 2 <xyz>  (returned if not null)
         real8 *f3    ,      ///< force on node 3 <xyz>  (returned if not null)
         real8 *f4    ,      ///< force on node 4 <xyz>  (returned if not null)
   const real8 *p1    ,      ///< position of node 1 <xyz>
   const real8 *p2    ,      ///< position of node 2 <xyz>
   const real8 *p3    ,      ///< position of node 3 <xyz>
   const real8 *p4    ,      ///< position of node 4 <xyz>
   const real8 *b1    ,      ///< burgers vector (node 1->2)
   const real8 *b3    ,      ///< burgers vector (node 3->4)
   const real8  a     ,      ///< core radius
   const real8  mu    ,      ///< shear modulus
   const real8  nu    ,      ///< poisson ratio
   const real8  ecrit        ///< critical angle for parallelism test
)
{
   // compute constants of integration...

   const real8 a2     = a*a         ;
   const real8 m4p    = 0.25*mu/M_PI;
   const real8 m4pn   = m4p/(1.0-nu);
   const real8 a2m4pn = a2*m4pn     ;
   const real8 a2m8p  = a2*m4p/2.0  ;

   real8 v1[3]     = { p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2] };                   // v1 = (p2-p1)
   real8 v3[3]     = { p4[0]-p3[0], p4[1]-p3[1], p4[2]-p3[2] };                   // v3 = (p4-p3)
                                                                                  //
   real8 v1dot     = (v1[0]*v1[0]) + (v1[1]*v1[1]) + (v1[2]*v1[2]);               // v1dot = dot(v1,v1)
   real8 v3dot     = (v3[0]*v3[0]) + (v3[1]*v3[1]) + (v3[2]*v3[2]);               // v3dot = dot(v3,v3)
                                                                                  //
   real8 oneoverL  = ( (v3dot<1.0e-20) ? 0.0 : 1.0/sqrt(v3dot) );                 // oneoverL  = 1/||p4-p3||
   real8 oneoverLp = ( (v1dot<1.0e-20) ? 0.0 : 1.0/sqrt(v1dot) );                 // oneoverLp = 1/||p2-p1||
   real8 t [3]     = { (v3[0]*oneoverL ), (v3[1]*oneoverL ), (v3[2]*oneoverL ) }; // t  = (p4-p3)/||p4-p3||
   real8 tp[3]     = { (v1[0]*oneoverLp), (v1[1]*oneoverLp), (v1[2]*oneoverLp) }; // tp = (p2-p1)/||p2-p1||
                                                                                  //
   real8 c         = (t[0]*tp[0]) + (t[1]*tp[1]) + (t[2]*tp[2]);                  // c         = dot(t,tp)
   real8 c2        = (c*c);                                                       // c2        = c*c
         c2        = ( (c2<1.0) ? c2 : 1.0 );                                     // c2        = min(c*c,1.0)
   real8 onemc2    = 1.0-c2;                                                      // onemc2    = (1-c*c)
   real8 onemc2inv = ( (onemc2>0.0) ? (1.0/onemc2) : 0.0 );                       // onemc2inv = 1/(1-c*c)
                                                                                  //
   real8 txtp[3]   = { ((t[1]*tp[2]) - (t[2]*tp[1])),                             // txtp = cross(t,tp)
                       ((t[2]*tp[0]) - (t[0]*tp[2])),                             //
                       ((t[0]*tp[1]) - (t[1]*tp[0])) };                           //
                                                                                  //
   real8 mtxtp[3]  = { -txtp[0], -txtp[1], -txtp[2] };                            // mtxtp = -cross(t,tp);
                                                                                  //
   real8 R[12]     = { (p3[0]-p1[0]), (p3[1]-p1[1]), (p3[2]-p1[2]),               //
                       (p3[0]-p2[0]), (p3[1]-p2[1]), (p3[2]-p2[2]),               //
                       (p4[0]-p1[0]), (p4[1]-p1[1]), (p4[2]-p1[2]),               //
                       (p4[0]-p2[0]), (p4[1]-p2[1]), (p4[2]-p2[2]) };             //

   real8 d         = 0.5 * (  (R[0]*txtp[0]) + (R[ 1]*txtp[1]) + (R[ 2]*txtp[2])
                            + (R[9]*txtp[0]) + (R[10]*txtp[1]) + (R[11]*txtp[2]) ) * onemc2inv;

   real8 fyi[16]   = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };
   real8 fzi[16]   = { 0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0 };
   real8   y[ 2]   = { 0,0 };
   real8   z[ 2]   = { 0,0 };

   SSF_Iso_Force_Integrals(fyi,fzi,y,z,a2,d,c,c2,onemc2,onemc2inv,R,t,tp);

   if (f1 || f2) { SSF_Iso_Forces(f1,f2,fzi,tp,t ,mtxtp,b1,b3,c,-d,m4p,m4pn,a2m4pn,a2m8p,z,oneoverLp); }
   if (f3 || f4) { SSF_Iso_Forces(f3,f4,fyi,t ,tp, txtp,b3,b1,c, d,m4p,m4pn,a2m4pn,a2m8p,y,oneoverL ); }
}

// SSF_Iso_Parallel()
//
// Given a pair of parallel segments and their associated burgers vectors, will compute the
// forces on each node. Basically, this routine will fudge one of the segments to be
// non-parallel, invoke the non-parallel version of the force routines and interpolate
// for the parallel case.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso_Parallel
(
         real8  *f1    ,    ///< resulting force on node 1 <xyz>
         real8  *f2    ,    ///< resulting force on node 2 <xyz>
         real8  *f3    ,    ///< resulting force on node 3 <xyz>
         real8  *f4    ,    ///< resulting force on node 4 <xyz>
   const real8  *p1    ,    ///< position of node 1 <xyz>
   const real8  *p2    ,    ///< position of node 2 <xyz>
   const real8  *p3    ,    ///< position of node 3 <xyz>
   const real8  *p4    ,    ///< position of node 4 <xyz>
   const real8  *b1    ,    ///< burgers vector for segment p1->p2
   const real8  *b3    ,    ///< burgers vector for segment p3->p4
   const real8   a     ,    ///< core radius
   const real8   mu    ,    ///< shear modulus
   const real8   nu    ,    ///< poisson ratio
   const real8   ecrit      ///< critical angle
)
{
   // make local copies of the node positions so that we can translate them
   // to the mean center of the two segments.

   real8 lp1[3] = {  p1[0],  p1[1],  p1[2] };                 // lp1 =  p1 (local copy)
   real8 lp2[3] = {  p2[0],  p2[1],  p2[2] };                 // lp2 =  p2 (local copy)
   real8 lp3[3] = {  p3[0],  p3[1],  p3[2] };                 // lp3 =  p3 (local copy)
   real8 lp4[3] = {  p4[0],  p4[1],  p4[2] };                 // lp4 =  p4 (local copy)
   real8 lb1[3] = {  b1[0],  b1[1],  b1[2] };                 // lb1 =  b1 (local copy)
   real8 lb3[3] = {  b3[0],  b3[1],  b3[2] };                 // lb3 =  b3 (local copy)
   real8 mb1[3] = { -b1[0], -b1[1], -b1[2] };                 // mb1 = -b1 (local copy)

   // compute the mean center of the two segments.

   real8  lpc[3] = { 0.25*(lp1[0]+lp2[0]+lp3[0]+lp4[0]),      // lpc[0] = mean X
                     0.25*(lp1[1]+lp2[1]+lp3[1]+lp4[1]),      // lpc[1] = mean Y
                     0.25*(lp1[2]+lp2[2]+lp3[2]+lp4[2]) };    // lpc[2] = mean Z

   // translate the nodes to the mean center...

   lp1[0]-=lpc[0];  lp1[1]-=lpc[1];  lp1[2]-=lpc[2];
   lp2[0]-=lpc[0];  lp2[1]-=lpc[1];  lp2[2]-=lpc[2];
   lp3[0]-=lpc[0];  lp3[1]-=lpc[1];  lp3[2]-=lpc[2];
   lp4[0]-=lpc[0];  lp4[1]-=lpc[1];  lp4[2]-=lpc[2];

   // This test is used to flip the directions of the segments so that they
   // are parallel (and not antiparallel) when they are used as arguments
   // to the SSF_Iso_Special_Remote_Node_Force_Half() routines.

   int flip = ( ((   (lp4[0]-lp3[0])*(lp2[0]-lp1[0])
                   + (lp4[1]-lp3[1])*(lp2[1]-lp1[1])
                   + (lp4[2]-lp3[2])*(lp2[2]-lp1[2]) )<0.0) ? 1 : 0 );

   // Setup the arguments to the SSF_Iso_Special_Remote_Node_Force_Half() based on results of flip test.
   // This should hopefully eliminate one source of warp-diversion stalls on the GPU.

   real8 *pf1=0, *pf2=0, *px1=0, *px2=0, *pb1=0;

   if (flip) { pf1=f2; pf2=f1; px1=lp2; px2=lp1; pb1=mb1; }
   else      { pf1=f1; pf2=f2; px1=lp1; px2=lp2; pb1=lb1; }

   if (f3 || f4) { SSF_Iso_Special_Remote_Node_Force_Half( f3, f4, px1,px2,lp3,lp4, pb1,lb3, a,mu,nu,ecrit); }
   if (f1 || f2) { SSF_Iso_Special_Remote_Node_Force_Half(pf1,pf2, lp3,lp4,px1,px2, lb3,pb1, a,mu,nu,ecrit); }
}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso
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
      if (parallel) { SSF_Iso_Parallel    (f1,f2,f3,f4, x1,x2,x3,x4, b1,b3, a,mu,nu,ecrit); }
      else          { SSF_Iso_Non_Parallel(f1,f2,f3,f4, x1,x2,x3,x4, b1,b3, a,mu,nu,ecrit); }
   }
}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso
(
          real8    *fv   ,  ///< resulting forces                     (f1,f2,f3,f4,...)
    const real8    *pv   ,  ///< source positions and burgers vectors (p1,p2,p3,p4,b12,b34,...)
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

   SSF_Iso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit, lx,ly,lz, sx,sy,sz);
}

//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso
(
          SSF_FV_t *fv   ,  ///< returned forces        (f1,f2,f3,f4,...)
    const real8    *nv   ,  ///< source node positions  (n1,n2,n3,n4,n5,n6,n7,...)
    const SSF_PV_t *pv   ,  ///< source pvec structures (pv1,pv2,pv3,pv4,pv5,...)
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

   SSF_Iso(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit, lx,ly,lz, sx,sy,sz);
}

// SSF_Iso()
//
// Alternate version that doesn't have precomputed constants of integration.
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso
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
   const real8   a         ,    ///< core radius
   const real8   mu        ,    ///< shear modulus
   const real8   nu        ,    ///< poisson ratio
   const real8   ecrit          ///< critical angle for parallelism test
)
{
   // set default result (zero)...

   if (f1) { V3_ZERO(f1); }
   if (f2) { V3_ZERO(f2); }
   if (f3) { V3_ZERO(f3); }
   if (f4) { V3_ZERO(f4); }

   // cope with zero-length segments...

   if ( SSF_Seg_Zero(p1,p2) || SSF_Seg_Zero(p3,p4) ) { return; }

   // compute forces based on whether the segment pair is parallel/non-parallel...

   const int   para   = SSF_Segs_Parallel(p1,p2,p3,p4,ecrit);

   if (para) SSF_Iso_Parallel    (f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
   else      SSF_Iso_Non_Parallel(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);

#ifdef USE_RATIONAL_SEG_FORCES
   SSF_Iso_Rational_Correction_Half(f1,f2, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
   SSF_Iso_Rational_Correction_Half(f3,f4, p3,p4,p1,p2, b3,b1, a,mu,nu,ecrit);
#endif
}

// SSF_Iso() (alternate)
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso
(
         real8  &f1x,       real8  &f1y,       real8  &f1z,    ///< resulting force on node 1 <xyz>
         real8  &f2x,       real8  &f2y,       real8  &f2z,    ///< resulting force on node 2 <xyz>
         real8  &f3x,       real8  &f3y,       real8  &f3z,    ///< resulting force on node 3 <xyz>
         real8  &f4x,       real8  &f4y,       real8  &f4z,    ///< resulting force on node 4 <xyz>

   const real8   p1x, const real8   p1y, const real8   p1z,    ///< position of node 1 <xyz>
   const real8   p2x, const real8   p2y, const real8   p2z,    ///< position of node 2 <xyz>
   const real8   p3x, const real8   p3y, const real8   p3z,    ///< position of node 3 <xyz>
   const real8   p4x, const real8   p4y, const real8   p4z,    ///< position of node 4 <xyz>

   const real8   b1x, const real8   b1y, const real8   b1z,    ///< burgers vector for segment p1->p2
   const real8   b3x, const real8   b3y, const real8   b3z,    ///< burgers vector for segment p3->p4

   const real8   a     ,                                       ///< core radius
   const real8   mu    ,                                       ///< shear modulus
   const real8   nu    ,                                       ///< poisson ratio
   const real8   ecrit                                         ///< critical angle for parallelism test
)
{
         real8 f1[3] = { 0.0 };
         real8 f2[3] = { 0.0 };
         real8 f3[3] = { 0.0 };
         real8 f4[3] = { 0.0 };

   const real8 p1[3] = { p1x, p1y, p1z };
   const real8 p2[3] = { p2x, p2y, p2z };
   const real8 p3[3] = { p3x, p3y, p3z };
   const real8 p4[3] = { p4x, p4y, p4z };

   const real8 b1[3] = { b1x, b1y, b1z };
   const real8 b3[3] = { b3x, b3y, b3z };

   // set default result (zero)...

   f1x = f1y = f1z = 0.0;
   f2x = f2y = f2z = 0.0;
   f3x = f3y = f3z = 0.0;
   f4x = f4y = f4z = 0.0;

   // cope with zero-length segments...

   if ( SSF_Seg_Zero(p1,p2) || SSF_Seg_Zero(p3,p4) ) { return; }

   // compute forces based on whether the segment pair is parallel/non-parallel...

   const int   para   = SSF_Segs_Parallel(p1,p2,p3,p4,ecrit);

   if (para) SSF_Iso_Parallel    (f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
   else      SSF_Iso_Non_Parallel(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);

#ifdef USE_RATIONAL_SEG_FORCES
   SSF_Iso_Rational_Correction_Half(f1,f2, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
   SSF_Iso_Rational_Correction_Half(f3,f4, p3,p4,p1,p2, b3,b1, a,mu,nu,ecrit);
#endif

   // transfer results to args...

   f1x = f1[0];  f1y = f1[1];  f1z = f1[2];
   f2x = f2[0];  f2y = f2[1];  f2z = f2[2];
   f3x = f3[0];  f3y = f3[1];  f3z = f3[2];
   f4x = f4[0];  f4y = f4[1];  f4z = f4[2];
}

// SSF_Iso() (alternate)
//-----------------------------------------------------------------------------------------------

__cuda_hdev__
void SSF_Iso
(
         real8  *f1x,       real8  *f1y,       real8  *f1z,    ///< resulting force on node 1 <xyz>
         real8  *f2x,       real8  *f2y,       real8  *f2z,    ///< resulting force on node 2 <xyz>
         real8  *f3x,       real8  *f3y,       real8  *f3z,    ///< resulting force on node 3 <xyz>
         real8  *f4x,       real8  *f4y,       real8  *f4z,    ///< resulting force on node 4 <xyz>

   const real8   p1x, const real8   p1y, const real8   p1z,    ///< position of node 1 <xyz>
   const real8   p2x, const real8   p2y, const real8   p2z,    ///< position of node 2 <xyz>
   const real8   p3x, const real8   p3y, const real8   p3z,    ///< position of node 3 <xyz>
   const real8   p4x, const real8   p4y, const real8   p4z,    ///< position of node 4 <xyz>

   const real8   b1x, const real8   b1y, const real8   b1z,    ///< burgers vector for segment p1->p2
   const real8   b3x, const real8   b3y, const real8   b3z,    ///< burgers vector for segment p3->p4

   const real8   a     ,                                       ///< core radius
   const real8   mu    ,                                       ///< shear modulus
   const real8   nu    ,                                       ///< poisson ratio
   const real8   ecrit                                         ///< critical angle for parallelism test
)
{
         real8 f1[3] = { 0.0 };
         real8 f2[3] = { 0.0 };
         real8 f3[3] = { 0.0 };
         real8 f4[3] = { 0.0 };

   const real8 p1[3] = { p1x, p1y, p1z };
   const real8 p2[3] = { p2x, p2y, p2z };
   const real8 p3[3] = { p3x, p3y, p3z };
   const real8 p4[3] = { p4x, p4y, p4z };

   const real8 b1[3] = { b1x, b1y, b1z };
   const real8 b3[3] = { b3x, b3y, b3z };

   // set default result (zero)...

   *f1x = *f1y = *f1z = 0.0;
   *f2x = *f2y = *f2z = 0.0;
   *f3x = *f3y = *f3z = 0.0;
   *f4x = *f4y = *f4z = 0.0;

   // cope with zero-length segments...

   if ( SSF_Seg_Zero(p1,p2) || SSF_Seg_Zero(p3,p4) ) { return; }

   // compute forces based on whether the segment pair is parallel/non-parallel...

   const int   para   = SSF_Segs_Parallel(p1,p2,p3,p4,ecrit);

   if (para) SSF_Iso_Parallel    (f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
   else      SSF_Iso_Non_Parallel(f1,f2,f3,f4, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);

   // apply rational correction (if active)...

#ifdef USE_RATIONAL_SEG_FORCES
   SSF_Iso_Rational_Correction_Half(f1,f2, p1,p2,p3,p4, b1,b3, a,mu,nu,ecrit);
   SSF_Iso_Rational_Correction_Half(f3,f4, p3,p4,p1,p2, b3,b1, a,mu,nu,ecrit);
#endif

   // transfer results to args...

   *f1x = f1[0];  *f1y = f1[1];  *f1z = f1[2];
   *f2x = f2[0];  *f2y = f2[1];  *f2z = f2[2];
   *f3x = f3[0];  *f3y = f3[1];  *f3z = f3[2];
   *f4x = f4[0];  *f4y = f4[1];  *f4z = f4[2];
}

