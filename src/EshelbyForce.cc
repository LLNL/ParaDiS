/*****************************************************************************
 *
 *      Module:      EshelbyForce.c
 *
 *      Description: Contains primary functions for calculating the
 *                   force on a displocation segment from an Eshelby
 *                   inclusion. This function assumes that the strain is such that
 *                   it is zero except for strain[0][0], strain[1][1], strain[2][2].
 *                   Also the these three components are equal: Volumetric strain only.
 *
 *      Includes public functions:
 *
 *          EshelbyForce()
 *
 *      Includes private functions:
 *          EshelbyIntegrals
 *          EshelbySpecialIntegrals
 *
 *****************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "V3.h"


#ifdef ESHELBY

#ifdef ESHELBYFORCE
static void EshelbySpecialIntegrals(real8 s0, real8 s1, real8 wIntSp[2])
{
   // [s0,s1] is inside the particle
        wIntSp[0] = s1 - s0;
        wIntSp[1] = 0.5 * (s1*s1 - s0*s0);

        return;
}

static void EshelbyIntegrals(real8 s0, real8 s1, real8 d2, real8 wInt[5])
{
   // [s0,s1] is outside the particle
        int   i;
        real8 temp, eps;
        real8 s[2], s2[2], d2v[2], d2Inv[2], mf[2];
        real8 Ra[2], RaInv[2], Ra3Inv[2], sRa3Inv[2];
        real8 s_03[2], s_13[2], s_05[2], s_15[2], s_25[2];


        eps = 1.0e-08;

        s[0] = s0;
        s[1] = s1;

        s2[0] = s[0] * s[0];
        s2[1] = s[1] * s[1];

        d2v[0] = d2;
        d2v[1] = d2;

        temp = 1.0 / d2;

        d2Inv[0] = temp;
        d2Inv[1] = temp;

        for (i = 0; i < 2; i++) {

            Ra[i]      = sqrt(d2v[i] + s[i]*s[i]);
            RaInv[i]   = 1.0 / Ra[i];
            Ra3Inv[i]  = RaInv[i]*RaInv[i]*RaInv[i];
            sRa3Inv[i] = s[i]*Ra3Inv[i];

/*
 *          If the tangent line of the segment hits the center of the
 *          inclusion, we have to do some special handling, otherwise
 *          deal with it as usual via the else clause.
 */
            if (d2 < ((s2[0] + s2[1]) * eps)) {
                s_03[i] = -0.5 * sRa3Inv[i];
                s_13[i] = -RaInv[i];
                s_05[i] = 0.0;
                s_15[i] = 0.0;
                s_25[i] = 0.0;
            } else {
                s_03[i] = s[i]*RaInv[i]*d2Inv[i];
                s_13[i] = -RaInv[i];
                s_05[i] = (2.0*s_03[i] + sRa3Inv[i]) * d2Inv[i];
                s_15[i] = -Ra3Inv[i];
                s_25[i] = s_03[i] - sRa3Inv[i];
            }
        }

        mf[0] = -1.0;
        mf[1] =  1.0;

        wInt[0] = s_03[0]*mf[0] + s_03[1]*mf[1];
        wInt[1] = s_13[0]*mf[0] + s_13[1]*mf[1];
        wInt[2] = s_05[0]*mf[0] + s_05[1]*mf[1];
        wInt[3] = s_15[0]*mf[0] + s_15[1]*mf[1];
        wInt[4] = s_25[0]*mf[0] + s_25[1]*mf[1];

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:        Tom Arsenlis
 *
 *      Function:      EshelbyForce
 *
 *      Description:   Caluclate the force on the endpoints of a dislocation
 *                     segment going from <pos1> to <pos2> resulting from
 *                     an Eshelby inclusion.
 *
 *      Arguments:
 *          EPos     Center position of the Eshelby inclusion
 *          radius   Radius of the Eshelby inclusion
 *          EStrain  Strain field of the Eshelby inclusion
 *          cellID   ID of the cell encompassing the inclusion
 *          pos1     1st endpoint position of dislocation segment
 *          pos2     2nd endpoint position of dislocation segment
 *          burg     burgers vectore of the dislocation segment
 *                   going from pos1 to pos2
 *          MU       Shear modulus
 *          NU       Poisson ratio
 *          f1,f2    Arrays in which to return the force components at
 *                   the dislocation endpoints (pos1 and pos2 respectively)
 *
 *      Last Modified: 01/28/2008 - original versionQ
 *                     07/01/2016 - Added intersection with an ellipse (S. Aubry).
 *
 *-------------------------------------------------------------------------*/
void EshelbyForce(Home_t *home, real8 EPos[3], real8 radius[3], 
                  real8 EStrain[6], real8 newPos1[3], real8 newPos2[3],
                  real8 burg[3], real8 f1[3], real8 f2[3])
{
        int   i;
        real8 d2, oneOverL, Ev, test;
        real8 factor, factorSP;
        real8 F_03, F_05, F_15, F_SP;
        real8 bdnd, bdt;
        real8 s[2], sinter[2];
        real8 MaxRadius, diff[3], t[3], nd[3], bct[3], ndct[3];
        real8 R[2][3];
        real8 I_03[3], I_05[3], I_15[3];
        real8 wInt[5], wIntTmp[5], wIntSp[2];

        Param_t *param;
        param = home->param;

        real8 MU, NU, MU2, NU2, factorSP2;
        MU = param->shearModulus;
        NU = param->pois;
        MU2 = param->shearModulus2;
        NU2 = param->pois2;

        f1[0] = 0.0;  f1[1] = 0.0;  f1[2] = 0.0;
        f2[0] = 0.0;  f2[1] = 0.0;  f2[2] = 0.0;

/*
 *      If the strain field for the particle is zero, then 
 *      there's nothing more to be done in this routine.
 */
        if ( ((EStrain[0] == 0.0) &&
              (EStrain[1] == 0.0) &&
              (EStrain[2] == 0.0)) ) 
        {
           return;
        }
        else
        {
           // If EStrain is non zero, then we are probably dealing with spherical particles.
           // The following is valid only for spherical particles.
           MaxRadius = V3_MAX(radius[0],radius[1],radius[2]);

           for (i = 0; i < 3; i++) {
              diff[i] = newPos2[i] - newPos1[i];
              R[0][i] = newPos1[i] - EPos[i];
              R[1][i] = newPos2[i] - EPos[i];
           }
           
           oneOverL = 1.0 / Normal(diff);
           
           for (i = 0; i < 3; i++) {
              t[i] = diff[i] * oneOverL;
           }
           
           s[0] = R[0][0]*t[0] + R[0][1]*t[1] + R[0][2]*t[2];
           s[1] = R[1][0]*t[0] + R[1][1]*t[1] + R[1][2]*t[2];
           
           for (i = 0, d2 = 0.0; i < 3; i++) {
              nd[i] = R[0][i] - s[0]*t[i];
              d2 += nd[i] * nd[i];
           }
           
           test = MaxRadius * MaxRadius - d2;
           

           // sinter are the two intersection points between the sphere and the dislocation segment.
           if (test > 0.0) {
              sinter[0] = -sqrt(test);
              sinter[1] = -sinter[0];
           } else {
              sinter[0] = 0.0;
              sinter[1] = 0.0;
           }
 
           if ((test < 0.0) | (sinter[0] > s[1]) | (sinter[1] < s[0])) 
              {
                 EshelbyIntegrals(s[0], s[1], d2, wInt);
                 wIntSp[0] = 0.0;
                 wIntSp[1] = 0.0;
              } 
           else if ((sinter[0] < s[0]) && (sinter[1] < s[1])) 
           {
              EshelbyIntegrals(sinter[1], s[1], d2, wInt);
              EshelbySpecialIntegrals(s[0], sinter[1], wIntSp);
           } 
           else if ((sinter[0] > s[0]) && (sinter[1] > s[1])) 
           {
              EshelbyIntegrals(s[0], sinter[0], d2, wInt);
              EshelbySpecialIntegrals(sinter[0], s[1], wIntSp);
           } 
           else if ((sinter[0] < s[0]) && (sinter[1] > s[1])) 
           {
              wInt[0] = 0.0;
              wInt[1] = 0.0;
              wInt[2] = 0.0;
              wInt[3] = 0.0;
              wInt[4] = 0.0;
              EshelbySpecialIntegrals(s[0], s[1], wIntSp);
           } 
           else 
           { /* ((sinter[0] > s[0]) && (sinter[1] < s[1])) */
              EshelbyIntegrals(s[0], sinter[0], d2, wInt);
              EshelbyIntegrals(sinter[1], s[1], d2, wIntTmp);
              wInt[0] += wIntTmp[0];
              wInt[1] += wIntTmp[1];
              wInt[2] += wIntTmp[2];
              wInt[3] += wIntTmp[3];
              wInt[4] += wIntTmp[4];
              EshelbySpecialIntegrals(sinter[0], sinter[1], wIntSp);
           }
           
           bdnd = DotProduct(burg, nd);
           bdt = DotProduct(burg, t);
           
           cross(burg, t, bct);
           cross(nd, t, ndct);
           
           for (i = 0; i < 3; i++) {
              I_03[i] = bct[i];
              I_05[i] = bdnd * ndct[i];
              I_15[i] = bdt * ndct[i];
           }
           
           F_03 = wInt[1] - s[0]*wInt[0];
           F_05 = wInt[3] - s[0]*wInt[2];
           F_15 = wInt[4] - s[0]*wInt[3];
           F_SP = wIntSp[1] - s[0]*wIntSp[0];

           Ev = EStrain[0] + EStrain[1] + EStrain[2];
           factorSP = 2.0 * MU * (1.0 + NU) * Ev / (1.0 - (2.0 * NU)) * oneOverL;
           factor = -factorSP * MaxRadius * MaxRadius * MaxRadius * 0.5;
           
           factorSP2 = 2.0 * MU2 * (1.0 + NU2) * Ev / (1.0 - (2.0 * NU2)) * oneOverL;


           for (i = 0; i < 3; i++) {
              f2[i] += factor * ((I_03[i]*F_03)+(I_05[i]*F_05)+(I_15[i]*F_15)) +
                 (factorSP2 * F_SP * I_03[i]);
           }
           
           F_03 = s[1]*wInt[0] - wInt[1];
           F_05 = s[1]*wInt[2] - wInt[3];
           F_15 = s[1]*wInt[3] - wInt[4];
           F_SP = s[1]*wIntSp[0] - wIntSp[1];
           
           for (i = 0; i < 3; i++) {
              f1[i] += factor * ((I_03[i]*F_03)+(I_05[i]*F_05)+(I_15[i]*F_15)) +
                 (factorSP2 * F_SP * I_03[i]);
           }
        }

        return;
}
#endif // ESHELBYFORCE
#endif // ESHELBY
