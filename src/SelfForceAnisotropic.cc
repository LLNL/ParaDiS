/***************************************************************************
 *
 *      Module:      SelfForceAnisotropic.c
 *
 *      Description: This module contains the functions needed for
 *                   calculating anisotropic self-force for a given
 *                   dislocation segment. The self force is composed
 *                   of the elastic self-force and the core force.
 *                   Both are computed using an expansion in spherical
 *                   harmonics.
 *
 *                   This routine uses only the single integral
 *
 *                       Jz_ijk = \int (R . c12)^i (R . e3)^j dz
 *                                    ----------------------
 *                                            Ra^k
 *
 *      Includes public functions:
 *          SelfForceAnisotropic()
 *          CoreForceAnisotropic()
 *          CoreForceandDerivativeAnisotropic()
 *
 *
 ***************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

#ifdef ANISOTROPIC

//---------------------------------------------------------------------------
// The IBM compiler continues to manifest errors when casting variable 
// length and multidimensional matrices.  The ptrcpy() macro copies pointers 
// via memcpy() to avoid the casting problems.
//---------------------------------------------------------------------------

#define ptrcpy(a,b)  memcpy(&(a),&(b),sizeof(void *))

/*---------------------------------------------------------------
 *
 *      Function:        CoreForceAnisotropic
 *      Description:     Calculate the core force for anisotropic
 *                       elasticity. 
 *
 *      Parameters:
 *          IN:  x1,y1,z1, x2,y,2,z2   Dislocation segment
 *          IN:  b                     Burgers vector
 *          IN:  qMax                  Number of expansion coeffs. for 
 *                                     spherical harmonics calculations of
 *                                     the core.
 *          IN:  Ecore                 Core energy
 *          OUT: fCore                 Force coming from core effects at node 2
 *
 *--------------------------------------------------------------------------*/
void CoreForceAnisotropic(Home_t *home, 
                          real8 x1, real8 y1, real8 z1, 
                          real8 x2, real8 y2, real8 z2,
                          real8 b[3], real8 Ecore, int qMax, 
                          real8 fCore[3])
{

    int i, p;
    int Lmax = 2*qMax;
    int qMaxProd = (qMax+2)*(qMax+1)-1;
    int ind, twoind,M,P,local;
    Param_t *param = home->param;

    real8 L;
    real8 t[3];
    real8 bv[6];
    real8 te3;
    real8 amzt[3];
    real8 qv[Lmax+1][2];
    real8 reA[6], reC[6];
    real8 *e3;
    real8 reCbv=0.0, reAbv=0.0;

    //
    // "complex" arrays below are represented as pairs of 
    // arrays for the real and imaginary components.
    //
    real8 (*FqReal)[qMaxProd]; /* dimensions are [6][qMaxProd] */
    real8 (*FqImag)[qMaxProd]; /* dimensions are [6][qMaxProd] */

    real8 te12Real;
    real8 te12Imag;

    real8 cmytReal[3];
    real8 cmytImag[3];

    real8 pvReal[Lmax+1][2];
    real8 pvImag[Lmax+1][2];

    real8 BReal[6];
    real8 BImag[6];

    real8 AtempReal, BtempReal, CtempReal;
    real8 AtempImag, BtempImag, CtempImag;

    real8 SReal, *e12Real;
    real8 SImag, *e12Imag;

    real8 BbvReal = 0.0;
    real8 BbvImag = 0.0;

/*
 *  For HCP Be, we have MD data that gives Ecore as a function of Burgers 
 *  vector type. Use it here.
 */
    if (home->param->materialType == MAT_TYPE_HCP) 
    {
        real8 bCryst[3];

/*
 *      If needed, rotate the burgers vector from the laboratory frame
 *      to the crystal frame.
 */
        if (param->useLabFrame) {
            real8 bLab[3];
            bLab[X] = b[0];
            bLab[Y] = b[1];
            bLab[Z] = b[2];
            Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[X] = b[0];
            bCryst[Y] = b[1];
            bCryst[Z] = b[2];
        }

        Ecore = HCPEcoreFactor(home, bCryst);

    }  /* end if (param->materialType == MAT_TYPE_HCP) */

/*
 *  Need some stuff from 'home'
 */
    ptrcpy(FqReal,home->anisoVars.coreFqReal);
    ptrcpy(FqImag,home->anisoVars.coreFqImag);

    e12Real = &home->anisoVars.anisoRotMatrix[0][0];
    e12Imag = &home->anisoVars.anisoRotMatrix[1][0];

    e3  = &home->anisoVars.anisoRotMatrix[2][0];

    GetUnitVector(1, x1, y1, z1, x2, y2, z2, &(t[0]), &(t[1]), &(t[2]), &L);

    te12Real = DotProduct(t, e12Real);
    te12Imag = DotProduct(t, e12Imag);

    te3  = DotProduct(t,  e3);

    for (i=0;i<3;i++) {

        cmytReal[i] = e12Real[i] - te12Real * t[i];
        cmytImag[i] = e12Imag[i] - te12Imag * t[i];

        amzt[i] = e3[i]  - te3*t[i];
    }

    pvReal[0][0]=1.0;
    pvImag[0][0]=0.0;

    pvReal[0][1]=0.0;
    pvImag[0][1]=0.0;

    qv[0][0]=1.0;
    qv[0][1]=0.0;

    for (i=1;i<=Lmax;i++) {
        pvReal[i][0] = (pvReal[i-1][0] * te12Real) -
                       (pvImag[i-1][0] * te12Imag);
        pvImag[i][0] = (pvReal[i-1][0] * te12Imag) +
                       (pvImag[i-1][0] * te12Real);

        pvReal[i][1] = (pvReal[i-1][0] * i);
        pvImag[i][1] = (pvImag[i-1][0] * i);

        qv[i][0]=qv[i-1][0]*te3;
        qv[i][1]=qv[i-1][0]*i;
    }
 
    for (p=0;p<6;p++) {
        reA[p]   = 0.0;
        reC[p]   = 0.0;
        BReal[p] = 0.0;
        BImag[p] = 0.0;
    }
 
    local = 0;

    for (ind=0;ind<qMax+1;ind++) {
        twoind = ind+ind;
        P = twoind;
        for (M=0;M<twoind+1;M++) { 
            AtempReal = qv[P][1] * pvReal[M][0];
            AtempImag = qv[P][1] * pvImag[M][0];

            BtempReal = qv[P][0] * pvReal[M][1];
            BtempImag = qv[P][0] * pvImag[M][1];

            CtempReal = qv[P][0] * pvReal[M][0];
            CtempImag = qv[P][0] * pvImag[M][0];

            for (p=0;p<6;p++)
            {
                SReal = FqReal[p][local];
                SImag = FqImag[p][local];
    
                reA[p] += (SReal*AtempReal) - (SImag*AtempImag);

                BReal[p] += (SReal * BtempReal) - 
                            (SImag * BtempImag);
                BImag[p] += (SReal * BtempImag) + 
                            (SImag * BtempReal);

                reC[p] += (SReal*CtempReal) - (SImag*CtempImag);
            }

            P--;
            local++;
        }
    }

    bv[0] = b[0]*b[0];
    bv[1] = b[1]*b[1];
    bv[2] = b[2]*b[2];
    bv[3] = 2*b[1]*b[2];
    bv[4] = 2*b[2]*b[0];
    bv[5] = 2*b[0]*b[1];

    for (p=0;p<6;p++) {

        reAbv += reA[p]*bv[p];

        BbvReal   +=   BReal[p] * bv[p];
        BbvImag   +=   BImag[p] * bv[p];

        reCbv += reC[p]*bv[p];
    }

    for (i=0;i<3;i++) {
        fCore[i] = -Ecore * (t[i]*reCbv +  amzt[i]*reAbv +
                             (cmytReal[i]*BbvReal) - (cmytImag[i]*BbvImag));  
    }

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:        CoreForceandDerivativeAnisotropic
 *      Description:     Calculate the core force and its derivative for
 *                       anisotropic elasticity.  The force derivative
 *                       is only computed if requested (see below).
 *
 *      Parameters:
 *          IN: x1,y1,z1 and x2,y2,z2  Positions of the two endpoints of the
 *                                     segment for which forces are being
 *                                     calculated
 *          IN:  b                     Burgers vector
 *          OUT: fCore                 Force coming from core effects at node
 *                                     x2, y2, z2
 *          OUT: dfCoredx              Derivative of Force coming from core
 *                                     effects at node x2, y2, z2
 *
 *--------------------------------------------------------------------------*/
void CoreForceandDerivativeAnisotropic(Home_t *home, 
                real8 x1, real8 y1, real8 z1, 
                real8 x2, real8 y2, real8 z2,
                real8 b[3], real8 Ecore, int qMax,
                real8 fCore[3], real8 dfCoredx[3][3])
{
    int   i, j, p;
    int   Lmax = 2*qMax;
    int   qMaxProd = (qMax+2)*(qMax+1)-1;
    int   ind, twoind, M, P, local;
    Param_t *param = home->param;

    real8 L, Linv;
    real8 t[3], Fac, temp;
    real8 bv[6];
    real8 te3;
    real8 amzt[3];
    real8 Dmat[3][3], Lmat[3][3];
    real8 qv[Lmax+1][3];
    real8 reA[6], reC[6], reD[6];
    real8 *e3;
    real8 reCbv=0.0, reAbv=0.0, reDbv=0.0;

    //
    // "complex" arrays below are represented as separate
    // arrays for the real and imaginary components.
    //
    real8 (*FqReal)[qMaxProd]; /* dimensions are [6][qMaxProd] */
    real8 (*FqImag)[qMaxProd]; /* dimensions are [6][qMaxProd] */

    real8 te12Real;
    real8 te12Imag;

    real8 cmytReal[3];
    real8 cmytImag[3];

    real8 EmatReal[3][3], GmatReal[3][3];
    real8 EmatImag[3][3], GmatImag[3][3];

    real8 pvReal[Lmax+1][3];
    real8 pvImag[Lmax+1][3];

    real8 BReal[6], EReal[6], GReal[6];
    real8 BImag[6], EImag[6], GImag[6];

    real8 AtempReal, BtempReal, CtempReal, DtempReal, EtempReal, GtempReal;
    real8 AtempImag, BtempImag, CtempImag, DtempImag, EtempImag, GtempImag;

    real8 SReal, *e12Real;
    real8 SImag, *e12Imag;

    real8 BbvReal = 0.0;
    real8 BbvImag = 0.0;

    real8 EbvReal = 0.0;
    real8 EbvImag = 0.0;

    real8 GbvReal = 0.0;
    real8 GbvImag = 0.0;

/*
 *  For HCP Be, we have MD data that gives Ecore as a function of Burgers 
 *  vector type. Use it here.
 */
    if (home->param->materialType == MAT_TYPE_HCP) 
    {
        real8 bCryst[3];

/*
 *      If needed, rotate the burgers vector from the laboratory frame
 *      to the crystal frame.
 */
        if (param->useLabFrame) {
            real8 bLab[3];
            bLab[X] = b[0];
            bLab[Y] = b[1];
            bLab[Z] = b[2];
            Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[X] = b[0];
            bCryst[Y] = b[1];
            bCryst[Z] = b[2];
        }

        Ecore = HCPEcoreFactor(home, bCryst);

    }  /* end if (param->materialType == MAT_TYPE_HCP) */

/*
 *  Need some stuff from 'home'
 */
    ptrcpy(FqReal,home->anisoVars.coreFqReal);
    ptrcpy(FqImag,home->anisoVars.coreFqImag);

    e12Real = &home->anisoVars.anisoRotMatrix[0][0];
    e12Imag = &home->anisoVars.anisoRotMatrix[1][0];
    e3      = &home->anisoVars.anisoRotMatrix[2][0];

    GetUnitVector(1, x1, y1, z1, x2, y2, z2, &(t[0]), &(t[1]), &(t[2]), &L);

/*
 *  In case of a zero length vector...
 */
    if (L == 0.0) {
        Linv = 0.0;
    } else {
        Linv= 1./L;
    }

    te12Real = DotProduct(t,  e12Real);
    te12Imag = DotProduct(t,  e12Imag);
    te3      = DotProduct(t,  e3);

    for (i=0;i<3;i++) {
        cmytReal[i] = e12Real[i] - te12Real * t[i];
        cmytImag[i] = e12Imag[i] - te12Imag * t[i];
        amzt[i] = e3[i]  - te3*t[i];
    }

    for (i=0;i<3;i++) {
        for (j=0;j<3;j++) {
            Dmat[i][j] = amzt[i]*amzt[j];
            Lmat[i][j] =-t[i]*t[j];

            EmatReal[i][j] = (cmytReal[i] * amzt[j]) +
                             (amzt[i] * cmytReal[j]);

            EmatImag[i][j] = (cmytImag[i] * amzt[j]) +
                             (amzt[i] * cmytImag[j]);

            GmatReal[i][j] = (cmytReal[i] * cmytReal[j]) -
                             (cmytImag[i] * cmytImag[j]);

            GmatImag[i][j] = (cmytReal[i] * cmytImag[j]) +
                             (cmytImag[i] * cmytReal[j]);
        }
    }

    Lmat[0][0] += 1.0;
    Lmat[1][1] += 1.0;
    Lmat[2][2] += 1.0;

    pvReal[0][0]=1.0;
    pvImag[0][0]=0.0;

    pvReal[0][1]=0.0;
    pvImag[0][1]=0.0;

    pvReal[0][2]=0.0;
    pvImag[0][2]=0.0;

    qv[0][0]=1.0;
    qv[0][1]=0.0;
    qv[0][2]=0.0;

    for (i=1;i<=Lmax;i++) {

        pvReal[i][0] = (pvReal[i-1][0] * te12Real) - 
                       (pvImag[i-1][0] * te12Imag);
        pvImag[i][0] = (pvReal[i-1][0] * te12Imag) + 
                       (pvImag[i-1][0] * te12Real);

        pvReal[i][1] = pvReal[i-1][0] * i;
        pvImag[i][1] = pvImag[i-1][0] * i;

        pvReal[i][2] = pvReal[i-1][1] * i;
        pvImag[i][2] = pvImag[i-1][1] * i;
     
        qv[i][0]=qv[i-1][0]*te3;
        qv[i][1]=qv[i-1][0]*i;
        qv[i][2]=qv[i-1][1]*i;
    }
 
    for (p=0;p<6;p++) {
        reA[p] = 0.0;
        reC[p] = 0.0;
        reD[p] = 0.0;

        BReal[p] = 0.0;
        BImag[p] = 0.0;

        EReal[p] = 0.0;
        EImag[p] = 0.0;

        GReal[p] = 0.0;
        GImag[p] = 0.0;
    }
 
    local = 0;

    for (ind=0;ind<qMax+1;ind++) {

        twoind = ind+ind;
        P = twoind;

        for (M=0;M<twoind+1;M++) { 
            AtempReal = pvReal[M][0] * qv[P][1];
            AtempImag = pvImag[M][0] * qv[P][1];

            BtempReal = pvReal[M][1] * qv[P][0];
            BtempImag = pvImag[M][1] * qv[P][0];

            CtempReal = pvReal[M][0] * qv[P][0];
            CtempImag = pvImag[M][0] * qv[P][0];
        
            DtempReal = pvReal[M][0] * qv[P][2]; 
            DtempImag = pvImag[M][0] * qv[P][2]; 

            EtempReal = pvReal[M][1] * qv[P][1]; 
            EtempImag = pvImag[M][1] * qv[P][1]; 

            GtempReal = pvReal[M][2] * qv[P][0];
            GtempImag = pvImag[M][2] * qv[P][0];

            for (p=0;p<6;p++) {
                SReal = FqReal[p][local];
                SImag = FqImag[p][local];

                reA[p] += SReal*AtempReal - SImag*AtempImag;
                reC[p] += SReal*CtempReal - SImag*CtempImag;
                reD[p] += SReal*DtempReal - SImag*DtempImag;

                BReal[p] += (SReal * BtempReal) - (SImag * BtempImag);
                BImag[p] += (SReal * BtempImag) + (SImag * BtempReal);

                EReal[p] += (SReal * EtempReal) - (SImag * EtempImag);
                EImag[p] += (SReal * EtempImag) + (SImag * EtempReal);

                GReal[p] += (SReal * GtempReal) - (SImag * GtempImag);
                GImag[p] += (SReal * GtempImag) + (SImag * GtempReal);
            }

            P--;
            local++;
        }
    }

    bv[0] = b[0]*b[0];
    bv[1] = b[1]*b[1];
    bv[2] = b[2]*b[2];
    bv[3] = 2*b[1]*b[2];
    bv[4] = 2*b[2]*b[0];
    bv[5] = 2*b[0]*b[1];

    for (p=0;p<6;p++) {

        reAbv += reA[p]*bv[p];
        reCbv += reC[p]*bv[p];
        reDbv += reD[p]*bv[p];

        BbvReal += BReal[p] * bv[p];
        BbvImag += BImag[p] * bv[p];

        EbvReal += EReal[p] * bv[p];
        EbvImag += EImag[p] * bv[p];

        GbvReal += GReal[p] * bv[p];
        GbvImag += GImag[p] * bv[p];
    }

/*
 *  This function computes the force and derivatives at node 2
 */
    Fac = -Ecore*Linv;

    temp = Fac*reCbv -
           te3*reAbv -
           ((te12Real * BbvReal) - (te12Imag * BbvImag));

    reDbv *= Fac;

    EbvReal *= Fac;
    EbvImag *= Fac;

    GbvReal *= Fac;
    GbvImag *= Fac;

    for (i=0;i<3;i++) {
        fCore[i] = -Ecore *
                   (t[i]*reCbv + amzt[i]*reAbv +
                    ((cmytReal[i]*BbvReal) - (cmytImag[i]*BbvImag)));
      
        for (j=0;j<3;j++) {
            real8 tmpReal;

            tmpReal = ((EmatReal[i][j] * EbvReal) -
                       (EmatImag[i][j] * EbvImag)) +
                      ((GmatReal[i][j] * GbvReal) -
                       (GmatImag[i][j] * GbvImag));

            dfCoredx[i][j] = Dmat[i][j]*reDbv + tmpReal + Lmat[i][j]*temp;
        }
    }
  
//
// STORE dFdx in a 6 matrix instead of 3x3??
//
    return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SelfForceAnisotropic()
 *      Description:  Calculates anisotropic self-force for the specified
 *                    dislocation segment
 *
 *      Parameters:
 *          IN:   b{xyz}        Burgers vector of segment 2 from p3 to p4
 *          IN:   {xyz}1       Position of first segment's first endpoint
 *          IN:   {xyz}2       Position of first segment's second endpoint
 *          IN:   a             Core radius
 *          IN:   qMax          Base factor for determining the number of terms
 *                              in the expansion for spherical harmonics.  The
 *                              total number of terms is 2 * qMax + 1.
 *          OUT:  f1{xyz}       Force on segment at position {xyz}1
 *          OUT:  f2{xyz}       Force on segment at position {xyz}2
 *
 *-------------------------------------------------------------------------*/
void SelfForceAnisotropic(Home_t *home, int coreOnly,
                          real8 bx, real8 by, real8 bz,
                          real8 x1, real8 y1, real8 z1,
                          real8 x2, real8 y2, real8 z2,
                          real8 a, real8 Ecore,
                          real8 f1[3], real8 f2[3])
{
        int           ii, jj, kk, iii;
        int           qMax = home->param->anisoHarmonicsNumTermsBase;
        int           qMaxp1 = qMax+1;
        int           qMaxp2 = qMax+2;
        int           twoqMaxp3 = 2*qMax+3;
        int           i, j, k;
        int           ind, indb;
        int           local, local2, local3, local4;
        int           offsetG;
        short         *RemArray, *ImmArray;
        short         *RelArray, *ImlArray;
        short         *Re18Array, *Im18Array;
        short         *ReiiArray, *ImiiArray;
        real8         *FqReal, *FqImag;
        real8         a2, coeff12;
        real8         L, L2, invL;
        real8         Ra2inv, Rasum;
        real8         delta;
        real8         bt11, bt12, bt13;
        real8         bt21, bt22, bt23;
        real8         bt31, bt32, bt33;
        real8         Ap12, Ap23, Ap31, Ap41, Ap43, Ap51, Ap52, Ap62, Ap63;
        real8         B12, B14, B16, B23, B24, B25, B31, B35, B36;
        real8         *rotMatrix3;
        real8         f2Core[3];
        real8         t[3];
        real8         Ra[4];
        real8         Ap[6][3];
        real8         B[3][6];
        real8         G[18];
        real8         ApG[6];
        real8         ApGaux[3][3];
        real8         ypzin[twoqMaxp3][2];
        real8         deltaind[twoqMaxp3];
        real8         factor[qMaxp2];
        real8         factorinv[qMaxp2];
        real8         Rainv[qMaxp1];
        real8         (*ecMatrix)[6];
        real8         reHy, imHy;
        complex8 *rotMatrix12;
        complex8 alpha, zJyzsum;
        complex8 alphaind[twoqMaxp3];
        complex8 Ht[2], Hyt[2];
        complex8 Jyzt00[2], Jyzt01[2], Jyzt02[2], Jyzt03[2];
        complex8 Jyzt10[2], Jyzt11[2], Jyzt12[2], Jyzt13[2];
	real8 b[3];


/* 
 *      Compute core effects.
 *      Note that Ecore might change in the definition of SelfForce/Core Force.
 *      Explicitly pass Ecore.
 */
	b[0]=bx; b[1]=by; b[2]=bz;

	/* To debug and compare to isotropic elasticity set Ecore = 1 */

        CoreForceAnisotropic(home, x1, y1, z1, x2, y2, z2, b, Ecore,
			     qMax,f2Core);


        if (coreOnly) {
/* 
 *          We're doing a core-only calculation, so just use the
 *          core effects calculated above and we're done.
 */
            f2[0] = f2Core[0];
            f2[1] = f2Core[1];
            f2[2] = f2Core[2];

            f1[0] = -f2[0];
            f1[1] = -f2[1];
            f1[2] = -f2[2];
            return;
        }

/*
 *      Define line direction
 */
        t[0] = x2 - x1;
        t[1] = y2 - y1;
        t[2] = z2 - z1;

        L2 = DotProduct(t, t);
        L = sqrt(L2);
        invL = 1.0 / L;

        coeff12 = invL * INV_4_PI_SQUARED;

        t[0] = t[0] * invL;
        t[1] = t[1] * invL;
        t[2] = t[2] * invL;

/*
 *      Need some of the stuff in 'home'
 */
        FqReal = home->anisoVars.FqReal;
        FqImag = home->anisoVars.FqImag;

        RemArray = home->anisoVars.sphericalHarmonicsRemKept;
        ImmArray = home->anisoVars.sphericalHarmonicsImmKept;
        RelArray = home->anisoVars.sphericalHarmonicsRelKept;
        ImlArray = home->anisoVars.sphericalHarmonicsImlKept;
        Re18Array= home->anisoVars.sphericalHarmonicsRe18Kept;
        Im18Array= home->anisoVars.sphericalHarmonicsIm18Kept;
        ReiiArray= home->anisoVars.sphericalHarmonicsReiKept;
        ImiiArray= home->anisoVars.sphericalHarmonicsImiKept;

        rotMatrix12 = home->anisoVars.anisoRotMatrix12;
        rotMatrix3  = &home->anisoVars.anisoRotMatrix[2][0];

        ecMatrix = home->param->elasticConstantMatrix;

/*
 *      Initialization
 */
        for (i = 0; i < 18; i++) {
            G[i] = 0.0;
        }

/*
 *      Compute integrals via recurrence
 */
        alpha = DotProduct(t,  rotMatrix12);
        delta = DotProduct(t,  rotMatrix3);

        ypzin[0][0] = 1.0;
        ypzin[0][1] = 1.0;

        ypzin[1][0] = -L;
        ypzin[1][1] =  L;

        alphaind[0] = 1.0;
        alphaind[1] = alpha;

        deltaind[0] = 1.0;
        deltaind[1] = delta;

        for (i = 2; i < twoqMaxp3; i++) {
            alphaind[i] = alphaind[i-1] * alpha;
            deltaind[i] = deltaind[i-1] * delta;
            ypzin[i][0] = ypzin[i-1][0] * ypzin[1][0];
            ypzin[i][1] = ypzin[i-1][1] * ypzin[1][1];
        }

        for (i = 0; i < qMaxp2; i++) {
            factorinv[i] = 2*i-1;
            factor[i] = 1.0 / factorinv[i];
        }

        a2 = a * a;
        Ra[0] = a;
        Ra[1] = sqrt(L2+a2);
        Ra[2] = Ra[1];
        Ra[3] = a;

        Rainv[0] = 1.0 / Ra[1];

        Ra2inv = Rainv[0] * Rainv[0];

        for (i = 1; i < qMaxp1; i++) {
            Rainv[i] = Rainv[i-1] * Ra2inv;
        }

        Jyzt00[0] = log(a);
        Jyzt01[0] = log(Ra[1] - L);
        Jyzt02[0] = log(Ra[1] + L);
        Jyzt03[0] = Jyzt00[0];

        Jyzt10[0] = Ra[0];
        Jyzt11[0] = Ra[1];
        Jyzt12[0] = Ra[2];
        Jyzt13[0] = Ra[3];

        Rasum = Ra[0] + Ra[0] - Ra[1] - Ra[1];

        zJyzsum = L * (Jyzt01[0] - Jyzt03[0]);

        Hyt[0] = zJyzsum - Rasum;
        Ht[0] = L * (Jyzt01[0] - Jyzt02[0]) - Rasum - Rasum;

        local3 = 0;
        local4 = 0;

        reHy = delta * creal(Hyt[0]);
        imHy = delta * cimag(Hyt[0]);

        for (i = 0; i < 18; i++) {
            G[i] += reHy * FqReal[local3];
            local3++;
        }

#if 0
//
//      TEMP FIX BEFORE FULL REMOVAL OF COMPLEX VARIABLES
//
//      The new FqImag array does not contain the first 18 elements
//      because they were all zero, so this code is obsolete?
//
        for (i = 0; i < 18; i++) {
            G[i] -= imHy * FqImag[local4];
            local4++;
        }
#endif

        reHy = creal(alpha * Hyt[0]);
        imHy = cimag(alpha * Hyt[0]);

        for (i = 0; i < 18; i++) {
            G[i] += reHy * FqReal[local3];
            local3++;
        }

        for (i = 0; i < 18; i++) {
            G[i] -= imHy * FqImag[local4];
            local4++;
        }

/*
 *      Recursive integrals
 */
        indb = 1;
        local = 2;
        local2= 2;

        for (i = 2; i < qMaxp2; i++) {
            int           inda;
            complex8 cTmp1, cTmp2;

            inda = indb;
            indb = 1 - indb;

            kk = 2*i;
            ind = kk-4;

            Jyzt02[inda] = factor[i-1] * (complex8(ind+1,0)*Jyzt02[indb] -
                                          ypzin[ind+1][1] * Rainv[i-2]);

            Jyzt03[inda] = factor[i-1] * (ind+1) * Jyzt03[indb];

            cTmp1 = factor[i-1] * (2 * (i+1) - 4);
            cTmp2 = factor[i-1] * Rainv[i-2];

            ind = kk-2;

            Jyzt10[inda] = cTmp1*Jyzt10[indb];
            Jyzt11[inda] = cTmp1*Jyzt11[indb] - cTmp2*ypzin[ind][0];
            Jyzt12[inda] = cTmp1*Jyzt12[indb] - cTmp2*ypzin[ind][1];
            Jyzt13[inda] = cTmp1*Jyzt13[indb];

            cTmp1 = L * (Jyzt03[inda] - Jyzt02[inda]);
            Hyt[inda] = factor[i] * (Ht[indb] + complex8(ind,0) *
                                     Hyt[indb] - cTmp1);

            cTmp2 = Jyzt10[inda] - Jyzt11[inda] -
                    Jyzt12[inda] + Jyzt13[inda];
            Ht[inda] = factor[i] * (complex8(ind+1,0) * Ht[indb] - cTmp2);

            for (jj = 1; jj < RemArray[i]; jj++) {

                j = RelArray[local];
                k = kk - j;
                reHy = creal(alphaind[j-1] * deltaind[k] * Hyt[inda]);

                for (iii = 0; iii < Re18Array[local]; iii++) {
                    ii = ReiiArray[local3];
                    G[ii] += reHy*FqReal[local3];
                    local3++;
                }

                local++;  
            }

            for (jj = 1; jj < ImmArray[i]; jj++) {

                j = ImlArray[local2];
                k = kk - j;
                imHy = cimag(alphaind[j-1] * deltaind[k] * Hyt[inda]);

                for (iii = 0; iii < Im18Array[local2]; iii++) {
                    ii = ImiiArray[local4];
                    G[ii] -= imHy*FqImag[local4];
                    local4++;
                }

                local2++;
            }

        }  /* end for (i = 2; i < qMaxp2; i++) */

/*
 *      Compute pre-factor eps C C b t b t in the form of
 *      matrices Ap, B, Bp, A
 */
        bt11 = bx * t[0];
        bt22 = by * t[1];
        bt33 = bz * t[2];

        bt12 = bx * t[1];
        bt23 = by * t[2];
        bt31 = bz * t[0];

        bt21 = by * t[0];
        bt32 = bz * t[1];
        bt13 = bx * t[2];

        Ap12 = bt13;
        Ap62 = bt23;
        Ap51 = bt12;
        Ap23 = bt21;
        Ap43 = bt31;
        Ap31 = bt32;

        Ap63 = bt11 - bt22;
        Ap41 = bt22 - bt33;
        Ap52 = bt33 - bt11;

        B12 = bt23;
        B35 = bt32;
        B31 = bt12;
        B16 = bt13;
        B23 = bt31;
        B24 = bt21;

        B36 = bt22 - bt11;
        B14 = bt33 - bt22;
        B25 = bt11 - bt33;

        switch(home->param->materialType) {
            case MAT_TYPE_BCC:
            case MAT_TYPE_FCC:
            case MAT_TYPE_HCP:
/*
 *              For hexagonal or cubic elastic constant matrix
 */
                for (i = 0; i < 3; i++) {
                    Ap[i][0] = -ecMatrix[i][1]*Ap62 + ecMatrix[i][2]*Ap31;
                    Ap[i][1] =  ecMatrix[i][0]*Ap12 - ecMatrix[i][2]*Ap43;
                    Ap[i][2] = -ecMatrix[i][0]*Ap51 + ecMatrix[i][1]*Ap23;
                }

                Ap[3][0] =  ecMatrix[3][3]*Ap41;
                Ap[3][1] = -ecMatrix[3][3]*Ap23;
                Ap[3][2] =  ecMatrix[3][3]*Ap43;

                Ap[4][0] =  ecMatrix[4][4]*Ap51;
                Ap[4][1] =  ecMatrix[4][4]*Ap52;
                Ap[4][2] = -ecMatrix[4][4]*Ap31;

                Ap[5][0] = -ecMatrix[5][5]*Ap12;
                Ap[5][1] =  ecMatrix[5][5]*Ap62;
                Ap[5][2] =  ecMatrix[5][5]*Ap63;

                for (i = 0; i < 6; i++) {
                    for (j = 0; j < 3; j++) {
                        B[j][i] = -Ap[i][j];
                    }
                }
                break;

            case MAT_TYPE_RHOMBOHEDRAL_VA:
/*
 *              For general elastic constant matrix
 */
                for (i = 0; i < 6; i++) {

                    Ap[i][0] = -ecMatrix[i][1]*Ap62 + ecMatrix[i][2]*Ap31 +
                                ecMatrix[i][3]*Ap41 + ecMatrix[i][4]*Ap51 -
                                ecMatrix[i][5]*Ap12;

                    Ap[i][1] =  ecMatrix[i][0]*Ap12 - ecMatrix[i][2]*Ap43 -
                                ecMatrix[i][3]*Ap23 + ecMatrix[i][4]*Ap52 +
                                ecMatrix[i][5]*Ap62;

                    Ap[i][2] = -ecMatrix[i][0]*Ap51 + ecMatrix[i][1]*Ap23 +
                                ecMatrix[i][3]*Ap43 - ecMatrix[i][4]*Ap31 +
                                ecMatrix[i][5]*Ap63;

                    B[0][i] =  ecMatrix[1][i]*B12 - ecMatrix[2][i]*B35 +
                               ecMatrix[3][i]*B14 - ecMatrix[4][i]*B31 +
                               ecMatrix[5][i]*B16;

                    B[1][i] = -ecMatrix[0][i]*B16 + ecMatrix[2][i]*B23 +
                               ecMatrix[3][i]*B24 + ecMatrix[4][i]*B25 -
                               ecMatrix[5][i]*B12;

                    B[2][i] =  ecMatrix[0][i]*B31 - ecMatrix[1][i]*B24 -
                               ecMatrix[3][i]*B23 + ecMatrix[4][i]*B35 +
                               ecMatrix[5][i]*B36;
                }

                break;
#if 0
            case XXX:
/*
 *              For symetric elastic constant matrix
 */
                for (i = 0; i < 6; i++) {
                    Ap[i][0] = -ecMatrix[i][1]*Ap62 + ecMatrix[i][2]*Ap31 +
                                ecMatrix[i][3]*Ap41 + ecMatrix[i][4]*Ap51 -
                                ecMatrix[i][5]*Ap12;

                    Ap[i][1] =  ecMatrix[i][0]*Ap12 - ecMatrix[i][2]*Ap43 -
                                ecMatrix[i][3]*Ap23 + ecMatrix[i][4]*Ap52 +
                                ecMatrix[i][5]*Ap62;

                    Ap[i][2] = -ecMatrix[i][0]*Ap51 + ecMatrix[i][1]*Ap23 +
                                ecMatrix[i][3]*Ap43 - ecMatrix[i][4]*Ap31 +
                                ecMatrix[i][5]*Ap63;
                }

                for (i = 0; i < 6; i++) {
                    for (j = 0; j < 3; j++) {
                        B[j][i] = -Ap[i][j];
                    }
                }

                break;
#endif  /* end for symetric elastic constant matrix */

        }  /* end switch(home->param->materialType) */




        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                ApGaux[i][j] = 0.0;
                for (k = 0; k < 6; k++) {
                    offsetG = 6*j+k;
                    ApGaux[i][j] += Ap[k][i] * G[offsetG];
                }
            }
        }

        ApG[0] = ApGaux[0][0];
        ApG[1] = ApGaux[1][1];
        ApG[2] = ApGaux[2][2];

        ApG[3] = ApGaux[1][2] + ApGaux[2][1];
        ApG[4] = ApGaux[2][0] + ApGaux[0][2];
        ApG[5] = ApGaux[0][1] + ApGaux[1][0];

        for (i = 0; i < 3; i++) {
            f2[i] = 0.0;
            for (j = 0; j < 6; j++) {
                f2[i] += B[i][j] * ApG[j];
            }
        }

        f2[0] *= coeff12;
        f2[1] *= coeff12;
        f2[2] *= coeff12;

/*
 *      Include the core contribution to the force
 */

        f2[0] += f2Core[0];
        f2[1] += f2Core[1];
        f2[2] += f2Core[2];

        f1[0] = -f2[0];
        f1[1] = -f2[1];
        f1[2] = -f2[2];

        return;
}
#endif  /* end ifdef ANISOTROPIC */
