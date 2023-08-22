/**************************************************************************
 *
 *      Module:       MobilityLaw_BCC_Va_nl.c
 *      Description:  Contains functions for calculating isotropic
 *                    non-linear mobility of nodes in BCC Vanadium.
 *                    Includes non-linear components for screw and edge, 
 *                    although the edge nonlinearity in the non-glide
 *                    directions has been removed.
 *
 *
 *      Includes public functions:
 *            Mobility_BCC_Va_nl()
 *
 *     Includes private functions:
 *            EdgeDrag()
 *            JunctionDrag()
 *            ScrewDrag()
 *            FscaleEdge()
 *            FscaleScrew()
 *                
 *      Last Modified: 04/18/2008 gh - Original version based on the
 *                                     Arsenlis matlab code mobbccVanl3c.m
 *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)


//static int dumpDebugData = 0;

/*
 *      Define a number of the mobility-specific constants that
 *      have been derived from information from atomistic siumulations
 *      and declare related global variables that are used in
 *      multiple functions.
 */
#define ALPHA_SCREW 0.14
#define N_SCREW     5
#define ALPHA_EDGE  0.0
#define BETA_EDGE   2.0
#define N_EDGE      3

static real8 C0;
static real8 EpsC0;
static real8 DeltaH0Screw;
static real8 PeierlsScrew;
static real8 BetaScrew;
static real8 V0Screw;
static real8 LinearScrew;
static real8 BJglide;
static real8 DeltaH0Edge;
static real8 PeierlsEdge;
static real8 V0Edge;
static real8 LinearEdge;


/*
 *      Some of the array/matrix sizes are dependent on the number of
 *      segments connected to a node or the number of different types
 *      of screw segments attached to the node.  By defining a large
 *      enough upper limit on these numbers a node might have, we can
 *      allocate the arrays/matrices to a maximum size at the start.
 *      Makes life much simpler.
 */
#define MAX_NSCREW MAX_NBRS
#define MAX_NEDGE  MAX_NBRS

/*
 *      The mobility depends on a number of values that are dependent
 *      on the input parameters.  This function is just called once
 *      during initialization to set those values appropriately.
 */
void InitMob_BCC_Va_nl(Home_t *home) 
{
        real8   G0, P, burgMag, tempK;

        burgMag = home->param->burgMag;
        tempK   = home->param->TempK;
        P       = home->param->pressure;
        G0      = home->param->MobG0;

        C0           = 5077 / burgMag;
        EpsC0        = 1.0e-06 * 5077 / burgMag;
        DeltaH0Screw = BOLTZMANNS_CONST * (7891+(11292-6700*(P/G0))*(P/G0));
        PeierlsScrew = G0*(6.939e-3-(1.8786e-3-3.53e-3*(P/G0))*(P/G0));
        BetaScrew    = 1.9426e-2+(9.115e-4-1.3776e-7*tempK)*tempK;
        V0Screw      = (3.1e-2-(3.65e-5-2.303e-8*tempK)*tempK)*5077/burgMag;
        LinearScrew  = (1.8788e-2-(3.484e-6+1.0695e-8*tempK)*tempK) *
                       5077/burgMag;
        BJglide      = 100 / LinearScrew * PeierlsScrew;
        DeltaH0Edge  = (696+(P/G0)*(579-630*(P/G0)))*BOLTZMANNS_CONST;
        PeierlsEdge  = G0*(1.265e-3+(5.81e-4+7.11e-4*(P/G0))*(P/G0));
        V0Edge       = 0.005*5077/burgMag;
        LinearEdge   = (1.45e-2 - (2.53e-5*tempK) + (1.33e-8 * tempK * tempK)) *
                       5077 / burgMag;

        return;
}


/*
 *      Just another function to dump information used for debugging
 */
static void DumpNodeAndParams(Home_t *home, Node_t *node)
{
        int j;
        Node_t *nbrNode;

        printf("rTol         = %.15e\n", home->param->rTol);
        printf("nextDT       = %.15e\n", home->param->deltaTT);

        printf("Pressure     = %.15e\n", home->param->pressure);
        printf("G0           = %.15e\n", home->param->MobG0);
        printf("C0           = %.15e\n", C0);
        printf("EpsC0        = %.15e\n", EpsC0);

        printf("SPeierls     = %.15e\n", PeierlsScrew);
        printf("SDeltaH0     = %.15e\n", DeltaH0Screw);
        printf("SAlpha       = %.15e\n", ALPHA_SCREW);
        printf("SBeta        = %.15e\n", BetaScrew);
        printf("SLinear      = %.15e\n", LinearScrew);
        printf("SV0          = %.15e\n", V0Screw);
        printf("SN           = %d\n",    N_SCREW);
        printf("EPeierls     = %.15e\n", PeierlsEdge);
        printf("EDeltaH0     = %.15e\n", DeltaH0Edge);
        printf("EAlpha       = %.15e\n", ALPHA_EDGE);
        printf("EBeta        = %.15e\n", BETA_EDGE);
        printf("ELinear      = %.15e\n", LinearEdge);
        printf("EV0          = %.15e\n", V0Edge);
        printf("EN           = %d\n",    N_EDGE);
        printf("BJglide      = %.15e\n", BJglide);

        printf("Boltz        = %.15e\n", BOLTZMANNS_CONST);
        printf("shearModulus = %.15e\n", home->param->shearModulus);
        printf("poisson ratio= %.15e\n", home->param->pois);
        printf("Burgers mag. = %.15e\n", home->param->burgMag);
        printf("TempK        = %.15e\n", home->param->TempK);

        PrintNode(node);
        for (j = 0; j < node->numNbrs; j++) {
            if ((nbrNode = GetNeighborNode(home, node, j)) == (Node_t *)NULL) {
                continue;
            }
            PrintNode(nbrNode);
        }

        return;
}


static void JunctionDrag(real8 BJglide, real8 vel[3], real8 burg[3],
                         real8 t[3], real8 fjdrag[3], real8 dfjdragdv[3][3])
{
        real8 BJline;
        real8 tmp[3][3];

        BJline = BJglide * 1.0e-06;
/*
 *      dfjdragdv = beclimb.*(eye(3)-(t'*t)) + beline.*(t'*t)
 */
        Vec3TransposeAndMult(t, tmp);
        
        dfjdragdv[0][0] = BJglide * (1.0 - tmp[0][0]) + BJline * tmp[0][0];
        dfjdragdv[0][1] = BJglide * (0.0 - tmp[0][1]) + BJline * tmp[0][1];
        dfjdragdv[0][2] = BJglide * (0.0 - tmp[0][2]) + BJline * tmp[0][2];

        dfjdragdv[1][0] = BJglide * (0.0 - tmp[1][0]) + BJline * tmp[1][0];
        dfjdragdv[1][1] = BJglide * (1.0 - tmp[1][1]) + BJline * tmp[1][1];
        dfjdragdv[1][2] = BJglide * (0.0 - tmp[1][2]) + BJline * tmp[1][2];

        dfjdragdv[2][0] = BJglide * (0.0 - tmp[2][0]) + BJline * tmp[2][0];
        dfjdragdv[2][1] = BJglide * (0.0 - tmp[2][1]) + BJline * tmp[2][1];
        dfjdragdv[2][2] = BJglide * (1.0 - tmp[2][2]) + BJline * tmp[2][2];

        Matrix33Vector3Multiply(dfjdragdv, vel, fjdrag);

        return;
}


static void FscaleEdge(Param_t *param, real8 feinit, real8 vin, real8 burg[3],
                       real8 *fout, real8 *dfdv, real8 *dfdvlin, real8 *d2fdv2)
{
        int   notConverged;
        real8 error, errormag, edgeerrortol;
        real8 fmag, vmag, vsign, hratio, test;
        real8 peierls, deltaH0, tempK, beta, v0, Ml, n;
        real8 Hrstar,burgnormal;
        real8 ninv, nlmob;
        real8 vnl, dvnldf, d2vnldf2;
        real8 vlin, dvlindf, d2vlindf2;
        real8 vout, dvoutdf, d2voutdf2, dvoutdvnl, dvoutdvlin;
        real8 correction, temp1, temp2, fnlmax, dvnldfnlmax;
        real8 d2voutdvlin2, d2voutdvnl2, d2voutdvnldvlin;
        real8 fmin, fmax;


        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);
        deltaH0 = DeltaH0Edge;
        tempK = param->TempK;
        peierls = PeierlsEdge * burgnormal;
        beta = BETA_EDGE;
        v0 = V0Edge;
        Ml = LinearEdge;
        Ml = Ml / peierls;
        n = N_EDGE;

        vmag = fabs(vin);
        vsign = (real8)Sign(vin);
        hratio = deltaH0 / (BOLTZMANNS_CONST*tempK);
        ninv = 1.0 / n;

        edgeerrortol = MAX(EpsC0 * 1.0e-6, vmag * 1.0e-12);
        Hrstar = beta * hratio;
        nlmob = v0 * exp(-Hrstar);
        fnlmax = asinh(C0/nlmob) * peierls / Hrstar;
        dvnldfnlmax = nlmob * cosh(Hrstar*fnlmax/peierls) * Hrstar / peierls;

        if (vin == 0.0) {
            fmag = 0.0;
            dvoutdf = nlmob * Hrstar / peierls;
            d2voutdf2 = 0.0;
            error = 0.0;
        } else {
            if (feinit == 0.0) {
                fmag = 0.5 * peierls;
            } else {
                fmag = 0.666666 * fabs(feinit);
            }
/*
 *          Bound the solution to within a factor of two
 */
            fmax = -1;
            fmin = -1;

            while ((fmax < 0.0) || (fmin < 0.0)) {
                vlin = Ml * fmag;
                if (fmag < fnlmax) {
                    vnl = nlmob * sinh(Hrstar*fmag/peierls);
                } else {
                    vnl = C0 + dvnldfnlmax * (fmag-fnlmax);
                }
                temp1 = 1.0 / (pow(vlin+EpsC0, n) + pow(vnl,n));
                temp2 = pow(temp1, ninv);
                vout = (vlin+EpsC0) * vnl * temp2;
                error = vmag - vout;
                if (error < 0.0) {
                    fmax = fmag;
                    fmag = 0.5 * fmag;
                } else {
                    fmin = fmag;
                    fmag = 2.0 * fmag;
                }
            }

            fmag = fabs(feinit);

            if ((fmag >= fmax) || (fmag <= fmin)) {
/*
 *              feinit was not a good guess
 */
                fmag = 0.5 * (fmax+fmin);
            }
/*
 *          fmag now ranges between fmin and fmax
 */
            notConverged = 1;

            while (notConverged) {
                vlin = Ml * fmag;
                dvlindf = Ml;
                d2vlindf2 = 0.0;

                if (fmag < fnlmax) {
                    vnl = nlmob * sinh(Hrstar*fmag/peierls);
                    dvnldf = nlmob*cosh(Hrstar*fmag/peierls)*Hrstar/peierls;
                    d2vnldf2 = nlmob * sinh(Hrstar*fmag/peierls) *
                               (Hrstar*Hrstar) / (peierls*peierls);
                } else {
                    vnl = C0 + dvnldfnlmax * (fmag-fnlmax);
                    dvnldf = dvnldfnlmax;
                    d2vnldf2 = 0.0;
                }

                temp1 = 1.0 / (pow(vlin+EpsC0, n) + pow(vnl,n));
                temp2 = pow(temp1, ninv);
                vout = (vlin+EpsC0) * vnl * temp2;
                dvoutdvlin = vnl * temp2 * (1.0 - pow(vlin+EpsC0,n) * temp1);
                dvoutdvnl = (vlin+EpsC0) * temp2 * (1.0 - pow(vnl,n) * temp1);

                error = vmag - vout;
                if (error < 0.0) {
                    fmax = fmag;
                } else {
                    fmin = fmag;
                }

                dvoutdf = dvoutdvlin * dvlindf + dvoutdvnl * dvnldf;
                d2voutdvlin2 = dvoutdvlin * (n+1.0) * pow(vlin,n-1.0) * temp1;
                d2voutdvnl2 = dvoutdvnl * (n+1.0) * pow(vnl,n-1.0) * temp1;
                d2voutdvnldvlin = temp2 * ((1.0 - pow(vlin+EpsC0,n) * temp1) *
                                           (1.0 - pow(vnl,n) * temp1) +
                                           n * pow(vnl,n) * pow(vlin+EpsC0,n) *
                                           temp1 * temp1);
                d2voutdf2 = 2.0 * (d2voutdvlin2 * dvlindf +
                                   d2voutdvnldvlin * dvlindf * dvnldf +
                                   d2voutdvnl2 * dvnldf);
                d2voutdf2 += dvoutdvlin*d2vlindf2 +
                             dvoutdvnl*d2vnldf2;

                errormag = fabs(error);

                if (errormag > edgeerrortol) {
                    test = dvoutdf * dvoutdf + (2 * error * d2voutdf2);
                    if ((test < 0.0) ||
                        (fabs(d2voutdf2) < edgeerrortol / (peierls*peierls))) {
                        correction = error / dvoutdf;
                    } else {
                        correction = -(dvoutdf - sqrt(test)) / d2voutdf2;
                        if (correction == 0.0) {
                            correction = error / dvoutdf;
                        }
                    }
/*
 *                  if correction has pushed fmag outside the known range, use
 *                  a bisection instead
 */
                    fmag += correction;
                    if ((fmag >= fmax) || (fmag <= fmin)) {
                        fmag = 0.5 * (fmax + fmin);
                    }
                } else {
                    notConverged = 0;
                }
            }
        }

        *fout = fmag * vsign;
        *dfdv = 1.0 / dvoutdf;
        *dfdvlin = 1.0 / Ml;
#ifdef NAN_CHECK
        if (isnan(*fout) || isinf(*fout)) {
            Fatal("fscale: fout = %lf", *fout);
        }

        if (isnan(*dfdv) || isinf(*dfdv)) {
            Fatal("fscale: dfdv = %lf = 1.0 / dvoutdf\n"
                  "  where dvoutdf = %lf", *dfdv, dvoutdf);
        }

        if (isnan(*dfdvlin) || isinf(*dfdvlin)) {
            Fatal("fscale: dfdvlin = %lf = 1.0 / Ml\n"
                  "  where Ml = %lf", *dfdvlin, Ml);
        }
#endif

#if 0
/*
 *      Use of the 2nd derivatives is still somewhat problematic
 *      so for the time being this value will be zero out.
 */
        *d2fdv2 = -d2voutdf2 * pow(*dfdv, 3.0) * vsign;
#else
        *d2fdv2 = 0.0;
#endif

        return;
}


static void EdgeDrag(Param_t *param, real8 vel[3], real8 burg[3],
                     real8 climbdir[3],  real8 fedrag[3],
                     real8 dfedragdv[3][3])
{
        int   m, n;
        real8 vmag, fout, dfdv, d2fdv2;
        real8 Beclimb, Beline;
        real8 dp_climbvel, dp_linevel;
        real8 ftest;
        real8 linedir[3], glidedir[3], temp31[3];
        real8 tmp1[3][3], tmp2[3][3], tmp3[3][3];
        real8 temp33B[3][3];
        real8 dfdvlin;


        glidedir[0] = burg[0];
        glidedir[1] = burg[1];
        glidedir[2] = burg[2];

        Normalize(&glidedir[0], &glidedir[1], &glidedir[2]);

        cross(glidedir, climbdir, linedir);
    
        vmag = DotProduct(vel, glidedir);

        ftest = DotProduct(fedrag, glidedir);

        FscaleEdge(param, ftest, vmag, burg, &fout, &dfdv, &dfdvlin, &d2fdv2);

        Beclimb = dfdvlin * 1.0e+03;
        Beline = dfdv * 1.0e-03;

        Vec3TransposeAndMult(glidedir, tmp1);
        Vec3TransposeAndMult(climbdir, tmp2);
        Vec3TransposeAndMult(linedir,  tmp3);
        
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] = Beclimb * tmp2[m][n] +
                                  Beline  * tmp3[m][n]; 
            }
        }

        for (m = 0; m < 3; m++) {
            Matrix33Vector3Multiply(dfedragdv, vel, temp31);
            fedrag[m] = temp31[m] + fout * glidedir[m];
        }

        dp_climbvel = DotProduct(climbdir, vel);
        dp_linevel = DotProduct(linedir, vel);

        for (m = 0; m < 3; m++) {
            temp31[m] = d2fdv2 *
                        (1.0e+03 * dp_climbvel * climbdir[m] +
                         1.0e-03 * dp_linevel * linedir[m]);
        }

        Matrix33_VVt_Mul(temp33B,temp31,glidedir);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] += dfdv * tmp1[m][n] + temp33B[m][n];
            }
        }

        return;
}

static void FscaleScrew(Param_t *param, real8 fsinit, real8 vin, real8 burg[3],
                   real8 *fout, real8 *dfdv, real8 *d2fdv2)
{
        int   notConverged;
        real8 error, errormag, screwerrortol;
        real8 fmag, vmag, vsign, hratio, test;
        real8 peierls, deltaH0, tempK, alpha, beta, v0, Ml, n;
        real8 limit, Hrstar,burgnormal;
        real8 ninv, nlmob, nlexp;
        real8 vnl, dvnldf, d2vnldf2;
        real8 vlin, dvlindf, d2vlindf2;
        real8 vout, dvoutdf, d2voutdf2, dvoutdvnl, dvoutdvlin;
        real8 correction, temp1, temp2;
        real8 d2voutdvlin2, d2voutdvnl2, d2voutdvnldvlin;
        real8 fmin, fmax;

        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);
        deltaH0 = DeltaH0Screw;
        tempK = param->TempK;
        peierls = PeierlsScrew * burgnormal;
        alpha = ALPHA_SCREW;
        beta = BetaScrew;
        v0 = V0Screw;
        Ml = LinearScrew;
        Ml = Ml / peierls;
        n = N_SCREW;

        vmag = fabs(vin);
        vsign = (real8)Sign(vin);
        fmag = 0.5 * peierls;
        hratio = deltaH0 / (BOLTZMANNS_CONST*tempK);
        limit = 0.95 * alpha * peierls;
        ninv = 1.0 / n;

        screwerrortol = MAX(EpsC0 * 1.0e-6, vmag * 1.0e-12);
        Hrstar = beta * hratio * log(1.0+alpha) / log((1.0+alpha)/alpha);
        nlexp = Hrstar / log(1.0+alpha);
        nlmob = v0 * exp(-Hrstar);

        if (vin == 0.0) {
            fmag = 0.0;
            dvoutdf = nlmob * nlexp * 2.0 * pow(alpha, nlexp-1.0) / peierls;
            d2voutdf2 = 0.0;
            error = 0.0;
        } else {
            if (fsinit == 0.0) {
                fmag = 0.5 * peierls;
            } else {
                fmag = 0.666666 * fabs(fsinit);
            }
/*
 *          Bound the solution to within a factor of two
 */
            fmax = -1;
            fmin = -1;

            while ((fmax < 0.0) || (fmin < 0.0)) {
                vlin = Ml * fmag;
                if (fmag > limit) {
                    vnl = nlmob * (pow(fmag/peierls+alpha,nlexp) -
                                   pow(alpha-limit/peierls,nlexp));
                } else {
                    vnl = nlmob * (pow(alpha+fmag/peierls,nlexp) -
                                   pow(alpha-fmag/peierls,nlexp));
                }

                temp1 = 1.0 / (pow(vlin+EpsC0, n) + pow(vnl,n));
                temp2 = pow(temp1, ninv);
                vout = (vlin+EpsC0) * vnl * temp2;
                error = vmag - vout;
                if (error < 0.0) {
                    fmax = fmag;
                    fmag = 0.5 * fmag;
                } else {
                    fmin = fmag;
                    fmag = 2.0 * fmag;
                }
            }

            fmag = fabs(fsinit);
/*
 *          if fsinit was not a good guess...
 */
            if ((fmag >= fmax) || (fmag <= fmin)) {
                fmag = 0.5 * (fmax+fmin);
            }

            notConverged = 1;

            while (notConverged) {
                vlin = Ml * fmag;
                dvlindf = Ml;
                d2vlindf2 = 0.0;
                if (fmag > limit) {
                    vnl = nlmob * (pow(fmag/peierls+alpha,nlexp) -
                                   pow(alpha-limit/peierls,nlexp));
                    dvnldf = nlmob * nlexp *
                             pow(fmag/peierls+alpha,nlexp-1.0) / peierls;
                    d2vnldf2 = nlmob * nlexp * (nlexp-1.0) *
                               pow(fmag/peierls+alpha,nlexp-2.0) /
                               (peierls*peierls);
                } else {
                    vnl = nlmob * (pow(fmag/peierls+alpha,nlexp) -
                                   pow(alpha-fmag/peierls,nlexp));
                    dvnldf = nlmob * nlexp * 
                             (pow(fmag/peierls+alpha,nlexp-1.0) +
                              pow(alpha-fmag/peierls,nlexp-1.0)) / peierls;
                    d2vnldf2 = nlmob * nlexp * (nlexp-1.0) *
                               (pow(fmag/peierls+alpha,nlexp-2.0) -
                                pow(alpha-fmag/peierls,nlexp-2.0)) /
                               (peierls*peierls);
                }

                temp1 = 1.0 / (pow(vlin+EpsC0, n) + pow(vnl,n));
                temp2 = pow(temp1, ninv);
                vout = (vlin+EpsC0) * vnl * temp2;
                dvoutdvlin = vnl * temp2 * (1.0 - pow(vlin+EpsC0,n) * temp1);
                dvoutdvnl = (vlin+EpsC0) * temp2 * (1.0 - pow(vnl,n) * temp1);

                error = vmag - vout;
                if (error < 0.0) {
                    fmax = fmag;
                } else {
                    fmin = fmag;
                }

                dvoutdf = dvoutdvlin * dvlindf + dvoutdvnl * dvnldf;
                d2voutdvlin2 = dvoutdvlin * (n+1.0) * pow(vlin,n-1.0) * temp1;
                d2voutdvnl2 = dvoutdvnl * (n+1.0) * pow(vnl,n-1.0) * temp1;
                d2voutdvnldvlin = temp2 * ((1.0 - pow(vlin+EpsC0,n) * temp1) *
                                           (1.0 - pow(vnl,n) * temp1) +
                                           n * pow(vnl,n) * pow(vlin+EpsC0,n) *
                                           temp1 * temp1);
                d2voutdf2 = 2.0 * (d2voutdvlin2 * dvlindf +
                                   d2voutdvnldvlin * dvlindf * dvnldf +
                                   d2voutdvnl2 * dvnldf);
                d2voutdf2 += dvoutdvlin*d2vlindf2 + dvoutdvnl*d2vnldf2;

                errormag = fabs(error);

                if (errormag > screwerrortol) {
                    test = dvoutdf*dvoutdf + 2.0 * error * d2voutdf2;
                    if ((test < 0.0) ||
                        (fabs(d2voutdf2) < screwerrortol / (peierls*peierls))) {
                        correction = error / dvoutdf;
                    } else {
                        correction = -(dvoutdf - sqrt(test)) / d2voutdf2;
                        if (correction == 0.0) {
                            correction = error / dvoutdf;
                        }
                    }
/*
 *                  if correction has pushed fmag outside the known range, use
 *                  a bisection instead
 */
                    fmag += correction;
                    if ((fmag >= fmax) || (fmag <= fmin)) {
                        fmag = 0.5 * (fmax + fmin);
                    }

                } else {
                    notConverged = 0;
                }
            }
        }

        *fout = fmag * vsign;
        *dfdv = 1.0 / dvoutdf;
#if 0
/*
 *      Use of the 2nd derivatives is still somewhat problematic
 *      so for the time being this value will be zero out.  It may
 *      be that it is just the 2nd derivatives for the edges that
 *      are troublesome, though, so we may be able to use the
 *      screw values... need to look into it.
 */
        *d2fdv2 = -d2voutdf2 * pow(*dfdv,3) * vsign;
#else
        *d2fdv2 = 0.0;
#endif

        return;
}


/*
 *    finit:  Initial guess... updated before return to caller.
 */
static void ScrewDrag(Param_t *param, real8 finit[3], real8 vel[3],
                      real8 burg[3], real8 dfsdragdv[3][3])
{
        int   i, m, n;
        real8 eps;
        real8 Bsline;
        real8 burgnormal;
        real8 tmp3[3], tmp33[3][3];
        real8 vproj[3];
        real8 linedir[3], linedirTM[3][3];
        real8 glidedirmat[3][3];
        real8 vMag, vMagInv, vNorm, vdir[3];
        real8 fout, dfdv, d2fdv2;
        real8 dvdirdv[3][3];
        real8 dp_linevel;
        real8 ftest, fproj[3];


        eps = 1.0e-12;

        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);

        linedir[0] = burg[0] / burgnormal;
        linedir[1] = burg[1] / burgnormal;
        linedir[2] = burg[2] / burgnormal;

        Vec3TransposeAndMult(linedir, linedirTM);

        glidedirmat[0][0] = (1.0 - linedirTM[0][0]);
        glidedirmat[0][1] = (0.0 - linedirTM[0][1]);
        glidedirmat[0][2] = (0.0 - linedirTM[0][2]);
        
        glidedirmat[1][0] = (0.0 - linedirTM[1][0]);
        glidedirmat[1][1] = (1.0 - linedirTM[1][1]);
        glidedirmat[1][2] = (0.0 - linedirTM[1][2]);
        
        glidedirmat[2][0] = (0.0 - linedirTM[2][0]);
        glidedirmat[2][1] = (0.0 - linedirTM[2][1]);
        glidedirmat[2][2] = (1.0 - linedirTM[2][2]);

        Matrix33Vector3Multiply(glidedirmat, vel, vproj);

        vMag = sqrt(DotProduct(vproj, vproj));
        vMagInv = 1.0 / (vMag + eps);

        for (i = 0; i < 3; i++) {
            vdir[i] = vproj[i] * vMagInv;
        }

        Vec3TransposeAndMult(vdir, tmp33);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dvdirdv[m][n] = (glidedirmat[m][n]-tmp33[m][n]) * vMagInv;
            }
        }

        fout = 0.0;
        dfdv = 0.0;
        d2fdv2 = 0.0;

        Matrix33Vector3Multiply(glidedirmat, finit, fproj);
        ftest = sqrt(DotProduct(fproj, fproj));

        FscaleScrew(param, ftest, vMag, burg, &fout, &dfdv, &d2fdv2);

        Bsline = dfdv * 1.0e-03;

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfsdragdv[m][n] = Bsline * linedirTM[m][n];
            }
        }

        Matrix33Vector3Multiply(dfsdragdv, vel, tmp3);

        for (i = 0; i < 3; i++) {
            finit[i] = tmp3[i] + fout*vdir[i];
        }

        vNorm = Normal(vel);

        if ((vNorm > eps) && ((vMag/vNorm) > eps)) {
            Vec3TransposeAndMult(vdir, tmp33);
            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dfsdragdv[m][n] +=  dfdv * tmp33[m][n] +
                                        fout * dvdirdv[m][n];
                }
            }
        } else {
            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dfsdragdv[m][n] +=  dfdv * ((m == n) - linedirTM[m][n]);
                }
            }
        }

        dp_linevel = DotProduct(linedir, vel);

        Matrix33_VVt_Mul(tmp33,linedir, vdir);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfsdragdv[m][n] += d2fdv2 * 1.0e-03 * dp_linevel * tmp33[m][n];
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     Mobility_BCC_Va_nl
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node         Pointer to the node for which to
 *                            calculate velocity
 *          IN/OUT: mobArgs   Structure containing additional
 *                            parameters for conveying information
 *                            to/from the mobility function.
 * 
 *      Returns: 0 if successful
 *               1 if function could not converge on a velocity
 *
 *************************************************************************/
int Mobility_BCC_Va_nl(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, j, m, n, iAll, numNbrs, numscrews, thisScrew;
        int     iterCnt, maxIter, maxIter2, count;
        int     numedges, thisEdge;
        int     zeroLenSeg[MAX_NBRS];
        int     edgeExists[MAX_NBRS], screwExists[MAX_NBRS];
        int     notConverged;
        int     numNonZeroLenSegs = 0;
        real8   velerrortol, velerror, minvelerror;
        real8   mag1, mag2;
        real8   eps, rt2, burgProd, mag, avgFactor, correctionAdjustment;
        real8   vtest[3], fn[3], rt[3], zero3[3], tmp3[3], ferror[3];
        real8   dfsdv[3][3], dferrordv[3][3];
        real8   dflindv[3][3], dfjdragdv[3][3];
        real8   dferrordvold[3][3], dfedv[3][3]; 
        real8   LTotal    [MAX_NBRS  ]   ={ 0.0 }, LScrew    [MAX_NBRS  ]   ={ 0.0 }, LEdge[MAX_NBRS]   ={ 0.0 }, L2;
        real8   linedir   [MAX_NBRS  ][3]={{0.0}}, ndir      [MAX_NBRS  ][3]={{0.0}}, burg [MAX_NBRS][3]={{0.0}};
        real8   screwb    [MAX_NSCREW][3]={{0.0}}, dfndfscrew[MAX_NSCREW]   ={ 0.0 };
        real8   fscrew    [MAX_NSCREW][3]={{0.0}}, edgeb     [MAX_NEDGE ][3]={{0.0}};
        real8   edgenorm  [MAX_NEDGE ][3]={{0.0}}, dfndfedge [MAX_NEDGE ]   ={ 0.0 };
        real8   fedge     [MAX_NEDGE ][3]={{0.0}};
        real8   correction[3]={0.0};
        real8   vtmp[3]={0.0, 0.0, 0.0}, fjdrag[3];
        real8   mu, lmax;
        real8   forceerror, ferrortol;
        real8   totLen, massMult=0.0, massMatrix[3][3], vOld[3], vDiff[3], tmp3B[3];
        Param_t *param;
        Node_t  *nbrNode;


        eps = 1.0e-12;
        iterCnt = 0;
        maxIter = 50;
        maxIter2 = 100;

        real8 (*invdferrordv)[3]     =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdferrordv);

/*
 *      This mobility module imposes no glide constraints on the node.
 */
        *numGlideConstraints = 0;

        param = home->param;
        mu = param->shearModulus;
        lmax = param->maxSeg;

        memset(zero3, 0, sizeof(zero3));
/*
 *      If we need to include inertial terms, need some more initializations
 */
        if (param->includeInertia != 0) {
            massMult = 0.25 * param->massDensity *
                       (param->burgMag * param->burgMag);
        }

/*
 *      If the node is pinned in ALL dimensions it cannot be moved
 *      so just zero the velocity and return
 */
        if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      If any of the arms of the node has a burgers vector that
 *      has explicitly been set to be sessile (via control file
 *      inputs), the node may not be moved.
 */
        if (NodeHasSessileBurg(home, node)) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

        numNbrs = node->numNbrs;

/*
 *      This function saves various pieces of info for each neighbor
 *      node using static arrays for speed rather than doing memory
 *      allocations.  If the arrays are too  small, abort with an error
 *      to let the caller increase the array sizes.
 */
        if (numNbrs > MAX_NBRS) {
            Fatal("Mobility_BCC_Va_nl: Node segment count (%d) exceeds\n"
                  "maximum allowed (%d) in static arrays.  Increase MAX_NBRS\n"
                  "and recompile\n", node->numNbrs, MAX_NBRS);
        }

/*
 *      Initialize some arrays for later...
 */
        vtest[0] = node->vX;
        vtest[1] = node->vY;
        vtest[2] = node->vZ;

        vOld[0] = vtest[0];
        vOld[1] = vtest[1];
        vOld[2] = vtest[2];

        memset(ndir, 0, sizeof(ndir));

/*
 *      Make a temporary copy of the total forces
 */
        fn[0] = node->fX;
        fn[1] = node->fY;
        fn[2] = node->fZ;

/*
 *      If necessary, rotate the force and velocity vectors from the
 *      laboratory frame to the crystal frame
 */
        if (param->useLabFrame) {
            real8 fnRot[3], vtestRot[3];
            Matrix33Vector3Multiply(param->rotMatrixInverse, fn, fnRot);
            Matrix33Vector3Multiply(param->rotMatrixInverse, vtest, vtestRot);
            VECTOR_COPY(fn, fnRot);
            VECTOR_COPY(vtest, vtestRot);
            VECTOR_COPY(vOld, vtest);
        }

#if 0
/*
 *      For debug only...
 */
{
        real8   x[4], y[4], z[4];

        home->param->deltaTT = 3.716023500000000e-14;

        numNbrs = 3;

        fn[0] = -1.583682898404573e+10;
        fn[1] = -8.048410323907104e+09;
        fn[2] = -5.252440014827087e+10;

        vtest[0] =  2.099470451717293e+11;
        vtest[1] =  1.438929709467170e+11;
        vtest[2] = -2.420183864050443e+11;

        x[0] = -2.144137264335070e+03;
        y[0] =  3.466244481031347e+03;
        z[0] =  2.907975450800293e+03;

        x[1] = -1.981221517269450e+03;
        y[1] =  3.623083948402581e+03;
        z[1] =  2.754331630921930e+03;

        x[2] = -2.142823470470344e+03;
        y[2] =  3.455974819109506e+03;
        z[2] =  2.837428798156832e+03;

        x[3] = -2.147730498877328e+03;
        y[3] =  3.463782347197094e+03;
        z[3] =  2.912119965662012e+03;

        burg[0][0] = -5.773503000000000e-01;
        burg[0][1] = -5.773503000000000e-01;
        burg[0][2] =  5.773503000000000e-01;

        burg[1][0] =  5.773503000000000e-01;
        burg[1][1] =  5.773503000000000e-01;
        burg[1][2] =  5.773503000000000e-01;

        burg[2][0] =  0.000000000000000e+00;
        burg[2][1] =  0.000000000000000e+00;
        burg[2][2] = -1.154700600000000e+00;


        for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {
            rt[0] = x[i+1] - x[0];
            rt[1] = y[i+1] - y[0];
            rt[2] = z[i+1] - z[0];
            ZImage(param, &rt[0], &rt[1], &rt[2]);
            if ((rt2 = DotProduct(rt, rt)) < eps) {
                i--;
                zeroLenSeg[iAll] = 1;
                continue;
            } else {
                zeroLenSeg[iAll] = 0;
            }
            
            LTotal[i] = sqrt(rt2);

            LScrew[i] = fabs(rt[0]*burg[i][0] +
                             rt[1]*burg[i][1] +
                             rt[2]*burg[i][2]) / Normal(burg[i]);

            LScrew[i] = MIN(LTotal[i], LScrew[i]);
            L2 = MAX(0.0, (LTotal[i]*LTotal[i]-LScrew[i]*LScrew[i]));

            if (L2 <= 0.0) {
                LEdge[i] = 0.0;
            } else {
                LEdge[i] = sqrt(L2);
            }

            linedir[i][0] = rt[0]/LTotal[i];
            linedir[i][1] = rt[1]/LTotal[i];
            linedir[i][2] = rt[2]/LTotal[i];

            edgeExists[i]  = ((LEdge[i]  / LTotal[i]) > 1.0e-06);
            screwExists[i] = ((LScrew[i] / LTotal[i]) > 1.0e-06);

            if (edgeExists[i]) {
                cross(burg[i], linedir[i], ndir[i]);
                Normalize(&ndir[i][0], &ndir[i][1], &ndir[i][2]);
            }
        }
}
#else
/*
 *      Loop over all arms attached to the node
 */
        totLen = 0.0;

        for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {

            nbrNode = GetNeighborNode(home, node, iAll);

            if (nbrNode == (Node_t *)NULL) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

            rt[0] = nbrNode->x - node->x;
            rt[1] = nbrNode->y - node->y;
            rt[2] = nbrNode->z - node->z;
            
            ZImage(param, &rt[0], &rt[1], &rt[2]);

/*
 *          If the segment is zero length (which may happen when the
 *          mobility function is called from SplitMultiNodes()), just
 *          skip the segment.
 */
            if ((rt2 = DotProduct(rt, rt)) < eps) {
                i--;
                zeroLenSeg[iAll] = 1;
                continue;
            } else {
                zeroLenSeg[iAll] = 0;
                numNonZeroLenSegs++;
            }

            burg[i][0] = node->burgX[iAll];
            burg[i][1] = node->burgY[iAll];
            burg[i][2] = node->burgZ[iAll];

/*
 *          If necessary, rotate the burgers vector and line sense from the
 *          laboratory frame to the crystal frame
 */
            if (param->useLabFrame) {
                real8 burgRot[3], rtRot[3];

                Matrix33Vector3Multiply(param->rotMatrixInverse,burg[i],burgRot);
                Matrix33Vector3Multiply(param->rotMatrixInverse,rt,rtRot);

                VECTOR_COPY(burg[i], burgRot);
                VECTOR_COPY(rt, rtRot);
            }

/*
 *          Calculate the segment lengths (total length, screw length,
 *          edge length)
 */
            LTotal[i] = sqrt(rt2);

            LScrew[i] = fabs(rt[0]*burg[i][0] +
                             rt[1]*burg[i][1] +
                             rt[2]*burg[i][2]) / Normal(burg[i]);

/*
 *          Sanity check.  Due to the different calculations of LTotal
 *          and LScrew above, it is possible (due to machine roundoff
 *          errors) for LScrew to be larger than LTotal.  If that
 *          happens, we need to correct the situation or all hell
 *          breaks loose with subsequent calculations... 
 */
            LScrew[i] = MIN(LTotal[i], LScrew[i]);
            L2 = MAX(0.0, (LTotal[i]*LTotal[i]-LScrew[i]*LScrew[i]));

            if (L2 <= 0.0) {
                LEdge[i] = 0.0;
            } else {
                LEdge[i] = sqrt(L2);
            }

/*
 *          Calculate tangents
 */
            linedir[i][0] = rt[0]/LTotal[i];
            linedir[i][1] = rt[1]/LTotal[i];
            linedir[i][2] = rt[2]/LTotal[i];

            edgeExists[i]  = ((LEdge[i]  / LTotal[i]) > 1.0e-06);
            screwExists[i] = ((LScrew[i] / LTotal[i]) > 1.0e-06);

            if (edgeExists[i]) {
                cross(burg[i], linedir[i], ndir[i]);
                Normalize(&ndir[i][0], &ndir[i][1], &ndir[i][2]);
            }

            totLen += LTotal[i];

        }  /* loop over neighbors */
#endif

/*
 *      It's possible this function was called for a node which only
 *      had zero length segments (during SplitSurfaceNodes() for example).
 *      If that is the case, just set the velocity to zero and return;
 */
        if (numNonZeroLenSegs == 0) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      Set up the mass matrix: if not needed, zero it out.
 */
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                massMatrix[m][n] = 0.0;
                if ((m == n) && (param->includeInertia)) {
                    massMatrix[m][n] = totLen * massMult / param->deltaTT;
                }
            }
        }

        memset(dflindv, 0, sizeof(dflindv));

        for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {

/*
 *          Skip any zero length segments...
 */
            if (zeroLenSeg[iAll]) {
                i--;
                continue;
            }

            burgProd = burg[i][0] * burg[i][1] * burg[i][2];

            if (fabs(burgProd) < eps) {
                JunctionDrag(BJglide, zero3, burg[i], linedir[i],
                             fjdrag, dfjdragdv);
                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        dflindv[m][n] += (0.5 * LTotal[i]) * dfjdragdv[m][n];
                    }
                }
            }

        }  /* Loop over all segments */

/*
 *      Loop over all segments (skipping zero length segments)
 */
        numscrews = 0;
        numedges = 0;

        for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {

            if (zeroLenSeg[iAll]) {
                i--;
                continue;
            }

            burgProd = burg[i][0] * burg[i][1] * burg[i][2];

            if (fabs(burgProd) > eps) {

                if (screwExists[i]) {
                    thisScrew = -1;
                    for (j = 0; j < numscrews; j++) {
                        mag = ((burg[i][0] * screwb[j][0] +
                                burg[i][1] * screwb[j][1] +
                                burg[i][2] * screwb[j][2]) /
                               (Normal(burg[i]) * Normal(screwb[j])));
                        if (fabs(fabs(mag) - 1.0) < eps) {
                            thisScrew = j;
                        }
                    }

                    if (thisScrew < 0) {
                        screwb[numscrews][0] = burg[i][0];
                        screwb[numscrews][1] = burg[i][1];
                        screwb[numscrews][2] = burg[i][2];
                        dfndfscrew[numscrews] = 0.5*LScrew[i];
                        numscrews++;
                    } else {
                        dfndfscrew[thisScrew] += 0.5*LScrew[i];
                    }
                }

                if (edgeExists[i]) {
                    thisEdge = -1;
                    for (j = 0; j < numedges; j++) {
                        mag1 = DotProduct(ndir[i], edgenorm[j]);
                        mag2 = ((burg[i][0] * edgeb[j][0] +
                                 burg[i][1] * edgeb[j][1] +
                                 burg[i][2] * edgeb[j][2]) /
                                (Normal(burg[i]) * Normal(edgeb[j])));
                        if ((fabs(fabs(mag1) - 1.0) < eps) &&
                            (fabs(fabs(mag2) - 1.0) < eps)) {
                            thisEdge = j;
                        }
                    }
                    if (thisEdge < 0) {
                        edgeb[numedges][0] = burg[i][0];
                        edgeb[numedges][1] = burg[i][1];
                        edgeb[numedges][2] = burg[i][2];
                        edgenorm[numedges][0] = ndir[i][0];
                        edgenorm[numedges][1] = ndir[i][1];
                        edgenorm[numedges][2] = ndir[i][2];
                        dfndfedge[numedges] = 0.5 * LEdge[i];
                        numedges++;
                    } else {
                        dfndfedge[thisEdge] += 0.5 * LEdge[i];
                    }
                }
            }
        }

        count = 0;
        memset(fscrew, 0, sizeof(fscrew));
        memset(fedge, 0, sizeof(fedge));
        memset(dferrordvold, 0, sizeof(dferrordvold));
        ferrortol = MAX((DotProduct(fn,fn)*1.0e-16),
                        (mu*mu*lmax*lmax*1.0e-24));
        velerrortol = MAX((1.0e-14*param->rTol*param->rTol/
                           (param->deltaTT*param->deltaTT)),
                           (EpsC0*EpsC0*1.0e-20));
        minvelerror = velerrortol;
        correctionAdjustment = 1.0;
        notConverged = 1;

        while (notConverged) {

            count++;

            Matrix33Vector3Multiply(dflindv, vtest, tmp3);

/*
 *          If we need to include inertial terms...
 */
            if (param->includeInertia) {
                vDiff[0] = vtest[0] - vOld[0];
                vDiff[1] = vtest[1] - vOld[1];
                vDiff[2] = vtest[2] - vOld[2];
                Matrix33Vector3Multiply(massMatrix, vDiff, tmp3B);
            } else {
                tmp3B[0] = 0.0;
                tmp3B[1] = 0.0;
                tmp3B[2] = 0.0;
            }

            for (i = 0; i < 3; i++) {
                ferror[i] = fn[i] - tmp3[i] - tmp3B[i];
            }


            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dferrordv[m][n] = -dflindv[m][n] - massMatrix[m][n];
                }
            }

            for (i = 0; i < numscrews; i++) {
                ScrewDrag(param, fscrew[i], vtest, screwb[i], dfsdv);
                for (m = 0; m < 3; m++) {
                    ferror[m] -= dfndfscrew[i] * fscrew[i][m];
                    for (n = 0; n < 3; n++) {
                        dferrordv[m][n] -= dfndfscrew[i] * dfsdv[m][n];
                    }
                }
            }

            for (i = 0; i < numedges; i++) {
                EdgeDrag(param, vtest, edgeb[i], edgenorm[i], fedge[i], dfedv);
                for (m = 0; m < 3; m++) {
                    ferror[m] -= dfndfedge[i] * fedge[i][m];
                    for (n = 0; n < 3; n++) {
                        dferrordv[m][n] -= dfndfedge[i] * dfedv[m][n];
                    }
                }
            }

/*
 *          Calculate the velocity correction.  Note: if we cannot
 *          invert dferrordv, just return an error
 */
            if (count > 1) {
                avgFactor = 0.50;
                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        dferrordv[m][n] = avgFactor * dferrordv[m][n] +
                                         (1.0-avgFactor) * dferrordvold[m][n];
                    }
                }
            }

            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dferrordvold[m][n] = dferrordv[m][n];
                }
            }

            if ( M33_INVERSE(invdferrordv,dferrordv) < 0 ) 
            {
                printf("%s::%s(%d) : Cannot invert matrix dferrordv!\n", __FILE__, __func__, __LINE__ );
                printf("iteration = %d\n", count);
                printf("numscrews = %d\n", numscrews);
                printf("dferrordv = \n    ");

                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        printf("%.13e ", dferrordv[m][n]);
                    }
                    printf("\n    ");
                }
                printf("\n");

                DumpNodeAndParams(home, node);

                return(1);
            }

            Matrix33Vector3Multiply(invdferrordv, ferror, correction);

/*
 *          Adjust the correction values by a factor.  This factor
 *          is initially 1.0, but in some cases, the function is 
 *          unable to converge on an answer because the corrrection
 *          oscillates between a couple values.  In those cases,
 *          the adjustment factor will be lowered, which can help
 *          achieve convergence.
 */
            for (i = 0; i < 3; i++) {
                correction[i] *= correctionAdjustment;
            }

            forceerror = DotProduct(ferror, ferror);
            velerror = DotProduct(correction, correction);
            if (forceerror > ferrortol) {
                for (i = 0; i < 3; i++) {
                    vtest[i] -= correction[i];
                }
            } else {
                notConverged = 0;
            }

/*
 *          If the initial velocity guess was too far off, it's
 *          better to just zero the vtest value.  So, after the
 *          first iteration, if the force error is too large (an
 *          indication the initial velocity guess was way off),
 *          just zero the vtest value.
 */
            if ((count == 1) && (forceerror > 1.0e+01 * ferrortol)) {
                for (i = 0; i < 3; i++) {
                    vtest[i] = 0.0;
                }
            }

/*
 *          Preserve the vtest value for the iteration that has the
 *          lowest absolute error.  If the function fails to converge
 *          we may be able to use this preserved velocity.
 */
            if (velerror < minvelerror) {
                minvelerror = velerror;
                vtmp[0] = vtest[0];
                vtmp[1] = vtest[1];
                vtmp[2] = vtest[2];
            }

            if (++iterCnt > maxIter) {
/*
 *              We didn't converge on an answer, but if there was an iteration
 *              of the above loop that resulted in an absolute error within
 *              the absolute error tolerance, go ahead and use the velocity
 *              calculated during that iteration.
 */
                if (minvelerror < velerrortol) {
                    vtest[0] = vtmp[0];
                    vtest[1] = vtmp[1];
                    vtest[2] = vtmp[2];
                    break;
                }

/*
 *              If the function failed to converge on a velocity in the
 *              first set of iterations, try a second set of iterations
 *              adjusting the correction values by a factor of 0.5.  This
 *              helps in some cases and although it increases the cost of
 *              calculating the mobility for a single node (on occassion),
 *              it can keep the timestep from being cut, which more than
 *              compensates.
 */
                if (iterCnt < maxIter2) {
                    maxIter = maxIter2;
                    count = 0;
                    memset(fscrew, 0, sizeof(fscrew));
                    memset(fedge, 0, sizeof(fedge));
                    memset(dferrordvold, 0, sizeof(dferrordvold));
                    velerror = 10.0 * velerrortol;
                    minvelerror = velerror;
                    correctionAdjustment = 0.5;
                    vtest[0] = node->vX;
                    vtest[1] = node->vY;
                    vtest[2] = node->vZ;
/*
 *                  If necessary, rotate the old velocity vector from the
 *                  laboratory frame to the crystal frame
 */
                    if (param->useLabFrame) {
                        real8 vtestRot[3];
                        Matrix33Vector3Multiply(param->rotMatrixInverse, vtest,
                                                vtestRot);
                        VECTOR_COPY(vtest, vtestRot);
                    }

                    continue;
                }
#if 0
/*
 *              Debug code only.
 */
                if ((home->param->deltaTT < 1.0e-16) && dumpDebugData) {
                    printf("mobility: (%d,%d) no convergence after %d iterations\n", 
                           node->myTag.domainID, node->myTag.index, iterCnt);
                    printf("forceerror = %e, correction = %e %e %e\n",
                           forceerror, tmp3[0], tmp3[1], tmp3[2]);
                    printf("velerror = %e, velerrortol = %e\n",
                           velerror, velerrortol);
                    printf("vtest = %e %e %e\n", vtest[0], vtest[1], vtest[2]);
  
                    DumpNodeAndParams(home, node);
                }
#endif
                return(1);
            }

        }  /* while (!converged) */

/*
 *      We were able to converge on a velocity or at least find one
 *      within the error tolerances, so set the nodal velocity
 *      to the new values.
 *
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame first
 */
        if (param->useLabFrame) {
            real8 vtestRot[3];
            Matrix33Vector3Multiply(param->rotMatrix, vtest, vtestRot);
            VECTOR_COPY(vtest, vtestRot);
        }

        node->vX = vtest[0];
        node->vY = vtest[1];
        node->vZ = vtest[2];

        return(0);
}
