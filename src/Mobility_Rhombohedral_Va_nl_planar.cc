/**************************************************************************
 *
 *      Module:       MobilityLaw_Rhombohedral_Va_nl_planar.c
 *      Description:  Based on the Arsenlis matlab code mobrhombVanlplanar.m
 *
 *                    There are three types of slip systems in this geometry:
 *                      - the long burgers vector on the  < 1  1 0> planes
 *                      - the short burgers vector on the < 1  1 0> planes
 *                      - the short burgers vector on the <10 10 1> planes
 *
 *                    There are also a number of junction burgers vectors of
 *                    the <20 1 1> type, <1 1 0> type and the <10 10 1> type.
 *                    The junctions will be allowed to move only along their
 *                    line direction.  The glide dislocations will only be
 *                    allowed to move along their glide planes.
 *
 *
 *      Includes public functions:
 *            MobilityLaw_Rhombohedral_Va_nl_planar()
 *
 *     Includes private functions:
 *            FastDrag()
 *            FscaleFast()
 *            EdgeDrag()
 *            JunctionDrag()
 *            ScrewDrag()
 *            FscaleEdge()
 *            FscaleScrew()
 *                
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
 *      Declare some variables that are specific to this module, but needed
 *      by various functions within the module.
 *
 *      There are three components of the 'screw' related constants
 *      PeierlsScrew, AlphaScrew, BetaScrew, GammaScrew, DeltaScrew,
 *      PhiScrew, EpsilonScrew, and EtaScrew.
 *
 *      There are also three components of PeierlsEdge, AlphaEdge,
 *      BetaEdge, GammaEdge, DeltaEdge, PhiEdge, EpsilonEdge, and
 *      EtaEdge.
 *
 *      The *Edge variables correspond to the non screw trench dislocations
 *
 *      The first element of these variables will apply to the long burgers
 *      vectors on the <1 1 0> planes
 *
 *      The second element of these variables will apply to the short burgers
 *      vectors on the <1 1 0> planes
 *
 *      The third element of these variables will apply to the short burgers
 *      vectors on the non <1 1 0> planes
 */
static real8 PeierlsScrew[3];
static real8 AlphaScrew[3];
static real8 BetaScrew[3];
static real8 GammaScrew[3];
static real8 DeltaScrew[3];
static real8 PhiScrew[3];
static real8 EpsilonScrew[3];
static real8 EtaScrew[3];

static real8 PeierlsEdge[3];
static real8 AlphaEdge[3];
static real8 BetaEdge[3];
static real8 GammaEdge[3];
static real8 DeltaEdge[3];
static real8 PhiEdge[3];
static real8 EpsilonEdge[3];
static real8 EtaEdge[3];

static real8 a, b, c, d, e, f;
static real8 thetaCrit[12][2];
static real8 C0;
static real8 epsC0;
static real8 BJglide;


/*
 *      Geometry of the rhombohedral slip systems burgers vectors in
 *      the rhombohedral geometry.  We'll normalize these during
 *      initialization.
 */
static real8 slipSystem[6][9];
static real8 burgRef[13][3] = {
     /* glide burgers vectors */
   { 1.45000000000e+00,  1.45000000000e+00,  1.45000000000e+00},
   { 1.31850000000e+00,  1.31850000000e+00, -1.18700000000e+00},
   { 1.31850000000e+00, -1.18700000000e+00,  1.31850000000e+00},
   {-1.18700000000e+00,  1.31850000000e+00,  1.31850000000e+00},
     /* likely junction burgers vectors */
   { 1.31500000000e-01,  1.31500000000e-01,  2.63700000000e+00},
   { 1.31500000000e-01,  2.63700000000e+00,  1.31500000000e-01},
   { 2.63700000000e+00,  1.31500000000e-01,  1.31500000000e-01},
     /* unlikely junction burgers vectors */
   {                 0,  2.50550000000e+00, -2.50550000000e+00},
   { 2.50550000000e+00,                  0, -2.50550000000e+00},
   { 2.50550000000e+00, -2.50550000000e+00,                  0},
   { 2.76850000000e+00,  2.76850000000e+00,  2.63000000000e-01},
   { 2.76850000000e+00,  2.63000000000e-01,  2.76850000000e+00},
   { 2.63000000000e-01,  2.76850000000e+00,  2.76850000000e+00}
};



/*
 *      The mobility depends on a number of values that are dependent
 *      on the input parameters.  This function is just called once
 *      during initialization to set those values appropriately.
 */
void InitMob_Rhombohedral_Va_nl_planar(Home_t *home)
{
        int     i;

        Param_t *param = home->param;

     // real8 MU          = param->shearModulus;
        real8 temperature = param->TempK;
        real8 burgmag     = param->burgMag;

/*
 *      Sound speed should be in unit of burgMag
 */ 
        C0    = 3300.0/burgmag;

/*
 *      Espilon value for use in covergence checking (in b/sec)
 */
        epsC0 = 1.0e-6 * C0;

/*
 *      Normalize the components of the burgRef[] array.
 */
        for (i = 0; i < 13; i++) {
            NormalizeVec(burgRef[i]);
        }

/*
 *      Initialize slip planes and close packed slip direction in the
 *      rhombohedral geometry.  Normalize slip planes.
 */
        slipSystem[0][0] = -3.63297500000e+00;
        slipSystem[0][1] =  3.63297500000e+00;
        slipSystem[0][2] =  0.00000000000e+00;
        slipSystem[0][3] =  burgRef[0][0];
        slipSystem[0][4] =  burgRef[0][1];
        slipSystem[0][5] =  burgRef[0][2];
        slipSystem[0][6] =  burgRef[1][0];
        slipSystem[0][7] =  burgRef[1][1];
        slipSystem[0][8] =  burgRef[1][2];

        slipSystem[1][0] =  3.63297500000e+00;
        slipSystem[1][1] =  0.00000000000e+00;
        slipSystem[1][2] = -3.63297500000e+00;
        slipSystem[1][3] =  burgRef[0][0];
        slipSystem[1][4] =  burgRef[0][1];
        slipSystem[1][5] =  burgRef[0][2];
        slipSystem[1][6] =  burgRef[2][0];
        slipSystem[1][7] =  burgRef[2][1];
        slipSystem[1][8] =  burgRef[2][2];

        slipSystem[2][0] =  0.00000000000e+00;
        slipSystem[2][1] = -3.63297500000e+00;
        slipSystem[2][2] =  3.63297500000e+00;
        slipSystem[2][3] =  burgRef[0][0];
        slipSystem[2][4] =  burgRef[0][1];
        slipSystem[2][5] =  burgRef[0][2];
        slipSystem[2][6] =  burgRef[3][0];
        slipSystem[2][7] =  burgRef[3][1];
        slipSystem[2][8] =  burgRef[3][2];

        slipSystem[3][0] =  3.29473250000e-01;
        slipSystem[3][1] = -3.30350175000e+00;
        slipSystem[3][2] = -3.30350175000e+00;
        slipSystem[3][3] =  burgRef[1][0];
        slipSystem[3][4] =  burgRef[1][1];
        slipSystem[3][5] =  burgRef[1][2];
        slipSystem[3][6] =  burgRef[2][0];
        slipSystem[3][7] =  burgRef[2][1];
        slipSystem[3][8] =  burgRef[2][2];

        slipSystem[4][0] = -3.30350175000e+00;
        slipSystem[4][1] =  3.29473250000e-01;
        slipSystem[4][2] = -3.30350175000e+00;
        slipSystem[4][3] =  burgRef[1][0];
        slipSystem[4][4] =  burgRef[1][1];
        slipSystem[4][5] =  burgRef[1][2];
        slipSystem[4][6] =  burgRef[3][0];
        slipSystem[4][7] =  burgRef[3][1];
        slipSystem[4][8] =  burgRef[3][2];

        slipSystem[5][0] = -3.30350175000e+00;
        slipSystem[5][1] = -3.30350175000e+00;
        slipSystem[5][2] =  3.29473250000e-01;
        slipSystem[5][3] =  burgRef[2][0];
        slipSystem[5][4] =  burgRef[2][1];
        slipSystem[5][5] =  burgRef[2][2];
        slipSystem[5][6] =  burgRef[3][0];
        slipSystem[5][7] =  burgRef[3][1];
        slipSystem[5][8] =  burgRef[3][2];

/*
 *      Only need to normalize first 3 components; other components are
 *      from <burgRef> which is already normalized.
 */
        for (i = 0; i < 6; i++) {
            Normalize(&slipSystem[i][0], &slipSystem[i][1], &slipSystem[i][2]);
        }

/*
 *         Specify fitted mobility parameters
 *
 *         The correspondence between the parameter names here and those in 
 *         Rich Becker's Ta-2 paper:
 *
 *             alpha = alpha_0/pielrs
 *             beta = alpha_t
 *             delta = beta_0
 *             epsilon = T/beta_t
 *             c0 = c0
 *             gamma = nu_0 = 1.e-9
 *
 *             phi = kai_0/pielrs
 *             eta = kai_0*kai_1/pielrs = phi*kai_1
 *
 *             ftherm = tao_therm (it's tao_stress/tao_p)
 *
 *         so fdrag and ftherm are in units of pierls; however, fmag is in 
 *         unit of Pascal*burgnormal. 
 *
 *         Note: in the fitting for R-V mobility, it's the actual stress values
 *         are used, not the reduced value. therefore the alpha value needs to
 *         be scaled with their corresponding pierls stress
 */

/*
 *      Screw parameters
 */
        GammaScrew[0]=1.e-9;
        GammaScrew[1]=1.e-9;
        GammaScrew[2]=1.e-9;

        PeierlsScrew[0] = 8.21e8;
        PeierlsScrew[1] = 8.28e8;
        PeierlsScrew[2] = 8.28e8;

/* 
 *      For long screws (always on (110) planes only
 */
        AlphaScrew[0]=11475.0*1.e6/PeierlsScrew[0];
        BetaScrew[0]=331.0;
        DeltaScrew[0]=0.49;
        EpsilonScrew[0]=temperature/1395.3;
        PhiScrew[0]=(6.933e4+163.0*temperature)*1.e6/PeierlsScrew[0];
        EtaScrew[0]=PhiScrew[0]*(-3.4e-3+1.e-5*temperature);

/* 
 *      For short screws on (1,1,0) planes
 */
        AlphaScrew[1]=1913.9*1.e6/PeierlsScrew[1];
        BetaScrew[1]=179.4;
        DeltaScrew[1]=0.24;
        EpsilonScrew[1]=temperature/769.2307;
        PhiScrew[1]=(1.5852e5+70.0*temperature)*1.e6/PeierlsScrew[1];
        EtaScrew[1]=PhiScrew[1]*(0.0205-8.3e-6*temperature);

/* 
 *      For short screws on (1, 10,10) planes 
 */
        AlphaScrew[2]=7772.0*1.e6/PeierlsScrew[2];
        BetaScrew[2]=-1837.0;
        DeltaScrew[2]=0.4;
        EpsilonScrew[2]=0.0;
        PhiScrew[2]=(8931.0+210.4*temperature)*1.e6/PeierlsScrew[2];
        EtaScrew[2]=PhiScrew[2]*(-1.19e-2 + 2.87e-5*temperature);

/* 
 *      Non screw parameters 
 */
        GammaEdge[0]=1.e-9;
        GammaEdge[1]=1.e-9;
        GammaEdge[2]=1.e-9;

        PeierlsEdge[0] = 3.02e8;
        PeierlsEdge[1] = 3.37e8;
        PeierlsEdge[2] = 3.37e8;

/* 
 *      For long non-screws (always on (110) planes only 
 */
        AlphaEdge[0]=13.2*1.e6/PeierlsEdge[0];
        BetaEdge[0]=123.2;
        DeltaEdge[0]=0.1;
        EpsilonEdge[0]=temperature/6172.5;
        PhiEdge[0]=(2.3072e4)*1.e6/PeierlsEdge[0];
        EtaEdge[0]=PhiEdge[0]*(0.0118);

/* 
 *      For short non-screws on (1,1,0) planes 
 */
        AlphaEdge[1]=2997.8*1.e6/PeierlsEdge[1];
        BetaEdge[1]=6791.7;
        DeltaEdge[1]=1.0;
        EpsilonEdge[1]=0.0;
        PhiEdge[1]=(-3526.0+56.1*temperature)*1.e6/PeierlsEdge[1];
        EtaEdge[1]=PhiEdge[1]*(0.1573+5.53e-5*temperature);

/*
 *      For short non-screws on (1, 10,10) planes 
 */
        AlphaEdge[2]=1640.3*1.e6/PeierlsEdge[2];
        BetaEdge[2]=463.5;
        DeltaEdge[2]=0.75;
        EpsilonEdge[2]=temperature/1200.0;
        PhiEdge[2]=(-11695+82.0*temperature)*1.e6/PeierlsEdge[2];
        EtaEdge[2]=PhiEdge[2]*(0.1203+1.28e-4*temperature);

/*
 *      For junction glide, set it to fastest mobilit slope
 *
 *      previous value: BJglide = 5.979902204668377e-01;
 */
        BJglide = PhiEdge[1] / 9.0 * PeierlsEdge[1] / C0;

/*
 *      Calculate the width of the trenches based on Vasily Bulatov's
 *      estimate of kink widths
 */
#if 0
        a = 0.7 * sqrt(PeierlsScrew[0] / MU);
        b = 0.7 * sqrt(PeierlsEdge[0] / MU);
        c = 0.7 * sqrt(PeierlsScrew[1] / MU);
        d = 0.7 * sqrt(PeierlsEdge[1] / MU);
        e = 0.7 * sqrt(PeierlsScrew[2] / MU);
        f = 0.7 * sqrt(PeierlsEdge[2] / MU);
#endif

/*
 *      Use 5 degrees on each side of trench rather than about 2
 */
        a = b = c = d = e = f = 0.0875;

        thetaCrit[0][0] =  atan(a);  thetaCrit[0][1]  = atan(b);
        thetaCrit[1][0] =  atan(c);  thetaCrit[1][1]  = atan(d);
        thetaCrit[2][0] =  atan(a);  thetaCrit[2][1]  = atan(b);
        thetaCrit[3][0] =  atan(c);  thetaCrit[3][1]  = atan(d);
        thetaCrit[4][0] =  atan(a);  thetaCrit[4][1]  = atan(b);
        thetaCrit[5][0] =  atan(c);  thetaCrit[5][1]  = atan(d);
        thetaCrit[6][0] =  atan(e);  thetaCrit[6][1]  = atan(f);
        thetaCrit[7][0] =  atan(e);  thetaCrit[7][1]  = atan(f);
        thetaCrit[8][0] =  atan(e);  thetaCrit[8][1]  = atan(f);
        thetaCrit[9][0] =  atan(e);  thetaCrit[9][1]  = atan(f);
        thetaCrit[10][0] = atan(e);  thetaCrit[10][1] = atan(f);
        thetaCrit[11][0] = atan(e);  thetaCrit[11][1] = atan(f);


        return;
}


static void JunctionDrag(real8 BJglide, real8 vel[3], real8 t[3],
                         real8 fjdrag[3], real8 dfjdragdv[3][3])
{
        real8 BJline;
        real8 tmp[3][3];

        BJline = BJglide;

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

static void FscaleFast(Param_t *param, real8 ftest, real8 vin, real8 burg[3],
                       int mobIndex, real8 *fout, real8 *dfdv, real8 *d2fdv2)
{
        real8  burgnormal,vmag,vsign,peierls;
        real8  phi; 
        real8  relativity, ratio, ratioMax, finter;

        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);	

#if 0
/*
 *      The factor of 6 choice below is based on short (110) 75 Deg 
 *      The slope for the thermal part is around 1, (the smaller, the larger
 *      the mobility) the slope for the drag part if around 4
 *      the choice of the slope for the fast is around 0.5 = 4/8
 */ 
        phi = PhiEdge[1]/6.0;
        peierls = PeierlsEdge[1] * burgnormal;
#endif

/*
 *      MRhee:
 *
 *      For now, the fastest one is "approximately" two times faster
 *      than the second fastest one.
 *
 *      With T=600, the slope becomes larger than the non-linear
 *      kinetic ones.  It was agreed that we use a constant value
 *      consistent with T=300K.
 */

        /* PhiEdge[1]=(-3526.0+56.1*temperature)*1.e6/PeierlsEdge[1]; */

        phi = (-3526.0+56.1 * 300.0 ) * 1.e6 / PeierlsEdge[1] / 9.0 ;
        peierls = PeierlsEdge[1] ;

        vmag = fabs(vin);
        vsign = (real8)Sign(vin);
        ratio = vmag / C0;
        ratioMax = 0.90;

/*
 *      This magic number is from Graph on page 1
 *      We want to slow down the fast ones for the long 111 burgers vector
 *      3.e6 Pa.sec is the slope on p.1
 */
        if (burgnormal > 0.95) {
/*
 *          These values give exact values of 2.e6 Pa.sec drag
 */
            if (ratio < ratioMax) {
                relativity= 1.0 / sqrt(1.0 - (ratio * ratio));
                *fout = 3.e6  * param->burgMag * vmag * vsign * relativity;
                *dfdv = 3.e6  * param->burgMag * relativity * relativity * relativity;
                *d2fdv2 = 0.0;  
            } else {
                relativity= 1.0 / sqrt(1.0 - (ratioMax * ratioMax));
                finter= 3.e6  * param->burgMag * ratioMax * C0 * vsign * relativity;
                *dfdv=3.e6  * param->burgMag * relativity * relativity * relativity;
                *fout=finter+ *dfdv * C0 * ( ratio - ratioMax ) * vsign;
                *d2fdv2= 0.0e0;
            }

        } else {
            /* PeierlsEdge[1] = 3.37e8; */
            /* PhiEdge[1]=(-3526.0+56.1*temperature)*1.e6/PeierlsEdge[1]; */
            phi = (-3526.0+56.1 * 300.0 ) * 1.e6 / PeierlsEdge[1] / 9.0 ;
            if (ratio < ratioMax) {
                relativity= 1.0 / sqrt(1.0 - (ratio * ratio));
                *fout = peierls * phi * ratio * vsign * relativity;
/*
 *              This value is 4.48e5 * burgmag Pa.s
 */
                *dfdv = peierls * phi / C0 * relativity * relativity * relativity;
                *d2fdv2 = 0.0; 
            } else {
                relativity= 1.0 / sqrt(1.0 - (ratioMax * ratioMax));
                finter= peierls * phi * ratioMax * vsign * relativity;
                *dfdv = peierls * phi / C0 * relativity * relativity * relativity;
                *fout=finter+ *dfdv * C0 * ( ratio - ratioMax ) * vsign;
                *d2fdv2= 0.0;
            }
        }

#ifdef NAN_CHECK
        if (isnan(*fout) || isinf(*fout)) {
            printf("FscaleFast: fout = %lf\n", *fout);
        }
        if (isnan(*dfdv) || isinf(*dfdv)) {
            printf("FscaleFast: dfdv = %lf\n", *dfdv);
        }
#endif

        return;
}


static void FscaleScrew(Param_t *param, real8 fsinit, real8 vin, real8 burg[3],
                        int mobIndex, real8 *fout, real8 *dfdv, real8 *d2fdv2)
{
        real8  alpha, beta, gamma, delta, phi, epsilon, eta, n, ninv;
        real8  peierls, burgnormal, temperature, vmag, vsign;
        real8  dfdragdv, dfthermdv, dfmagdfdrag, dfmagdftherm, dfmagdv;
        real8  relativity, drelativitydratio;
        real8  fmag, fdrag, ftherm, finter, power, test;
        real8  ratio, ratioMax;


        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);

        temperature = param->TempK;

        peierls     = PeierlsScrew[mobIndex] * burgnormal;
        alpha       = AlphaScrew[mobIndex];
        beta        = BetaScrew[mobIndex];
        gamma       = GammaScrew[mobIndex];
        delta       = DeltaScrew[mobIndex];
        phi         = PhiScrew[mobIndex];
        epsilon     = EpsilonScrew[mobIndex];
        eta         = EtaScrew[mobIndex];
        n           = 15.0;

        vmag = fabs(vin);
        vsign = (real8)Sign(vin);
        ninv = 1.0 / n;

        ratio = vmag / C0;
        ratioMax = 0.90;
        power = delta + epsilon;

        ftherm = alpha * exp(temperature / beta) *
                 (pow(ratio+gamma, power) - pow(gamma, power));
        dfthermdv = alpha * power / C0 * exp(temperature / beta) *
                    pow(ratio+gamma, power-1.0);

        test = phi * ratio - eta;

        if (vmag < 1.0) {
            *fout = 0.0e0;
            *dfdv = peierls * alpha * power / C0 * exp(temperature/beta) *
                    pow(gamma, power-1.0);
            *d2fdv2 = 0.0;
#ifdef NAN_CHECK
            if (isnan(*dfdv) || isinf(*dfdv)) {
                printf("FscaleScrew: vmag < 1.0, dfdv = %lf\n", *dfdv);
            }
#endif
            return;
        } else if (test < 0.0) {
            fdrag = 0.0;
            dfdragdv = 0.0;
        } else if (ratio > ratioMax) {
            relativity = 1.0 / sqrt(1.0 - (ratioMax * ratioMax));
            drelativitydratio = ratioMax * relativity * relativity * relativity;
            finter = (phi * ratioMax - eta) * relativity;
            dfdragdv = 1.0 / C0 * (phi * relativity +
                                   (phi * ratioMax -eta) * drelativitydratio);
            fdrag = finter + dfdragdv * C0 * (ratio - ratioMax);
        } else {
            relativity = 1.0 / sqrt(1.0 - (ratio * ratio));
            drelativitydratio = ratio * relativity * relativity * relativity;
            fdrag = (phi * ratio - eta) * relativity;
            dfdragdv = 1.0 / C0 *
                       (phi * relativity +
                        (phi * ratio - eta) * drelativitydratio);
        }

        fmag = peierls * pow(pow(fdrag, n)+pow(ftherm, n), ninv);
        dfmagdfdrag = peierls * pow(fdrag, n-1.0) *
                      pow(pow(fdrag, n)+pow(ftherm, n), ninv-1.0);
        dfmagdftherm = peierls * pow(ftherm, n-1.0) *
                       pow(pow(fdrag, n)+pow(ftherm, n), ninv-1.0);
        dfmagdv = dfmagdfdrag * dfdragdv + dfmagdftherm * dfthermdv;

        *fout = fmag * vsign;
        *dfdv = dfmagdv;
        *d2fdv2 = 0.0;

#ifdef NAN_CHECK
        if (isnan(*fout) || isinf(*fout)) {
            printf("FscaleScrew: fout = %lf\n", *fout);
        }
        if (isnan(*dfdv) || isinf(*dfdv)) {
            printf("FscaleScrew: dfdv = %lf\n", *dfdv);
        }
#endif
        return;
}


static void FastDrag(Param_t *param, real8 vel[3], real8 burg[3],
                     real8 climbDir[3],  int slipSystemIndex, int mobIndex,
                     real8 lineDir[3], real8 fastDrag[3],
                     real8 dFastDragdv[3][3])
{
        int   m, n;
        real8 ftest, fout, vMag;
        real8 dfdv, d2fdv2;
        real8 glideDir[3], finit[3];
        real8 climbDirTM[3][3], glideDirTM[3][3], lineDirTM[3][3];

/*
 *      Save a copy of the initial value of fastDrag for later
 */
        for (m = 0; m < 3; m++) {
            finit[m] = fastDrag[m];
        }

        cross(climbDir, lineDir, glideDir);
        vMag = DotProduct(vel, glideDir);
        ftest = DotProduct(finit, glideDir);

        FscaleFast(param, ftest, vMag, burg, mobIndex, &fout, &dfdv, &d2fdv2);


        for (m = 0; m < 3; m++) {
            fastDrag[m] = fout * glideDir[m];
        }

        Vec3TransposeAndMult(glideDir, glideDirTM);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dFastDragdv[m][n] = dfdv * glideDirTM[m][n];
            }
        }

        vMag = DotProduct(vel, lineDir);
        ftest = DotProduct(finit, lineDir);

        FscaleFast(param, ftest, vMag, burg, mobIndex, &fout, &dfdv, &d2fdv2);

        for (m = 0; m < 3; m++) {
            fastDrag[m] += fout * lineDir[m];
        }

        Vec3TransposeAndMult(lineDir, lineDirTM);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dFastDragdv[m][n] += dfdv * lineDirTM[m][n];
            }
        }

        vMag = DotProduct(vel, climbDir);
        ftest = DotProduct(finit, climbDir);

        FscaleScrew(param, ftest, vMag, burg, mobIndex, &fout, &dfdv, &d2fdv2);

        Vec3TransposeAndMult(climbDir, climbDirTM);

        for (m = 0; m < 3; m++) {
            fastDrag[m] += fout * climbDir[m];
            for (n = 0; n < 3; n++) {
                dFastDragdv[m][n] += dfdv * climbDirTM[m][n];
            }
        }

        return;
}


static void FscaleEdge(Param_t *param, real8 feinit, real8 vin, real8 burg[3],
                       int mobIndex, real8 *fout, real8 *dfdv, real8 *d2fdv2)
{
        real8  alpha, beta, gamma, delta, phi, epsilon, eta, n, ninv;
        real8  peierls, burgnormal, temperature, vmag, vsign;
        real8  dfdragdv, dfthermdv, dfmagdfdrag, dfmagdftherm, dfmagdv;
        real8  relativity, drelativitydratio;
        real8  fmag, fdrag, ftherm, finter, power, test;
        real8  ratio, ratioMax;


        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);

        temperature = param->TempK;

        peierls     = PeierlsEdge[mobIndex] * burgnormal; 
        alpha       = AlphaEdge[mobIndex];
        beta        = BetaEdge[mobIndex];
        gamma       = GammaEdge[mobIndex];
        delta       = DeltaEdge[mobIndex];
        phi         = PhiEdge[mobIndex];
        epsilon     = EpsilonEdge[mobIndex];
        eta         = EtaEdge[mobIndex];
        n           = 15.0;

        vmag = fabs(vin);
        vsign = (real8)Sign(vin);
        ninv = 1.0 / n;

        ratio = vmag / C0;
        ratioMax = 0.90;
        power = delta + epsilon;

        ftherm = alpha * exp(temperature / beta) *
                 (pow(ratio+gamma, power) - pow(gamma, power));
        dfthermdv = alpha * power / C0 * exp(temperature / beta) *
                    pow(ratio+gamma, power-1.0);

        test = phi * ratio - eta;

        if (vmag < 1.0) {
            *fout = 0.0e0;
            *dfdv = peierls * alpha * power / C0 * exp(temperature/beta) *
                    pow(gamma, power-1.0);
            *d2fdv2 = 0.0;
#ifdef NAN_CHECK
            if (isnan(*dfdv) || isinf(*dfdv)) {
                printf("FscaleEdge: vmag < 1.0, dfdv = %lf\n", *dfdv);
            }
#endif
            return;
        } else if (test < 0.0) {
            fdrag = 0.0;
            dfdragdv = 0.0;
        } else if (ratio > ratioMax) {
            relativity = 1.0 / sqrt(1.0 - (ratioMax * ratioMax));
            drelativitydratio = ratioMax * relativity * relativity * relativity;
            finter = (phi * ratioMax - eta) * relativity;
            dfdragdv = 1.0 / C0 * (phi * relativity +
                                   (phi * ratioMax -eta) * drelativitydratio);
            fdrag = finter + dfdragdv * C0 * (ratio - ratioMax);
        } else {
            relativity = 1.0 / sqrt(1.0 - (ratio * ratio));
            drelativitydratio = ratio * relativity * relativity * relativity;
            fdrag = (phi * ratio - eta) * relativity;
            dfdragdv = 1.0 / C0 *
                       (phi * relativity +
                        (phi * ratio - eta) * drelativitydratio);
        }

        fmag = peierls * pow(pow(fdrag, n)+pow(ftherm, n), ninv);
        dfmagdfdrag = peierls * pow(fdrag, n-1.0) *
                      pow(pow(fdrag, n)+pow(ftherm, n), ninv-1.0);
        dfmagdftherm = peierls * pow(ftherm, n-1.0) *
                       pow(pow(fdrag, n)+pow(ftherm, n), ninv-1.0);
        dfmagdv = dfmagdfdrag * dfdragdv + dfmagdftherm * dfthermdv;

        *fout = fmag * vsign;
        *dfdv = dfmagdv;
        *d2fdv2 = 0.0;

#ifdef NAN_CHECK
        if (isnan(*fout) || isinf(*fout)) {
            printf("FscaleEdge: fout = %lf\n", *fout);
        }
        if (isnan(*dfdv) || isinf(*dfdv)) {
            printf("FscaleEdge: dfdv = %lf\n", *dfdv);
        }
#endif
        return;
}


static void EdgeDrag(Param_t *param, real8 vel[3], real8 burg[3],
                     real8 climbDir[3],  int slipSystemIndex, int mobIndex,
                     real8 fedrag[3], real8 dfedragdv[3][3], int isLongBurg)
{
        int   m, n;
        int   tmpInt, planeID;
        real8 vMag, fout, ftest;
        real8 dfdv, d2fdv2;
        real8 finit[3];
        real8 lineDir[3], glideDir[3];
        real8 lineDirTM[3][3], glideDirTM[3][3], climbDirTM[3][3];

/*
 *      Save a copy of the initial value of fedrag for later
 */
        for (m = 0; m < 3; m++) {
            finit[m] = fedrag[m];
        }

        planeID = slipSystemIndex / 2;
        tmpInt = slipSystemIndex % 2;

        if (tmpInt == 0) {
            lineDir[0] = slipSystem[planeID][6];
            lineDir[1] = slipSystem[planeID][7];
            lineDir[2] = slipSystem[planeID][8];
        } else {
            lineDir[0] = slipSystem[planeID][3];
            lineDir[1] = slipSystem[planeID][4];
            lineDir[2] = slipSystem[planeID][5];
        }

        cross(climbDir, lineDir, glideDir);
        
        vMag = DotProduct(vel, lineDir);
        ftest = DotProduct(finit, lineDir);

        FscaleFast(param, ftest, vMag, burg, mobIndex, &fout, &dfdv, &d2fdv2);

        for (m = 0; m < 3; m++) {
            fedrag[m] = fout * lineDir[m];
        }

        Vec3TransposeAndMult(lineDir,  lineDirTM);
        
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] = dfdv  * lineDirTM[m][n]; 
            }
        }

        vMag = DotProduct(vel, glideDir);
        ftest = DotProduct(fedrag, glideDir);

        FscaleEdge(param, ftest, vMag, burg, mobIndex, &fout, &dfdv, &d2fdv2);

        for (m = 0; m < 3; m++) {
            fedrag[m] +=  fout * glideDir[m];
        }

        Vec3TransposeAndMult(glideDir, glideDirTM);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] += dfdv * glideDirTM[m][n];
            }
        }

        vMag = DotProduct(vel, climbDir);
        ftest = DotProduct(finit, climbDir);

/*
 *      Note: for now, the 2nd derivative (d2fdv2) from FscaleEdge()
 *            is zero on return from the function call!
 */

        FscaleScrew(param, ftest, vMag, burg, mobIndex, &fout, &dfdv, &d2fdv2);

        Vec3TransposeAndMult(climbDir, climbDirTM);

        for (m = 0; m < 3; m++) {
            fedrag[m] += fout * climbDir[m];
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] += dfdv * climbDirTM[m][n];
            }
        }

        return;
}



/*
 *    finit:  Initial guess... updated before return to caller.
 */
static void ScrewDrag(Param_t *param, real8 fsdrag[3], real8 vel[3],
                      real8 burg[3], real8 climbDir[3], int slipSysID,
                      int mobIndex, real8 dfsdragdv[3][3], int isLongBurg)
{
        int   m, n;
        real8 burgNorm;
        real8 vMag, ftest;
        real8 fout1, dfdv1, d2fdv21;
        real8 fout2, dfdv2, d2fdv22;
        real8 lineDir[3], glideDir[3], finit[3];
        real8 lineDirTM[3][3], glideDirTM[3][3], climbDirTM[3][3];

/*
 *      Save a copy of initial fsdrag for later
 */
        for (m = 0; m < 3; m++) {
            finit[m] = fsdrag[m];
        }

        burgNorm = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);

        lineDir[0] = burg[0] / burgNorm;
        lineDir[1] = burg[1] / burgNorm;
        lineDir[2] = burg[2] / burgNorm;

        cross(climbDir, lineDir, glideDir);

        Vec3TransposeAndMult(lineDir,  lineDirTM);
        Vec3TransposeAndMult(glideDir, glideDirTM);
        Vec3TransposeAndMult(climbDir, climbDirTM);
        
/*
 *      Note: for now, the 2nd derivative (d2fdv21) from FscaleScrew()
 *            is zero on return from the function call!
 */
        vMag = DotProduct(vel, lineDir);
        ftest = DotProduct(finit, lineDir);

        FscaleFast(param, ftest, vMag, burg, mobIndex, &fout1, &dfdv1, &d2fdv21);

        for (m = 0; m < 3; m++) {
            fsdrag[m] = fout1 * lineDir[m];
        }

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfsdragdv[m][n] = dfdv1   * lineDirTM[m][n]; 
            }
        }

        vMag = DotProduct(vel, glideDir);
        ftest = DotProduct(finit, glideDir);

        FscaleScrew(param, ftest, vMag, burg, mobIndex, &fout1, &dfdv1,
                    &d2fdv21);

        for (m = 0; m < 3; m++) {
            fsdrag[m] += fout1 * glideDir[m];
        }

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfsdragdv[m][n] += dfdv1   * glideDirTM[m][n]; 
            }
        }

        vMag = DotProduct(vel, climbDir);
        ftest = DotProduct(finit, climbDir);
/*
 *      Note: for now, the 2nd derivative (d2fdv22) from FscaleScrew()
 *            is zero on return from the function call!
 */

        FscaleScrew(param, ftest, vMag, burg, mobIndex, &fout2, &dfdv2,
                    &d2fdv22);

        for (m = 0; m < 3; m++) {
            fsdrag[m] += fout2 * climbDir[m];
        }

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfsdragdv[m][n] += dfdv2 * climbDirTM[m][n];
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     MobilityLaw_Rhombohedral_Va_nl_planar
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node                Pointer to the node for which to
 *                                   calculate velocity
 *          OUT: invdferrordv        Array in which to return the inverse of
 *                                   the drag matrix used in the function.
 *          OUT: glideConstraints    Array in which to return the normals
 *                                   to the glide planes to which this node
 *                                   is restricted.
 *          OUT: numGlideConstraints Number of glide plane normals returned
 *
 *      Returns:  0 on success
 *                1 if velocity could not be determined
 *
 *************************************************************************/
int Mobility_Rhombohedral_Va_nl_planar(Home_t *home, Node_t *node,
                                       MobArgs_t *mobArgs)
{
        int     i, j, m, n, iAll, numNbrs;
        int     iterCnt, maxIter, maxIter2, count;
        int     zeroLenSeg[MAX_NBRS];
        int     notConverged;
        int     numNonZeroLenSegs = 0;
        int     numConstraints;
        int     index1, index2, index3;
        int     mobIndex;
        int     planeID;
        real8   velerrortol, velerror, minvelerror;
        real8   eps, avgFactor, correctionAdjustment;
        real8   vtest[3], fn[3], rt[3], zero3[3], ferror[3];
        real8   dferrordv[3][3];
        real8   dfjdragdv[3][3];
        real8   dferrordvold[3][3];
        real8   ndir[MAX_NBRS][3], burg[MAX_NBRS][3];
        real8   fedge[MAX_NEDGE][3], ffast[MAX_NBRS][3], fscrew[MAX_NSCREW][3];
        real8   dfedv[3][3], dffdv[3][3], dfsdv[3][3];
        real8   correction[3];
        real8   vtmp  [3] = { 0.0, 0.0, 0.0 };
        real8   fjdrag[3] = { 0.0, 0.0, 0.0 };
        real8   mu, lmax;
        real8   forceerror, ferrortol;
        real8   totLen, massMult=0.0, massMatrix[3][3], vOld[3], vDiff[3], tmp3B[3];
        Param_t *param;
        Node_t  *nbrNode;
        real8   S1, S2;
        real8   M2[2];
        real8   maxVal;
        real8   burgRefXburg[13];
        real8   slipXndir[6];
        real8   dir2[3];
        real8   theta, thetaTest, sinThetaTest, cosThetaTest;
        real8   dirCrit[3];
        real8   L[MAX_NBRS][9];  /* See comment below */
        real8   nonSpecialMinLength;
        real8   minLimit[3], maxLimit[3];
        int     minLimitSet[3], maxLimitSet[3];
        int     offset;
        real8   shortSegCutoff, shortSegRatio = 1.0;
        int     hasShortSeg = 0;

/*
 *      For each arm, we need a number of data items.  These will
 *      be stored in L[i][j] where 'i' corresponds to the arm index.
 *      (The array name is useless, but matches the array used in
 *      the matlab code, which makes for easier debugging/comparisons).
 *      Components of L are:
 *
 *          L[i][0] = burgers vector index of the dislocation as ordered
 *                    in burgref
 *          L[i][1] = slip system of the dislocation as ordered in slipsys
 *                    rows increase id by 2, burg increase id by 1
 *          L[i][2] = total line length of the dislocation
 *          L[i][3] = if segment is a glide dislocation, this is the
 *                    length of the screw dislocation if it is in a
 *                    mobility trench
 *          L[i][4] = if segment is a glide dislocation this is the
 *                    length of the non-screw dislocation if it is in a
 *                    mobility trench
 *          L[i][5] = if segment is a glide dislocation this is the
 *                    length of the non-special dislocation character
 *          L[i][6-8] = if segment is not a junction, this is the direction
 *                    of the non-special dislocation length.  If segment
 *                    is a junctions, this is the dislocation line direction
 *
 *      We'll also define some constants that can be used as indices for
 *      the second dimension of the "L" array.  Makes references a little
 *      clearer
 */

#define BURG_INDEX      0  /* L[i][1] in matlab */
#define SLIP_SYS_ID     1  /* L[i][2] in matlab */
#define TOT_LEN         2  /* L[i][3] in matlab */
#define SCREW_LEN       3  /* L[i][4] in matlab */
#define NON_SCREW_LEN   4  /* L[i][5] in matlab */
#define NON_SPECIAL_LEN 5  /* L[i][6] in matlab */
#define LINE_DIR        6  /* denotes L[i][7:9] in matlab */


        real8 (*invdferrordv)[3]     =  mobArgs->invDragMatrix;
        real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdferrordv);
       *numGlideConstraints = 0;


        eps = 1.0e-12;
        iterCnt = 0;
        maxIter = 100;
        maxIter2 = 200;

        param = home->param;
        mu = param->shearModulus;
        lmax = param->maxSeg;
        shortSegCutoff = 0.50 * param->minSeg;

        memset(zero3, 0, sizeof(zero3));


/*
 *      If we need to include inertial terms, need some more initializations
 */
        if (param->includeInertia != 0) {
            massMult = 0.25 * param->massDensity * (param->burgMag * param->burgMag);
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
            Fatal("Mobility_Rhombohedral_Va_nl_planar: Node segment count (%d) "
                  "exceeds\nmaximum allowed (%d) in static arrays.  Increase "
                  "MAX_NBRS\nand recompile\n", node->numNbrs, MAX_NBRS);
        }

/*
 *      Initialize some arrays for later...
 */
        vtest[0] = node->vX;
        vtest[1] = node->vY;
        vtest[2] = node->vZ;

        VECTOR_COPY(vOld, vtest);

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
            if (DotProduct(rt, rt) < eps) {
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

            MatrixMult((real8 *)burgRef, 13, 3, 3,
                       burg[i], 1, 1,
                       burgRefXburg, 1);

            FindAbsMax(burgRefXburg, 13, &maxVal, &index1);
            
            L[i][0]    = index1;
            L[i][2]    = Normal(rt);

/*
 *          If the node has a very short segment, we need to flag it for later
 */
            if (L[i][2] < shortSegCutoff) {
                hasShortSeg = 1;
                shortSegRatio = MIN(shortSegRatio, L[i][2]/param->minSeg);
            }

/*
 *          Initialize dislocation to be outside of vicinal ranges
 */
            L[i][3] = 0.0;
            L[i][4] = 0.0;
            L[i][5] = L[i][2];
            L[i][6] = rt[0] / L[i][2];
            L[i][7] = rt[1] / L[i][2];
            L[i][8] = rt[2] / L[i][2];

/*
 *          For the non-special length, set a minimum length under
 *          which we treat the length as zero.
 */
            nonSpecialMinLength = 1.0e-3 * L[i][2];


            if (L[i][0] < 4) {
                real8 rtDotSlip, rtDotDir, slipDotDir, slipDotDirSq;
/*
 *              Dislocation is a glide dislocation
 */
                ndir[i][0] = node->nx[iAll];
                ndir[i][1] = node->ny[iAll];
                ndir[i][2] = node->nz[iAll];

/*
 *              If necessary, rotate the glide plane normal from the
 *              laboratory frame to the crystal frame
 */
                if (param->useLabFrame) {
                    real8 ndirRot[3];
                    Matrix33Vector3Multiply(param->rotMatrixInverse,
                                            ndir[i], ndirRot);
                    VECTOR_COPY(ndir[i], ndirRot);
                }

/*
 *              Identify the slip system
 */
                MatrixMult((real8 *)slipSystem, 6, 3, 9,
                           ndir[i], 1, 1,
                           slipXndir, 1);

                FindAbsMax(slipXndir, 6, &maxVal, &index1);

                S1 = fabs(DotProduct((&slipSystem[index1][3]),
                                     burgRef[(int)L[i][0]]));
                S2 = fabs(DotProduct((&slipSystem[index1][6]),
                                     burgRef[(int)L[i][0]]));

                index2 = (S1 > S2 ? 0 : 1);
               
                L[i][1] = 2 * index1 + index2;

/*
 *              Now check if dislocation is within a trench
 */
                M2[0] = fabs(DotProduct((&slipSystem[index1][3]),
                                        (&L[i][6])));
                M2[1] = fabs(DotProduct((&slipSystem[index1][6]),
                                        (&L[i][6])));

                FindAbsMax(M2, 2, &maxVal, &index3);

                theta = acos(maxVal);

                if (index2 == index3) {
                    offset = 3;
                    thetaTest = thetaCrit[(int)L[i][1]][0];
                } else {
                    offset = 4;
                    thetaTest = thetaCrit[(int)L[i][1]][1];
                }

                if (theta < thetaTest) {
/*
 *                  Dislocation might be in a screw trench
 */
                    if (index3 == 0) {

                        cross(&slipSystem[index1][3], slipSystem[index1], dir2);
                        S1 = DotProduct(dir2, (&L[i][6]));
                        dir2[0] *= Sign(S1);
                        dir2[1] *= Sign(S1);
                        dir2[2] *= Sign(S1);

                        S1 = DotProduct((&slipSystem[index1][3]),
                                        (&L[i][6]));

                        cosThetaTest = cos(thetaTest);
                        sinThetaTest = sin(thetaTest);

                        for (j = 0; j < 3; j++) {
                            dirCrit[j] = (cosThetaTest * Sign(S1) *
                                          slipSystem[index1][j+3]) +
                                         (sinThetaTest * dir2[j]);
                            L[i][6+j] = dirCrit[j];
                        }

                        rtDotSlip = DotProduct(rt, (&slipSystem[index1][3]));
                        rtDotDir  = DotProduct(rt, dirCrit);
                        slipDotDir = DotProduct((&slipSystem[index1][3]),
                                                 dirCrit);
                        slipDotDirSq = slipDotDir * slipDotDir;

                        L[i][offset] = (rtDotSlip - (rtDotDir * slipDotDir)) /
                                  (1.0 - slipDotDirSq);
                        L[i][offset] = fabs(L[i][offset]);
                        L[i][5] = ((-rtDotSlip * slipDotDir) + rtDotDir) /
                                  (1.0 - slipDotDirSq);
                        
/*
 *                      The non-special length should be positive.  Just add a
 *                      check to make sure.  Abort or set to zero?
 */
                        if (L[i][5] < nonSpecialMinLength) {
                            L[i][5] = 0.0;
                        }
                    } else { /* index3 != 0 */
/*
 *                      Dislocation might be in a non-screw trench
 */
                        cross(&slipSystem[index1][6], slipSystem[index1], dir2);
                        S1 = DotProduct(dir2, (&L[i][6]));
                        dir2[0] *= Sign(S1);
                        dir2[1] *= Sign(S1);
                        dir2[2] *= Sign(S1);

                        S1 = DotProduct((&slipSystem[index1][6]),
                                        (&L[i][6]));

                        cosThetaTest = cos(thetaTest);
                        sinThetaTest = sin(thetaTest);

                        for (j = 0; j < 3; j++) {
                            dirCrit[j] = (cosThetaTest * Sign(S1) *
                                          slipSystem[index1][j+6]) +
                                         (sinThetaTest * dir2[j]);
                            L[i][6+j] = dirCrit[j];
                        }

                        rtDotSlip = DotProduct(rt, (&slipSystem[index1][6]));
                        rtDotDir  = DotProduct(rt, dirCrit);
                        slipDotDir = DotProduct((&slipSystem[index1][6]),
                                                 dirCrit);
                        slipDotDirSq = slipDotDir * slipDotDir;

                        L[i][offset] = (rtDotSlip - (rtDotDir * slipDotDir)) /
                                  (1.0 - slipDotDirSq);
                        L[i][offset] = fabs(L[i][offset]);
                        L[i][5] = ((-rtDotSlip * slipDotDir) + rtDotDir) /
                                  (1.0 - slipDotDirSq);

/*
 *                      The non-special length should be positive.  Just add a
 *                      check to make sure.  Abort or set to zero?
 */
                        if (L[i][5] < nonSpecialMinLength) {
                            L[i][5] = 0.0;
                        }
                    }
                }
            }  /* end if (L[i][0] < 4) */

            totLen += L[i][2];

        }  /* end loop over all neighbors */

/*
 *      Create a glide constraint matrix
 */
        numConstraints = 0;

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                glideConstraints[m][n] = 0.0;
            }
        }


        for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {

/*
 *          Skip any zero length segments...
 */
            if (zeroLenSeg[iAll]) {
                i--;
                continue;
            }

            if ((L[i][0] < 4) && (numConstraints != 3)) {

                planeID = (int)L[i][1] / 2;

/*
 *              numConstraints:
 *                  1 allows movement in plane
 *                  2 allows movement along line
 *                  3 is immobile
 */
                switch(numConstraints) {
                    case 0:
                        numConstraints = 1;
                        VECTOR_COPY(glideConstraints[0], slipSystem[planeID]);
                        break;

                    case 1:
                        S1 = DotProduct(glideConstraints[0],
                                        slipSystem[planeID]);

                        if ((1.0-fabs(S1)) > eps) {
                            numConstraints = 2;

                            glideConstraints[1][0] = slipSystem[planeID][0] - 
                                                  S1*glideConstraints[0][0];
                            glideConstraints[1][1] = slipSystem[planeID][1] - 
                                                  S1*glideConstraints[0][1];
                            glideConstraints[1][2] = slipSystem[planeID][2] - 
                                                  S1*glideConstraints[0][2];

                            NormalizeVec(glideConstraints[1]);

                            cross(glideConstraints[0], glideConstraints[1],
                                  glideConstraints[2]);
                        }
                        break;

                    case 2:
                        S1 = DotProduct(glideConstraints[2],
                                        slipSystem[planeID]);

                        if (fabs(S1) > eps) {
                            numConstraints = 3;
                        }
                        break;

                }  /* end switch(numConstraints) */

            } else if ((L[i][0] > 3) && (numConstraints != 3)) {
/*
 *              Dislocation is a junction
 */
                switch(numConstraints) {
                    case 0:
                        numConstraints = 2;
                        glideConstraints[0][0] = 1.0 - (L[i][6] * L[i][6]);
                        glideConstraints[0][1] = 0.0 - (L[i][6] * L[i][7]);
                        glideConstraints[0][2] = 0.0 - (L[i][6] * L[i][8]);
                        S1 = Normal(glideConstraints[0]);

                        if (S1 > eps) {
                            glideConstraints[0][0] /= S1;
                            glideConstraints[0][1] /= S1;
                            glideConstraints[0][2] /= S1;
                        } else {
                            glideConstraints[0][0] = 0.0 - (L[i][7] * L[i][6]);
                            glideConstraints[0][1] = 1.0 - (L[i][7] * L[i][7]);
                            glideConstraints[0][2] = 0.0 - (L[i][7] * L[i][8]);
                            S1 = Normal(glideConstraints[0]);

                            if (S1 > eps) {
                                glideConstraints[0][0] /= S1;
                                glideConstraints[0][1] /= S1;
                                glideConstraints[0][2] /= S1;
                            } else {
                                glideConstraints[0][0] = 0.0 - (L[i][8] * L[i][6]);
                                glideConstraints[0][1] = 0.0 - (L[i][8] * L[i][7]);
                                glideConstraints[0][2] = 1.0 - (L[i][8] * L[i][8]);

                                S1 = Normal(glideConstraints[0]);

                                glideConstraints[0][0] /= S1;
                                glideConstraints[0][1] /= S1;
                                glideConstraints[0][2] /= S1;
                            }
                        }

                        cross(glideConstraints[0], (&L[i][6]),
                              glideConstraints[1]);

                        glideConstraints[2][0] = L[i][6];
                        glideConstraints[2][1] = L[i][7];
                        glideConstraints[2][2] = L[i][8];

                        break;

                    case 1:
                        S1 = DotProduct(glideConstraints[0],
                                        (&L[i][6]));

                        if (fabs(S1) > eps) {
                            numConstraints = 3;
                        } else {
                            cross(glideConstraints[0], (&L[i][6]),
                                  glideConstraints[1]);

                            glideConstraints[2][0] = L[i][6];
                            glideConstraints[2][1] = L[i][7];
                            glideConstraints[2][2] = L[i][8];
                        }

                        break;

                    case 2:
                        S1 = DotProduct(glideConstraints[2],
                                        (&L[i][6]));
                        if ((1.0-fabs(S1)) > eps) {
                            numConstraints = 3;
                        }

                        break;
                }
            }
        }  /* end loop creating glide constraint matrix */


	*numGlideConstraints = numConstraints;

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

        count = 0;
        memset(fscrew, 0, sizeof(fscrew));
        memset(fedge, 0, sizeof(fedge));
        memset(ffast, 0, sizeof(ffast));
        memset(dferrordvold, 0, sizeof(dferrordvold));
        ferrortol = MAX((DotProduct(fn,fn)*1.0e-16),
                        (mu*mu*lmax*lmax*1.0e-24));
        velerrortol = MAX((1.0e-14*param->rTol*param->rTol/
                           (param->deltaTT*param->deltaTT)),
                           (epsC0*epsC0*1.0e-20));
        minvelerror = velerrortol;
        correctionAdjustment = 1.0;
        notConverged = 1;

/*
 *      Set up some initial stuff for limiting the velocity corrections
 */
        minLimit[0] = minLimit[1] = minLimit[2] = -1e+100;
        minLimitSet[0] = minLimitSet[1] = minLimitSet[2] = 0;

        maxLimit[0] = maxLimit[1] = maxLimit[2] =  1e+100;
        maxLimitSet[0] = maxLimitSet[1] = maxLimitSet[2] = 0;

/*
 *      Loop over all segments (skipping zero length segments)
 */
        while ((notConverged) && (numConstraints != 3)) {

            count++;

/*
 *          If we need to include inertial terms...
 */
            if (param->includeInertia) {
                vDiff[0] = vtest[0] - vOld[0];
                vDiff[1] = vtest[1] - vOld[1];
                vDiff[2] = vtest[2] - vOld[2];
                Matrix33Vector3Multiply(massMatrix, vDiff, tmp3B);
            } else {
                VECTOR_ZERO(tmp3B);
            }

            for (i = 0; i < 3; i++) {
                ferror[i] = fn[i] - tmp3B[i];
            }

            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dferrordv[m][n] = -massMatrix[m][n];
                }
            }


            for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {
/*
 *              Skip any zero length segments...
 */
                if (zeroLenSeg[iAll]) {
                    i--;
                    continue;
                }

                if (L[i][0] > 3) {
/*
 *                  Dislocation is a junction.
 */
                    JunctionDrag(BJglide, vtest, &L[i][6],
                                 fjdrag, dfjdragdv);

                    S1 = 0.5 * L[i][2];

                    for (m = 0; m < 3; m++) {
                        ferror[m] -= S1 * fjdrag[m];
                        for (n = 0; n < 3; n++) {
                            dferrordv[m][n] -= S1 * dfjdragdv[m][n];
                        }
                    }

                } else {
                    int isLongBurg;
/*
 *                  Dislocation is a glide dislocation.  Decide which of
 *                  the three glide types and variables are appropriate.
 */
                    if (L[i][0] == 0) {
                        mobIndex = 0;
                    } else if (L[i][1] < 6) {
                        mobIndex = 1;
                    } else {
                        mobIndex = 2;
                    }

/*
 *                  A few of the functions need to know if they're dealing
 *                  with a long or short burgers vector, so set a flag here.
 */
                    isLongBurg = (L[i][0] == 0 ? 1 : 0);


                    if (L[i][3] > 0.0) {
/*
 *                      Dislocation in in the screw trench. Use the 
 *                      'screw' drag function.
 */
                        ScrewDrag(param, fscrew[i], vtest, burg[i], ndir[i],
                                  L[i][1], mobIndex, dfsdv, isLongBurg);

                        S1 = 0.5 * L[i][3];

                        for (m = 0; m < 3; m++) {
                            ferror[m] -= S1 * fscrew[i][m];
                            for (n = 0; n < 3; n++) {
                                dferrordv[m][n] -= S1 * dfsdv[m][n];
                            }
                        }

                    } else if (L[i][4] > 0.0) {
/*
 *                      Dislocation is in the non-screw trench.  Use the
 *                      'edge' drag function.
 */
                        EdgeDrag(param, vtest, burg[i], ndir[i],
                                 L[i][1], mobIndex,
                                 fedge[i], dfedv, isLongBurg);

                        S1 = 0.5 * L[i][4];

                        for (m = 0; m < 3; m++) {
                            ferror[m] -= S1 * fedge[i][m];
                            for (n = 0; n < 3; n++) {
                                dferrordv[m][n] -= S1 * dfedv[m][n];
                            }
                        }
                    }

                    if (L[i][5] > 0.0) {
/*
 *                      Dislocation has non-trench components.  Use the
 *                      'fast' mobility function.
 */
                        FastDrag(param, vtest, burg[i], ndir[i],
                                 L[i][1], mobIndex, &L[i][6],
                                 ffast[i], dffdv);

                        S1 = 0.5 * L[i][5];

                        for (m = 0; m < 3; m++) {
                            ferror[m] -= S1 * ffast[i][m];
                            for (n = 0; n < 3; n++) {
                                dferrordv[m][n] -= S1 * dffdv[m][n];
                            }
                        }
                    }
                }
            }

/*
 *          Calculate the velocity correction.  Note: if we cannot
 *          invert dferrordv, just return an error
 */
            if (count > 1) {
#if 0
                avgFactor = 0.50;
#else
                avgFactor = 1.00;
#endif
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
                printf("dferrordv = \n    ");

                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        printf("%.13e ", dferrordv[m][n]);
                    }
                    printf("\n    ");
                }
                printf("\n");

                for (i = 0; i < numNonZeroLenSegs; i++) {
                    printf("L[%d][0]   = %e\n", i, L[i][0]);
                    printf("L[%d][1]   = %e\n", i, L[i][1]);
                    printf("L[%d][2]   = %e\n", i, L[i][2]);
                    printf("L[%d][3]   = %e\n", i, L[i][3]);
                    printf("L[%d][4]   = %e\n", i, L[i][4]);
                    printf("L[%d][5]   = %e\n", i, L[i][5]);
                    printf("L[%d][6-8] = %e %e %e\n", i,
                           L[i][6], L[i][7], L[i][8]);
                }

                PrintNode(node);
#if 0
                DumpNodeAndParams(home, node);
#endif

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

/*
 *          Allow a few loop iterations in which to get past any potentially
 *          trash velocity estimates, then begin putting boundaries on the
 *          test velocities based on the correction values.  This can help
 *          prevent cases where the test velocity simply oscillates between
 *          a couple values without actually converging.
 *
 *          Need to do this before vtest gets updated!
 */
            if (count > 5) {
                real8 velNew;

                for (i = 0; i < 3; i++) {
/*
 *                  This looks backwards, but is okay since correction is
 *                  subtracted from vtest.
 *
 *                  If the correction will decrement the test velocity, use
 *                  the current test velocity as the new upper bound.  if the
 *                  correction will increment the test velocity, use the
 *                  current test velocity as the new lower bound.
 */
                    if (correction[i] > 0.0) {
                        maxLimit[i] = vtest[i];
                        maxLimitSet[i] = 1;
                    } else {
                        minLimit[i] = vtest[i];
                        minLimitSet[i] = 1;
                    }

                    velNew = vtest[i] - correction[i];

/*
 *                  If the correction will move the test velocity outside of
 *                  the current boundaries, modify the correction value to
 *                  generate a test velocity within the range.  (If both an
 *                  upper bound and lower bound have been defined, choose a
 *                  correction that will result in a test velocity at the
 *                  midpoint of the permitted velocity range.  Otherwise
 *                  set the correction to place the new test velocity at
 *                  the currently defined boundary.)
 */
                    if ((velNew < minLimit[i]) || (velNew > maxLimit[i])) {
                        if (!maxLimitSet[i]) {
                            correction[i] = vtest[i] - minLimit[i];
                        } else if (!minLimitSet[i]) {
                            correction[i] = vtest[i] - maxLimit[i];
                        } else {
                            correction[i] = vtest[i] - 
                                            0.5 * (maxLimit[i] + minLimit[i]);
                        }
                    }
                }
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
 *          Preserve the vtest value for the iteration that has the
 *          lowest absolute error.  If the function fails to converge
 *          we may be able to use this preserved velocity.
 */
            if (velerror < minvelerror) {
                minvelerror = velerror;
                VECTOR_COPY(vtmp, vtest);
            }

/*
 *          If the initial velocity guess was too far off, it's
 *          better to just zero the vtest value.  So, after the
 *          first iteration, if the force error is too large (an
 *          indication the initial velocity guess was way off),
 *          just zero the vtest value.
 */
            if ((count == 1) && (forceerror > 1.0e+01 * ferrortol)) {
                VECTOR_ZERO(vtest);
            }


            if (++iterCnt > maxIter) {
/*
 *              We didn't converge on an answer, but if there was an iteration
 *              of the above loop that resulted in an absolute error within
 *              the absolute error tolerance, go ahead and use the velocity
 *              calculated during that iteration.
 */
                if (minvelerror < velerrortol) {
                    VECTOR_COPY(vtest, vtmp);
                    break;
                }

/*
 *              On occassion the mobility function doesn't converge, even
 *              though the correction values are only a fraction of the
 *              test velocity.  In those cases, go ahead and accept the
 *              magnitude of the correction is less than 1.5% of the
 *              magnitude of the test velocity.
 */
                if ((velerror / DotProduct(vtest, vtest)) < 2.25e-4) {
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
                    memset(ffast, 0, sizeof(ffast));
                    memset(dferrordvold, 0, sizeof(dferrordvold));
                    velerror = 10.0 * velerrortol;
                    minvelerror = velerror;
                    correctionAdjustment = 0.5;
                    vtest[0] = node->vX;
                    vtest[1] = node->vY;
                    vtest[2] = node->vZ;

                    minLimit[0] = minLimit[1] = minLimit[2] = -1e+100;
                    minLimitSet[0] = minLimitSet[1] = minLimitSet[2] = 0;

                    maxLimit[0] = maxLimit[1] = maxLimit[2] =  1e-100;
                    maxLimitSet[0] = maxLimitSet[1] = maxLimitSet[2] = 0;

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

                return(1);
            }

        }  /* while ((notConverged) && (numConstraints != 3)) */

/*
 *      We were able to converge sufficiently on a velocity, so
 *      apply the glide constraints and reset the nodal velocity
 *      to the new values.
 *
 *      If there is one planar constraint, the first row in the glide
 *      constraint matrix contains the plane normal direction.  Any velocity
 *      projection in this direction is removed.  If there are two planar
 *      constraints, the third row in the constrainmat contains the only
 *      direction that the dislocation may move in and the vector projection
 *      in that direction is returned.  If there are three planar constraints
 *      the node simply can't move.
 */
        switch(numConstraints) {
            case 1:
                S1 = DotProduct(vtest, glideConstraints[0]);
                vtest[0] -= S1 * glideConstraints[0][0];
                vtest[1] -= S1 * glideConstraints[0][1];
                vtest[2] -= S1 * glideConstraints[0][2];
                break;
            case 2:
                S1 = DotProduct(vtest, glideConstraints[2]);
                vtest[0] = S1 * glideConstraints[2][0];
                vtest[1] = S1 * glideConstraints[2][1];
                vtest[2] = S1 * glideConstraints[2][2];
                break;
            case 3:
                VECTOR_ZERO(vtest);
        }
        
/*
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame.
 */
        if (param->useLabFrame) {
            real8 vtestRot[3];
            Matrix33Vector3Multiply(param->rotMatrix, vtest, vtestRot);
            VECTOR_COPY(vtest, vtestRot);
        }

        node->vX = vtest[0];
        node->vY = vtest[1];
        node->vZ = vtest[2];

/*
 *      Sometimes nodes oscillate at high frequency causing timestep to drop.
 *      Usually, these are 2-nodes that have short segments attached.  So, if
 *      we find a 2-node that has a least 1 short segment, and the velocity
 *      of the node has changed direction since the previous cycle, slow
 *      the node down for this cycle.  Hopefully this will help dampen the
 *      oscillation and improve the timestep.
 */
        if ((numNbrs == 2) && (hasShortSeg)) {
            if (((node->vX * node->oldvX) < 0.0) ||
                ((node->vY * node->oldvY) < 0.0) ||
                ((node->vZ * node->oldvZ) < 0.0)) {
                node->vX *= shortSegRatio;
                node->vY *= shortSegRatio;
                node->vZ *= shortSegRatio;
            }
        }

#ifdef _ARLFEM
/*
 *      The velocity of any surface node along the negative surface
 *      normal direction should be zero to prevent the node moving into
 *      the box.  Make a call to adjsut the velocity accordingly.
 *
 *      Note: If the node moves outside the box, it's allowed here, then
 *      position is adjusted in AdjustNodePosition().
 */
        AdjustSurfaceNodeVel(home, node);
#endif

#ifdef NAN_CHECK
        if ((isnan(node->vX) || isinf(node->vX)) ||
            (isnan(node->vY) || isinf(node->vY)) ||
            (isnan(node->vZ) || isinf(node->vZ))) {
            PrintNode(node);
            for (i = 0; i < node->numNbrs; i++) {
                nbrNode = GetNeighborNode(home, node, i);
                if (nbrNode == (Node_t *)NULL) {
                    printf("WARNING: Neighbor not found at %s line %d\n",
                           __FILE__, __LINE__);
                    continue;
                }
                PrintNode(nbrNode);
            }
            Fatal("(%d, %d) velocity: %lf %lf %lf",
                  node->myTag.domainID, node->myTag.index,
                  node->vX, node->vY, node->vZ);
        }
#endif

/*
 *      Adjust the node's velocity if the node has any constraints that
 *      affect its velocity.
 */
        ApplyNodeConstraintsToVelocity(home, node);

        return(0);
}
