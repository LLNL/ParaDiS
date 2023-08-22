/**************************************************************************
 *
 *      Module:       Mobility_BCC_nl.c
 *      Description:  Contains functions for calculating isotropic
 *                    non-linear mobility of nodes in BCC metals.
 *                    Based on the Arsenlis matlab code mobbccnl6b.m
 *
 *      Includes public functions:
 *            Mobility_BCC_nl()
 *
 *     Includes private functions:
 *            EdgeDrag()
 *            JunctionDrag()
 *            ScrewDrag()
 *            ScrewMobility()
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

//static int doDebug = 0;
//static int dumpDebugData = 0;

static real8 Linear;
static real8 NonLinear;
static real8 V0;

/*
 *      Some of the array/matrix sizes are dependent on the number of
 *      segments connected to a node or the number of different types
 *      of screw segments attached to the node.  By defining a large
 *      enough upper limit on these numbers a node might have, we can
 *      allocate the arrays/matrices to a maximum size at the start.
 *      Makes life much simpler.
 */
#define MAX_NSCREW MAX_NBRS



/*
 *      The mobility depends on a number of values that are dependent
 *      on the input parameters.  This function is just called once
 *      during initialization to set those values appropriately.
 */
void InitMob_BCC_nl(Home_t *home)
{
        real8   p, q, alpha, y, hratio, pressureInGPa;
        Param_t *param;

        param = home->param;

        p = param->MobExpP;
        q = param->MobExpQ;
        alpha = param->MobAlpha;

        pressureInGPa = param->pressure / 1.0e+09;

        y = pow(alpha, p);
        hratio = param->MobDeltaH0 / (BOLTZMANNS_CONST * param->TempK);
        Linear = param->MobCoeffA * param->MobCoeffC0 *
                 (1.0 - param->MobCoeffC * param->TempK /
                 param->meltTemp + param->MobCoeffD *
                 pressureInGPa);
        NonLinear = Linear / exp(-hratio*pow((1.0-y),q)) /
                    (1.0+p*q*hratio*pow((1.0-y),(q-1.0))*y);
        V0 = Linear*alpha-NonLinear*alpha*exp(-hratio*pow((1.0-y),q));

        return;
}


/*
 *      Just another function to dump information used for debugging
 */
static void DumpNodeAndParams(Home_t *home, Node_t *node)
{
        int j;
        Node_t *nbrNode;

        printf("rTol     = %.15e\n", home->param->rTol);
        printf("nextDT   = %.15e\n", home->param->deltaTT);

        printf("MobLine  = %.15e\n", home->param->MobLine);
        printf("MobGlide = %.15e\n", home->param->MobGlide);
        printf("MobClimb = %.15e\n", home->param->MobClimb);
        printf("Peierls  = %.15e\n", home->param->MobPeierls);
        printf("DeltaH0  = %.15e\n", home->param->MobDeltaH0);
        printf("Boltz    = %.15e\n", BOLTZMANNS_CONST);
        printf("TempK    = %.15e\n", home->param->TempK);
        printf("P        = %.15e\n", home->param->MobExpP);
        printf("Q        = %.15e\n", home->param->MobExpQ);
        printf("MobAlpha = %.15e\n", home->param->MobAlpha);
        printf("MobLinear    = %.15e\n", Linear);
        printf("MobNonLinear = %.15e\n", NonLinear);
        printf("MobV0        = %.15e\n", V0);

        printf("shearModulus = %.15e\n", home->param->shearModulus);
        printf("poisson ratio= %.15e\n", home->param->pois);
        printf("Burgers mag. = %.15e\n", home->param->burgMag);

        PrintNode(node);
        for (j = 0; j < node->numNbrs; j++) {
            if ((nbrNode = GetNeighborNode(home, node, j)) == (Node_t *)NULL) {
                continue;
            }
            PrintNode(nbrNode);
        }

        return;
}

static void JunctionDrag(real8 Beclimb, real8 Beline, real8 vel[3],
                         real8 burg[3], real8 t[3], /* real8 fjdrag[3],*/
                         real8 dfjdragdv[3][3])
{
        real8 tmp[3][3];

/*
 *      dfjdragdv = beclimb.*(eye(3)-(t'*t)) + beline.*(t'*t)
 */
        Vec3TransposeAndMult(t, tmp);
        
        dfjdragdv[0][0] = Beclimb * (1.0 - tmp[0][0]) + Beline * tmp[0][0];
        dfjdragdv[0][1] = Beclimb * (0.0 - tmp[0][1]) + Beline * tmp[0][1];
        dfjdragdv[0][2] = Beclimb * (0.0 - tmp[0][2]) + Beline * tmp[0][2];

        dfjdragdv[1][0] = Beclimb * (0.0 - tmp[1][0]) + Beline * tmp[1][0];
        dfjdragdv[1][1] = Beclimb * (1.0 - tmp[1][1]) + Beline * tmp[1][1];
        dfjdragdv[1][2] = Beclimb * (0.0 - tmp[1][2]) + Beline * tmp[1][2];

        dfjdragdv[2][0] = Beclimb * (0.0 - tmp[2][0]) + Beline * tmp[2][0];
        dfjdragdv[2][1] = Beclimb * (0.0 - tmp[2][1]) + Beline * tmp[2][1];
        dfjdragdv[2][2] = Beclimb * (1.0 - tmp[2][2]) + Beline * tmp[2][2];

        return;
}


static void EdgeDrag(real8 Beglide, real8 Beclimb, real8 Beline, real8 vel[3],
                     real8 burg[3], real8 climbdir[3], /* real8 fedrag[3],*/
                     real8 dfedragdv[3][3])
{
        int   m, n;
        real8 linedir[3], glidedir[3];
        real8 tmp1[3][3], tmp2[3][3], tmp3[3][3];

        glidedir[0] = burg[0];
        glidedir[1] = burg[1];
        glidedir[2] = burg[2];

        Normalize(&glidedir[0], &glidedir[1], &glidedir[2]);

        cross(glidedir, climbdir, linedir);
    
        Vec3TransposeAndMult(glidedir, tmp1);
        Vec3TransposeAndMult(climbdir, tmp2);
        Vec3TransposeAndMult(linedir,  tmp3);
        
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] = Beglide * tmp1[m][n] +
                                  Beclimb * tmp2[m][n] +
                                  Beline  * tmp3[m][n]; 
            }
        }

        return;
}


static void fscale(Param_t *param, real8 vin, real8 burg[3],
                   real8 *fout, real8 *dfdv)
{
        int   path;
        real8 err, derrdf, screwerrortol, pfactor1, pfactor2;
        real8 fmag, vmag, vtrans, vsign, hratio, fratio, lratio, dvdf;
        real8 peierls, deltaH0, tempK, omega, p, q;
        real8 burgnormal, eps;
        real8 alpha, v0, mobLinear, mobNonLinear;

        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);
        deltaH0 = param->MobDeltaH0;
        p = param->MobExpP;
        q = param->MobExpQ;
        tempK = param->TempK;

        mobLinear = Linear / param->burgMag;
	mobNonLinear = NonLinear / param->burgMag; 
        v0 = V0 / param->burgMag;
        peierls = param->MobPeierls * burgnormal;
        alpha = param->MobAlpha;

        omega = mobNonLinear / peierls;

        fmag = peierls;
/*
 *      vin may already be an absolute value now... this may be superfluous
 */
        vmag = fabs(vin);
        vsign = (real8)Sign(vin);
        vtrans = alpha * mobLinear - v0;
        hratio = deltaH0 / (BOLTZMANNS_CONST*tempK);
        screwerrortol = 1.0e-10;
        eps = 1.0e-12;

        if (vmag >= vtrans) {
/*
 *          When the velocity is above a specific threshold, use
 *          linear mobility
 */
            *fout = (vmag + v0) * peierls / mobLinear * vsign;
            *dfdv = peierls / mobLinear;
        } else if (vmag < eps) {
/*
 *          We need this check for vmag == 0 or we end up with a 
 *          divide by zero calculating lratio, which causes havoc
 */
            *fout = 0.0;
            fratio = 0.0;
            dvdf = omega * exp(-hratio);
            *dfdv = 1.0 / dvdf;
        } else {
            err = 1.0;
            derrdf = 0.0;
            path = 1;

            pfactor1 = 0.50 * peierls;
            pfactor2 = 0.60 * peierls;

            while (err > screwerrortol) {
                if ((((fmag > pfactor1) && (path == 1)) ||
                     ((fmag > pfactor2) && (path == 2))) &&
                    (fmag > (vmag / omega))) {
                    lratio = log((omega * fmag) / vmag);

                    err = fmag - peierls *
                          pow(1.0-pow(lratio/hratio,1.0/q),1.0/p);
                    derrdf = 1.0 + peierls / (p*q*hratio*fmag) * 
                             pow(1.0-pow(lratio/hratio,1.0/q),1.0/p-1.0) *
                             pow(lratio/hratio,(1.0/q)-1.0);
                    path = 1;
                } else {
                    fratio = fmag / peierls;
                    err = fmag - vmag /
                          (omega*exp(-hratio*pow(1.0-pow(fratio,p),q)));
                    derrdf = 1.0 + vmag /
                          (omega*exp(-hratio*pow(1.0-pow(fratio,p),q))) *
                          hratio * pow(1.0-pow(fratio,p),q-1.0) *
                          q * p * pow(fratio,p-1.0) / peierls;
                    path = 2;
                }

                fmag -= err / derrdf;

                if (fmag < 0.0) {
                    fmag = 1.0e-16;
                }

                err = fabs(err/(fmag+eps));

            }  /* while (error > screwerrortol) */

            err = 1.0;

            while (err > screwerrortol) {
                fratio = fmag / peierls;
                err = vmag - omega * fmag *
                      exp(-hratio * pow(1.0-pow(fratio,p),q));
                derrdf = -omega * exp(-hratio * pow(1.0-pow(fratio,p),q)) *
                         (1.0 + p*q*hratio*pow(1.0-pow(fratio,p),(q-1.0) *
                          pow(fratio,p)));
                fmag -= err / derrdf;
                err = fabs(err / (vmag+eps));
            }

            *fout = fmag * vsign;
            fratio = fmag / peierls;
            dvdf = omega * exp(-hratio*pow(1.0-pow(fratio,p),q)) *
                   (1.0 + p * q * hratio * pow(1.0-pow(fratio,p),q-1.0) *
                    pow(fratio,p));
            *dfdv = 1.0 / dvdf;

        }  /* vmag < vtrans && vmag != 0.0 */

#ifdef NAN_CHECK
        if (isnan(*fout) || isinf(*fout)) {
            Fatal("fscale: fout = %lf", *fout);
        }

        if (isnan(*dfdv) || isinf(*dfdv)) {
            Fatal("fscale: dfdv = %lf = 1.0 / dvdf\n"
                  "  where dvdf = %lf", *dfdv, dvdf);
        }
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
        real8 dvdirdv[3][3];
        real8 vMag, vMagInv, velNorm, vdir[3];
        real8 fout, dfdv;


        Bsline = 1.0 / param->MobLine;

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

        fscale(param, vMag, burg, &fout, &dfdv);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfsdragdv[m][n] = Bsline * linedirTM[m][n];
            }
        }

        Matrix33Vector3Multiply(dfsdragdv, vel, tmp3);

        for (i = 0; i < 3; i++) {
            finit[i] = tmp3[i] + fout*vdir[i];
        }

        velNorm = Normal(vel);

        if ((velNorm > eps) && (vMag/velNorm > eps)) {
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
                    dfsdragdv[m][n] +=  dfdv * glidedirmat[m][n];
                }
            }
        }

        return;
}


/**************************************************************************
 *
 *      Function:     Mobility_BCC_nl
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
 *      Returns: 0 on success
 *               1 if function could not converge on a velocity
 *
 *************************************************************************/
int Mobility_BCC_nl(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, j, m, n, iAll, numNbrs, numscrews, thisScrew;
        int     iterCnt, maxIter, maxIter2, count;
        int     isJunction;
        int     zeroLenSeg[MAX_NBRS];
        int     edgeExists[MAX_NBRS], screwExists[MAX_NBRS];
        int     numNonZeroLenSegs = 0;
        real8   abserrortol, abserror, minabserror;
        real8   normerrortol, normerror, vline;
        real8   Beclimb, Beline, Beglide;
        real8   eps, rt2, burgProd, mag, avgFactor, correctionAdjustment;
        real8   vtest[3], fn[3], rt[3], zero3[3], tmp3[3], ferror[3];
        real8   dfsdv[3][3], dferrordv[3][3];
        real8   dflindv[3][3], dfedragdv[3][3], dfjdragdv[3][3];
        real8   dferrordvold[3][3]; 
        real8   LTotal[MAX_NBRS], LScrew[MAX_NBRS], LEdge[MAX_NBRS], L2;
        real8   linedir[MAX_NBRS][3], ndir[MAX_NBRS][3], burg[MAX_NBRS][3];
        real8   screwb[MAX_NSCREW][3], dfndfscrew[MAX_NSCREW];
        real8   fscrew[MAX_NSCREW][3];
        real8   tmp, scalar, vdirNorm, vtestNorm, screwNorm;
        real8   rhs[2], tmp2[2], tmp22[2][2], tmp22Inv[2][2];
        real8   vdir[3], perpdir[3], screwlinedir[3];
        real8   corr[3], correction[3];
        real8   vtmp[3]={0.0, 0.0, 0.0}, fnstar[3];
        Param_t *param;
        Node_t  *nbrNode;

        real8  (*invdferrordv)[3]     =  mobArgs->invDragMatrix;
     // real8  (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int     *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdferrordv);

        *numGlideConstraints = 0;

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

        param = home->param;
        eps = 1.0e-12;
        iterCnt = 0;
        maxIter = 50;
        maxIter2 = 100;
        isJunction = 0;

#if 0
/*
 *      Temporary debug code... 
 */
        if (doDebug) {
            DumpNodeAndParams(home, node);
        }
#endif

        abserrortol = 1.0e-14 * (param->rTol * param->rTol) /
                      (param->deltaTT * param->deltaTT);
        normerrortol = 1.0e-12;

        numscrews = 0;

        Beclimb  = 1.0 / param->MobClimb;
        Beline   = 1.0 / param->MobLine;
        Beglide  = 1.0 / param->MobGlide;

        memset(zero3, 0, sizeof(zero3));

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

        numNbrs = node->numNbrs;

/*
 *      This function saves various pieces of info for each neighbor
 *      node using static arrays for speed rather than doing memory
 *      allocations.  If the array is too  small, abort with an error
 *      to let the caller increase the array size.
 */
        if (numNbrs > MAX_NBRS) {
            Fatal("Mobility_BCC_nl: Node segment count (%d) exceeds\n"
                  "maximum allowed (%d) in static arrays.  Increase MAX_NBRS\n"
                  "and recompile\n", node->numNbrs, MAX_NBRS);
        }

/*
 *      Initialize some arrays for later...
 */
        vtest[0] = node->vX;
        vtest[1] = node->vY;
        vtest[2] = node->vZ;

        memset(ndir, 0, sizeof(ndir));

/*
 *      Make a temporary copy of the total forces
 */
        fn[0] = node->fX;
        fn[1] = node->fY;
        fn[2] = node->fZ;

/*
 *      If necessary, rotate the force and current velocity
 *      vectors from the laboratory frame to the crystal frame
 */
        if (param->useLabFrame) {
            real8 forceRot[3], velRot[3];
            Matrix33Vector3Multiply(param->rotMatrixInverse, fn, forceRot);
            Matrix33Vector3Multiply(param->rotMatrixInverse, vtest, velRot);
            VECTOR_COPY(fn, forceRot);
            VECTOR_COPY(vtest, velRot);
        }

/*
 *      Loop over all arms attached to the node
 */
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
 *          If necessary, rotate the line sense from the laboratory
 *          frame to the crystal frame
 */
            if (param->useLabFrame) {
                real8 lineRot[3];
                Matrix33Vector3Multiply(param->rotMatrixInverse, rt, lineRot);
                VECTOR_COPY(rt, lineRot);
            }

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
 *          If necessary, rotate the burgers vector from the laboratory
 *          frame to the crystal frame
 */
            if (param->useLabFrame) {
                real8 burgRot[3];
                Matrix33Vector3Multiply(param->rotMatrixInverse,burg[i],burgRot);
                VECTOR_COPY(burg[i], burgRot);
            }

/*
 *          If the node is a junction node (i.e. it has any
 *          non 111 type burgers vectors) set a flag.  We'll
 *          need to know this later to determine how we iterate
 *          to find the velocity.
 */
            if (fabs(burg[i][0]*burg[i][1]*burg[i][2]) < eps) {
                isJunction = 1;
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

        }  /* loop over neighbors */

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
                JunctionDrag(Beclimb, Beline, zero3, burg[i], linedir[i],
                                 /* fjdrag,*/ dfjdragdv);
                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        dflindv[m][n] += (0.5 * LTotal[i]) * dfjdragdv[m][n];
                    }
                }
            } else if (edgeExists[i]) {
                EdgeDrag(Beglide, Beclimb, Beline, zero3, burg[i],
                         ndir[i], /* fedrag,*/ dfedragdv);
                for (m = 0; m < 3; m++) {
                    for (n = 0; n < 3; n++) {
                        dflindv[m][n] += (0.5 * LEdge[i]) * dfedragdv[m][n];
                    }
                }
            }
        }  /* Loop over all segments */

/*
 *      Loop over all segments (skipping zero length segments)
 */
        numscrews = 0;

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
            }
        }

/*
 *      For non-junction nodes with only 1 screw burgers vector,
 *      the velocity component along the line can be computed
 *      directly without iterating.  So calculate this value now
 *      and temporarily subtract the corresponding force from
 *      the node.  Later after we've iterated to converge on
 *      the non-linear components of the velocity, the linear
 *      velocity will be added back. 
 */
#if 0
        if ((numscrews == 1) && (isJunction == 0)) {
            tmp = 0.0;
            for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {
                if (zeroLenSeg[iAll]) {
                    i--;
                    continue;
                }
                tmp += LEdge[i];
            }
            scalar = -dfndfscrew[0] * 0.5 * Beline - 0.5 * Beglide * tmp;
            screwNorm = Normal(screwb[0]);
            for (i = 0; i < 3; i++) {
                screwlinedir[i] = screwb[0][i] / screwNorm;
                fnstar[i] = fn[i];
            }
#endif
#if 1
        if (numscrews == 1) {
            screwNorm = Normal(screwb[0]);
            for (i = 0; i < 3; i++) {
                screwlinedir[i] = screwb[0][i] / screwNorm;
                fnstar[i] = fn[i];
            }
            Matrix33_Vmul(tmp3, screwlinedir, dflindv );   // tmp3 = screwlinedir * dflindv
            tmp = DotProduct(tmp3, screwlinedir);
            scalar = -dfndfscrew[0] * 0.5 * Beline - tmp;
#endif
            tmp = DotProduct(fnstar, screwlinedir);
            for (i = 0; i < 3; i++) {
                fn[i] = fnstar[i] - tmp * screwlinedir[i];
            }
            tmp = DotProduct(vtest, screwlinedir);
            for (i = 0; i < 3; i++) {
                vtest[i] = vtest[i] - screwlinedir[i] * tmp;
            }
            vline = -(DotProduct(fnstar, screwlinedir)) / scalar;
        }

        count = 0;
        normerror = 1.0;
        memset(fscrew, 0, sizeof(fscrew));
        memset(dferrordvold, 0, sizeof(dferrordvold));
        abserror = 10.0 * abserrortol;
        minabserror = abserror;
        correctionAdjustment = 1.0;

        while (normerror > normerrortol) {

            count++;

            Matrix33Vector3Multiply(dflindv, vtest, tmp3);

            for (i = 0; i < 3; i++) {
                ferror[i] = fn[i] - tmp3[i];
            }


            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dferrordv[m][n] = -dflindv[m][n];
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

/*
 *          For non-junction nodes with only 1 screw burgers vector,
 *          the velocity component along the line can be computed
 *          directly without iterating.  As such, we only need to iterate
 *          to converge on the non-linear portion of velocity for these
 *          nodes, and only after convergence will the velocity along
 *          the line be added back into the final calculated velocity.
 */
#if 0
            if ((numscrews == 1) && (isJunction == 0)) {
#endif
#if 1
            if (numscrews == 1) {
#endif

                tmp = DotProduct(vtest, screwlinedir);
                for (i = 0; i < 3; i++) {
                    vdir[i] = vtest[i]  - tmp * screwlinedir[i];
                }
                vdirNorm = Normal(vdir);
                vtestNorm = Normal(vtest);
                if ((vtestNorm < eps) || (vdirNorm/vtestNorm < eps)) {
                    vdir[0] =  screwb[0][0];
                    vdir[1] = -screwb[0][1];
                    vdir[2] =  0.0;
                    vdirNorm = Normal(vdir);
                }
                for (i = 0; i < 3; i++) {
                    vdir[i] = vdir[i] / (vdirNorm+eps);
                }
                cross(vdir, screwlinedir, perpdir);

                Matrix33_Vmul(tmp3, vdir   , dferrordv );   // tmp3 = vdir    * dferrordv
                tmp22[0][0] = DotProduct(tmp3, vdir);
                tmp22[0][1] = DotProduct(tmp3, perpdir);

                Matrix33_Vmul(tmp3, perpdir, dferrordv );   // tmp3 = perpdir * dferrordv
                tmp22[1][0] = DotProduct(tmp3, vdir);
                tmp22[1][1] = DotProduct(tmp3, perpdir);

                rhs[0] = DotProduct(ferror, vdir);
                rhs[1] = DotProduct(ferror, perpdir);

                Matrix22_Inverse(tmp22Inv, tmp22);
                Matrix22_Vmul   (tmp2,tmp22Inv,rhs);

                corr[0] = tmp2[0];
                corr[1] = tmp2[1];

                for (i = 0; i < 3; i++) {
                    correction[i] = corr[0]*vdir[i] +
                                    corr[1]*perpdir[i];
                }

            } else {
/*
 *              Calculate the velocity correction for junction
 *              nodes and any node with 2 or more screw burgers vectors.
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
 *              Adjust the correction values by a factor.  This factor
 *              is initially 1.0, but in some cases, the function is 
 *              unable to converge on an answer because the corrrection
 *              oscillates between a couple values.  In those cases,
 *              the adjustment factor will be lowered, which can help
 *              achieve convergence.
 */
                for (i = 0; i < 3; i++) {
                    correction[i] *= correctionAdjustment;
                }

            }

            for (i = 0; i < 3; i++) {
                vtest[i] -= correction[i];
            }

            normerror = DotProduct(ferror,ferror) / (DotProduct(fn,fn)+eps);
            abserror = DotProduct(correction,correction);

/*
 *          If the initial velocity guess was too far off, it's
 *          better to just zero the vtest value.  So, after the
 *          first iteration, if the normerror is too large (an
 *          indication the initial velocity guess was way off),
 *          just zero the vtest value.
 */
            if ((count == 1) && (normerror > 1.0e-10)) {
                for (i = 0; i < 3; i++) {
                    vtest[i] = 0.0;
                }
            }

/*
 *          Preserve the vtest value for the iteration that has the
 *          lowest absolute error.  If the function fails to converge
 *          we may be able to use this preserved velocity.
 */
            if (abserror < minabserror) {
                minabserror = abserror;
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
                if (minabserror < abserrortol) {
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
                    normerror = 1.0;
                    count = 0;
                    memset(fscrew, 0, sizeof(fscrew));
                    memset(dferrordvold, 0, sizeof(dferrordvold));
                    abserror = 10.0 * abserrortol;
                    minabserror = abserror;
                    correctionAdjustment = 0.5;
                    vtest[0] = node->vX;
                    vtest[1] = node->vY;
                    vtest[2] = node->vZ;
/*
 *                  If necessary, rotate the current velocity vector from
 *                  the laboratory frame to the crystal frame
 */
                    if (param->useLabFrame) {
                        real8 velRot[3];
                        Matrix33Vector3Multiply(param->rotMatrixInverse,
                                                vtest, velRot);
                        VECTOR_COPY(vtest, velRot);
                    }

                    continue;
                }
#if 0
/*
 *              Temporary debug code to be removed in near future.
 */
                if ((home->param->deltaTT < 5.0e-13) && dumpDebugData) {
                    printf("mobility: (%d,%d) no convergence after %d iterations\n", 
                           node->myTag.domainID, node->myTag.index, iterCnt);
                    printf("normerror = %e, correction = %e %e %e\n",
                           normerror, tmp3[0], tmp3[1], tmp3[2]);
                    printf("abserror = %e, abserrortol = %e\n",
                           abserror, abserrortol);
                    printf("vtest = %e %e %e\n", vtest[0], vtest[1], vtest[2]);
  
                    DumpNodeAndParams(home, node);
                }
#endif
                return(1);
            }

        }  /* while (normerror > normerrortol) */

/*
 *      We were able to converge on a velocity or at least find one
 *      within the error tolerances, so set the nodal velocity
 *      to the new values and (if necessary) add the velocity along
 *      the line back in if necessary.
 */
        if ((numscrews == 1) && (isJunction == 0)) {
            for (i = 0; i < 3; i++) {
                vtest[i] = vtest[i] + screwlinedir[i] * vline;
            }
        }

/*
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame
 */
        if (param->useLabFrame) {
            real8 vRot[3];
            Matrix33Vector3Multiply(param->rotMatrix, vtest, vRot);
            VECTOR_COPY(vtest, vRot);
        }

        node->vX = vtest[0];
        node->vY = vtest[1];
        node->vZ = vtest[2];

        return(0);
}
