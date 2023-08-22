/**************************************************************************
 *
 *      Description:  Contains functions for calculating mobility of nodes
 *                    in BCC Fe.  Includes non-linear components for screw,
 *                    linear components for edge, and optional inertial components.
 *
 *                    Note: In this version of mobility, nodes with 4 or
 *                    more segments are pinned, and junction nodes are
 *                    allowed to move only along the junction line to
 *                    allow the junction to expand or contract.
 *      
 *
 *                    Note: This mobility function is based on
 *                          Mobility_BCC_Ta_nl.c with primary modifications
 *                          to the fscalescrew() and fscaleedge() functions.
 *
 *      Includes public functions:
 *            InitMob_BCC_Fe_nl()
 *            Mobility_BCC_Fe_nl()
 *
 *      Includes private functions:
 *            EdgeDrag()
 *            ScrewDrag()
 *            FscaleEdge()
 *            FscaleScrew()
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
 *      Define a number of the mobility-specific constants that
 *      have been derived from information from atomistic siumulations
 *      and declare a few related global variables that are used in
 *      multiple functions.
 */
#define ALPHA_SCREW       3.3
#define BETA_SCREW        2.0e-4
#define B_SCREW_PRECONST  2.48
#define PEIERLS_SCREW     1.2e9
#define N_SCREW           10

#define PEIERLS_EDGE110   1.0e8
#define N_EDGE            5

static real8 C0;
static real8 V0Edge;

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
void InitMob_BCC_Fe_nl(Home_t *home)
{
        C0 = 2000 / home->param->burgMag;
        V0Edge = C0;

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

        printf("C0           = %.15e\n", C0);
        printf("V0Edge       = %.15e\n", V0Edge);

        printf("SPeierls     = %.15e\n", PEIERLS_SCREW);
        printf("SAlpha       = %.15e\n", ALPHA_SCREW);
        printf("SBeta        = %.15e\n", BETA_SCREW);
        printf("SPreconstB   = %.15e\n", B_SCREW_PRECONST);
        printf("SN           = %d\n",    N_SCREW);
        printf("EPeierls110  = %.15e\n", PEIERLS_EDGE110);
        printf("EN           = %d\n",    N_EDGE);
        printf("Boltz        = %.15e\n", BOLTZMANNS_CONST);
        printf("shearModulus = %.15e\n", home->param->shearModulus);
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


static void FscaleEdge(Param_t *param, real8 feinit, real8 vin, real8 burg[3],
                       real8 *fout, real8 *dfdv, real8 *dfdvlin, real8 *d2fdv2)
{
        real8  temp;
        real8  peierls, burgnormal;
        real8  vmag, vsign;
        real8  b, bt; 
        real8  fdrag, dfdragdv;
        real8  fmag, dfmagdv;

        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);
        temp = param->TempK;

        peierls = PEIERLS_EDGE110*burgnormal;
        b = 370.1/1e6;
        bt = temp/b/peierls; 

        vmag = fabs(vin); 
        vsign = (real8)Sign(vin);

        fdrag = bt*vmag*param->burgMag;
        dfdragdv = bt*param->burgMag; 

        fmag = peierls * fdrag; 
        dfmagdv = peierls*dfdragdv; 

        *fout = fmag * vsign;
        *dfdv = dfmagdv;
        *dfdvlin = *dfdv; 
        *d2fdv2 = 0.0;

        return;
}


static void EdgeDrag(Param_t *param, real8 vel[3], real8 burg[3],
                     real8 climbdir[3],  real8 fedrag[3],
                     real8 dfedragdv[3][3])
{
        int   m, n;
        real8 vmag, fout, dfdv, dfdvlin, d2fdv2;
        real8 Beclimb, Beline;
        real8 ftest;
        real8 linedir[3], glidedir[3], temp31[3];
        real8 tmp1[3][3], tmp2[3][3], tmp3[3][3];


        glidedir[0] = burg[0];
        glidedir[1] = burg[1];
        glidedir[2] = burg[2];

        Normalize(&glidedir[0], &glidedir[1], &glidedir[2]);

        cross(glidedir, climbdir, linedir);
    
        vmag = DotProduct(vel, glidedir);

        ftest = DotProduct(fedrag, glidedir);

        FscaleEdge(param, ftest, vmag, burg, &fout, &dfdv, &dfdvlin, &d2fdv2);

        Beclimb = dfdvlin * 1.0e+03;
        Beline = dfdvlin * 1.0e-03;

        Vec3TransposeAndMult(glidedir, tmp1);
        Vec3TransposeAndMult(climbdir, tmp2);
        Vec3TransposeAndMult(linedir,  tmp3);
        
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] = Beclimb * tmp2[m][n] +
                                  Beline  * tmp3[m][n]; 
            }
        }

        Matrix33Vector3Multiply(dfedragdv, vel, temp31);

        for (m = 0; m < 3; m++) {
            fedrag[m] = temp31[m] + fout * glidedir[m];
        }

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfedragdv[m][n] += dfdv * tmp1[m][n];
            }
        }

        return;
}


static void FscaleScrew(Param_t *param, real8 fsinit, real8 vin, real8 burg[3],
                   real8 *fout, real8 *dfdv, real8 *dfdvlin, real8 *d2fdv2)
{
        real8  alpha, beta, n, ninv;
        real8  burgnormal, temp, vmag, vsign;
        real8  dfdragdv, dfthermdv, dfmagdfdrag, dfmagdftherm, dfmagdv;
        real8  fmag, fdrag, ftherm;
        real8  ratio;
        real8  b, v0; 
        real8  taus, c0,c0p, vs;
        real8  factor110; 
        real8  bzero; 
        real8  test1; 


        factor110 = sqrt(3)/2.0; 
        burgnormal = sqrt(burg[0]*burg[0] + burg[1]*burg[1] + burg[2]*burg[2]);
        temp = param->TempK;
        taus = PEIERLS_SCREW*(1.0-1.e-3*temp)*(1.0-1.e-3*temp)*burgnormal;
        v0 = 1.33318333e-10*temp*temp*temp - 1.4985e-8*temp*temp -
             2.3379833333e-6*temp+2.5045e-4;

        vs = 0.0; 
        c0 = 3710.0/sqrt(temp);
        c0p = c0/param->burgMag;

        alpha       = ALPHA_SCREW;
        beta        = BETA_SCREW*temp;
        b           = B_SCREW_PRECONST;
        b = b/1e6; 
        n           = N_SCREW;
        ninv = 1.0 / n;

/*
 *      bzero = 1.68 at 300K.  It's the new coefficient that simply sets
 *      tau=0 at v=0. This should replace bt/b
 */ 
        bzero = 1.61*1e6; 

        vmag = fabs(vin); 
        vmag = vmag*factor110;
        vsign = (real8)Sign(vin);

        test1 = (vmag < 1);
        ratio = vmag / c0p; 

        if (test1 == 1) {
/* 
 *          below a very small velocity value, same as taus*dfthermdv but
 *          set the vmag to be zero
 *
 *          note: dfdvlin = taus * dfdragdv
 */ 
            *fout = 0.0e0;
            *dfdv = taus * alpha * beta /c0p * pow(v0,beta-1.0);
            *dfdvlin = bzero*param->burgMag;  
            *d2fdv2 = 0.0; 
            return;  
        } else {
            ftherm = alpha * (pow(ratio+v0, beta) - pow(v0, beta)); 
            dfthermdv = alpha * beta / c0p * pow(ratio+v0, beta-1.0); 
            fdrag = bzero*(vmag*param->burgMag + vs)/taus; 
            dfdragdv = bzero*param->burgMag/taus;
        }

        fmag = taus * pow(pow(fdrag, n)+pow(ftherm, n), ninv);
        dfmagdfdrag = taus * pow(fdrag, n-1.0) *
                      pow(pow(fdrag, n)+pow(ftherm, n), ninv-1.0);
        dfmagdftherm = taus * pow(ftherm, n-1.0) *
                       pow(pow(fdrag, n)+pow(ftherm, n), ninv-1.0);
        dfmagdv = dfmagdfdrag * dfdragdv + dfmagdftherm * dfthermdv;

        *fout = fmag * vsign;
        *dfdv = dfmagdv;
        *dfdvlin = bzero*param->burgMag;
        *d2fdv2 = 0.0; 

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
        real8 ftest, fproj[3];
        real8 dfdvlin;


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

        FscaleScrew(param, ftest, vMag, burg, &fout, &dfdv, &dfdvlin, &d2fdv2);

        Bsline = dfdvlin * 1.0e-03;

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

        return;
}


/**************************************************************************
 *
 *      Function:     Mobility_BCC_Fe_nl
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node          Pointer to the node for which to
 *                             calculate velocity
 *          IN/OUT: mobArgs    Structure containing additional
 *                             parameters for conveying information
 *                             to/from the mobility function.
 *
 *      Returns: 0 if successful
 *               1 if function could not converge on a velocity
 *
 *************************************************************************/
int Mobility_BCC_Fe_nl(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int     i, j, m, n, iAll, numNbrs, numscrews, thisScrew;
        int     numjunct, thisjunct;
        int     iterCnt, maxIter, maxIter2, count;
        int     numedges, thisEdge;
        int     zeroLenSeg[MAX_NBRS];
        int     edgeExists[MAX_NBRS], screwExists[MAX_NBRS];
        int     notConverged;
        int     numNonZeroLenSegs = 0;
        real8   velerrortol, velerror, minvelerror;
        real8   BJline, fscaleError, dfdvscale;
        real8   mag1, mag2;
        real8   eps, rt2, burgProd, mag, correctionAdjustment;
        real8   vtest[3], fn[3], rt[3], zero3[3], tmp3[3], ferror[3];
        real8   dfsdv[3][3], dferrordv[3][3];
        real8   dflindv[3][3];
        real8   dferrordvold[3][3], dfedv[3][3]; 
        real8   LTotal  [MAX_NBRS  ]   = {0.0} , LScrew    [MAX_NBRS  ]   = {0.0} , LEdge[MAX_NBRS]   = {0.0} , L2;
        real8   linedir [MAX_NBRS  ][3]={{0.0}}, ndir      [MAX_NBRS  ][3]={{0.0}}, burg [MAX_NBRS][3]={{0.0}};
        real8   screwb  [MAX_NSCREW][3]={{0.0}}, dfndfscrew[MAX_NSCREW]   = {0.0} ;
        real8   fscrew  [MAX_NSCREW][3]={{0.0}}, edgeb     [MAX_NEDGE ][3]={{0.0}};
        real8   junctb  [MAX_NBRS  ][3]={{0.0}}, junctdir  [MAX_NBRS  ][3]={{0.0}}, dfndfjunct[MAX_NBRS]={0.0};
        real8   edgenorm[MAX_NEDGE ][3]={{0.0}}, dfndfedge [MAX_NEDGE ]   = {0.0} ;
        real8   fedge   [MAX_NEDGE ][3]={{0.0}};
        real8   tmp;
        real8   tmp33[3][3];
        real8   correction[3];
        real8   vtmp[3] = { 0.0, 0.0, 0.0 };
        real8   mu, lmax;
        real8   forceerror, ferrortol;
        real8   ferrorold[3], vtestold[3], vtestnew[3], linmin;
        real8   totLen, massMult=0.0, massMatrix[3][3], vOld[3], vDiff[3], tmp3B[3];
        real8   EpsC0;
        Param_t *param;
        Node_t  *nbrNode;


        real8 (*invdferrordv)[3]     =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdferrordv);

/*
 *      This mobility module imposes no glide constraints on the node.
 */
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

        eps = 1.0e-12;
        iterCnt = 0;
        maxIter = 100;
        maxIter2 = 200;

        param = home->param;
        mu = param->shearModulus;
        lmax = param->maxSeg;

        EpsC0  = 1.0e-06 * C0;
        BJline = 1.e-3; 

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

        numNbrs = node->numNbrs;

/*
 *      For now, zero the velocity for nodes with more than
 *      3 connections...  a bit of a kludge!
 */
        if (numNbrs > 3) {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      This function saves various pieces of info for each neighbor
 *      node using static arrays for speed rather than doing memory
 *      allocations.  If the arrays are too  small, abort with an error
 *      to let the caller increase the array sizes.
 */
        if (numNbrs > MAX_NBRS) {
            Fatal("Mobility_BCC_Ta_nl: Node segment count (%d) exceeds\n"
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

/*
 *      Loop over all segments (skipping zero length segments)
 */
        numscrews = 0;
        numedges = 0;
        numjunct = 0;

        for (i = 0, iAll = 0; iAll < numNbrs; i++, iAll++) {

            if (zeroLenSeg[iAll]) {
                i--;
                continue;
            }

            burgProd = burg[i][0] * burg[i][1] * burg[i][2];

            if (fabs(burgProd) <= eps) {

                thisjunct = -1;

                for (j = 0; j < numjunct; j++) {
                    mag1 = ((burg[i][0] * junctb[j][0] +
                             burg[i][1] * junctb[j][1] +
                             burg[i][2] * junctb[j][2]) /
                            (Normal(burg[i]) * Normal(junctb[j])));
                    mag2 = ((linedir[i][0] * junctdir[j][0] +
                             linedir[i][1] * junctdir[j][1] +
                             linedir[i][2] * junctdir[j][2]) /
                            (Normal(linedir[i]) * Normal(junctdir[j])));
                    if ((fabs(fabs(mag1)-1.0) < eps) &&
                        (fabs(fabs(mag2)-1.0) < eps)) {
                        thisjunct = j;
                    }
                }

                if (thisjunct < 0) {
                    junctb[numjunct][0] = burg[i][0];
                    junctb[numjunct][1] = burg[i][1];
                    junctb[numjunct][2] = burg[i][2];
                    junctdir[numjunct][0] = linedir[i][0];
                    junctdir[numjunct][1] = linedir[i][1];
                    junctdir[numjunct][2] = linedir[i][2];
                    dfndfjunct[numjunct]  = 0.5 * LTotal[i];
                    numjunct++;
                } else {
                    dfndfjunct[thisjunct] += 0.5 * LTotal[i];
                }
            }

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

        ferrorold[0] = 0.0;
        ferrorold[1] = 0.0;
        ferrorold[2] = 0.0;

        vtestold[0] = 0.0;
        vtestold[1] = 0.0;
        vtestold[2] = 0.0;

        if (numjunct > 1) {
            vtest[0] = 0.0;
            vtest[1] = 0.0;
            vtest[2] = 0.0;
            notConverged = 0;
        } else {
            if (numjunct == 1) {
                tmp = DotProduct(vtest, junctdir[0]);
                vtest[0] = junctdir[0][0] * tmp;
                vtest[1] = junctdir[0][1] * tmp;
                vtest[2] = junctdir[0][2] * tmp;
            }
        }

        while (notConverged) {

            count++;

            correction[0] = 0.0;
            correction[1] = 0.0;
            correction[2] = 0.0;

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
                ferror[i] = fn[i] - tmp3B[i];
            }

            for (m = 0; m < 3; m++) {
                for (n = 0; n < 3; n++) {
                    dferrordv[m][n] = -massMatrix[m][n];
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

            for (i = 0; i < numjunct; i++) {
                Vec3TransposeAndMult(junctdir[i], tmp33);
                for (m = 0; m < 3; m++) {
                    ferror[m] -= BJline * vtest[m] * dfndfjunct[i];
                    for (n = 0; n < 3; n++) {
                        dferrordv[m][n] -= (BJline*dfndfjunct[i]) * tmp33[m][n];
                    }
                }
            }

/*
 *          Look for oscillatory behavior
 */
            if (numjunct == 1) {
                fscaleError = DotProduct(ferror, junctdir[0]);
                forceerror = fscaleError * fscaleError;
            } else {
                forceerror = DotProduct(ferror, ferror);
            }

            if (forceerror > ferrortol) {
                tmp3[0] = ferror[0] - ferrorold[0];
                tmp3[1] = ferror[1] - ferrorold[1];
                tmp3[2] = ferror[2] - ferrorold[2];

                if (DotProduct(tmp3, ferrorold) != 0.0) {
                    linmin = -(DotProduct(tmp3, ferrorold) /
                               DotProduct(tmp3, tmp3));
                } else {
                    linmin = 0.0;
                }
                if (fabs(linmin - 0.5) < 1.0e-5) {
                    for (m = 0; m < 3; m++) {
                        vtestnew[m] = linmin * vtest[m] +
                                      (1.0 - linmin) * vtestold[m];
                        tmp3[m] = vtestnew[m] - vtest[m];
                    }
                    velerror = DotProduct(tmp3, tmp3);
                } else {
                    if (numjunct == 1) {
                        Matrix33Vector3Multiply(dferrordv, junctdir[0], tmp3);
                        dfdvscale = DotProduct(junctdir[0], tmp3);
                        for (m = 0; m < 3; m++) {
                            correction[m] = (fscaleError / dfdvscale) *
                                            junctdir[0][m];
                        }
                        velerror = DotProduct(correction, correction);
                    } else {
                        forceerror = DotProduct(ferror, ferror);
                        if ( M33_INVERSE(invdferrordv,dferrordv) < 0 )
                        {
                            printf("%s::%s(%d) : Cannot invert matrix dferrordv!\n", __FILE__, __func__, __LINE__ );
                            printf("iteration = %d\n", count);
                            printf("numscrews = %d\n", numscrews);
                            printf("dferrordv = \n    ");

                            for (m = 0; m < 3; m++) {
                                for (n = 0; n < 3; n++) {
                                    printf("%.15e ", dferrordv[m][n]);
                                }
                                printf("\n    ");
                            }
                            printf("\n");
    
                            DumpNodeAndParams(home, node);

                            return(1);
                        }

                        Matrix33Vector3Multiply(invdferrordv, ferror,
                                                correction);
                        velerror = DotProduct(correction, correction);
                    }

                    for (i = 0; i < 3; i++) {
                        correction[i] *= correctionAdjustment;
                        vtestnew[i] = vtest[i] - correction[i];
                    }
                }

                for (m = 0; m < 3; m++) {
                    ferrorold[m] = ferror[m];
                    vtestold[m] = vtest[m];
                    vtest[m] = vtestnew[m];
                    for (n = 0; n < 3; n++) {
                        dferrordvold[m][n] = dferrordv[m][n];
                    }
                }
            } else {
                notConverged = 0;
                break;
            }

/*
 *          If the initial velocity guess was too far off, it's
 *          better to just zero the vtest value.  So, after the
 *          first iteration, if the force error is too large (an
 *          indication the initial velocity guess was way off),
 *          just zero the vtest value.
 */
            if ((count == 1) && (forceerror > 1.0e+02 * ferrortol)) {
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
                    correctionAdjustment = 0.50;
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

                return(1);
            }

        }  /* while (notConverged) */

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

#ifdef NAN_CHECK
        if ((isnan(node->vX) || isinf(node->vX)) ||
            (isnan(node->vY) || isinf(node->vY)) ||
            (isnan(node->vZ) || isinf(node->vZ))) {
            Fatal("(%d, %d) velocity: %lf %lf %lf",
                  node->myTag.domainID, node->myTag.index,
                  node->vX, node->vY, node->vZ);
        }
#endif

        return(0);
}
