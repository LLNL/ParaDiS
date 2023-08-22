/**************************************************************************
 *
 *      Module:  Functions common to the Chebyshev and the uniform grid methods.
 *
 *************************************************************************/

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include <math.h>
#include "Home.h"

#if defined (BBFMM) || defined (UNIFORMFMM)
#include "FM.h"

#if defined ANISOTROPIC
#include "Anisotropic.h"
#endif

#define MFREE(a)  if (a) { free(a); a=0; }

//---------------------------------------------------------------------------
// The IBM compiler continues to manifest errors when casting variable
// length and multidimensional matrices.  The ptrcpy() macro copies pointers
// via memcpy() to avoid the casting problems.
//---------------------------------------------------------------------------

#define ptrcpy(a,b)  memcpy(&(a),&(b),sizeof(void *))

//---------------------------------------------------------------------------

static int   ctab_exporder = 0;
static int   ctab_levels   = 0;
static real8 ctab_MU       = 0.0;
static real8 ctab_NU       = 0.0;

// summation using Clenshaw's recurrence relation
static double ClenshawSum(int n, double x, double *Tk)
{
    int j;
    double d0, d1, d2;
    d0 = d1 = 0;
    for (j = n - 1; j > 0; j--)
    {
        d2 = d0;
        d0 = 2.0 * x * d0 - d1 + Tk[j];
        d1 = d2;
    }
    return x * d0 - d1 + 0.5 * Tk[0];
}

void ComputeSn2(double *point[3], double *Tkmat, int n, int N, double *Sn[3])
{
    int i, k, m;
    double pfac = 2./(1.0*n);
    double *Tk;

    // Loop over Chebyshev node m
    for (m=0; m<n; m++)
    {
        k  = m*N;
        Tk = Tkmat + m*n;
        for (i=0; i<N; i++)
        {
            // Compute S_n for each direction using Clenshaw
            Sn[0][k+i] = pfac * ClenshawSum(n, point[0][i], Tk);
            Sn[1][k+i] = pfac * ClenshawSum(n, point[1][i], Tk);
            Sn[2][k+i] = pfac * ClenshawSum(n, point[2][i], Tk);
        }
    }
}

/*
 * Function: ComputeSn
 * ------------------------------------------------------------------
 * Computes S_n(x_m,x_i) for all Chebychev node-point pairs using
 * Clenshaw's recurrence relation.
 */
void ComputeSn(double *point[3], double *Tkz, int n, int N, double *Sn[3])
{
    int i, j, k, m;
    double vec[n], d[n+2], x;

    for (m=0; m<n; m++)
    {
        // Extract T_k for the Chebychev node x_m
        k = m*n;
        for (j=0; j<n; j++)
            vec[j] = Tkz[k+j];

        // Compute S_n for each direction independently using Clenshaw
        k = m*N;
        for (i=0; i<N; i++)
        {
            x = point[0][i];
            d[n] = d[n+1] = 0.;
            for (j=n-1; j>0; j--)
                d[j] = 2.0*x*d[j+1] - d[j+2] + vec[j];
            Sn[0][k+i] = x*d[1] - d[2] + 0.5*vec[0];

            x = point[1][i];
            d[n] = d[n+1] = 0.;
            for (j=n-1; j>0; j--)
                d[j] = 2.0*x*d[j+1] - d[j+2] + vec[j];
            Sn[1][k+i] = x*d[1] - d[2] + 0.5*vec[0];

            x = point[2][i];
            d[n] = d[n+1] = 0.;
            for (j=n-1; j>0; j--)
                d[j] = 2.0*x*d[j+1] - d[j+2] + vec[j];
            Sn[2][k+i] = x*d[1] - d[2] + 0.5*vec[0];
        }
    }
}

#ifdef ANISOTROPIC
void DerivativeOfGreensFunction
(
    Home_t *home,
    real8 x, real8 y, real8 z,
    real8 px, real8 py, real8 pz,
    int qMax, real8 dGdx[3][6]
)
{
    real8     te3 [2*qMax+3];
    complex8  te12[2*qMax+3];

    // Need some stuff from 'home'

    AnisotropicVars_t *avars = &home->anisoVars;

    real8    *fqr = (real8 *) avars->FqReal_v3;
    real8    *fqi = (real8 *) avars->FqImag_v3;
    real8   (*FqReal)[(qMax+2)*(qMax+1)]=0;  ptrcpy(FqReal,fqr);
    real8   (*FqImag)[(qMax+2)*(qMax+1)]=0;  ptrcpy(FqImag,fqi);

    complex8 *rotMatrix12   =   avars->anisoRotMatrix12;
    real8    *rotMatrixRow3 = &(avars->anisoRotMatrix[2][0]);

    // Define line direction

    real8 dX = px - x;
    real8 dY = py - y;
    real8 dZ = pz - z;

    real8 temp = dX*dX + dY*dY + dZ*dZ;
    real8 invL = 1.0 / sqrt(temp);

    real8 t[3] = { dX * invL,
                   dY * invL,
                   dZ * invL
                 };

    real8 coeff = INV_4_PI_SQUARED*invL*invL;

    te12[0] = 1.0;
    te12[1] = DotProduct(t, rotMatrix12);

    te3[0]  = 1.0;
    te3[1]  = DotProduct(t, rotMatrixRow3);

    for (int i = 2; i < (2*qMax+3); i++)
    {
        te12[i] = te12[i-1] * te12[1];
        te3 [i] = te3 [i-1] * te3 [1];
    }

    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 6; j++)
            dGdx[i][j] = 0.0;

    int local=0;
    for (int ind = 0; ind < qMax+1; ind++)
        for (int ind1 = 0; ind1 < 2*ind+2; ind1++)
        {
            complex8 J = te12[ind1]*te3[2*ind+1-ind1];
            for (int i = 0; i < 3; i++)
            {
                for (int j = 0; j < 6; j++)
                {
                    dGdx[i][j] += (creal(J)*FqReal[6*i+j][local] -
                                   cimag(J)*FqImag[6*i+j][local])*coeff;
                }
            }
            local++;
        }
}
#else
static void IsotropicKernel
(
    real8 x, real8 y, real8 z,
    real8 px, real8 py, real8 pz,
    real8 a, real8 mu, real8 nu,
    real8 *K
)
{
    int i,j,k;
    real8 RRx, RRy, RRz;
    real8 RR,RR3,RR5;
    real8 OneOverRR3, OneOverRR5;
    real8 Fac, Fac2, TwoFac;
    real8 RxOverRR3, RyOverRR3, RzOverRR3;
    real8 d3Rdx3, d3Rdy3, d3Rdz3;
    real8 d3Rdxdydz, d3Rdxdy2, d3Rdydx2, d3Rdzdx2, d3Rdxdz2, d3Rdydz2, d3Rdzdy2;
    real8 ThreeRx, ThreeRy, ThreeRz;
    real8 d3Ri[3], R[3];

    R[0] = x - px;
    R[1] = y - py;
    R[2] = z - pz;

    RRx  = R[0]*R[0];
    RRy  = R[1]*R[1];
    RRz  = R[2]*R[2];

    ThreeRx = 3*R[0];
    ThreeRy = 3*R[1];
    ThreeRz = 3*R[2];

    RR  = sqrt(RRx + RRy + RRz + a*a);
    RR3 = RR *RR*RR;
    RR5 = RR3*RR*RR;

    OneOverRR3 = 1/RR3;
    OneOverRR5 = 1/RR5;

    RxOverRR3 = R[0]*OneOverRR3;
    RyOverRR3 = R[1]*OneOverRR3;
    RzOverRR3 = R[2]*OneOverRR3;

    d3Rdx3 = -3*RxOverRR3 + ThreeRx*RRx*OneOverRR5;
    d3Rdy3 = -3*RyOverRR3 + ThreeRy*RRy*OneOverRR5;
    d3Rdz3 = -3*RzOverRR3 + ThreeRz*RRz*OneOverRR5;

    d3Rdxdy2 = -RxOverRR3 + ThreeRx*RRy*OneOverRR5;
    d3Rdydx2 = -RyOverRR3 + ThreeRy*RRx*OneOverRR5;
    d3Rdzdx2 = -RzOverRR3 + ThreeRz*RRx*OneOverRR5;

    d3Rdxdz2 = -RxOverRR3 + ThreeRx*RRz*OneOverRR5;
    d3Rdydz2 = -RyOverRR3 + ThreeRy*RRz*OneOverRR5;
    d3Rdzdy2 = -RzOverRR3 + ThreeRz*RRy*OneOverRR5;

    d3Rdxdydz = 3*R[0]*R[1]*R[2]*OneOverRR5;

    Fac  = mu/8.0/M_PI;
    Fac2 = Fac*2.0/(1-nu);
    TwoFac = 2*Fac;

    for (i=0; i<54; i++) K[i] = 0.0;

    d3Ri[0] = d3Rdx3+d3Rdxdy2+d3Rdxdz2;
    d3Ri[1] = d3Rdydx2+d3Rdy3+d3Rdydz2;
    d3Ri[2] = d3Rdzdx2+d3Rdzdy2+d3Rdz3;

    int s1=0,s2=1,s3=2,s4=3,s5=4,s6=5;
    int XX = 0, YY=1, ZZ=2;
    //int X=0, Y=1, Z=2;

    K[s6+6*YY+18*YY]=      -Fac * d3Ri[Z];
    K[s5+6*ZZ+18*ZZ]=       Fac * d3Ri[Y];
    K[s4+6*ZZ+18*ZZ]=      -Fac * d3Ri[X];
    K[s6+6*XX+18*XX]=       Fac * d3Ri[Z];
    K[s5+6*XX+18*XX]=      -Fac * d3Ri[Y];
    K[s4+6*YY+18*YY]=       Fac * d3Ri[X];

    K[s1+6*YY+18*XX]= - 2 * Fac * d3Ri[Z] + Fac2 * (d3Rdzdy2+d3Rdz3);
    K[s1+6*ZZ+18*XX]=   2 * Fac * d3Ri[Y] - Fac2 * (d3Rdy3+d3Rdydz2);
    K[s1+6*XX+18*YY]=                     - Fac2 * (d3Rdzdy2+d3Rdz3);
    K[s1+6*ZZ+18*YY]=                       Fac2 * (d3Rdxdy2+d3Rdxdz2);
    K[s1+6*XX+18*ZZ]=                       Fac2 * (d3Rdy3+d3Rdydz2);
    K[s1+6*YY+18*ZZ]=                     - Fac2 * (d3Rdxdy2+d3Rdxdz2);
    K[s6+6*YY+18*XX]=                     - Fac2 * d3Rdxdydz;
    K[s6+6*ZZ+18*XX]=      -Fac * d3Ri[X] + Fac2 * d3Rdxdy2;
    K[s6+6*XX+18*YY]=                       Fac2 * d3Rdxdydz;
    K[s6+6*ZZ+18*YY]=       Fac * d3Ri[Y] - Fac2 * d3Rdydx2;
    K[s6+6*XX+18*ZZ]=                     - Fac2 * d3Rdxdy2;
    K[s6+6*YY+18*ZZ]=                       Fac2 * d3Rdydx2;
    K[s5+6*YY+18*XX]=       Fac * d3Ri[X] - Fac2 * d3Rdxdz2;
    K[s5+6*ZZ+18*XX]=                       Fac2 * d3Rdxdydz;
    K[s5+6*XX+18*YY]=                       Fac2 * d3Rdxdz2;
    K[s5+6*ZZ+18*YY]=                     - Fac2 * d3Rdzdx2;
    K[s5+6*XX+18*ZZ]=                     - Fac2 * d3Rdxdydz;
    K[s5+6*YY+18*ZZ]=     - Fac * d3Ri[Z] + Fac2 * d3Rdzdx2;

    K[s2+6*YY+18*XX]=                       Fac2 * (d3Rdzdx2+d3Rdz3);
    K[s2+6*ZZ+18*XX]=                     - Fac2 * (d3Rdydx2+d3Rdydz2);
    K[s2+6*XX+18*YY]=   2 * Fac * d3Ri[Z] - Fac2 * (d3Rdzdx2+d3Rdz3);
    K[s2+6*ZZ+18*YY]= - 2 * Fac * d3Ri[X] + Fac2 * (d3Rdx3+d3Rdxdz2);
    K[s2+6*XX+18*ZZ]=                       Fac2 * (d3Rdydx2+d3Rdydz2);
    K[s2+6*YY+18*ZZ]=                     - Fac2 * (d3Rdx3+d3Rdxdz2);
    K[s4+6*YY+18*XX]=                     - Fac2 * d3Rdydz2;
    K[s4+6*ZZ+18*XX]=                       Fac2 * d3Rdzdy2;
    K[s4+6*XX+18*YY]=      -Fac * d3Ri[Y] + Fac2 * d3Rdydz2;
    K[s4+6*ZZ+18*YY]=                     - Fac2 * d3Rdxdydz;
    K[s4+6*XX+18*ZZ]=       Fac * d3Ri[Z] - Fac2 * d3Rdzdy2;
    K[s4+6*YY+18*ZZ]=                       Fac2 * d3Rdxdydz;

    K[s3+6*YY+18*XX]=                       Fac2 * (d3Rdzdx2+d3Rdzdy2);
    K[s3+6*ZZ+18*XX]=                     - Fac2 * (d3Rdydx2+d3Rdy3);
    K[s3+6*XX+18*YY]=                     - Fac2 * (d3Rdzdx2+d3Rdzdy2);
    K[s3+6*ZZ+18*YY]=                       Fac2 * (d3Rdx3+d3Rdxdy2);
    K[s3+6*XX+18*ZZ]= - 2 * Fac * d3Ri[Y] + Fac2 * (d3Rdydx2+d3Rdy3);
    K[s3+6*YY+18*ZZ]=   2 * Fac * d3Ri[X] - Fac2 * (d3Rdx3+d3Rdxdy2);
    return;
}
#endif

/*
 * Function: EvaluateKernel
 * -------------------------------------------------------------------
 * Evaluates the kernel given a source and a field point.
 * Anisotropic stress function implemented here.
 *          IN:     (x,y,z)      : A field point
 *          IN:     (px,py,pz)   : A point on the segment
 */
void EvaluateKernel(Home_t *home,
                    real8 x, real8 y, real8 z,
                    real8 px, real8 py, real8 pz, real8 *K)
{

    int i, j, k;
    Param_t *param;
    param= home->param;

#ifdef ANISOTROPIC
    real8 dGdx[3][6];
    real8 C633[6][3][3];

    int qMax= param->anisoHarmonicsNumTermsBase;;

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            C633[0][i][j] = home->anisoVars.elasticConstantMatrix4D[0][0][i][j];
            C633[1][i][j] = home->anisoVars.elasticConstantMatrix4D[1][1][i][j];
            C633[2][i][j] = home->anisoVars.elasticConstantMatrix4D[2][2][i][j];
            C633[3][i][j] = home->anisoVars.elasticConstantMatrix4D[1][2][i][j];
            C633[4][i][j] = home->anisoVars.elasticConstantMatrix4D[2][0][i][j];
            C633[5][i][j] = home->anisoVars.elasticConstantMatrix4D[1][0][i][j];
        }
    }

    DerivativeOfGreensFunction(home, x,y,z, px,py,pz, qMax,dGdx);

    int v,w,n,js;
    int alpha;
    real8 D[3][3][3];

    for (v=0; v<3; v++)
    {
        for (w=0; w<3; w++)
        {
            for (n=0; n<3; n++)
            {
                D[v][w][n] = 0.0;
                for (alpha=0; alpha<6; alpha++)
                {
                    D[v][w][n] += C633[alpha][w][n]*dGdx[v][alpha];
                }
            }
        }
    }


    // (j,s) indices for stress components and (w,r) indices for b_w and \xi_r components
    for (w=0; w<3; ++w)
    {
        for (js=0; js<6; ++js)
        {
            K[js+w*6     ] = 0.0;
            K[js+w*6+  18] = 0.0;
            K[js+w*6+2*18] = 0.0;

            for (v=0; v<3; v++)
            {
                K[js+w*6]      += C633[js][v][2] * D[v][w][1]-C633[js][v][1] * D[v][w][2];
                K[js+w*6+18]   += C633[js][v][0] * D[v][w][2]-C633[js][v][2] * D[v][w][0];
                K[js+w*6+2*18] += C633[js][v][1] * D[v][w][0]-C633[js][v][0] * D[v][w][1];
            }
        }
    }
#else
    real8 a  = param->rc;
    real8 mu = param->shearModulus;
    real8 nu = param->pois;

    IsotropicKernel(x,y,z, px,py,pz, a,mu,nu,K);
#endif
}

/*
 * Function below are for PBC corrections for the Chebyshev and the uniform grid
 * methods.
 */

static real8 *correctionTbl = NULL;


/*
 * Function: ComputeWeightsPBC
 * ------------------------------------------------------------------
 * Computes the weights for the Chebyshev nodes of all children cells
 * (identical for all cells and all levels so just compute once and
 * store in memory) (for PBC calculation each cell has 27 children
 * instead of 8).
 */
static void ComputeWeightsPBC(int n, real8 *Tkz,real8 *UpMat)
{
    int i, j, k, l, m, l1, l2, l3, count1, count2;
    int n3, n6, Nc;
    real8 fac, prefac, prefac3;
    real8 vtmp[3];

    real8 *nodes, *vec;
    real8 *Sxyz, *W;
    real8 *fieldt[3], *Sn[3];

    nodes = (real8 *)malloc(n * sizeof(real8));
    vec   = (real8 *)malloc(n * sizeof(real8));

    n3 = n*n*n;
    n6 = n3*n3;

    prefac = 2.0/(real8)n;
    prefac3 = prefac*prefac*prefac;
    Nc = 27*n3;                     // Number of child Chebyshev nodes

    for (i=0; i < 3; i++)
    {
        fieldt[i] = (real8 *)malloc(Nc * sizeof(real8)); // Chebyshev-transformed coordinates
        Sn[i]     = (real8 *)malloc(n * Nc * sizeof(real8));
    }

    Sxyz   = (real8 *)malloc(n3 * Nc * sizeof(real8));
    W      = (real8 *)malloc(n6 * sizeof(real8));

    fac=1.0/3.0;

    // Compute the n Chebyshev nodes of T_n(x)
    for (m=0; m<n; m++)
        nodes[m] = cos(M_PI*((real8)m+0.5)/(real8)n);

    // Map all Chebyshev nodes from the children cells to the parent
    k = 0;
    for (i=0; i<27; i++)
    {
        // Determine the mapping function for the specific child cell
        if (i<9)
        {
            vtmp[0] = -2;
            if      (i<3) vtmp[1] = -2;
            else if (i<6) vtmp[1] =  0;
            else          vtmp[1] =  2;
        }
        else if (i<18)
        {
            vtmp[0] = 0;
            if      (i<12) vtmp[1] = -2;
            else if (i<15) vtmp[1] =  0;
            else           vtmp[1] =  2;
        }
        else
        {
            vtmp[0] = 2;
            if      (i<21) vtmp[1] = -2;
            else if (i<24) vtmp[1] =  0;
            else           vtmp[1] =  2;
        }

        if      (i%3 == 0) vtmp[2] = -2;
        else if (i%3 == 1) vtmp[2] =  0;
        else               vtmp[2] =  2;

        for (l1=0; l1<n; l1++)
        {
            for (l2=0; l2<n; l2++)
            {
                for (l3=0; l3<n; l3++)
                {
                    fieldt[0][k] = fac*(nodes[l1] + vtmp[0]);
                    fieldt[1][k] = fac*(nodes[l2] + vtmp[1]);
                    fieldt[2][k] = fac*(nodes[l3] + vtmp[2]);
                    k++;
                }
            }
        }
    }

    // Compute Sc, the mapping function for the field points
    ComputeSn(fieldt,Tkz,n,Nc,Sn);

    // Compute Sxyz, the weights for the sources
    // Sxyz[Nc][n3]: Nc = 27n3
    // 27 L2L matrices n3 x n3 in column

    for (k=0; k<Nc; k++)
    {
        l = 0;
        for (l1=0; l1<n; l1++)
        {
            for (l2=0; l2<n; l2++)
            {
                for (l3=0; l3<n; l3++)
                {
                    Sxyz[l*Nc+k] = Sn[0][l1*Nc+k]*Sn[1][l2*Nc+k]*Sn[2][l3*Nc+k];
                    l++;
                }
            }
        }
    }

    // Accumulate the weights into a single weight matrix W
    for (i=0; i<n6; i++) { W[i] = 0; }

    for (i=0; i<27; i++)
    {
        count1 = 0;
        count2 = i*n3;
        for (j=0; j<n3; j++)
        {
            for (k=0; k<n3; k++)
            {
                W[count1] += Sxyz[k*Nc+count2]; // Notice the transpose here
                count1++;
            }
            count2++;
        }
    }

    count1 = 0;
    for (i=0; i<n3; i++)
        for (j=0; j<n3; j++, count1++)
            W[count1] *= prefac3;

    // Output the Upward Pass Matrix
    for (i=0; i<n6; i++) { UpMat[i] = W[i]; }

    MFREE(nodes);
    MFREE(vec);

    for (i=0; i<3; i++)
    {
        MFREE(fieldt[i]);
        MFREE(Sn    [i]);
    }
    MFREE(Sxyz);
    MFREE(W);
}


static void FORNOWSetElasticConstantMatrix4D
(
    real8 ecMatrix[6][6],
    real8 ecMatrix4D[3][3][3][3]
)
{
    int i, j, k, l;
    int i3to6[3][3] = { {0, 3, 4}, {3, 1, 5}, {4, 5, 2} };

    for (i = 0; i < 3; i++)
    {
        for (j = 0; j < 3; j++)
        {
            for (k = 0; k < 3; k++)
            {
                for (l = 0; l < 3; l++)
                {
                    ecMatrix4D[i][j][k][l] = ecMatrix[i3to6[i][j]][i3to6[k][l]];
                }
            }
        }
    }

    return;
}

/*
 * Read the FMM PBC Table
 */

void Stanford_CorrectionTableInit(Home_t *home)
{
    int i, n, n3, dofn3source, dofn3field, dof2n6;
    FILE    *fp;
    real8 c3[3];
    char    *fileName;
    Param_t *param;
    real8   eps = 1.0e-06;
    char    line[256];

    param    = home->param;
    fileName = param->fmCorrectionTbl;
    n        = param->fmExpansionOrder;

    // Just a sanity check to verify we have period boundaries
    // in at least one dimension.

    if (    (param->xBoundType != Periodic)
            && (param->yBoundType != Periodic)
            && (param->zBoundType != Periodic) )
    {
        return;
    }

    // Only the domain owning the single FM cell at the coarsest
    // layer has anything to do here.

    if (home->fmLayer[0].ownedCnt == 0) return;

    if (fileName[0] == 0)
    {
        Fatal("Error: The file containing the PBC image correction\n"
              "    table was not specified.  Please set the \n"
              "    <fmCorrectionTbl> value in the control file\n"
              "    to the name of the file containing the table.\n");
    }

    if ((fp = fopen(fileName, "r")) == (FILE *)NULL)
    {
        Fatal("Error %d opening PBC correction table file %s", errno, fileName);
    }

    // Do a few sanity checks before reading the whole FMM table

    printf("\nReading the FMM correction table : %s for PBC\n",fileName);

    // First, read some parameters used to set up the FMM table
#ifdef ANISOTROPIC
    real8 C[6][6];
    for (int i = 0; i < 6; i++)
    {
        for (int j = 0; j < 6; j++)
        {
            Getline(line, sizeof(line)-1, fp);
            sscanf(line,"%lf", &C[i][j]);
        }
    }
#else
    Getline(line, sizeof(line)-1, fp);
    sscanf(line,"%lf", &ctab_MU);

    Getline(line, sizeof(line)-1, fp);
    sscanf(line,"%lf", &ctab_NU);
#endif

    Getline(line, sizeof(line)-1, fp);
    sscanf(line,"%d", &ctab_exporder);

    Getline(line, sizeof(line)-1, fp);
    sscanf(line,"%d", &ctab_levels);

// Now, check that these parameters correspond to the one we have in param
#ifdef ANISOTROPIC
    if (!Matrix66_Near(param->elasticConstantMatrix, C, eps))
    {
        Fatal("%s: C matrix from table does not match C matrix\n"
              "for current simulation", __func__);
    }

#else
    real8 NU = param->pois;
    if ((NU < 0.0) || ((NU > 0.0) && fabs((NU-ctab_NU)/NU) > eps))
    {
        Fatal("Mismatch between configured value "
              "for pois\n    and value used to build FMM PBC correction "
              "table:\n    Configured pois = %e    Table pois = %e"
              "    Use ctablegen to create a new correction table",
              NU, ctab_NU);
    }
#endif

    if (param->fmExpansionOrder != ctab_exporder)
    {
        Fatal("Error: Current expansion order "
              "= %d, but\n    %s is built for "
              "expansion order = %d",
              param->fmExpansionOrder, fileName, ctab_exporder);
    }

    // Reading the actual FMM table coefficients

    n3 = n*n*n;
    dofn3source = n3 * 9;
    dofn3field  = n3 * 6;
    dof2n6 = dofn3source * dofn3field;

    correctionTbl = (real8 *) malloc( dof2n6 * sizeof(real8) );
    for (i=0; i<dof2n6; i++)
    {
        fscanf(fp, "%lf", correctionTbl+i);
    }
    fclose(fp);
    printf("\nCorrection table for FMM read\n");

}

void Stanford_FreeCorrectionTable(void)
{
    MFREE(correctionTbl);
}

void Stanford_CreateTableCorrection(Home_t *home, int n, int numLevels, int pbc[3])
{
    int i,j,k,m,lll;
    int n2, n3, n6, l, l1, l2, l3, count;
    int doff, dof2;
    int dofn3source, dofn3field, dof2n6;
    int incr=1;
    int ipbc, idof, isource, ifield;
    real8 scale;
    real8 sourcepos1[3], sourcepos2[3], origin[3];
    real8 alpha=1, beta=0;
    real8 boxl, cshift[3], vtmp[3];
    real8 *UpMat, *Ulevel, *Tkz, *nodes, *vec;
    real8 Kij1[54], Kij2[54], Kij3[54], Kij4[54];
    real8 *Ilevel, *MMresult, *KPBC, *KPBC_j;
    real8 *nodepos[3];
    char transa='n';
    Param_t *param;
    param = home->param;

#ifdef ANISOTROPIC
    // If anisotropic elasticity is on, need the coefficients
    // Compute them here
    AnisotropicInit(home);
#endif

    boxl = param->maxSideX - param->minSideX;

    n2 = n*n;
    n3 = n2*n;
    n6 = n3*n3;

    doff =  6;
    dof2 = 54;
    dofn3source = 9*n3;
    dofn3field  = 6*n3;
    dof2n6 = dofn3source * dofn3field;

    nodes = (double *)malloc(n * sizeof(double));
    vec   = (double *)malloc(n * sizeof(double));
    Tkz = (double *)malloc(n2 * sizeof(double));
    for (m=0; m<3; m++)
        nodepos[m] = (double *)malloc(n3 * sizeof(double));
    UpMat = (real8 *) malloc(n6 *sizeof(real8) );

    // Ulevel is U^i at the ith level,
    Ulevel   = (real8 *) calloc(1, n6*sizeof(real8) );
    MMresult = (real8 *) malloc(n6 * sizeof(real8) );

    // Allocate and initialize KPBC
    KPBC = (real8 *) calloc(1, dof2n6*sizeof(real8) );
    KPBC_j = (real8 *) malloc( dof2n6 * sizeof(real8) );

    // Compute the Chebyshev nodes of T_n(x)
    for (m=0; m<n; m++)
        nodes[m] = cos(M_PI*((real8)m+0.5)/(real8)n);

    // Compute the locations of the field points
    count = 0;
    for (l1=0; l1<n; l1++)
    {
        vtmp[0] = 0.5*nodes[l1];
        for (l2=0; l2<n; l2++)
        {
            vtmp[1] = 0.5*nodes[l2];
            for (l3=0; l3<n; l3++)
            {
                nodepos[0][count] = vtmp[0];
                nodepos[1][count] = vtmp[1];
                nodepos[2][count] = 0.5*nodes[l3];
                count++;
            }
        }
    }

    // Compute Tkz
    for (m=0; m<n; m++)
    {
        ComputeTk(nodes[m],n,vec);
        i = m*n;
        for (k=0; k<n; k++) Tkz[i+k] = vec[k];
    }

    // Compute the upward pass matrix
    ComputeWeightsPBC(n, Tkz, UpMat);

    // Set the initial scaling factor
    scale = 1.0;
    origin[0] = 0, origin[1] = 0, origin[2] = 0;

    // Initialize Ulevel to be identity matrix
    for (k=0; k<n3; k++)
        Ulevel[k*(n3+1)] = 1.;

    // Compute the field values due to periodic sources in different shells
    for (ipbc=1; ipbc<numLevels; ipbc++)
    {
        // Compute the jth column that M(j) = KPBC_j * U^i(j)
        for (j=0; j<n3; j++)
        {
            // Initialize KPBC_j
            for (i=0; i<dof2n6; i++) KPBC_j[i] = 0.;

            // Compute KPBC_j
            for (l=0; l<n3; l++)
            {

                for (l1=-4; l1<5; l1++)
                {
                    cshift[0] = (real8)l1;
                    for (l2=-4; l2<5; l2++)
                    {
                        cshift[1] = (real8)l2;
                        for (l3=-4; l3<5; l3++)
                        {
                            cshift[2] = (real8)l3;
                            if (abs(l1) > 1 || abs(l2) > 1 || abs(l3) > 1)
                            {

                                sourcepos1[0] = (cshift[0]+nodepos[0][l])*scale;
                                sourcepos1[1] = (cshift[1]+nodepos[1][l])*scale;
                                sourcepos1[2] = (cshift[2]+nodepos[2][l])*scale;

                                sourcepos2[0] = sourcepos1[0] - nodepos[0][j];
                                sourcepos2[1] = sourcepos1[1] - nodepos[1][j];
                                sourcepos2[2] = sourcepos1[2] - nodepos[2][j];

                                EvaluateKernel(home,     origin[0],     origin[1],     origin[2],
                                               sourcepos1[0], sourcepos1[1], sourcepos1[2], Kij1);

                                EvaluateKernel(home,     origin[0],     origin[1],     origin[2],
                                               sourcepos2[0], sourcepos2[1], sourcepos2[2], Kij2);

                                for (k=0; k<n3; k++)
                                {

                                    EvaluateKernel(home,    nodepos[0][k], nodepos[1][k], nodepos[2][k],
                                                   sourcepos1[0], sourcepos1[1], sourcepos1[2], Kij3);
                                    EvaluateKernel(home, nodepos[0][k], nodepos[1][k], nodepos[2][k],
                                                   sourcepos2[0], sourcepos2[1], sourcepos2[2], Kij4);

                                    for (idof=0; idof < dof2; idof++)
                                    {
                                        KPBC_j[idof*n6 + l*n3 + k] += Kij3[idof] + Kij2[idof] -
                                                                      Kij4[idof] - Kij1[idof];
                                    }
                                }
                            }
                        }
                    }
                }
            }

            // Compute M_j
            for (isource=0; isource < 9; isource++)
                for (ifield=0; ifield < 6; ifield++)
                    dgemv_(&transa, &n3, &n3, &alpha, KPBC_j+(isource*6 + ifield)*n6, &n3, Ulevel+j*n3,
                           &incr, &alpha, KPBC+j*dofn3field*9 + ifield + dofn3field*isource, &doff);

        }  // j


        // Update 'Ulevel'
        dgemm_(&transa, &transa, &n3, &n3, &n3, &alpha, Ulevel, &n3, UpMat, &n3, &beta, MMresult, &n3);

        for (i=0; i<n6; i++)
            Ulevel[i] = MMresult[i];

        // Increase scaling factor
        scale *= 3.0;
    } // ipbc

    MFREE(nodes);
    MFREE(vec  );
    MFREE(Tkz  );
    MFREE(UpMat);

    MFREE(Ulevel  );
    MFREE(MMresult);
    MFREE(KPBC_j  );

    /* Write PBC table */
    FILE *fp;
    fp = fopen(param->fmCorrectionTbl,"w");
    if (fp == (FILE *)NULL) Fatal("Error %d opening output file %s", errno, param->fmCorrectionTbl);


    // Write the elastic constants
#ifdef ANISOTROPIC
    /*
     *    When support for anisotropy is compiled in and enabled, we have to
     *    dump to the file the elastic constants matrix used when creating the
     *    correction table
     */
    for (i = 0; i < 6; i++)
    {
        for (j = 0; j < 6; j++)
        {
            fprintf(fp, "%-25.15lf   # C[%d][%d]\n",
                    param->elasticConstantMatrix[i][j], i, j);
        }
    }
#else
    double mu, nu;
    mu = param->shearModulus;
    nu = param->pois;

    fprintf(fp, "%-25.15lf   # MU\n", mu);
    fprintf(fp, "%-25.15lf   # NU\n", nu);
#endif

    fprintf(fp, "%d   # Expansion order\n", n);
    fprintf(fp, "%d   # # of levels of PBC images\n", numLevels);

    // Write pbc kernel
    for (i=0; i<dof2n6; i++)
        fprintf(fp, "%.15e\n", KPBC[i]);
    fclose(fp);

    MFREE(KPBC);
}


void Stanford_DoTableCorrection(Home_t *home, int n, real8 homogen)
{
    Param_t *param = home->param;

    // Only the domain owning the single FM cell at the highest(coarsest)
    // layer has anything to do here.

    FMLayer_t *layer = &home->fmLayer[0];
    if (layer->ownedCnt == 0) return;

    FMCell_t  *cell = LookupFMCell(layer->cellTable, 0);

    // Scale for the box size, since pbc kernel is created corresponding
    // to unit boxes. Len is the Length.
    real8 Len   = param->Lx;
    real8 scale = pow(1/Len, homogen);

    int  n3          = n*n*n;
    int  dofn3field  = n3 * 6;
    int  dofn3source = n3 * 9;
    int  dof2n6      = dofn3source * dofn3field;

    // Compute influence from PBC
    real8 beta   = 0;
    int   incr   = 1;
    char  transb = 'n';
    dgemv_(&transb, &dofn3field, &dofn3source, &scale, correctionTbl,
           &dofn3field, cell->mpCoeff, &incr, &beta,
           cell->expansionCoeff, &incr);
}


/*
 * Function: ComputeTk
 * ------------------------------------------------------------------
 * Computes T_k(x) for k between 0 and n-1 inclusive.
 */
void ComputeTk(real8 x, int n, real8 *vec)
{
    vec[0] = 1;
    vec[1] = x;

    for (int k=2; k<n; k++)
        vec[k] = 2.0*x*vec[k-1] - vec[k-2];
}

static  double LagrangeWeight(int n, double x, int m)
{
    double num   = 1.0;
    double denom = 1.0;
    for (int j=0; j<n; j++)
    {
        if(m!=j)
        {
            num   *= x - (1-2*(double)j/(double)(n-1));
            denom *= -2*((double)m-(double)j)/(double)(n-1);
        }
    }
    return num/denom;
}

/*
   in : expansion coeff
  out : mean stress
*/

static void Print3x32(const char *msg, real8 A[3][3])
{
    printf("%s\n", msg);
    printf("%.15e %.15e %.15e\n", A[0][0], A[0][1], A[0][2]);
    printf("%.15e %.15e %.15e\n", A[1][0], A[1][1], A[1][2]);
    printf("%.15e %.15e %.15e\n\n", A[2][0], A[2][1], A[2][2]);
}

static void LocalCorrectionStress
(
    Home_t *home, int n,  real8 *Tkz,
    real8 MeanStress[3][3])
{
    int i,j,m;
    int l1, l2, l3, idof;
    int count, cellID;
    real8 meanWeight, Volume;
    real8 sigma[6], stress[3][3];

    Param_t   *param;
    FMLayer_t *layer;
    FMCell_t  *cell;

    param = home->param;

    Volume = 8.0;

    // Integrate Sn
    real8 SnInt[n];

    for (i=0; i<n; i++)
    {
        SnInt[i] = 1;
        for (j=2; j<n; j+=2)
        {
            // Int_{-1}^1 Tk(x) dx
            SnInt[i] -= 2.0/(j*j*1.0-1.0) * Tkz[i*n+j];
        }

        SnInt[i] *= 2.0/(n*1.0);
    }

    layer = &home->fmLayer[0];
    //if (layer->ownedCnt == 0) return;

    if (layer->ownedCnt != 0)
    {
        cell = LookupFMCell(layer->cellTable, 0);

        // Sum up the stress for all cells at the layer

        for (j = 0; j < 6; j++) sigma[j] = 0;

        for (l1=0; l1<n; l1++)
            for (l2=0; l2<n; l2++)
                for (l3=0; l3<n; l3++)
                {
                    count = (l1 * n*n + l2*n + l3) * 6;

                    meanWeight = SnInt[l1] * SnInt[l2] * SnInt[l3];

                    for (idof=0; idof< 6; idof++)
                    {
                        sigma[idof] += cell->expansionCoeff[count + idof] * meanWeight;
                    }
                }

        stress[0][0] = sigma[0];
        stress[1][1] = sigma[1];
        stress[2][2] = sigma[2];

        stress[0][1] = sigma[5];
        stress[1][2] = sigma[3];
        stress[0][2] = sigma[4];


        stress[1][0] = stress[0][1];
        stress[2][0] = stress[0][2];
        stress[2][1] = stress[1][2];
    }

#ifdef PARALLEL
    MPI_Allreduce(stress, MeanStress, 9, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

    for (j = 0; j < 3; j++)
    {
        for (m = 0; m < 3; m++)
        {
            MeanStress[j][m] /= Volume;
        }
    }

#else
    for (m = 0; m < 3; m++)
    {
        for (j = 0; j < 3; j++)
        {
            MeanStress[m][j] = stress[m][j]/Volume;
        }
    }
#endif

    return;
}


/*---------------------------------------------------------------------------
 *
 *      Author:      Chao Chen, Eric Darve and Sylvie Aubry
 *      Function:    bbfmm_MeanStressCorrrection
 *      Description: Subtract the average stress computed in BB_MeanL2P()
 *                   from the first order coefficients of the taylor
 *                   coefficients of each cell this domain intersects.
 *
 *-------------------------------------------------------------------------*/
void Stanford_MeanStressCorrection(Home_t *home, int n, real8 *Tkz)
{
    Param_t   *param =  home->param;
    FMLayer_t *layer = &home->fmLayer[param->fmNumLayers-1];

    real8 MeanStress[3][3];
    LocalCorrectionStress(home, n, Tkz, MeanStress);

    for (int i=0; i < layer->ownedCnt; i++)
    {
        int       cellID = layer->ownedMin + i;
        FMCell_t *cell   = LookupFMCell(layer->cellTable, cellID);

        for (int j=0; j<3; j++)
            for (int m=0; m<3; m++)
                cell->expansionCoeff[j*3+m] -= MeanStress[j][m];
    }
}

#endif // BBFMM or UNIFORMFMM
