/***************************************************************************
 *
 *      Module:      AnisotropicInitFuncs.c
 *
 *      Description: This module contains functions to initialize certain
 *                   tables needed for use with anisotropic elasticity.
 *
 *      Includes public functions:
 *          AnisotropicInit()
 *
 *      Includes private functions:
 *          CalcFqTable()
 *          CalcCoreFqTable()
 *          GetGreensFuncDerivatives()
 *          GetBmatrix()
 *          ComputeLegendre()
 *          CalcGLMandBinomialFactCoeffs()
 *          FindNonZeroFqComponents()
 *          SetElasticConstantMatrix()
 *          SetElasticConstantMatrix4D()
 *          SetCubicElasticModuli()
 *          SetHexagonalElasticModuli()
 *          AnisotropicInitElasticConstants()
 *          MultVecElastConstVec()
 *
 **************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

#ifdef ANISOTROPIC

#define GREENS_FUNC_NUM_POINTS 100
#define NUM_INTEGRATION_POINTS 100

void CalcCoreFqTable(int qMax, int numTheta, int numPhi,
                     real8 *phi, real8 *BgridIn,
                     real8 *gaussQuadPoints, real8 *gaussQuadWeights,
                     real8 *FqRealPtrIn, real8 *FqImagPtrIn);

//---------------------------------------------------------------------------
// The IBM compiler continues to manifest errors when casting variable
// length and multidimensional matrices.  The ptrcpy() macro copies pointers
// via memcpy() to avoid the casting problems.
//---------------------------------------------------------------------------

#define ptrcpy(a,b)  memcpy(&(a),&(b),sizeof(void *))

/*---------------------------------------------------------------------------
 *
 *      Function:        MultVecElastConstVec
 *      Description:     Calculates UV = U * C * V from provided vectors
 *                       U, V and the elastic constant matrix C.  In
 *                       other words: (UV_jk = U_i * Cijkl * V_l)
 *
 *                       Note: this function could be simplified if
 *                             C is cubic.
 *
 *      Parameters:
 *          IN:  U   First vector
 *          IN:  V   Second vector
 *          OUT: UV  Resultant matrix
 *          IN:  C   Elastic constants provided in the 3x3x3x3 matrix form
 *
 *--------------------------------------------------------------------------*/
static void MultVecElastConstVec(real8 U[3], real8 V[3], real8 UV[3][3],
                                 real8 C[3][3][3][3])
{
        int i, j, k, l;

        for (j = 0; j < 3; j++) {
            for (k = 0; k < 3; k++) {

                UV[j][k] = 0.0;

                for (i = 0; i < 3; i++) {
                    for (l = 0; l < 3; l++) {

                        UV[j][k] += U[i] * C[i][j][k][l] * V[l];
                    }
                }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    Computes B and dBdt
 *      Description: B and dBdt are elements of the line tension force and
 *                   its derivative.
 *
 *      Parameters:
 *          IN:     nMax
 *          IN:     C
 *          OUT:    dGdx
 *
 *-------------------------------------------------------------------------*/
void GetBmatrix(real8 t[3], real8 C[3][3][3][3], real8 B6[6])
{
        int   i, j, k;
        int   alpha;
        real8 omega, domega;
        real8 MaxNormL, NormLmati;
        real8 Lmat[3][3];
        real8 temp1[3][3], temp2[3][3];
        real8 n[3], m[3], M[3], N[3], MM[3][3], NN[3][3];
        real8 nn[3][3],nninv[3][3],nm[3][3],mn[3][3],mm[3][3];
        real8 Fac, Fac2;
        real8 cosomega, sinomega, B[3][3];

/*
 *      Initialization
 */
        omega = 0.0;

/*
 *      Define N and M (two vectors orthonormal to the line direction <t>)
 */
        for (i = 0; i < 3; i++) {

            for (j = 0; j < 3; j++) {
                Lmat[i][j] = (i==j) - t[i]*t[j];
            }
        }


        alpha = 0;
        MaxNormL = 0.0;

        for (i = 0; i < 3; i++) {

            NormLmati = sqrt(Lmat[0][i]*Lmat[0][i] +
                             Lmat[1][i]*Lmat[1][i] +
                             Lmat[2][i]*Lmat[2][i]);

            if (NormLmati > MaxNormL) {

                MaxNormL = NormLmati;
                alpha = i;
            }
        }

        N[0] = Lmat[0][alpha];
        N[1] = Lmat[1][alpha];
        N[2] = Lmat[2][alpha];

        NormalizeVec(N);
        cross(N,t,M);

/*
 *      Compute B as in Eq. (4.1.39) of Bacon, Barnett and
 *      Scattergood p. 149
 */
        MultVecElastConstVec(M,M,MM,C);
        MultVecElastConstVec(N,N,NN,C);

/*
 *      integration starts here: integrate 0 to pi, NOT 2*pi
 */
        domega = M_PI/(1.0*(NUM_INTEGRATION_POINTS));

        for (i = 0; i < 3; i++) {
            VECTOR_ZERO(B[i]);
        }

        for (i = 1; i <= NUM_INTEGRATION_POINTS; i++) {

            omega = i*domega;
            cosomega = cos(omega);
            sinomega = sin(omega);

/*
 *
 *          Rotating integration vectors n and m around N and M
 *          in a plane orthogonal to t
 */
            m[0]=cosomega*M[0] + sinomega*N[0];
            m[1]=cosomega*M[1] + sinomega*N[1];
            m[2]=cosomega*M[2] + sinomega*N[2];

            n[0]=cosomega*N[0] - sinomega*M[0];
            n[1]=cosomega*N[1] - sinomega*M[1];
            n[2]=cosomega*N[2] - sinomega*M[2];

/*
 *          Compute B
 */
            MultVecElastConstVec(n,m,nm,C);
            MultVecElastConstVec(m,n,mn,C);
            MultVecElastConstVec(n,n,nn,C);

            Matrix33_Inverse(nninv,nn);

            Matrix33_Mul(temp1,mn,nninv);   // temp1=mn*nninv
            Matrix33_Mul(temp2,temp1,nm);   // temp2=mn*nninv*nm

            Matrix33_Sub(B,B,temp2);        // B = B - (mn*nninv*nm)

        }  /* end loop over integration points */

/*
 *      Compute analyticaly \int_0^pi mm dt
 */
        MultVecElastConstVec(M, M, temp1, C);
        MultVecElastConstVec(N, N, temp2, C);

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                mm[i][j] = temp1[i][j] + temp2[i][j];
            }
        }


        Fac = domega/(4.0*M_PI*M_PI);
        Fac2= 1./(8.0*M_PI);  /* Ov8pi */

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                B[i][j] *= Fac;
                B[i][j] += mm[i][j] * Fac2;
	    }
	}


	/* change the 3x3 into 6 */
	B6[0] = B[0][0];
	B6[1] = B[1][1];
	B6[2] = B[2][2];
	B6[3] = B[1][2];
	B6[4] = B[2][0];
	B6[5] = B[0][1];
	
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    CalcFqTable()
 *
 *      Description: For a pair of dislocation segments, the force and
 *                   stress's spherical harmonics decomposition can be
 *                   written as the sum of the product of two
 *                   functions: Fq and H, where Fq can be precomputed
 *                   independently of the two dislocation segments. This
 *                   function computes the pre-factor Fq.
 *
 *      Parameters:
 *          IN:   qMax       Base factor for determining the number of terms
 *                           in the expansion for spherical harmonics.  The
 *                           total number of terms is 2 * qMax + 1.
 *          IN:   binFactTblPtr Pointer to the table of binomial factorial
 *                           coefficients needed for calculating the Fq table.
 *                           Internally, this table is cast to a 3 dimensional
 *                           array of size: [2*qMax+2][2*qMax+2][2*qMax+2]
 *          IN:   glmPtr     Pointer to expansion coefficient array. This
 *                           array is cast internally into a multi-dimensional
 *                           array of dimensions: [2*qMax+2][2*qMax+2][3][6]
 *          OUT:  FqTblReal  Pointer the the first element of the array
 *                           of real components of the Fqtable.  Dimensions
 *                           of the table are [6][3][(qMax+2)*(qMax+1)]
 *          OUT:  FqTblImag  Pointer the the first element of the array
 *                           of imaginary components of the Fqtable.
 *                           What is the Fq table? Dimensions of the
 *                           table are [6][3][(qMax+2)*(qMax+1)]?
 *                           This table will be updated by this function.
 *          OUT:  FqRealMax  Contents will be set to the absolute value
 *                           of the element in the FqReal array with the
 *                           largest absolute value.
 *          OUT:  FqImagMax  Contents will be set to the absolute value
 *                           of the element in the FqImag array with the
 *                           largest absolute value.
 *
 *-------------------------------------------------------------------------*/
static void CalcFqTable(int qMax, real8 *binFactTblPtr,
                        real8 *glmRePtr, real8 *glmImPtr,
                        real8 *FqRealPtr, real8 *FqImagPtr,
                        real8 *FqRealMax, real8 *FqImagMax)
{
        int   i, j, k, l, m, n;
        real8 fac;
        real8 (*binomialFactCoeffs)[2*qMax+2][2*qMax+2];  ptrcpy(binomialFactCoeffs,binFactTblPtr);    // binomialFactCoeffs = binFactTblPtr (pointer cast)
        real8 (*FqReal)[3][(qMax+2)*(qMax+1)];            ptrcpy(FqReal,FqRealPtr);                    // FqReal = FqRealPtr (pointer cast)
        real8 (*FqImag)[3][(qMax+2)*(qMax+1)];            ptrcpy(FqImag,FqImagPtr);                    // FqImag = FqImagPtr (pointer cast)
        real8 (*glmRe )[2*qMax+2][3][6];                  ptrcpy(glmRe ,glmRePtr );                    // glmRe  = glmRePtr  (pointer cast)
        real8 (*glmIm )[2*qMax+2][3][6];                  ptrcpy(glmIm ,glmImPtr );                    // glmIm  = glmImPtr  (pointer cast)

/*
 *      Initially zero the arrays for the real and imaginary portions of
 *      the Fq table.
 */
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 6; n++) {
                for (i = 0; i < (qMax+1)*(qMax+2); i++) {
                    FqReal[n][m][i] = 0.0;
                    FqImag[n][m][i] = 0.0;
                }
            }
        }

/*
 *      ???
 */
        for (m = 0; m < 3; m++) {
            for (n = 0; n < 6; n++) {
                for (i = 0; i <= qMax; i++) {

                    l = 2 * i + 1;

                    for (j = 0; j <= l; j++) {
                        int intlmj;

                        if (j == 0) {
                            fac = 1;
                        } else {
                            fac = 2;
                        }

                        if ((l-j) % 2 == 0) {
                            intlmj = (l-j) / 2;
                        } else {
                            intlmj = (l-j-1) / 2;
                        }

                        for (k = 0; k <= intlmj; k++) {
                            int   ind, ind1, local;

                            ind1 = j;
                            ind = i-k;
                            local = ind*(ind+1) + ind1;

                            FqReal[n][m][local] += fac *
                                                   binomialFactCoeffs[j][l][k] *
                                                   glmRe[l][j][m][n];

                            FqImag[n][m][local] += fac *
                                                   binomialFactCoeffs[j][l][k] *
                                                   glmIm[l][j][m][n];
                        }
                    }
                }
            }
        }

/*
 *      Find the real and imaginary components of the Fq table with
 *      the largest magnitude.
 */
        *FqRealMax = 0.0;
        *FqImagMax = 0.0;

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 6; n++) {
                for (i = 0; i < (qMax+2)*(qMax+1); i++) {

                    if (fabs(FqReal[n][m][i]) > *FqRealMax) {
                        *FqRealMax = fabs(FqReal[n][m][i]);
                    }

                    if (fabs(FqImag[n][m][i]) > *FqImagMax) {
                        *FqImagMax = fabs(FqImag[n][m][i]);
                    }
                }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    GetGreensFuncDerivatives()
 *      Description: The angular part of the derivative of the
 *                   Green's function is given by an integral.  This
 *                   function computes that integral and returns the
 *                   derivative.
 *
 *      Parameters:
 *          IN:  numPoints  Number of points used to approximate the integral
 *          IN:  lineDir    Direction vector for sperical harmonics.
 *          IN:  ecMatrix   Elastic constant matrix in it's
 *                          [3][3][3][3] form
 *          OUT: dGdx       Derivative values of the angular part of the
 *                          Green's function.
 *
 *-------------------------------------------------------------------------*/
static void GetGreensFuncDerivatives(int numPoints, real8 lineDir[3],
                                     real8 ecMatrix[3][3][3][3],
                                     real8 dGdx[3][3][3])
{
        int   i, j, k, l, m, n, p;
        int   dirMinIndex;
        real8 tmp1;
        real8 xHat[3], yHat[3];

/*
 *      Build the basis (lineDir, z, y)
 */
        NormalizeVec(lineDir);

        dirMinIndex = (fabs(lineDir[1]) < fabs(lineDir[0]));
        dirMinIndex = (fabs(lineDir[2]) < fabs(lineDir[dirMinIndex])) ?
                      2 : dirMinIndex;

        xHat[0] = 0.0;
        xHat[1] = 0.0;
        xHat[2] = 0.0;

        xHat[dirMinIndex] = 1.0;

        tmp1 = DotProduct(lineDir, xHat);

        for (i = 0; i < 3; i++) {
            xHat[i] -= tmp1 * lineDir[i];
        }

        NormalizeVec(xHat);
        cross(lineDir, xHat, yHat);

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    dGdx[i][j][k] = 0.0;
                }
            }
        }

        for (i = 0; i < numPoints; i++) {
            real8 psi, sinPsi, cosPsi;
            real8 z[3];
            real8 mZ[3][3];
            real8 fZ[3][3];
            real8 mStar[3][3];

            psi = M_PI * (i+1) / numPoints;

            sinPsi = sin(psi);
            cosPsi = cos(psi);

            for (j = 0; j < 3; j++) {
                z[j] = cosPsi * xHat[j] + sinPsi * yHat[j];
            }

/*
 *          mStar = (zz)^({-1}
 */
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    mZ[j][k] = 0.0;
                }
            }

            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    for (l = 0; l < 3; l++) {
                        for (m = 0; m < 3; m++) {
                            mZ[l][k] += ecMatrix[l][j][k][m] * z[j] * z[m];
                        }
                    }
                }
            }  /* end for (j = 0; j < 3; j++) */

            if ( Matrix33_Inverse(mStar,mZ) < 0) {
                printf("ERROR: Unabled to invert mZ\n");
                exit(0);
            }

            for (j = 0; j < 3; j++) {
                fZ[j][0] = 0.0;
                fZ[j][1] = 0.0;
                fZ[j][2] = 0.0;
            }

            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    for (l = 0; l < 3; l++) {
                        for (m = 0; m < 3; m++) {
                            for (n = 0; n < 3; n++) {
                                for (p = 0; p < 3; p++) {
                                    fZ[j][k] += ecMatrix[l][m][n][p] *
                                                mStar[j][l] *
                                                mStar[n][k] *
                                                (z[m]*lineDir[p] + z[p]*lineDir[m]);
                                }
                            }
                        }
                    }
                }
            }

/*
 *          Integrand
 */
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    for (l = 0; l < 3; l++) {
                        dGdx[j][k][l] = dGdx[j][k][l] -
                                        lineDir[l] * mStar[j][k] +
                                        z[l] * fZ[j][k];
                    }
                }
            }

        }  /* end for (i = 0; i < numPoints; i++) */

        tmp1 = M_PI / (real8)numPoints;

        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {
                    dGdx[i][j][k] *= tmp1;
                }
            }
        }

        return;
}


/*
 *      For use when calculating the Legendre polynomials; to
 *      access P_l,m in an array.
 *
 *      {0,0}{1,0}{1,1}{2,0}{2,1}{2,2} ...
 *
 *      summation series over l + m => (l*(l+1))/2 + m
 */
#define ATLM(l,m) ( (((l)*((l)+1))>>1) + (m) )

/*---------------------------------------------------------------------------
 *
 *      Function:    ComputeLegendre
 *      Description: Compute the legendre polynomial by recurrence.
 *
 *      Parameters:
 *          IN:  P         Degree of the polynomial
 *          IN:  x         Value for which the polynomial is being calculated
 *          OUT: legendre  Array of length P+1 in which to return the
 *                         legendre function for <x>.
 *
 *-------------------------------------------------------------------------*/
void ComputeLegendre(int P, real8 x, real8 *legendre){
        int   l, m;
        real8 factor = -sqrt(1.0-pow(x,2));
        real8 lTemp[((P+2)*(P+1))/2+2];

/*
 *      Don't remove this check.  Seems to be an xlc compiler bug
 *      on the BG/Q system, but having this unnecessary check gets
 *      us past the problem...
 */
        if (ATLM(P,0) < 0) {
            Fatal("Bad value in ComputeLegendre()");
        }

/*
 *      Init legendre
 */
        lTemp[ATLM(0,0)] = 1.0;        /* P_0,0(x) = 1 */
        lTemp[ATLM(1,0)] = x;          /* P_1,0(x) = x */
        lTemp[ATLM(1,1)] = factor;     /* P_1,1(x) = -sqrt(1 - x^2) */

        for (l = 2; l <= P ; ++l ){
            for (m = 0; m < l - 1 ; ++m ){
                /* P_l,m = (2l-1)*x*P_l-1,m - (l+m-1)*x*P_l-2,m / (l-k) */
                lTemp[ATLM(l,m)] = ((real8)(2*l-1) * x * lTemp[ATLM(l-1,m)] -
                                    (real8)(l + m - 1) * lTemp[ATLM(l-2,m)]) /
                                   (real8)(l-m);
            }

            /* P_l,l-1 = (2l-1)*x*P_l-1,l-1 */
            lTemp[ATLM(l,l-1)] = (real8)(2*l-1) * x * lTemp[ATLM(l-1,l-1)];

            /* P_l,l = (2l-1)*factor*P_l-1,l-1 */
            lTemp[ATLM(l,l)] = (real8)(2*l-1) * factor * lTemp[ATLM(l-1,l-1)];
        }

        for (m = 0 ; m <= P ; ++m){
            legendre[m] = lTemp[ATLM(P,m)];
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    CalcGLMandBinomialFactCoeffs()
 *
 *      Description: The Fq table is a sum of the product of a set of
 *                   expansion coefficients <glm> and a factor Q that is
 *                   a product of binomial coefficients.  This function
 *                   computes the glm coefficients using the double integral
 *                   of the angular part of the derivative of the Green's
 *                   function evaluated at grid points and computes Q using
 *                   recurrences.
 *
 *                   Note: The glm expansion coefficients are defined by
 *                   a double integral over two angles, theta and phi.
 *                   These angles are discretized into a grid of
 *                   numTheta X numPhi points.
 *
 *      Parameters:
 *          IN:  qMax        Base factor for determining the number of terms
 *                           in the expansion for spherical harmonics.  The
 *                           total number of terms is 2 * qMax + 1.
 *          IN:  numTheta    Number of theta points in grid
 *          IN:  numPhi      Number of phi points in grid
 *          IN:  phi         Array of <numPhi> angles.
 *          IN:  dGdxgridIn  Array of grid points at which the glm expansion
 *                           coefficients are computed.  Within this function
 *                           the single dimension input array is cast to a
 *                           multi-dimensional array of dimensions:
 *                           [numTheta][numPhi][3][3][3];
 *          IN:  gaussQuadPoints  Array of <numTheta> points for
 *                           gauss quadrature
 *          IN:  gaussQuadWeights Array of <numTheta> weights for
 *                           gauss quadrature.
 *          OUT: glmReIn     Pointer to the real components of the
 *                           expansion coeffients needed for
 *                           calculating the Fq table.  Within this function
 *                           the single dimension input array is cast to a
 *                           multi-dimensional array of dimensions:
 *                           [2*qMax+2][2*qMax+2][3][6]
 *          OUT: glmImIn     Pointer to the real components of the
 *                           expansion coeffients needed for
 *                           calculating the Fq table.  Within this function
 *                           the single dimension input array is cast to a
 *                           multi-dimensional array of dimensions:
 *                           [2*qMax+2][2*qMax+2][3][6]
 *          OUT: binomialFactCoeffsIn Pointer to table of binomial factorial
 *                           coefficients needed for calculating the Fq table.
 *                           Within this function the single dimension input
 *                           array is cast to a multi-dimensional array of
 *                           dimensions: [2*qMax+2][2*qMax+2][2*qMax+2]
 *
 *
 *-------------------------------------------------------------------------*/
static void CalcGLMandBinomialFactCoeffs(int qMax, int numTheta, int numPhi,
                real8 *phi, real8 *dGdxgridIn,
                real8 *gaussQuadPoints, real8 *gaussQuadWeights,
                real8 *glmReIn, real8 *glmImIn,
                real8 *binomialFactCoeffsIn)
{
        int   a, b;
        int   i, j, k, l, m, n;
        int   L, M, intLmM;
        int   lMax = 2*qMax+1;
        real8 factor;
        real8 plm[numTheta][lMax+1];
        real8 plmCol[numTheta];
        real8 q2[lMax+1][lMax+1];
        real8 (*dGdxgrid)[numPhi][3][3][3];                 ptrcpy(dGdxgrid,dGdxgridIn);
        real8 (*binomialFactCoeffs)[2*qMax+2][2*qMax+2];    ptrcpy(binomialFactCoeffs,binomialFactCoeffsIn);

        // "complex" arrays below are represented as separate
        // arrays for the real and imaginary components.

        real8 (*glmRe)[2*qMax+2][3][6];  ptrcpy(glmRe,glmReIn);
        real8 (*glmIm)[2*qMax+2][3][6];  ptrcpy(glmIm,glmImIn);

        real8 glmauxRe[lMax+1][lMax+1][3][3][3];
        real8 glmauxIm[lMax+1][lMax+1][3][3][3];

        real8 GRe[numTheta][numPhi];
        real8 GIm[numTheta][numPhi];

        real8 YlmRe[numTheta][numPhi];
        real8 YlmIm[numTheta][numPhi];

        real8 YlmbarRe[numTheta][numPhi];
        real8 YlmbarIm[numTheta][numPhi];

        real8 expPhiRe[numPhi];
        real8 expPhiIm[numPhi];

        real8 sumYlmbarGRe[numTheta];
        real8 sumYlmbarGIm[numTheta];

        real8  twoPIovernumPhi = 2 * M_PI / numPhi;

        for (i = 0; i < 2*qMax+2; i++) {
            for (j = 0; j < 2*qMax+2; j++) {
                for (k = 0; k < 2*qMax+2; k++) {
                    binomialFactCoeffs[i][j][k] = 0.0;
                }
            }
        }

        binomialFactCoeffs[0][0][0] = 1.0/sqrt(4*M_PI);

        for (L = 1; L <= lMax; L++) {
            factor = sqrt((2.0*L+1.0)*(2.0*L-1.0))/(L*1.0);
            binomialFactCoeffs[0][L][0] = binomialFactCoeffs[0][L-1][0]*factor;

            for (k = 1; k <= L; k++) {
                factor = -(L-2.0*k+2.0)*(L-2.0*k+1.0) /
                          (2.0*k*(2.0*L-2.0*k+1.0));

                binomialFactCoeffs[0][L][k] =
                        binomialFactCoeffs[0][L][k-1] * factor;
            }

            for (M = 1; M <= L; M++) {
                factor = -sqrt((2.0*L-1.0)*(2.0*L+1.0) /
                          (1.0*(L+M)*(L+M-1.0)));

                binomialFactCoeffs[M][L][0] =
                        binomialFactCoeffs[M-1][L-1][0] * factor;

                if ((L-M) % 2 == 0)
                    intLmM = (L-M)/2;
                else
                    intLmM = (L-M-1)/2;

                for (k = 1; k <= intLmM; k++) {
                    factor = sqrt((2.0*L+1.0)/(2.0*L-1.0)) /
                             sqrt(1.0*(L+M)*(L+M-1.0)) *
                             (L-2.0*k-M+2.0)*(L-2.0*k-M+1.0)/(k*2.0);

                    binomialFactCoeffs[M][L][k] =
                            binomialFactCoeffs[M-1][L-1][k-1]*factor;
                }
            }

        }  /* end for (L = 1; L <= lMax; L++) */

/*
 *      Insure all elements of <q2> are zeroed before going any further
 */
        for (i = 0; i < lMax+1; i++) {
            for (j = 0; j < lMax+1; j++) {
               q2[i][j] = 0.0;
            }
        }

        for (i = 0; i <= qMax; i++) {

            l = 2*i+1;

            q2[l][0] = 1.0;

            for (j = 1; j <= l; j++) {

                q2[l][j] = q2[l][j-1] / sqrt((real8)(l-j+1) * (real8)(l+j));
            }
        }


        for (i = 0; i < 3; i++) {
            for (j = 0; j < 3; j++) {
                for (k = 0; k < 3; k++) {

/*
 *                  Copy portion of dGdxgrid into G.  Initially
 *                  only the 'real' portion of G values are set.
 *                  imaginary portions are zeroed for now but will
 *                  be updated later.
 */
                    for (m = 0; m < numTheta; m++) {
                        for (n = 0; n < numPhi; n++) {
                            GRe[m][n] = dGdxgrid[m][n][i][j][k];
                            GIm[m][n] = 0.0;
                        }
                    }

                    for (m = 0; m <= lMax; m++) {
                        real8 nFac0;
/*
 *                      Plm is a numTheta X m+1 matrix (to max of
 *                      numTheta X lMax+1)
 */
                        for (n = 0; n < numTheta; n++) {
                            ComputeLegendre(m, gaussQuadPoints[n], plm[n]);
                        }

                        nFac0 = sqrt((2.0 * (real8)m + 1.0) / FOUR_PI);

                        for (n = 0; n <= m; n++) {
                            real8 nFac;

                            nFac = nFac0 * q2[m][abs(n)];

                            for (a = 0; a < numTheta; a++) {
                                plmCol[a] = plm[a][abs(n)];
                            }

                            for (b = 0; b < numPhi; b++) {
                                expPhiRe[b] = cos(n*phi[b]);
                                expPhiIm[b] = sin(n*phi[b]);
                            }

/*
 *                          Have to init YlmPtr and use that as the arg to
 *                          MatrixMult() because at least one
 *                          version of the IBM compiler on BGQ has an
 *                          internal bug that causes the compiler to fail
 *                          without this.
 */

                            real8 *YlmRePtr = (real8 *) &YlmRe[0][0];
                            real8 *YlmImPtr = (real8 *) &YlmIm[0][0];

                            MatrixMult(plmCol, numTheta, 1, 1,
                                       expPhiRe, numPhi, numPhi,
                                       YlmRePtr, numPhi);

                            MatrixMult(plmCol, numTheta, 1, 1,
                                       expPhiIm, numPhi, numPhi,
                                       YlmImPtr, numPhi);

                            for (a = 0; a < numTheta; a++) {
                                for (b = 0; b < numPhi; b++) {

                                    YlmRe[a][b] *= nFac;
                                    YlmIm[a][b] *= nFac;

                                    YlmbarRe[a][b] = YlmRe[a][b];
                                    YlmbarIm[a][b] = -YlmIm[a][b];
                                }
                            }

                            glmauxRe[m][n][i][j][k] = 0.0;
                            glmauxIm[m][n][i][j][k] = 0.0;

                            for (a = 0; a < numTheta; a++) {

                                sumYlmbarGRe[a] = 0.0;
                                sumYlmbarGIm[a] = 0.0;

                                for (b = 0; b < numPhi; b++) {
                                    sumYlmbarGRe[a] += twoPIovernumPhi *
                                        ((YlmbarRe[a][b] * GRe[a][b]) -
                                         (YlmbarIm[a][b] * GIm[a][b]));

                                    sumYlmbarGIm[a] += twoPIovernumPhi *
                                        ((YlmbarRe[a][b] * GIm[a][b]) +
                                         (YlmbarIm[a][b] * GRe[a][b]));
                                }

                                glmauxRe[m][n][i][j][k] += sumYlmbarGRe[a] *
                                                           gaussQuadWeights[a];
                                glmauxIm[m][n][i][j][k] += sumYlmbarGIm[a] *
                                                           gaussQuadWeights[a];

                            }

                            for (a = 0; a < numTheta; a++) {
                                for (b = 0; b < numPhi; b++) {
                                    GRe[a][b] -=
                                        (glmauxRe[m][n][i][j][k]*YlmRe[a][b]) -
                                        (glmauxIm[m][n][i][j][k]*YlmIm[a][b]);

                                    GIm[a][b] -=
                                        (glmauxRe[m][n][i][j][k]*YlmIm[a][b]) +
                                        (glmauxIm[m][n][i][j][k]*YlmRe[a][b]);
                                }
                            }

                        }  /* end for n ... */
                    }  /* end for m ... */
                }  /* end for k ... */
            }  /* end for j ... */
        }  /* end for i ... */

        for (i = 0; i < 3; i++) {
            for (j = 0; j <= lMax; j++) {
                for (k = 0; k <= j; k++) {
                    glmRe[j][k][i][0] = glmauxRe[j][k][i][0][0];
                    glmRe[j][k][i][1] = glmauxRe[j][k][i][1][1];
                    glmRe[j][k][i][2] = glmauxRe[j][k][i][2][2];
                    glmRe[j][k][i][3] = glmauxRe[j][k][i][1][2] +
                                        glmauxRe[j][k][i][2][1];
                    glmRe[j][k][i][4] = glmauxRe[j][k][i][2][0] +
                                        glmauxRe[j][k][i][0][2];
                    glmRe[j][k][i][5] = glmauxRe[j][k][i][0][1] +
                                        glmauxRe[j][k][i][1][0];

                    glmIm[j][k][i][0] = glmauxIm[j][k][i][0][0];
                    glmIm[j][k][i][1] = glmauxIm[j][k][i][1][1];
                    glmIm[j][k][i][2] = glmauxIm[j][k][i][2][2];
                    glmIm[j][k][i][3] = glmauxIm[j][k][i][1][2] +
                                        glmauxIm[j][k][i][2][1];
                    glmIm[j][k][i][4] = glmauxIm[j][k][i][2][0] +
                                        glmauxIm[j][k][i][0][2];
                    glmIm[j][k][i][5] = glmauxIm[j][k][i][0][1] +
                                        glmauxIm[j][k][i][1][0];
                }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    FindNonZeroFqComponents()
 *
 *      Description: This function locates all components of the specified
 *                   <FqTable> that are non-zero and returns a 1D array
 *                   containing all such components.
 *
 *                   Within the force routine there are many zero entries in
 *                   S(i,j,ii) (See the "Use of spherical harmonics for
 *                   disloction dynamics in anisotropic elastic media" paper
 *                   in the docs directory).  Rather than do many zero
 *                   multiplies and have many zeros stored in the matrix the
 *                   array is stored in a column vector with the zeros
 *                   removed.  Additional arrays are needed to specify which
 *                   entries remain and to make effective use of the
 *                   analytical integrals used to compute forces.
 *
 *                   This function also calculates that set of arrays and
 *                   returns them to the caller.
 *
 *                   For example, if indmArray[] is the three-element array
 *                   indmAray = {1, 0, 3}, it indicates there is 1 (j,ii)
 *                   element with non-zero components for i=0, none for
 *                   i=1 and 3 for i=2.  The indlArray would then consist
 *                   of 4 elements.  The first being the j index of the
 *                   single (j,ii) element with a non-zero component for
 *                   i=1, the remaining three being the j indices of the
 *                   (j,ii) elements with non-zero components for i=2.
 *                   The ind18Array then contains 4 elements (1 for each
 *                   index in indlArray) each giving the count of non-zero
 *                   elements within the ii block (each ii block represents
 *                   18 components) for the corresponding j index.  Lastly,
 *                   indiiArray consists of the indices (range 0 - 17) of
 *                   the non-zero components of each ii block specified
 *                   in ind18Array (similar to the relationship between
 *                   indlArray and indmArray).
 *
 *      Parameters:
 *          IN:   qMax       Base factor for determining the number of terms
 *                           in the expansion for spherical harmonics.  The
 *                           total number of terms is 2 * qMax + 1.
 *          IN:  maxVal      Maximum absolute value of any component from the
 *                           array <FqTable>.
 *          IN:  FqTable     Pointer to the full Fq table.  Internally, the
 *                           function casts the 1D array to a 3D table of
 *                           dimensions [6][3][(qMax+2)*(qMax+1)]
 *          OUT: indmArray   Based on S(i,j,ii) (see function description
 *                           above), array in which to return the number of
 *                           j elements for each i that have at least one
 *                           (j,ii) element that is non-zero.
 *          OUT: indlArray   Based on S(i,j,ii) (see function description
 *                           above), array in which to return the j indices
 *                           associated with the non-zero entries in indmArray.
 *                           The number of elements in indlArray equals
 *                           the sum of all elements in indmArray.
 *          OUT: ind18Array  Based on S(i,j,ii) (see function description
 *                           above), array in which to return the number
 *                           of ii elements associated with an (i,j) that
 *                           are non-zero.
 *          OUT: indiiArray  Based on S(i,j,ii) (see function description
 *                           above), array in which to return the indices
 *                           (between 0 and 17) associated with the non-zero
 *                           entries in ind18Array.  The number of elements
 *                           in indiiArray equals the sum of all elements in
 *                           ind18Array.
 *          OUT: numElem     Number of elements in  <FqOut> and <indiiArray>
 *                           arrays.  Maximum size of these arrays is
 *                           (qMax+1)*(qMax+2)*18
 *          OUT: numElem2    Number of elements in each of the <indlArray>
 *                           and <ind18Array> arrays.  Maximum size of
 *                           these arrays is (qMax+1)*(qMax+2)
 *
 *
 *-------------------------------------------------------------------------*/
static void FindNonZeroFqComponents(int qMax, real8 maxF, real8 *FqIn,
                                    short *indmArray, short *indlArray,
                                    short *ind18Array, short *indiiArray,
                                    real8 *FqOut,
                                    int *numElem,int *numElem2)
{
        int   i,j,ii;
        int   indl,indi;
        int   local, localRe;
        int   local3;
        int   qMaxProd = (qMax+2)*(qMax+1);
        real8 *FqArray[18];

        real8 (*Fq)[3][qMaxProd]; ptrcpy(Fq,FqIn);

        for (i = 0; i < 18; i++) {
            FqArray[i] = (real8 *)calloc(qMaxProd, sizeof(real8));
        }

        for (ii = 0; ii < 3; ii++) {
            for (j = 0; j < 6; j++) {
                for (local = 0; local < (qMax+2)*(qMax+1); local++) {
                    FqArray[6*ii+j][local] = Fq[j][ii][local];
                }
            }
        }

/*
 *      local is the general counter for the Fq table. It counts the zeros
 *      and the nonzeros elements among the (qMax+1)*(qMax+2) elements of
 *      the Fqtable.
 *
 *      local3 counts the number of elements among the (qMax+1)*(qMax+2)*18
 *      elements in the Fq table that are non zeros.
 *
 *      localRe counts the number of elements among the (qMax+1)*(qMax+2)
 *      elements in the Fq table that are non zeros.
 */
        local3 = 0;
        local  = 0;

        for (ii = 0; ii < 18; ii++) {
            FqOut[local3] = FqArray[ii][local];
            indiiArray[local3] = 0.0;
            local3++;
        }

        local = 1;

        for (ii = 0; ii < 18; ii++) {
            FqOut[local3] = FqArray[ii][local];
            indiiArray[local3] = 0.0;
            local3++;
        }

        local = 2;
        localRe = 2;

        for (i = 2; i < qMax+2; i++) {

            indl=1;

            for (j = 1; j < 2*i+1; j++) {
                int isNonZero = 0;

                for (ii = 0; ii < 18; ii++) {
                    if (fabs(FqArray[ii][local]) >  maxF*1e-15) {
                        isNonZero = 1;
                    }
                }

                if (isNonZero == 1) {
                    indi = 0;

                    for (ii = 0; ii < 18; ii++) {
                        if (fabs(FqArray[ii][local]) >  maxF*1e-15) {
                            FqOut[local3] = FqArray[ii][local];
                            indiiArray[local3] = ii;
                            local3++;
                            indi++;
                        }
                    }

                    ind18Array[localRe] = indi;
                    indlArray[localRe] = j;

                    localRe++;
                    indl++;
                }

                local++;
            }

            indmArray[i] = indl;
        }

        *numElem  = local3;
        *numElem2 = localRe;

        for (i = 0; i < 18; i++) {
            free(FqArray[i]);
        }

        return;
}


static void FindNonZeroFqComponentsIm(int qMax, real8 maxF, real8 *FqIn,
                                      short *indmArray, short *indlArray,
                                      short *ind18Array, short *indiiArray,
                                      real8 *FqOut,
                                      int *numElem,int *numElem2)
{
        int   i,j,ii;
        int   indl,indi;
        int   local, localRe;
        int   local3;
        int   qMaxProd = (qMax+2)*(qMax+1);
        real8 *FqArray[18];

        real8 (*Fq)[3][qMaxProd]; ptrcpy(Fq,FqIn);

        for (i = 0; i < 18; i++) {
            FqArray[i] = (real8 *)calloc(qMaxProd, sizeof(real8));
        }

        for (ii = 0; ii < 3; ii++) {
            for (j = 0; j < 6; j++) {
                for (local = 0; local < (qMax+2)*(qMax+1); local++) {
                    FqArray[6*ii+j][local] = Fq[j][ii][local];
                }
            }
        }

/*
 *      local is the general counter for the Fq table. It counts the zeros
 *      and the nonzeros elements among the (qMax+1)*(qMax+2) elements of
 *      the Fqtable.
 *
 *      local3 counts the number of elements among the (qMax+1)*(qMax+2)*18
 *      elements in the Fq table that are non zeros.
 *
 *      localRe counts the number of elements among the (qMax+1)*(qMax+2)
 *      elements in the Fq table that are non zeros.
 */
        local  = 1;
        local3 = 0;

        for (ii = 0; ii < 18; ii++) {
            FqOut[local3] = FqArray[ii][local];
            indiiArray[local3] = 0.0;
            local3++;
        }

        local = 2;
        localRe = 2;

        for (i = 2; i < qMax+2; i++) {

            indl=1;

            for (j = 1; j < 2*i+1; j++) {
                int isNonZero = 0;

                for (ii = 0; ii < 18; ii++) {
                    if (fabs(FqArray[ii][local]) >  maxF*1e-15) {
                        isNonZero = 1;
                    }
                }

                if (isNonZero == 1) {
                    indi = 0;

                    for (ii = 0; ii < 18; ii++) {
                        if (fabs(FqArray[ii][local]) >  maxF*1e-15) {
                            FqOut[local3] = FqArray[ii][local];
                            indiiArray[local3] = ii;
                            local3++;
                            indi++;
                        }
                    }

                    ind18Array[localRe] = indi;
                    indlArray[localRe] = j;

                    localRe++;
                    indl++;
                }

                local++;
            }

            indmArray[i] = indl;
        }

        *numElem  = local3;
        *numElem2 = localRe;

        for (i = 0; i < 18; i++) {
            free(FqArray[i]);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    SetElasticConstantMatrix()
 *
 *      Description: Using the individual elastic constants provided by
 *                   the user, initialize the full 6X6 elastic constant
 *                   matrix appropriate to the material being simulated.
 *
 *                   NOTE: This is currently on valid for BCC, FCC and
 *                         HCP materials.  For materials with a rhombohedral
 *                         structure, the user is required to provide the
 *                         full elastic constant matrix explicitly.
 *
 *-------------------------------------------------------------------------*/
void SetElasticConstantMatrix(Home_t *home)
{
        int     i, j;
        real8   C66;
        Param_t *param;

        param = home->param;

/*
 *      First zero out the whole array
 */
        for (i = 0; i < 6; i++) {
            for (j = 0; j < 6; j++) {
                param->elasticConstantMatrix[i][j] = 0.0;
            }
        }

/*
 *      Fill in the proper constants based on the type of material
 *      in use.
 *
 *      NOTE: Rhombohedral is not handled here because the user
 *            is currently required to provide the entire elastic
 *            constamt matrix in the control file for those materials
 */
        switch (home->param->materialType) {
            case MAT_TYPE_BCC:
            case MAT_TYPE_FCC:

                param->elasticConstantMatrix[0][0] = param->C11;
                param->elasticConstantMatrix[0][1] = param->C12;
                param->elasticConstantMatrix[0][2] = param->C12;

                param->elasticConstantMatrix[1][0] = param->C12;
                param->elasticConstantMatrix[1][1] = param->C11;
                param->elasticConstantMatrix[1][2] = param->C12;

                param->elasticConstantMatrix[2][0] = param->C12;
                param->elasticConstantMatrix[2][1] = param->C12;
                param->elasticConstantMatrix[2][2] = param->C11;

                param->elasticConstantMatrix[3][3] = param->C44;
                param->elasticConstantMatrix[4][4] = param->C44;
                param->elasticConstantMatrix[5][5] = param->C44;

                break;

            case MAT_TYPE_HCP:

                C66 = 0.5 * (param->C11 - param->C12);

                param->elasticConstantMatrix[0][0] = param->C11;
                param->elasticConstantMatrix[0][1] = param->C12;
                param->elasticConstantMatrix[0][2] = param->C13;

                param->elasticConstantMatrix[1][0] = param->C12;
                param->elasticConstantMatrix[1][1] = param->C11;
                param->elasticConstantMatrix[1][2] = param->C13;

                param->elasticConstantMatrix[2][0] = param->C13;
                param->elasticConstantMatrix[2][1] = param->C13;
                param->elasticConstantMatrix[2][2] = param->C33;

                param->elasticConstantMatrix[3][3] = param->C44;
                param->elasticConstantMatrix[4][4] = param->C44;
                param->elasticConstantMatrix[5][5] = C66;

                break;
        }

        return;
}



/*---------------------------------------------------------------------------
 *
 *      Function:    SetElasticConstantMatrix4D()
 *
 *      Description: Uses the full 6x6 elastic constant matrix for
 *                   the material to initialize the 3X3X3X3 matrix
 *                   of elastic constants required by some portions
 *                   of the code.
 *
 *      Parameters:
 *          IN:     ecMatrix    Elastic constants as a 6x6 matrix
 *          OUT:    ecMatrix4D  Elastic constants from <ecMatrix> transformed
 *                              into a 3x3x3x3 matrix.
 *
 *-------------------------------------------------------------------------*/
void SetElasticConstantMatrix4D(real8 ecMatrix[6][6],
                                real8 ecMatrix4D[3][3][3][3])
{
        int i3to6[3][3] = { {0, 5, 4}, {5, 1, 3}, {4, 3, 2} };

        for (int i = 0; (i<3); i++)
        for (int j = 0; (j<3); j++)
        for (int k = 0; (k<3); k++)
        for (int l = 0; (l<3); l++)
           ecMatrix4D[i][j][k][l] = ecMatrix[i3to6[i][j]][i3to6[k][l]];

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    SetCubicElasticModuli()
 *
 *      Description: Uses provided elastic constants to calculates
 *                   elastic moduli (MU, NU, etc) for cubic crystals
 *                   using Voight averaging.
 *
 *-------------------------------------------------------------------------*/
static void SetCubicElasticModuli(Home_t *home)
{
        Param_t *param;
        real8 C12, C44, Cpr;

        param = home->param;

        C12 = param->C12;
        C44 = param->C44;
        Cpr = param->Cpr;

/*
 *      Use VOIGT averaging
 */
        param->shearModulus  = (2.00*Cpr+3.00*C44)/5.00;

        param->pois          = 0.50*(2.00*Cpr+5.00*C12-2.00*C44)/
                               (4.00*Cpr+5.00*C12+C44);

        param->YoungsModulus = (4.00*Cpr*Cpr + 6.00*Cpr*C12 +
                                6.00*Cpr*C44 + 9.00*C12*C44) /
                               (4.00*Cpr+5.00*C12+C44);

        if (home->myDomain == 0) {
            printf("Recalculating constants for cubic crystal anisotropy:\n");
            printf("    Shear Modulus = %.8e\n", param->shearModulus);
            printf("    Poisson Ratio = %.8e\n", param->pois);
            printf("    Young Modulus = %.8e\n", param->YoungsModulus);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    SetHexagonalElasticModuli()
 *
 *      Description: Uses provided elastic constants to calculates
 *                   elastic moduli (MU, NU, etc) for hexagonal crystals
 *                   using Voight averaging.
 *
 *-------------------------------------------------------------------------*/
static void SetHexagonalElasticModuli(Home_t *home)
{
        Param_t *param;
        real8 C11, C12, C13, C33, C44;

        param = home->param;

        C11 = param->C11;
        C12 = param->C12;
        C13 = param->C13;
        C33 = param->C33;
        C44 = param->C44;

/*
 *      Use VOIGT averaging: Hirth and Lothe p. 436
 */
        param->shearModulus = (1.00/30.00) *
                              ( 7.00*C11 - 5.00*C12 + 2.00*C33 +
                               12.00*C44 - 4.00*C13);

        param->pois = (C11 +      C33 + 5.00*C12 - 4.00*C44 +  8.00*C13 )/
                      (9.00*C11 + 4.00*C33 + 5.00*C12 + 4.00*C44 + 12.00*C13);

        param->YoungsModulus = 2 * param->shearModulus * (1+ param->pois);

        if (home->myDomain == 0) {
            printf("Recalculating constants for HCP crystal anisotropy\n");
            printf("    Shear Modulus = %.8e\n", param->shearModulus);
            printf("    Poisson Ratio = %.8e\n", param->pois);
            printf("    Young Modulus = %.8e\n", param->YoungsModulus);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    AnisotropicInitElasticConstants()
 *
 *      Description: This function determines whether the user has provided
 *                   the individual elastic constants and/or the elastic
 *                   constant matrix.  Sanity checks are done to verify
 *                   that the minimum necessary information has been
 *                   provided, then the function uses the available values
 *                   to insure all neede elastic constants, the elastic
 *                   constant matrix (in both its 6x6 and 3x3x3x3 form),
 *                   and the elastic moduli such as MU, Nu and Youngs
 *                   modulus are initialize appropriately.
 *
 *-------------------------------------------------------------------------*/
static void AnisotropicInitElasticConstants(Home_t *home)
{
        int     haveElasticConstantMatrix = 0;;
        int     haveElasticConstants = 0;;
        int     haveElasticModulii = 0;;
        int     haveCubicConstants = 0;;
        int     haveHexagonalConstants = 0;;
        real8   old_NU=0.0, old_MU=0.0;
        real8   epsC = 1e-5;
        Param_t *param;
        AnisotropicVars_t *anisoVars;



        param = home->param;
        anisoVars = &home->anisoVars;

/*
 *      Check if the user provided the full elastic constant matrix
 *      for the material
 */
        if (fabs(param->elasticConstantMatrix[1][1]) > epsC ||
            fabs(param->elasticConstantMatrix[1][2]) > epsC ||
            fabs(param->elasticConstantMatrix[1][3]) > epsC ||
            fabs(param->elasticConstantMatrix[3][3]) > epsC) {

            haveElasticConstantMatrix = 1;
        }

/*
 *      Check if the user provided crystal specific elastic constants
 */
        if (((param->materialType == MAT_TYPE_BCC) ||
             (param->materialType == MAT_TYPE_FCC)) &&
            (fabs(param->C11) * fabs(param->C12) * fabs(param->C44) > epsC)) {
            haveCubicConstants = 1;
            haveElasticConstants = 1;
        }

        if ((param->materialType == MAT_TYPE_HCP) &&
            (fabs(param->C11) * fabs(param->C12) * fabs(param->C13) *
             fabs(param->C33) * fabs(param->C44) > epsC)) {
            haveHexagonalConstants = 1;
            haveElasticConstants = 1;
        }

/*
 *      Check if the elastic modulii were provided
 */
        if ((fabs(param->shearModulus) > epsC) && (fabs(param->pois) > epsC)) {
            haveElasticModulii = 1;
            old_NU = param->shearModulus;
            old_MU = param->pois;
        }

/*
 *      If neither the elastic constants matrix nor the crystal specific
 *      elastic constants have been provided, we have a problem
 */
        if (!haveElasticConstantMatrix && !haveElasticConstants) {
            if (home->myDomain == 0) {
                Fatal("Either the <elasticConstants> matrix or the "
                      "individual\nelastic constants (C11, C12, etc) must "
                      "be provided when\nusing anisotropic elasticity.");
            }
        }

/*
 *      For rhombohedral crystals, we currently require the user to
 *      provide the full elastic constants matrix plus the elastic
 *      modulii
 */
        if (param->materialType == MAT_TYPE_RHOMBOHEDRAL_VA) {

            if ((home->myDomain == 0) && (haveElasticConstantMatrix == 0)) {
                Fatal("The <elasticConstants> must be defined when "
                      "using anisotropic\nelasticity with rhombohedral "
                      "material types");
            }

            if ((home->myDomain == 0) && (haveElasticModulii == 0)) {
                Fatal("<shearModulus> and <pois> must be defined when using "
                      "anisotropic\nelasticity with rhombohedral material "
                      "types");
            }

            anisoVars->Set_EC_Matrices(param->elasticConstantMatrix);

        } else {
/*
 *          For material structures other than rhombohedral, if the user
 *          provided the full elastic constant matrix it will take
 *          precedence over any crystal specific elastic constants
 *          provided.
 */
            if (haveElasticConstantMatrix) {

                anisoVars->Set_EC_Matrices(param->elasticConstantMatrix);

/*
 *              Explicity define the individual elastic constants from
 *              the user-specified matrix.  Needed so the elastic modulii
 *              can be recalculated.
 */
                if ((param->materialType == MAT_TYPE_BCC) ||
                    (param->materialType == MAT_TYPE_FCC)) {

                    param->C11 = anisoVars->elasticConstantMatrix4D[0][0][0][0];
                    param->C12 = anisoVars->elasticConstantMatrix4D[0][0][1][1];
                    param->Cpr = 0.5 * (anisoVars->elasticConstantMatrix4D[0][0][0][0] -
                                        anisoVars->elasticConstantMatrix4D[0][0][1][1]);

                    param->C44 = anisoVars->elasticConstantMatrix4D[0][1][0][1];

                } else if (param->materialType == MAT_TYPE_HCP) {

                    param->C11 = anisoVars->elasticConstantMatrix4D[0][0][0][0];
                    param->C12 = anisoVars->elasticConstantMatrix4D[0][0][1][1];
                    param->C13 = anisoVars->elasticConstantMatrix4D[2][2][0][0];
                    param->C33 = anisoVars->elasticConstantMatrix4D[2][2][2][2];
                    param->C44 = anisoVars->elasticConstantMatrix4D[0][1][0][1];
                }

            } else {
/*
 *              We don't have the full elastic constant matrix, so the
 *              2D and 4D versions of the elastic constant matrix must
 *              be built from the individual crystal specific elastic
 *              constants.
 */
                switch (param->materialType) {

                    case MAT_TYPE_BCC:
                    case MAT_TYPE_FCC:

                        if ((home->myDomain == 0) && (!haveCubicConstants)) {
                            Fatal("One of the elastic constants C12, C44, "
                                  "or Cpr has not been provided");
                        }

                        SetElasticConstantMatrix(home);

                        break;

                    case MAT_TYPE_HCP:

                        if ((home->myDomain == 0)&&(!haveHexagonalConstants)) {
                            Fatal("One of the elastic constants C11, C44, "
                                  "C13, C33, or C44 has not been provided");
                        }

                        SetElasticConstantMatrix(home);

                        break;

                }  /* end of switch(materialType) */

                anisoVars->Set_EC_Matrices(param->elasticConstantMatrix);
            }

/*
 *          Use the elastic constants to explicitly calculate
 *          the proper elastic modulii.
 */
            switch (param->materialType) {

                case MAT_TYPE_BCC:
                case MAT_TYPE_FCC:
/*
 *                  We need the Cpr valuye when calculating the elastic
 *                  moduli, but it *may* not have been set above...
 */
                    param->Cpr = 0.5 * (param->C11 - param->C12);

                    SetCubicElasticModuli(home);
                    break;

                case MAT_TYPE_HCP:
                    SetHexagonalElasticModuli(home);
                    break;
            }

/*
 *          If the user provided values for the elastic modulii that are
 *          not consistant with the elastic constants, issue a warning.
 */
            if ((home->myDomain == 0) && haveElasticModulii) {

                if (fabs((old_NU-param->shearModulus)/old_NU) > epsC) {
                    printf("WARNING: <shearModulus> reset to %e from %e\n",
                           param->shearModulus, old_NU);
                }

                if (fabs((old_MU-param->pois)/old_MU) > epsC) {
                    printf("WARNING: <pois> reset to %e from %e\n",
                           param->pois, old_MU);
                }
            }

        }  /* end code for non-rhombohedral materials */

/*
 *      If the user has not explicitly provided an Ecore value, we
 *      need to set it now based on the newly calculated shear modulus.
 */
        if (param->Ecore < 0.0) {
            param->Ecore = log(param->rc/0.1);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    CalcCoreFqTable()
 *
 *      Description: The spherical harmonics decomposition of the line
 *                   tension core force can be
 *                   written as the sum of the product of two
 *                   functions: FqCore and B, where Fq can be precomputed
 *                   independently of the dislocation segment. This
 *                   function computes the pre-factor FqCore.
 *
 *
 *      Parameters:
 *          IN:  qMax        Base factor for determining the number of terms
 *                           in the expansion for spherical harmonics.  The
 *                           total number of terms is 2 * qMax + 1.
 *          IN:  numTheta    Number of theta points in grid
 *          IN:  numPhi      Number of phi points in grid
 *          IN:  phi         Array of <numPhi> angles.
 *          IN:  BgridIn     Array of grid points at which the Blm expansion
 *                           coefficients are computed.  Within this function
 *                           the single dimension input array is cast to a
 *                           multi-dimensional array of dimensions:
 *                           [numTheta][numPhi][6];
 *          IN:  gaussQuadPoints  Array of <numTheta> points for
 *                           gauss quadrature
 *          IN:  gaussQuadWeights Array of <numTheta> weights for
 *                           gauss quadrature.
 *          OUT:  FqRealPtrIn Pointer the first element of the array
 *                           of the real components of the Fqtable.  Dimensions
 *                           of the table are [6][(qMax+2)*(qMax+1)-1]
 *          OUT:  FqImagPtrIn Pointer the first element of the array
 *                           of the real components of the Fqtable.  Dimensions
 *                           of the table are [6][(qMax+2)*(qMax+1)-1]
 *
 *-------------------------------------------------------------------------*/
void CalcCoreFqTable(int qMax, int numTheta, int numPhi,
                     real8 *phi, real8 *BgridIn,
                     real8 *gaussQuadPoints, real8 *gaussQuadWeights,
                     real8 *FqRealPtrIn, real8 *FqImagPtrIn)
{
	
        int   a, b;
        int   i, j, k, l, m, n;
        int   L, M, intLmM;
        int   lMax = 2*qMax;  // for line tension force, it is 2*qMax.
        int   qMaxProd = (qMax+2)*(qMax+1)-1;
        real8 factor;
        real8 twoPIovernumPhi;
        real8 plm[numTheta][lMax+1];
        real8 plmCol[numTheta];
        real8 q2[lMax+1][lMax+1];
        real8 binomialFactCoeffs[2*qMax+2][2*qMax+2][2*qMax+2];

        //
        // "complex" arrays below are represented as pairs of
        // arrays for the real and imaginary components.
        //
        real8   BlmReal[lMax+1][lMax+1][6];
        real8   BlmImag[lMax+1][lMax+1][6];

        real8   GReal[numTheta][numPhi];
        real8   GImag[numTheta][numPhi];

        real8   YlmReal[numTheta][numPhi];
        real8   YlmImag[numTheta][numPhi];

        real8   YlmbarReal[numTheta][numPhi];
        real8   YlmbarImag[numTheta][numPhi];

        real8   expPhiReal[numPhi];
        real8   expPhiImag[numPhi];

        real8   sumYlmbarGReal[numTheta];
        real8   sumYlmbarGImag[numTheta];

        real8  *YlmRealPtr;
        real8  *YlmImagPtr;

        real8 (*FqReal)[qMaxProd];  ptrcpy(FqReal,FqRealPtrIn);
        real8 (*FqImag)[qMaxProd];  ptrcpy(FqImag,FqImagPtrIn);

        real8 (*Bgrid)[numPhi][6];  ptrcpy(Bgrid,BgridIn);

/*
 *      Compute binomial terms
 */
        twoPIovernumPhi = 2 * M_PI / numPhi;

        for (i = 0; i < 2*qMax+2; i++) {
            for (j = 0; j < 2*qMax+2; j++) {
                for (k = 0; k < 2*qMax+2; k++) {
                    binomialFactCoeffs[i][j][k] = 0.0;
                }
            }
        }

        binomialFactCoeffs[0][0][0] = sqrt(INV_4_OVER_PI);

        for (L = 1; L <= lMax; L++) {
            factor = sqrt((2.0*L+1.0)*(2.0*L-1.0))/(L*1.0);
            binomialFactCoeffs[0][L][0] = binomialFactCoeffs[0][L-1][0]*factor;

            for (k = 1; k <= L; k++) {
                factor = -(L-2.0*k+2.0)*(L-2.0*k+1.0) /
                          (2.0*k*(2.0*L-2.0*k+1.0));

                binomialFactCoeffs[0][L][k] =
                        binomialFactCoeffs[0][L][k-1] * factor;
            }

            for (M = 1; M <= L; M++) {
                factor = -sqrt((2.0*L-1.0)*(2.0*L+1.0) /
                          (1.0*(L+M)*(L+M-1.0)));

                binomialFactCoeffs[M][L][0] =
                        binomialFactCoeffs[M-1][L-1][0] * factor;

                if ((L-M) % 2 == 0)
                    intLmM = (L-M)/2;
                else
                    intLmM = (L-M-1)/2;

                for (k = 1; k <= intLmM; k++) {
                    factor = sqrt((2.0*L+1.0)/(2.0*L-1.0)) /
                             sqrt(1.0*(L+M)*(L+M-1.0)) *
                             (L-2.0*k-M+2.0)*(L-2.0*k-M+1.0)/(k*2.0);

                    binomialFactCoeffs[M][L][k] =
                            binomialFactCoeffs[M-1][L-1][k-1]*factor;
                }
            }

        }  /* end for (L = 1; L <= lMax; L++) */

/*
 *      Insure all elements of <q2> are zeroed before going any further
 */
        for (i = 0; i < lMax+1; i++) {
            for (j = 0; j < lMax+1; j++) {
               q2[i][j] = 0.0;
            }
        }

        for (i = 0; i <= qMax; i++) {

            l = 2*i;  // 2* l here

            q2[l][0] = 1.0;

            for (j = 1; j <= l; j++) {

                q2[l][j] = q2[l][j-1] / sqrt((real8)(l-j+1) * (real8)(l+j));
            }
        }

/*
 *      Compute Blm coefficients
 */
        for (i = 0; i < 6; i++) {

/*
 *          Copy portion of Bgrid into G.  Initially
 *          only the 'real' portion of G values are set.
 *          further on the imaginary portions may be updated.
 */
            for (m = 0; m < numTheta; m++) {
                for (n = 0; n < numPhi; n++) {
                    GReal[m][n] = Bgrid[m][n][i];
                    GImag[m][n] = 0.0;
                }
            }

            for (m = 0; m <= lMax; m++) {
                real8 nFac0;
/*
 *              Plm is a numTheta X m+1 matrix (to max of
 *              numTheta X lMax+1)
 */
                for (n = 0; n < numTheta; n++) {
                    ComputeLegendre(m, gaussQuadPoints[n], plm[n]);
                    //printf("plm[%d][%d]=%f\n",n,m,plm[n][m]);
                }

                nFac0 = sqrt((2.0 * (real8)m + 1.0) * INV_4_OVER_PI);

                for (n = 0; n <= m; n++) {
                    real8 nFac;

                    nFac = nFac0 * q2[m][abs(n)];

                    for (a = 0; a < numTheta; a++) {
                        plmCol[a] = plm[a][abs(n)];
                    }

                    for (b = 0; b < numPhi; b++) {
                        expPhiReal[b] = cos(n*phi[b]);
                        expPhiImag[b] = sin(n*phi[b]);
                    }

/*
 *                  Have to init YlmPtr and use that as the arg to
 *                  MatrixDblMultComplex() because at least one version of
 *                  the IBM compiler on BGQ has an internal bug that causes
 *                  the compiler to fail without this.
 */
                    YlmRealPtr = &YlmReal[0][0];
                    YlmImagPtr = &YlmImag[0][0];

                    MatrixMult(plmCol, numTheta, 1, 1,
                               expPhiReal, numPhi, numPhi,
                               YlmRealPtr, numPhi);

                    MatrixMult(plmCol, numTheta, 1, 1,
                               expPhiImag, numPhi, numPhi,
                               YlmImagPtr, numPhi);

                    for (a = 0; a < numTheta; a++) {
                        for (b = 0; b < numPhi; b++) {

                            YlmReal[a][b] *= nFac;
                            YlmImag[a][b] *= nFac;

                            YlmbarReal[a][b] =  YlmReal[a][b];
                            YlmbarImag[a][b] = -YlmImag[a][b];
                        }
                    }

                    BlmReal[m][n][i] = 0.0;
                    BlmImag[m][n][i] = 0.0;

                    for (a = 0; a < numTheta; a++) {

                        sumYlmbarGReal[a] = 0.0;
                        sumYlmbarGImag[a] = 0.0;

                        for (b = 0; b < numPhi; b++) {

                            sumYlmbarGReal[a] += twoPIovernumPhi *
                                        ((YlmbarReal[a][b] * GReal[a][b]) -
                                         (YlmbarImag[a][b] * GImag[a][b]));

                            sumYlmbarGImag[a] += twoPIovernumPhi *
                                        ((YlmbarReal[a][b] * GImag[a][b]) +
                                         (YlmbarImag[a][b] * GReal[a][b]));
                        }

                        BlmReal[m][n][i] += sumYlmbarGReal[a] *
                                            gaussQuadWeights[a];

                        BlmImag[m][n][i] += sumYlmbarGImag[a] *
                                            gaussQuadWeights[a];

                    }

                    for (a = 0; a < numTheta; a++) {
                        for (b = 0; b < numPhi; b++) {

                            GReal[a][b] -= (BlmReal[m][n][i] * YlmReal[a][b]) -
                                           (BlmImag[m][n][i] * YlmImag[a][b]);

                            GImag[a][b] -= (BlmReal[m][n][i] * YlmImag[a][b]) +
                                           (BlmImag[m][n][i] * YlmReal[a][b]);
                    }
                }

            }  /* end for n ... */
        }  /* end for m ... */
    }  /* end for i ... */

/*
 *  Compute the product of binomial and Blm function in prefactor
 *  arrays.
 *
 *  Initially zero the arrays for the real and imaginary portions of
 *  the Fq table.
 */
    for (i = 0; i < 6; i++) {
        for (m = 0; m < qMaxProd; m++) {
            FqReal[i][m] = 0.0;
            FqImag[i][m] = 0.0;
        }
    }

/*
 *  Product of binomial and Blm coefficients in prefactor
 */
    int ind1;

    for (i = 0; i <= qMax; i++) {
        real8 fac;

        l = 2 * i; // different here than segseg force prefactor

        for (j = 0; j <= l; j++) {
            int intlmj;

            if (j == 0) {
                fac = 1;
            } else {
                fac = 2;
            }

            if ((l-j) % 2 == 0) {
                intlmj = (l-j) / 2;
            } else {
                intlmj = (l-j-1) / 2;
            }

            ind1 = j;
            for (k = 0; k <= intlmj; k++) {
                int   ind, local;
                real8 val;

                ind = i-k;
                local = ind*ind + ind1;

                val = fac * binomialFactCoeffs[j][l][k];

                for (m = 0; m < 6; m++) {

                    FqReal[m][local] += val * BlmReal[l][j][m];
                    FqImag[m][local] += val * BlmImag[l][j][m];

                }

            }  // end for (k = 0; k <= intlmj; k++)

        }  // end for (j = 0; j <= l; j++)

    }  // end for (i = 0; i <= qMax; i++)

    return;
}

//---------------------------------------------------------------------------
// CalcFqTableV3()
//
// Converts an original, uncompressed <Fq> table to a version using by the
// GPU enabled aniso forces.
//---------------------------------------------------------------------------

static void CalcFqTableV3(real8 *fq_dst, const real8 *fq_src, const int qmax)
{
    if (fq_dst && fq_src)
    {
        real8 (*fq  )[3][(qmax+2)*(qmax+1)]=0;   ptrcpy(fq  ,fq_src);
        real8 (*fqv3)   [(qmax+2)*(qmax+1)]=0;   ptrcpy(fqv3,fq_dst);

        for (int i=0; i<3; i++)
        for (int j=0; j<6; j++)
        for (int k=0; k < (qmax+2)*(qmax+1); k++)
           fqv3[6*i+j][k]=fq[j][i][k];
    }
}

/*---------------------------------------------------------------------------
 *
 *      Function:    AnisotropicInit()
 *
 *      Description: This function is called only one during initialization
 *                   to insure proper setup of elastic constants, the
 *                   elastic constant matrix, elastic moduli, and to
 *                   create certain static tables required during the
 *                   anisotropic force calculations.
 *
 *-------------------------------------------------------------------------*/
void AnisotropicInit(Home_t *home)
{
        Param_t *param          = home->param;
        int      numPhiPoints   = param->anisoNumPhiPoints;
        int      numThetaPoints = param->anisoNumThetaPoints;
        int      qMax           = param->anisoHarmonicsNumTermsBase;
        int      i, j;
        real8    FqRealMax, FqImagMax;
        real8   *binomialFactCoeffsPtr;
        real8   dGdx[3][3][3];
        /* full dimensions are dGdxgrid[numTheta][numPhi][3][3][3] */
        /* full dimensions are Bgrid[numTheta][numPhi][6] */
        real8   binomialFactCoeffs[2*qMax+2][2*qMax+2][2*qMax+2];
        real8   glmRe[2*qMax+2][2*qMax+2][3][6];
        real8   glmIm[2*qMax+2][2*qMax+2][3][6];

        if (home->myDomain == 0) { printf("Initializing anisotropy tables...\n"); }

        AnisotropicVars_t *anisoVars = &home->anisoVars;

/*
 *      Since the elastic constants have been provided in the laboratory
 *      frame, the basis for the rotation in the anisotropy code will
 *      always be {1 0 0, 0 1 0, 0 0 1}, so initialize it here.
 */
        anisoVars->anisoRotMatrix[0][0] = 1.0;
        anisoVars->anisoRotMatrix[0][1] = 0.0;
        anisoVars->anisoRotMatrix[0][2] = 0.0;

        anisoVars->anisoRotMatrix[1][0] = 0.0;
        anisoVars->anisoRotMatrix[1][1] = 1.0;
        anisoVars->anisoRotMatrix[1][2] = 0.0;

        anisoVars->anisoRotMatrix[2][0] = 0.0;
        anisoVars->anisoRotMatrix[2][1] = 0.0;
        anisoVars->anisoRotMatrix[2][2] = 1.0;

        anisoVars->anisoRotMatrix12[0] = complex8(anisoVars->anisoRotMatrix[0][0],
                                                  anisoVars->anisoRotMatrix[1][0]);

        anisoVars->anisoRotMatrix12[1] = complex8(anisoVars->anisoRotMatrix[0][1],
                                                  anisoVars->anisoRotMatrix[1][1]);

        anisoVars->anisoRotMatrix12[2] = complex8(anisoVars->anisoRotMatrix[0][2],
                                                  anisoVars->anisoRotMatrix[1][2]);

/*
 *      Check the user provide control parameters and insure that
 *      the elastic constant matrix and elastic modulii are defined
 *      appropriately.
 */
        AnisotropicInitElasticConstants(home);

/*
 *      Allocate a bunch of temporary arrays/tables
 */
        real8 *phi              = (real8 *)malloc(numPhiPoints   * sizeof(real8));
        real8 *theta            = (real8 *)malloc(numThetaPoints * sizeof(real8));
        real8 *gaussQuadPoints  = (real8 *)malloc(numThetaPoints * sizeof(real8));
        real8 *gaussQuadWeights = (real8 *)malloc(numThetaPoints * sizeof(real8));

/*
 *      We need to set up a simple pointer to the dGdx grid as well as
 *      the full grid because at least one version of the IBM compiler on
 *      BGQ has an internal bug that causes the compiler to fail when
 *      when we try to use <dGdxgrid> as a function argument.
 */
        real8  *dGdxGridPtr     = (real8 *)calloc(1, 27* numThetaPoints * numPhiPoints * sizeof(real8));

        real8 (*dGdxgrid)[numPhiPoints][3][3][3]; ptrcpy(dGdxgrid,dGdxGridPtr);

        real8  *FqTableReal     = (real8 *)calloc(1, 18* (qMax+2)*(qMax+1)   *sizeof(real8));
        real8  *FqTableImag     = (real8 *)calloc(1, 18* (qMax+2)*(qMax+1)   *sizeof(real8));

        real8  *FqReal_v3       = (real8 *)calloc(1, 18* (qMax+2)*(qMax+1)   *sizeof(real8));
        real8  *FqImag_v3       = (real8 *)calloc(1, 18* (qMax+2)*(qMax+1)   *sizeof(real8));

        real8  *coreFqTableReal = (real8 *)calloc(1,  6*((qMax+1)*(qMax+2)-1)*sizeof(real8));
        real8  *coreFqTableImag = (real8 *)calloc(1,  6*((qMax+1)*(qMax+2)-1)*sizeof(real8));

/*
 *      Calculate the evaluation points and weights for an <numThetaPoints>
 *      Gaussian quadrature.
 */
        GaussQuadraturePointsAndWgts(numThetaPoints, gaussQuadPoints, gaussQuadWeights);

        for (i = 0; i < numPhiPoints  ; i++) { phi  [i] = (real8) i / numPhiPoints * 2 * M_PI; }
        for (i = 0; i < numThetaPoints; i++) { theta[i] = acos(gaussQuadPoints[i]); }

/*
 *      Define the derivative of the Green's function
 *      Decompose derivative of the Green's function into a grid of
 *      theta and phi angles
 */
        for (i = 0; i < numThetaPoints; i++) {
            for (j = 0; j < numPhiPoints; j++) {
                real8 lineDir[3];

                lineDir[0] = sin(theta[i]) * cos(phi[j]);
                lineDir[1] = sin(theta[i]) * sin(phi[j]);
                lineDir[2] = cos(theta[i]);

                GetGreensFuncDerivatives(GREENS_FUNC_NUM_POINTS, lineDir,
                                         anisoVars->elasticConstantMatrix4D,
                                         dGdx);

                for (int a=0; a<3; a++)
                for (int b=0; b<3; b++)
                for (int c=0; c<3; c++)
                { dGdxgrid[i][j][a][b][c] = dGdx[a][b][c]; }
            }
        }

/*
 *      Compute expansion coefficients (g^{lm}_{vpg}) of the derivative
 *      of the Green's function and compute the table of binomial
 *      factorial coefficients Q^{lm}(k)
 *
 *      NOTE: We need to set up simple pointers to the glm and
 *      binomialFactCoeffs arrays because at least one version of the
 *      IBM compiler on BGQ has an internal bug that causes the
 *      compiler to fail when when we try to use <glm> or
 *      <binomialFactCoeffs> as function arguments.
 */
        binomialFactCoeffsPtr = &binomialFactCoeffs[0][0][0];

        real8 *glmRePtr = &glmRe[0][0][0][0];
        real8 *glmImPtr = &glmIm[0][0][0][0];

        CalcGLMandBinomialFactCoeffs(qMax, numThetaPoints, numPhiPoints, phi,
                                     dGdxGridPtr, gaussQuadPoints,
                                     gaussQuadWeights, glmRePtr, glmImPtr,
                                     binomialFactCoeffsPtr);

/*
 *      Compute the initial Fq table (F = Q^{lm}(k) * g^{lm}_{vpg})
 */
        CalcFqTable(qMax, binomialFactCoeffsPtr, glmRePtr, glmImPtr,
                    FqTableReal, FqTableImag, &FqRealMax, &FqImagMax);

        // create the <Fq> tables used by the GPU-enabled forces...

        CalcFqTableV3(FqReal_v3, FqTableReal, qMax);
        CalcFqTableV3(FqImag_v3, FqTableImag, qMax);

        anisoVars->FqReal_v3 = FqReal_v3;
        anisoVars->FqImag_v3 = FqImag_v3;

/*
 *      Define B matrix: decompose B onto a grid of theta and phi angles.
 */
        real8   *BgridPtr  = (real8 *) calloc(1, numThetaPoints * numPhiPoints * 6 * sizeof(real8));
        real8  (*Bgrid)[numPhiPoints][6];  ptrcpy(Bgrid,BgridPtr);

        for (i = 0; i < numThetaPoints; i++) {
            for (j = 0; j < numPhiPoints; j++) {
                int   a;
                real8 B[6];
                real8 sphericalVecNorm[3];

                sphericalVecNorm[0] = sin(theta[i]) * cos(phi[j]);
                sphericalVecNorm[1] = sin(theta[i]) * sin(phi[j]);
                sphericalVecNorm[2] = cos(theta[i]);

                GetBmatrix(sphericalVecNorm, anisoVars->elasticConstantMatrix4D, B);

                for (a = 0; a < 6; a++) {
                    Bgrid[i][j][a] = B[a];
                }
            }
        }

/*
 *      Now calulate the inital Fq table for core forces
 *      Define B matrix: decompose B onto a grid of theta and phi angles.
 */
        CalcCoreFqTable(qMax, numThetaPoints, numPhiPoints,
                        phi, BgridPtr, gaussQuadPoints, gaussQuadWeights,
                        coreFqTableReal, coreFqTableImag);

/*
 *      Now remove the near-zero components from the Fq tables, flatten
 *      them into 1D arrays and generate the corresponding indices lists
 *      needed further on, and store the resulting arrays for later.
 *
 *      First the 'real' portion of the Fq table...
 */
        int    numElemMax0 =    (qMax+1)*(qMax+2);
        int    numElemMax  = 18*(qMax+1)*(qMax+2);
        int    numElem0    = 0;
        int    numElem     = 0;

        real8 *FqArrayReal = (real8 *)malloc(numElemMax * sizeof(real8));

        short *RemArray    = (short *)calloc(1,    (qMax+2) * sizeof(short));
        short *RelArray    = (short *)calloc(1, numElemMax0 * sizeof(short));
        short *Re18Array   = (short *)calloc(1, numElemMax0 * sizeof(short));
        short *ReiiArray   = (short *)calloc(1, numElemMax  * sizeof(short));

        FindNonZeroFqComponents(qMax, FqRealMax, FqTableReal,
                                RemArray, RelArray, Re18Array, ReiiArray,
                                FqArrayReal, &numElem, &numElem0);

        anisoVars->sphericalHarmonicsRemKept  = (short *) realloc(RemArray   , (qMax+2) * sizeof(short) );
        anisoVars->sphericalHarmonicsRelKept  = (short *) realloc(RelArray   , numElem0 * sizeof(short) );
        anisoVars->sphericalHarmonicsRe18Kept = (short *) realloc(Re18Array  , numElem0 * sizeof(short) );
        anisoVars->sphericalHarmonicsReiKept  = (short *) realloc(ReiiArray  , numElem  * sizeof(short) );
        anisoVars->FqReal                     = (real8 *) realloc(FqArrayReal, numElem  * sizeof(real8) );

        anisoVars->nz_elem_re0 = numElem0;
        anisoVars->nz_elem_re  = numElem ;

/*
 *      Now the 'imaginary' portion of the Fq table...
 */
        real8 *FqArrayImag = (real8 *)malloc(numElemMax * sizeof(real8));

        short *ImmArray    = (short *)calloc(1,    (qMax+2) * sizeof(short));
        short *ImlArray    = (short *)calloc(1, numElemMax0 * sizeof(short));
        short *Im18Array   = (short *)calloc(1, numElemMax0 * sizeof(short));
        short *ImiiArray   = (short *)calloc(1, numElemMax  * sizeof(short));

        FindNonZeroFqComponentsIm(qMax, FqImagMax, FqTableImag,
                                  ImmArray, ImlArray, Im18Array, ImiiArray,
                                  FqArrayImag, &numElem, &numElem0);

        anisoVars->sphericalHarmonicsImmKept  = (short *) realloc(ImmArray   , (qMax+2) * sizeof(short) );
        anisoVars->sphericalHarmonicsImlKept  = (short *) realloc(ImlArray   , numElem0 * sizeof(short) );
        anisoVars->sphericalHarmonicsIm18Kept = (short *) realloc(Im18Array  , numElem0 * sizeof(short) );
        anisoVars->sphericalHarmonicsImiKept  = (short *) realloc(ImiiArray  , numElem  * sizeof(short) );
        anisoVars->FqImag                     = (real8 *) realloc(FqArrayImag, numElem  * sizeof(real8) );

        anisoVars->nz_elem_im0 = numElem0;
        anisoVars->nz_elem_im  = numElem ;
        anisoVars->qMax        = qMax    ;
/*
 *      Save the <Fq> table for core forces.  These arrays will be freed
 *      at the end of the simulation.
 */
        anisoVars->coreFqReal = coreFqTableReal;
        anisoVars->coreFqImag = coreFqTableImag;

/*
 *      Now free up all the temporary arrays
 */
        free(phi);
        free(theta);
        free(gaussQuadPoints);
        free(gaussQuadWeights);
        free(dGdxgrid);
        free(FqTableReal);
        free(FqTableImag);
        free(Bgrid);

#ifdef TAYLORFMM
/*
 *      If the Taylor FMM is enabled and the user has requested anisotropic FMM,
 *      initialize the FMM anisotropy table. THIS METHOD DOES NOT WORK YET.
 */
        if ((param->fmEnabled) && (param->fmEnableAnisotropy)) {
            FMAnisotropicInit(home);
        }
#endif

        return;
}
#endif  /* end ifdef ANISOTROPIC */
