/********************************************************************************
 *
 *      Module:       MobilityLaw_Ta_pencil.c
 *		Author:       Nicolas Bertin
 *      Description:  Contains functions for calculating the twinning/anti-
 *                    twinning screw non-linear mobility based on MD simulations.
 *		              The force-velocities relations are defined explicitely in
 *		              polar coordinate system of the screw dislocation.
 *                    This mobility law uses a decomposition of the line segment
 *                    onto screw, and edge components.
 *		              The code must be linked to LAPACK when the pseudo-inverse
 *                    calculation is used to eliminate the Bline coefficient.
 *
 *      ========================================================================
 *
 *      The non-linear mobility law computes the velocity of a node
 *      solution of
 *
 *          F^{elastic} (v) + F^{drag}(v) = 0.
 *
 *      The elastic forces F^{elastic} (v) at the node are given as
 *      input. The function F^{drag}(v) is constructed using molecular
 *      dynamics input data. A Newton-Raphson procedure is used to
 *      calculate v solution of that equation.
 *
 *      Includes public functions:
 *                    MobilityLaw_Ta_pencil.cc()
 *
 ***************************************************************************/
#include <stdio.h>
#include "Home.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
#define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)

#define PSEUDO_INVERSE 0

#if 0
#include "dgesvd.h"

/**************************************************************************
 *
 *      Function:     Matrix33SVD
 *      Description:  This function computes the SVD decomposition of a
 *                    3x3 matrix using LAPACK subroutines.
 *
 *      In C order, this computes SVD of At = A'; with
 *		A' = U' * [s1 0 0 ; 0 s2 0 ; 0 0 s3] * V
 *
 *************************************************************************/
static int Matrix33SVD(double A[3][3], double U[3][3], double S[3], double V[3][3])
{
	char JOBU = 'A', JOBVT = 'A';
	double WORK[100];
	int M = 3, N = 3, LDA = 3, LDU = 3, LDVT = 3;
	int INFO = 0, LWORK = sizeof(WORK)/sizeof(*WORK);

	double AT[3][3], UT[3][3];
	Matrix33_Transpose(AT,A);

	dgesvd_(&JOBU, &JOBVT, &M, &N, AT[0], &LDA, S, UT[0],
		    &LDU, V[0], &LDVT, WORK, &LWORK, &INFO);

	Matrix33_Transpose(U,UT);

	return INFO;
}

/**************************************************************************
 *
 *      Function:     Matrix33_PseudoInverse
 *
 *************************************************************************/
static int Matrix33PseudoInverse(real8 Ainv[3][3], real8 A[3][3], int numNbrs)
{
    int stat;
    real8 U[3][3], S[3], V[3][3];
    stat = Matrix33SVD(A, U, S, V);

    real8 meps = 2.22e-16; // machine precision for double
    real8 tol = 3.0*meps*fmax(S[0], fmax(S[1], S[2]));
    int rank = 0;
    for (int i = 0; i < 3; i++)
        if (S[i] > tol) rank++;

    //if (numNbrs == 2 && rank == 3 && S[1] > 1e2*S[2]) rank--;

    if (rank == 0) {
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Ainv[i][j] = 0.0;
    } else {
        real8 Vs[3][rank];
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < rank; j++)
                Vs[i][j] = V[i][j]*1.0/S[j];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                Ainv[i][j] = 0.0;
                for (int k = 0; k < rank; k++)
                    Ainv[i][j] += Vs[i][k]*U[j][k];
            }
        }
    }

    return stat;
}
#endif


/*****************************************************************************
 *
 *      Function:    BubbleSortIndices
 *
 *****************************************************************************/
void BubbleSortIndices(int id[], double val[], int n)
{
	for (int i = 0; i < n-1; i++) {
		for (int j = 0; j < n-i-1; j++) {
			if (val[j] > val[j+1]) {
				double vtmp = val[j];
				val[j] = val[j+1];
				val[j+1] = vtmp;
				int itmp = id[j];
				id[j] = id[j+1];
				id[j+1] = itmp;
			}
		}
	}
}


/**************************************************************************
 *
 *      Function:     ScrewMDFit
 *      Description:  MD fit for (110), (112) T and AT  planes based on MD data.
 *                    This functional form fits all MD data. Parameters B, KmA, alpha
 *                    vary for each plane.
 *
 *************************************************************************/
static void ScrewMDFit(real8 burgMag, int TATdir, real8 v, real8 *f, real8 *dfdv)
{
	real8 B, KmA, alpha;
	real8 Bv, f0, expmBv, Fac, Fac2;

	v = v * burgMag;

	if (TATdir == 0) {
		// Twinning
		KmA = 2200;
		B = 0.12;
		alpha = 1.7;
	} else {
		// Anti-twinning
		KmA = 10600;
		B = 0.5;
		alpha = 3.0;
	}

    Bv = B*v;
    expmBv = exp(-Bv);
    Fac = (1 + expmBv);

    Fac2 = KmA*B;
    f0 = KmA / 2.0;

    *f = KmA / Fac - f0;
    *f +=  alpha * v;

    *dfdv = Fac2 * expmBv / pow(Fac,2) + alpha;

    // Non-linear velocity at high forces
    real8 ea, eb;
    ea = 0.012;
    eb = 800.0;
    *f += exp(ea*(v-eb));
    *dfdv += ea*exp(ea*(v-eb));

    // Avoid divergence in NR
    real8 fmax = 10000.0;
    if (*f > fmax) {
        real8 fm, dfmdv;
        real8 vmax = 0.0;
        for (int iter = 0; iter < 1000; iter++) {
            expmBv = exp(-B*vmax);
            fm = alpha*vmax + KmA/(1+expmBv) - f0 + exp(ea*(vmax-eb)) - fmax;
            dfmdv = alpha + KmA*B*expmBv/(1.0+expmBv)/(1.0+expmBv) + ea*exp(ea*(vmax-eb));
            vmax -= fm/dfmdv;
            if (fabs(fm) < 1e-6) break;
        }
        *dfdv = 100.0*alpha;
        *f = fmax + (*dfdv) * (v-vmax);
    }

    *f = (*f) * 1e6;  // f from MPa to Pa
    *dfdv = (*dfdv) * 1e6 * burgMag;  // df from MPa / (m/s) to Pa / (|b|/s)
}


/*****************************************************************************
 *
 *      Function:    ForceScrewPencil
 *      Description: This function computes force required to move a screw
 *                   dislocation of Burgers vector b with velocity v.
 *                   The input velocity vector is assumed to be perpendicular
 *                   to the line direction (and b) and the output force vector
 *                   is also perpendicular to b.
 *
 *   Inputs:
 *
 *       v   Screw velocity, in [|b|/s]
 *       b   The Burgers vector, in [arbitrary units]
 *     ksi   The dislocation line (tangent) vector, in [arbitrary units]
 *
 *   Output:
 *
 *       f   PK force per unit length of the screw dislocation, in [Pa*b]
 *      dfdv Derivative of PK force wrt velocity per unit length of the screw
 *           dislocation, in [Pa*s]
 *
 *****************************************************************************/
static void ForceScrewPencil(real8 burgMag, real8 vel[3], real8 b[3],
                             real8 ksi[3], real8 Fscrew[3], real8 dFscrewdv[3][3])
{
	int   id[3];
	real8 blen, klen, Coeff;
	real8 t10[3], t20[3], t30[3];
	real8 tt[3][3], vt[3], t[3], at[3];
	real8 dx[3], dy[3], dz[3];
	real8 dproj, vproj[3];
	real8 psi, vmag, chi, fmag, k, A, u;
	real8 fT, fAT, dfTdv, dfATdv, F[2];
	real8 dchidpsi, dchidv, dfdpsi, dfdv, dudpsi;
	real8 dpsidx, dpsidy, dvdx, dvdy;
	real8 dXdD[3][2], dVdX[2][3], dF[2][2], dfdV[2][2], dFdV[3][2];

	// Define three twinning motion directions perpendicular to b
	blen = Normal(b);
	klen = sqrt(DotProduct(ksi, ksi));
	Coeff = 1.0/blen/sqrt(2);

	t10[X] = 2*Coeff*b[X];
	t10[Y] =  -Coeff*b[Y];
	t10[Z] =  -Coeff*b[Z];

	t20[X] =  -Coeff*b[X];
	t20[Y] = 2*Coeff*b[Y];
	t20[Z] =  -Coeff*b[Z];

	t30[X] =  -Coeff*b[X];
	t30[Y] =  -Coeff*b[Y];
	t30[Z] = 2*Coeff*b[Z];

	cross(t10, ksi, tt[0]);
	cross(t20, ksi, tt[1]);
	cross(t30, ksi, tt[2]);

	for (int i = 0; i < 3; i++) {
		for (int j = 0; j < 3; j++)
        	tt[i][j] /= klen;
		id[i] = i;
		vt[i] = DotProduct(vel, tt[i]);
	}

	// Find which of three T vectors is most aligned with the velocity
	// The bounding AT direction is obtained by summing the 1st and the 2nd most aligned tt
	BubbleSortIndices(id, vt, 3);

	for (int i = 0; i < 3; i++) {
		t[i] = tt[id[2]][i];
	       at[i] = tt[id[2]][i] + tt[id[1]][i];
	}

	// (110) center plane and psi angle
	dx[X] = t[X] + at[X];
	dx[Y] = t[Y] + at[Y];
	dx[Z] = t[Z] + at[Z];
	NormalizeVec(dx);
	NormalizedCrossVector(dx, at, dz);
	NormalizedCrossVector(dz, dx, dy);


	// The rest uses force-velocity functions defined in polar coordinates
	psi = atan2(DotProduct(vel, dy), DotProduct(vel, dx));
	psi = fmax(fmin(psi, M_PI/6.0), -M_PI/6.0);

	dproj = DotProduct(vel, dz);
	vproj[X] = vel[X] - dproj*dz[X];
	vproj[Y] = vel[Y] - dproj*dz[Y];
	vproj[Z] = vel[Z] - dproj*dz[Z];

	vmag = Normal(vproj);
	if (vmag < 1e-10) {
		vmag = 1.0/burgMag;
		//psi = 0.0;
	}

	//printf("vmag = %e\n",vmag);


	k = 15.0;
	A = M_PI/3.0 + 0.00081;

	chi = -1.0/k*log(A/(psi+A/2.0)-1.0);
	chi = fmax(fmin(chi, M_PI/6.0), -M_PI/6.0);

	dchidpsi = 4.0*A/(A*A*k-4.0*k*psi*psi);
	dchidv = 0.0;

	ScrewMDFit(burgMag, 0, vmag, &fT,  &dfTdv);
	ScrewMDFit(burgMag, 1, vmag, &fAT, &dfATdv);

	u = sin(psi)+0.5;
	fmag = (1.0-u)*fT + u*fAT;

	dudpsi = cos(psi);
	dfdpsi = -dudpsi*fT + dudpsi*fAT;
	dfdv = (1.0-u)*dfTdv + u*dfATdv;

	// 2D polar to 3D Cartesian transforms
	dpsidx = -1.0/vmag*sin(psi);
	dpsidy = 1.0/vmag*cos(psi);
	dvdx = cos(psi);
	dvdy = sin(psi);

	for (int i = 0; i < 3; i ++) {
		dXdD[i][0] = dx[i];
		dXdD[i][1] = dy[i];
	}

	F[0] = fmag*cos(chi);
	F[1] = fmag*sin(chi);
	Matrix32_Vmul(Fscrew, dXdD, F );

	dF[0][0] = -fmag*sin(chi);
	dF[0][1] = cos(chi);
	dF[1][0] = fmag*cos(chi);
	dF[1][1] = sin(chi);


	dVdX[0][0] = dpsidx*dXdD[0][0] + dpsidy*dXdD[0][1];
	dVdX[0][1] = dpsidx*dXdD[1][0] + dpsidy*dXdD[1][1];
	dVdX[0][2] = dpsidx*dXdD[2][0] + dpsidy*dXdD[2][1];
	dVdX[1][0] = dvdx*dXdD[0][0] + dvdy*dXdD[0][1];
	dVdX[1][1] = dvdx*dXdD[1][0] + dvdy*dXdD[1][1];
	dVdX[1][2] = dvdx*dXdD[2][0] + dvdy*dXdD[2][1];

	dfdV[0][0] = dF[0][0]*dchidpsi + dF[0][1]*dfdpsi;
	dfdV[0][1] = dF[0][0]*dchidv   + dF[0][1]*dfdv;
	dfdV[1][0] = dF[1][0]*dchidpsi + dF[1][1]*dfdpsi;
	dfdV[1][1] = dF[1][0]*dchidv   + dF[1][1]*dfdv;

	dFdV[0][0] = dXdD[0][0]*dfdV[0][0] + dXdD[0][1]*dfdV[1][0];
	dFdV[0][1] = dXdD[0][0]*dfdV[0][1] + dXdD[0][1]*dfdV[1][1];
	dFdV[1][0] = dXdD[1][0]*dfdV[0][0] + dXdD[1][1]*dfdV[1][0];
	dFdV[1][1] = dXdD[1][0]*dfdV[0][1] + dXdD[1][1]*dfdV[1][1];
	dFdV[2][0] = dXdD[2][0]*dfdV[0][0] + dXdD[2][1]*dfdV[1][0];
	dFdV[2][1] = dXdD[2][0]*dfdV[0][1] + dXdD[2][1]*dfdV[1][1];

	dFscrewdv[0][0] = dFdV[0][0]*dVdX[0][0] + dFdV[0][1]*dVdX[1][0];
	dFscrewdv[0][1] = dFdV[0][0]*dVdX[0][1] + dFdV[0][1]*dVdX[1][1];
	dFscrewdv[0][2] = dFdV[0][0]*dVdX[0][2] + dFdV[0][1]*dVdX[1][2];
	dFscrewdv[1][0] = dFdV[1][0]*dVdX[0][0] + dFdV[1][1]*dVdX[1][0];
	dFscrewdv[1][1] = dFdV[1][0]*dVdX[0][1] + dFdV[1][1]*dVdX[1][1];
	dFscrewdv[1][2] = dFdV[1][0]*dVdX[0][2] + dFdV[1][1]*dVdX[1][2];
	dFscrewdv[2][0] = dFdV[2][0]*dVdX[0][0] + dFdV[2][1]*dVdX[1][0];
	dFscrewdv[2][1] = dFdV[2][0]*dVdX[0][1] + dFdV[2][1]*dVdX[1][1];
	dFscrewdv[2][2] = dFdV[2][0]*dVdX[0][2] + dFdV[2][1]*dVdX[1][2];

	return;
}

/**************************************************************************
 *
 *      Function:     ScrewDragPencil
 *      Description:  This function returns the non-linear mobility of a screw
 *                    dislocation in tantalum based on MD simulations
 *
 *************************************************************************/
static void ScrewDragPencil(real8 burgMag, real8 vel[3], real8 burg[3], real8 t[3],
                            real8 Bline, real8 f[3],real8 dfdv[3][3])
{
    int   m, n;
    real8 lineDir[3];
    real8 lineDirMat[3][3];
    real8 tdotb;
    real8 Fscrew[3], dFscrewdv[3][3];

    // Project the line direction along the Burgers vector
    // while keeping the sign of the line direction.
    tdotb = DotProduct(t,burg);
    lineDir[X] = tdotb * burg[X];
    lineDir[Y] = tdotb * burg[Y];
    lineDir[Z] = tdotb * burg[Z];
    NormalizeVec(lineDir);
    Vec3TransposeAndMult(lineDir,  lineDirMat);

	ForceScrewPencil(burgMag, vel, burg, lineDir, Fscrew, dFscrewdv);

    // First add the linear line drag component.
    // On the derivative
    for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
            dfdv[m][n] =  Bline * lineDirMat[m][n];

    // On the force
    Matrix33Vector3Multiply(dfdv, vel, f);

    // Second add the new non-linear mobility law with components
    // along glide and climb.
    // On the force
    for (m = 0; m < 3; m++) f[m] += Fscrew[m];

    // On the derivative
    for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
            dfdv[m][n] += dFscrewdv[m][n];

    return;
}


/**************************************************************************
 *
 *      Function:     EdgeMDFit
 *
 *************************************************************************/
static void EdgeMDFit(real8 burgMag, real8 v, real8 *f, real8 *dfdv)
{
    v = v * burgMag;

    real8 vsign = (v > 0.0) ? 1.0 : -1.0;
    v = fabs(v); // velocity magnitude

    real8 A = 0.9196;
    real8 B = 22.2860;
    real8 C = 1.1050e3;
    real8 D = 28.1015;

    real8 fmax = 5000.0;
    real8 vmax = 1.3265e3;

    *f = A*v + B*pow(v/C, D);
    *dfdv = A + B*D/C*pow(v/C, D-1);

    // Avoid divergence in NR
    if (*f > fmax) {
        *dfdv = 100.0*A;
        *f = fmax + (*dfdv) * (v-vmax);
    }

    *f = vsign * (*f) * 1e6;  // f from MPa to Pa
    *dfdv = (*dfdv) * 1e6 * burgMag;  // df from MPa / (m/s) to Pa / (|b|/s)
}


/**************************************************************************
 *
 *      Function:     TAEdgeDragNonLinear
 *      Description:  This function returns the drag function g_edge
 *                    along the line, the glide directions and the climb
 *                    directions as well as its derivative Dg_edge/Dv.
 *                    EdgeDrag is non-linear fitted from MD data.
 *
 *************************************************************************/
static void TAEdgeDragNonLinear(real8 burgMag, real8 vel[3], real8 burg[3],
                                real8 lineDir[3], real8 dragClimb, real8 dragLine,
                                real8 fEdgeDrag[3], real8 dfEdgeDragdv[3][3])
{
    int   m, n;
    real8 climbDir[3], glideDir[3];
    real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];

    // glide direction is the Burgers vector normalized.
    VECTOR_COPY(glideDir, burg);
    NormalizeVec(glideDir);

    NormalizedCrossVector(lineDir, glideDir, climbDir);

    real8 projv, projf=0, dprojfdv=0;
    projv = DotProduct(vel, glideDir);

    EdgeMDFit(burgMag, projv, &projf, &dprojfdv);

    Vec3TransposeAndMult(glideDir, glideDirMat);
    Vec3TransposeAndMult(climbDir, climbDirMat);
    Vec3TransposeAndMult( lineDir,  lineDirMat);

    for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
            dfEdgeDragdv[m][n] = dragClimb * climbDirMat[m][n] +
                                 dragLine  * lineDirMat [m][n];

    Matrix33Vector3Multiply(dfEdgeDragdv, vel, fEdgeDrag);

    for (m = 0; m < 3; m++)
        for (n = 0; n < 3; n++)
            dfEdgeDragdv[m][n] += dprojfdv * glideDirMat[m][n];

    for (m = 0; m < 3; m++)
        fEdgeDrag[m] += projf * glideDir[m];

    return;
}


/**************************************************************************
 *
 *      Function:     ForceErrorCalc
 *
 *************************************************************************/
static void ForceErrorCalc(Param_t *param, real8 vTest[3], real8 fn[3], int numNbrs,
                           int eEx[MAX_NBRS], int sEx[MAX_NBRS], int jctEx[MAX_NBRS],
                           real8 Ltot[MAX_NBRS],real8 Lscrew[MAX_NBRS],real8 Ledge[MAX_NBRS],
                           real8 lineDir[MAX_NBRS][3], real8 burg[MAX_NBRS][3], real8 edgeDir[MAX_NBRS][3],
                           real8 burgMag, real8 dragLine, real8 dragClimb, real8 dragEdge,
                           real8 fError[3], real8 dfErrordv[3][3])
{
    int m,n,i;
    real8 eps = 1.0e-12;


    for (i = 0; i < numNbrs; i++)
    {
        real8 fJunction[3], dfJunctiondv[3][3];
        real8 fScrew[3], dfScrewdv[3][3];
        real8 fEdge[3], dfEdgedv[3][3];

        if (fabs(Ltot[i]) < eps) continue;

/*
 *      Junction contribution
 */
        if (jctEx[i])
        {
			// Effective climb component from vacancy concentration
			real8 invbeMag = 1.0;
			real8 dragClimbEff = invbeMag*invbeMag*dragClimb;

			TAJunctionDrag(vTest, lineDir[i], dragLine,
                           dragClimbEff, fJunction, dfJunctiondv);

            for (m = 0; m < 3; m++) {
                fError[m] -= 0.5 * Ltot[i] * fJunction[m];
                for (n = 0; n < 3; n++) {
                    dfErrordv[m][n] -= 0.5 * Ltot[i] * dfJunctiondv[m][n];
                }
            }

        }
        else
        {  /* not a junction */
/*
 *          Screw contribution
 */
            if (sEx[i])
            {
                ScrewDragPencil(burgMag, vTest, burg[i], lineDir[i], dragLine,
                                fScrew, dfScrewdv);

                for (m = 0; m < 3; m++) {
                    fError[m] -= 0.5 * Lscrew[i] * fScrew[m];
                    for (n = 0; n < 3; n++) {
                        dfErrordv[m][n] -= 0.5 * Lscrew[i] * dfScrewdv[m][n];
                    }
                }

            }  /* end if (sEx[i]) */

/*
 *                  Edge contribution
 */
            if (eEx[i])
            {
                TAEdgeDragNonLinear(burgMag, vTest, burg[i],
                                    edgeDir[i], dragClimb, dragLine,
                                    fEdge, dfEdgedv);

                for (m = 0; m < 3; m++) {
                    fError[m] -= 0.5 * Ledge[i] * fEdge[m];
                    for (n = 0; n < 3; n++) {
                        dfErrordv[m][n] -= 0.5 * Ledge[i] * dfEdgedv[m][n];
                    }
                }

            }  /* end if (eEx[i]) */

        }

    } /* end loop over segments */

}

/**************************************************************************
 *
 *      Function:     MobilityLaw_Ta_pencil.cc
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *************************************************************************/
int Mobility_Ta_pencil(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int   i, m, n;
        int   count, TotalCount;
        int   numNbrs;
        int   numNonZeroLenSegs = 0;

        int   edgeExists[MAX_NBRS];
        int   screwExists[MAX_NBRS];
        int   junctionExists[MAX_NBRS];

        real8 mu, maxSegLen, burgMag, totLen;
        real8 massMult, massMatrix[3][3], eps, adjustment = 1;
        real8 vTest[3], vOld[3], fn[3], rt[3], fmag;

        real8 segTotLength[MAX_NBRS],segScrewLength[MAX_NBRS];
        real8 segEdgeLength[MAX_NBRS];
        real8 burg[MAX_NBRS][3], lineDir[MAX_NBRS][3], edgeDir[MAX_NBRS][3];

        real8 dragEdge, dragLine, dragClimb;

        real8 fError[3], dfErrordv[3][3];
        real8 fErrorTol, forceError, forceErrorOld;
        real8 correction[3], inertia[3];
        int exitNR;
        Param_t *param;
        Node_t  *nbrNode;

        int forceisgood = 0;
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

        param = home->param;

        numNbrs = node->numNbrs;

        mu = param->shearModulus;
        maxSegLen = param->maxSeg;

        burgMag = param->burgMag;

        real8 (*invdfErrordv)[3]     =  mobArgs->invDragMatrix;
        real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdfErrordv);
        *numGlideConstraints = 0;

        eps = 1.0e-12;

        // Drag edge is (\sigma . |b| ) / |v| in units of GPa * |b| / (m/s)= GPa * |b| * s / m.
        // dragEdge converted to Pa * s
        // burgMag is such that : b = sqrt(3)/2 * (2*V)^(1./3.). V= 17.8806 A^3;
        // burgMag = 2.853217687613880
        dragEdge   = 0.9e6*burgMag;
        //dragClimb  = dragEdge * 1.0e5;

	dragClimb  = 1.0 / param->MobClimb;

#if 0
	/* Climb mobility for junction nodes */
        if (param->MobClimbJunc > 0.0) {
	  if (node->numNbrs > 2)
	    dragClimb = 1.0 / param->MobClimbJunc;
        }
#else
	/* Use BCC_0b for junction nodes */
	if (param->MobClimbJunc > 0.0 && node->numNbrs > 2) {
	  Mobility_BCC_0b(home, node, mobArgs);
	  return(0);
	}
#endif


#if PSEUDO_INVERSE
        dragLine = 0.0;
#else
        //dragLine = dragEdge * 1.0e-2;
        dragLine = 2.0e-5;
#endif

/*
 *      If we need to include inertial terms, need some more initializations
 */
        if (param->includeInertia != 0)
        {
           massMult = 0.25 * param->massDensity *
              (burgMag * burgMag);
        }
        else
           massMult = 0.0;

/*
 *      This function saves various pieces of info for each neighbor
 *      node using static arrays for speed rather than doing memory
 *      allocations.  If the arrays are too  small, abort with an error
 *      to let the caller increase the array sizes.
 */
        if (numNbrs > MAX_NBRS) {
            Fatal("Mobility_Ta_pencil: Node segment count (%d) exceeds\n"
                  "maximum allowed (%d) in static arrays.  Increase MAX_NBRS\n"
                  "and recompile\n", node->numNbrs, MAX_NBRS);
        }

/*
 *      Loop over all arms attached to the node to get some information
 *      for the Newton-Raphson procedure
 */
        totLen = 0.0;

        for (i = 0; i < numNbrs; i++) {
            real8 rtSquared;

            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL) {
                printf("WARNING: Neighbor not found at %s line %d\n",
                       __FILE__, __LINE__);
                continue;
            }

/*
 *          Get the segment line
 */
            rt[X] = nbrNode->x - node->x;
            rt[Y] = nbrNode->y - node->y;
            rt[Z] = nbrNode->z - node->z;

            ZImage(param, &rt[X], &rt[Y], &rt[Z]);

/*
 *          If the segment is zero length (which may happen when the
 *          mobility function is called from SplitMultiNodes()), just
 *          skip the segment.
 */
            if ((rtSquared = DotProduct(rt, rt)) < eps) {
                segTotLength[i] = 0.0;
                continue;
            }

            numNonZeroLenSegs++;

/*
 *          Save the burgers vector
 */
            burg[i][X] = node->burgX[i];
            burg[i][Y] = node->burgY[i];
            burg[i][Z] = node->burgZ[i];

/*
 *          If necessary, rotate the burgers vector and line sense from the
 *          laboratory frame to the crystal frame
 */
            if (param->useLabFrame) {
               real8 burgRot[3], rtRot[3];

                Matrix33Vector3Multiply(param->rotMatrixInverse, burg[i],
                                        burgRot);
                Matrix33Vector3Multiply(param->rotMatrixInverse, rt,
                                        rtRot);

                VECTOR_COPY(burg[i], burgRot);
                VECTOR_COPY(rt, rtRot);
            }

/*
 *          Calculate the segment total length
 */
            segTotLength[i] = sqrt(rtSquared);
            totLen += segTotLength[i];

            int bIndex = -1;
            real8 (*burgList)[3];
            burgList = home->burgData.burgList;
            bIndex = GetBurgIndexBCC(burg[i], BCC_NUM_TOT_BURG, burgList);

            if (bIndex >= BCC_NUM_GLIDE_BURG)
            {
/*
 *              Junction
 */
                junctionExists[i] = 1;
                edgeExists[i]     = 0;
                screwExists[i]    = 0;
            }
            else
            {
/*
 *              Not a junction
 */
                real8 L2, nburg, intFac,ls;
                junctionExists[i] = 0;

                nburg = Normal(burg[i]);

                ls = DotProduct(rt, burg[i])/nburg;

                segScrewLength[i] = fabs(ls);
                segScrewLength[i] = MIN(segScrewLength[i],segTotLength[i]);

                L2 = MAX(0.0, (segTotLength[i]*segTotLength[i]-segScrewLength[i]*segScrewLength[i]));

                if (L2 <= 0.0)
                {
                    segEdgeLength[i] = 0.0;
                    VECTOR_ZERO(edgeDir[i]);
                }
                else
                {
                    segEdgeLength[i] = sqrt(L2);
                    intFac = ls/nburg;
                    edgeDir[i][X] = (rt[X] - intFac*burg[i][X])/segEdgeLength[i];
                    edgeDir[i][Y] = (rt[Y] - intFac*burg[i][Y])/segEdgeLength[i];
                    edgeDir[i][Z] = (rt[Z] - intFac*burg[i][Z])/segEdgeLength[i];
                }

                edgeExists[i]  = ((segEdgeLength[i]/segTotLength[i]) > 1.0e-04);
                screwExists[i] = ((segScrewLength[i]/segTotLength[i])> 1.0e-04);

            } /* end glide Burgers vector */

/*
 *          Calculate tangents
 */
            lineDir[i][X] = rt[X] / segTotLength[i];
            lineDir[i][Y] = rt[Y] / segTotLength[i];
            lineDir[i][Z] = rt[Z] / segTotLength[i];

        }  /* loop over neighbors i */

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

/*
 *      Setup for the non-linear Newton-Raphson loop
 *      vTest has been initialized to the linear solution.
 */
        count = 0;
        TotalCount = 5000;

        // ParaDiS force in Pa.
        fn[0] = node->fX;
        fn[1] = node->fY;
        fn[2] = node->fZ;

        // Initial velocity guess
        vTest[0] = node->vX;
    	vTest[1] = node->vY;
        vTest[2] = node->vZ;

        VECTOR_COPY(vOld, vTest);

        forceisgood = 0;

        fErrorTol = MAX((DotProduct(fn, fn) * 1.0e-16),
                        (mu * mu * maxSegLen * maxSegLen * 1.0e-24));

        fmag = sqrt(DotProduct(fn, fn));

/*
 *      If necessary, rotate the force and velocity vectors from the
 *      laboratory frame to the crystal frame
 */
        if (param->useLabFrame) {
            real8 fnRot[3], vTestRot[3];
            Matrix33Vector3Multiply(param->rotMatrixInverse, fn, fnRot);
            Matrix33Vector3Multiply(param->rotMatrixInverse, vTest, vTestRot);
            VECTOR_COPY(fn, fnRot);
            VECTOR_COPY(vTest, vTestRot);
            VECTOR_COPY(vOld, vTest);
        }

/*
 *      If we need to include inertial terms...
 */
        if (param->includeInertia) {
            real8 vDiff[3];
            vDiff[0] = vTest[0] - vOld[0];
            vDiff[1] = vTest[1] - vOld[1];
            vDiff[2] = vTest[2] - vOld[2];
            Matrix33Vector3Multiply(massMatrix, vDiff, inertia);
        } else {
            VECTOR_ZERO(inertia);
        }

        forceError = fErrorTol + 1;
        forceErrorOld = forceError;
        exitNR = 0;

/*
 *      Begin the Newton-Raphson loop
 */
        while ((count <= TotalCount) && (forceError > fErrorTol) && (!exitNR)) {


           if (count < 2) forceErrorOld = forceError;

           for (m = 0; m < 3; m++) fError[m] = fn[m] - inertia[m];

           for (m = 0; m < 3; m++)
              for (n = 0; n < 3; n++)
                 dfErrordv[m][n] = -massMatrix[m][n];

/*
 *          Loop over all segments (skipping zero length segments)
 *          adding up screw, edge and junction contributions to fError
 *          and dfErrordv.
 */
           ForceErrorCalc( param, vTest,  fn,  numNbrs,
                           edgeExists,    screwExists,    junctionExists,
                           segTotLength, segScrewLength,  segEdgeLength,
                           lineDir,  burg,  edgeDir,
                           burgMag,  dragLine, dragClimb,  dragEdge,
                           fError,  dfErrordv);

#if PSEUDO_INVERSE
            if (Matrix33_PseudoInverse(invdfErrordv, dfErrordv) != 0) {
                Fatal("%s::%s(%d) : Cannot pseudoinvert dfErrordv!", __FILE__, __func__, __LINE__ );
            }
#else
            if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 ) 
            {
                printf("det of derivative = %g\n", Matrix33_Det(dfErrordv) );
                Print3x3("dfErrordv", dfErrordv);
                Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
            }
#endif

            Matrix33Vector3Multiply(invdfErrordv, fError, correction);
            forceError = DotProduct(fError, fError);

            vTest[0] -= adjustment*correction[0];
            vTest[1] -= adjustment*correction[1];
            vTest[2] -= adjustment*correction[2];

#if PSEUDO_INVERSE
            if (fabs(forceError-forceErrorOld)/fmag < 1e-6 || adjustment < 1e-20) exitNR = 1;
#endif

            // This forceError is good. Continue Newton
            if (forceError <= forceErrorOld + 1e-6)
            {
               forceErrorOld = forceError;
               if (forceisgood > 100)
               {
                  adjustment = MIN(adjustment*2.0,1.0);
                  forceisgood = 0;
                  //printf("forceisgood = %d force is good adjustment = %g\n",forceisgood, adjustment);
               }
               forceisgood += 1;
            }

            // If the current step has led to a worst error than before. Modify the Newton step.
            while (forceError > forceErrorOld && count > 2 && adjustment > 0.0156 /*0.1*adjustmentOld*/)
              {
                 forceisgood = 0;
                 // Adjust Newton step
                 adjustment *= 0.5;

                 // Re-calculate error
                 vTest[0] -= adjustment*correction[0];
                 vTest[1] -= adjustment*correction[1];
                 vTest[2] -= adjustment*correction[2];

                 for (m = 0; m < 3; m++) fError[m] = fn[m] - inertia[m];

                 for (m = 0; m < 3; m++)
                    for (n = 0; n < 3; n++)
                       dfErrordv[m][n] = -massMatrix[m][n];

                 ForceErrorCalc( param, vTest,  fn,  numNbrs,
                                 edgeExists,    screwExists,    junctionExists,
                                 segTotLength, segScrewLength,  segEdgeLength,
                                 lineDir,  burg,  edgeDir,
                                 burgMag,  dragLine, dragClimb,  dragEdge,
                                 fError,  dfErrordv);

#if PSEUDO_INVERSE
                if (Matrix33_PseudoInverse(invdfErrordv, dfErrordv) != 0) {
                    Fatal("%s::%s(%d) : Cannot pseudoinvert dfErrordv!", __FILE__, __func__, __LINE__ );
                }
#else
                 if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 ) 
                 {
                    printf("det of derivative = %g\n", Matrix33_Det(dfErrordv) );
                    Print3x3("dfErrordv", dfErrordv);
                    Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
                 }
#endif

                 Matrix33Vector3Multiply(invdfErrordv, fError, correction);
                 forceError = DotProduct(fError, fError);
              }

              count++;
        } /* end Newton-Raphson while */



#if 0 // !PSEUDO_INVERSE
        if ((count > TotalCount) || (forceError > fErrorTol)) {
            printf("WARNING: MobilityLaw_Ta_pencil.cc has not converged for node id %d in %d iterations\n",
                   node->myTag.index,count);
        }
#endif


/*
 *      We were able to converge on a velocity
 *
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame.
 */
        if (param->useLabFrame) {
            real8 vTestRot[3];
            Matrix33Vector3Multiply(param->rotMatrix, vTest, vTestRot);
            VECTOR_COPY(vTest, vTestRot);
        }

        node->vX = vTest[0];
        node->vY = vTest[1];
        node->vZ = vTest[2];

        return(0);
}
