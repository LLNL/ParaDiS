/**************************************************************************
 *
 *      Module:       MobilityLaw_Ta.cc
 *		Author:       Nicolas Bertin, based on Sylvie Aubry's MobilityLaw_Ta.cc function
 *      Description:  BCC mobility law with non-linear functions
 *                    for edge and screw dislocations in tantalum
 * 		              fitted from MD data.
 *		              The code must be linked to LAPACK when the pseudo-inverse
 *                    calculation is used to eliminate the Bline coefficient.
 *
 *      ==============================================================
 *
 ***************************************************************************/
#include <stdio.h>
#include <math.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Mobility.h"

#ifdef _ARLFEM
#include "FEM.h"
#endif

#define NON_LINEAR_EDGE 1
#define PSEUDO_INVERSE  0

//#define M33_INVERSE(a,b) Matrix33_Inverse(a,b)
//#define M33_INVERSE(a,b) Matrix33_SVD_Inverse(a,b)
  #define M33_INVERSE(a,b) Matrix33_PseudoInverse(a,b)

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
 *      Function:     Matrix33PseudoInverse
 *
 *************************************************************************/
static int Matrix33PseudoInverse(real8 A[3][3], real8 Ainv[3][3])
{
    int stat;
    real8 U[3][3], S[3], V[3][3];
    stat = Matrix33SVD(A, U, S, V);

    real8 meps = 2.22e-16; // machine precision for double
    real8 tol = 3.0*meps*fmax(S[0], fmax(S[1], S[2]));
    int rank = 0;
    for (int i = 0; i < 3; i++)
        if (S[i] > tol) rank++;

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


/**************************************************************************
 *
 *      Function:     PrintNodeDebug
 *
 *************************************************************************/
static void PrintNodeDebug(Home_t *home, Node_t *node)
{
    printf("Node %d\n",node->myTag.index);
    for (int i = 0; i < node->numNbrs; i++) {
        real8 rt[3];
        Node_t *nbrNode = GetNeighborNode(home, node, i);
        rt[X] = nbrNode->x - node->x;
        rt[Y] = nbrNode->y - node->y;
        rt[Z] = nbrNode->z - node->z;
        ZImage(home->param, &rt[X], &rt[Y], &rt[Z]);
        printf("  arm %d l = %e %e %e\n",i,rt[0],rt[1],rt[2]);
        printf("         b = %e %e %e\n",node->burgX[i],node->burgY[i],node->burgZ[i]);
    }
    printf("  f   = %e %e %e\n",node->fX,node->fY,node->fZ);
}


/**************************************************************************
 *
 *      Function:     ScrewDrag
 *      Description:  This function returns the drag function g_screw
 *                    along the line, glide and climb directions as well
 *                    as its derivative Dg_screw/Dv.  ScrewDrag is linear.
 *
 *************************************************************************/
static void LinScrewDrag(real8 vel[3], real8 ksi[3], real8 burg[3],
                         real8 climbDir[3], real8 glideDir[3],
                         real8 dragLine, real8 dragClimb,
                         real8 B112T, real8 B112AT, real8 B110,
                         real8 fScrewDrag[3], real8 dfScrewDragdv[3][3])
{
        int   m, n;
        real8 dragGlide = 0.0;
        real8 lineDir[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];

	// For screws, line direction is the Burgers vector normalized.
        lineDir[X] = burg[X];
        lineDir[Y] = burg[Y];
        lineDir[Z] = burg[Z];
        NormalizeVec(lineDir);

        if (fabs(climbDir[X]*climbDir[Y]*climbDir[Z]) < 1e-10)
        {
           // {110} plane
           dragGlide = B110;
        }
         else
        {
           // {112} plane
           real8 tdotb, projv;
           projv = DotProduct(vel,glideDir);
           tdotb = DotProduct(ksi,burg);

           if ( (tdotb >= 0 && projv >= 0.0) || (tdotb <0 && projv <0) )
           {
              // twinning
              dragGlide = B112T;
           }
           else
           {
              // anti-twinning
              dragGlide = B112AT;
           }
        }

        Vec3TransposeAndMult(glideDir, glideDirMat);
        Vec3TransposeAndMult(climbDir, climbDirMat);
        Vec3TransposeAndMult(lineDir,  lineDirMat);

        for (m = 0; m < 3; m++)
           for (n = 0; n < 3; n++)
              dfScrewDragdv[m][n] = dragGlide * glideDirMat[m][n] +
                                    dragClimb * climbDirMat[m][n] +
                                    dragLine  * lineDirMat [m][n];

        Matrix33Vector3Multiply(dfScrewDragdv, vel, fScrewDrag);
        return;
}


/**************************************************************************
 *
 *      Function:     ScrewMDFit
 *      Description:  MD fit for (110), (112) T and AT  planes based on MD data.
 *                    This functional form fits all MD data. Parameters B, KmA, alpha
 *                    vary for each plane.
 *
 *************************************************************************/
static void ScrewMDFit(real8 burgMag, real8 B, real8 KmA, real8 alpha,
                       real8 v, real8 *f, real8 *dfdv)
{
    v = v * burgMag;

    real8 vsign = (v > 0.0) ? 1.0 : -1.0;
    v = fabs(v); // velocity magnitude

    real8 Bv, f0, expmBv,Fac,Fac2;
    Bv = B*v;

    expmBv = exp(-Bv);
    Fac = (1 + expmBv);

    Fac2 = KmA*B;
    f0 = KmA / 2.0;

    *f = KmA / Fac - f0;
    *f +=  alpha * v;

    *dfdv = Fac2 * expmBv / pow(Fac,2) + alpha;

#if NON_LINEAR_EDGE
    // Non-linear velocity at high forces
    real8 ea, eb;
    ea = 0.005;
    eb = 1000.0;
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
#endif

    *f = vsign * (*f) * 1e6;  // f from MPa to Pa
    *dfdv = (*dfdv) * 1e6 * burgMag;  // df from MPa / (m/s) to Pa / (|b|/s)
}


/**************************************************************************
 *
 *      Function:     ScrewDrag
 *      Description:  This function returns the drag function g_screw
 *                    along the line, glide and climb directions as well
 *                    as its derivative Dg_screw/Dv.  ScrewDrag is linear.
 *
 *************************************************************************/
static void ScrewDrag(real8 burgMag, real8 vel[3], real8 ksi[3], real8 burg[3],
                      real8 climbDir[3], real8 glideDir[3],
                      real8 dragLine, real8 dragClimb,
                      real8 K110, real8 B110, real8 a110,
                      real8 K112T, real8 B112T, real8 a112T,
                      real8 K112AT, real8 B112AT, real8 a112AT,
                      real8 fScrewDrag[3], real8 dfScrewDragdv[3][3])
{
        int   m, n;
        real8 lineDir[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];

	// For screws, line direction is the Burgers vector normalized.
        lineDir[X] = burg[X];
        lineDir[Y] = burg[Y];
        lineDir[Z] = burg[Z];
        NormalizeVec(lineDir);

        real8 projv, projf=0, dprojfdv=0;
        projv = DotProduct(vel,glideDir);
        if (fabs(climbDir[X]*climbDir[Y]*climbDir[Z]) < 1e-10)
        {
           // {110} plane
           ScrewMDFit(burgMag, B110, K110, a110, projv, &projf, &dprojfdv);
        }
         else
        {
           // {112} plane
           real8 tdotb;
           tdotb = DotProduct(ksi,burg);

           if ( (tdotb >= 0 && projv >= 0.0) || (tdotb <0 && projv <0) )
           {
              // twinning
              ScrewMDFit(burgMag, B112T, K112T, a112T, projv, &projf, &dprojfdv);
           }
           else
           {
              // anti-twinning
              ScrewMDFit(burgMag, B112AT, K112AT, a112AT, projv, &projf, &dprojfdv);
           }
        }

        Vec3TransposeAndMult(glideDir, glideDirMat);
        Vec3TransposeAndMult(climbDir, climbDirMat);
        Vec3TransposeAndMult(lineDir,  lineDirMat);

        for (m = 0; m < 3; m++)
           for (n = 0; n < 3; n++)
              dfScrewDragdv[m][n] = dragClimb * climbDirMat[m][n] +
                                    dragLine  * lineDirMat [m][n];

        Matrix33Vector3Multiply(dfScrewDragdv, vel, fScrewDrag);

        for (m = 0; m < 3; m++)
           for (n = 0; n < 3; n++)
              dfScrewDragdv[m][n] += dprojfdv * glideDirMat[m][n];


        for (m = 0; m < 3; m++)
           fScrewDrag[m] += projf * glideDir[m];

        return;
}


/**************************************************************************
 *
 *      Function:     LinEdgeDrag
 *      Description:  This function returns the drag function g_edge
 *                    along the line, the glide directions and the climb
 *                    directions as well as its derivative Dg_edge/Dv.
 *                    EdgeDrag is linear.
 *
 *************************************************************************/
static void LinEdgeDrag(real8 vel[3], real8 burg[3], real8 lineDir[3],
                     real8 dragLine, real8 dragClimb, real8 Be110g, real8 Be112g,
                     real8 fEdgeDrag[3], real8 dfEdgeDragdv[3][3])
{
        int   m, n;
        real8 dragGlide;
        real8 climbDir[3], glideDir[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];

	// For edges, glide direction is the Burgers vector normalized.
        VECTOR_COPY(glideDir, burg);
        NormalizeVec(glideDir);

        NormalizedCrossVector(lineDir, glideDir,climbDir);

        if ( fabs(climbDir[X]*climbDir[Y]*climbDir[Z]) < 1e-10)
        {
           // {110} plane
           dragGlide = Be110g;
        }
        else
        {
           // {112} plane
           dragGlide = Be112g;
        }

        Vec3TransposeAndMult(glideDir, glideDirMat);
        Vec3TransposeAndMult(climbDir, climbDirMat);
        Vec3TransposeAndMult( lineDir,  lineDirMat);

        for (m = 0; m < 3; m++) {
            for (n = 0; n < 3; n++) {
                dfEdgeDragdv[m][n] = dragClimb * climbDirMat[m][n] +
                                     dragLine  *  lineDirMat[m][n] +
                                     dragGlide * glideDirMat[m][n];
            }
        }

        Matrix33Vector3Multiply(dfEdgeDragdv, vel, fEdgeDrag);

        return;
}


#if NON_LINEAR_EDGE
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
 *      Function:     EdgeDrag
 *      Description:  This function returns the drag function g_edge
 *                    along the line, the glide directions and the climb
 *                    directions as well as its derivative Dg_edge/Dv.
 *                    EdgeDrag is linear.
 *
 *************************************************************************/
static void EdgeDrag(real8 burgMag, real8 vel[3], real8 burg[3], real8 lineDir[3],
                     real8 dragLine, real8 dragClimb,
                     real8 fEdgeDrag[3], real8 dfEdgeDragdv[3][3])
{
        int   m, n;
        real8 climbDir[3], glideDir[3];
        real8 glideDirMat[3][3], climbDirMat[3][3], lineDirMat[3][3];

        // For edges, glide direction is the Burgers vector normalized.
        VECTOR_COPY(glideDir, burg);
        NormalizeVec(glideDir);

        NormalizedCrossVector(lineDir, glideDir,climbDir);

        real8 projv, projf=0, dprojfdv=0;
        projv = DotProduct(vel, glideDir);

        if ( fabs(climbDir[X]*climbDir[Y]*climbDir[Z]) < 1e-10)
        {
            // {110} plane
            EdgeMDFit(burgMag, projv, &projf, &dprojfdv);
        }
        else
        {
            // {112} plane
            EdgeMDFit(burgMag, projv, &projf, &dprojfdv);
        }

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
#endif


/**************************************************************************
 *
 *      Function:     DecomposeArm
 *      Description:  This function decomposes a dislocation segment into
 *                    screw and edges components.  There are two edge
 *                    components which are the highest projections of the
 *                    dislocation segment onto the different edge directions.
 *                    This function assumes the first edge component index
 *                    is known. The second one is computed within this function.
 *
 *                    Decomposition of a segment is given by:
 *
 *                        seg = fabs(screwLength * burgNorm) +
 *                              fabs(edge1Length * refPlanes[edge1Index]) +
 *                              fabs(edge2Length * refPlanes[edge2Index])
 *      Arguments:
 *
 *          IN:   seg        Dislocation segment vector specified in the
 *                           crystal frame.
 *
 *          IN:   numEdgePlanes Number of reference planes associated
 *                           with the burgers vector.
 *
 *          IN:   burgNorm   Normalized burgers vector specified in the
 *                           crystal frame.
 *
 *          IN:   refEdges   Edge reference planes for this burgers vector
 *                           specified in the crystal frame.
 *
 *          OUT:  screwLen   Length of the segment along the screw direction
 *
 *          OUT:  edge1Len   Length of the segment along the first edge
 *
 *          OUT:  edge2Len   Length of the segment along the second edge
 *
 *          IN:   edge1Index Index in <refPlanes> of the edge direction
 *                           with the highest projection
 *
 *          OUT:  edge1Index Index in <refPlanes> of the edge direction
 *                           with the 2nd highest projection
 *
 *************************************************************************/
static void DecomposeArm(real8 seg[3], int numEdgePlanes,
                         real8 burgNorm[3],
                         real8 refEdges[6][3],
                         real8 *screwLength, real8 *edge1Length,
                         real8 *edge2Length, int edge1Index, int *edge2Index)
{
   real8 ltot, ls, segLeft[3], tempLength[6];

/*
 *      Calculate screw length
 */
        ls = DotProduct(seg, burgNorm);

        ltot = Normal(seg);
        ls = MIN(ls,ltot);

/*
 *      Now get the lengths for both edges
 */
        segLeft[0] = seg[0] - ls * burgNorm[0];
        segLeft[1] = seg[1] - ls * burgNorm[1];
        segLeft[2] = seg[2] - ls * burgNorm[2];

        for (int i=0; (i<numEdgePlanes); i++)
            tempLength[i] = fabs(DotProduct(refEdges[i], segLeft));

/*
 *      We know the index of the highest edge value, now find the index
 *      of the second highest.
 */
        *edge2Index = Get2ndMaxIndex(tempLength, numEdgePlanes, edge1Index);

        DecompVec(segLeft, refEdges[edge1Index], refEdges[*edge2Index], tempLength);

        *edge1Length = MIN(ltot,fabs(tempLength[0]));
        *edge2Length = MIN(ltot,fabs(tempLength[1]));
        *screwLength = fabs(ls);

#if 0
        // Check seg decomposition

        Print3("seg",seg);
        printf("seg = [%.15e %.15e %.15e];\n\n",
               ls*burgNorm[X] + *edge1Length * refEdges[edge1Index][X] + *edge2Length * refEdges[*edge2Index][X],
               ls*burgNorm[Y] + *edge1Length * refEdges[edge1Index][Y] + *edge2Length * refEdges[*edge2Index][Y],
               ls*burgNorm[Z] + *edge1Length * refEdges[edge1Index][Z] + *edge2Length * refEdges[*edge2Index][Z]);
#endif

        return;
}

/**************************************************************************
 *
 *      Function:     ForceErrorCalc
 *
 *************************************************************************/
static void ForceErrorCalc(Param_t *param, real8 vTest[3], int numNbrs,
                           int eEx[MAX_NBRS], int sEx[MAX_NBRS], int jctEx[MAX_NBRS],
                           real8 Ltot[MAX_NBRS],real8 Lscrew[MAX_NBRS],
                           real8 Ledge1[MAX_NBRS],real8 Ledge2[MAX_NBRS],
                           real8 lineDir[MAX_NBRS][3], real8 burg[MAX_NBRS][3],
                           real8 edgeDir1[MAX_NBRS][3], real8 edgeDir2[MAX_NBRS][3],
                           real8 glideDirRef[MAX_NBRS][3], real8 ndir[MAX_NBRS][3],
                           real8 fError[3], real8 dfErrordv[3][3],
                           int iLin)
{
   int m,n,i;
   real8 eps = 1.0e-12;
   real8 burgMag = param->burgMag;

   // Drag parameters for linear screw case
   real8 drags[8];
   drags[0] = param->TaDragEdge110;
   drags[1] = param->TaDragEdge112;
   drags[2] = param->TaDragScrew112T;
   drags[3] = param->TaDragScrew112AT;
   drags[4] = param->TaDragScrew110;
   drags[5] = param->TaDragLine;
   drags[6] = param->TaDragClimbEdge;
   drags[7] = param->TaDragClimbScrew;

   real8 dragJunction = param->TaDragLine;

#if PSEUDO_INVERSE
    drags[5]     = 0.0;
    dragJunction = 0.0;
#endif

   // Drag parameters for non-linear screw case
   real8 K110   = param->TaDragParam110K;
   real8 B110   = param->TaDragParam110B;
   real8 a110   = param->TaDragParam110alpha;

   real8 K112T  = param->TaDragParam112TK;
   real8 B112T  = param->TaDragParam112TB;
   real8 a112T  = param->TaDragParam112Talpha;

   real8 K112AT = param->TaDragParam112ATK;
   real8 B112AT = param->TaDragParam112ATB;
   real8 a112AT = param->TaDragParam112ATalpha;

   for (i=0; i<numNbrs; i++)
   {
      real8 fJunction[3], dfJunctiondv[3][3];
      real8 fScrew[3], dfScrewdv[3][3];
      real8 fEdge[3], dfEdgedv[3][3];

      if (fabs(Ltot[i]) < eps) continue;

/*
 *       Junction contribution
 */
      if (jctEx[i])
      {
         TAJunctionDrag(vTest, lineDir[i], dragJunction, drags[7], fJunction, dfJunctiondv);

         for (m = 0; m < 3; m++)
         {
            fError[m] -= 0.5 * Ltot[i] * fJunction[m];
            for (n = 0; n < 3; n++)
               dfErrordv[m][n] -= 0.5 * Ltot[i] * dfJunctiondv[m][n];
         }
      }
      else
      {  /* not a junction */
/*
 *       Screw contribution
 */
         if (sEx[i])
         {
            if (iLin == 1)
            {
               LinScrewDrag(vTest, lineDir[i], burg[i], ndir[i], glideDirRef[i],
                            drags[5], drags[7], drags[2], drags[3], drags[4],
                            fScrew, dfScrewdv);
            }
            else
            {
               ScrewDrag(burgMag, vTest, lineDir[i], burg[i], ndir[i], glideDirRef[i],
                         drags[5], drags[7], K110, B110, a110, K112T, B112T, a112T,
                         K112AT, B112AT, a112AT, fScrew, dfScrewdv);
            }


            for (m = 0; m < 3; m++)
            {
               fError[m] -= 0.5 * Lscrew[i] * fScrew[m];
               for (n = 0; n < 3; n++)
                  dfErrordv[m][n] -= 0.5 * Lscrew[i] * dfScrewdv[m][n];
            }
         }  /* end if (sEx[i]) */

/*
 *       Edge contribution
 */

         if (eEx[i])
         {
/*
 *          First edge
 */
#if NON_LINEAR_EDGE
            if (iLin == 1)
#endif
            {
                LinEdgeDrag(vTest, burg[i], edgeDir1[i],
                            drags[5], drags[6], drags[0], drags[1],
                            fEdge, dfEdgedv);
            }
#if NON_LINEAR_EDGE
            else
            {
                EdgeDrag(burgMag, vTest, burg[i], edgeDir1[i],
                         drags[5], drags[6], fEdge, dfEdgedv);
            }
#endif

            for (m=0; m<3; m++)
            {
               fError[m] -= 0.5 * Ledge1[i] * fEdge[m];
               for (n=0; n<3; n++)
                  dfErrordv[m][n] -= 0.5 * Ledge1[i] * dfEdgedv[m][n];
            }
/*
 *          Second edge
 */
#if NON_LINEAR_EDGE
            if (iLin == 1)
#endif
            {
                LinEdgeDrag(vTest, burg[i], edgeDir2[i],
                         drags[5], drags[6], drags[0], drags[1],
                         fEdge, dfEdgedv);
            }
#if NON_LINEAR_EDGE
            else
            {
                EdgeDrag(burgMag, vTest, burg[i], edgeDir2[i],
                         drags[5], drags[6], fEdge, dfEdgedv);
            }
#endif

            for (m=0; m<3; m++)
            {
               fError[m] -= 0.5 * Ledge2[i] * fEdge[m];
               for (n=0; n<3; n++)
                  dfErrordv[m][n] -= 0.5 * Ledge2[i] * dfEdgedv[m][n];
            }
         } /* end if (eEx[i]) */

      }

   } /* end loop over segments */

}


/**************************************************************************
 *
 *      Function:     MobilityLaw_Ta
 *      Description:  This function calculates the velocity for a single
 *                    specified node.
 *
 *      Arguments:
 *          IN:  node                Pointer to the node for which to
 *                                   calculate velocity
 *          OUT: invdfErrordv        Array in which to return the inverse of
 *                                   the drag matrix used in the function.
 *          OUT: glideConstraints    Array in which to return the normals
 *                                   to the glide planes to which this node
 *                                   is restricted.
 *          OUT: numGlideConstraints Number of glide plane normals returned
 *
 *      Returns: 0 if successful
 *               1 if function could not converge on a velocity
 *
 *************************************************************************/
int Mobility_Ta(Home_t *home, Node_t *node, MobArgs_t *mobArgs)
{
        int   i, j, m, n;
        int   count;
        int   TotalCount;
        int   numNbrs;
        int   numNonZeroLenSegs = 0;
        int   numPlanes=6, iLin=1;
        int  *burgFirstPlaneIndex;
        int   edgeExists[MAX_NBRS];
        int   screwExists[MAX_NBRS];
        int   junctionExists[MAX_NBRS];
        real8 mu, maxSegLen,  adjustment=1.0;
        real8 massMult=0.0, massMatrix[3][3], eps;
        real8 vTest[3], vOld[3], vLin[3], fn[3], rt[3];
        real8 (*planeList)[3];
        real8 segTotLength[MAX_NBRS],segScrewLength[MAX_NBRS];
        real8 segEdge1Length[MAX_NBRS], segEdge2Length[MAX_NBRS];

        real8 burg[MAX_NBRS][3], ndir[MAX_NBRS][3], lineDir[MAX_NBRS][3];
        real8 glideDirRef[MAX_NBRS][3];
        real8 edgeDir1[MAX_NBRS][3], edgeDir2[MAX_NBRS][3];
        real8 edgeLength;

        real8 fError[3], dfErrordv[3][3];
        real8 fErrorTol;
        real8 forceError, totLen,forceErrorOld;
        real8 correction[3], inertia[3];
        Param_t *param;
        Node_t  *nbrNode;

        int forceisgood = 0;

/*
 *      If any of the arms of the node has a burgers vector that
 *      has explicitly been set to be sessile (via control file
 *      inputs), the node may not be moved.
 */
        if (NodeHasSessileBurg(home, node))
        {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      If the node is pinned in ALL dimensions it cannot be moved
 *      so just zero the velocity and return
 */
        if (HAS_ALL_OF_CONSTRAINTS(node->constraint, PINNED_NODE))
        {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

        param = home->param;

        mu = param->shearModulus;
        maxSegLen = param->maxSeg;

        real8 (*invdfErrordv)[3]     =  mobArgs->invDragMatrix;
     // real8 (*glideConstraints)[3] =  mobArgs->glideConstraints;
        int    *numGlideConstraints  = &mobArgs->numGlideConstraints;

        Matrix33_Zero(invdfErrordv);
        *numGlideConstraints = 0;

        eps = 1.0e-12;

/*
 *      If we need to include inertial terms, need some more initializations
 */
        if (param->includeInertia != 0)
            massMult = 0.25 * param->massDensity * (param->burgMag * param->burgMag);

        numNbrs = node->numNbrs;

/*
 *      This function saves various pieces of info for each neighbor
 *      node using static arrays for speed rather than doing memory
 *      allocations.  If the arrays are too  small, abort with an error
 *      to let the caller increase the array sizes.
 */
        if (numNbrs > MAX_NBRS)
        {
            Fatal("Mobility_Ta: Node segment count (%d) exceeds\n"
                  "maximum allowed (%d) in static arrays.  Increase MAX_NBRS\n"
                  "and recompile\n", node->numNbrs, MAX_NBRS);
        }

/*
 *      Set some pointers to the reference burgers vector and glide plane
 *      lists that were created during initialization.
 */
        //burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;

/*
 *      Loop over all arms attached to the node to get some information
 *      for the Newton-Raphson procedure
 */
        totLen = 0.0;

        for (i=0; (i<numNbrs); i++)
        {
            real8 rtSquared;

            nbrNode = GetNeighborNode(home, node, i);

            if (nbrNode == (Node_t *)NULL)
            {
                printf("WARNING: Neighbor not found at %s line %d\n", __FILE__, __LINE__);
                continue;
            }

/*
 *          Get the segment line direction
 */
            rt[0] = nbrNode->x - node->x;
            rt[1] = nbrNode->y - node->y;
            rt[2] = nbrNode->z - node->z;
            ZImage(param, &rt[0], &rt[1], &rt[2]);

/*
 *          If the segment is zero length (which may happen when the
 *          mobility function is called from SplitMultiNodes()), just
 *          skip the segment.
 */
            if ((rtSquared = DotProduct(rt, rt)) < eps)
            {
                segTotLength[i] = 0.0;
                continue;
            }

/*
 *          Calculate the segment total length
 */
            segTotLength[i] = sqrt(rtSquared);
            totLen += segTotLength[i];
            numNonZeroLenSegs++;

/*
 *          Save the burgers vector and glide plane of the segment
 */
            burg[i][X] = node->burgX[i];
            burg[i][Y] = node->burgY[i];
            burg[i][Z] = node->burgZ[i];

            ndir[i][X] = node->nx[i];
            ndir[i][Y] = node->ny[i];
            ndir[i][Z] = node->nz[i];
            NormalizeVec(ndir[i]);

/*
 *          If necessary, rotate the burgers vector and line sense from the
 *          laboratory frame to the crystal frame
 */
            if (param->useLabFrame)
            {
                real8 burgRot[3], rtRot[3], planeRot[3];

                Matrix33Vector3Multiply(param->rotMatrixInverse, burg[i], burgRot );
                Matrix33Vector3Multiply(param->rotMatrixInverse, rt     , rtRot   );
                Matrix33Vector3Multiply(param->rotMatrixInverse, ndir[i], planeRot);

                VECTOR_COPY(burg[i], burgRot);
                VECTOR_COPY(rt, rtRot);
                VECTOR_COPY(ndir[i], planeRot);
            }
/*
 *          Get the index in the reference Burgers vectors list
 *          of current Burgers vector, glide direction and plane
 *          direction.
 */
            int bIndex=-1, pIndex = -1, gIndex = -1;
            GetBCCAllIndices(home,burg[i],ndir[i],&bIndex, &pIndex, &gIndex);

/*
 *          A dislocation segment is a junction if both its
 *          Burgers vector and its planes are pre-defined as
 *          junction Burgers vector and sessile plane.
 *          A bIndex > 0 can have a sessile plane too.
 */

            if (bIndex >= BCC_NUM_GLIDE_BURG)
            {
/*
 *             Junction Burgers vectors
 */
               junctionExists[i] = 1;
               edgeExists    [i] = 0;
               screwExists   [i] = 0;
            }
            else
            {
/*
 *             Not a junction
 */
               junctionExists[i] = 0;
               real8 edgeRef[6][3];

/*
 *             define the reference glideDir
 */
               VECTOR_COPY(glideDirRef[i],planeList[gIndex]);

/*
 *             Get the corresponding planes and edges directions.
 *             Can look among glide planes only. There should not
 *             be more than 6 glide planes.
 */
               int edgeIndex1=0;
               int edgeIndex2=0;

/*
 *             The second edge direction is found
 *             in the group of glissile planes only.
 */
               int start = burgFirstPlaneIndex[bIndex];
               int end   = burgFirstPlaneIndex[bIndex] + numPlanes;

               for (m=0, j=start; j<end; j++, m++)
                  VECTOR_COPY(edgeRef[m], planeList[j]);

/*
 *             The first edge direction is the glide direction
 *             corresponding to gIndex.
 */
               edgeIndex1 = gIndex -start;

/*
 *             Decompose segment into screw and edges components.
 *             Compute length of screw and edges as well as their
 *             directions.
 */
               DecomposeArm(rt, numPlanes, burg[i],  edgeRef,
                            &(segScrewLength[i]), &(segEdge1Length[i]),
                            &(segEdge2Length[i]), edgeIndex1, &edgeIndex2);

               VECTOR_COPY(edgeDir1[i], edgeRef[edgeIndex1]);
               VECTOR_COPY(edgeDir2[i], edgeRef[edgeIndex2]);

               edgeLength = sqrt((segEdge1Length[i] * segEdge1Length[i]) +
                                 (segEdge2Length[i] * segEdge2Length[i]));

               edgeExists[i]  = ((edgeLength / segTotLength[i]) > 1.0e-04);
               screwExists[i] = ((segScrewLength[i]/segTotLength[i])>1.0e-04);
           } /* end glide Burgers vector */

/*
 *          Calculate tangents
 */
            lineDir[i][X] = rt[X] / segTotLength[i];
            lineDir[i][Y] = rt[Y] / segTotLength[i];
            lineDir[i][Z] = rt[Z] / segTotLength[i];
        }  /* loop i over neighbors */

/*
 *      It's possible this function was called for a node which only
 *      had zero length segments (during SplitSurfaceNodes() for example).
 *      If that is the case, just set the velocity to zero and return;
 */
        if (numNonZeroLenSegs == 0)
        {
            node->vX = 0.0;
            node->vY = 0.0;
            node->vZ = 0.0;
            return(0);
        }

/*
 *      Set up the mass matrix: if not needed, zero it out.
 */
        for (m=0; (m<3); m++)
        {
            for (n=0; (n<3); n++)
            {
                massMatrix[m][n] = 0.0;

                if ((m == n) && (param->includeInertia))
                    massMatrix[m][n] = totLen * massMult / param->deltaTT;
            }
        }

#if !PSEUDO_INVERSE
/*
 *     Linear solution in case of non-convergence of
 *     non-linear solution
 */
        count = 0;
        TotalCount   = 5000;
        iLin = 1;

        // ParaDiS force in Pa.
        fn[0] = node->fX;
        fn[1] = node->fY;
        fn[2] = node->fZ;

        // Initial velocity guess
        vTest[0] = node->vX;
        vTest[1] = node->vY;
        vTest[2] = node->vZ;

        VECTOR_COPY(vOld, vTest);
        VECTOR_ZERO(vLin);

        forceisgood = 0;

        fErrorTol = MAX((DotProduct(fn, fn) * 1.0e-16),
                        (mu * mu * maxSegLen * maxSegLen * 1.0e-24));

/*
 *      If necessary, rotate the force and velocity vectors from the
 *      laboratory frame to the crystal frame
 */
        if (param->useLabFrame)
        {
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
        if (param->includeInertia)
        {
            real8 vDiff[3];
            vDiff[0] = vTest[0] - vOld[0];
            vDiff[1] = vTest[1] - vOld[1];
            vDiff[2] = vTest[2] - vOld[2];
            Matrix33Vector3Multiply(massMatrix, vDiff, inertia);
        }
        else
        {
            VECTOR_ZERO(inertia);
        }

        forceError = fErrorTol + 1;
        forceErrorOld = forceError;

/*
 *      Begin the Newton-Raphson loop
 */
        while ((count < TotalCount) && (forceError > fErrorTol))
        {
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
           ForceErrorCalc(param, vTest, numNbrs, edgeExists, screwExists, junctionExists,
                          segTotLength, segScrewLength,  segEdge1Length, segEdge2Length,
                          lineDir, burg, edgeDir1, edgeDir2, glideDirRef, ndir,
                          fError, dfErrordv, iLin);


           if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 )
           {
               printf("det of derivative = %g\n", Matrix33_Det(dfErrordv) );
               Print3x3("dfErrordv", dfErrordv);
               PrintNodeDebug(home, node);
               Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
            }

            Matrix33Vector3Multiply(invdfErrordv, fError, correction);
            forceError = DotProduct(fError, fError);

            vTest[0] -= adjustment*correction[0];
            vTest[1] -= adjustment*correction[1];
            vTest[2] -= adjustment*correction[2];

            // This forceError is good. Continue Newton
            if (forceError <= forceErrorOld + 1e-6)
            {
               forceErrorOld = forceError;
               if (forceisgood > 100)
               {
                  adjustment = MIN(adjustment*2.0,1.0);
                  forceisgood = 0;
               }
               forceisgood += 1;
            }

            // If the current step has led to a worst error than before. Modify the Newton step.
            while (forceError > forceErrorOld && count > 2 && adjustment > 0.0156)
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

                 ForceErrorCalc(param, vTest, numNbrs, edgeExists, screwExists, junctionExists,
                                segTotLength, segScrewLength,  segEdge1Length, segEdge2Length,
                                lineDir, burg, edgeDir1, edgeDir2, glideDirRef, ndir,
                                fError, dfErrordv, iLin);


                 if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 )
                 {
                    printf("det of derivative = %g\n",Matrix33_Det(dfErrordv));
                    Print3x3("dfErrordv", dfErrordv);
                    PrintNodeDebug(home, node);
                    Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
                 }

                 Matrix33Vector3Multiply(invdfErrordv, fError, correction);
                 forceError = DotProduct(fError, fError);
              }

              count++;
        } /* end Newton-Raphson while */


        if ((count > TotalCount) || (forceError > fErrorTol))
        {
           printf("WARNING: MobilityLaw_Ta has not converged for node id %d in %d iterations\n",
                  node->myTag.index,count);
           return(1);
        }

/*
 *      A linear solution has been found. Keep it in case the
 *      non-linear solution does not converge.
 */
        VECTOR_COPY(vLin,vTest);
#endif

/*
 *      Setup for the Newton-Raphson loop
 */
        count = 0;
        TotalCount = 5000;
        iLin = 0;

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

/*
 *      If necessary, rotate the force and velocity vectors from the
 *      laboratory frame to the crystal frame
 */
        if (param->useLabFrame)
        {
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
        if (param->includeInertia)
        {
            real8 vDiff[3];
            vDiff[0] = vTest[0] - vOld[0];
            vDiff[1] = vTest[1] - vOld[1];
            vDiff[2] = vTest[2] - vOld[2];
            Matrix33Vector3Multiply(massMatrix, vDiff, inertia);
        }
        else
        {
            VECTOR_ZERO(inertia);
        }

        forceError = fErrorTol + 1;
        forceErrorOld = forceError;

#if PSEUDO_INVERSE
        real8 fmag = sqrt(DotProduct(fn,fn));
#endif

/*
 *      Begin the Newton-Raphson loop
 */
        int exitNR=0;
        while ((count < TotalCount) && (forceError > fErrorTol) && (!exitNR))
        {
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
           ForceErrorCalc(param, vTest, numNbrs, edgeExists, screwExists, junctionExists,
                          segTotLength, segScrewLength,  segEdge1Length, segEdge2Length,
                          lineDir, burg, edgeDir1, edgeDir2, glideDirRef, ndir,
                          fError, dfErrordv, iLin);

#if PSEUDO_INVERSE
            if (Matrix33_PseudoInverse(invdfErrordv, dfErrordv) != 0) {
                Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
            }
#else
           if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 ) 
           {
               printf("det of derivative = %g\n",Matrix33_Det(dfErrordv));
               Print3x3("dfErrordv", dfErrordv);
               PrintNodeDebug(home, node);
               Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
            }
#endif

            Matrix33Vector3Multiply(invdfErrordv, fError, correction);
            forceError = DotProduct(fError, fError);

            vTest[0] -= adjustment*correction[0];
            vTest[1] -= adjustment*correction[1];
            vTest[2] -= adjustment*correction[2];

#if PSEUDO_INVERSE
            if (fabs(forceError-forceErrorOld)/fmag < 1e-6) exitNR = 1;
#endif

            // This forceError is good. Continue Newton
            if (forceError <= forceErrorOld + 1e-6)
            {
               forceErrorOld = forceError;
               if (forceisgood > 100)
               {
                  adjustment = MIN(adjustment*2.0,1.0);
                  forceisgood = 0;
               }
               forceisgood += 1;
            }

            // If the current step has led to a worst error than before. Modify the Newton step.
            while (forceError > forceErrorOld && count > 2 && adjustment > 0.0156)
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

                 ForceErrorCalc(param, vTest, numNbrs, edgeExists, screwExists, junctionExists,
                                segTotLength, segScrewLength,  segEdge1Length, segEdge2Length,
                                lineDir, burg, edgeDir1, edgeDir2, glideDirRef, ndir,fError, dfErrordv, iLin);

#if PSEUDO_INVERSE
                if (Matrix33_PseudoInverse(invdfErrordv, dfErrordv) != 0) {
                    Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
                }
#else
                 if ( M33_INVERSE(invdfErrordv,dfErrordv) < 0 ) 
                 {
                    printf("det of derivative = %g\n",Matrix33_Det(dfErrordv));
                    Print3x3("dfErrordv", dfErrordv);
                    PrintNodeDebug(home, node);
                    Fatal("%s::%s(%d) : Cannot invert dfErrordv!", __FILE__, __func__, __LINE__ );
                 }
#endif

                 Matrix33Vector3Multiply(invdfErrordv, fError, correction);
                 forceError = DotProduct(fError, fError);
              }

              count++;
        } /* end Newton-Raphson while */


#if !PSEUDO_INVERSE
        if ((count > TotalCount) || (forceError > fErrorTol))
        {
           printf("WARNING: MobilityLaw_Ta has not converged for node id %d in %d iterations\n",
                  node->myTag.index,count);
           VECTOR_COPY(vTest,vLin);
        }
#endif

/*
 *      We were able to converge on a velocity
 *
 *      If needed, rotate the velocity vector back to the laboratory frame
 *      from the crystal frame.
 */
        if (param->useLabFrame)
        {
            real8 vTestRot[3];
            Matrix33Vector3Multiply(param->rotMatrix, vTest, vTestRot);
            VECTOR_COPY(vTest, vTestRot);
        }


        node->vX = vTest[0];
        node->vY = vTest[1];
        node->vZ = vTest[2];



#ifdef NAN_CHECK
        if ( (isnan(node->vX) || isinf(node->vX)) ||
             (isnan(node->vY) || isinf(node->vY)) ||
             (isnan(node->vZ) || isinf(node->vZ)))
        {
            PrintNode(node);
            Fatal("Mobility_Ta: (%d, %d) velocity error.", node->myTag.domainID, node->myTag.index);
        }
#endif

        return(0);
}
