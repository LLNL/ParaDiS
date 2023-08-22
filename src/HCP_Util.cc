/****************************************************************************
 *
 *      Module:         HCP_Util.c
 *
 *      Author:         S. Aubry March 2011
 *
 *      Description:    This module contains various functions required
 *                      to calculate the burgers vectors and associated
 *                      glide planes for an HCP material with the specified
 *                      c/a parameter.
 *
 *      Includes public functions:
 *
 *          FindJctPlanes()
 *          GetBurgIndexHCP()
 *          GetBurgList_HCP()
 *          FindGlidePlaneHCP()
 *
 *      Includes private functions:
 *
 *          AssignWeightToPlane()
 *          CheckPlanes()
 *          RotateSystem()
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

extern int  CheckPlanes (real8 b[3],real8 p[3]);

extern void RotateSystem(real8 XX[6], real8 YY[6], real8 ZZ[6],
                         real8 r1[3], real8 r2[3], real8 r3[3],
                         real8 Xp[6], real8 Yp[6], real8 Zp[6]);

/*-------------------------------------------------------------------------
 *
 *      Function:    CheckPlanes
 *      Description: Verify that the direction belongs to the plane
 *
 *      Arguments:
 *          IN:  b
 *          IN:  p
 *
 *      Returns:  1 if OK, zero if not OK
 *
 *------------------------------------------------------------------------*/
int CheckPlanes(real8 b[3],real8 p[3])
{
        real8 dot = DotProduct(b,p);

        if (fabs(dot) > 1e-5) {
            printf("b = %f %f %f\n", b[0], b[1], b[2]);
            printf("p = %f %f %f\n", p[0], p[1], p[2]);
            printf("b . p= %f\n", dot);
            return(1);
        }

        return(0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    MergePlaneLists
 *      Description: Given two lists of glide planes, merge them into
 *                   a single array in which there are no duplicate
 *                   planes.
 *
 *      Arguments:
 *          IN:  planeList1,         The two lists of planes to be merged
 *               planeList2          together
 *
 *          IN:  planeList1Len,      The number of planes on each of the
 *               planeList2Len       respective incoming plane lists.
 *
 *          IN:  newPlaneListMaxLen  Maximum number of planes that can be
 *                                   returned on <newPlaneList>
 *          OUT: newPlaneList        Array in whcih to return the merged
 *                                   plane list
 *          OUT: newPlaneListLen     Number of planes returned in
 *                                   <newPlaneList>
 *
 *------------------------------------------------------------------------*/
static void MergePlaneLists(real8 (*planeList1)[3],
                            real8 (*planeList2)[3],
                            int planeList1Len,
                            int planeList2Len,
                            int newPlaneListMaxLen,
                            real8 (*newPlaneList)[3],
                            int *newPlaneListLen)
{
        int   i, j;
        real8 eps = 1e-5;


/*
 *      Copy the first plane list into the new plane list
 */
        *newPlaneListLen = planeList1Len;

        if (*newPlaneListLen >= newPlaneListMaxLen) {
            Fatal("MergePlaneLists: newPlaneList is too small to "
                  "accomodate\n    the merged plane list!");
        }

        for (i = 0; i < planeList1Len; i++) {
            VECTOR_COPY(newPlaneList[i], planeList1[i]);
        }


/*
 *      Now find any planes in the second list that are not in
 *      the first list and append those to the new list
 */
        for (i = 0; i < planeList2Len; i++) {
            int isDuplicate;

            isDuplicate = 0;

            for (j = 0; j < planeList1Len; j++) {
                real8 p1Dotp2;
                real8 signVal;

                p1Dotp2 = DotProduct(planeList1[j], planeList2[i]);
                signVal = ((p1Dotp2 < 0.0) ?  -1.0 : 1.0);

/*
 *              If this plane from list 2 is also included in list 1, there's
 *              no need to check any more for the plane, so set a flag and
 *              break out of the loop.
 */
                if ((fabs(planeList1[j][0] - signVal*planeList2[i][0]) < eps) &&
                    (fabs(planeList1[j][1] - signVal*planeList2[i][1]) < eps) &&
                    (fabs(planeList1[j][2] - signVal*planeList2[i][2]) < eps)) {

                    isDuplicate = 1;

                    break;
                }
            }

/*
 *          If the plane from list 2 was not a duplicate of any from list
 *          1, add the plane to the new list.
 */
            if (!isDuplicate) {

                if (*newPlaneListLen >= newPlaneListMaxLen) {
                    Fatal("MergePlaneLists: newPlaneList is too small to "
                          "accomodate\n    the merged plane list!");
                }

                VECTOR_COPY(newPlaneList[*newPlaneListLen], planeList2[i]);
                *newPlaneListLen += 1;
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    FindJctPlanes
 *
 *      Description: Calculate the list of possible planes for
 *                   a junction with burgers vector <burg> where
 *                   <burg> is the the sum of two glide burgers
 *                   vectors (we'll refer to them as burgP and burgQ).
 *
 *      Arguments:
 *          IN:  sizeP  Number of planes in plane list <P>
 *          IN:  sizeQ  Number of planes in plane list <Q>
 *          IN:  P      List of possible planes for the first glide
 *                      burgers vector ( i.e. burgP).
 *          IN:  Q      List of possible planes for the second glide
 *                      burgers vector ( i.e. burgQ).
 *          IN:  burg   Burgers vector of the junction
 *          OUT: PQ     Array of plane normals returned to caller
 *          OUT: sizePQ Number of plane normals returned in <PQ>
 *
 *------------------------------------------------------------------------*/
static void FindJctPlanes(int sizeP, int sizeQ,
                          real8 (*P)[3], real8 (*Q)[3],
                          real8 burg[3], real8 (*PQ)[3], int *sizePQ)
{
        int   i, j, k, m;
        real8 paux[3], qaux[3];
        real8 normq;


        m = 0;

        for (i = 0; i < sizeP; i++) {

            for (j = 0; j < sizeQ; j++) {

                cross(P[i], Q[j], paux);
                cross(paux, burg, qaux);

                normq = Normal(qaux);

/*
 *              In the case where P == Q, the two dislocations are on the
 *              same plane and the cross product is ill-defined.
 *              In that particular case, qaux = P = Q.
 */ 
                if (normq <= 1e-5) { 
                    // Check that the planes are really the same 
                    real8 p1Dotp2;
                    real8 signVal;
                    real8 eps = 1e-5;

                    p1Dotp2 = DotProduct(P[i],Q[j]);
                    signVal = ((p1Dotp2 < 0.0) ?  -1.0 : 1.0);

                    if ((fabs(Q[j][0] - signVal*P[i][0]) > eps) ||
                        (fabs(Q[j][1] - signVal*P[i][1]) > eps) ||
                        (fabs(Q[j][2] - signVal*P[i][2]) > eps)) {

                        printf("P=%f %f %f\n", P[i][0], P[i][1], P[i][2]);
                        printf("Q=%f %f %f\n", Q[j][0], Q[j][1], Q[j][2]);

                        Fatal("FindJctPlanes: Planes are different but "
                              "their cross product is zero?\n");
                    } 

                    qaux[0] = P[i][0];
                    qaux[1] = P[i][1];
                    qaux[2] = P[i][2];
                    normq = Normal(qaux);
                }

                NormalizeVec(qaux);

/*
 *              Need to find all the planes but without repetitions, so
 *              check if the plane is already on the <PQ> array before
 *              adding it in.
 */
                if (normq > 1e-5) {

                    for  (k = 0; k < m; k++) {
/*
 *                      PQ[k][*] and qaux are already normalized and their
 *                      normals are greater than zero
 */
                        if (fabs(DotProduct(PQ[k], qaux)) > 0.99995) {
                            break;
                        }
                    }

/*
 *                  If the plane was not found on the <PQ> array, add it.
 *                  Note: If a new vector has a near-zero component, the
 *                        component will be zeroed and the vector
 *                        renormalized
 */
                    if (k >= m) {
                        int n;

                        VECTOR_COPY(PQ[m], qaux);

                        for (n = 0; n < 3; n++) {
                            if (fabs(PQ[m][n]) < 1e-10) {
                                PQ[m][n] = 0.0;
                                NormalizeVec(PQ[m]);
                            }
                        }

                        m = m + 1;
                    }

                } // end if (normq > 1e-5)

            } // end for (j = 0; j < sizeQ; j++)

        } // end for (i = 0; i < sizeP; i++)

        *sizePQ = m;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetBurgIndexHCP
 *
 *      Description: Find the index in the reference burgers vector list
 *                   of the burgers vector matching the one provided by
 *                   the caller.  Criteria for matching burgers vectors
 *                   are 1) burgers vector normals are the same and vectors
 *                   are equivalent up to the sign.
 *
 *      Parameters:
 *          IN:  burgVec    Burgers vector to look up.  This burgers vector
 *                          is NOT normalized and MUST be specified in the
 *                          crystalographic frame.
 *          IN:  useClosest Flag indicating if the function should return
 *                          the index of the closest matching burgers vector
 *                          if an exact match is not found. (0 == NO, 1 == YES)
 *          IN:  numBurg    Number of burgers vectors in the reference list
 *          IN:  burgList   Reference burgers vectors.  These burgers vectors
 *                          are NOT normalized and are specified in the
 *                          crystalographic frame.
 *          OUT: burgIndex  Location in which to return the index of the
 *                          reference burgers vector matching the one specified
 *                          by the caller.  If <useClosest> == 0 and no exact
 *                          match is found, this value will be set to -1.
 *
 *-----------------------------------------------------------------------*/
void GetBurgIndexHCP(real8 burgVec[3], int useClosest, int numBurg,
                     real8 (*burgList)[3], int *burgIndex)
{
        int   i, indexOfMin;
        real8 minVal;

/*
 *      The comparison is done this way because there are burgers vectors
 *      in the reference list whose normalized values are the same with
 *      a difference only in magnitude.  For example, for HCP, it is
 *      possible to have a [1 0 0] glide burgers vector and a [2 0 0]
 *      junction burgers vector.  The only difference is the magnitude
 *      but they need to be handled very differently.
 */
        indexOfMin = -1;
        minVal = 1e+10;
        
        for (i = 0; i < numBurg; i++) {
            real8 sb, tmp;
            
            real8 tmpVec[3];

            sb = (DotProduct(burgList[i], burgVec) < 0.0 ? -1.0 : 1.0);

            if ((fabs(Normal(burgList[i]) - Normal(burgVec)) < 1e-4) &&
                ((fabs(burgList[i][X] - sb * burgVec[X]) < 1.e-4) &&
                 (fabs(burgList[i][Y] - sb * burgVec[Y]) < 1.e-4) &&
                 (fabs(burgList[i][Z] - sb * burgVec[Z]) < 1.e-4))) {
                *burgIndex = i;
                return;
            }

/*
 *          Just a safety catch.  In case we don't find an exact match
 *          (perhaps due to precision differences between teh internally
 *          calcukated reference burges vectors and the user-supplied 
 *          burgers vectors in the restart file), keep track of the
 *          closest match.
 */
            if (useClosest) {
                cross(burgList[i], burgVec, tmpVec);
                tmp = DotProduct(tmpVec, tmpVec);

                if (tmp < minVal) {
                    minVal = tmp;
                    indexOfMin = i;
                }
            }
        }

/*
 *      We didn't find an exact match on the burgers vector, so return
 *      the index of the closest match.
 */
        *burgIndex = indexOfMin;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetBurgList_HCP()
 *
 *      Description: Calculate all the unique burgers vectors for an HCP
 *                   material with the specified cOa value, plus the
 *                   glide plane normals for each of those burgers vectors.
 *
 *                   WARNING!  The arrays returned from this routine are
 *                             dynamically calculated, but the number of
 *                             elements of many of these arrays have been
 *                             pre-computed.  As such, there are no checks
 *                             to verify that the array bounds are not
 *                             being exceeded.  Should you alter the
 *                             content of these arrays resulting in arrays
 *                             of different lengths, you must take care
 *                             to resize the arrays appropriately!
 *
 *      Arguments:
 *
 *          IN:  cOa
 *
 *          OUT: burgList  Array of <numBurg> burgers vectors.  Caller is 
 *                         responsible for deallocating this memory.
 *                         The burgers vectors in this list are NOT normalized
 *                         and are specified in the crystalographic frame
 *
 *          OUT: planeList Array of <numPlanes> plane normals.  Caller is
 *                         responsible for deallocating this memory.  The
 *                         plane normals in this list are specified in
 *                         the crystalographic frame.
 *
 *          OUT: splinterableBurgList Array of <numSplinterableBurgs>
 *                         structs specifying the burgers vectors that
 *                         may be splintered and the corresponding burgers
 *                         vectors into which they may be splintered.
 *                         Caller is responsible for deallocating this memory.
 *                         The burgers vectors in this list are NOT normalized
 *                         and are specified in the crystalographic frame
 *
 *          OUT: numPlanesPerBurg  Array of <numBurg> integers indicating the
 *                         number of planes associated with the corresponding
 *                         burgers vector.  Caller is responsible for
 *                         deallocating this memory.
 *
 *          OUT: burgFirstPlaneIndex  Array of <numBurg> integers indicating
 *                         the index in <planeList> of the first plane
 *                         associated with the corresponding burgers vector.
 *                         Caller is responsible for deallocating this memory.
 *
 *          OUT: numGlissilePlanes  Array of <numBurg> integers indicating
 *                         the number of glissile planes associated with
 *                         the corresponding burgers vectors.
 *
 *          OUT: planeType Array of <numPlanes> integers indicating the type
 *                         associated with each plane.  See comments in
 *                         GetPlaneType() function for list of valid values.
 *
 *          OUT: numBurg   Number of burgers vectors returned in <burgList>
 *
 *          OUT: numPlanes Number of planes returned in <planeList>
 *
 *          OUT: numSplinterableBurgs Number of brugers vectors returned
 *                         in <splinterableBurgList.>
 *
 *------------------------------------------------------------------------*/
void GetBurgList_HCP(real8 cOa, real8 (**burgList)[3], real8 (**planeList)[3],
                     SplinterableBurg_t **splinterableBurgList,
                     int **numPlanesPerBurg, int **burgFirstPlaneIndex,
                     int **numGlissilePlanes, int **planeType,
                     int *numBurg, int *numPlanes, int *numSplinterableBurgs)   
{

        int   i, j, k;
        int   planeOffset, bufSize;
        int   burgCnt, splinterableBurgCnt;
        int   splinterGlideBurg[4], splinterJunctBurg[4];
        int   *typeOfPlane, *glissilePlanesPerBurg;
        int   *planesPerBurg, *burgPlaneIndex;
        real8 plane[3];
        real8 tmpDir[3], tmpBurg[3];
        real8 (*bList)[3];
        real8 (*pList)[3];
        SplinterableBurg_t *splinterableBList;

        /* Points */
        real8 A[3] = {0.0, 0.0, 0.0};
        real8 B[3] = {-0.5, sqrt(3)*0.5, 0};
        real8 C[3] = { 0.5, sqrt(3)*0.5, 0};
  
        real8  delta[3] = {0, sqrt(3)/3., 0.5*cOa};
  
        real8 Ap[3] = {0, 0, cOa};
        real8 Bp[3] = {-0.5, sqrt(3)*0.5, cOa};
        real8 Cp[3] = { 0.5, sqrt(3)*0.5, cOa};
/*
 *      The number of HCP burgers vectors and the number of glide planes
 *      associated with each has been precalculated.  There are a total
 *      of 41 burgers vectors, (the first 10 are glide burgers vectors,
 *      the rest are junction burgers vectors).  There are a total
 *      of 422 planes, and the burgPlaneIndex array defines the index
 *      in the plane list of the first plane associated with each of
 *      the 41 burgers vectors.  For any given burgers vector, the
 *      glissile planes will be listed first, followed by the sessile
 *      planes.  The <glissilePlanesPerBurg> array will be populated
 *      with the number of glissile planes associated with each burgers
 *      vector.
 */
        static int burgStartPlaneIndex[HCP_NUM_TOT_BURG] = {
                0,   12,  24,  36,  48,  60,  72,  84,  96,  108,
                126, 140, 154, 162, 172, 182, 190, 200, 210, 224,
                232, 242, 250, 260, 268, 276, 290, 304, 310, 322,
                328, 334, 348, 362, 368, 374, 380, 394, 400, 406, 412
        };

/*
 *      Allocate some of the dynamic arrays... sizes are known
 *      ahead of time.
 */
        bufSize = HCP_NUM_TOT_BURG * sizeof(double [3]);
        bList = (double (*)[3])calloc(1, bufSize);

        bufSize = HCP_NUM_TOT_PLANES * sizeof(double [3]);
        pList = (double (*)[3])calloc(1, bufSize);

        bufSize =  4 * sizeof(SplinterableBurg_t);
        splinterableBList = (SplinterableBurg_t *)calloc(1, bufSize);

        bufSize = HCP_NUM_TOT_BURG * sizeof(double [3]);
        planesPerBurg = (int *)calloc(1, bufSize * sizeof(int));
        burgPlaneIndex = (int *)calloc(1, bufSize * sizeof(int));

        bufSize = HCP_NUM_TOT_PLANES * sizeof(double [3]);
        glissilePlanesPerBurg = (int *)calloc(1, bufSize * sizeof(int));
        typeOfPlane = (int *)calloc(1, bufSize * sizeof(int));

/* 
 *      Initialize the type of plane to be sessile. 
 *      Only glissile planes will be assigned later on.
 */
	for (k = 0; k < HCP_NUM_TOT_PLANES; k++) {
	  typeOfPlane[k] = 4;
	} 


/*
 *      Since we precalculated the number of planes for every burgers
 *      vector, we can set the index (in the plane array) of the first
 *      plane for each burgers vector right now.
 */
        for (k = 0; k < HCP_NUM_TOT_BURG; k++) {
            burgPlaneIndex[k] = burgStartPlaneIndex[k];
        }

/* 
 *      First HCP_NUM_GLIDE_BURG Burgers vectors are the glide burgers
 *      vectors.  Initialize them here.
 *
 *      WARNING!  These glide burgers vectors MUST be the first ones in
 *                the burgers vector list.  Other portions of the code
 *                depend on this being so!
 */

/*
 *      Set the plane type (character) for every plane in the
 *      list.  Type will be set to one of the following:
 *
 *          0 if plane is basal
 *          1 if plane is prismatic
 *          2 if plane is first pyramidal
 *          3 if plane is second pyramidal
 *          4 if plane is sessile
 */

        planeOffset = burgPlaneIndex[0];

        bList[0][0] = B[0] - A[0];
        bList[0][1] = B[1] - A[1];
        bList[0][2] = B[2] - A[2];

        GetPlaneNormFromPoints(C,A,B,plane);       
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=0; /* basal */

        GetPlaneNormFromPoints(Cp,A,B,plane);      
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2; /* 1st pyramidal */

        GetPlaneNormFromPoints(A,B,Bp,plane);
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1; /* prismatic */

        GetPlaneNormFromPoints(C,Bp,Ap,plane);
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2; /* 1st pyramidal */

        glissilePlanesPerBurg[0] = 4;
        planesPerBurg[0] = 4;

/*
 *      BC
 */
        planeOffset = burgPlaneIndex[1];

        bList[1][0] = C[0] - B[0];
        bList[1][1] = C[1] - B[1];
        bList[1][2] = C[2] - B[2];

        GetPlaneNormFromPoints(C,A,B,plane);
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=0; /* basal */

        GetPlaneNormFromPoints(A,Bp,Cp,plane);
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2; /* 1st pyramidal */

        GetPlaneNormFromPoints(Cp,C,B,plane);
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1; /* prismatic */

        GetPlaneNormFromPoints(Ap,C,B,plane);   
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2;  /* 1st pyramidal */

        glissilePlanesPerBurg[1] = 4;
        planesPerBurg[1] = 4;

/*
 *      CA
 */
        planeOffset = burgPlaneIndex[2];

        bList[2][0] = A[0] - C[0];
        bList[2][1] = A[1] - C[1];
        bList[2][2] = A[2] - C[2];

        GetPlaneNormFromPoints(A,C,B,plane); 
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=0; /* basal */

        GetPlaneNormFromPoints(Bp,A,C,plane);
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2; /* 1st pyramidal */

        GetPlaneNormFromPoints(A,C,Ap,plane);
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1; /* prismatic */

        GetPlaneNormFromPoints(B,Cp,Ap,plane);
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2;  /* 1st pyramidal */

        glissilePlanesPerBurg[2] = 4;
        planesPerBurg[2] = 4;

/*
 *      ABp
 */
        planeOffset = burgPlaneIndex[3];

        bList[3][0] = Bp[0] - A[0];
        bList[3][1] = Bp[1] - A[1];
        bList[3][2] = Bp[2] - A[2];

        GetPlaneNormFromPoints(delta,A,Bp,plane);  
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=3; /* 2nd pyramidal */

        GetPlaneNormFromPoints(A,Bp,Cp,plane);     
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2; /* 1st pyramidal */

        GetPlaneNormFromPoints(A,B,Bp,plane);   
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1; /* prismatic */

        GetPlaneNormFromPoints(Bp,A,C,plane);   
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2; /* 1st pyramidal */

        glissilePlanesPerBurg[3] = 4;
        planesPerBurg[3] = 4;

/*
 *      BCp
 */
        planeOffset = burgPlaneIndex[4];

        bList[4][0] = Cp[0] - B[0];
        bList[4][1] = Cp[1] - B[1];
        bList[4][2] = Cp[2] - B[2];

        GetPlaneNormFromPoints(C,B,Cp,plane); 
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=1; /* prismatic */

        GetPlaneNormFromPoints(Cp,A,B,plane);    
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2; /* 1st pyramidal */

        GetPlaneNormFromPoints(delta,B,Cp,plane);  
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=3; /* 2nd pyramidal */

        GetPlaneNormFromPoints(B,Cp,Ap,plane);
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2; /* 1st pyramidal */

        glissilePlanesPerBurg[4] = 4;
        planesPerBurg[4] = 4;

/*
 *      CAp
 */
        planeOffset = burgPlaneIndex[5];
 
        bList[5][0] = Ap[0] - C[0];
        bList[5][1] = Ap[1] - C[1];
        bList[5][2] = Ap[2] - C[2];

        GetPlaneNormFromPoints(C,delta,Ap,plane);  
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=3;/* 2nd pyramidal */

        GetPlaneNormFromPoints(Ap,C,B,plane);      
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2;/* 1st pyramidal */

        GetPlaneNormFromPoints(C,A,Ap,plane);      
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1;/* prismatic */

        GetPlaneNormFromPoints(Bp,C,Ap,plane);     
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2;/* 1st pyramidal */

        glissilePlanesPerBurg[5] = 4;
        planesPerBurg[5] = 4;

/*
 *      ApB
 */
        planeOffset = burgPlaneIndex[6];

        bList[6][0] = B[0] - Ap[0];
        bList[6][1] = B[1] - Ap[1];
        bList[6][2] = B[2] - Ap[2];

        GetPlaneNormFromPoints(Ap,delta,B,plane);  
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=3;/* 2nd pyramidal */

        GetPlaneNormFromPoints(Ap,C,B,plane);      
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2;/* 1st pyramidal */

        GetPlaneNormFromPoints(Ap,B,Bp,plane);     
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1;/* prismatic */

        GetPlaneNormFromPoints(B,Cp,Ap,plane);     
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2;/* 1st pyramidal */

        glissilePlanesPerBurg[6] = 4;
        planesPerBurg[6] = 4;

/*
 *      BpC
 */
        planeOffset = burgPlaneIndex[7];

        bList[7][0] = C[0] - Bp[0];
        bList[7][1] = C[1] - Bp[1];
        bList[7][2] = C[2] - Bp[2];

        GetPlaneNormFromPoints(Bp,delta,C,plane);  
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=3;/* 2nd pyramidal */

        GetPlaneNormFromPoints(C,Bp,Ap,plane);     
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2;/* 1st pyramidal */

        GetPlaneNormFromPoints(Bp,Cp,C,plane);     
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1;/* prismatic */

        GetPlaneNormFromPoints(Bp,C,A,plane);      
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2;/* 1st pyramidal */

        glissilePlanesPerBurg[7] = 4;
        planesPerBurg[7] = 4;

/*
 *      CpA
 */
        planeOffset = burgPlaneIndex[8];

        bList[8][0] = A[0] - Cp[0];
        bList[8][1] = A[1] - Cp[1];
        bList[8][2] = A[2] - Cp[2];

        GetPlaneNormFromPoints(Cp,A,delta,plane);  
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=3;/* 2nd pyramidal */

        GetPlaneNormFromPoints(A,Bp,Cp,plane);     
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=2;/* 1st pyramidal */

        GetPlaneNormFromPoints(Cp,A,Ap,plane);     
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1;/* prismatic */

        GetPlaneNormFromPoints(Cp,B,A,plane);      
        VECTOR_COPY(pList[planeOffset+3], plane);
	typeOfPlane[planeOffset+3]=2;/* 1st pyramidal */

        glissilePlanesPerBurg[8] = 4;
        planesPerBurg[8] = 4;

/*
 *      AAp
 */
        planeOffset = burgPlaneIndex[9];

        bList[9][0] = Ap[0] - A[0];
        bList[9][1] = Ap[1] - A[1];
        bList[9][2] = Ap[2] - A[2];

        GetPlaneNormFromPoints(C,A,Ap,plane);      
        VECTOR_COPY(pList[planeOffset+0], plane);
	typeOfPlane[planeOffset+0]=1;/* prismatic */

        GetPlaneNormFromPoints(C,B,Cp,plane);   
        VECTOR_COPY(pList[planeOffset+1], plane);
	typeOfPlane[planeOffset+1]=1;/* prismatic */

        GetPlaneNormFromPoints(Ap,A,Bp,plane);  
        VECTOR_COPY(pList[planeOffset+2], plane);
	typeOfPlane[planeOffset+2]=1;/* prismatic */

        glissilePlanesPerBurg[9] = 3;
        planesPerBurg[9] = 3;



/*
 *      Set initial count for the known burgers vectors that have been
 *      added to the arrays above.
 */
        burgCnt = HCP_NUM_GLIDE_BURG;
        splinterableBurgCnt = 0;

/*
 *      Calculate all the possible junction burgers vectors that might
 *      be formed by interactions between the glide burgers vectors
 *      and find the associated planes.
 */
        for (i = 0; i < HCP_NUM_GLIDE_BURG; i++) {

            for (j = i + 1; j < HCP_NUM_GLIDE_BURG; j++) {
                int   bIndex;
                real8 dir[3];
/*
 *              Possible direction 1
 */ 
                dir[0] = bList[i][0] + bList[j][0];
                dir[1] = bList[i][1] + bList[j][1];
                dir[2] = bList[i][2] + bList[j][2];

/*
 *              Check if this Burgers vector is already in the list
 *              Criteria are: burgers vector normals are the same
 *              and vectors are equivalent up to a sign.
 */
                GetBurgIndexHCP(dir, 0, burgCnt, bList, &bIndex);

/*
 *              If the junction burgers vector is already in the list,
 *              this interaction may result in some new planes.  Add the
 *              new planes (if any) in now.
 */
                if (bIndex >= 0) {
                    int   numJunctPlanes, newListLen;
                    real8 junctPlaneList[HCP_MAX_PLANES_PER_BURG][3];
                    real8 newList[HCP_MAX_PLANES_PER_BURG][3];

                    numJunctPlanes = 0;

/*
 *                  Find all the junction planes that may be formed by
 *                  interactions of these two glide burgers vectors
 */
                    FindJctPlanes(glissilePlanesPerBurg[i],
                                  glissilePlanesPerBurg[j],
                                  &pList[burgPlaneIndex[i]],
                                  &pList[burgPlaneIndex[j]],
                                  dir, junctPlaneList, &numJunctPlanes);

/*
 *                  Merge the list of junction planes with the current
 *                  list of known planes for this dir
 */
                    MergePlaneLists(&pList[burgPlaneIndex[bIndex]],
                                    junctPlaneList, planesPerBurg[bIndex],
                                    numJunctPlanes, HCP_MAX_PLANES_PER_BURG,
                                    newList, &newListLen);

/*
 *                  Copy the merged plane list back into the section
 *                  of the full plane list associated with this burgers
 *                  vector
 */
                    planesPerBurg[bIndex] = newListLen;

                    for (k = 0; k < planesPerBurg[bIndex]; k++) {
                        VECTOR_COPY(pList[burgPlaneIndex[bIndex]+k],
                                    newList[k]);
                    }

                }  /* if (bIndex >= 0) */

/*
 *              If we did not find the burgers vector on the existing list,
 *              add it in.  
 */
                if (bIndex < 0) {
                    int   numJunctPlanes;
                    real8 junctPlaneList[HCP_MAX_PLANES_PER_BURG][3];

                    VECTOR_COPY(bList[burgCnt], dir);

/*
 *                  It may be possible for segments with this burgers vector
 *                  to be splintered into two identical segments with a
 *                  different glide burgers vector in a manner which is
 *                  energetically favorable.  If this is the case, add it
 *                  to the list for use when we check for splinterable
 *                  segments elsewhere in the code.
 */
                    VECTOR_COPY(tmpDir, dir);
                    NormalizeVec(tmpDir);

                    for (k = 0; k < HCP_NUM_GLIDE_BURG; k++) {

                        VECTOR_COPY(tmpBurg, bList[k]);
                        NormalizeVec(tmpBurg);

                        if (fabs(DotProduct(tmpDir, tmpBurg)) > 0.99995) {
                            real8              burgSign;
                            SplinterableBurg_t *sBurg;

                            sBurg = &splinterableBList[splinterableBurgCnt];

                            VECTOR_COPY(sBurg->refBurg, dir);

                            burgSign = (DotProduct(bList[k], dir) < 0.0) ?
                                       -1.0 : 1.0;

                            sBurg->newBurg[X] = burgSign * bList[k][X];
                            sBurg->newBurg[Y] = burgSign * bList[k][Y];
                            sBurg->newBurg[Z] = burgSign * bList[k][Z];


                            splinterGlideBurg[splinterableBurgCnt] = k;
                            splinterJunctBurg[splinterableBurgCnt] = burgCnt;

                            splinterableBurgCnt +=1;
                        }
                    }

/*
 *                  Find the possible planes associated with this direction
 *                  and add them to the full plane list.
 */
                    FindJctPlanes(glissilePlanesPerBurg[i],
                                  glissilePlanesPerBurg[j],
                                  &pList[burgPlaneIndex[i]],
                                  &pList[burgPlaneIndex[j]],
                                  dir, junctPlaneList, &numJunctPlanes);

                    planesPerBurg[burgCnt]  = numJunctPlanes;

                    for (k = 0; k < numJunctPlanes; k++) {
                        VECTOR_COPY(pList[burgPlaneIndex[burgCnt]+k],
                                    junctPlaneList[k]);
                    }

                    burgCnt += 1;

                }  /* if (bIndex < 0) */

/*
 *              Possible direction 2
 */ 
                dir[0] = bList[i][0] - bList[j][0];
                dir[1] = bList[i][1] - bList[j][1];
                dir[2] = bList[i][2] - bList[j][2];

/*              Check if this Burgers vector is already in the list
 *              Criteria are: burgers vector normals are the same
 *              and vectors are equivalent up to a sign.
 */
                GetBurgIndexHCP(dir, 0, burgCnt, bList, &bIndex);

/*
 *              If the burgers vector is already in the list,
 *              this interaction may bring new planes. Add these planes here.
 */
                if (bIndex >= 0) {
                    int   numJunctPlanes, newListLen;
                    real8 junctPlaneList[HCP_MAX_PLANES_PER_BURG][3];
                    real8 newList[HCP_MAX_PLANES_PER_BURG][3];

                    numJunctPlanes = 0;

/*
 *                  Find all the junction planes that may be formed by
 *                  interactions of these two glide burgers vectors
 */
                    FindJctPlanes(glissilePlanesPerBurg[i],
                                  glissilePlanesPerBurg[j],
                                  &pList[burgPlaneIndex[i]],
                                  &pList[burgPlaneIndex[j]],
                                  dir, junctPlaneList, &numJunctPlanes);

/*
 *                  Merge the list of junction planes with the current
 *                  list of known planes for this dir
 */
                    MergePlaneLists(&pList[burgPlaneIndex[bIndex]],
                                    junctPlaneList, planesPerBurg[bIndex],
                                    numJunctPlanes, HCP_MAX_PLANES_PER_BURG,
                                    newList, &newListLen);

/*
 *                  Copy the merged plane list back into the section
 *                  of the full plane list associated with this burgers
 *                  vector
 */
                    planesPerBurg[bIndex] = newListLen;

                    for (k = 0; k < planesPerBurg[bIndex]; k++) {
                        VECTOR_COPY(pList[burgPlaneIndex[bIndex]+k],
                                    newList[k]);
                    }

                }  /* if (bIndex >= 0) */

/*
 *              If we did not find the burgers vector on the existing list,
 *              add it in.  
 */
                if (bIndex < 0) {
                    int   numJunctPlanes;
                    real8 junctPlaneList[HCP_MAX_PLANES_PER_BURG][3];

                    VECTOR_COPY(bList[burgCnt], dir);

/*
 *                  It may be possible for segments with this burgers vector
 *                  to be splintered into two identical segments with a
 *                  different glide burgers vector in a manner which is
 *                  energetically favorable.  If this is the case, add it
 *                  to the list for use when we check for splinterable
 *                  segments elsewhere in the code.
 */
                    VECTOR_COPY(tmpDir, dir);
                    NormalizeVec(tmpDir);

                    for (k = 0; k < HCP_NUM_GLIDE_BURG; k++) {

                        VECTOR_COPY(tmpBurg, bList[k]);
                        NormalizeVec(tmpBurg);

                        if (fabs(DotProduct(tmpDir, tmpBurg)) > 0.99995) {
                            real8              burgSign;
                            SplinterableBurg_t *sBurg;

                            sBurg = &splinterableBList[splinterableBurgCnt];

                            VECTOR_COPY(sBurg->refBurg, dir);

                            burgSign = (DotProduct(bList[k], dir) < 0.0) ?
                                       -1.0 : 1.0;

                            sBurg->newBurg[X] = burgSign * bList[k][X];
                            sBurg->newBurg[Y] = burgSign * bList[k][Y];
                            sBurg->newBurg[Z] = burgSign * bList[k][Z];

                            splinterGlideBurg[splinterableBurgCnt] = k;
                            splinterJunctBurg[splinterableBurgCnt] = burgCnt;

                            splinterableBurgCnt +=1;
                        }
                    }
/*
 *                  Find the possible planes associated with this direction
 *                  and add them to the full plane list.
 */
                    FindJctPlanes(glissilePlanesPerBurg[i],
                                  glissilePlanesPerBurg[j],
                                  &pList[burgPlaneIndex[i]],
                                  &pList[burgPlaneIndex[j]],
                                  dir, junctPlaneList, &numJunctPlanes);

                    planesPerBurg[burgCnt]  = numJunctPlanes;

                    for (k = 0; k < numJunctPlanes; k++) {
                        VECTOR_COPY(pList[burgPlaneIndex[burgCnt]+k],
                                    junctPlaneList[k]);
                    }

		    burgCnt += 1;
			
		}  /* if (bIndex < 0) */

            }  /* for ( j = i + 1; ... ) */

        }  /* for (i = 0; ... ); */

/* 
 *      For any of the junction burgers vectors that are splinterable,
 *      we have to ensure that all the planes for that burgers vector
 *      are included in the plane list for the burgers vector into which
 *      it may be splintered.  (i.e. all planes for junction burgers 
 *      vector [2a 0 0] must also be on the list for glide burgers vector
 *      [a 0 0].
 */
        for (k = 0; k < splinterableBurgCnt; k++) {
            int   numGlidePlanes, numJunctPlanes, newPlaneListLen;
            int   firstGlidePlane, firstJunctPlane;
            real8 newPlaneList[HCP_MAX_PLANES_PER_BURG][3];

            numGlidePlanes = planesPerBurg[splinterGlideBurg[k]];
            numJunctPlanes = planesPerBurg[splinterJunctBurg[k]];

            firstGlidePlane = burgPlaneIndex[splinterGlideBurg[k]];
            firstJunctPlane = burgPlaneIndex[splinterJunctBurg[k]];

            MergePlaneLists(&pList[firstGlidePlane], &pList[firstJunctPlane],
                            numGlidePlanes, numJunctPlanes,
                            HCP_MAX_PLANES_PER_BURG, newPlaneList,
                            &newPlaneListLen);

            for (j = 0; j < newPlaneListLen; j++) {
                VECTOR_COPY(pList[firstGlidePlane+j], newPlaneList[j]);
            }

            planesPerBurg[splinterGlideBurg[k]] = newPlaneListLen;
        }


/* 
 *      Initialize arrays
 */
        *numBurg = HCP_NUM_TOT_BURG;
        *numPlanes = HCP_NUM_TOT_PLANES;
        *numSplinterableBurgs = splinterableBurgCnt;

        *burgList = bList;
        *planeList = pList;
        *splinterableBurgList = splinterableBList;

        *numPlanesPerBurg = planesPerBurg;
        *burgFirstPlaneIndex = burgPlaneIndex;
        *numGlissilePlanes = glissilePlanesPerBurg;
        *planeType = typeOfPlane;


#if 0
/*
 *      The following can be used for debugging:
 *      Print the Burgers vectors and planes for HCP
 *      crystals
 */
	printf("List of Burgers vectors and plane for c/a=%f\n",cOa);
        int n,p, firstplane;
        for (p = 0; p < *numBurg ; p++) 
	  {
	    printf("\nBurgers[%d]=%f %f %f\n",p,bList[p][0],bList[p][1],bList[p][2]);
	    firstplane = burgPlaneIndex[p];
	    for (n = firstplane; n < firstplane+planesPerBurg[p]; n++) 
	      {
		printf("Plane[%d]=%f %f %f [type=%d]\n",n,pList[n][0],pList[n][1],pList[n][2],typeOfPlane[n]);
	      }
	  }
	exit(0);
#endif

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    RotateSystem
 *      Description: 
 *
 *      Arguments:
 *          IN:     XX
 *          IN:     YY
 *          IN:     ZZ
 *          IN/OUT: r1
 *          IN/OUT: r2
 *          IN/OUT: r3
 *          OUT:    Xp
 *          OUT:    Yp
 *          OUT:    Zp
 *
 *------------------------------------------------------------------------*/
void RotateSystem(real8 XX[6], real8 YY[6], real8 ZZ[6],
                  real8 r1[3], real8 r2[3], real8 r3[3],
                  real8 Xp[6], real8 Yp[6], real8 Zp[6])
{
        int   i, n;
        real8 A[3][3];

        n = 6;

        NormalizeVec(r1);
        NormalizeVec(r2);
        NormalizeVec(r3);

        VECTOR_COPY(A[0], r1);
        VECTOR_COPY(A[1], r2);
        VECTOR_COPY(A[2], r3);

        for (i = 0; i < n; i++) {
            Xp[i] = A[0][0] * XX[i] + A[0][1] * YY[i] + A[0][2] * ZZ[i];
            Yp[i] = A[1][0] * XX[i] + A[1][1] * YY[i] + A[1][2] * ZZ[i];
            Zp[i] = A[2][0] * XX[i] + A[2][1] * YY[i] + A[2][2] * ZZ[i];
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    FindGlidePlaneHCP()
 *      Description: Given the index of a burgers vector in the list of
 *                   HCP burgers vectors, find the index in the plane list
 *                   of the plane associed with the burgers vector and
 *                   matching the one provided by the caller.
 *
 *                   Even glide burgers vectors have sessile planes.
 *
 *                   NOTE: This function can be used in mobility_law_hcp_*
 *
 *      Arguments:
 *          IN:   bIndex  Index into the list of HCP burgers vector
 *                        stored in <home->burgData.burgList>
 *          IN:   plane   Plane to look up
 *          OUT:  pIndex  Index into <home->burgData.planeList> which
 *                        is returned to caller.
 *
 *------------------------------------------------------------------------*/
void FindGlidePlaneHCP(Home_t *home, int bIndex, real8 plane[3], int *pIndex)
{
    int     numGlidePlanes, firstGlidePlaneIndex;
    int    *numPlanesPerBurg, *burgFirstPlaneIndex;
    real8   maxVal;
    real8  *planeListXplane;
    real8 (*planeList)[3];

    planeList = home->burgData.planeList;
    numPlanesPerBurg = home->burgData.numPlanesPerBurg;
    burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;

    numGlidePlanes  = numPlanesPerBurg[bIndex];
    planeListXplane = (real8 *)calloc(1, numGlidePlanes * sizeof(real8));

    firstGlidePlaneIndex = burgFirstPlaneIndex[bIndex];
    MatrixMult(&planeList[firstGlidePlaneIndex][0],
               numGlidePlanes, 3, 3,
               plane, 1, 1,
               planeListXplane, 1);

    FindAbsMax(planeListXplane, numGlidePlanes, &maxVal, pIndex);

    free(planeListXplane);

    return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    HCPEcoreFactor()
 *      Description: For HCP simulations only, determine the Burgers vector
 *                   index and return the corresponding Ecore prefactor value.
 *
 *      Arguments:
 *          IN:   burg  burgers vector specified in crystalographic frame
 *
 *------------------------------------------------------------------------*/
real8 HCPEcoreFactor(Home_t *home, real8 burg[3])
{
    int   bIndex;
    real8 (*burgList)[3];
    real8 Ecore;
    Param_t *param;

    param = home->param;
  
    burgList = home->burgData.burgList;
  
    GetBurgIndexHCP(burg, 0, HCP_NUM_TOT_BURG , burgList, &bIndex);
  
    if (bIndex < 3) {
        // <a> type burgers vector
        Ecore = param->HCPEcoreA;
    } else if (bIndex == 9) {
        // <c> type burgers vector
        Ecore = param->HCPEcoreC;
    } else if (bIndex < 9) {
        // Any burgers vector with an index > 3 but less than 9 is
        // a <c+a> type burgers vector.
        Ecore = param->HCPEcoreCpA;
    } else {
        // All other types of burgers vectors are junction Burgers vectors.
        // For now these are assigned the same Ecore as <a> types.  Is this
        // correct???
        Ecore = param->HCPEcoreA;
    }
    
    return(Ecore);
}
