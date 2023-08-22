/****************************************************************************
 *
 *      Module:         FCC_Util.c
 *
 *      Author:         S. Aubry June 2016
 *
 *      Description:    This module contains various functions required
 *                      to list the burgers vectors and associated
 *                      glide planes for a FCC material
 *
 *      Includes public functions:
 *
 *          GetBurgList_FCC()
 *          GetBurgIndexFCC()
 *          FindGlidePlaneFCC()
 *
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Matrix.h"

/*-------------------------------------------------------------------------
 *
 *      Function:    GetBurgList_FCC()
 *
 *      Description: Calculate all the unique burgers vectors for an FCC
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
 *------------------------------------------------------------------------*/
void GetBurgList_FCC(real8 (**burgList)[3], real8 (**planeList)[3],
                     int **numPlanesPerBurg, int **burgFirstPlaneIndex,
                     int **numGlissilePlanes, int **planeType,
                     int *numBurg, int *numPlanes)
{

        int   bufSize;
        int   *typeOfPlane, *glissilePlanesPerBurg;
        int   *planesPerBurg, *burgPlaneIndex;
        real8 (*bList)[3];
        real8 (*pList)[3];

        static int pPerBurg   [21] = { 3, 3, 3, 3, 3, 3, 4, 4, 4, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2, 2 };
        static int bPlaneIndex[21] = { 0, 3, 6, 9, 12, 15, 18, 22, 26, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52 };

        // The values in the reference burgers vector array are
        // normalized burgers vectors.

        static real8 burgRef[][3] =
        {
           // glide and junction burgers vectors
           {   0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0]
           {   0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0]
           {   0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1]
           {   0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1]
           {   0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1]
           {   0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1]

           // other likely junction burgers vectors
           {   1.0               ,   0.0               ,   0.0                },   // [ 1  0  0]
           {   0.0               ,   1.0               ,   0.0                },   // [ 0  1  0]
           {   0.0               ,   0.0               ,   1.0                },   // [ 0  0  1]

           /* unlikely junction burgers vectors */
           {   0.8164965809277261,   0.4082482904638631,   0.4082482904638631 },   // [ 2  1  1]
           {  -0.8164965809277261,   0.4082482904638631,   0.4082482904638631 },   // [-2  1  1]
           {   0.8164965809277261,  -0.4082482904638631,   0.4082482904638631 },   // [ 2 -1  1]
           {   0.8164965809277261,   0.4082482904638631,  -0.4082482904638631 },   // [ 2  1 -1]
           {   0.4082482904638631,   0.8164965809277261,   0.4082482904638631 },   // [ 1  2  1]
           {  -0.4082482904638631,   0.8164965809277261,   0.4082482904638631 },   // [-1  2  1]
           {   0.4082482904638631,  -0.8164965809277261,   0.4082482904638631 },   // [ 1 -2  1]
           {   0.4082482904638631,   0.8164965809277261,  -0.4082482904638631 },   // [ 1  2 -1]
           {   0.4082482904638631,   0.4082482904638631,   0.8164965809277261 },   // [ 1  1  2]
           {  -0.4082482904638631,   0.4082482904638631,   0.8164965809277261 },   // [-1  1  2]
           {   0.4082482904638631,  -0.4082482904638631,   0.8164965809277261 },   // [ 1 -1  2]
           {   0.4082482904638631,   0.4082482904638631,  -0.8164965809277261 },   // [ 1  1 -2]
        };

        static real8 glideRef[][3] =
        {
            // burgRef[0][*]...
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1] (common glide plane)
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1] (common glide plane)
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1] (junction plane)

            // burgRef[1][*]...
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1] (common glide plane)
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1] (common glide plane)
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1] (junction plane)

            // burgRef[2][*]...
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1] (common glide plane)
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1] (common glide plane)
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0] (junction plane)

            // burgRef[3][*]...
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1] (common glide plane)
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1] (common glide plane)
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0] (junction plane)

            // burgRef[4][*]...
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1] (common glide plane)
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1] (common glide plane)
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0] (junction plane)

            // burgRef[5][*]...
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1] (common glide plane)
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1] (common glide plane)
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0] (junction plane)

            // burgRef[6][*]...
            {  0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1] (junction plane)
            {  0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1] (junction plane)
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0] (junction plane)
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1] (junction plane)

            // burgRef[7][*]...
            {  0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1] (junction plane)
            {  0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1] (junction plane)
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0] (junction plane)
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1] (junction plane)

            // burgRef[8][*]...
            {  0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0] (junction plane)
            {  0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0] (junction plane)
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0] (junction plane)
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0] (junction plane)

            // burgRef[9][*]...
            {  0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1]
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1]

            // burgRef[10][*]...
            {  0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1]
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1]

            // burgRef[11][*]...
            {  0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1]
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1]

            // burgRef[12][*]...
            {  0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1]
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1]

            // burgRef[13][*]...
            {  0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1]
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1]

            // burgRef[14][*]...
            {  0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1]
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1]

            // burgRef[15][*]...
            {  0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1]
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1]

            // burgRef[16][*]...
            {  0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1]
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1]

            // burgRef[17][*]...
            {  0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0]
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1]

            // burgRef[18][*]...
            {  0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0]
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1]

            // burgRef[19][*]...
            {  0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0]
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1]

            // burgRef[20][*]...
            {  0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0]
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 }    // [ 1  1  1]
        };

        // Allocate some of the dynamic arrays... sizes are known
        // ahead of time.

        bufSize               = FCC_NUM_TOT_BURG   * sizeof(double [3]);
        bList                 = (double (*)[3]) calloc(1, bufSize);

        bufSize               = FCC_NUM_TOT_PLANES * sizeof(double [3]);
        pList                 = (double (*)[3]) calloc(1, bufSize);

        bufSize               = FCC_NUM_TOT_BURG   * sizeof(double [3]);
        planesPerBurg         = (int *) calloc(1, bufSize * sizeof(int));
        burgPlaneIndex        = (int *) calloc(1, bufSize * sizeof(int));

        bufSize               = FCC_NUM_TOT_PLANES * sizeof(double [3]);
        glissilePlanesPerBurg = (int *) calloc(1, bufSize * sizeof(int));
        typeOfPlane           = (int *) calloc(1, bufSize * sizeof(int));

        for (int k=0; k < FCC_NUM_TOT_BURG; k++)
        {
            VECTOR_COPY(bList[k],burgRef[k]);
            burgPlaneIndex[k] = bPlaneIndex[k];
            planesPerBurg[k] = pPerBurg[k];
        }

        for (int k=0; k < FCC_NUM_TOT_PLANES; k++)
        {
            VECTOR_COPY(pList[k],glideRef[k]);
            typeOfPlane[k] = 0;
            glissilePlanesPerBurg[k] = 0;
        }

        for (int k=0; k < 6; k++)
        {
            glissilePlanesPerBurg[k] = 2;
        }

        // Set the plane type (character) for every plane in the
        // list.  Type will be set to one of the following:
        //
        //     1 if plane is glissile
        //     0 if plane is sessile

        typeOfPlane[ 0] = 1;
        typeOfPlane[ 1] = 1;

        typeOfPlane[ 3] = 1;
        typeOfPlane[ 4] = 1;

        typeOfPlane[ 6] = 1;
        typeOfPlane[ 7] = 1;

        typeOfPlane[ 9] = 1;
        typeOfPlane[10] = 1;

        typeOfPlane[12] = 1;
        typeOfPlane[13] = 1;

        typeOfPlane[15] = 1;
        typeOfPlane[16] = 1;

        // Initialize arrays

        *numBurg   = FCC_NUM_TOT_BURG;
        *numPlanes = FCC_NUM_TOT_PLANES;

        *burgList  = bList;
        *planeList = pList;

        *numPlanesPerBurg    = planesPerBurg;
        *burgFirstPlaneIndex = burgPlaneIndex;
        *numGlissilePlanes   = glissilePlanesPerBurg;
        *planeType           = typeOfPlane;

#if 0
        // The following can be used for debugging:
        // Print the Burgers vectors and planes for FCC crystals

        printf("List of Burgers vectors and plane for FCC\n");
        int n,p, firstplane;
        for (p=0; p < *numBurg ; p++)
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
 *      Function:    GetBurgIndexFCC
 *
 *      Description: Find the index in the reference burgers vector list
 *                   of the burgers vector matching the one provided by
 *                   the caller.  Criteria for matching burgers vectors
 *                   are 1) burgers vector normals are the same and 2) vectors
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
void GetBurgIndexFCC(real8 burgVec[3], int numBurg,
                     real8 (*burgList)[3], int *burgIndex)
{
        real8 maxVal, nb;
        real8  b[3], *burgListX;

        nb   = Normal(burgVec);
        b[X] = burgVec[X]/nb;
        b[Y] = burgVec[Y]/nb;
        b[Z] = burgVec[Z]/nb;

        burgListX = (real8 *)calloc(1, numBurg * sizeof(real8));

        MatrixMult(&burgList[0][0],
        	        numBurg, 3, 3,
        	        b, 1, 1,
        	        burgListX, 1);

        FindAbsMax(burgListX, numBurg, &maxVal,burgIndex);

        free(burgListX);
        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    FindGlidePlaneFCC()
 *      Description: Given the index of a burgers vector in the list of
 *                   FCC burgers vectors, find the index in the plane list
 *                   of the plane associed with the burgers vector and
 *                   matching the one provided by the caller.
 *
 *      Arguments:
 *          IN:   bIndex  Index into the list of FCC burgers vector
 *                        stored in <home->burgData.burgList>
 *          IN:   plane   Plane to look up
 *          OUT:  pIndex  Index into <home->burgData.planeList> which
 *                        is returned to caller.
 *
 *------------------------------------------------------------------------*/
void FindGlidePlaneFCC(Home_t *home, int bIndex, real8 plane[3], int *pIndex)
{
    int     numGlidePlanes, firstGlidePlaneIndex;
    int    *numPlanesPerBurg, *burgFirstPlaneIndex;
    real8   maxVal;
    real8  *planeListXplane;
    real8 (*planeList)[3];

    planeList           = home->burgData.planeList;
    numPlanesPerBurg    = home->burgData.numPlanesPerBurg;
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

