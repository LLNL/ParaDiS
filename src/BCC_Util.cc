/****************************************************************************
 *
 *      Module:         BCC_Util.c
 *
 *      Author:         S. Aubry January 2017
 *
 *      Description:    This module contains various functions required
 *                      to list the burgers vectors and associated
 *                      glide planes for a BCC material
 *
 *      Includes public functions:
 *
 *          GetBurgList_BCC()
 *          GetBurgIndexBCC()
 *          FindGlidePlaneBCC()
 *
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Matrix.h"

/*-------------------------------------------------------------------------
 *
 *      Function:    GetBurgList_BCC()
 *
 *      Description: Calculate all the unique burgers vectors for an BCC
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
void GetBurgList_BCC(real8 (**burgList)[3], real8 (**planeList)[3],
                     int **numPlanesPerBurg, int **burgFirstPlaneIndex,
                     int **numGlissilePlanes, int **planeType,
                     int *numBurg, int *numPlanes)
{
        int   bufSize;
        int   *typeOfPlane, *glissilePlanesPerBurg;
        int   *planesPerBurg, *burgPlaneIndex;
        real8 (*bList)[3];
        real8 (*pList)[3];

        // The values in the reference burgers vector array are
        // normalized burgers vectors.

        static int pPerBurg   [] = { 6, 6, 6, 6, 16, 16, 16, 12, 12, 12, 12, 12, 12};
        static int bPlaneIndex[] = { 0, 6, 12, 18, 24, 40, 56, 72, 84, 96, 108, 120, 132};

        static real8 burgRef[][3] = 
        {
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1]
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1]
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1]
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1]
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1]
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0]
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0]
            {  0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1]
            {  0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1]
            {  0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0]
            {  0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1]
            {  0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1]
            {  0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0]
        };

        static real8 planeRef[][3] = 
        {
            // burgRef[0][*]...
            {  0.0               ,  -0.7071067811865475,   0.7071067811865475 },   // [ 0 -1  1]
            {  0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1]
            { -0.7071067811865475,   0.7071067811865475,   0.0                },   // [-1  1  0]
            {  0.8164965809277261,  -0.4082482904638631,  -0.4082482904638631 },   // [ 2 -1 -1]
            { -0.4082482904638631,   0.8164965809277261,  -0.4082482904638631 },   // [-1  2 -1]
            { -0.4082482904638631,  -0.4082482904638631,   0.8164965809277261 },   // [-1 -1  2]

            // burgRef[1][*]...
            {  0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1]
            { -0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [-1  0 -1]
            {  0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0]
            {  0.8164965809277261,  -0.4082482904638631,   0.4082482904638631 },   // [ 2 -1  1]
            { -0.4082482904638631,   0.8164965809277261,   0.4082482904638631 },   // [-1  2  1]
            { -0.4082482904638631,  -0.4082482904638631,  -0.8164965809277261 },   // [-1 -1 -2]

            // burgRef[2][*]...
            {  0.0               ,  -0.7071067811865475,  -0.7071067811865475 },   // [ 0 -1 -1]
            { -0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [-1  0  1]
            {  0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0]
            {  0.8164965809277261,   0.4082482904638631,  -0.4082482904638631 },   // [ 2  1 -1]
            { -0.4082482904638631,  -0.8164965809277261,  -0.4082482904638631 },   // [-1 -2 -1]
            { -0.4082482904638631,   0.4082482904638631,   0.8164965809277261 },   // [-1  1  2]

            // burgRef[3][*]...
            {  0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1]
            {  0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1]
            { -0.7071067811865475,  -0.7071067811865475,   0.0                },   // [-1 -1  0]
            { -0.8164965809277261,  -0.4082482904638631,  -0.4082482904638631 },   // [-2 -1 -1]
            {  0.4082482904638631,   0.8164965809277261,  -0.4082482904638631 },   // [ 1  2 -1]
            {  0.4082482904638631,  -0.4082482904638631,   0.8164965809277261 },   // [ 1 -1  2]

            // burgRef[4][*]...
            {  0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0]
            {  0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0]
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0]
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0]
            {  0.3162277660168379,  -0.9486832980505138,   0.0                },   // [ 1 -3  0]
            {  0.3162277660168379,   0.9486832980505138,   0.0                },   // [ 1  3  0]
            { -0.9486832980505138,   0.3162277660168379,   0.0                },   // [-3  1  0]
            {  0.9486832980505138,   0.3162277660168379,   0.0                },   // [ 3  1  0]
            {  0.8944271909999159,  -0.4472135954999579,   0.0                },   // [ 2 -1  0]
            { -0.4472135954999579,   0.8944271909999159,   0.0                },   // [-1  2  0]
            {  0.9805806756909202,  -0.1961161351381840,   0.0                },   // [ 5 -1  0]
            { -0.1961161351381840,   0.9805806756909202,   0.0                },   // [-1  5  0]
            {  0.8944271909999159,   0.4472135954999579,   0.0                },   // [ 2  1  0]
            {  0.4472135954999579,   0.8944271909999159,   0.0                },   // [ 1  2  0]
            {  0.9805806756909202,   0.1961161351381840,   0.0                },   // [ 5  1  0]
            {  0.1961161351381840,   0.9805806756909202,   0.0                },   // [ 1  5  0]

            // burgRef[5][*]...
            {  0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1]
            {  0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1]
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1]
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0]
            {  0.3162277660168379,   0.0               ,   0.9486832980505138 },   // [ 1  0  3]
            { -0.3162277660168379,   0.0               ,   0.9486832980505138 },   // [-1  0  3]
            { -0.9486832980505138,   0.0               ,   0.3162277660168379 },   // [-3  0  1]
            {  0.9486832980505138,   0.0               ,   0.3162277660168379 },   // [ 3  0  1]
            {  0.8944271909999159,   0.0               ,  -0.4472135954999579 },   // [ 2  0 -1]
            { -0.4472135954999579,   0.0               ,   0.8944271909999159 },   // [-1  0  2]
            {  0.9805806756909202,   0.0               ,  -0.1961161351381840 },   // [ 5  0 -1]
            {  0.1961161351381840,   0.0               ,  -0.9805806756909202 },   // [ 1  0 -5]
            {  0.8944271909999159,   0.0               ,   0.4472135954999579 },   // [ 2  0  1]
            {  0.4472135954999579,   0.0               ,   0.8944271909999159 },   // [ 1  0  2]
            {  0.9805806756909202,   0.0               ,   0.1961161351381840 },   // [ 5  0  1]
            {  0.1961161351381840,   0.0               ,   0.9805806756909202 },   // [ 1  0  5]

            // burgRef[6][*]...
            {  0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1]
            {  0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1]
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1]
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0]
            {  0.0               ,   0.3162277660168379,   0.9486832980505138 },   // [ 0  1  3]
            {  0.0               ,  -0.3162277660168379,   0.9486832980505138 },   // [ 0 -1  3]
            {  0.0               ,   0.9486832980505138,   0.3162277660168379 },   // [ 0  3  1]
            {  0.0               ,   0.9486832980505138,  -0.3162277660168379 },   // [ 0  3 -1]
            {  0.0               ,  -0.9805806756909202,   0.1961161351381840 },   // [ 0 -5  1]
            {  0.0               ,   0.1961161351381840,  -0.9805806756909202 },   // [ 0  1 -5]
            {  0.0               ,   0.9805806756909202,   0.1961161351381840 },   // [ 0  5  1]
            {  0.0               ,   0.1961161351381840,   0.9805806756909202 },   // [ 0  1  5]
            {  0.0               ,   0.8944271909999159,  -0.4472135954999579 },   // [ 0  2 -1]
            {  0.0               ,  -0.4472135954999579,   0.8944271909999159 },   // [ 0 -1  2]
            {  0.0               ,   0.8944271909999159,   0.4472135954999579 },   // [ 0  2  1]
            {  0.0               ,   0.4472135954999579,   0.8944271909999159 },   // [ 0  1  2]

            // burgRef[7][*]...
            {  0.0               ,   0.7071067811865475,   0.7071067811865475 },   // [ 0  1  1]
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0]
            {  0.8164965809277261,   0.4082482904638631,   0.4082482904638631 },   // [ 2  1  1]
            { -0.8164965809277261,   0.4082482904638631,   0.4082482904638631 },   // [-2  1  1]
            {  0.9428090415820635,   0.2357022603955159,   0.2357022603955159 },   // [ 4  1  1]
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1]
            { -0.9428090415820635,   0.2357022603955159,   0.2357022603955159 },   // [-4  1  1]
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1]
            {  0.6859943405700353,   0.5144957554275265,   0.5144957554275265 },   // [ 4  3  3]
            { -0.6859943405700353,   0.5144957554275265,   0.5144957554275265 },   // [-4  3  3]
            {  0.4264014327112208,   0.6396021490668312,   0.6396021490668312 },   // [ 4  6  6]
            { -0.4264014327112208,   0.6396021490668312,   0.6396021490668312 },   // [-4  6  6]

            // burgRef[8][*]...
            {  0.7071067811865475,   0.0               ,   0.7071067811865475 },   // [ 1  0  1]
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0]
            {  0.4082482904638631,   0.8164965809277261,   0.4082482904638631 },   // [ 1  2  1]
            {  0.4082482904638631,  -0.8164965809277261,   0.4082482904638631 },   // [ 1 -2  1]
            {  0.2357022603955159,   0.9428090415820635,   0.2357022603955159 },   // [ 1  4  1]
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1]
            {  0.2357022603955159,  -0.9428090415820635,   0.2357022603955159 },   // [ 1 -4  1]
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1]
            {  0.5144957554275265,   0.6859943405700353,   0.5144957554275265 },   // [ 3  4  3]
            {  0.5144957554275265,  -0.6859943405700353,   0.5144957554275265 },   // [ 3 -4  3]
            {  0.6396021490668312,   0.4264014327112208,   0.6396021490668312 },   // [ 6  4  6]
            {  0.6396021490668312,  -0.4264014327112208,   0.6396021490668312 },   // [ 6 -4  6]

            // burgRef[9][*]...
            {  0.7071067811865475,   0.7071067811865475,   0.0                },   // [ 1  1  0]
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1]
            {  0.4082482904638631,   0.4082482904638631,   0.8164965809277261 },   // [ 1  1  2]
            {  0.4082482904638631,   0.4082482904638631,  -0.8164965809277261 },   // [ 1  1 -2]
            {  0.2357022603955159,   0.2357022603955159,   0.9428090415820635 },   // [ 1  1  4]
            {  0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [ 1  1  1]
            {  0.2357022603955159,   0.2357022603955159,  -0.9428090415820635 },   // [ 1  1 -4]
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1]
            {  0.5144957554275265,   0.5144957554275265,   0.6859943405700353 },   // [ 3  3  4]
            {  0.5144957554275265,   0.5144957554275265,  -0.6859943405700353 },   // [ 3  3 -4]
            {  0.6396021490668312,   0.6396021490668312,   0.4264014327112208 },   // [ 6  6  4]
            {  0.6396021490668312,   0.6396021490668312,  -0.4264014327112208 },   // [ 6  6 -4]

            // burgRef[10][*]...
            {  0.0               ,   0.7071067811865475,  -0.7071067811865475 },   // [ 0  1 -1]
            {  1.0               ,   0.0               ,   0.0                },   // [ 1  0  0]
            {  0.8164965809277261,   0.4082482904638631,  -0.4082482904638631 },   // [ 2  1 -1]
            {  0.8164965809277261,  -0.4082482904638631,   0.4082482904638631 },   // [ 2 -1  1]
            {  0.9428090415820635,  -0.2357022603955159,   0.2357022603955159 },   // [ 4 -1  1]
            {  0.9428090415820635,   0.2357022603955159,  -0.2357022603955159 },   // [ 4  1 -1]
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1]
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1]
            {  0.6859943405700353,   0.5144957554275265,  -0.5144957554275265 },   // [ 4  3 -3]
            {  0.6859943405700353,  -0.5144957554275265,   0.5144957554275265 },   // [ 4 -3  3]
            {  0.4264014327112208,  -0.6396021490668312,   0.6396021490668312 },   // [ 4 -6  6]
            {  0.4264014327112208,   0.6396021490668312,  -0.6396021490668312 },   // [ 4  6 -6]

            // burgRef[11][*]...
            {  0.7071067811865475,   0.0               ,  -0.7071067811865475 },   // [ 1  0 -1]
            {  0.0               ,   1.0               ,   0.0                },   // [ 0  1  0]
            {  0.4082482904638631,   0.8164965809277261,  -0.4082482904638631 },   // [ 1  2 -1]
            { -0.4082482904638631,   0.8164965809277261,   0.4082482904638631 },   // [-1  2  1]
            {  0.5773502691896258,   0.5773502691896258,  -0.5773502691896258 },   // [ 1  1 -1]
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1]
            { -0.2357022603955159,   0.9428090415820635,   0.2357022603955159 },   // [-1  4  1]
            {  0.2357022603955159,   0.9428090415820635,  -0.2357022603955159 },   // [ 1  4 -1]
            {  0.6396021490668312,   0.4264014327112208,  -0.6396021490668312 },   // [ 6  4 -6]
            { -0.6396021490668312,   0.4264014327112208,   0.6396021490668312 },   // [-6  4  6]
            {  0.5144957554275265,   0.6859943405700353,  -0.5144957554275265 },   // [ 3  4 -3]
            { -0.5144957554275265,   0.6859943405700353,   0.5144957554275265 },   // [-3  4  3]

            // burgRef[12][*]...
            {  0.7071067811865475,  -0.7071067811865475,   0.0                },   // [ 1 -1  0]
            {  0.0               ,   0.0               ,   1.0                },   // [ 0  0  1]
            {  0.4082482904638631,  -0.4082482904638631,   0.8164965809277261 },   // [ 1 -1  2]
            { -0.4082482904638631,   0.4082482904638631,   0.8164965809277261 },   // [-1  1  2]
            {  0.5773502691896258,  -0.5773502691896258,   0.5773502691896258 },   // [ 1 -1  1]
            { -0.5773502691896258,   0.5773502691896258,   0.5773502691896258 },   // [-1  1  1]
            {  0.2357022603955159,  -0.2357022603955159,   0.9428090415820635 },   // [ 1 -1  4]
            { -0.2357022603955159,   0.2357022603955159,   0.9428090415820635 },   // [-1  1  4]
            { -0.6396021490668312,   0.6396021490668312,   0.4264014327112208 },   // [-6  6  4]
            {  0.6396021490668312,  -0.6396021490668312,   0.4264014327112208 },   // [ 6 -6  4]
            { -0.5144957554275265,   0.5144957554275265,   0.6859943405700353 },   // [-3  3  4]
            {  0.5144957554275265,  -0.5144957554275265,   0.6859943405700353 }    // [ 3 -3  4]

        };

        // Allocate some of the dynamic arrays... sizes are known
        // ahead of time.

        bufSize        = BCC_NUM_TOT_BURG * sizeof(double [3]);
        bList          = (double (*)[3])calloc(1, bufSize);

        bufSize        = BCC_NUM_TOT_PLANES * sizeof(double [3]);
        pList          = (double (*)[3])calloc(1, bufSize);

        bufSize        = BCC_NUM_TOT_BURG * sizeof(double [3]);
        planesPerBurg  = (int *)calloc(1, bufSize * sizeof(int));
        burgPlaneIndex = (int *)calloc(1, bufSize * sizeof(int));

        bufSize               = BCC_NUM_TOT_PLANES * sizeof(double [3]);
        glissilePlanesPerBurg = (int *)calloc(1, bufSize * sizeof(int));
        typeOfPlane           = (int *)calloc(1, bufSize * sizeof(int));

        for (int k=0; k < BCC_NUM_TOT_BURG; k++) 
        {
            VECTOR_COPY(bList[k],burgRef[k]);
            burgPlaneIndex[k] = bPlaneIndex[k];
            planesPerBurg[k] = pPerBurg[k];
        }

        for (int k=0; k < BCC_NUM_TOT_PLANES; k++) 
        {
            VECTOR_COPY(pList[k],planeRef[k]);
            typeOfPlane[k] = 0;
            glissilePlanesPerBurg[k] = 0;
        }

        for (int k=0; k < 4; k++) 
        {   glissilePlanesPerBurg[k] = 6; }

        // Set the plane type (character) for every plane in the
        // list.  Type will be set to one of the following:
        //
        //     1 if plane is glissile
        //     0 if plane is sessile

        for (int k=0; k < 24; k++)
        {   typeOfPlane[k] = 1; }

        // Initialize arrays

        *numBurg = BCC_NUM_TOT_BURG;
        *numPlanes = BCC_NUM_TOT_PLANES;

        *burgList = bList;
        *planeList = pList;

        *numPlanesPerBurg = planesPerBurg;
        *burgFirstPlaneIndex = burgPlaneIndex;
        *numGlissilePlanes = glissilePlanesPerBurg;
        *planeType = typeOfPlane;

#if 0
        // The following can be used for debugging:
        // Print the Burgers vectors and planes for BCC crystals

        printf("List of Burgers vectors and plane for BCC\n");
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
 *      Function:    GetBurgIndexBCC
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
int GetBurgIndexBCC(real8 burgVec[3], int numBurg, real8 (*burgList)[3])
{
        int burgIndex;
        real8 maxVal=-1.0, nb;
        real8  b[3], *burgListX;

        nb = Normal(burgVec);
        b[X] = burgVec[X]/nb;
        b[Y] = burgVec[Y]/nb;
        b[Z] = burgVec[Z]/nb;

        burgListX = (real8 *)calloc(1, numBurg * sizeof(real8));

        MatrixMult(&burgList[0][0],
                    numBurg, 3, 3,
                    b, 1, 1,
                    burgListX, 1);

        FindAbsMax(burgListX,numBurg,&maxVal,&burgIndex);
        free(burgListX);

        return burgIndex;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    FindGlidePlaneBCC()
 *      Description: Given the index of a burgers vector in the list of
 *                   BCC burgers vectors, find the index in the plane list
 *                   of the plane associed with the burgers vector and
 *                   matching the one provided by the caller.
 *
 *      Arguments:
 *          IN:   bIndex  Index into the list of BCC burgers vector
 *                        stored in <home->burgData.burgList>
 *          IN:   plane   Plane to look up
 *          OUT:  pIndex  Index into <home->burgData.planeList> which
 *                        is returned to caller.
 *
 *------------------------------------------------------------------------*/
void FindGlidePlaneBCC(Home_t *home, int bIndex, real8 plane[3], int *pIndex)
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

void GetBCCAllIndices(Home_t *home, real8 b[3], real8 p[3],
                      int *bIndex, int *pIndex, int *gIndex)
{
   *bIndex = -1;
   *pIndex = -1;
   *gIndex = -1;

   int firstPlaneIndex;
   real8 maxVal;
   int  *burgFirstPlaneIndex;
   real8 (*burgList)[3], (*planeList)[3];
   burgList            = home->burgData.burgList;
   planeList           = home->burgData.planeList;
   burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;

   *bIndex = GetBurgIndexBCC(b, BCC_NUM_TOT_BURG , burgList);

   if (*bIndex < BCC_NUM_GLIDE_BURG)
   {
      int numPlanes = 6;
      real8 *planeListXplane;
      firstPlaneIndex = burgFirstPlaneIndex[*bIndex];
      planeListXplane = (real8 *)calloc(1, numPlanes * sizeof(real8));

      MatrixMult(&planeList[firstPlaneIndex][0],
                 numPlanes, 3, 3,
                 p    , 1, 1,
                 planeListXplane, 1);
      FindAbsMax(planeListXplane, numPlanes, &maxVal, pIndex);

      *pIndex += firstPlaneIndex;

      real8 g[3];
      NormalizedCrossVector(p,b,g);

      MatrixMult(&planeList[firstPlaneIndex][0],
                 numPlanes, 3, 3,
                 g    , 1, 1,
                 planeListXplane, 1);
      FindAbsMax(planeListXplane, numPlanes, &maxVal, gIndex);

      *gIndex += firstPlaneIndex;


      free(planeListXplane);
   }

   // NEW: identify plane index for [100] binary junctions (bertin1)

   else if (*bIndex >= 4 && *bIndex <= 6)
   {
      int numPlanes = 4;        // only consider [010] and [011] planes
      real8 *planeListXplane;
      firstPlaneIndex = burgFirstPlaneIndex[*bIndex];
      planeListXplane = (real8 *)calloc(1, numPlanes * sizeof(real8));

      MatrixMult(&planeList[firstPlaneIndex][0],
                 numPlanes, 3, 3,
                 p    , 1, 1,
                 planeListXplane, 1);
      FindAbsMax(planeListXplane, numPlanes, &maxVal, pIndex);

      *pIndex += firstPlaneIndex;

      free(planeListXplane);
   }

}
