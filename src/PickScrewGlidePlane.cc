/***************************************************************************
 *
 *      Module:       PickScrewGlidePlane.c
 *      Description:  Contains functions needed to select an appropriate
 *                    glide plane for a newly created screw dislocation
 *                    based on the burgers vector and crystal structure.
 *
 *      Includes:
 *
 *          PickBCCScrewGlidePlane() static func
 *          PickFCCScrewGlidePlane() static func
 *          PickHCPScrewGlidePlane() static func
 *          PickRhomboVaScrewGlidePlane() static func
 *          PickScrewGlidePlane()
 *
 **************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

/*-------------------------------------------------------------------------
 *
 *      Function:     PickBCCScrewGlidePlane
 *      Description:  Selects, at random, an appropriate glide plane for
 *                    BCC screw dislocations based on the burgers vector.
 *    
 *      Changed S.A. Jan 31, 2017 : 6 glide planes instead of 3 for BCC
 *
 *      Arguments:
 *          burgVec      Input array containing the 3 components of the
 *                       burgers vector. This burgers vector should 
 *                       already have been rotated from the laboratory
 *                       frame to crystal frame by the caller.
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.  Returned glide plane is in the 
 *                       crystalographic frame.
 *
 *------------------------------------------------------------------------*/
static void PickBCCScrewGlidePlane(Home_t *home, real8 burgVec[3], real8 glidePlane[3])
{
        static int   seed = 8917346;
        int          numBurg, burgIndex, planeIndex, numGlidePlanes;
        int          *numGlissilePlanesPerBurg, *burgFirstPlaneIndex;
        real8        randVal;
        real8        (*burgList)[3], (*planeList)[3];
#ifdef _OPENMP
#pragma omp threadprivate (seed)
#endif


        VECTOR_ZERO(glidePlane);

/*
 *      Set some pointers to the burgers vector and glide plane lists
 *      that were created during initialization.
 *
 *      The burgers vectors in the pre-computed list are already normalized.
 */
        burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;
        numBurg = home->burgData.numBurgVectors;
        numGlissilePlanesPerBurg = home->burgData.numGlissilePlanesPerBurg;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;

/*
 *      Find the index in the reference burgers vector list of the
 *      current burgers vector, and find out how many glide planes
 *      are associated with the burgers vector.
 */
        burgIndex = GetBurgIndexBCC(burgVec, numBurg, burgList);

/*
 *      Select one of the burger's glide planes at random.
 */
        randVal = randm(&seed);

        numGlidePlanes  = numGlissilePlanesPerBurg[burgIndex];

	planeIndex = MIN(((int)(randVal*numGlidePlanes)), numGlidePlanes-1);
	planeIndex += burgFirstPlaneIndex[burgIndex];
	    
	VECTOR_COPY(glidePlane, planeList[planeIndex]);

 
        return;
}

#if 0
/*-------------------------------------------------------------------------
 *
 *      Function:     PickBCCScrewGlidePlaneOld
 *      Description:  Selects, at random, an appropriate glide plane for
 *                    BCC screw dislocations based on the burgers vector.
 *
 *      Arguments:
 *          burgVec      Input array containing the 3 components of the
 *                       burgers vector. This burgers vector should 
 *                       already have been rotated from the laboratory
 *                       frame to crystal frame by the caller.
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.  Returned glide plane is in the 
 *                       crystalographic frame.
 *
 *------------------------------------------------------------------------*/
static void PickBCCScrewGlidePlaneOld(real8 burgVec[3], real8 glidePlane[3])
{
        real8        randVal;
        static int   seed = 8917346;
#ifdef _OPENMP
#pragma omp threadprivate (seed)
#endif

        randVal = randm(&seed);

        if (randVal < 0.3333) {
            glidePlane[0] = 0.0;
            glidePlane[1] = 1.0 / sqrt(2.0);
            glidePlane[2] = -Sign(burgVec[1] * burgVec[2]) * glidePlane[1];
        } else if ((randVal >= 0.3333) && (randVal < 0.6666)) {
            glidePlane[0] = 1.0 / sqrt(2.0);
            glidePlane[1] = 0.0;
            glidePlane[2] = -Sign(burgVec[0] * burgVec[2]) * glidePlane[0];
        } else {
            glidePlane[0] = 1.0 / sqrt(2.0);
            glidePlane[1] = -Sign(burgVec[0] * burgVec[1]) * glidePlane[0];
            glidePlane[2] = 0.0;
        }

        return;
}
#endif


/*-------------------------------------------------------------------------
 *
 *      Function:     PickFCCScrewGlidePlane
 *      Description:  Selects, at random, an appropriate glide plane for
 *                    FCC screw dislocations based on the burgers vector.
 *
 *      Arguments:
 *          burgVec      Input array containing the 3 components of the
 *                       burgers vector. This burgers vector should 
 *                       already have been rotated from the laboratory
 *                       frame to crystal frame by the caller.
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.  Returned glide plane is in the 
 *                       crystalographic frame.
 *
 *------------------------------------------------------------------------*/
static void PickFCCScrewGlidePlane(real8 burgVec[3], real8 glidePlane[3])
{
        real8        randVal;
        static int   seed = 8917346;
#ifdef _OPENMP
#pragma omp threadprivate (seed)
#endif

        randVal = randm(&seed);

        if (burgVec[0] == 0.0e0) {
            glidePlane[0] =  1 / sqrt(3.0);
            glidePlane[1] =  glidePlane[0];
            glidePlane[2] = -Sign(burgVec[1] * burgVec[2]) * glidePlane[1];
            if (randVal < 0.5){
                glidePlane[0] = -glidePlane[1];
            }
        } else if (burgVec[1] == 0.0e0) {
            glidePlane[0] =  1 / sqrt(3.0);
            glidePlane[1] =  glidePlane[0];
            glidePlane[2] = -Sign(burgVec[0] * burgVec[2]) * glidePlane[0];
            if (randVal < 0.5){
                glidePlane[1] = -glidePlane[0];
            }
        } else {
            glidePlane[0] =  1 / sqrt(3.0);
            glidePlane[1] = -Sign(burgVec[0] * burgVec[1]) * glidePlane[0];
            glidePlane[2] =  glidePlane[0];
            if (randVal < 0.5){
                glidePlane[2] = -glidePlane[0];
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     PickHCPScrewGlidePlane
 *      Description:  Selects, at random, an appropriate glide plane for
 *                    HCP screw dislocations based on the burgers vector.
 *
 *      Arguments:
 *          burgVec      Input array containing the 3 components of the
 *                       burgers vector. This burgers vector should 
 *                       already have been rotated from the laboratory
 *                       frame to crystal frame by the caller.
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.  Returned glide plane is in the 
 *                       crystalographic frame.
 *
 *------------------------------------------------------------------------*/
static void PickHCPScrewGlidePlane(Home_t *home, real8 burgVec[3],
                                   real8 glidePlane[3])
{
        int          numBurg, burgIndex, planeIndex, numGlidePlanes;
        int          *numGlissilePlanesPerBurg, *burgFirstPlaneIndex;
        real8        randVal;
        real8        (*burgList)[3], (*planeList)[3];
        static int   seed = 8917346;
#ifdef _OPENMP
#pragma omp threadprivate (seed)
#endif


        VECTOR_ZERO(glidePlane);

/*
 *      Set some pointers to the burgers vector and glide plane lists
 *      that were created during initialization.
 *
 *      The burgers vectors in the pre-computed list are already normalized.
 */
        burgList = home->burgData.burgList;
        planeList = home->burgData.planeList;
        numBurg = home->burgData.numBurgVectors;
        numGlissilePlanesPerBurg = home->burgData.numGlissilePlanesPerBurg;
        burgFirstPlaneIndex = home->burgData.burgFirstPlaneIndex;

/*
 *      Find the index in the reference burgers vector list of the
 *      current burgers vector, and find out how many glide planes
 *      are associated with the burgers vector.
 */
        GetBurgIndexHCP(burgVec, 1, numBurg, burgList, &burgIndex);

/*
 *      Select one of the burger's glide planes at random
 *      Some burgers vectors have one glide plane only.
 */
        randVal = randm(&seed);

        numGlidePlanes  = numGlissilePlanesPerBurg[burgIndex];

	planeIndex = MIN(((int)(randVal*numGlidePlanes)), numGlidePlanes-1);
	planeIndex += burgFirstPlaneIndex[burgIndex];
	    
	VECTOR_COPY(glidePlane, planeList[planeIndex]);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     PickRhombohedralVaScrewGlidePlane
 *      Description:  Selects, at random, an appropriate glide plane for
 *                    screw dislocations in rhombohedral Vanadium based
 *                    on the burgers vector.
 *
 *      Arguments:
 *          burgVec      Input array containing the 3 components of the
 *                       burgers vector. This burgers vector should 
 *                       already have been rotated from the laboratory
 *                       frame to crystal frame by the caller.
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.  Returned glide plane is in the 
 *                       crystalographic frame.
 *
 *      NOTE: The burgers vectors and glide planes used in the function
 *            below are normalized values.  The non-normalized values
 *            are as follows:
 *
 *            glide burgers vectors
 *
 *             1.45000000000e+00,  1.45000000000e+00,  1.45000000000e+00
 *             1.31850000000e+00,  1.31850000000e+00, -1.18700000000e+00
 *             1.31850000000e+00, -1.18700000000e+00,  1.31850000000e+00
 *            -1.18700000000e+00,  1.31850000000e+00,  1.31850000000e+00
 *
 *            glide planes
 *
 *            -3.63297500000e+00,  3.63297500000e+00,                  0
 *             3.63297500000e+00,                  0, -3.63297500000e+00
 *                             0, -3.63297500000e+00,  3.63297500000e+00
 *            -3.63297500000e+00,  3.63297500000e+00,                  0
 *             3.29473250000e-01, -3.30350175000e+00, -3.30350175000e+00
 *            -3.30350175000e+00,  3.29473250000e-01, -3.30350175000e+00
 *             3.63297500000e+00,                  0, -3.63297500000e+00
 *             3.29473250000e-01, -3.30350175000e+00, -3.30350175000e+00
 *            -3.30350175000e+00, -3.30350175000e+00,  3.29473250000e-01
 *                             0, -3.63297500000e+00,  3.63297500000e+00
 *             3.29473250000e-01, -3.30350175000e+00, -3.30350175000e+00
 *            -3.30350175000e+00, -3.30350175000e+00,  3.29473250000e-01
 *
 *------------------------------------------------------------------------*/
static void PickRhomboVaScrewGlidePlane(real8 burgVec[3], real8 glidePlane[3])
{
        int          i, j, indexOfMax;
        real8        randVal, maxVal;
        real8        burgRefXburg[4];
        static int   seed = 8917346;
#ifdef _OPENMP
#pragma omp threadprivate (seed)
#endif
        static real8 burgRef[4][3] = {
        /* normalized glide burgers vectors */
        { 5.773502691896258e-01, 5.773502691896258e-01, 5.773502691896258e-01},
        { 5.964992662313060e-01, 5.964992662313060e-01,-5.370076822271977e-01},
        { 5.964992662313060e-01,-5.370076822271977e-01, 5.964992662313060e-01},
        {-5.370076822271977e-01, 5.964992662313060e-01, 5.964992662313060e-01}};
        /* normalized glide planes */
        static real8 glidePlanes[12][3]={
            { -7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
            {  7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
            {  0.0000000000e+00, -7.0710678119e-01,  7.0710678119e-01},

            { -7.0710678119e-01,  7.0710678119e-01,  0.0000000000e+00},
            {  7.0348253195e-02, -7.0535491891e-01, -7.0535491891e-01},
            {  7.0535491891e-01, -7.0348253195e-02,  7.0535491891e-01},

            {  7.0710678119e-01,  0.0000000000e+00, -7.0710678119e-01},
            {  7.0348253195e-02, -7.0535491891e-01, -7.0535491891e-01},
            { -7.0535491891e-01, -7.0535491891e-01,  7.0348253195e-02},

            {  0.0000000000e+00, -7.0710678119e-01,  7.0710678119e-01},
            {  7.0535491891e-01, -7.0348253195e-02,  7.0535491891e-01},
            { -7.0535491891e-01, -7.0535491891e-01,  7.0348253195e-02}};

        MatrixMult((real8 *)burgRef, 4, 3, 3,
                   burgVec, 1, 1,
                   burgRefXburg, 1);
            
        FindAbsMax(burgRefXburg, 4, &maxVal, &indexOfMax);

        randVal = randm(&seed);

        if (randVal < 0.3333) {
            j=0;
        } else if ((randVal >= 0.3333) && (randVal < 0.6666)) {
            j=1;
        } else {
            j=2;
        }

        for (i = 0; i < 3; i++) {
            glidePlane[i] = glidePlanes[3*indexOfMax+j][i];
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     PickScrewGlidePlane
 *      Description:  This is a generic dispatch function which calls an
 *                    appropriate material-specific routine to select, at
 *                    random, an appropriate glide plane for screw
 *                    dislocations in the specified type of material
 *                    crystal structure based on the burgers vector.
 *
 *      Arguments:
 *          burgVecIn    Input array containing the 3 components of the
 *                       burgers vector
 *          glidePlane   Output array in which the 3 components of the
 *                       selected glide plane will be returned to the
 *                       caller.
 *
 *------------------------------------------------------------------------*/
void PickScrewGlidePlane(Home_t *home, real8 burgVecIn[3], real8 glidePlane[3])
{
        real8   *burgVec = burgVecIn;
        real8    burgVecRotated[3];
        Param_t *param = home->param;

/*
 *      If necessary, rotate the burgers vector from the laboratory frame
 *      to the crystal frame.
 */
        if (param->useLabFrame) {
            Matrix33Vector3Multiply(param->rotMatrixInverse, burgVecIn,
                                    burgVecRotated);
            burgVec = burgVecRotated;
        }

        switch(param->materialType) {
            case MAT_TYPE_BCC:
               PickBCCScrewGlidePlane(home, burgVec, glidePlane);
                break;

            case MAT_TYPE_FCC:
                PickFCCScrewGlidePlane(burgVec, glidePlane);
                break;

            case MAT_TYPE_HCP:
                PickHCPScrewGlidePlane(home, burgVec, glidePlane);
                break;

            case MAT_TYPE_RHOMBOHEDRAL_VA:
                PickRhomboVaScrewGlidePlane(burgVec, glidePlane);
                break;
        }

/*
 *      If necessary, rotate the plane normal from the crystal frame
 *      back to the laboratory frame.
 */
        if (param->useLabFrame) {
            real8 plane[3];
            Matrix33Vector3Multiply(param->rotMatrix, glidePlane, plane);
            glidePlane[0] = plane[0];
            glidePlane[1] = plane[1];
            glidePlane[2] = plane[2];
        }

        return;
}
