/*****************************************************************************
 *
 *      Module:         WriteVisitBurgID.c
 *      Description:    This module contains a generic dispatch function
 *                      and material-specific functions to associate a
 *                      burgers vector 'ID' with a specific segment which
 *                      can be included in the VisIt 'segment' file.
 *
 *      Includes public functions:
 *          GetBurgID()
 *
 *      Includes private functions:
 *          GetBurgID_BCC()
 *          GetBurgID_FCC()
 *          GetBurgID_HCP()
 *          GetBurgID_RhombohedralV()
 *
 *
 ****************************************************************************/

#include "mpi_portability.h"

#include "Home.h"


/*---------------------------------------------------------------------------
 *
 *      Function:    GetBurgID_BCC
 *      Description: Return the type of burgers vector (identified by an
 *                   integer value) of the specified segment for the given
 *                   node.
 *
 *      Arguments:
 *          node   Pointer to node
 *          segID  Index of segment
 *
 *      Returns: Integer value indicating the type of burgers vector
 *               as defined below:
 *
 *                     value       burgers vector types
 *                       0         [ 1  1  1]
 *                       1         [-1  1  1]
 *                       2         [ 1 -1  1]
 *                       3         [ 1  1 -1]
 *                       4         [ 1  0  0]
 *                       5         [ 0  1  0]
 *                       6         [ 0  0  1]
 *                       7         [ 1  1  0] [ 1 -1  0]
 *                                 [ 1  0  1] [ 1  0 -1]
 *                                 [ 0  1  1] [ 0  1 -1]
 *                       8         All other types
 *
 ****************************************************************************/
static int GetBurgID_BCC(Home_t *home, Node_t *node, int segID)
{
        int          burgIndex, burgID;
        real8        bCryst[3], crossProd[3];
        Param_t     *param = home->param;
        static int   numBurgs = 13;
/*
 *      Define the IDs for each type (or grouping) of burgers vector
 *      we're interested in (plus a value with which to associate any
 *      burgers vector for which we don't explicitly look.
 */
        static int   burg2IDMap[14] =
                { 0,                /* [ 1  1  1] */
                  1,                /* [-1  1  1] */
                  2,                /* [ 1 -1  1] */
                  3,                /* [ 1  1 -1] */
                  4,                /* [ 1  0  0] */
                  5,                /* [ 0  1  0] */
                  6,                /* [ 0  0  1] */
                  7, 7, 7, 7, 7, 7, /* all 6 [1 1 0] types */
                  8 };              /* all other burgers vectors */
/*
 *      Define the list of burgers vectors for which we are explicitly
 *      looking.
 */
        static real8 burgList[14][3] = {
                /* glide burgers vectors */
                { 1.0,  1.0,  1.0},
                {-1.0,  1.0,  1.0},
                { 1.0, -1.0,  1.0},
                { 1.0,  1.0, -1.0},
                /* [100] type junction burgers vectors */
                { 1.0,  0.0,  0.0},
                { 0.0,  1.0,  0.0},
                { 0.0,  0.0,  1.0},
                /* [110] type junction burgers vectors */
                { 1.0,  1.0,  0.0},
                { 1.0, -1.0,  0.0},
                { 1.0,  0.0,  1.0},
                { 1.0,  0.0, -1.0},
                { 0.0,  1.0,  1.0},
                { 0.0,  1.0, -1.0} };
                

/*
 *      If needed, rotate the provided burgers vector from the lab frame
 *      into the crystal frame.
 */
        if (param->useLabFrame) {
            real8 bLab[3];
            bLab[X] = node->burgX[segID];
            bLab[Y] = node->burgY[segID];
            bLab[Z] = node->burgZ[segID];
            Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[X] = node->burgX[segID];
            bCryst[Y] = node->burgY[segID];
            bCryst[Z] = node->burgZ[segID];
        }

/*
 *      Look through the list of reference burgers vectors for one
 *      that matches the segment burgers vector.  If we find a match
 *      return the ID associated with that burgers vector.  If no
 *      match is found, the ID defaults to the final one in the list
 *      to indicate it was some *other* burgers vector type.
 */
        burgID = burg2IDMap[numBurgs];

        for (burgIndex = 0; burgIndex < numBurgs; burgIndex++) {
            cross(bCryst, burgList[burgIndex], crossProd);
            if (fabs(DotProduct(crossProd, crossProd)) < 1.0e-3) {
                burgID = burg2IDMap[burgIndex];
                break;
            }
        }

        return(burgID);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    GetBurgID_FCC
 *      Description: Return the type of burgers vector (identified by an
 *                   integer value) of the specified segment for the given
 *                   node.
 *
 *      Arguments:
 *          node   Pointer to node
 *          segID  Index of segment
 *
 *      Returns: Integer value indicating the type of burgers vector
 *               as defined below:
 *
 *                     value       burgers vector types
 *                       0         [ 1  1  0]
 *                       1         [ 1 -1  0]
 *                       2         [ 1  0  1]
 *                       3         [ 1  0 -1]
 *                       4         [ 0  1  1]
 *                       5         [ 0  1 -1]
 *                       6         [ 1  0  0]
 *                       7         [ 0  1  0]
 *                       8         [ 0  0  1]
 *                       9         All [ 2  1  1] types
 *                       10        All other types
 *
 ****************************************************************************/
static int GetBurgID_FCC(Home_t *home, Node_t *node, int segID)
{
        int          burgIndex, burgID;
        real8        bCryst[3], crossProd[3];
        Param_t     *param = home->param;
        static int   numBurgs = 21;
/*
 *      Define the IDs for each type (or grouping) of burgers vector
 *      we're interested in (plus a value with which to associate any
 *      burgers vector for which we don't explicitly look.
 */
        static int   burg2IDMap[22] =
                { 0,                /* [ 1  1  0] */
                  1,                /* [ 1 -1  0] */
                  2,                /* [ 1  0  1] */
                  3,                /* [ 1  0 -1] */
                  4,                /* [ 0  1  1] */
                  5,                /* [ 0  1 -1] */
                  6,                /* [ 1  0  0] */
                  7,                /* [ 0  1  0] */
                  8,                /* [ 0  0  1] */
                  9, 9, 9, 9, 9, 9, /* all 12 [2 1 1] types */
                  9, 9, 9, 9, 9, 9, 
                  10 };             /* all other burgers vectors */
/*
 *      Define the list of burgers vectors for which we are explicitly
 *      looking.
 */
        static real8 burgList[21][3] = {
                /* glide burgers vectors */
                { 1.0,  1.0,  0.0},
                { 1.0, -1.0,  0.0},
                { 1.0,  0.0,  1.0},
                { 1.0,  0.0, -1.0},
                { 0.0,  1.0,  1.0},
                { 0.0,  1.0, -1.0},
                /* junction burgers vectors */
                { 1.0,  0.0,  0.0},
                { 0.0,  1.0,  0.0},
                { 0.0,  0.0,  1.0},
                { 2.0,  1.0,  1.0},
                {-2.0,  1.0,  1.0},
                { 2.0, -1.0,  1.0},
                { 2.0,  1.0, -1.0},
                { 1.0,  2.0,  1.0},
                {-1.0,  2.0,  1.0},
                { 1.0, -2.0,  1.0},
                { 1.0,  2.0, -1.0},
                { 1.0,  1.0,  2.0},
                {-1.0,  1.0,  2.0},
                { 1.0, -1.0,  2.0},
                { 1.0,  1.0, -2.0} };

/*
 *      If needed, rotate the provided burgers vector from the lab frame
 *      into the crystal frame.
 */
        if (param->useLabFrame) {
            real8 bLab[3];
            bLab[X] = node->burgX[segID];
            bLab[Y] = node->burgY[segID];
            bLab[Z] = node->burgZ[segID];
            Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[X] = node->burgX[segID];
            bCryst[Y] = node->burgY[segID];
            bCryst[Z] = node->burgZ[segID];
        }

/*
 *      Look through the list of reference burgers vectors for one
 *      that matches the segment burgers vector.  If we find a match
 *      return the ID associated with that burgers vector.  If no
 *      match is found, the ID defaults to the final one in the list
 *      to indicate it was some *other* burgers vector type.
 */
        burgID = burg2IDMap[numBurgs];

        for (burgIndex = 0; burgIndex < numBurgs; burgIndex++) {
            cross(bCryst, burgList[burgIndex], crossProd);
            if (fabs(DotProduct(crossProd, crossProd)) < 1.0e-3) {
                burgID = burg2IDMap[burgIndex];
                break;
            }
        }

        return(burgID);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    GetBurgID_HCP
 *      Description: Return the type of burgers vector (identified by an
 *                   integer value) of the specified segment for the given
 *                   node.
 *
 *      Arguments:
 *          node   Pointer to node
 *          segID  Index of segment
 *
 *      Returns: Integer value indicating the type of burgers vector.
 *               ID values 0 thru 9 correspond to the first 10 glide
 *               burgers vectors in the pre-computed HCP burgers vector
 *               list.  A value of 10 indicates the burgers vector was
 *               some other type.
 *
 ****************************************************************************/
static int GetBurgID_HCP(Home_t *home, Node_t *node, int segID)
{
        int   burgIndex;
        real8 bCryst[3];
        real8 (*burgList)[3];
        Param_t *param = home->param;
        static int numBurgs = 10;

/*
 *      Use the reference burgers vector list that was already pre-computed.
 *      We're only interested in the first 10 glide-burgers vectors for
 *      this routine.  All other types of burgers vector are grouped
 *      together with a single ID value.
 */
        burgList = home->burgData.burgList;

/*
 *      If needed, rotate the provided burgers vector from the lab frame
 *      into the crystal frame.
 */
        if (param->useLabFrame) {
            real8 bLab[3];
            bLab[X] = node->burgX[segID];
            bLab[Y] = node->burgY[segID];
            bLab[Z] = node->burgZ[segID];
            Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[X] = node->burgX[segID];
            bCryst[Y] = node->burgY[segID];
            bCryst[Z] = node->burgZ[segID];
        }

/*
 *      Look through the list of reference burgers vectors for one
 *      that matches the segment burgers vector.  If we find a match
 *      return the ID associated with that burgers vector.  If no
 *      match is found, the ID defaults to the final one in the list
 *      to indicate it was some *other* burgers vector type.
 */
        GetBurgIndexHCP(bCryst, 0, numBurgs, burgList, &burgIndex);

        if (burgIndex < 0) {
            burgIndex = numBurgs;
        }

        return(burgIndex);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    GetBurgID_RhombohedralV
 *      Description: Return the type of burgers vector (identified by an
 *                   integer value) of the specified segment for the given
 *                   node.
 *
 *      Arguments:
 *          node   Pointer to node
 *          segID  Index of segment
 *
 *      Returns: Integer value indicating the type of burgers vector
 *               as defined below:
 *
 *                     value       burgers vector types
 *                       0         [ 1.4500  1.4500  1.4500]
 *                       1         [-1.1870  1.3185  1.3185]
 *                       2         [ 1.3185 -1.1870  1.3185]
 *                       3         [ 1.3185  1.3185 -1.1870]
 *                       4         [ 1.315e-01 1.315e-01 2.637e+00]
 *                       5         [ 1.315e-01 2.637e+00 1.315e-01]
 *                       6         [ 2.637e+00 1.315e-01 1.315e-01]
 *                       7         3 [0 2.5055 -2.5055] types
 *                                 3 [2.7685 2.7685 0.263] types
 *                       8         All other types
 *
 ****************************************************************************/
static int GetBurgID_RhombohedralV(Home_t *home, Node_t *node, int segID)
{
        int          burgIndex, burgID;
        real8        bCryst[3], crossProd[3];
        Param_t     *param = home->param;
        static int   numBurgs = 13;
/*
 *      Define the IDs for each type (or grouping) of burgers vector
 *      we're interested in (plus a value with which to associate any
 *      burgers vector for which we don't explicitly look.
 */
        static int   burg2IDMap[14] =
                { 0,                /* [ 1.4500  1.4500  1.4500] */
                  1,                /* [-1.1870  1.3185  1.3185] */
                  2,                /* [ 1.3185 -1.1870  1.3185] */
                  3,                /* [ 1.3185  1.3185 -1.1870] */
                  4,                /* [ 1.315e-01 1.315e-01 2.637e+00] */ 
                  5,                /* [ 1.315e-01 2.637e+00 1.315e-01] */
                  6,                /* [ 2.637e+00 1.315e-01 1.315e-01] */
                  7, 7, 7,          /* 3 [0 2.5055 -2.5055] types */
                  7, 7, 7,          /* 3 [2.7685 2.7685 0.263] types */
                  8 };              /* all other burgers vectors */
/*
 *      Define the list of burgers vectors for which we are explicitly
 *      looking.
 */
        static real8 burgList[13][3] = {
                /* glide burgers vectors */
                { 1.4500,  1.4500,  1.4500},
                {-1.1870,  1.3185,  1.3185},
                { 1.3185, -1.1870,  1.3185},
                { 1.3185,  1.3185, -1.1870},
                /* junction burgers vectors */
                { 1.3150e-01,  1.3150e-01,  2.6370e+00},
                { 1.3150e-01,  2.6370e+00,  1.3150e-01},
                { 2.6370e+00,  1.3150e-01,  1.3150e-01},
                { 0.0000e+00,  2.5055e+00, -2.5055e+00},
                { 2.5055e+00,  0.0000e+00, -2.5055e+00},
                { 2.5055e+00, -2.5055e+00,  0.0000e+00},
                { 2.7685e+00,  2.7685e+00,  2.6300e-01},
                { 2.7685e+00,  2.6300e-01,  2.7685e+00},
                { 2.6300e-01,  2.7685e+00,  2.7685e+00} };

/*
 *      If needed, rotate the provided burgers vector from the lab frame
 *      into the crystal frame.
 */
        if (param->useLabFrame) {
            real8 bLab[3];
            bLab[X] = node->burgX[segID];
            bLab[Y] = node->burgY[segID];
            bLab[Z] = node->burgZ[segID];
            Matrix33Vector3Multiply(param->rotMatrixInverse, bLab, bCryst);
        } else {
            bCryst[X] = node->burgX[segID];
            bCryst[Y] = node->burgY[segID];
            bCryst[Z] = node->burgZ[segID];
        }

/*
 *      Look through the list of reference burgers vectors for one
 *      that matches the segment burgers vector.  If we find a match
 *      return the ID associated with that burgers vector.  If no
 *      match is found, the ID defaults to the final one in the list
 *      to indicate it was some *other* burgers vector type.
 */
        burgID = burg2IDMap[numBurgs];

        for (burgIndex = 0; burgIndex < numBurgs; burgIndex++) {
            cross(bCryst, burgList[burgIndex], crossProd);
            if (fabs(DotProduct(crossProd, crossProd)) < 1.0e-3) {
                burgID = burg2IDMap[burgIndex];
                break;
            }
        }

        return(burgID);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    GetBurgID
 *      Description: Generic dispatch function to invoke a material-specific
 *                   function to return the type of burgers vector 
 *                   (identified by an integer value) of the specified
 *                   segment for the given node.
 *
 *      Arguments:
 *          node   Pointer to node
 *          segID  Index of segment
 *
 *      Returns:  Integer value indicating the type of burgers vector.
 *                See material specific functions for mapping of values
 *                to burgers vectors.
 *          
 ****************************************************************************/
int GetBurgID(Home_t *home, Node_t *node, int segID)
{
        int burgID;

        switch (home->param->materialType) {

            case MAT_TYPE_BCC:
                burgID = GetBurgID_BCC(home, node, segID);
                break;

            case MAT_TYPE_FCC:
                burgID = GetBurgID_FCC(home, node, segID);
                break;

            case MAT_TYPE_HCP:
                burgID = GetBurgID_HCP(home, node, segID);
                break;

            case MAT_TYPE_RHOMBOHEDRAL_VA:
                burgID = GetBurgID_RhombohedralV(home, node, segID);
                break;

            default:
                burgID = 0;
                break;
        }

        return(burgID);
}
