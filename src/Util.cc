/***************************************************************************
 *
 *      Module:       Util.c
 *      Description:  Contains various and sundry small utility routines
 *
 *      Includes:
 *
 *          3D to linear index mapping functions:
 *
 *              DecodeCellIdx()
 *              DecodeCell2Idx()
 *              DecodeDomainIdx()
 *              EncodeCellIdx()
 *              EncodeCell2Idx()
 *              EncodeDomainIdx()
 *
 *          Numerical operations:
 *
 *              cross()
 *              CSpline()
 *              CSplint()
 *              DecompVec()
 *              FindAbsMax()
 *              FindAbsMin()
 *              FindMax()
 *              FindMin()
 *              GetPlaneNormalFromPoints()
 *              GetUnitVector()
 *              InterpolateL()
 *              InterpolateBL()
 *              Normal()
 *              Normalize()
 *              NormalizeVec()
 *              NormalizedCrossVector()
 *              Orthogonalize()
 *              xvector()
 *
 *          Nodal operations:
 *
 *              AllocNodeArms()
 *              FreeNode()
 *              FreeNodeArms()
 *              GetNeighborNode()
 *              GetNodeFromIndex()
 *              GetNodeFromTag()
 *              InitNodeArm()        static func
 *              NodeHasSessileBurg()
 *              PrintNode()
 *              ReallocNodeArms()
 *
 *          Recycled node handling functions:
 *
 *              GetRecycledNodeTag()
 *              GetFreeNodeTag()
 *              RecycleNodeTag()
 *
 *          PBC image reconciliation functions:
 *
 *              FoldBox()
 *              PBCPOSITION()
 *              ZImage()
 *
 *          Topology change support functions:
 *
 *              ChangeArmBurg()
 *              ChangeConnection()
 *              CompressArmLists()
 *              DomainOwnsSeg()
 *              InsertArm()
 *              MarkNodeForceObsolete()
 *              NodeDistance()       static func
 *              RepositionNode()
 *              ResetGlidePlane()
 *              ResetSegForces2()
 *              ResetSegForces()
 *              SubtractSegForce()
 *
 *          Node sorting and ordering functions:
 *
 *              CollisonNodeOrder()
 *              IntCompare()
 *              OrderNodes()
 *              OrderTags()
 *
 *          Miscellaneous functions:
 *
 *              ApplyConstraintsToVelocity()
 *              Connected()
 *              Fatal()
 *              Warning()
 *              GetArmID()
 *              GetBurgersVectorNormal()
 *              GetGlideConstraintList()
 *              Getline()
 *              FindCellCenter()
 *              LocateCell()
 *              NodePinned()
 *              randm()
 *              ReadTabulatedData()
 *              Sign()
 *              StrEquiv();
 *              Uniq()
 *              Print3()
 *              Print2x3()
 *              Print3x3()
 *
 **************************************************************************/
#include <stdio.h>
#include <memory.h>
#include <ctype.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Abort.h"
#include "InData.h"
#include "QueueOps.h"

/*-------------------------------------------------------------------------
 *
 *      Function:     DecodeCellIdx
 *      Description:  Given the ID of a cell, return the indices of
 *                    the cell in the X, Y and Z dimensions.
 *
 *                    IMPORTANT!  The base set of cells is padded
 *                    by 1 cell on each side of each dimension to
 *                    allow for periodic boundaries.   The indices
 *                    returned to the caller will have been
 *                    adjusted to take this layer of ghost cells
 *                    into account
 *
 *                    NOTE: It is assumed in the cell index encoding/
 *                    decoding functions that cells are assigned in
 *                    the 3D space varying first in Z and last in the
 *                    X dimension.
 *
 *      Arguments:
 *          cellID    Cell Id as returned by EncodeCellIdx()
 *          xIndex    Pointer to locaton in which to return to the
 *                    caller the index of the cell in the X dimension
 *          yIndex    Pointer to locaton in which to return to the
 *                    caller the index of the cell in the Y dimension
 *          zIndex    Pointer to locaton in which to return to the
 *                    caller the index of the cell in the Z dimension
 *
 *------------------------------------------------------------------------*/
void DecodeCellIdx(Home_t *home, int cellID, int *xIndex,
                   int *yIndex, int *zIndex)
{
        int rem, fullPlane, fullLine;

        fullLine  = home->param->nZcells + 2;
        fullPlane = (home->param->nZcells + 2) * (home->param->nYcells + 2);

        *xIndex = cellID / fullPlane;
        rem = cellID - (*xIndex) * fullPlane;
        *yIndex = rem / fullLine;
        *zIndex = rem - (*yIndex) * fullLine;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     DecodeCell2Idx
 *      Description:  Given a cell2 id, calculate the indices of the
 *                    cell2 in the X, Y and Z dimensions.
 *
 *                    NOTE: It is assumed in the cell2 index encoding/
 *                    decoding functions that cell2s are assigned in
 *                    the 3D space varying first in Z and last in the
 *                    X dimension.
 *
 *      Arguments:
 *          cell2ID   Cell2 identifier as returned by EncodeCell2Idx().
 *          xIndex    Pointer to location in which to return to the
 *                    caller the cell2 index in the X dimension
 *          yIndex    Pointer to location in which to return to the
 *                    caller the cell2 index in the Y dimension
 *          zIndex    Pointer to location in which to return to the
 *                    caller the cell2 index in the Z dimension
 *
 *------------------------------------------------------------------------*/
void DecodeCell2Idx(Home_t *home, int cell2ID, int *xIndex,
                    int *yIndex, int *zIndex)
{
        int rem, fullPlane, fullLine;

        fullLine  = home->cell2nz;
        fullPlane = home->cell2nz * home->cell2ny;

        *xIndex = cell2ID / fullPlane;
        rem = cell2ID - (*xIndex) * fullPlane;
        *yIndex = rem / fullLine;
        *zIndex = rem - (*yIndex) * fullLine;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     DecodeDomainIdx
 *      Description:  Given a domain ID, calculate the indices of the
 *                    domain in the X, Y and Z dimensions.
 *
 *                    NOTE: It is assumed in the domain index encoding/
 *                    decoding functions that domains are assigned in
 *                    the 3D space varying first in Z and last in the
 *                    X dimension.
 *
 *      Arguments:
 *          domID     domain identifier as calculated by EncodeDomainIdx()
 *          xIndex    Pointer to location in which to return the
 *                    domain's index in the X dimension.
 *          yIndex    Pointer to location in which to return the
 *                    domain's index in the Y dimension.
 *          zIndex    Pointer to location in which to return the
 *                    domain's index in the Z dimension.
 *
 *------------------------------------------------------------------------*/
void DecodeDomainIdx(Home_t *home, int domID, int *xIndex,
                     int *yIndex, int *zIndex)
{
        int rem, fullPlane, fullLine;

        fullLine  = home->param->nZdoms;
        fullPlane = home->param->nZdoms * home->param->nYdoms;

        *xIndex = domID / fullPlane;
        rem = domID - (*xIndex) * fullPlane;
        *yIndex = rem / fullLine;
        *zIndex = rem - (*yIndex) * fullLine;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     EncodeCellIdx
 *      Description:  Given the indices of a cell in each of the X, Y
 *                    and Z dimensions, return its ID encoded as
 *                    a single number.
 *
 *                    IMPORTANT!  The base set of cells is padded
 *                    by 1 cell on each side of each dimension to
 *                    allow for periodic boundaries.  The encoded
 *                    index returned to the caller will have been
 *                    adjusted to take this layer of ghost cells
 *                    into account
 *
 *                    NOTE: It is assumed in the cell index encoding/
 *                    decoding functions that cells are assigned in
 *                    the 3D space varying first in Z and last in the
 *                    X dimension.
 *
 *      Arguments:
 *          xIndex    Base index of the cell in the X dimension
 *          yIndex    Base index of the cell in the Y dimension
 *          zIndex    Base index of the cell in the Z dimension
 *
 *      Returns:  Cell ID as single index into a 1-dimensional array.
 *
 *------------------------------------------------------------------------*/
int EncodeCellIdx(Home_t *home, int xIndex, int yIndex, int zIndex)
{
        int cellID;

        cellID = zIndex +
                 yIndex * (home->param->nZcells+2) +
                 xIndex * (home->param->nZcells+2) * (home->param->nYcells+2);

        return(cellID);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     EncodeCell2Idx
 *      Description:  Given the indices of a cell2 in each of the X, Y
 *                    and Z dimensions, return its index encoded as
 *                    a single number.
 *
 *                    NOTE: It is assumed in the cell2 index encoding/
 *                    decoding functions that cell2s are assigned in
 *                    the 3D space varying first in Z and last in the
 *                    X dimension.
 *
 *      Arguments:
 *          xIndex    Index of the cell2 in the X dimension
 *          yIndex    Index of the cell2 in the Y dimension
 *          zIndex    Index of the cell2 in the Z dimension
 *
 *------------------------------------------------------------------------*/
int EncodeCell2Idx(Home_t *home, int xIndex, int yIndex, int zIndex)
{
        int index;

        index = zIndex +
                yIndex * home->cell2nz +
                xIndex * home->cell2nz * home->cell2ny;

        return(index);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     EncodeDomainIdx
 *      Description:  Given the indices of a domain in each of the X, Y
 *                    and Z dimensions, return its domainID encoded as
 *                    a single number.
 *
 *                    NOTE: It is assumed in the domain index encoding/
 *                    decoding functions that domains are assigned in
 *                    the 3D space varying first in Z and last in the
 *                    X dimension.
 *
 *      Arguments:
 *          xIndex    Index of the domain in the X dimension
 *          yIndex    Index of the domain in the Y dimension
 *          zIndex    Index of the domain in the Z dimension
 *
 *------------------------------------------------------------------------*/
int EncodeDomainIdx(Home_t *home, int xIndex, int yIndex, int zIndex)
{
        int domID;

        domID = zIndex +
                yIndex * home->param->nZdoms +
                xIndex * home->param->nZdoms * home->param->nYdoms;

        return(domID);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     cross
 *      Description:  Calculates the cross product of two vector
 *                    supplied as arrays by the caller.
 *
 *      Arguments:
 *          a    3 element array containing components of the first vector
 *          b    3 element array containing components of the second vector
 *          c    3 element array in which to return to the caller the
 *               components of the cross product of vectors <a> and <b>.
 *
 *------------------------------------------------------------------------*/
void cross(real8 a[3], real8 b[3], real8 c[3])
{
    c[0]=a[1]*b[2]-a[2]*b[1];
    c[1]=a[2]*b[0]-a[0]*b[2];
    c[2]=a[0]*b[1]-a[1]*b[0];
}


/*------------------------------------------------------------------------
 *
 *      Function:    CSpline
 *      Description: Construct the cubic spline interpolant y2
 *                   for a function f(), defined at the numbers
 *                   x(0) < x(1) < ... < x(n), satisfying y2''(x(0)) = 0,
 *                   and y2''(x(n)) = 0:
 *
 *      Arguments
 *          x,y        Contains tabulated function, i.e.  y = f(x[n]) with
 *                     x[0] < x[1] < x[2] ... < x[numPoints-1]
 *          y2         array in which to return the 2nd derivatives
 *                     of the interpolating function at the tabulated
 *                     points x.
 *          numPoints  number of elements in arrays <x>, <y>, and <y2>
 *
 *-----------------------------------------------------------------------*/
void CSpline(real8 *x, real8 *y, real8 *y2, int numPoints)
{
        int   i, n;
        real8 p, qn, sig, un;
        real8 *u;

        u = (real8 *)calloc(1, sizeof(real8) * numPoints);
        n = numPoints-1;

/*
 *      Set derivatives to zero at lower boundary for "natural cubic spline"
 */
        y2[0] = 0.0;
        u[0]  = 0.0;

/*
 *      Decomposition loop of the tridiagonal algorithm.  y2 and u
 *      are used for temporary storage of the decomposed factors
 */
        for (i = 1; i < n; i++) {
            sig = (x[i] - x[i-1]) / (x[i+1] - x[i-1]);
            p = sig * y2[i-1] + 2.0;
            y2[i] = (sig - 1.0) / p;
            u[i] = (6.0 * ((y[i+1]-y[i]) / (x[i+1]-x[i]) - (y[i] - y[i-1]) /
                    (x[i]-x[i-1])) / (x[i+1]-x[i-1]) - (sig * u[i-1])) / p;
        }

/*
 *      Set derivatives to zero at upper boundary for "natural cubic spline"
 */
        qn = 0.0;
        un = 0.0;

        y2[n] = (un - (qn * u[n-1])) / ((qn * y2[n-1]) + 1.0);

/*
 *      Backsubstitution loop of the tridiagonal algorithm
 */
        for (i = n-1; i >= 0; i--) {
            y2[i] = y2[i] * y2[i+1] + u[i];
        }

        free(u);

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     CSplint
 *      Description:  Perform cubic spline interpolation to obtain
 *                    the value <y>
 *
 *      Arguments
 *          xa, ya     Contain tabulated function, i.e.  ya = f(xa[n]) with
 *                     xa[0] < xa[1] < xa[2] ... < xa[numPoints-1]
 *          y2         array containing the 2nd derivatives of the
 *                     interpolating function at the tabulated points x
 *                     as returned from the CSpline() function.
 *          numPoints  number of elements in arrays <x>, <y>, and <y2>
 *          x          Value for which to obtain the value of the
 *                     interpolated function.
 *          y          location n which to return the value of
 *                     the interpolated function at <x>.
 *
 *-----------------------------------------------------------------------*/
void CSplint(real8 *xa, real8 *ya, real8 *y2, int numPoints,
            real8 x, real8 *y)
{
        int   high, mid, low;
        real8 a, b, h;

/*
 *      Just a sanity check on the input value.
 */
        if ((x < xa[0]) || (x > xa[numPoints-1])) {
            Fatal("CSplint(): specified X outside interpolation range");
        }

/*
 *      Find the indices of the values in <xa> that bracket the
 *      specified <x> value.
 */
        low = 0;
        high = numPoints-1;

        while ((high - low) > 1) {
            mid = (high + low) >> 1;
            if (xa[mid] > x) {
                high = mid;
            } else {
                low = mid;
            }
        }

/*
 *      Evaluate the cubic spline polynomial
 */
        h = xa[high] - xa[low];
        a = (xa[high] - x) / h;
        b = (x - xa[low]) / h;

        *y = (a * ya[low]) + (b * ya[high]) +
             ((pow(a,3.0)-a) * y2[low] + (pow(b,3)-b) * y2[high]) *
             (h*h) / 6.0;

        return;
}


/***************************************************************************
 *
 *      Function:    DecompVec
 *
 *      Description: Decompose a vector (inVec) into the two directions given
 *                   by <vec1> and <vec2>.  This should work with any
 *                   coordinate system.
 *
 *      Arguments:
 *          IN:  inVec       Vector to be decomposed
 *          IN:  vec1, vec2  Direction vectors
 *          OUT: ansVec      Resulting decomposition
 *
 ***************************************************************************/
void DecompVec(real8 inVec[3], real8 vec1[3], real8 vec2[3], real8 ansVec[2])
{
        real8 cc, dec1, dec2;

        cc   = DotProduct(vec1,  vec2);
        dec1 = DotProduct(inVec, vec1);
        dec2 = DotProduct(inVec, vec2);

        ansVec[0] = (dec1 - (dec2 * cc)) / (1.0 - (cc * cc));
        ansVec[1] = (dec2 - (dec1 * cc)) / (1.0 - (cc * cc));

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FindAbsMax
 *      Description:  Given an array of real8's, find the element with the
 *                    highest absolute value, return its absolute value
 *                    and its position in the array.
 *
 *                    Note: This function assumes there is at least
 *                    one element in the array!
 *
 *      Arguments:
 *          array       Array to be searched
 *          numElements Number of elements in array to search
 *          maxVal      Location in which to return the absolute value
 *                      of the element with the highest absolute value
 *          indexOfMax  Index of the element with the highest absolute value
 *
 *------------------------------------------------------------------------*/
void FindAbsMax(real8 *array, int numElements, real8 *maxVal, int *indexOfMax)
{
        int i;

        *maxVal     = fabs(array[0]);
        *indexOfMax = 0;

        for (i = 1; i < numElements; i++) {
            if (fabs(array[i]) > *maxVal) {
                *maxVal = fabs(array[i]);
                *indexOfMax = i;
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FindMax
 *      Description:  Given an array of real8's, find the element with the
 *                    highest value, return its value and position in the
 *                    array.
 *
 *                    Note: This function assumes there is at least
 *                    one element in the array!
 *
 *      Arguments:
 *          array       Array to be searched
 *          numElements Number of elements in array to search
 *          maxVal      Location in which to return the value
 *                      of the element with the highest value
 *          indexOfMax  Index of the element with the highest value
 *
 *------------------------------------------------------------------------*/
void FindMax(real8 *array, int numElements, real8 *maxVal, int *indexOfMax)
{
        int i;

        *maxVal     = array[0];
        *indexOfMax = 0;

        for (i = 1; i < numElements; i++) {
            if (array[i] > *maxVal) {
                *maxVal = array[i];
                *indexOfMax = i;
            }
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     FindMaxExcluding
 *      Description:  Given an array of real8's, find the element with the
 *                    highest value, return its value and position in the
 *                    array.
 *
 *                    Note: This function assumes there is at least
 *                    one element in the array!
 *
 *      Arguments:
 *          array       Array to be searched
 *          numElements Number of elements in array to search
 *          maxVal      Location in which to return the value
 *                      of the element with the highest value
 *          indexOfMax  Index of the element with the highest value
 *
 *------------------------------------------------------------------------*/

void FindMaxExcluding(real8 *array, int numElements, real8 *maxVal, int *indexOfMax, int nonIndex1, int nonIndex2)
{
        *maxVal     = -1.0e+15;
        *indexOfMax = -1;

        for (int i=0; i<numElements; i++)
        {
            if( (i != nonIndex1) && (i != nonIndex2) )
            {
                if (array[i] > *maxVal)
                {
                    *maxVal     = array[i];
                    *indexOfMax = i;
                }
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FindAbsMin
 *      Description:  Given an array of real8's, find the element with the
 *                    lowest absolute value, return its absolute value
 *                    and its position in the array.
 *
 *                    Note: This function assumes there is at least
 *                    one element in the array!
 *
 *      Arguments:
 *          array       Array to be searched
 *          numElements Number of elements in array to search
 *          minVal      Location in which to return the absolute value
 *                      of the element with the lowest absolute value
 *          indexOfMin  Index of the element with the lowest absolute value
 *
 *------------------------------------------------------------------------*/
void FindAbsMin(real8 *array, int numElements, real8 *minVal, int *indexOfMin)
{
        int i;

        *minVal = fabs(array[0]);
        *indexOfMin = 0;

        for (i = 1; i < numElements; i++) {
            if (fabs(array[i]) < *minVal) {
                *minVal = fabs(array[i]);
                *indexOfMin = i;
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FindMin
 *      Description:  Given an array of real8's, find the element with the
 *                    lowest value, return its value and position in the array.
 *
 *                    Note: This function assumes there is at least
 *                    one element in the array!
 *
 *      Arguments:
 *          array       Array to be searched
 *          numElements Number of elements in array to search
 *          minVal      Location in which to return the value
 *                      of the element with the lowest value
 *          indexOfMin  Index of the element with the lowest value
 *
 *------------------------------------------------------------------------*/
void FindMin(real8 *array, int numElements, real8 *minVal, int *indexOfMin)
{
        int i;

        *minVal = array[0];
        *indexOfMin = 0;

        for (i = 1; i < numElements; i++) {
            if (array[i] < *minVal) {
                *minVal = array[i];
                *indexOfMin = i;
            }
        }

        return;
}



/*-------------------------------------------------------------------------
 *
 *      Function:    GetPlaneNormFromPoints
 *      Description: Calculate the normal vector for a plane containing
 *                   three specified points.
 *
 *      Arguments:
 *          IN:  p1, p2, p3 Three element arrays containing coordinates
 *                          (x, y, z) of the three points respectively.
 *
 *          OUT: normalVec  normal vector of the plane containing p1,
 *                          p2 and p3
 *
 *------------------------------------------------------------------------*/
void GetPlaneNormFromPoints(real8 p1[3], real8 p2[3], real8 p3[3],
                            real8 normalVec[3])
{
        real8 x1, x2, x3;
        real8 y1, y2, y3;
        real8 z1, z2, z3;
        real8 A[3][3];

        x1 = p1[0]; y1 = p1[1]; z1 = p1[2];
        x2 = p2[0]; y2 = p2[1]; z2 = p2[2];
        x3 = p3[0]; y3 = p3[1]; z3 = p3[2];

        A[0][0] = 1; A[0][1] = y1; A[0][2] = z1;
        A[1][0] = 1; A[1][1] = y2; A[1][2] = z2;
        A[2][0] = 1; A[2][1] = y3; A[2][2] = z3;

        normalVec[0] = Matrix33_Det(A);

        A[0][0] = x1; A[0][1] = 1; A[0][2] = z1;
        A[1][0] = x2; A[1][1] = 1; A[1][2] = z2;
        A[2][0] = x3; A[2][1] = 1; A[2][2] = z3;

        normalVec[1] = Matrix33_Det(A);

        A[0][0] = x1; A[0][1] = y1; A[0][2] = 1;
        A[1][0] = x2; A[1][1] = y2; A[1][2] = 1;
        A[2][0] = x3; A[2][1] = y3; A[2][2] = 1;

        normalVec[2] = Matrix33_Det(A);

        NormalizeVec(normalVec);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     GetUnitVector
 *      Description:  Given the coordinates of two points, find
 *                    the vector pointing from point 1 to point 2.
 *
 *      Arguments:
 *          unitFlag     Flag indicating if the returned vector should
 *                       be a unit vector.  1 if yes, 0 if no.
 *          x1, y1, z1   Coordinates of point 1
 *          x2, y2, z2   Coordinates of point 2
 *          ux, uy, uz   Pointers to components of (unit) vector
 *                       returned to caller.
 *          disMag       Magnitude of the vector from point 1 to
 *                       point 2.
 *
 *------------------------------------------------------------------------*/
void GetUnitVector(int unitFlag,
                   real8 x1, real8 y1, real8 z1,
                   real8 x2, real8 y2, real8 z2,
                   real8 *ux, real8 *uy, real8 *uz,
                   real8 *disMag)
{
        *ux = x2 - x1;
        *uy = y2 - y1;
        *uz = z2 - z1;

        *disMag = sqrt(*ux * *ux + *uy * *uy + *uz * *uz);

/*
 *      Assumption is that the points involved are node
 *      coordinates.  If the distance between the two points
 *      is too small, it means we have two nodes basically on
 *      top of each other which should not be.  However, if
 *      it does happen, just print a warning, treat it like
 *      a zero-length vector and let the caller deal with it.
 */
        if (*disMag < 1.0e-10) {
#if 0
           printf("GetUnitVector:  disMag < min:\n"
                  "    pointA = (%20.15f %20.15f %20.15f)\n"
                  "    pointB = (%20.15f %20.15f %20.15f)",
                  x1, y1, z1, x2, y2, z2);
#endif
           *ux = 0.0;
           *uy = 0.0;
           *uz = 0.0;

           *disMag = 0.0;

        } else {

           if (unitFlag == 1) {
              *ux = *ux / *disMag;
              *uy = *uy / *disMag;
              *uz = *uz / *disMag;
           }

        }

        return;
}


/*************************************************************************
 *
 *      Function:     InterpolateL
 *      Description:  Use linear interpolation to calculate the
 *                    value <y> for f(x) given a set of values for
 *                    the function f().
 *
 *      Arguments
 *          xa, ya     Contains tabulated function, i.e.  ya = f(xa[n]) with
 *                     xa[0] < xa[1] < xa[2] ... < xa[numPoints-1]
 *          numPoints  number of elements in arrays <x>, <y>, and <y2>
 *          x          Value for which to obtain the value of the
 *                     interpolated function.
 *          y          location n which to return the value of
 *                     the interpolated function at <x>.
 *
 *************************************************************************/
void InterpolateL(real8 *xa, real8 *ya, int numPoints, real8 x, real8 *y)
{
        int   low, mid, high;
        real8 a, b;

/*
 *      Just a sanity check on the input value.
 */
        if ((x < xa[0]) || (x > xa[numPoints-1])) {
            Fatal("InterpolateL(): specified X outside interpolation range");
        }

/*
 *      Find the indices of the values in <xa> that bracket the
 *      specified <x> value.
 */
        low = 0;
        high = numPoints-1;

        while ((high - low) > 1) {
            mid = (high + low) >> 1;
            if (xa[mid] > x) {
                high = mid;
            } else {
                low = mid;
            }
        }

/*
 *      low and high now bracket the specified <x> value
 */
        a = (xa[high] - x) / (xa[high] - xa[low]);
        b = 1.0 - a;

        *y = (a * ya[low]) + (b * ya[high]);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     InterpolateBL
 *      Description:  Use bilinear interpolation to calculate
 *                    the value of a function f(x,y) based on
 *                    the values of the function provided in
 *                    array <val> such that val[N] = f(x[N],y[N])
 *
 *
 *                    NOTE: This function assumes that for
 *                    each unique value in <x> there will
 *                    be an identical number of <y> values
 *                    and that the same set of <y> values are
 *                    used for each unique <x> value as in
 *                    the following:
 *
 *                      x[0] = xval1,  y[0] = yval1,  val[0] = value1
 *                      x[1] = xval1,  y[1] = yval2,  val[0] = value2
 *                      x[2] = xval1,  y[2] = yval3,  val[0] = value3
 *                      x[3] = xval2,  y[3] = yval1,  val[0] = value4
 *                      x[4] = xval2,  y[4] = yval2,  val[0] = value5
 *                      x[5] = xval2,  y[5] = yval3,  val[0] = value6
 *                      x[6] = xval3,  y[6] = yval1,  val[0] = value7
 *                      x[7] = xval3,  y[7] = yval2,  val[0] = value8
 *                      ...
 *
 *      Arguments
 *          x         Array of length <numElem> such that
 *                    x[0] <= x[1] <= x[2] <= ... <= x[numElem-1]
 *          y         Array of length <numElem> as defined in
 *                    description above
 *          val       Array of <numElem> values specifying the
 *                    value of the desired function at the various
 *                    combinations of <x> and <y> values.
 *          numElem   Number of elements in arrays <x>, <y> and <val>
 *          targetX   X value at which the function should be evaluated
 *          targetY   Y value at which the function should be evaluated
 *          result    Pointer to location in which to return
 *                    value of f(targetX,TargetY) to caller.
 *
 *      Returns: 1 if successful
 *               0 on failure
 *
 *------------------------------------------------------------------------*/
int InterpolateBL(real8 *x, real8 *y, real8 *val, int numElem,
                  real8 targetX, real8 targetY, real8 *result)
{
        int    i;
        real8  t, u;
        real8  lowX, highX, lowY, highY;
        real8  val1, val2, val3, val4, tmpVal;


/*
 *      If the target X value is outside the range in the
 *      table, return an error.
 */
        if ((targetX < x[0]) || (targetX > x[numElem-1])) {
            return(0);
        }

/*
 *      Find the min and max indices of the samples matching the x
 *      values bracketing the target x value on both the high and
 *      low ends.
 */
        int xLowMin  = 0;
        int xHighMin = 0;

        for (tmpVal = x[0], i = 0; i < numElem; i++) {
            if (tmpVal != x[i]) {
                tmpVal = x[i];
                if (tmpVal < targetX) {
                    xLowMin = i;
                } else if (tmpVal >= targetX) {
                    xHighMin = i;
                    break;
                }
            }
        }

        int xLowMax  = xHighMin-1;
        int xHighMax = xHighMin;

        for (i = xHighMin; i < numElem; i++) {
            if (x[xHighMin] == x[i]) xHighMax = i;
            else break;
        }

/*
 *      Find the indices of the samples bracketing the target y
 *      value within the samples at the lower x value.  If the
 *      target Y value is outside the range provided for the
 *      low X value, return an error.
 */
        if ((targetY < y[xLowMin]) || (targetY > y[xLowMax])) {
            return(0);
        }

        int yLow = xLowMin;

        for (yLow = xLowMin, i = xLowMin; i <= xLowMax; i++) {
            if (targetY <= y[i]) break;
            yLow = i;
        }

        int yHigh = yLow + 1;

/*
 *      Save the sample values bracketing the target y at the
 *      lower x value.
 */
        lowX  = x[yLow];
        lowY  = y[yLow];
        val1  = val[yLow];
        val2  = val[yHigh];

/*
 *      Find the indices of the samples bracketing the target y
 *      value within the samples at the upper x value.  If the
 *      target Y value is outside the range provided for the
 *      high X value, return an error.
 */
        if ((targetY < y[xHighMin]) || (targetY > y[xHighMax])) {
            return(0);
        }

        for (yLow = xHighMin, i = xHighMin; i <= xHighMax; i++) {
            if (targetY <= y[i]) break;
            yLow = i;
        }

        yHigh = yLow + 1;

/*
 *      Save the sample values bracketing the target y value at the
 *      upper x value.
 */
        highX = x[yHigh];
        highY = y[yHigh];
        val3  = val[yLow];
        val4  = val[yHigh];

/*
 *      And do the interpolation...
 */
        t = (targetX - lowX) / (highX - lowX);
        u = (targetY - lowY) / (highY - lowY);

        *result = (val1 * (1.0-t) * (1.0-u)) +
                  (val3 * t * (1.0-u))       +
                  (val4 * t * u)             +
                  (val2 * (1.0-t) * u);

        return(1);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Normal
 *      Description:  Calculate and return to the caller the normal of
 *                    vector <a>
 *
 *      Arguments:
 *          a   Array containing the three element vector.
 *
 *------------------------------------------------------------------------*/
real8 Normal(real8 a[3])
{
        return( sqrt(a[0]*a[0]+a[1]*a[1]+a[2]*a[2]) );
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Normalize
 *      Description:  Normalize vector a
 *
 *      Arguments:
 *          ax, ay, az  Pointers to the three elements of vector a.
 *                      The normalized values will be returned to the
 *                      caller via these same pointers.
 *
 *------------------------------------------------------------------------*/
void Normalize(real8 *ax, real8 *ay, real8 *az)
{
        real8 a2, a;

        a2 = ((*ax)*(*ax) + (*ay)*(*ay) + (*az)*(*az));

        if (a2 > 0.0) {
            a = sqrt(a2);
            *ax /= a;
            *ay /= a;
            *az /= a;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     NormalizeVec
 *      Description:  Normalize vector a
 *
 *      Arguments:
 *          vec   Three-element vector to be normalized.  Contents
 *                are updated before control returned to the caller.
 *
 *------------------------------------------------------------------------*/
void NormalizeVec(real8 vec[3])
{
        real8 a2, a;

        a2 = (vec[0]*vec[0] + vec[1]*vec[1] + vec[2]*vec[2]);

        if (a2 > 0.0) {
            a = sqrt(a2);
            vec[0] /= a;
            vec[1] /= a;
            vec[2] /= a;
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     CrossVector
 *      Description:  Calculate the cross product of two vectors.
 *
 *      Arguments:
 *          a  3 element array containing components of first vector
 *          b  3 element array containing components of second vector
 *          c  3 element array in which the resulting cross product
 *
 *------------------------------------------------------------------------*/
void CrossVector(real8 a[3], real8 b[3], real8 c[3])
{
        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     NormalizedCrossVector
 *      Description:  Calculate the normalized cross product of
 *                    two vectors.
 *
 *      Arguments:
 *          a  3 element array containing components of first vector
 *          b  3 element array containing components of second vector
 *          c  3 element array in which the resulting cross product
 *
 *------------------------------------------------------------------------*/
void NormalizedCrossVector(real8 a[3], real8 b[3], real8 c[3])
{
        real8 mag;

        c[0] = a[1]*b[2] - a[2]*b[1];
        c[1] = a[2]*b[0] - a[0]*b[2];
        c[2] = a[0]*b[1] - a[1]*b[0];

        mag = sqrt(c[0]*c[0] + c[1]*c[1] + c[2]*c[2]);

        if (mag < 1.0e-10) {
            printf("Warning! NormalizedCrossVector() All zeros.\n");
            c[0] = 0.0;
            c[1] = 0.0;
            c[2] = 0.0;
        } else {
            c[0] = c[0] / mag;
            c[1] = c[1] / mag;
            c[2] = c[2] / mag;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Orthogonalize
 *      Description:  Orthogonalize vector a against b with the
 *                    results overwriting the original values
 *                    for vector a.  i.e  a = a - b(a.b)/(b.b)
 *
 *      Arguments:
 *          ax, ay, az   Pointers to the three elements of vector a.
 *                       The results of this function will be returned
 *                       to the caller via these same pointers.
 *          bx, by, bz   Three elements of vector b
 *
 *------------------------------------------------------------------------*/
void Orthogonalize(real8 *ax, real8 *ay, real8 *az,
                   real8 bx, real8 by, real8 bz)
{
        real8 b2, p;

        b2 = (bx*bx + by*by + bz*bz);

        if (b2 <= 0.0) {
            return;
        } else {
            p = ((*ax)*bx + (*ay)*by + (*az)*bz) / b2;
            *ax -= p*bx;
            *ay -= p*by;
            *az -= p*bz;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     xvector
 *      Description:  Calculates the cross product of two vectors.
 *
 *      Arguments:
 *          ax, ay, az  Components of first vector
 *          bx, by, bz  Components of second vector
 *          cx, cy, cz  Pointers to locations in which to return to
 *                      the caller the components of the cross product
 *                      of the two vectors.
 *
 *------------------------------------------------------------------------*/
void xvector(real8 ax, real8 ay, real8 az,
             real8 bx, real8 by, real8 bz,
             real8 *cx, real8 *cy, real8 *cz)
{
        *cx =  ay*bz - az*by;
        *cy =  az*bx - ax*bz;
        *cz =  ax*by - ay*bx;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     InitNodeArm
 *      Description:  Initialize all components of the specified arm
 *                    (segment) of a node.
 *
 *      Arguments:
 *          node   Pointer to the node
 *          armID  Index (0 offset) of the arm/segment to be initialized.
 *
 *------------------------------------------------------------------------*/
static void InitNodeArm(Node_t *node, int armID)
{
        node->nbrTag[armID].domainID = -1;
        node->nbrTag[armID].index = -1;
        node->burgX[armID] = 0.0;
        node->burgY[armID] = 0.0;
        node->burgZ[armID] = 0.0;
        node->armfx[armID] = 0.0;
        node->armfy[armID] = 0.0;
        node->armfz[armID] = 0.0;
        node->nx[armID] = 0.0;
        node->ny[armID] = 0.0;
        node->nz[armID] = 0.0;
        node->sigbLoc[3*armID  ] = 0.0;
        node->sigbLoc[3*armID+1] = 0.0;
        node->sigbLoc[3*armID+2] = 0.0;
        node->sigbRem[3*armID  ] = 0.0;
        node->sigbRem[3*armID+1] = 0.0;
        node->sigbRem[3*armID+2] = 0.0;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    AllocNodeArms
 *      Description: Allocate memory for 'n' arms in the specified
 *                   node structure.
 *
 *      Arguments:
 *          node       Pointer to node for which to allocate arms
 *
 *------------------------------------------------------------------------*/
void AllocNodeArms(Node_t *node, int n)
{

/*
 *      The node structure *may* have some arm arrays already allocated.
 *      If this is the case but the number of currently allocated
 *      arms is not what we need, we have to free the previously
 *      allocated arrays before allocating new ones.
 */
        if (node->numNbrs != n) {

            if (node->numNbrs > 0) {

                free(node->nbrTag);
                free(node->burgX);
                free(node->burgY);
                free(node->burgZ);
                free(node->armfx);
                free(node->armfy);
                free(node->armfz);
                free(node->nx);
                free(node->ny);
                free(node->nz);
                free(node->sigbLoc);
                free(node->sigbRem);
            }

            node->numNbrs = n;

            node->nbrTag = (Tag_t *)malloc(n*sizeof(Tag_t));
            node->burgX = (real8 *)malloc(n*sizeof(real8));
            node->burgY = (real8 *)malloc(n*sizeof(real8));
            node->burgZ = (real8 *)malloc(n*sizeof(real8));
            node->armfx = (real8 *)malloc(n*sizeof(real8));
            node->armfy = (real8 *)malloc(n*sizeof(real8));
            node->armfz = (real8 *)malloc(n*sizeof(real8));
            node->nx = (real8 *)malloc(n*sizeof(real8));
            node->ny = (real8 *)malloc(n*sizeof(real8));
            node->nz = (real8 *)malloc(n*sizeof(real8));

            node->sigbLoc = (real8 *)malloc(n*3*sizeof(real8));
            node->sigbRem = (real8 *)malloc(n*3*sizeof(real8));
        }

/*
 *      And just make sure the arm specific arrays are initialized
 *      as necessary before returning to the caller.
 */
        for (int i=0; i < node->numNbrs; i++) {
            InitNodeArm(node, i);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FreeNode
 *      Description:  Removes a local node from various node queues,
 *                    recycles the tag, adds the node to the free node
 *                    queue and  re-initializes the node's index in the
 *                    <nodeKeys> array.
 *
 *                    NOTE: We don't free the arm allocations since
 *                    there's a good chance they'll be the same the
 *                    next time.  If later it is determined the arm
 *                    count differs, the arm arrays will be reallocated
 *                    at that time.
 *      Arguments:
 *          index   Index in the local <nodeKeys> array of the pointer
 *                  to the node to be freed.
 *
 *------------------------------------------------------------------------*/
void FreeNode(Home_t *home, int index)
{
	RemoveNodeFromCellQ(home, home->nodeKeys[index]);
	RemoveNodeFromCell2Q(home, home->nodeKeys[index]);
	RecycleNodeTag(home, index);
	PushFreeNodeQ (home, home->nodeKeys[index]);
	home->nodeKeys[index] = (Node_t *)NULL;

	return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     FreeNodeArms
 *      Description:  Release all the specified node's arm related
 *                    arrays, zero out the pointers, and set the arm.
 *                    count to zero.
 *      Arguments:
 *          node    Pointer to node structure in which to free
 *                  all arm related arrays.
 *
 *------------------------------------------------------------------------*/
void FreeNodeArms(Node_t *node)
{
        if (node->numNbrs == 0) {
            return;
        }

        free(node->nbrTag);      node->nbrTag = (Tag_t *)NULL;
        free(node->burgX);       node->burgX = (real8 *)NULL;
        free(node->burgY);       node->burgY = (real8 *)NULL;
        free(node->burgZ);       node->burgZ = (real8 *)NULL;
        free(node->armfx);       node->armfx = (real8 *)NULL;
        free(node->armfy);       node->armfy = (real8 *)NULL;
        free(node->armfz);       node->armfz = (real8 *)NULL;
        free(node->nx);          node->nx = (real8 *)NULL;
        free(node->ny);          node->ny = (real8 *)NULL;
        free(node->nz);          node->nz = (real8 *)NULL;
        free(node->sigbLoc);     node->sigbLoc = (real8 *)NULL;
        free(node->sigbRem);     node->sigbRem = (real8 *)NULL;

        node->numNbrs = 0;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     GetNeighborNode
 *      Description:  Given a node pointer, return the pointer to
 *                    the n'th neighbor of the node.
 *
 *      Arguments:
 *          node    Pointer to a node.
 *          n       Index in the node's neighbor list of the
 *                  neighbor to return to the caller.
 *
 *      Returns:    Pointer to the requested neighbor if found.
 *                  NULL in all other cases.
 *
 * FIX ME!  Check if this function will ever be called with sparse arm list
 *
 *------------------------------------------------------------------------*/
Node_t *GetNeighborNode(Home_t *home, Node_t *node, int n)
{
        Node_t *neighbor=0;
#if 0
/*
 *      New version which assumes no invalid arms in arm list
 */
        if (n >= node->numNbrs) {
            printf("GetNeighborNode: Error finding neighbor %d\n", n);
            PrintNode(node);
            return((Node_t *)NULL);
        }

        neighbor = GetNodeFromTag(home, node->nbrTag[n]);

        return(neighbor);
#else
/*
 *      Old version which assumes the arm list may be sparsely
 *      populated and returns the n'th valid neighbor, which may
 *      not be at index n.
 */
        int j = -1;

        for (int i=0; i < node->numNbrs; i++) {

            if (node->nbrTag[i].domainID >= 0) j++;

            if (j == n) {
                neighbor = GetNodeFromTag(home, node->nbrTag[i]);
                return(neighbor);
            }
        }

        printf("GetNeighborNode returning NULL for node (%d,%d) nbr %d\n",
               node->myTag.domainID, node->myTag.index, n);
        PrintNode(node);

        return((Node_t *)NULL);
#endif
}


/*-------------------------------------------------------------------------
 *
 *      Function:     GetNodeFromIndex
 *      Description:  Given the domain ID and index return
 *                    a pointer to the corresponding node.
 *
 *      Arguments:
 *          domID   ID of the domain owning the node in question
 *          index   Index in the <nodeKeys> array for domain <domID>
 *                  of the node in question.
 *
 *      Returns:    Pointer to the requested node if found.
 *                  NULL in all other cases.
 *
 *------------------------------------------------------------------------*/
Node_t *GetNodeFromIndex(Home_t *home, int domID, int index)
{
        Node_t         *node;
        RemoteDomain_t *remDom;


        if ((domID < 0) || (index < 0)) {
            return((Node_t *)NULL);
        }

/*
 *      If the node in question is owned by the current domain,
 *      look up the node in the local <nodeKeys> array.
 */
        if (domID == home->myDomain) {

            if (index >= home->newNodeKeyPtr) {
                return((Node_t *)NULL);
            }

            if ((node = home->nodeKeys[index]) == (Node_t *)NULL) {
                return((Node_t *)NULL);
            }

            return(node);

        } else {
/*
 *          Node is owned by a remote domain, so look up the
 *          node in the appropriate <remoteDomainKeys> array.
 */
            remDom = home->remoteDomainKeys[domID];

            if (remDom == (RemoteDomain_t *)NULL) {
                return((Node_t *)NULL);
            }

            if (index >= remDom->maxTagIndex) {
                return((Node_t *)NULL);
            }

            if ((node = remDom->nodeKeys[index]) == (Node_t *)NULL) {
                return((Node_t *)NULL);
            }

            return(node);
        }
}


/*-------------------------------------------------------------------------
 *
 *      Function:     GetNodeFromTag
 *      Description:  Given a node tag, returns a pointer to the
 *                    corresponding node structure.
 *
 *                    NOTE: If the specified tag is for a local node
 *                    and the local node is not found, the function
 *                    will force a code abort with a fatal error.
 *
 *      Arguments:
 *          tag    Tag identifying the desired node
 *
 *      Returns:    Pointer to the requested node if found.
 *                  NULL in all other cases.
 *
 *------------------------------------------------------------------------*/
Node_t *GetNodeFromTag (Home_t *home, Tag_t tag)
{
        Node_t         *node;
        RemoteDomain_t *remDom;

/*
 *      If the tag is for a local domain, look up the node in
 *      the local <nodeKeys> array.
 */
        if (tag.domainID == home->myDomain) {

            if (tag.index >= home->newNodeKeyPtr) {
                return((Node_t *)NULL);
            }

            node = home->nodeKeys[tag.index];
            return(node);

        } else {
/*
 *          If the node is in a remote domain, there are valid situations
 *          in which the current domain does not have information on
 *          either the remote doamin or the remote node.  Hence, it
 *          is not an error to return a NULL pointer.
 */
            remDom = home->remoteDomainKeys[tag.domainID];

            if (remDom == NULL) {
                return((Node_t *)NULL);
            }

            if (tag.index >= remDom->maxTagIndex) {
                return((Node_t *)NULL);
            }

            node = remDom->nodeKeys[tag.index];
            return(node);
        }
}


/*-------------------------------------------------------------------------
 *
 *      Function:     PrintNode
 *      Description:  For the specified node, print some interesting
 *                    items of data.
 *
 *------------------------------------------------------------------------*/
void PrintNode(Node_t *node)
{
        if (node == (Node_t *)NULL) return;

#if 1
        printf("  node(%d,%d) arms %d, cst=%d",
               node->myTag.domainID, node->myTag.index,
               node->numNbrs,node->constraint);
#if 0
        for (int i=0; i < node->numNbrs; i++) {
            printf("(%d,%d) ", node->nbrTag[i].domainID,
                   node->nbrTag[i].index);
        }
#endif
        printf("\n");
#endif

#if 1
/*
 *      Print the nodal position
 */
        printf("  node(%d,%d) position = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index,
               node->x, node->y, node->z);
#endif

#if 0
/*
 *      Print the nodal velocity and total node force
 */
        printf("  node(%d,%d) v = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index,
               node->vX, node->vY, node->vZ);
        printf("  node(%d,%d) f = (%.15e %.15e %.15e)\n",
               node->myTag.domainID, node->myTag.index,
               node->fX, node->fY, node->fZ);
#endif

#if 0
/*
 *      Print the arm specific forces
 */
        for (int i=0; i < node->numNbrs; i++) {
            printf("  node(%d,%d) arm[%d]-> (%d %d) f = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID, node->myTag.index, i,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->armfx[i], node->armfy[i], node->armfz[i]);
        }
#endif

#if 1
/*
 *      Print the burger's vector for each arm of the node
 */
        for (int i=0; i < node->numNbrs; i++) {
            printf("  node(%d,%d) arm[%d]-> (%d %d) b = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID, node->myTag.index, i,
                   node->nbrTag[i].domainID, node->nbrTag[i].index,
                   node->burgX[i], node->burgY[i], node->burgZ[i]);
        }
#endif

#if 0
/*
 *      Print the glide plane normal for each arm of the node
 */
        for (int i=0; i < node->numNbrs; i++) {
            printf("  node(%d,%d) arm[%d]-> (%d %d) n = (%.15e %.15e %.15e)\n",
                   node->myTag.domainID,node->myTag.index, i,
                   node->nbrTag[i].domainID,node->nbrTag[i].index,
                   node->nx[i],node->ny[i],node->nz[i]);
        }
#endif

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReallocNodeArms
 *      Description: Reallocate memory for 'n' arms in the specified
 *                   node structure preserving any previously
 *                   existing arm data.
 *
 *      Arguments:
 *          node       Pointer to node for which to allocate arms
 *
 *------------------------------------------------------------------------*/
void ReallocNodeArms(Node_t *node, int n)
{
        int i, origNbrCnt;

        origNbrCnt = node->numNbrs;

        node->numNbrs = n;

        node->nbrTag = (Tag_t *) realloc(node->nbrTag, sizeof(Tag_t)*n);
        node->burgX = (real8 *) realloc(node->burgX, sizeof(real8  )*n);
        node->burgY = (real8 *) realloc(node->burgY, sizeof(real8  )*n);
        node->burgZ = (real8 *) realloc(node->burgZ, sizeof(real8  )*n);
        node->armfx = (real8 *) realloc(node->armfx, sizeof(real8  )*n);
        node->armfy = (real8 *) realloc(node->armfy, sizeof(real8  )*n);
        node->armfz = (real8 *) realloc(node->armfz, sizeof(real8  )*n);
        node->nx = (real8 *) realloc(node->nx, sizeof(real8  )*n);
        node->ny = (real8 *) realloc(node->ny, sizeof(real8  )*n);
        node->nz = (real8 *) realloc(node->nz, sizeof(real8  )*n);
        node->sigbLoc  = (real8 *) realloc(node->sigbLoc, n*3*sizeof(real8));
        node->sigbRem  = (real8 *) realloc(node->sigbRem, n*3*sizeof(real8));

/*
 *      And just initialize the newly allocated arms only, leaving
 *      the previously existing arms as they were.
 */
        for (i = origNbrCnt; i < node->numNbrs; i++) {
            InitNodeArm(node, i);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     GetFreeNodeTag
 *      Description:  Return a free node tag.  If possible, return a
 *                    recycled tag. If there are no recycled tags available
 *                    return a new tag at the end of the nodeKeys array.
 *                    If there are no tags available there either, expand
 *                    nodeKeys first, then return new tag.
 *
 *------------------------------------------------------------------------*/
int GetFreeNodeTag(Home_t *home)
{
        int tagIndex;

        tagIndex = GetRecycledNodeTag(home);

        if (tagIndex < 0) {

            if (home->newNodeKeyPtr == home->newNodeKeyMax) {
                home->newNodeKeyMax += NEW_NODEKEY_INC;
                home->nodeKeys = (Node_t **) realloc(home->nodeKeys,
                home->newNodeKeyMax * sizeof(Node_t *));
            }

            tagIndex = home->newNodeKeyPtr++;
        } else if (tagIndex >= home->newNodeKeyPtr) {
            home->newNodeKeyPtr = tagIndex + 1;
        }

        return(tagIndex);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetRecycledNodeTag
 *      Description: Return the lowest available recycled node tag
 *                   from the recycled tag heap.
 *
 *      Returns:  -1 if no recycled tag is available, otherwise returns the
 *                lowest available recycled tag index.
 *
 *------------------------------------------------------------------------*/
int GetRecycledNodeTag(Home_t *home)
{
        int tagIndex;

        tagIndex = HeapRemove(home->recycledNodeHeap,
                              &home->recycledNodeHeapEnts);

        return(tagIndex);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    RecycleNodeTag
 *      Description: Add the specified tag to the list of available
 *                   recycled nodes.  List is maintained as a heap.
 *
 *      Arguments:
 *          tagIndex  Index of the local node tag to be recycled.
 *
 *------------------------------------------------------------------------*/
void RecycleNodeTag(Home_t *home, int tagIndex)
{
/*
 *      Add the indicated tag to the recycle heap.  If the heap is
 *      not large enough it's size will automatically be increased
 */
        HeapAdd(&home->recycledNodeHeap, &home->recycledNodeHeapSize,
                &home->recycledNodeHeapEnts, tagIndex);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     NodeDistance
 *      Description:  Compute the distance between two nodes
 *
 *      Arguments:
 *          node1   Pointer to first node structure.
 *          node2   Pointer to second node structure.
 *
 *      Returns: Distance between the nodes
 *
 *------------------------------------------------------------------------*/
static real8 NodeDistance(Home_t *home, Node_t *node1, Node_t *node2)
{
        real8   dx, dy, dz, dr;
        Param_t *param;

        param = home->param;

        dx = node1->x - node2->x;
        dy = node1->y - node2->y;
        dz = node1->z - node2->z;

        ZImage(param, &dx, &dy, &dz);

        dr = sqrt(dx*dx + dy*dy + dz*dz);

        return(dr);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ChangeArmBurg
 *      Description:  This function performs 1 of two tasks; it either
 *                    resets the burgers vector for an arm between
 *                    the two specified nodes (removing redundant arms
 *                    while doing so), or removes the arm completely.
 *                    The behavior depends on the provided burgers vector.
 *
 *                    NOTE: This only removes/modifies arms for node1.
 *                    A second call to the function is required to fix
 *                    the arm in the other direction.
 *
 *      Arguments:
 *          node1        Pointers to the two nodes we're interested
 *                       in.  Only arms of node1 will be modified
 *                       by this call.
 *          tag2         Pointer to tag indicating the node terminating
 *                       the arms of node1 that we're interested in.
 *          bx, by, bz   New burgers vector for the arm from
 *                       node1 to node2.  If zero, any arms of node1
 *                       terminating at node2 will be removed.  If non-
 *                       zero, the burgers vector for the first arm
 *                       from node1 to node2 will be rest to the
 *                       provided value, and any redundant arms from
 *                       node1 to node2 will be removed.
 *          nx, ny, nz   Glide plane normal for the new segment.  If
 *                       the provided burgers vector is zero, these
 *                       are ignored.  If the burgers vector is
 *                       non-zero, and the glide plane values are
 *                       zero, the function will explicitly re-
 *                       calculate the correct glide plane, otherwise
 *                       the provided values are used.
 *          log          Set to 1 if this operation is an operation
 *                       that will have to be passed on to remote
 *                       domains
 *          del_seg_factor  Factor specifying the portion of the
 *                       deleted segment that should be accumulated
 *                       in the running total of deleted segments.
 *
 *------------------------------------------------------------------------*/
void ChangeArmBurg(Home_t *home, Node_t *node1, Tag_t *tag2,
                   real8 bx, real8 by, real8 bz, real8 nx, real8 ny, real8 nz,
                   int log, real8 del_seg_factor)
{
        int     i, nc;
        real8   gpDotgp, eps = 1.0e-12;
        real8   burg[3], vel[3], dir[3], glidePlane[3];
        Node_t  *node2 = NULL;
        Param_t *param;

        if (node1 == NULL) return;

        param = home->param;

        if (log) {
            node2 = GetNodeFromIndex(home, tag2->domainID, tag2->index);
        }

/*
 *      If we're not removing the link (i.e. provided burgers vector
 *      is non-zero) and the glide plane is zero, calculate the glide
 *      plane for the node1->node2 segment. (When remote domains
 *      perfrom this operation, the correct glide plane will have
 *      been provided via the function arguments and will not need
 *      to be recomputed.)
 */
        if (((bx != 0) || (by != 0) || (bz != 0)) &&
            ((nx == 0.0) && (ny == 0.0) && (nz == 0.0))) {

            if (!log) {
                node2 = GetNodeFromIndex(home, tag2->domainID, tag2->index);
            }

            if (node2 == (Node_t *)NULL) {
/*
                printf("WARNING: Task %d -- ChangeArmBurg() "
                       "Can't compute new\nglide plane for segment "
                       "(%d,%d)->(%d,%d); 2nd node not found!\n",
                       home->myDomain, node1->myTag.domainID,
                       node1->myTag.index, tag2->domainID,
                       tag2->index);
*/
            } else {

                dir[X] = node1->x - node2->x;
                dir[Y] = node1->y - node2->y;
                dir[Z] = node1->z - node2->z;

                ZImage(param, &dir[X], &dir[Y], &dir[Z]);

                burg[X] = bx;
                burg[Y] = by;
                burg[Z] = bz;

                FindPreciseGlidePlane(home, burg, dir, glidePlane,
                                      home->param->allowFuzzyGlidePlanes);

                gpDotgp = DotProduct(glidePlane, glidePlane);

                if (gpDotgp < eps) {
                    vel[X] = (node1->vX + node2->vX) * 0.5;
                    vel[Y] = (node1->vY + node2->vY) * 0.5;
                    vel[Z] = (node1->vZ + node2->vZ) * 0.5;

                    FindPreciseGlidePlane(home, vel, dir, glidePlane,
                                          home->param->allowFuzzyGlidePlanes);

                    gpDotgp = DotProduct(glidePlane, glidePlane);

                    if (gpDotgp < eps) {
                        VECTOR_ZERO(glidePlane);
                    }

/*
 *                  If segment is screw, glide plane will still be undefined
 *                  so pick a plane appropriate to the burgers vector
 */
                    if (DotProduct(glidePlane, glidePlane) < 1.0e-03) {
                        PickScrewGlidePlane(home, burg, glidePlane);
                    }
                }

                nx = glidePlane[X];
                ny = glidePlane[Y];
                nz = glidePlane[Z];
            }
        }

/*
 *      If this operation needs to be sent to remote domains (ie. the
 *      <log> variable is set), add the operation to the list to be
 *      sent to the remote domains.
 */
        if (log) {

            glidePlane[X] = nx;
            glidePlane[Y] = ny;
            glidePlane[Z] = nz;

            burg[X] = bx;
            burg[Y] = by;
            burg[Z] = bz;

            AddOpChangeBurg(home, &node1->myTag, tag2, burg, glidePlane);
        }

        nc = 0;

/*
 *      Loop over node1's arms looking for all arms terminating at node2.
 */
        for (i=0; i < node1->numNbrs; i++) {

            if ((node1->nbrTag[i].domainID != tag2->domainID) ||
                (node1->nbrTag[i].index != tag2->index)) {
                continue;
            }

            nc++;

/*
 *          If the provided burgers vector is zero, remove the connection
 *          from node1 to node2.
 */
            if ((bx == 0) && (by == 0) && (bz == 0)) {

                node1->nbrTag[i].domainID = -1;
                node1->nbrTag[i].index = -1;
                node1->burgX[i] = 0;
                node1->burgY[i] = 0;
                node1->burgZ[i] = 0;
/*
 *              Accumulate half of deleted segment (since this routine
 *              is called twice - once for each direction). Only on
 *              original domain.
 */
                if (log) {
                    home->param->delSegLength +=
                            del_seg_factor * NodeDistance(home, node1, node2);
                }

            } else {
/*
 *              Burger's vector is not zero...
 *
 *              If this is the first link we found between the two
 *              nodes, update the burgers vector (and other necessary
 *              values).  Subsequent links are redundant and will be
 *              removed.
 */
                if (nc == 1) {
                    node1->burgX[i] = bx;
                    node1->burgY[i] = by;
                    node1->burgZ[i] = bz;

                    node1->sigbLoc[3*i]   = 0;
                    node1->sigbLoc[3*i+1] = 0;
                    node1->sigbLoc[3*i+2] = 0;
                    node1->sigbRem[3*i]   = 0;
                    node1->sigbRem[3*i+1] = 0;
                    node1->sigbRem[3*i+2] = 0;

                    node1->nx[i] = nx;
                    node1->ny[i] = ny;
                    node1->nz[i] = nz;
                } else {
/*
 *                  This is not the first link between the two
 *                  nodes, so this link is redundant; remove it.
 *                  (Note: this removes the link in one direction;
 *                  the link in the other direction is be removed
 *                  by a separate call to this function.)
 */
                    node1->nbrTag[i].domainID = -1;
                    node1->nbrTag[i].index = -1;
                    node1->burgX[i] = 0;
                    node1->burgY[i] = 0;
                    node1->burgZ[i] = 0;

/*
 *                  Accumulate half of deleted segment (since
 *                  this routine is called twice - once for each
 *                  direction). Only on original domain
 */
                    if (log) {
                        home->param->delSegLength +=
                                del_seg_factor *
                                NodeDistance(home,node1,node2);
                    }
                }
            }
        }

/*
 *      Compress out any null (i.e. domainID == -1) arms in the
 *      node's arm lists
 */
        CompressArmLists(node1);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ChangeConnection
 *      Description:  Looks for a node1->node2 arm.  If found, the
 *                    terminating end of the arm in node1 is changed from
 *                    node2 to node3.
 *
 *      Arguments:
 *          node1    Pointer to node to be modified
 *          tag2     Tag indicating node from which to disassociate
 *                   the node1 arm.
 *          tag3     Tag indicating target node to which to reattach
 *                   the node1 arm.
 *          log      Set to 1 if this operation is an operation
 *                   that will have to be passed on to remote
 *                   domains
 *
 *      Returns:   0 on success, -1 if the domain has no knowledge
 *                 of one or more of the specified nodes, or no
 *                 matching arm is found.
 *
 *------------------------------------------------------------------------*/
int ChangeConnection(Home_t *home, Node_t *node1, Tag_t *tag2,
                     Tag_t *tag3, int log)
{
        if (node1 == NULL) { return(-1); }

/*
 *      Given the current rules for topology changes, there
 *      are certain valid sequences of operations that can
 *      result in the formation of doubly-linked nodes (which
 *      is a bad thing).  So, we flag here any local nodes
 *      involved in this operation so we can later check
 *      if it has double links to any neighbors.
 */
        if (node1->myTag.domainID == home->myDomain) {
            node1->flags |= NODE_CHK_DBL_LINK;
        }

        if (tag3->domainID == home->myDomain) {
            Node_t *node3 = GetNodeFromTag(home, *tag3);
            if (node3 != (Node_t *)NULL) {
                node3->flags |= NODE_CHK_DBL_LINK;
            }
        }

        if (log) {
            AddOpChangeConn(home, &node1->myTag, tag2, tag3);
        }

        int domID = tag2->domainID;
        int index = tag2->index;

        for (int i=0; i < node1->numNbrs; i++) {

            if ((domID == node1->nbrTag[i].domainID) &&
                (index == node1->nbrTag[i].index)) {

                node1->nbrTag[i].domainID = tag3->domainID;
                node1->nbrTag[i].index    = tag3->index;

                node1->sigbLoc[3*i]   = 0;
                node1->sigbLoc[3*i+1] = 0;
                node1->sigbLoc[3*i+2] = 0;
                node1->sigbRem[3*i]   = 0;
                node1->sigbRem[3*i+1] = 0;
                node1->sigbRem[3*i+2] = 0;

                return(0);
            }
        }

        return(-1);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    CompressArmLists
 *      Description: During certain topological operations, multiple
 *                   arms of a node may be combined (or annihilate
 *                   each other).  During this process, one or more
 *                   of the node's arms are flagged as such, but
 *                   are left in place during the operation.  Once
 *                   the operation is complete this function is invoked
 *                   to remove any now obsolete arms from a node.
 *
 *      Arguments:
 *          node       Pointer to node for which to compress the
 *                     arm list.
 *
 *------------------------------------------------------------------------*/
void CompressArmLists (Node_t *node)
{
        int i, j, found;

        int nulls=0;

        for (i=0; i < node->numNbrs; i++) {
/*
 *          Deal with any null arms (identified by a negative
 *          domain ID)
 */
            if (node->nbrTag[i].domainID == -1) {

                nulls++;
/*
 *              find the next non-null entry, and move it to i
 */
                found = 0;

                for (j = i+1; j < node->numNbrs; j++) {

                    if (node->nbrTag[j].domainID != -1) {
                        found = 1;
                        node->nbrTag[i].domainID = node->nbrTag[j].domainID;
                        node->nbrTag[i].index = node->nbrTag[j].index;
                        node->burgX[i] = node->burgX[j];
                        node->burgY[i] = node->burgY[j];
                        node->burgZ[i] = node->burgZ[j];
                        node->nx[i] = node->nx[j];
                        node->ny[i] = node->ny[j];
                        node->nz[i] = node->nz[j];
                        node->armfx[i] = node->armfx[j];
                        node->armfy[i] = node->armfy[j];
                        node->armfz[i] = node->armfz[j];
                        node->sigbLoc[3*i]  =node->sigbLoc[3*j];
                        node->sigbLoc[3*i+1]=node->sigbLoc[3*j+1];
                        node->sigbLoc[3*i+2]=node->sigbLoc[3*j+2];
                        node->sigbRem[3*i]  =node->sigbRem[3*j];
                        node->sigbRem[3*i+1]=node->sigbRem[3*j+1];
                        node->sigbRem[3*i+2]=node->sigbRem[3*j+2];
                        node->nbrTag[j].domainID = -1;
                        node->nbrTag[j].index = -1;

                        break;
                    }
                }

                if (!found) break;
            }
        }

/*
 *      If no more non-null arms were found, then arm i-1 is the
 *      last non-null, and there are i non-nulls all together.
 *      If there were no null arms, we're done; otherwise, we
 *      have to reallocate the node's arm arrays to the proper
 *      size and free the now unused arms.
 */
        if (!nulls) {
            return;
        }

        ReallocNodeArms(node, i);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       DomainOwnsSeg
 *      Description:    Determines if the specified domain owns
 *                      the segment beginning in the current domain
 *                      and terminating at the domain indicated by
 *                      <endTag>.
 *
 *                      Note:  Ownership of segments crossing
 *                             domain boundaries alternates on
 *                             and even/odd cycles.  Additionally,
 *                             the ownership is dependent on the
 *                             class of topological operation being
 *                             considered (i.e. during remesh operations
 *                             the ownership is the exact reverse of
 *                             ownership during collision handling.)
 *
 *      Arguments:
 *              opClass         Class of topological operation during which
 *                              this function is being invoked.  Valid
 *                              values are:
 *                                  OPCLASS_SEPARATION
 *                                  OPCLASS_COLLISION
 *                                  OPCLASS_REMESH
 *              thisDomain      domain containing the first endpoint of
 *                              the segment in question.
 *              endTag          pointer to the tag of the second endpoint
 *                              of the segment.
 *
 *      Returns:  1 if <thisDomain> owns the segment
 *                0 in all other cases.
 *
 *-------------------------------------------------------------------------*/
int DomainOwnsSeg(Home_t *home, int opClass, int thisDomain, Tag_t *endTag)
{
        int ownsSeg = 1;

/*
 *      If both endpoints are in the same domain, the domain owns
 *      the segment.
 */
        if (thisDomain == endTag->domainID) {
            return(1);
        }

/*
 *      For collision handling and node separations, ownership
 *      of segments crossing domain boundaries is the lower
 *      numbered domain on even numbered cycles and the higher
 *      numbered domain for odd numbered cycles.
 *
 *      For remesh operations, ownership rules are the opposite
 *      of those used for collision handling.
 */
        switch (opClass) {
        case OPCLASS_SEPARATION:
        case OPCLASS_COLLISION:
            if (home->cycle & 0x01)
                ownsSeg = (thisDomain > endTag->domainID);
            else
                ownsSeg = (thisDomain < endTag->domainID);
            break;
        case OPCLASS_REMESH:
            if (home->cycle & 0x01)
                ownsSeg = (thisDomain < endTag->domainID);
            else
                ownsSeg = (thisDomain > endTag->domainID);
            break;
        default:
            Fatal("Invalid opClass %d in DomainOwnsSeg()", opClass);
            break;
        }

        return(ownsSeg);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     InsertArm
 *      Description:  Add a new arm in nodeA pointing to nodeB with given
 *                    burgers vector and glide plane.
 *
 *      Arguments:
 *          nodeA        pointer to node to which to add a a new arm
 *          nodeBtag     pointer to the tag of the node at which the
 *                       new arm terminates
 *          bx, by, bz   burgers vector for the new arm
 *          nx, ny, nz   glide plane for the new arm
 *          log          set to 1 if this operation is an operation
 *                       that will have to be passed on to remote
 *                       domains
 *
 *------------------------------------------------------------------------*/
void InsertArm (Home_t *home, Node_t *nodeA, Tag_t *nodeBtag,
                real8 bx, real8 by, real8 bz,
                real8 nx, real8 ny, real8 nz, int log)
{
        int      i;
        Node_t  *nodeB;

        if (log) {
            real8 burg[3], plane[3];

            burg[0] = bx;
            burg[1] = by;
            burg[2] = bz;

            plane[0] = nx;
            plane[1] = ny;
            plane[2] = nz;

            AddOpInsertArm(home, &nodeA->myTag, nodeBtag, burg, plane);
        }

/*
 *      Allocate another arm and set the burgers vector and glide plane
 */
        i = nodeA->numNbrs;

        ReallocNodeArms(nodeA, i+1);

        nodeA->nbrTag[i].domainID = nodeBtag->domainID;
        nodeA->nbrTag[i].index    = nodeBtag->index;

        nodeA->burgX[i] = bx;
        nodeA->burgY[i] = by;
        nodeA->burgZ[i] = bz;

        nodeA->nx[i] = nx;
        nodeA->ny[i] = ny;
        nodeA->nz[i] = nz;

/*
 *      Given the current rules for topology changes, there
 *      are certain valid sequences of operations that can
 *      result in the formation of doubly-linked nodes (which
 *      is a bad thing).  So, we flag here any local nodes
 *      involved in this operation so we can later check
 *      if it has double links to any neighbors.
 */
        if (nodeA->myTag.domainID == home->myDomain) {
            nodeA->flags |= NODE_CHK_DBL_LINK;
        }

        if (nodeBtag->domainID == home->myDomain) {
            nodeB = GetNodeFromTag(home, *nodeBtag);
            if (nodeB != (Node_t *)NULL) {
                nodeB->flags |= NODE_CHK_DBL_LINK;
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     MarkNodeForceObsolete
 *      Description:  Sets a flag for the specified node indicating that
 *                    the nodal force and velocity values are obsolete and
 *                    must be recalculated.  If the node is not local to
 *                    the current domain, the operation list going to
 *                    the remote domains will be modified to include a
 *                    operation to force the owning domain to recalculate
 *                    the node's force and velocity.
 *
 *      Arguments:
 *          node    pointer to node for which the force and velocity
 *                  values must be recomputed
 *
 *------------------------------------------------------------------------*/
void MarkNodeForceObsolete(Home_t *home, Node_t *node)
{

#ifndef _BGQ
#ifdef _OPENMP
#pragma omp atomic
#endif
        node->flags |= NODE_RESET_FORCES;
#else
#ifdef _OPENMP
#pragma omp critical (MARK_NODE_FORCE_OBSOLETE)
#endif
        {
            node->flags |= NODE_RESET_FORCES;
        }
#endif

/*
 *      If the node is locally owned, we're done.  If not, we
 *      need to let the owning domain know it needs to recompute
 *      the force/velocity.
 */
        if (node->myTag.domainID == home->myDomain) {
            return;
        }

        AddOpMarkForceObsolete(home, &node->myTag);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     RecalcSegGlidePlane
 *      Description:  Calculate and reset (if necessary) the glide plane
 *                    normal for a segment.  This function will only update
 *                    the information locally.  It is assumed the information
 *                    will be passed to remote domains elsewhere as needed.
 *
 *      Arguments:
 *          node1, node2   pointers to nodal endpoints of the segment
 *          ignoreIfScrew  Toggle.  If set to 1, the function will
 *                         leave the segment glide plane as is if the
 *                         segment is screw.  Otherwise it will pick
 *                         an appropriate glide plane for the bugers vector
 *
 *------------------------------------------------------------------------*/
void RecalcSegGlidePlane(Home_t *home, Node_t *node1, Node_t *node2,
                         int ignoreIfScrew)
{
        int node1SegID, node2SegID;
        real8 burg[3], lineDir[3], newPlane[3];

/*
 *      Just a couple quick sanity checks.
 */
        if ((node1 == (Node_t *)NULL) ||
            (node2 == (Node_t *)NULL) ||
            (node1 == node2)) {
            return;
        }

/*
 *      It is possible that the two nodes are not really connected.
 *      This can happen in cases such as MeshCoarsen() where a node is
 *      removed leaving two nodes doubly linked, and when the double links
 *      get reconciled they annihilate each other, leaving the two nodes
 *      unconnected.  In situations like that, there's nothing for this
 *      function to do...
 */
        if (!Connected(node1, node2, &node1SegID)) {
            return;
        }

        node2SegID = GetArmID(node2, node1);

        burg[X] = node1->burgX[node1SegID];
        burg[Y] = node1->burgY[node1SegID];
        burg[Z] = node1->burgZ[node1SegID];

        lineDir[X] = node2->x - node1->x;
        lineDir[Y] = node2->y - node1->y;
        lineDir[Z] = node2->z - node1->z;

        ZImage(home->param, &lineDir[X], &lineDir[Y], &lineDir[Z]);
        NormalizeVec(lineDir);

        FindPreciseGlidePlane(home, burg, lineDir, newPlane,
                              home->param->allowFuzzyGlidePlanes);

        if (DotProduct(newPlane, newPlane) < 1.0e-03) {
            if (ignoreIfScrew) return;
            PickScrewGlidePlane(home, burg, newPlane);
        }

        Normalize(&newPlane[X], &newPlane[Y], &newPlane[Z]);

        node1->nx[node1SegID] = newPlane[X];
        node1->ny[node1SegID] = newPlane[Y];
        node1->nz[node1SegID] = newPlane[Z];

        node2->nx[node2SegID] = newPlane[X];
        node2->ny[node2SegID] = newPlane[Y];
        node2->nz[node2SegID] = newPlane[Z];

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ResetSegForces
 *      Description:  Reset the segment forces on nodeA for the
 *                    segment terminating at nodeB, then re-sum the
 *                    total nodal forces for nodeA
 *
 *      Arguments:
 *          nodeA       pointer to node for which to reset seg forces
 *          nodeBtag    pointer to the tag of the node at which the
 *                      segment in question terminates
 *          fx,fy,fz    segment forces
 *          globalOp    set to 1 if this operation is an operation
 *                      that will have to be passed on to remote
 *                      domains
 *
 *------------------------------------------------------------------------*/
void ResetSegForces(Home_t *home, Node_t *nodeA, Tag_t *nodeBtag,
                    real8 fx, real8 fy, real8 fz, int globalOp)
{
/*
 *      If other domains need to be notified of this operation, add the
 *      action to the operation list
 */
        if (globalOp) {
            real8 f1[3];

            f1[0] = fx;
            f1[1] = fy;
            f1[2] = fz;

            AddOpResetSegForce1(home, &nodeA->myTag, nodeBtag, f1);
        }

/*
 *      Locate the segment of nodeA terminating at nodeb and update
 *      forces for that segment.
 */
        for (int i=0; i < nodeA->numNbrs; i++) {
            if ((nodeA->nbrTag[i].domainID == nodeBtag->domainID) &&
                (nodeA->nbrTag[i].index == nodeBtag->index)) {
                nodeA->armfx[i] = fx;
                nodeA->armfy[i] = fy;
                nodeA->armfz[i] = fz;
                break;
            }
        }

/*
 *      Reset the total forces for nodeA based on all its segment forces
 */
        nodeA->fX = 0;
        nodeA->fY = 0;
        nodeA->fZ = 0;

        for (int i=0; i < nodeA->numNbrs; i++) {
            nodeA->fX += nodeA->armfx[i];
            nodeA->fY += nodeA->armfy[i];
            nodeA->fZ += nodeA->armfz[i];
        }

#ifndef _BGQ
#ifdef _OPENMP
#pragma omp atomic
#endif
        nodeA->flags |= NODE_RESET_FORCES;
#else
#ifdef _OPENMP
#pragma omp critical (RESET_SEG_FORCES)
#endif
        {
            nodeA->flags |= NODE_RESET_FORCES;
        }
#endif

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ResetSegForces2
 *      Description:  Reset the segment forces for segment nodeA--nodeB
 *                    in (potentially) both nodeA and nodeB, resumming
 *                    total nodal forces after the segment forces
 *                    are updated.
 *
 *     Arguments:
 *         nodeA         pointer to first nodal endpoint of the segment.
 *         nodeBtag      pointer to the tag of the second nodal
 *                       endpoint for the segment.
 *         fx1,fy1,f1z   Segment forces at nodeA
 *         fx2,fy2,f2z   Segment forces at nodeB
 *         globalOp      set to 1 if this operation is an operation
 *                       that will have to be passed on to remote
 *                       domains
 *
 *------------------------------------------------------------------------*/
void ResetSegForces2(Home_t *home, Node_t *nodeA, Tag_t *nodeBtag,
                   real8 f1x, real8 f1y, real8 f1z,
                   real8 f2x, real8 f2y, real8 f2z, int globalOp)
{
        int    i, nodeALocal, nodeBLocal;
        Node_t *nodeB;

/*
 *      If this is a global operation being sent to remote domains
 *      for processing, add the action to the operation list
 */
        if (globalOp) {
            real8 f1[3], f2[3];

            f1[0] = f1x;
            f1[1] = f1y;
            f1[2] = f1z;

            f2[0] = f2x;
            f2[1] = f2y;
            f2[2] = f2z;

            AddOpResetSegForce2(home, &nodeA->myTag, nodeBtag, f1, f2);

            return;
        }

        nodeALocal = (nodeA->myTag.domainID == home->myDomain);
        nodeBLocal = (nodeBtag->domainID == home->myDomain);

/*
 *      If both the nodes are local, it's because the current domain
 *      modified forces on both nodes, and this is a local operation
 *      only, hence there's no need to do anything else.
 *
 *      In all other cases, requests to update forces on remote
 *      nodes will be processed, and any request to update the
 *      forces for a node local to this domain will force the
 *      node to be flagged for a complete reevaluation of all
 *      its segment forces and its velocity.
 */
        if (nodeALocal && nodeBLocal) {
            return;
        }

/*
 *      Handle updates on nodeA forces if necessary.
 */
        if (nodeALocal) {
            nodeA->flags |= NODE_RESET_FORCES;
        } else {
            for (i = 0; i < nodeA->numNbrs; i++) {
                if ((nodeA->nbrTag[i].domainID == nodeBtag->domainID) &&
                    (nodeA->nbrTag[i].index == nodeBtag->index)) {
                    nodeA->armfx[i] = f1x;
                    nodeA->armfy[i] = f1y;
                    nodeA->armfz[i] = f1z;
                    break;
                }
            }
/*
 *          Reset the total forces for nodeA based on all its segment forces
 */
            nodeA->fX = 0;
            nodeA->fY = 0;
            nodeA->fZ = 0;

            for (i = 0; i < nodeA->numNbrs; i++) {
                nodeA->fX += nodeA->armfx[i];
                nodeA->fY += nodeA->armfy[i];
                nodeA->fZ += nodeA->armfz[i];
            }
        }

/*
 *      Handle updates on nodeB forces if necessary.  Note: it is
 *      possible that the local domain owned nodeB and deleted it
 *      via its own topological changes, so make sure it really
 *      still exists before messing with it.  :-)
 */
        if ((nodeB = GetNodeFromTag(home, *nodeBtag)) == (Node_t *)NULL) {
            return;
        }

        if (nodeBLocal) {
            nodeB->flags |= NODE_RESET_FORCES;
        } else {
            if (nodeB == (Node_t *)NULL) return;
            for (i = 0; i < nodeB->numNbrs; i++) {
                if ((nodeB->nbrTag[i].domainID == nodeA->myTag.domainID) &&
                    (nodeB->nbrTag[i].index == nodeA->myTag.index)) {
                    nodeB->armfx[i] = f2x;
                    nodeB->armfy[i] = f2y;
                    nodeB->armfz[i] = f2z;
                    break;
                }
            }
/*
 *          Reset the total forces for nodeB based on all its segment forces
 */
            nodeB->fX = 0;
            nodeB->fY = 0;
            nodeB->fZ = 0;

            for (i = 0; i < nodeB->numNbrs; i++) {
                nodeB->fX += nodeB->armfx[i];
                nodeB->fY += nodeB->armfy[i];
                nodeB->fZ += nodeB->armfz[i];
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     SubtractSegForce
 *      Description:  Subtract from the total force on node1, the force
 *                    contribution from segment node1/node2 and zero
 *                    out the segment specific force at node1 for the
 *                    node1/node2 segment.
 *
 *      Arguments:
 *          node1   pointer to node for which force values
 *                  are to be updated
 *          node2   pointer to the node terminating the arm
 *                  from node1
 *
 *      TBD:  Some functions calling this routine already know
 *            the ID of the arm connecting the two segments, so
 *            change <armID> to a parameter and only use
 *            GetArmID() if the provided ID is < zero.
 *
 *------------------------------------------------------------------------*/
void SubtractSegForce(Node_t *node1, Node_t *node2)
{
        int  armID;

        if ((armID = GetArmID(node1, node2)) < 0) return;

        node1->fX -= node1->armfx[armID];
        node1->fY -= node1->armfy[armID];
        node1->fZ -= node1->armfz[armID];

        node1->armfx[armID] = 0.0;
        node1->armfy[armID] = 0.0;
        node1->armfz[armID] = 0.0;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ResetGlidePlane
 *      Description:  Change the glide plane for the dislocation segment
 *                    between two nodes to the specified plane.
 *
 *      Arguments:
 *          newPlane  Array containing the X,Y and Z components of the
 *                    new glide plane
 *
 *------------------------------------------------------------------------*/
void ResetGlidePlane(Home_t *home, real8 newPlane[3], Tag_t *tag1,
                     Tag_t *tag2, int globalOp)
{
        int    armID;
        Node_t *node;

/*
 *      If this operation has to be communicated for implentation in
 *      remote domains, add it to the remote operation queue.
 */
        if (globalOp) {
            AddOpResetPlane(home, tag1, tag2, newPlane);
        }

/*
 *      The local domain may not have information about both endpoints
 *      of the segments, but we'll update the glide plane info in
 *      whichever nodal endpoints it does know about.
 */
        node = GetNodeFromTag(home, *tag1);

        if (node != (Node_t *)NULL) {
            for (armID = 0; armID < node->numNbrs; armID++) {
                if ((node->nbrTag[armID].domainID == tag2->domainID) &&
                    (node->nbrTag[armID].index == tag2->index)) {
                    node->nx[armID] = newPlane[X];
                    node->ny[armID] = newPlane[Y];
                    node->nz[armID] = newPlane[Z];
                }
            }
        }

        node = GetNodeFromTag(home, *tag2);

        if (node != (Node_t *)NULL) {
            for (armID = 0; armID < node->numNbrs; armID++) {
                if ((node->nbrTag[armID].domainID == tag1->domainID) &&
                    (node->nbrTag[armID].index == tag1->index)) {
                    node->nx[armID] = newPlane[X];
                    node->ny[armID] = newPlane[Y];
                    node->nz[armID] = newPlane[Z];
                }
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     RepositionNode
 *      Description:  Change the position of a node and if necessary send
 *                    the operation to remote domains.
 *
 *      Arguments:
 *          newPlane  Array containing the X,Y and Z components of the
 *                    new glide plane
 *
 *------------------------------------------------------------------------*/
void RepositionNode(Home_t *home, real8 newPos[3], Tag_t *tag, int globalOp)
{
        Node_t *node;

/*
 *      If this operation has to be communicated for implementation in
 *      remote domains, add it to the remote operation queue.
 */
        if (globalOp) {
            AddOpResetCoord(home, tag, newPos);
        }

/*
 *      If this domain knows about the specified node, update its position.
 */
        node = GetNodeFromTag(home, *tag);

        if (node != (Node_t *)NULL) {
            node->x = newPos[X];
            node->y = newPos[Y];
            node->z = newPos[Z];
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     CollisionNodeOrder
 *      Description:  Compares two node tags and returns a
 *                    value corresponding to their ordering.
 *
 *                    NOTE: This function differs from OrderNodes() in
 *                    that the criteria for determining whether one node
 *                    tag is lower than another reverses on even/odd
 *                    timesteps.  This behavior is needed during collision
 *                    handling because node/segment ownership plays a
 *                    part in determining whether certain topological
 *                    operations are permitted.  By altering ownership,
 *                    a topological change prevented during one timestep
 *                    due to ownership restrictions is likely to be
 *                    permitted the next timestep.
 *
 *      Arguments:
 *          tagA   Pointer to first tag
 *          tagB   Pointer to second tag
 *
 *      Returns:  -1 if tagA is lower than tagB
 *                 0 if tagA is the same as tagB
 *                 1 if tagA is higher than tagB
 *
 *------------------------------------------------------------------------*/
int CollisionNodeOrder(Home_t *home, Tag_t *tagA, Tag_t *tagB)
{
/*
 *      If the tags belong to nodes int the same domain, standard
 *      ordering applies (i.e. the node with lower index has the
 *      lowest order.
 */
        if (tagA->domainID == tagB->domainID) {
            if (tagA->index < tagB->index) return(-1);
            return(tagA->index > tagB->index);
        }
/*
 *      For nodes in different domains, the criteria for ordering
 *      is the domainID, but the ordering reverses with each
 *      timestep.  This reversal allows the collision handling
 *      code to reverse ownership of segments crossing domain
 *      boundaries in an attempt to prevent nodes in a higher
 *      number domain with arms terminating in lower domains
 *      from being permanently exempt from various collisions.
 */
        if (home->cycle & 0x01) {
            if (tagA->domainID > tagB->domainID) return(-1);
            else return(1);
        } else {
            if (tagA->domainID > tagB->domainID) return(1);
            else return(-1);
        }
}


/*-------------------------------------------------------------------------
 *
 *      Function:     IntCompare
 *      Description:  Comparison function for two integers suitable
 *                    for use by system sorting and searching functions
 *                    such as qsort(), bsearch(), etc.
 *
 *      Arguments:
 *          a   Pointer to first integer
 *          b   Pointer to second integer.
 *
 *      Returns:  a negative value if  a <  b
 *                0                if  a == b
 *                a positive value if  a >  b
 *
 *------------------------------------------------------------------------*/
int IntCompare(const void *a, const void *b)
{
        int *val1 = (int *)a;
        int *val2 = (int *)b;

        return(*val1 - *val2);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     OrderNodes
 *      Description:  Compares two Node_t structures based on the
 *                    tag values.  This function is compatible
 *                    for use by system sorting and searching functions
 *                    such as qsort(), bsearch(), etc.
 *
 *                    IMPORTANT!  This behavior of this function is duplicated
 *                                by the OrderTags() function below which uses
 *                                tag pointers rather than node pointers.
 *                                These two functions must be kept in sync
 *                                so any behavior change in one should be
 *                                reflected in the other.
 *
 *      Arguments:
 *          a   Pointer to first node structure.
 *          b   Pointer to second node structure.
 *
 *      Returns:  -1 if  a <  b
 *                 0 if  a == b
 *                 1 if  a >  b
 *
 *------------------------------------------------------------------------*/
int OrderNodes(const void *a, const void *b)
{
        Node_t *node1 = (Node_t *)a;
        Node_t *node2 = (Node_t *)b;

        if (node1->myTag.domainID < node2->myTag.domainID) return(-1);
        if (node1->myTag.domainID > node2->myTag.domainID) return(1);
        if (node1->myTag.index < node2->myTag.index) return(-1);

        return(node1->myTag.index > node2->myTag.index);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     OrderTags
 *      Description:  Comparison function for two Tag_t structures.  This
 *                    function is compatible for use by system sorting and
 *                    searching functions such as qsort(), bsearch(), etc.
 *
 *                    IMPORTANT!  This function duplicates the behavior
 *                                of the OrderNodes() function above using
 *                                tag pointers rather than node pointers.
 *                                These two functions must be kept in sync
 *                                so any behavior change in one should be
 *                                reflected in the other.
 *
 *      Arguments:
 *          a   Pointer to first tag structure.
 *          b   Pointer to second tag structure.
 *
 *      Returns:  -1 if  a <  b
 *                 0 if  a == b
 *                 1 if  a >  b
 *
 *------------------------------------------------------------------------*/
int OrderTags(const void *a, const void *b)
{
        Tag_t *tag1 = (Tag_t *)a;
        Tag_t *tag2 = (Tag_t *)b;

        if (tag1->domainID < tag2->domainID) return(-1);
        if (tag1->domainID > tag2->domainID) return(1);
        if (tag1->index < tag2->index) return(-1);

        return(tag1->index > tag2->index);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Connected
 *      Description:  Determines whether or not two nodes are connected
 *                    by a dislocation segment.
 *
 *      Arguments:
 *          node1     Pointer to first node
 *          node2     Pointer to second node
 *          armID     Pointer to location in which to return to the
 *                    caller the index (for node1) of the arm connecting
 *                    the two nodes.  (If the nodes are not connected
 *                    the contents of this pointer will be set to -1)
 *
 *      Returns:   0 if the nodes are not connected.
 *                 1 if the nodes are connected.
 *
 *------------------------------------------------------------------------*/
int Connected(Node_t *node1, Node_t *node2, int *armID)
{
        int i, dom2, idx2;

        *armID = -1;

        if ((node1 == NULL) || (node2 == NULL)) {
            return(0);
        }

        dom2 = node2->myTag.domainID;
        idx2 = node2->myTag.index;

        for (i = 0; i < node1->numNbrs; i++) {
            if ((dom2 == node1->nbrTag[i].domainID) &&
                (idx2==node1->nbrTag[i].index)) {
                *armID = i;
                return(1);
            }
        }

        return 0;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Fatal
 *      Description:  Prints a user specified message, aborts all
 *                    other parallel tasks (if any) and self-terminates
 *
 *      Arguments:    This function accepts a variable number of arguments
 *                    in the same fashion as printf(), with the first
 *                    option being the format string which determines
 *                    how the remainder of the arguments are to be
 *                    interpreted.
 *
 *------------------------------------------------------------------------*/

void Fatal(const char *format, ...)
{
        char    msg[512];
        va_list args;

        va_start(args, format);
        vsnprintf(msg, sizeof(msg)-1, format, args);
        msg[sizeof(msg)-1] = 0;
        va_end(args);
        printf("Fatal: %s\n", msg);

#ifdef DEBUG_ABORT
        ParadisAbort_Set(1);
#else
        MPI_Abort(MPI_COMM_WORLD, -1);
        exit(1);
#endif
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Warning
 *      Description:  Prints a user specified warning message.
 *
 *      Arguments:    This function accepts a variable number of arguments
 *                    in the same fashion as printf(), with the first
 *                    option being the format string which determines
 *                    how the remainder of the arguments are to be
 *                    interpreted.
 *
 *------------------------------------------------------------------------*/
void Warning(const char *format, ...)
{
        char    msg[512];
        va_list args;

        va_start(args, format);
        vsnprintf(msg, sizeof(msg)-1, format, args);
        msg[sizeof(msg)-1] = 0;
        va_end(args);
        printf("Warning: %s\n", msg);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetArmID
 *      Description: Given two node pointers, return arm ID for node1
 *                   if node2 is its neighbor.
 *
 *      Arguments:
 *          node1    Pointer to first node structure
 *          node2    Pointer to second node structure
 *
 *      Returns:  Non-negative index of the arm of node1 terminating
 *                at node2 if the nodes are connected.
 *                -1 if the two nodes are not connected.
 *
 *------------------------------------------------------------------------*/
int GetArmID (Node_t *node1, Node_t *node2)
{
        int i, domID, index;

        if ((node1 == NULL) || (node2 == NULL)) {
            return(-1);
        }

        domID = node2->myTag.domainID;
        index = node2->myTag.index;

        for(i = 0; i < node1->numNbrs; i++) {
            if ((domID == node1->nbrTag[i].domainID) &&
                (index == node1->nbrTag[i].index)) {
                return(i);
            }
        }

        return(-1);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetBurgersVectorNormal
 *      Description: Given two node pointers, return the burgers vector
 *                   (from node1 to node2) and glide-plane normal for
 *                   the segment.  If the two nodes are not connected
 *                   all values returned to the caller will be zero'ed
 *
 *                   QUESTION: Should this function ever be called
 *                             with a mismatched pair of nodes?
 *                             We might want to just simplify the
 *                             function and abort if that happens...
 *
 *      Arguments:
 *          node1       Pointer to first node structure
 *          node2       Pointer to second node structure
 *          bx, by, bz  Pointers to the three components of the
 *                      burgers vector to be returned to the caller.
 *          nx, ny, nz  Pointers to the three components of the
 *                      glide plane normal to be returned to the caller.
 *
 *------------------------------------------------------------------------*/
void GetBurgersVectorNormal(Node_t *node1, Node_t *node2,
                            real8 *bx, real8 *by, real8 *bz,
                            real8 *nx, real8 *ny, real8 *nz)
{
        int i, domID, index;

/*
 *      If we don't have two valid node pointers, just return zeros
 */
        if ((node1 == NULL) || (node2 == NULL)) {

            *bx = 0.0;
            *by = 0.0;
            *bz = 0.0;

            *nx = 0.0;
            *ny = 0.0;
            *nz = 0.0;

            return;
        }

/*
 *      Got valid pointers, so find node2 in node1's neighbor list
 *      and return the associated burgers vector and normal.
 */
        domID = node2->myTag.domainID;
        index = node2->myTag.index;

        for (i = 0; i < node1->numNbrs; i++) {

            if ((domID == node1->nbrTag[i].domainID) &&
                (index == node1->nbrTag[i].index)) {

                *bx = node1->burgX[i];
                *by = node1->burgY[i];
                *bz = node1->burgZ[i];

                *nx = node1->nx[i];
                *ny = node1->ny[i];
                *nz = node1->nz[i];

                return;
            }
        }

        printf("GetBurgersVectorNormal failed: "
               "node (%d,%d) no neighbor match (%d,%d)\n",
               node1->myTag.domainID, node1->myTag.index,
               node2->myTag.domainID, node2->myTag.index);

        *bx = 0.0;
        *by = 0.0;
        *bz = 0.0;

        *nx = 0.0;
        *ny = 0.0;
        *nz = 0.0;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       Getline
 *      Description:    Read the next non-empty, non-comment line
 *                      from the specified file stream and return
 *                      it to the caller.
 *
 *      Arguments:
 *              string          pointer to location at which to return
 *                              the line to the caller
 *              len             size (in bytes) of <string>
 *              fp              stream from which to read the data
 *
 *-----------------------------------------------------------------------*/
void Getline(char *string, int len, FILE *fp)
{
        char    *s;

        while (1) {

                string[0] = 0;
                if ((s = fgets(string, len, fp)) == (char *)NULL) break;
                string[len-1] = 0;

/*
 *              If the line is a coment line, skip it. If
 *              it has any non-whitespace characters, return it to
 *              the caller, otherwise move to the next line.
 */
                if (*s =='#') continue;

                for ( ; *s != 0; s++) {
                        if ((*s != ' ') && (*s != '\t') && (*s != '\n')) break;
                }

                if (*s != 0) break;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       FindCellCenter
 *      Description:    Given a set of coordinates within a cell, or the
 *                      natural indices of a cell (i.e. not adjusted to
 *                      allow for ghost cells), return the coordinates
 *                      of the cell center to the caller.
 *
 *      Arguments:
 *          x,y,z    Either the coordinates of a point within a cell, or
 *                   the indices of the cell in question.  The <type>
 *                   argument determiens how these will be interpreted.
 *          type     Integer indicating the type of values contained
 *                   in the <x>, <y> and <z> parameters.  A value of
 *                   1 indicates they represent coordinates of a point
 *                   within a cell, a value of 2 menas they represent
 *                   the indices of a cell (indices must be in the
 *                   range 0 thru numcells for each dimension respectively
 *          xCenter,
 *          yCenter,
 *          zCenter  Locations in which to store the coordinates of the
 *                   cell center for return to the caller.
 *
 *-----------------------------------------------------------------------*/
void FindCellCenter(Param_t *param, real8 x, real8 y, real8 z, int type,
                    real8 *xCenter, real8 *yCenter, real8 *zCenter)
{
        int     cellIndex[3];
        real8   xStart, yStart, zStart;
        real8   minCoord[3], cellSize[3];

        cellIndex[0] = 0;
        cellIndex[1] = 0;
        cellIndex[2] = 0;


        cellSize[X] = (param->maxSideX - param->minSideX) / param->nXcells;
        cellSize[Y] = (param->maxSideY - param->minSideY) / param->nYcells;
        cellSize[Z] = (param->maxSideZ - param->minSideZ) / param->nZcells;

        minCoord[X] = param->minSideX;
        minCoord[Y] = param->minSideY;
        minCoord[Z] = param->minSideZ;

        if (type == 1) {
/*
 *          x, y and z are positions.  Find the cell encompassing the
 *          point and return the center point of the cell
 */
            cellIndex[X] = (int)((x - minCoord[X]) / cellSize[X]);
            cellIndex[Y] = (int)((y - minCoord[Y]) / cellSize[Y]);
            cellIndex[Z] = (int)((z - minCoord[Z]) / cellSize[Z]);


        } else if (type == 2) {
/*
 *          x, y and z are the indices of a cell. Range is 0 thru num{XYZ}Cells
 *          i.e. not adjusted for ghost cell
 */
            cellIndex[X] = (int)x;
            cellIndex[Y] = (int)y;
            cellIndex[Z] = (int)z;

        } else {
            Fatal("Unknown arugment type in FindCellCenter()!");
        }

        xStart = minCoord[X] + cellSize[X]*0.5;
        yStart = minCoord[Y] + cellSize[Y]*0.5;
        zStart = minCoord[Z] + cellSize[Z]*0.5;

        *xCenter = xStart + (cellIndex[X] * cellSize[X]);
        *yCenter = yStart + (cellIndex[Y] * cellSize[Y]);
        *zCenter = zStart + (cellIndex[Z] * cellSize[Z]);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       LocateCell
 *      Description:    Given a set of coordinates, return to the
 *                      caller the ID of the cell encompassing them
 *
 *                      Warning: <position> MUST be located within
 *                               the problem space!
 *
 *      Arguments:
 *          cellID   location in which to return the cellID to the caller
 *          position Coordinates to locate
 *
 *-----------------------------------------------------------------------*/
void LocateCell(Home_t *home, int *cellID, real8 coord[3])
{
        int     cIndex;
        int     numCellsX, numCellsY, numCellsZ;
        int     cellIndexX, cellIndexY, cellIndexZ;
        real8   cellMax, cellMin;
        real8   minCoordX, minCoordY, minCoordZ;
        real8   cellSizeX, cellSizeY, cellSizeZ;
        Param_t *param;

/*
 *      Set the lower limit of the base area (excluding possible
 *      periodic cells) and the size of each cell
 */
        param = home->param;

        minCoordX = param->minSideX;
        minCoordY = param->minSideY;
        minCoordZ = param->minSideZ;

        numCellsX = param->nXcells;
        numCellsY = param->nYcells;
        numCellsZ = param->nZcells;

        cellSizeX = (param->maxSideX - minCoordX) / numCellsX;
        cellSizeY = (param->maxSideY - minCoordY) / numCellsY;
        cellSizeZ = (param->maxSideZ - minCoordZ) / numCellsZ;

/*
 *      Find the indices of the cell encompassing the specified
 *      coordinates (adjusted for periodic cells), and
 *      calculate the cell ID.
 *
 *      First find the cell index in X
 */
        cellIndexX = numCellsX;

        for (cIndex = 0; cIndex <= numCellsX; cIndex++) {

            cellMax = rint((minCoordX) + (cIndex * cellSizeX));
            cellMin = rint((minCoordX) + (cIndex-1) * cellSizeX);

            if ((coord[X] >= cellMin) && (coord[X] < cellMax)) {
                cellIndexX = cIndex;
                break;
            }
        }

/*
 *      Find the cell index in Y
 */
        cellIndexY = numCellsY;

        for (cIndex = 0; cIndex <= numCellsY; cIndex++) {

            cellMax = rint((minCoordY) + (cIndex * cellSizeY));
            cellMin = rint((minCoordY) + (cIndex-1) * cellSizeY);

            if ((coord[Y] >= cellMin) && (coord[Y] < cellMax)) {
                cellIndexY = cIndex;
                break;
            }
        }

/*
 *      Find the cell index in Z
 */
        cellIndexZ = numCellsZ;

        for (cIndex = 0; cIndex <= numCellsZ; cIndex++) {

            cellMax = rint((minCoordZ) + (cIndex * cellSizeZ));
            cellMin = rint((minCoordZ) + (cIndex-1) * cellSizeZ);

            if ((coord[Z] >= cellMin) && (coord[Z] < cellMax)) {
                cellIndexZ = cIndex;
                break;
            }
        }

        *cellID = EncodeCellIdx(home, cellIndexX, cellIndexY, cellIndexZ);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       NodePinned
 *
 *      Description:    Determine if a node should be considered 'pinned'
 *                      during the a cross slip procedure.  Nodes are
 *                      to be pinned if any of the following are true:
 *
 *                        - the node is pinned in ANY dimension
 *                        - the node is owned by another domain
 *                        - the node has any segments in a plane other
 *                          than the one indicated by <planeIndex>.
 *                        - the node is attached to any segments owned
 *                          by a remote domain.
 *
 *      Arguments:
 *         node         pointer to node
 *         planeIndex   row index of a glide plane in <glidedir>
 *         glideDir     <numGlideDirs> X 3 matrix of possible glide planes.
 *                      Each row representing a single plane.
 *         numGlideDirs Number of glide planes in <glideDir>.
 *
 *      Returns:   1 if the node should be considered pinned and unmovable
 *                 0 in all other cases
 *
 *-------------------------------------------------------------------------*/
int NodePinned(Home_t *home, Node_t *node, int planeIndex,
               real8 (*glideDir)[3], int numGlideDir)
{
        int    i, k;
        int    planetest;
        real8  planetestmin, ptest, ptest2;
        real8  segplane[3];

/*
 *      If the node is pinned due to constraints no need to
 *      check anything else.
 *
 *      NOTE: For now, a node that is pinned in ANY dimension is treated
 *            here as if it is pinned in ALL dimensions.
 */
        if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
            return(1);
        }

/*
 *      If the node is not owned by the current domain, it
 *      may not be repositioned.
 */
        if (node->myTag.domainID != home->myDomain) {
            return(1);
        }

/*
 *      Loop over all segments attached to the node.
 */
        for (i = 0; i < node->numNbrs; i++) {

/*
 *          If the segment is owned by another domain, treat this node
 *          as 'pinned'... otherwise, the remote domain *could* use
 *          the segment in a collision but might not have the correct
 *          nodal position since no communication is done between
 *          the cross slip and collision handling procedures.
 *
 *          NOTE: during cross slip, we'll use the same segment
 *                ownership rules we use during collision handling.
 */
            if (!DomainOwnsSeg(home, OPCLASS_COLLISION, home->myDomain,
                               &node->nbrTag[i])) {
                return(1);
            }

/*
 *          Check the glide plane for the segment, and if the segment
 *          has a different glide plane index than <planeIndex>, the
 *          node should be considered 'pinned'
 */
            segplane[X] = node->nx[i];
            segplane[Y] = node->ny[i];
            segplane[Z] = node->nz[i];

            planetest = 0;
            planetestmin = 10.0;

            for (k = 0; k < numGlideDir; k++) {

                ptest = fabs(DotProduct(glideDir[k], segplane));
                ptest2 = ptest * ptest;

                if (ptest2 < planetestmin) {
                    planetest = k;
                    planetestmin = ptest2;
                }

            }

            if (planeIndex != planetest) {
                return(1);
            }
        }

        return(0);
}


/*---------------------------------------------------------------------------
 *
 *      Function:     randm
 *      Description:  This is a special function for random number
 *                    generation on 32-bit machines that do not
 *                    support long integer multiplication and truncation.
 *                    The technique used is to do the multiplication and
 *                    addition in parts, by splitting all integers into a
 *                    'high' and a 'low' part.  The algorithm is exact,
 *                    and should give machine-independent results.
 *
 *                    The algorithm implemented is (following D.E. Knuth):
 *                    seed = seed*1592653589 + 453816691
 *                    if (seed.lt.0) seed = seed + 1 + 2147483647
 *                    Note that 1592653589 = 48603*2**15 + 30485
 *
 *      Returns:      real8 value in the range  0.0 <= value <= 1.0
 *
 *-------------------------------------------------------------------------*/
real8 randm(int *seed)
{
        int ia, ib, i1, i2, i3;

        ia = *seed / 32768;
        ib = (*seed % 32768);
        i1 = ia * 30485;
        i2 = ib * 30485;
        i3 = ib * 48603;
        i1 = (i1 % 65536);
        i3 = (i3 % 65536);
        i1 = i1 + i3 + 13849 + i2 / 32768 + (ia % 2) * 32768;
        i2 = (i2 % 32768) + 12659;
        i1 = i1 + i2 / 32768;
        i2 = (i2 % 32768);
        i1 = (i1 % 65536);

        *seed = i1 * 32768 + i2;

        return(*seed * 4.65661287308E-10);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ReadTabulatedData
 *      Description:  Reads columns of data (interpreted as real8
 *                    values) from the specified file and returns
 *                    the data to the caller.  Lines containing
 *                    nothing but white-space character and lines
 *                    beginning with a "#" are ignored by this function.
 *
 *                    Note: Caller is responsible for freeing memory
 *                    allocated by this function.
 *
 *      Arguments:
 *          fileName  Name of file containing data
 *          numCols   Number of colums of data to read and return to
 *                    the caller.
 *          colData   Location in which to return to the caller an array
 *                    of <numCols> pointers, each pointer to an array
 *                    containing the values for the corresponding column.
 *                    The caller is responsible for freeing the pointers
 *                    in colData and colData itself.  i.e. free(colData[0]),
 *                    free(colData[1]), ... free(colData[numRows-1]),
 *                    free(colData)
 *          numRows   Location in which to return to the caller the number
 *                    or rows of data read (i.e. the length of each array
 *                    pointed to by colData).
 *
 *------------------------------------------------------------------------*/
void ReadTabulatedData(char *fileName, int numCols, real8 ***colData,
                       int *numRows)
{
        int    col, row, rowsAllocated;
        real8  **column;
        char   *token;
        char   line[512];
        FILE   *fp;

        if ((fp = fopen(fileName, "r")) == (FILE *)NULL) {
            Fatal("ReadTabulatedData(): open error %d on %s",
                   errno, fileName);
        }

        row = 0;
        rowsAllocated = 0;
        column = (real8 **)calloc(1, sizeof(real8 *) * numCols);

/*
 *      Loop until there is no more data to be read.  The Getline()
 *      function will skip any line containing nothing but white-space
 *      characters and any line beginning with a "#".
 */
        while (1) {

            Getline(line, sizeof(line) - 1, fp);

            if (line[0] == 0) break;  /* no more data */

/*
 *          If we need to add more rows of data, increase the
 *          array sizes by a static amount each time.
 */
            if (row >= rowsAllocated) {
                rowsAllocated += 200;
                for (col = 0; col < numCols; col++) {
                    column[col] = (real8 *)realloc(column[col],
                                                   rowsAllocated *
                                                   sizeof(real8));
                }
            }

/*
 *          Read the next row of data into the respective column arrays.
 *          Multiple items on a line are delimited by a <space> or <tab>.
 */
            token = strtok(line, " 	");

            for (col = 0; col < numCols; col++) {
                if ((token == (char *)NULL) || (*token == '\n')) {
                    Fatal("ReadTabulatedData(): file %s, too few columns of "
                          "data on line\n", fileName);
                }
                column[col][row] = atof(token);
                token = strtok(NULL, " 	");
            }

            row++;
        }

/*
 *      Truncate any extra elements from the column arrays
 */
        if (rowsAllocated > row) {
            for (col = 0; col < numCols; col++) {
                column[col] = (real8 *)realloc(column[col],
                                               sizeof(real8) * row);
            }
        }

        *colData = column;
        *numRows = row;

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:     Sign
 *      Description:  Similar to the MATLAB sign() function.
 *
 *      Returns:     -1 if x <  zero
 *                    0 if x == zero
 *                   +1 if x >  zero
 *
 *----------------------------------------------------------------------*/
real8 Sign(const real8 x)
{
   return ( (x<0.0) ? -1.0 : (x>0.0) ? 1.0 : 0.0 );
}


/*------------------------------------------------------------------------
 *
 *      Function:     StrEquiv
 *      Description:  Do a case-insensitive comparison of two character
 *                    strings.
 *
 *      Arguments:
 *          s1, s2   pointers to the NULL-terminated strings to be compared
 *
 *      Returns:     1 if the only difference between the strings is the
 *                     case of 1 or more of the characters.  NOTE: if
 *                     both  string pointers passed to the function are
 *                     NULL, the function returns 1 indicating equivalence.
 *                   0 in all other cases.
 *
 *----------------------------------------------------------------------*/
int StrEquiv(const char *ps1, const char *ps2)
{
   char *s1 = (char *) ps1;
   char *s2 = (char *) ps2;

/*
 *      If both pointers are NULL treat as equivalent.  If only 1 of the
 *      pointers is NULL, treat as non-equivalent.
 */
        if ((s1 == (char *)NULL) && (s2 == (char *)NULL)) {
            return(1);
        } else if ((s1 == (char *)NULL) || (s2 == (char *)NULL)) {
            return(0);
        }

/*
 *      Loop though each character in s1 and do a case-insensitive
 *      comparison to the corresponding characfter in s2.  If there's
 *      a mismatch, return 0.
 */
        while (*s1 != 0) {
            if ((*s1 == toupper(*s2)) || (*s1 == tolower(*s2))) {
                s1++;
                s2++;
                continue;
            }

            return(0);
        }

/*
 *      If there are any remaining unchecked characters in s2, the
 *      strings are not equivalent, otherwise they are.
 */
        return(*s2 == 0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     Uniq
 *      Description:  Compress a sorted list of integers, removing all
 *                    duplicate numbers
 *
 *      Arguments:
 *          list    Array of integers to be compressed
 *          count   Pointer to integer indicating the number of
 *                  elements in the <list> array.  On return
 *                  to the caller, this value will have been
 *                  adjusted to allow for any elements removed
 *                  from the array.
 *
 *------------------------------------------------------------------------*/
void Uniq (int *list, int *count)
{
        int ptrIn, ptrOut, cntIn;

        ptrIn = 0;
        ptrOut = 0;
        cntIn = *count;

        while (ptrIn < cntIn) {
/*
 *          Copy from input pointer to output pointer
 */
            list[ptrOut] = list[ptrIn];

/*
 *          Skip redundant elements
 */
            while ((ptrIn < cntIn) && (list[ptrIn] == list[ptrOut])) {
                ptrIn++;
            }

            ptrOut++;
        }

        *count = ptrOut;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     NodeHasSessileBurg
 *      Description:  Determine if any segment of a node has a burgers
 *                    vector/line sense pair that has been specified by
 *                    by the user as sessile (i.e. via the control file
 *                    parameters)
 *
 *                    NOTE: The control file specification of sessile
 *                    burgers vector/line sense is done such that for
 *                    each specified burgers vector there is a corresponding
 *                    line sense.  If all components of the line sense
 *                    are zero, any segment with the specified burgers
 *                    vector is sessile.  If the line sense has non-zero
 *                    components, the segment is only considered sessile
 *                    if the burgers vector and the line sense match one
 *                    of the user-specified burgers vector/line sense pairs
 *
 *      Arguments:
 *          node    Pointer to node to be checked.
 *
 *      Returns:  0 if the node has no segments with sessile burgers vectors
 *                1 if it does.
 *
 *------------------------------------------------------------------------*/
int NodeHasSessileBurg(Home_t *home, Node_t *node)
{
        int     i, nbrIndex, numNbrs, nss;
        real8   eps = 1.0e-08;
        real8   mag, invMag;
        real8   bx, by, bz;
        real8   dx, dy, dz;
        real8   ssbx, ssby, ssbz, sslx, ssly, sslz;
        Node_t  *nbrNode;
        Param_t *param;

        param = home->param;
        nss = (int)param->sessileburgspec[0];

/*
 *      If there are no burgers vector/line sense pairs flagged
 *      as sessile, just return.  Otherwise loop through all the
 *      pairs and check if any of the segments of the node are of
 *      the specified burgers vector/line sense
 */
        if (nss <= 0) {
            return(0);
        }
        numNbrs = node->numNbrs;

        for (i = 0; i < nss; i++) {

            ssbx = param->sessileburgspec[i*3+1];
            ssby = param->sessileburgspec[i*3+2];
            ssbz = param->sessileburgspec[i*3+3];

            sslx = param->sessilelinespec[i*3+1];
            ssly = param->sessilelinespec[i*3+2];
            sslz = param->sessilelinespec[i*3+3];

            for (nbrIndex = 0; nbrIndex < numNbrs; nbrIndex++) {

                bx = node->burgX[nbrIndex];
                by = node->burgY[nbrIndex];
                bz = node->burgZ[nbrIndex];

                if (((bx ==  ssbx) && (by ==  ssby) && (bz ==  ssbz)) ||
                    ((bx == -ssbx) && (by == -ssby) && (bz == -ssbz))) {
/*
 *                  Burgers vector is a match.  If the corresponding line
 *                  sense is unspecified, all segments with this burgers
 *                  vector are sessile, so return true.  Otherwise, check
 *                  the line sense of the segment against the specification
 */
                    if ((sslx == 0) && (ssly == 0) && (sslz == 0)) {

                        return(1);

                    } else {

                        nbrNode = GetNeighborNode(home, node, nbrIndex);

                        if (nbrNode == (Node_t *)NULL) {
                            Fatal("Neighbor not found at %s line %d\n",
                                   __FILE__, __LINE__);
                        }

                        dx = nbrNode->x - node->x;
                        dy = nbrNode->y - node->y;
                        dz = nbrNode->z - node->z;

                        ZImage(param, &dx, &dy, &dz) ;

                        mag = sqrt(dx*dx + dy*dy + dz*dz);

/*
 *                      It's possible we have a zero-length segment.  If
 *                      so, skip the segment.
 */
                        if (mag < eps) {
                            continue;
                        }

                        invMag = 1.0 / mag;

                        dx *= invMag;
                        dy *= invMag;
                        dz *= invMag;

                        if (((fabs(dx-sslx) < 1e-2) &&
                             (fabs(dy-ssly) < 1e-2) &&
                             (fabs(dz-sslz) < 1e-2)) ||
                            ((fabs(dx+sslx) < 1e-2) &&
                             (fabs(dy+ssly) < 1e-2) &&
                             (fabs(dz+sslz) < 1e-2))) {
                            return(1);
                        }
                    }
                }
            }
        }

        return(0);
}


/*-------------------------------------------------------------------------
 *
 *      Function:     GetGlideConstraintList
 *      Description:  Find all the independent glide constraints for
 *                    a node based on the glide plane normals for all
 *                    segments attached to the node.
 *
 *      Arguments:
 *          IN:  node                Pointer to node for which glide
 *                                   constraints are being computed
 *          OUT: numGlideConstraints Number of glide constraints (plane
 *                                   normal vectors) returned in
 *                                   <glideConstraints>.  Maximum number
 *                                   of constraints is 3.
 *          OUT: glideConstraints    Array in which the glide constraints
 *                                   are returned to the caller.
 *
 *------------------------------------------------------------------------*/
void GetGlideConstraintList(Node_t *node,
                            int *numGlideConstraints,
                            real8 (*glideConstraints)[3])
{
        int     i;
        int     numConstraints;
        real8   eps = 1.0e-06;

        numConstraints = 0;

/*
 *      Define the linearly independent glide constraints
 */
        for (i = 0; i < node->numNbrs; i++) {
            real8   glidePlane[3], gpNormal;
            real8   temp, testVal, dpnnm;

            glidePlane[X] = node->nx[i];
            glidePlane[Y] = node->ny[i];
            glidePlane[Z] = node->nz[i];

            gpNormal = Normal(glidePlane);

            switch(numConstraints) {
                case 0:
                    glideConstraints[0][X] = glidePlane[X] / gpNormal;
                    glideConstraints[0][Y] = glidePlane[Y] / gpNormal;
                    glideConstraints[0][Z] = glidePlane[Z] / gpNormal;
                    numConstraints = 1;
                    break;

                case 1:
                    temp = glideConstraints[0][X] * glidePlane[X] +
                           glideConstraints[0][Y] * glidePlane[Y] +
                           glideConstraints[0][Z] * glidePlane[Z];

                    testVal = 1.0 - fabs(temp) / gpNormal;

                    if (testVal > eps) {
                        dpnnm = (glidePlane[X] * glideConstraints[0][X] +
                                 glidePlane[Y] * glideConstraints[0][Y] +
                                 glidePlane[Z] * glideConstraints[0][Z]);

                        glideConstraints[1][X] = glidePlane[X] -
                                                 dpnnm * glideConstraints[0][X];
                        glideConstraints[1][Y] = glidePlane[Y] -
                                                 dpnnm * glideConstraints[0][Y];
                        glideConstraints[1][Z] = glidePlane[Z] -
                                                 dpnnm * glideConstraints[0][Z];

                        Normalize(&glideConstraints[1][X],
                                  &glideConstraints[1][Y],
                                  &glideConstraints[1][Z]);

                        cross(glideConstraints[0], glideConstraints[1],
                              glideConstraints[2]);
                        numConstraints = 2;
                    }
                    break;

                case 2:
                    temp = glideConstraints[2][X] * glidePlane[X] +
                           glideConstraints[2][Y] * glidePlane[Y] +
                           glideConstraints[2][Z] * glidePlane[Z];

                    testVal = fabs(temp) / gpNormal;

                    if (testVal > eps) {
                        numConstraints = 3;
                    }

                    break;

            }  /* end switch(numConstraints) */

        }

        *numGlideConstraints = numConstraints;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ApplyConstraintsToVelocity
 *
 *      Description:  Apply the provided glide constraints to the
 *                    specified velocity vector.  If there are more
 *                    than 2 glide constraints, the velocity is
 *                    explicitly zeroed.
 *
 *      Arguments:
 *          IN:     numGlideConstraints Number of glide constraints (plane
 *                                      normal vectors) specified in
 *                                      <glideConstraints>.  Maximum number
 *                                      of constraints is 3.
 *
 *          IN:     glideConstraints    Array of <numGlideConstraints> glide
 *                                      constraints
 *
 *          IN/OUT: velocity            Caller-supplied velocity to which
 *                                      the glide constraints are to be
 *                                      applied.
 *
 *------------------------------------------------------------------------*/
void ApplyConstraintsToVelocity(int numGlideConstraints,
                                real8 (*glideConstraints)[3],
                                real8 velocity[3])
{
        real8 S1;

/*
 *      Apply glide constraints to node's velocity.  Note: assumes
 *      glideConstraints matrix is in the same frame (laboratory or
 *      crystalographic) as the velocity.
 */
        switch(numGlideConstraints) {
            case 1:
                S1 = DotProduct(velocity, glideConstraints[0]);
                velocity[0] -= S1 * glideConstraints[0][0];
                velocity[1] -= S1 * glideConstraints[0][1];
                velocity[2] -= S1 * glideConstraints[0][2];
                break;

            case 2:
                S1 = DotProduct(velocity, glideConstraints[2]);
                velocity[0] = S1 * glideConstraints[2][0];
                velocity[1] = S1 * glideConstraints[2][1];
                velocity[2] = S1 * glideConstraints[2][2];
                break;

            case 3:
                VECTOR_ZERO(velocity);
                break;
        }

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     ApplyShearVelocityCap
 *
 *      Description:  Cap the velocity to the shear sound velocity
 *
 *------------------------------------------------------------------------*/
void ApplyShearVelocityCap(Home_t *home, real8 velocity[3])
{
        if (home->param->shearVelocity <= 0.0) return;

        real8 alpha = 10.0;
        real8 vmax = home->param->shearVelocity; // m/s

        real8 vmag = Normal(velocity) * home->param->burgMag; // m/s
        if (vmag < 1e-5) return;

        real8 vcap = vmag / pow(1.0 + pow(vmag/vmax, alpha), 1.0/alpha);

        velocity[0] *= (vcap / vmag);
        velocity[1] *= (vcap / vmag);
        velocity[2] *= (vcap / vmag);

        return;
}

/*-------------------------------------------------------------------------
 *
 *      Function:    Get2ndMaxIndex
 *
 *      Description: Given a list of values and the index of the highest
 *                   value, find the index of the 2nd highest value.
 *
 *      Arguments:
 *          IN: A            Array of values to be searched
 *          IN: numElements  Number of elements in <A>
 *          IN: indexOfMax   Index in <A> of the largest value
 *
 *      Returns:
 *          -1 if <A> has only 1 element, otherwise returns the
 *          index in <A> of the largest value other than A[indexOfMax].
 *
 *------------------------------------------------------------------------*/
int Get2ndMaxIndex(real8 *A, int numElements, int indexOfMax)
{
        int   indexOf2ndMax, m;

        indexOf2ndMax = (numElements < 2) ? -1 : ((indexOfMax == 0) ? 1 : 0);

        if (indexOf2ndMax >= 0) {
            for (m = 0; m < numElements; m++) {
                if ((A[m] > A[indexOf2ndMax]) && (m != indexOfMax)) {
                    indexOf2ndMax = m;
                }
            }
        }

        return indexOf2ndMax;
}

int Get2ndMaxIndex(real8 *A, int numElements, int indexOfMax1, int indexOfMax2, int indexOfMax3)
{
        bool  found = false;
        int   indexOf2ndMax, m;
        int   enough_data;
        real8 tmp_A = 0.0;

        enough_data = (numElements < 2) ? -1 : 1; //((indexOfMax == 0) ? 1 : 0);

        if (enough_data >= 0)
        {
            for (m = 0; m < numElements; m++)
            {
                real8 Am = A[m];

                if (     (fabs(Am) > tmp_A)
                      && (m != indexOfMax1)
                      && (m != indexOfMax2)
                      && (m != indexOfMax3) )
                {
                    found         = true;
                    tmp_A         = fabs(Am);
                    indexOf2ndMax = m;
                }
            }
        }

        if(!found)
        indexOf2ndMax = -1;

        return indexOf2ndMax;
}

/*-------------------------------------------------------------------------
 *
 *      Function:     PrintAllNodes
 *      Description:  For all nodes present in the simulation, print some
 *                    interesting items of data.
 *
 *------------------------------------------------------------------------*/
void PrintAllNodes(Home_t *home)
{
  int     inode,numNodes;
  int     i;
  Node_t *node;


  numNodes = home->newNodeKeyPtr;


  for (inode = 0; inode < numNodes; inode++)
    {
      if ((node = home->nodeKeys[inode]) == (Node_t *)NULL) {
	continue;
      }

      for (i = 0; i < node->numNbrs; i++)
	{
	  if (fabs(node->burgX[i] * node->burgY[i]) < 1e-5)
	    {
	      printf("  node(%d,%d) b = (%.15e %.15e %.15e) n = (%.15e %.15e %.15e)\n",
		     node->myTag.domainID, node->myTag.index,
		     node->burgX[i], node->burgY[i], node->burgZ[i],
		     node->nx[i],node->ny[i],node->nz[i]);

	    }
	}
    }
}


/*-------------------------------------------------------------------------
 *
 *      Function:     ThreeExtremaIndex
 *      Description:  Find min, max indices of a three indice vector
 *                    Is there a more elegant way of doing this?
 *
 *------------------------------------------------------------------------*/
void ThreeExtremaIndex(real8 x, real8 y, real8 z, int *imax,int *imin)
{
  *imin = ( (x<y) ? ( (x<z) ? 0 : 2 ) : ( (y<z) ? 1 : 2 ) );
  *imax = ( (x>y) ? ( (x>z) ? 0 : 2 ) : ( (y>z) ? 1 : 2 ) );
}

/*-------------------------------------------------------------------------
 *
 *      Function:     Sort3
 *      Description:  three values in ascending order
 *
 *------------------------------------------------------------------------*/
static void Sort3_Swap(real8 *a, real8 *b, int *ia, int *ib)
{
    if(*a > *b)
    {
        real8 tmp  = *a ; *a  = *b ; *b  = tmp ;
        int   itmp = *ia; *ia = *ib; *ib = itmp;
    }
}

void Sort3(real8 values[3], int ind[3])
{
    ind[0]=0; ind[1]=1; ind[2]=2;

    Sort3_Swap(&values[0], &values[1], &ind[0], &ind[1]);
    Sort3_Swap(&values[1], &values[2], &ind[1], &ind[2]);
    Sort3_Swap(&values[0], &values[1], &ind[0], &ind[1]);
}


/**************************************************************************
 *
 *      Function:     Print3
 *      Description:  Print the components of a 3 element vector
 *
 *      Arguments:
 *          msg  text to be displayed with the vector data.
 *          A    3 element vector to be displayed
 *
 *************************************************************************/
void Print3(const char *msg, real8 A[3])
{
        printf("%s = ", msg);
        printf("[%.15e %.15e %.15e];\n", A[0], A[1], A[2]);

        return;
}

/**************************************************************************
 *
 *      Function:     PrintVec
 *      Description:  Print the components of a 3 element vector
 *
 *      Arguments:
 *          msg  text to be displayed with the vector data.
 *          A    3 element vector to be displayed
 *
 *************************************************************************/
void PrintVec(const char *msg, real8 A0, real8 A1, real8 A2, int i)
{
   real8 norm = A0*A0 + A1*A1 + A2*A2;
   norm = sqrt(norm);

   printf("%s = ", msg);

   if (norm > 1e-12)
      if (i == 0)
         printf("[%g %g %g]; %g %g\n", A0/norm, A1/norm, A2/norm, norm*2.85e-10,norm);
      else
         printf("[%g %g %g]; %g %g\n", A0/norm, A1/norm, A2/norm, norm*1e-9,norm);
   else
      printf("[0 0 0];\n");

   return;
}


/**************************************************************************
 *
 *      Function:     Print3x3
 *      Description:  Print the components of a 3 X 3 matrix.
 *
 *      Arguments:
 *          msg  text to be displayed with the matrix data.
 *          A    3 X 3 matrix to be displayed
 *
 *************************************************************************/
void Print3x3(const char *msg, real8 A[3][3])
{
        printf("%s= ...\n", msg);

        printf("[%.15e %.15e %.15e;\n"  , A[0][0], A[0][1], A[0][2]);
        printf(" %.15e %.15e %.15e;\n"  , A[1][0], A[1][1], A[1][2]);
        printf(" %.15e %.15e %.15e;];\n\n", A[2][0], A[2][1], A[2][2]);

        return;
}

/**************************************************************************
 *
 *      Function:     Print2x3
 *      Description:  Print the components of a 3 X 3 matrix.
 *
 *      Arguments:
 *          msg  text to be displayed with the matrix data.
 *          A    3 X 3 matrix to be displayed
 *
 *************************************************************************/
void Print2x3(const char *msg, real8 A[2][3])
{
        printf("%s= ...\n", msg);

        printf("[%.15e %.15e %.15e;\n"      , A[0][0], A[0][1], A[0][2]);
        printf(" %.15e %.15e %.15e;];\n\n"  , A[1][0], A[1][1], A[1][2]);
        return;
}

/**************************************************************************
 *
 *      Function:     Print6to3x3
 *      Description:  Print the components of a 6 tensor into a 3 X 3 matrix.
 *
 *      Arguments:
 *          msg  text to be displayed with the matrix data.
 *          A    6 tensor to be displayed
 *
 *************************************************************************/
void Print6to3x3(const char *msg, real8 A[6])
{
   real8 B[3][3];

        printf("%s= ...\n", msg);

        B[0][0] = A[0];
        B[1][1] = A[1];
        B[2][2] = A[2];
        B[1][2] = A[3];
        B[2][0] = A[4];
        B[0][1] = A[5];
        B[2][1] = B[1][2];
        B[0][2] = B[2][0];
        B[1][0] = B[0][1];



        printf("[%.15e %.15e %.15e;\n"  , B[0][0], B[0][1], B[0][2]);
        printf(" %.15e %.15e %.15e;\n"  , B[1][0], B[1][1], B[1][2]);
        printf(" %.15e %.15e %.15e;];\n\n", B[2][0], B[2][1], B[2][2]);

        return;
}

#ifdef CALCENERGY
void Write_Node_Energy(Home_t *home)
{
   int domID = home->myDomain;
   if (domID == 0) printf("Total energy = %f \n",home->param->TotalEnergy);
}
#endif
/*-------------------------------------------------------------------------
 *
 *      Function:       InitState()
 *
 *      Description:    Initialize the original node's state.
 *
 *-----------------------------------------------------------------------*/
void InitState(NodeState_t *origNodeState)
{
  origNodeState->numNbrs = 0;
  origNodeState->x = 0.0;  origNodeState->y = 0.0;  origNodeState->z = 0.0;
  origNodeState->fx= 0.0;  origNodeState->fy= 0.0;  origNodeState->fz= 0.0;
  origNodeState->vx= 0.0;  origNodeState->vy= 0.0;  origNodeState->vz= 0.0;	
  origNodeState->segState = (SegState_t *) NULL;
}
/*-------------------------------------------------------------------------
 *
 *      Function:       PreserveState()
 *
 *      Description:    Preserve all information pertaining to the specified
 *                      node which will be required to restore the node to
 *                      its original state after we mess with it trying
 *                      to determine if it should be splintered and cross-
 *                      slipped.  NOTE: This includes the force on BOTH ends
 *                      of each segment attached to this node.
 *
 *     Parameters:
 *         IN:  node          Pointer to node info to be preserved
 *
 *         OUT: origNodeState Pointer to structure in which the node (and
 *                            segment) data will be preserved.  This
 *                            structure will contain dynamically allocated
 *                            memory that must be freed by the caller when
 *                            it is no longer needed.
 *-----------------------------------------------------------------------*/
void PreserveState(Home_t *home, Node_t *node,
                   NodeState_t *origNodeState)
{
        int i, numNbrs;

        numNbrs = node->numNbrs;

        origNodeState->origNodeTag.domainID = node->myTag.domainID;
        origNodeState->origNodeTag.index = node->myTag.index;

        origNodeState->x = node->x;
        origNodeState->y = node->y;
        origNodeState->z = node->z;

        origNodeState->fx = node->fX;
        origNodeState->fy = node->fY;
        origNodeState->fz = node->fZ;

        origNodeState->vx = node->vX;
        origNodeState->vy = node->vY;
        origNodeState->vz = node->vZ;

        origNodeState->numNbrs = node->numNbrs;

/*
 *      Allocate an array of structures to hold segment info then loop
 *      over all segments attached to this node and save some of the
 *      segment-specific information.
 */
        origNodeState->segState = (SegState_t *)calloc(1, numNbrs *
                                                          sizeof(SegState_t));

        for (i = 0; i < numNbrs; i++) {
            int        nbrSegID;
            Node_t     *tmpNbr;
            SegState_t *segment = &origNodeState->segState[i];

            segment->nbrTag.domainID = node->nbrTag[i].domainID;
            segment->nbrTag.index = node->nbrTag[i].index;

            segment->burg[X] = node->burgX[i];
            segment->burg[Y] = node->burgY[i];
            segment->burg[Z] = node->burgZ[i];

            segment->planeNorm[X] = node->nx[i];
            segment->planeNorm[Y] = node->ny[i];
            segment->planeNorm[Z] = node->nz[i];

            segment->fNode[X] = node->armfx[i];
            segment->fNode[Y] = node->armfy[i];
            segment->fNode[Z] = node->armfz[i];

/*
 *          Save the force associated with this segment at the neighbor's
 *          end of the segment.
 */
            tmpNbr = GetNeighborNode(home, node, i);
            nbrSegID = GetArmID(tmpNbr, node);

            segment->fNbr[X] = tmpNbr->armfx[nbrSegID];
            segment->fNbr[Y] = tmpNbr->armfy[nbrSegID];
            segment->fNbr[Z] = tmpNbr->armfz[nbrSegID];
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       RestoreState()
 *
 *      Description:    Restore the specified node/segments to the original
 *                      state
 *
 *     Parameters:
 *         IN:  node          Pointer to node info to be restored
 *         OUT: origNodeState Pointer to structure in which the node (and
 *                            segment) data has been preserved.
 *
 *-----------------------------------------------------------------------*/
void  RestoreState(Home_t *home, Node_t *node,
                   NodeState_t *origNodeState)
{
        int i;

        node->x = origNodeState->x;
        node->y = origNodeState->y;
        node->z = origNodeState->z;

        node->fX = origNodeState->fx;
        node->fY = origNodeState->fy;
        node->fZ = origNodeState->fz;

        node->vX = origNodeState->vx;
        node->vY = origNodeState->vy;
        node->vZ = origNodeState->vz;

        node->numNbrs = origNodeState->numNbrs;

        for (i = 0; i < node->numNbrs; i++) {
            int        nbrSegID;
            Node_t     *tmpNbr;
            SegState_t *segment = &origNodeState->segState[i];

            node->nbrTag[i].domainID = segment->nbrTag.domainID;
            node->nbrTag[i].index = segment->nbrTag.index;

            node->burgX[i] = segment->burg[X];
            node->burgY[i] = segment->burg[Y];
            node->burgZ[i] = segment->burg[Z];

            node->nx[i] = segment->planeNorm[X];
            node->ny[i] = segment->planeNorm[Y];
            node->nz[i] = segment->planeNorm[Z];

            node->armfx[i] = segment->fNode[X];
            node->armfy[i] = segment->fNode[Y];
            node->armfz[i] = segment->fNode[Z];

/*
 *          Restore the force and burgers vector associated with this
 *          segment at the neighbor's end of the segment.
 */
            tmpNbr = GetNeighborNode(home, node, i);
            nbrSegID = GetArmID(tmpNbr, node);

            tmpNbr->armfx[nbrSegID] = segment->fNbr[X];
            tmpNbr->armfy[nbrSegID] = segment->fNbr[Y];
            tmpNbr->armfz[nbrSegID] = segment->fNbr[Z];

            tmpNbr->burgX[nbrSegID] = -segment->burg[X];
            tmpNbr->burgY[nbrSegID] = -segment->burg[Y];
            tmpNbr->burgZ[nbrSegID] = -segment->burg[Z];
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       DiscardState()
 *
 *      Description:    Frees any dynamically allocated memory used to
 *                      preserve the original node's state.
 *
 *     Parameters:
 *         IN/OUT: origNodeState Pointer to structure in which the node (and
 *                               segment) data was preserved.  Dynamically
 *                               alloctated memory is freed and pointers
 *                               NULL'ed before return to caller.
 *
 *-----------------------------------------------------------------------*/
void DiscardState(NodeState_t *origNodeState)
{

        if (origNodeState->segState != (SegState_t *)NULL) {
            free(origNodeState->segState);
            origNodeState->segState = (SegState_t *)NULL;
        }

        return;
}


#if ESHELBY
void PrintInclusion(EInclusion_t       *inclusion)
{
  printf("Inclusion ID = %d cell ID =%d\n",inclusion->id,inclusion->cellID);
  Print3("Position",inclusion->position);

}
#endif
