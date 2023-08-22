/****************************************************************************
 *
 *      Function:     InitCellNatives
 *      Description:  Find all the cells that intersect this domain.
 *                    Allocate and initialize a cell struct for each, add
 *                    the cell struct to the cell hash table, and
 *                    and add the decoded index to the domain's cellList
 *
 ****************************************************************************/

#include "mpi_portability.h"

#include "Home.h"

void InitCellNatives(Home_t *home)
{
        int     maxExtendedCells, cIndex;
        int     xIndex, yIndex, zIndex;
        int     xCellMin, yCellMin, zCellMin;
        int     xCellMax, yCellMax, zCellMax;
        real8   xMin, yMin, zMin;
        real8   xMax, yMax, zMax;
        real8   xCellSize, yCellSize, zCellSize;
        real8   cellMax, cellMin;
        Param_t *param;

        param = home->param;

/*
 *      Determine cell dimensions
 */
        xCellSize = (param->maxSideX - param->minSideX) / param->nXcells;
        yCellSize = (param->maxSideY - param->minSideY) / param->nYcells;
        zCellSize = (param->maxSideZ - param->minSideZ) / param->nZcells;

/*
 *      Make a temporary copy of the local domain boundaries.
 */
        xMin = home->domXmin;
        xMax = home->domXmax;

        yMin = home->domYmin;
        yMax = home->domYmax;

        zMin = home->domZmin;
        zMax = home->domZmax;

/*
 *      Find the min and max cell indices (in each dimension) of cells
 *      intersecting this domain
 *
 *      FIX ME!  We should be able to replace the loops below with
 *      direct calculations of the cell indices.
 */
        xCellMin = 0;
        yCellMin = 0;
        zCellMin = 0;

        xCellMax = param->nXcells-1;
        yCellMax = param->nYcells-1;
        zCellMax = param->nZcells-1;

/*
 *      Get min/max cell indices in X
 */
        for (cIndex = 0; cIndex <= param->nXcells; cIndex++) {

            cellMax = rint((param->minSideX) + (cIndex * xCellSize));
            cellMin = rint((param->minSideX) + (cIndex-1) * xCellSize);

            if ((xMin >= cellMin) && (xMin < cellMax)) xCellMin = cIndex;
            if ((xMax <= cellMax) && (xMax > cellMin)) xCellMax = cIndex;
        }

/*
 *      Get min/max cell indices in Y
 */
        for (cIndex = 0; cIndex <= param->nYcells; cIndex++) {
    
            cellMax = rint((param->minSideY) + (cIndex * yCellSize));
            cellMin = rint((param->minSideY) + (cIndex-1) * yCellSize);

            if ((yMin >= cellMin) && (yMin < cellMax)) yCellMin = cIndex;
            if ((yMax <= cellMax) && (yMax > cellMin)) yCellMax = cIndex;
    }

/*
 *      Get min/max cell indices in Z
 */
        for (cIndex = 0; cIndex <= param->nZcells; cIndex++) {

            cellMax = rint((param->minSideZ) + (cIndex * zCellSize));
            cellMin = rint((param->minSideZ) + (cIndex-1) * zCellSize);

            if ((zMin >= cellMin) && (zMin < cellMax)) zCellMin = cIndex;
            if ((zMax <= cellMax) && (zMax > cellMin)) zCellMax = cIndex;
    }


/*
 *      Save the minumum and maximum cells in each dimension in param, for
 *      later use in SortNativeNodes. Save the "natural" (i.e. pre-shifted)
 *      indices. They are easier to work with in [i,j,k] space.  Shift cell
 *      indices by one to make room for periodic images at index 0
 */
        param->iCellNatMin = xCellMin-1;
        param->iCellNatMax = xCellMax-1;

        param->jCellNatMin = yCellMin-1;
        param->jCellNatMax = yCellMax-1;

        param->kCellNatMin = zCellMin-1;
        param->kCellNatMax = zCellMax-1;


/*
 *      Calculate the max number of cells either native to (i.e. intersecting) 
 *      this domain, or neighbor to a native cell.
 */
        maxExtendedCells = (xCellMax-xCellMin+3) *
                           (yCellMax-yCellMin+3) *
                           (zCellMax-zCellMin+3);

/*
 *      The cellList is a list of the encoded indices of each native and
 *      neighbor cell for this domain. The encoded index of the cell
 *      (X, Y, Z) is Z + (nZCells * Y) + (nZcells * nYcells * X).
 *      The native cells are always listed first in cellList; neighbor
 *      cells will be added later.
 */
        home->cellList = (int *) malloc(maxExtendedCells * sizeof(int));
        home->cellCount = 0;

/*
 *      Allocate a Cell_t struct for each native cell, add it to the
 *      cell hash table, and add its encoded index into cellList.
 */
        for (xIndex = xCellMin; xIndex <= xCellMax; xIndex++) {
            for (yIndex = yCellMin; yIndex <= yCellMax; yIndex++) {
                for (zIndex = zCellMin; zIndex <= zCellMax; zIndex++) {
                    Cell_t *pCell;

                    cIndex = EncodeCellIdx(home, xIndex, yIndex, zIndex);

                    home->cellList[home->cellCount++] = cIndex;

/*
 *                  Check if the cell has already been allocated.  If not,
 *                  allocate it and add it to the call table.
 */
                    pCell = LookupCell(home, cIndex);

                    if (pCell == (Cell_t *)NULL) {

                        pCell = AddCellToTable(home, cIndex);

/*
 *                      Cell struct does not yet exist, so add it to the
 *                      cell table.  Some of the values associated with a
 *                      cell are static and will not change over time.
 *                      These values will be set now and never changed.
 *
 *
 *                      Since native cells cannot be periodic images,
 *                      initialize periodic shifts to 0 and index of the
 *                      base cell to -1 (i.e. none)
 */
                        pCell->baseIdx = -1;

                        pCell->xShift = 0.0;
                        pCell->yShift = 0.0;
                        pCell->zShift = 0.0;

/*
 *                      There are a few places we need to find the center
 *                      of a cell, so do the calculation here and save the
 *                      results for later on.
 */
                        FindCellCenter(param, xIndex-1, yIndex-1, zIndex-1, 2,
                                       &pCell->center[X], &pCell->center[Y],
                                       &pCell->center[Z]);

/*
 *                      Also save the per-dimension indices of the cell so
 *                      that we don't have to repeatedly call DecodeCellIdx()
 *                      later on.
 */
                        pCell->xIndex = xIndex;
                        pCell->yIndex = yIndex;
                        pCell->zIndex = zIndex;

/*
 *                      If threading has been enabled, we'll also need 
 *                      to create a lock by which a thread can gain exclusive
 *                      access to the cell-specific data.
 */
                        INIT_LOCK(&pCell->cellLock);

                    }

/*
 *                  Flag the cell structure as 'in use'
 */
                    pCell->inUse = 1;
                }
            }
        }

        home->nativeCellCount = home->cellCount;

/*
 *      Last thing to be done is (re)initialize some of the
 *      FM cell layers and associated info.
 */
#ifdef ESHELBY
        if ((param->fmEnabled) || (param->eshelbyfmEnabled)) { FMInit(home); }
#else
        if  (param->fmEnabled) { FMInit(home); }
#endif

        return;
}
