/****************************************************************************
 *
 *      Function:    InitCellNeighbors
 *      Description: For each native cell in domain, find all its neighbor
 *                   cells. If a neighbor cell has already been encountered
 *                   and allocated, just add its cell index to the native
 *                   cell's list of neighbors. If this is the first time,
 *                   allocate it, add it to the cell hash table, and add
 *                   it to home->cellList as well as to the native node's
 *                   neighbor list.
 *
 *                   If periodic boundaries are turned on and a neighbor
 *                   cell lies outside the base problem space, a cell is
 *                   allocated for the periodic image cell, but the cell
 *                   only contains offsets to apply to the corresponding
 *                   base cell, and a pointer to the base cell, which is
 *                   also allocated even if it isn't used by domain except
 *                   in shifted form. Only the base cell is added to
 *                   home->cellList
 *
 ****************************************************************************/
#include <float.h>

#include "mpi_portability.h"

#include "Home.h"


void InitCellNeighbors(Home_t *home)
{
    int     iCell, ntvIdx, iX, jY, kZ, iNbr, jNbr, kNbr, nbrIdx, isPeriodic;
    int     iBase, jBase, kBase, baseIdx;
    int     nXcells, nYcells, nZcells;
    real8   xShift, yShift, zShift;
    real8   sideX, sideY, sideZ;
    Cell_t  *ntvCell, *nbrCell, *baseCell;
    Param_t *param;

    param = home->param;
   
    nXcells = param->nXcells;
    nYcells = param->nYcells;
    nZcells = param->nZcells;

    sideX = (param->maxSideX - param->minSideX);
    sideY = (param->maxSideY - param->minSideY);
    sideZ = (param->maxSideZ - param->minSideZ);

/*
 *  Loop over all cells native to the domain
 */
    for (iCell = 0; iCell < home->nativeCellCount; iCell++) {
        int updateNbrList;

        ntvIdx = home->cellList[iCell];
        ntvCell = LookupCell(home, ntvIdx);

/*
 *      Since the list of cell neighboring a givebn cell does not 
 *      change throughout the simulation, we only have to allocate
 *      and initialize the list once after which it will be preserved
 *      from timestep to timestep.
 */
        if (ntvCell->nbrList == (int *)NULL) {
            ntvCell->nbrList = (int *) malloc(26 * sizeof(int));
            ntvCell->nbrCount = 0;
            updateNbrList = 1;
        } else {
            updateNbrList = 0;
        }

/*
 *      Include all cells surrounding the native cell.  Base cells
 *      are in the range [1,nXcells] X [1,nYcells] X [1,nZcells].
 *      Cells outside this range (i.e. with any cell index 0 or
 *      n[X|Y|Z]cells) are "periodic" cells.
 */
        iX = ntvCell->xIndex;
        jY = ntvCell->yIndex;
        kZ = ntvCell->zIndex;

        for (iNbr = iX-1; iNbr <= iX+1; iNbr++) {
            for (jNbr = jY-1; jNbr <= jY+1; jNbr++) {
                for (kNbr = kZ-1; kNbr <= kZ+1; kNbr++) {

/*
 *                  Skip the current native cell itself
 */
                    if (iNbr == iX && jNbr == jY && kNbr == kZ) continue;

/*
 *                  In the following, check for cell indices outside
 *                  the base cell region. If so, and if problem is
 *                  periodic in that dimension, determine how much
 *                  the base cell must be shifted in the periodic
 *                  image, what the corresponding base cell index
 *                  is, and flag the cell as periodic. If problem
 *                  not periodic in that dimension, just skip - there
 *                  is no neighbor in that direction.
 */
                    isPeriodic = 0;        /* assume not periodic */

                    xShift = 0.0;
                    yShift = 0.0;
                    zShift = 0.0;

                    iBase = iNbr;
                    jBase = jNbr;
                    kBase = kNbr;

                    if (iNbr == 0) {    /* periodic cell */
                        if (param->xBoundType == Periodic) {
    
                            /* i-index of corresponding base cell */
                            iBase = nXcells;

                            /* correction to base X coordinate */
                            xShift = -sideX;

                            /* flag cell as "periodic" */
                            isPeriodic = 1;

                        } else continue; /* free boundary; no neighbor */
                    }

                    if (iNbr > nXcells) {    /* periodic cell */
                        if (param->xBoundType == Periodic) {

                            iBase = 1;
                            xShift = sideX;
                            isPeriodic = 1;

                        } else continue;
                    }

                    if (jNbr == 0) {    /* periodic cell */
                        if (param->yBoundType == Periodic) {

                            /* j-index of corresponding base cell */
                            jBase = nYcells;

                            /* correction to base Y coordinate */
                            yShift = -sideY;

                            /* flag cell as "periodic" */
                            isPeriodic = 1;

                        } else continue; /* free boundary; no neighbor */
                    }

                    if (jNbr > nYcells) {    /* periodic cell */
                        if (param->yBoundType == Periodic) {

                            jBase = 1;
                            yShift = sideY;
                            isPeriodic = 1;

                        } else continue;
                    }

                    if (kNbr == 0) {    /* periodic cell */
                        if (param->zBoundType == Periodic) {

                            /* k-index of corres' base cell */
                            kBase = nZcells;

                            /* correction to base Z coordinate */
                            zShift = -sideZ;

                            /* flag cell as "periodic" */
                            isPeriodic = 1;

                        } else continue; /* free boundary; no neighbor */
                    }

                    if (kNbr > nZcells) {    /* periodic cell */
                        if (param->zBoundType == Periodic) {

                            kBase = 1;
                            zShift = sideZ;
                            isPeriodic = 1;

                        } else continue;
                    }

/*
 *                  If neighbor cell is already allocated, add its
 *                  index to the native cell's neighbor list, mark
 *                  the apporpriate cells 'in use' and update the
 *                  main cell list if necessary.
 */
                    nbrIdx = EncodeCellIdx(home, iNbr, jNbr, kNbr);
                    nbrCell = LookupCell(home, nbrIdx);

                    if (nbrCell != (Cell_t *)NULL) {

                        if (updateNbrList) {
                            ntvCell->nbrList[ntvCell->nbrCount++] = nbrIdx;
                        }

/*
 *                      If the cell is a periodic cell, flag both it and the
 *                      corresponding base cell as in use.  Also add the 
 *                      base cell to the main cell list if it had not been
 *                      flagged in use before now.
 */
                        if (isPeriodic) {

                            nbrCell->inUse = 1;

                            baseIdx  = EncodeCellIdx(home, iBase, jBase, kBase);
                            baseCell = LookupCell(home, baseIdx);

                            if (baseCell->inUse == 0) {
                                home->cellList[home->cellCount++] = baseIdx;
                            }

                            baseCell->inUse = 1;

                        } else {
/*
 *                          Cell is a base cell, so if it has not yet been
 *                          been flagged as in use, flag it and add it to
 *                          the main cell list.
 */
                            if (nbrCell->inUse == 0) {
                                nbrCell->inUse = 1;
                                home->cellList[home->cellCount++] = nbrIdx;
                            }
                        }

                        continue;
                    }

/*
 *                  If current neighbor cell is a periodic cell,
 *                  see if base period cell has been alloc'd. If 
 *                  not, alloc now. Then alloc periodic cell and
 *                  set its shifts and a pointer to its base cell.
 *                  If neighbor is not a periodic cell, just alloc
 *                  it. In any case, add the new neighbor to the
 *                  native cell's neighbor list.
 */
                    if (isPeriodic == 0) {

/*
 *                      This is a base cell, not a periodic cell.  If it
 *                      has not yet been allocated allocate it, add it to
 *                      the cell has table, add it to the domain cell list
 *                      and to the neighbor list for this native cell
 */
                        nbrCell = LookupCell(home, nbrIdx);

                        if (nbrCell == (Cell_t *)NULL) {

/*
 *                          Cell struct does not yet exist, so add it to the
 *                          cell table.  Some of the values associated with a
 *                          cell are static and will not change over time.
 *                          These values will be set now and never changed.
 *
 *
 *                          This cell is a base cell (not a periodic image) 
 *                          so initialize periodic shifts to 0 and index of
 *                          the base cell to -1 (i.e. none)
 */
                            nbrCell = AddCellToTable(home, nbrIdx);

                            nbrCell->baseIdx = -1;
   
                            nbrCell->xShift = xShift;
                            nbrCell->yShift = yShift;
                            nbrCell->zShift = zShift;

/*
 *                          There are a few places we need to find the center
 *                          of a cell, so do the caluclation here and save
 *                          the results for later on.
 */
                            FindCellCenter(param, iNbr-1, jNbr-1, kNbr-1, 2,
                                           &nbrCell->center[X],
                                           &nbrCell->center[Y],
                                           &nbrCell->center[Z]);

/*
 *                          For base cells we need to save the per-
 *                          dimension cell indices so we don't have to make
 *                          repeated calls to DecodeCellIdx() later to figure
 *                          them out.
 */
                            nbrCell->xIndex = iNbr;
                            nbrCell->yIndex = jNbr;
                            nbrCell->zIndex = kNbr;

/*
 *                          If threading has been enabled, we'll also need 
 *                          to create a lock by which a thread can gain
 *                          exclusive access to the cell-specific data.
 */
                            INIT_LOCK(&nbrCell->cellLock);
                        }

                        nbrCell->inUse = 1;
                        home->cellList[home->cellCount++] = nbrIdx;

                        if (updateNbrList) {
                            ntvCell->nbrList[ntvCell->nbrCount++] = nbrIdx;
                        }

                    } else {    /*  periodic cell  */

/*
 *                      If the base cell for periodic neighbor has not been
 *                      allocated, allocate it now.
 */
                        baseIdx  = EncodeCellIdx(home, iBase, jBase, kBase);
                        baseCell = LookupCell(home, baseIdx);

                        if (baseCell == (Cell_t *)NULL) {

/*
 *                          Cell struct does not yet exist, so add it to the
 *                          cell table.  Some of the values associated with a
 *                          cell are static and will not change over time.
 *                          These values will be set now and never changed.
 *
 *
 *                          Since this is not a periodic cell, initialize
 *                          periodic shifts to 0 and index of the base cell
 *                          to -1 (i.e. none)
 */
                            baseCell = AddCellToTable(home, baseIdx);

                            baseCell->xShift = 0.0;
                            baseCell->yShift = 0.0;
                            baseCell->zShift = 0.0;

                            baseCell->nbrList = 0;
                            baseCell->baseIdx = -1;

/*
 *                          For base cells we need to save the per-
 *                          dimension cell indices so we don't have to make
 *                          repeated calls to DecodeCellIdx() later to
 *                          figure them out.
 */
                            baseCell->xIndex = iBase;
                            baseCell->yIndex = jBase;
                            baseCell->zIndex = kBase;

/*
 *                          There are a few places we need to find the
 *                          center of a cell, so do the calculation
 *                          here and save the results for later on.
 *
 *                          NOTE: The center of the cell is the center
 *                          base cell, not the ghost cell.
 */
                            FindCellCenter(param,
                                           iBase-1, jBase-1, kBase-1, 2,
                                           &baseCell->center[X],
                                           &baseCell->center[Y],
                                           &baseCell->center[Z]);

/*
 *                          If threading has been enabled, we'll also need 
 *                          to create a lock by which a thread can gain
 *                          exclusive access to the cell-specific data.
 */
                            INIT_LOCK(&baseCell->cellLock);

                        }

/*
 *                      If the base cell has not yet been marked as in use 
 *                      this timestep, do so and add it to the main cell list
 */
                        if (baseCell->inUse == 0) {
                            baseCell->inUse = 1;
                            home->cellList[home->cellCount++] = baseIdx;
                        }

/*
 *                      The neighboring cell is a periodic cell that has not
 *                      yet been allocated, so do that now.  NOTE: periodic
 *                      cells are not added to the domain cell list, nor do
 *                      they get a "lock" for use when threading is enabled.
 */
                        nbrCell = AddCellToTable(home, nbrIdx);

                        nbrCell->xShift = xShift;
                        nbrCell->yShift = yShift;
                        nbrCell->zShift = zShift;

                        nbrCell->baseIdx = baseIdx;

/*
 *                      For periodic(ghost) cells we *shouldn't* need to
 *                      look up the per-dimension cell indices, so initialize
 *                      them to invalid indices.
 */
                        nbrCell->xIndex = -1;
                        nbrCell->yIndex = -1;
                        nbrCell->zIndex = -1;

/*
 *                      We also shouldn't need to look up the center
 *                      of the cell, so initialize the center to 
 *                      something outrageous...
 */
                        nbrCell->center[X] = DBL_MAX;
                        nbrCell->center[Y] = DBL_MAX;
                        nbrCell->center[Z] = DBL_MAX;

                        if (updateNbrList) {
                            ntvCell->nbrList[ntvCell->nbrCount++] = nbrIdx;
                        }

                    }  /* end if (isPeriodic == 0)  */

                }  /* end for (kNbr ...) */

            }  /* end for (jNbr ...)  */

        }  /* end for (iNbr ...) */

    }  /* end for (iCell ...) */

    return;
}
