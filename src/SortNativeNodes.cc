/**************************************************************************
 *
 *  Function    : SortNativeNodes 
 *  Description : Queue each native node onto the cell it currently falls into.
 *
 *    7/31/01  Replaced sideX with minSideX, maxSideX, etc     t.p.
 *
 **************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "Node.h"
#include "Cell.h"
#include "Util.h"

void AssignNodeToCell(Home_t *home, Node_t *node)
{
        Param_t *param;
        int iCell, jCell, kCell, cellIdx;
        real8 probXmin, probYmin, probZmin;
        real8 cellXsize, cellYsize, cellZsize;
        Cell_t *cell;

        if (node == (Node_t *)NULL) return;

/*
 *      set the lower limit of the base area (excluding possible periodic cells)
 *      and the size of each cell
 */
        param = home->param;
        probXmin = param->minSideX;
        probYmin = param->minSideY;
        probZmin = param->minSideZ;

        cellXsize = (param->maxSideX - param->minSideX) / param->nXcells;
        cellYsize = (param->maxSideY - param->minSideY) / param->nYcells;
        cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells;

/*
 *      put the node on its proper cell. If the index exceeds this domains
 *      range of native cells, put in the nearest native cell
 */
        iCell = (int)((node->x - probXmin) / cellXsize);
        if (iCell < param->iCellNatMin) iCell = param->iCellNatMin;
        if (iCell > param->iCellNatMax) iCell = param->iCellNatMax;
        iCell++;   /* compensate for periodic cells */

        jCell = (int)((node->y - probYmin) / cellYsize);
        if (jCell < param->jCellNatMin) jCell = param->jCellNatMin;
        if (jCell > param->jCellNatMax) jCell = param->jCellNatMax;
        jCell++;   /* compensate for periodic cells */

        kCell = (int)((node->z - probZmin) / cellZsize);
        if (kCell < param->kCellNatMin) kCell = param->kCellNatMin;
        if (kCell > param->kCellNatMax) kCell = param->kCellNatMax;
        kCell++;   /* compensate for periodic cells */

        cellIdx = EncodeCellIdx (home, iCell, jCell, kCell);
        cell = LookupCell(home, cellIdx);
        node->nextInCell = cell->nodeQ;
        cell->nodeQ = node;
        cell->nodeCount++;

        node->cellIdx = cellIdx;

        return;
}


void SortNativeNodes (Home_t *home)
{
        TimerStart(home, SORT_NATIVE_NODES);

        Param_t *param = home->param;

/*
 *      Set the lower limit of the base area (excluding possible
 *      periodic cells) and the size of each cell
 */
        real8 probXmin  = param->minSideX;
        real8 probYmin  = param->minSideY;
        real8 probZmin  = param->minSideZ;

        real8 cellXsize = (param->maxSideX - param->minSideX) / param->nXcells;
        real8 cellYsize = (param->maxSideY - param->minSideY) / param->nYcells;
        real8 cellZsize = (param->maxSideZ - param->minSideZ) / param->nZcells;

/*
 *      Loop through cells and reinitialize their node and eshelby
 *      inclusion queues to empty.
 */
        for (int i=0; (i<home->cellCount); i++) 
        {
            Cell_t *cell = LookupCell(home, home->cellList[i]);
                    cell->nodeQ          =  0;
                    cell->nodeCount      =  0;
                    cell->inclusionQ     = -1;
                    cell->inclusionCount =  0;
        }

/*
 *      Loop thru active nodes, putting them in their proper cell.
 *      If the index exceeds this domains range of native cells, put
 *      in the nearest native cell
 */
        for (int i=0; (i<home->newNodeKeyPtr); i++) 
        {
            Node_t *node = home->nodeKeys[i];

            if (!node) continue;

            int iCell = (int)((node->x - probXmin) / cellXsize);

            if (iCell < param->iCellNatMin) { iCell = param->iCellNatMin; } 
            if (iCell > param->iCellNatMax) { iCell = param->iCellNatMax; }

            iCell++;   /* compensate for periodic cells */

            int jCell = (int)((node->y - probYmin) / cellYsize);

            if (jCell < param->jCellNatMin) { jCell = param->jCellNatMin; } 
            if (jCell > param->jCellNatMax) { jCell = param->jCellNatMax; }

            jCell++;   /* compensate for periodic cells */

            int kCell = (int)((node->z - probZmin) / cellZsize);

            if (kCell < param->kCellNatMin) { kCell = param->kCellNatMin; } 
            if (kCell > param->kCellNatMax) { kCell = param->kCellNatMax; }

            kCell++;   /* compensate for periodic cells */

            int     cellIdx = EncodeCellIdx (home, iCell, jCell, kCell);
            Cell_t *cell    = LookupCell(home, cellIdx);

            node->nextInCell = cell->nodeQ;
            cell->nodeQ = node;
            cell->nodeCount++;

            node->cellIdx = cellIdx;
            node->native  = 1;
        }

#ifdef ESHELBY
/*
 *      Loop through all the native Eshelby inclusions and
 *      place them on the correct cell queues
 */
        for (int i=0; i<home->locInclusionCount; i++) 
        {
            EInclusion_t *inclusion = &home->eshelbyInclusions[i];

            int     cellID = inclusion->cellID;
            Cell_t *cell   = LookupCell(home, cellID);

            inclusion->nextInCell = cell->inclusionQ;
            cell->inclusionQ = i;
            cell->inclusionCount++;
        }
#endif

        TimerStop(home, SORT_NATIVE_NODES);

        return;
}
