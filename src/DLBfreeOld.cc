/***************************************************************************
 *
 *  function    : DLBfreeOld
 *  description : Remove the cell and remoteDomain structs, in preparation
 *                for reinitialization
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Decomp.h"


void DLBfreeOld(Home_t *home)
{
        int i, iRemDom, totRemDomCount;
        RemoteDomain_t *remDom;


/*
 *      To save time, we save the cell structures across timesteps, only
 *      allocating new ones as needed.  Some of the cell data, though,
 *      is dynamic in nature, so we need to clear those particular
 *      data items from all non-periodic cell structures each time
 *      the domain bounds shift...
 */
        for (i = 0; i < CELL_HASH_TABLE_SIZE; i++) {
            Cell_t *cell;

/*
 *          Check every cell in the linked list for this hash table entry
 */
            cell = home->cellTable[i];

            while (cell != (Cell_t *)NULL) {

/*
 *              Since cell structs persist across timesteps, we need to
 *              set a flag indicating if the cell is in use during a 
 *              cycle (needed when initializing the cell neighbor lists).
 *              Reset the flag to zero now.
 */
                cell->inUse = 0;

/*
 *              Only the base cells (non-periodic) have dynamic data that
 *              needs to be cleared.  If baseIdx >= zero, cell is periodic.
 */
                if (cell->baseIdx < 0) {

                    if (cell->domains) {
                        free(cell->domains);
                    }

                    cell->domains = (int *)NULL;
                    cell->domCount = 0;

                    cell->nodeQ = (Node_t *)NULL;
                    cell->nodeCount = 0;

                    cell->inclusionQ = -1;
                    cell->inclusionCount = 0;
                }

                cell = cell->next;
            }
        }

        if (home->cellList != (int *)NULL) {
            free(home->cellList);
        }

        home->cellList = (int *)NULL;

/*
 *      Clean out remote domain structures (including some MPI stuff)
 *      Need to do this for both primary and secondary remote domains.
 */
        totRemDomCount = home->remoteDomainCount +
                         home->secondaryRemoteDomainCount;

        for (i = 0; i < totRemDomCount; i++) {

            iRemDom = home->remoteDomains[i];
            remDom = home->remoteDomainKeys[iRemDom];

            if (!remDom) {
                printf ("DLBfreeOld: unexpected result!\n");
                continue;
            }

/*
 *          Secondary remote domains will not have exported cell lists
 */
            if (remDom->expCells != (int *)NULL) {
                free(remDom->expCells);
            }

            free(remDom->nodeKeys);
            free(remDom);
        }

        free(home->remoteDomains);
        free(home->remoteDomainKeys);

/*
 *      All cell queues have been wiped out, nodes have been moved
 *      and need to be assigned to new cells, but until that's done,
 *      we need to wipe out the cell index associated with each node.
 *      Wipe out cell2 references while we're at it.
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            Node_t *node;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            node->cellIdx      = -1;
            node->cell2Idx     = -1;
            node->cell2QentIdx = -1;
            node->nextInCell   = (Node_t *)NULL;
        }

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->ghostNodeCount; i++) {
            Node_t *node;

            node = home->ghostNodeList[i];

            node->cellIdx      = -1;
            node->cell2Idx     = -1;
            node->cell2QentIdx = -1;
            node->nextInCell   = (Node_t *)NULL;
        }

        return;
}
