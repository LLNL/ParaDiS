/**************************************************************************
 *
 *      Module:       InitRemoteDomains.c
 *      Description:  Contains functions needed for determining which
 *                    remote domains this domain talks to (i.e. all
 *                    domains intersecting cells intersected by this
 *                    domain or in cells neighboring cells intersected by
 *                    this domain).  For each remote domain, a list is
 *                    built indicating which native cells are to be sent
 *                    to the remote domain during communications.
 *
 *      Last Modified: 03/22/2008 gh - Added init of
 *                                     secondaryRemoteDomainCount.
 *
 **************************************************************************/

#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Cell.h"
#include "RemoteDomain.h"

#define NUM_EXPCELL_INC  50

void InitRemoteDomains(Home_t *home)
{
        int            rDomMax, i, remDomCount, myDom, iCell, cellIdx;
        int            iDom, domIdx, natCellIdx, nCell, nbrCellIdx;
        int            found, iExpCell, newCount;
        int            *remDomList;
        RemoteDomain_t *remDom;
        RemoteDomain_t **rDomKeys;
        Cell_t         *cell, *natCell, *nbrCell;

/*
 *      Alloc and initialize the array of remote domain keys; 1 pointer for
 *      each potential remote domain.
 */

        rDomMax = home->param->nXdoms *
                  home->param->nYdoms *
                  home->param->nZdoms;

        rDomKeys = (RemoteDomain_t **)malloc(rDomMax*sizeof(RemoteDomain_t *));
        
        for (i = 0; i < rDomMax; i++) {
            rDomKeys[i] = 0;
        }
        
        remDomList = (int *) malloc(rDomMax * sizeof(int));
        remDomCount = 0;
        myDom = home->myDomain;
        
/*
 *      Find the domains that intersect any of the cells in home->cellList
 */
        for (iCell = 0; iCell < home->cellCount; iCell++) {
        
            cellIdx = home->cellList[iCell];
            cell = LookupCell(home, cellIdx);
        
            for (iDom = 0; iDom < cell->domCount; iDom++) {
/*
 *              Don't count the local domain as remote.
 */
                if ((domIdx = cell->domains[iDom]) == myDom) {
                    continue;
                }
        
/*
 *              If this remote domain has not yet been encountered,
 *              allocate a structure for it, put its address in rDomKeys,
 *              and its index in remDomList.  Also, initialize count of
 *              exported cells to 0, and allocate the out migrator
 *              list and list size
 */
                if (rDomKeys[domIdx] == 0) {

                    remDom = (RemoteDomain_t *)malloc(sizeof(RemoteDomain_t));
                    rDomKeys[domIdx] = remDom;

                    remDomList[remDomCount++] = domIdx;
                    remDom->domainIdx = domIdx;
                    remDom->numExpCells = 0;
                    remDom->expCells = (int *)NULL;
                    remDom->maxTagIndex = 0;
                }
        
            } /* end for (iDom = 0; ...) */
        }  /* end for (iCell = 0; ...) */
        

/*
 *      Loop through the native cells. For each native cell, loop through
 *      all its neighbor cells. For each remote domain that intersects any
 *      of these neighbor cells (or the native cell itself) add the native
 *      cell to the list of cells exported to that remote domain.
 */
        for ( iCell = 0; iCell < home->nativeCellCount; iCell++) {
        
            natCellIdx = home->cellList[iCell];
            natCell = LookupCell(home, natCellIdx);
        
            for (nCell = 0; nCell <= natCell->nbrCount; nCell++) {
        
/*
 *              Splice in native cell itself
 */
                if (nCell == natCell->nbrCount) {
                    nbrCellIdx = natCellIdx;
                    nbrCell = natCell;
                } else {
                    nbrCellIdx = natCell->nbrList[nCell];
                    nbrCell = LookupCell(home, nbrCellIdx);
                }
        
/*
 *              If the neighbor cell is a periodic image cell, use
 *              its corresponding base cell instead
 */
                if (nbrCell->baseIdx != -1) {
                    nbrCellIdx = nbrCell->baseIdx;
                    nbrCell = LookupCell(home, nbrCellIdx);
                }
        
/*
 *              For each intersecting domain of the current nbrCell
 *              (except myDom), if native cell is not already in that
 *              domain's export list, put it in there now.
 */
                for (iDom = 0; iDom < nbrCell->domCount; iDom++) {
        
                    domIdx = nbrCell->domains[iDom];
                    if (domIdx == myDom) continue;   /* don't count self */
                    remDom = rDomKeys[domIdx];
        
/*
 *                  Check remote domain's current list of export cells
 *                  for native cell before adding it.
 */
                    found = 0;

                    for (iExpCell=0; iExpCell<remDom->numExpCells; iExpCell++) {
        
                        if (remDom->expCells[iExpCell] == natCellIdx) {
                            found = 1;
                            break;
                        }
        
                    }  /* end for (iExpCell = 0 ...) */
        
                    if (found == 0) {
/*
 *                      If the export cell array has not yet been allocated
 *                      (or is too small) reallocate it with a small number
 *                      of spare entries.
 */
                        if ((remDom->numExpCells % NUM_EXPCELL_INC) == 0) {
                            newCount = remDom->numExpCells + NUM_EXPCELL_INC;
                            remDom->expCells = (int *)realloc(remDom->expCells,
                                               newCount * sizeof(int));
                        }

                        remDom->expCells[remDom->numExpCells++] = natCellIdx;
                    }
        
                }  /* end for (iDom = 0 ...) */
            }  /* end for (nCell = 0 ...) */
        }  /* end for (iCell = 0 ...) */
        
/*
 *      Trim back the exported cell lists in each remote to domain to
 *      just the number actually used
 */
        for (iDom = 0; iDom < remDomCount; iDom++) {
            domIdx = remDomList[iDom];
            remDom = rDomKeys[domIdx];
            remDom->expCells = (int *)realloc(remDom->expCells, 
                               remDom->numExpCells * sizeof(int));
        }
        
/*
 *      Move the remote domain structs into home
 */
        home->remoteDomains = (int *)realloc(remDomList,
                              remDomCount*sizeof(int));

        home->remoteDomainCount = remDomCount;
        home->remoteDomainKeys = rDomKeys;

/*
 *      Initially, there are no secondary remote domains.  Those
 *      will be added as necessary when sending the secondary
 *      ghost nodes.
 */
        home->secondaryRemoteDomainCount = 0;
        

#ifdef PARALLEL
/*
 *      Allocate MPI structures in home to handle asynchronous msg's
 *      if they have not yet been allocated.
 */
        if (home->inRequests == ((MPI_Request *)NULL)) {
            home->inRequests  = (MPI_Request *) malloc(home->numDomains *
                                                       sizeof(MPI_Request));
        }

        if (home->outRequests == ((MPI_Request *)NULL)) {
            home->outRequests = (MPI_Request *) malloc(home->numDomains *
                                                       sizeof(MPI_Request));
        }

        if (home->inStatus == ((MPI_Status *)NULL)) {
            home->inStatus    = (MPI_Status *) malloc(home->numDomains *
                                                      sizeof(MPI_Status));
        }

        if (home->outStatus == ((MPI_Status *)NULL)) {
            home->outStatus   = (MPI_Status *) malloc(home->numDomains *
                                                      sizeof(MPI_Status));
        }
#endif
        
        return;
}
