/*---------------------------------------------------------------------------
 *
 *      Function:     InitCellDomains
 *      Description:  For each cell in home->cellList, determine which
 *                    domains intersect the cell and save the list in
 *                    the corresponding cell->domains list.
 *
 *--------------------------------------------------------------------------*/

#include "mpi_portability.h"

#include "Home.h"
#include "Cell.h"
#include "Decomp.h"


void InitCellDomains(Home_t *home)
{

#ifdef _OPENMP
#pragma omp parallel
#endif
        {
            int     i;
            int     threadID, threadIterStart, threadIterEnd;

            GetThreadIterationIndices(home->cellCount, &threadID,
                                      &threadIterStart, &threadIterEnd);

            for (i = threadIterStart; i < threadIterEnd; i++) {
                int     cellID;
                Cell_t  *cell;
    
                cellID = home->cellList[i];
                cell = LookupCell(home, cellID);

                GetCellDomainList(home, cellID, &cell->domCount, &cell->domains);
            }

        }  /* end omp parallel section */

        return;
}
