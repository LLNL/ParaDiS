/**************************************************************************
 *
 *  Module:      CellTable.c
 *  Description: Contains functions for managing the hash table
 *               for the standard cell structures.
 *
 *  Includes functions:
 *      AddCellToTable()
 *      FreeCellTable()
 *      LookupCell()
 *
 ***************************************************************************/

#include "mpi_portability.h"

#include <Home.h>


/*---------------------------------------------------------------------------
 *
 *  Function:     AddCellToTable
 *  Description:  Allocate a new cell structure, assign it the
 *                specified ID, and add it to the cell hash table
 *                in the home struct.
 *
 *                Note: If the cell has already been allocated, no
 *                      action is taken and a pointer to the
 *                      previously allocated cell is returned.
 *
 *  Arguments:
 *      IN: cellID    Cell ID of the cell to be allocated.  Note: this
 *                    cell id is assumed to have been adjusted to 
 *                    accommodate ghost cells.
 *
 *  Returns: Pointer to the cell structure associted with cell <cellID>
 *
 *--------------------------------------------------------------------------*/
Cell_t *AddCellToTable(Home_t *home, int cellID)
{
    int      hashVal;
    Cell_t   *cell, *tmpCell, **cellTable;

/*
 *  If the cell is already on the table, just mark it as in use
 *  and return.
 */
    cellTable = home->cellTable;

    cell = LookupCell(home, cellID);

    if (cell != (Cell_t *)NULL) {
        cell->inUse = 1;
        return(cell);
    }

/*
 *  Insert the new cell into the table
 */
    hashVal = cellID % CELL_HASH_TABLE_SIZE;
    cell = (Cell_t *)calloc(1, sizeof(Cell_t));

    cell->cellID = cellID;

    if (cellTable[hashVal] == (Cell_t *)NULL) {
        cellTable[hashVal] = cell;
    } else {
        tmpCell = cellTable[hashVal];
        cell->next = tmpCell;
        cellTable[hashVal] = cell;
    }

    return(cell);
}


/*---------------------------------------------------------------------------
 *
 *  Function:     FreeCellTable
 *  Description:  Free all dynamic memory associated with each cell in the
 *                cell hash table, and remove all cells from the table.
 *
 *--------------------------------------------------------------------------*/
void FreeCellTable(Home_t *home)
{
    int i;

/*
 *  Loop through all the hash table entries and free the linked
 *  list of cells for each.
 */
    for (i = 0; i < CELL_HASH_TABLE_SIZE; i++) {
        Cell_t *cell, *currCell;

        cell = home->cellTable[i];

        while (cell != (Cell_t *)NULL) {

            currCell = cell;

            if (cell->nbrList) {
                free(cell->nbrList);
            }

            if (cell->domains) {
                free(cell->domains);
            }

/*
 *          For base cells (non-periodic), we have to destroy any
 *          cell lock that was allocated.  Periodic cells have
 *          baseIdx values >= zero.
 */
            if (cell->baseIdx < 0) {
                DESTROY_LOCK(&cell->cellLock);
            }

            cell = cell->next;

            free(currCell);
        }

        home->cellTable[i] = (Cell_t *)NULL;
    }

    home->cellCount = 0;

    return;
}


/*---------------------------------------------------------------------------
 *
 *  Function:     LookupCell
 *  Description:  Find a specified cell in the cell hash table and return
 *                to the caller a pointer to the cell.
 *
 *  Arguments:
 *      cellID    ID number of the cell to locate.
 *
 *  Returns:  Pointer to the specified cell if found, NULL in all
 *            other cases.
 *
 *--------------------------------------------------------------------------*/
Cell_t *LookupCell(Home_t *home, int cellID)
{
        int    hashVal;
        Cell_t *cell = (Cell_t *)NULL;

        hashVal = cellID % CELL_HASH_TABLE_SIZE;
        cell = home->cellTable[hashVal];

        while (cell != (Cell_t *)NULL) {

            if (cell->cellID == cellID) {
                break;
            }

            cell = cell->next;
        }

        return(cell);
}
