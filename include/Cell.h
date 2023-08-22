#pragma once

#ifndef _PDS_CELL_H
#define _PDS_CELL_H

#include <stdio.h>

#include "Typedefs.h"
#include "Node.h"

class Cell_t
{
    public:
        int     cellID;                     // cell ID encoded from xyz indices to a single int

        int     inUse;                      // toggle indicating if the associated cell struct
                                            // is in use for a specific timestep

        Cell_t *next;                       // next cell in the hash table list

        Node_t *nodeQ;                      // queue head of nodes currently in this cell
        int     nodeCount;                  // number of nodes on nodeQ

        int     inclusionQ;                 // index of first inclusions in the inclusion queue for this cell
        int     inclusionCount;             // number of inclusions in this cell

        int    *nbrList;                    // list of neighbor cell encoded indices
        int     nbrCount;                   // number of neighbor cells

        int    *domains;                    // domains that intersect cell (encoded indices)
        int     domCount;                   // number of intersecting domains

        int     baseIdx;                    // encoded index of corresp' base cell (-1 if not  periodic)

        real8   xShift, yShift, zShift;     // if periodic, amount to shift corresponding base coordinates

        int     xIndex, yIndex, zIndex;     // indices of the cell in each dimension.
                                            // Note that these indices are in the range 0 <= xIndex <= nXcells
                                            // and likewise for Y and Z to account for the ghost cells

        double  center[3];                  // Defines the center position of the cell.
                                            // Note: This is curently only set for 'base' cells (i.e. cells in
                                            // the primary image), not the the extra layer of ghost cells.

#ifdef _OPENMP
        omp_lock_t cellLock;                // For use with threading to gain exclusive access to cell-specific data.
                                            // WARNING: This is only available for cells native to the domain for now!
#endif

    public:
        Cell_t(void);
       ~Cell_t();

       void Print (FILE *fd=stdout) const;
};

// Prototype functions needed for managing the standard cell hash table

extern Cell_t   *AddCellToTable   (Home_t *home, int cellID);
extern void      FreeCellTable    (Home_t *home);
extern Cell_t   *LookupCell       (Home_t *home, int cellID);

extern void      Cell_Print       (FILE *fd, const Cell_t  *cell);
extern void      Cell_Table_Print (FILE *fd, const Cell_t **ctbl, const int n);

#endif
