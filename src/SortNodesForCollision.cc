/*-------------------------------------------------------------------------
 *
 *      Function:    SortNodesForCollision
 *      Description: Impose a secondary cell grid (of cell2s) over
 *                   the standard cells encompassed by this domain
 *                   and assign each native node into one of these
 *                   secondary cells for a fast screening of well-
 *                   separated nodes in collision handling.
 *
 *      The cell2 queues are maintained in a single array (cell2QentArray)
 *      where each element of the array contains both a pointer to a node
 *      in the associated cell2, and the index in the array of the next
 *      node in the cell2.  The cell2 queue heads contain the index
 *      into this array of the first element associated with a node in
 *      the corresponding cell2 (an index of -1 is interpreted as the
 *      end of the queue).
 *
 *      Note: Nodes may be deleted during collision handling.  When this
 *      occurs, the node's entry in the cell2 array is modified by
 *      zeroing out that node pointer, but the associated index to
 *      the next array element for the cell2 is left intact.  The 
 *      basic nodal structure has been modified to include the 
 *      index of the nodes entry in the cell2QentArray in order to 
 *      facilitate node removal.
 *
 *------------------------------------------------------------------------*/

#include <memory.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Node.h"
#include "Cell.h"
#include "Util.h"

/*
 *      Limits the maximum size of the cell2 grid used during
 *      collision handling
 */
#define MAXCELL2PERCELL 20


/***************************************************************************
 *
 *      Function:    InitCell2Values
 *
 *      Description: A number of values are needed for determining the 'cell2'
 *                   with which a node is to be associated during a cycle.
 *                   Some of these values remain the same throughout the
 *                   entire simulation.  Rather than recalculate these
 *                   'static' values each cycle, we'll calculate them once
 *                   in this function which is only invoked the first time
 *                   nodes must be sorted onto their 'cell2' queues.
 *
 ***************************************************************************/
static void InitCell2Values(Home_t *home,
                            int *cell2PerCellX,
                            int *cell2PerCellY,
                            int *cell2PerCellZ,
                            real8 *cell2SizeX,
                            real8 *cell2SizeY,
                            real8 *cell2SizeZ,
                            real8 *cell2MinX,
                            real8 *cell2MinY,
                            real8 *cell2MinZ,
                            int *numBodyCell2X,
                            int *numBodyCell2Y,
                            int *numBodyCell2Z)
{
        real8   maxSeparation;
        real8   Lx, Ly, Lz;
        real8   probMinX, probMinY, probMinZ;
        real8   probMaxX, probMaxY, probMaxZ;
        real8   cellSizeX, cellSizeY, cellSizeZ;
        Param_t *param;


        param = home->param;

/*
 *      Get the problem dimensions
 */
        probMinX = param->minSideX; probMaxX = param->maxSideX;
        probMinY = param->minSideY; probMaxY = param->maxSideY;
        probMinZ = param->minSideZ; probMaxZ = param->maxSideZ;
        
        Lx = probMaxX - probMinX;
        Ly = probMaxY - probMinY;
        Lz = probMaxZ - probMinZ;
        
        cellSizeX = Lx / param->nXcells;
        cellSizeY = Ly / param->nYcells;
        cellSizeZ = Lz / param->nZcells;
        
        *cell2MinX = probMinX - cellSizeX; 
        *cell2MinY = probMinY - cellSizeY; 
        *cell2MinZ = probMinZ - cellSizeZ; 
        
/*
 *      Find the maximum separation between colliding node pairs
 */   
        maxSeparation = (1.1 * param->maxSeg) + param->rann;

/*
 *      Compute the size of cell2's, such that an integral number fall within
 *      each cell, and the total number of cell2's in the minimum image problem
 *      space
 */
        *cell2PerCellX = (int) floor(cellSizeX/maxSeparation);
        *cell2PerCellY = (int) floor(cellSizeY/maxSeparation);
        *cell2PerCellZ = (int) floor(cellSizeZ/maxSeparation);

        *cell2PerCellX = MIN(*cell2PerCellX, MAXCELL2PERCELL);
        *cell2PerCellY = MIN(*cell2PerCellY, MAXCELL2PERCELL);
        *cell2PerCellZ = MIN(*cell2PerCellZ, MAXCELL2PERCELL);

        *cell2PerCellX = MAX(*cell2PerCellX, 1);
        *cell2PerCellY = MAX(*cell2PerCellY, 1);
        *cell2PerCellZ = MAX(*cell2PerCellZ, 1);

        *cell2SizeX = cellSizeX / *cell2PerCellX;
        *cell2SizeY = cellSizeY / *cell2PerCellY;
        *cell2SizeZ = cellSizeZ / *cell2PerCellZ;

        *numBodyCell2X = *cell2PerCellX * param->nXcells;
        *numBodyCell2Y = *cell2PerCellY * param->nYcells;
        *numBodyCell2Z = *cell2PerCellZ * param->nZcells;

        return;
}


void SortNodesForCollision(Home_t *home)
{
        int            i, j, k;
        int            minDomCellX, minDomCellY, minDomCellZ;
        int            maxDomCellX, maxDomCellY, maxDomCellZ;
        int            minDomCell2X, minDomCell2Y, minDomCell2Z;
        int            maxDomCell2X, maxDomCell2Y, maxDomCell2Z;
        int            maxNativeCellX, maxNativeCellY, maxNativeCellZ;
        int            minNativeCellX, minNativeCellY, minNativeCellZ;
        int            maxNativeCell2X, maxNativeCell2Y, maxNativeCell2Z;
        int            minNativeCell2X, minNativeCell2Y, minNativeCell2Z;
        int            numCell2X, numCell2Y, numCell2Z, numCell2s;
        int            nextQent = 0;
        Param_t        *param;
        static int     initDone = 0;
        static int     cell2PerCellX, cell2PerCellY, cell2PerCellZ;
        static int     numBodyCell2X, numBodyCell2Y, numBodyCell2Z;
        static real8   cell2SizeX, cell2SizeY, cell2SizeZ;
        static real8   cell2MinX, cell2MinY, cell2MinZ;
        static real8   centerX, centerY, centerZ;
        

        param = home->param;

/*
 *      The first time into this routine, set some values that won't change
 *      for the remainder of the simulation
 */
        if (!initDone++) {
            InitCell2Values(home,
                            &cell2PerCellX, &cell2PerCellY, &cell2PerCellZ,
                            &cell2SizeX, &cell2SizeY, &cell2SizeZ,
                            &cell2MinX, &cell2MinY, &cell2MinZ,
                            &numBodyCell2X, &numBodyCell2Y, &numBodyCell2Z);
        }

/*
 *      Find the center of the current domain
 */
        centerX = 0.5 * (home->domXmin + home->domXmax);
        centerY = 0.5 * (home->domYmin + home->domYmax);
        centerZ = 0.5 * (home->domZmin + home->domZmax);
        
/*
 *      Find the current min and max cells for this domain,
 *      including ghost cells
 */
        maxDomCellX = 0;
        maxDomCellY = 0;
        maxDomCellZ = 0;

        minDomCellX = param->nXcells + 2;
        minDomCellY = param->nYcells + 2;
        minDomCellZ = param->nZcells + 2;

/*
 *      NOTE: we *could* thread this loop, but the number of iterations
 *            is likely to be small in most cases AND the overhead for
 *            synchronizing inter-thread access to the values being set
 *            would likely cost more than we would gain from threading
 *            the loop.
 */
        for (i = 0; i < home->nativeCellCount; i++) {
            int     ix, iy, iz;
            int     cellIndex;
            Cell_t  *cell;

            cellIndex = home->cellList[i];
            cell = LookupCell(home, cellIndex);

            ix = cell->xIndex;
            iy = cell->yIndex;
            iz = cell->zIndex;

            if (minDomCellX > ix) minDomCellX = ix;
            if (maxDomCellX < ix) maxDomCellX = ix;
            if (minDomCellY > iy) minDomCellY = iy;
            if (maxDomCellY < iy) maxDomCellY = iy;
            if (minDomCellZ > iz) minDomCellZ = iz;
            if (maxDomCellZ < iz) maxDomCellZ = iz;
        }

        minNativeCellX = minDomCellX;
        minNativeCellY = minDomCellY;
        minNativeCellZ = minDomCellZ;

        maxNativeCellX = maxDomCellX;
        maxNativeCellY = maxDomCellY;
        maxNativeCellZ = maxDomCellZ;

        minDomCellX--;
        minDomCellY--;
        minDomCellZ--;

        maxDomCellX++;
        maxDomCellY++;
        maxDomCellZ++;

/*
 *      Determine the min and max cell2's for this domain
 */
        minDomCell2X = minDomCellX * cell2PerCellX;
        maxDomCell2X = (maxDomCellX + 1) * cell2PerCellX - 1;
        minDomCell2Y = minDomCellY * cell2PerCellY;
        maxDomCell2Y = (maxDomCellY + 1) * cell2PerCellY - 1;
        minDomCell2Z = minDomCellZ * cell2PerCellZ;
        maxDomCell2Z = (maxDomCellZ + 1) * cell2PerCellZ - 1;

        minNativeCell2X = minNativeCellX * cell2PerCellX;
        maxNativeCell2X = (maxNativeCellX+1) * cell2PerCellX - 1;
        minNativeCell2Y = minNativeCellY * cell2PerCellY;
        maxNativeCell2Y = (maxNativeCellY+1) * cell2PerCellY - 1;
        minNativeCell2Z = minNativeCellZ * cell2PerCellZ;
        maxNativeCell2Z = (maxNativeCellZ+1) * cell2PerCellZ - 1;

        numCell2X = maxDomCell2X - minDomCell2X + 1;
        numCell2Y = maxDomCell2Y - minDomCell2Y + 1;
        numCell2Z = maxDomCell2Z - minDomCell2Z + 1;

        numCell2s = numCell2X * numCell2Y * numCell2Z;

/*
 *      Need to save the number cell2 geometry for collision handling
 *      which may need to convert cell2 IDs into the corresponding
 *      X, Y and Z indices of the cell2 in the cell2 grid.
 */
        home->cell2nx = numCell2X;
        home->cell2ny = numCell2Y;
        home->cell2nz = numCell2Z;

/*
 *      Allocate and initialize the array of cell2 queue heads for this domain
 *      and the single array containing the queues of node pointer/array index
 *      pairs.
 */
        home->cell2 = (int *)realloc(home->cell2, numCell2s * sizeof(int));

#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < numCell2s; i++) {
            home->cell2[i] = -1;
        }

        home->cell2QentArray = (C2Qent_t *)realloc(home->cell2QentArray,
                                                   home->newNodeKeyPtr *
                                                   sizeof(C2Qent_t));
/*
 *      Loop through all nodes. Queue each node onto one of the cell2's. A node
 *      may fall into more than one cell2 if the domain contains the whole
 *      problem in one or more directions. In this case, queue onto the lowest
 *      numbered image
 */
#ifdef _OPENMP
#pragma omp parallel for
#endif
        for (i = 0; i < home->newNodeKeyPtr; i++) {
            int    baseX, baseY, baseZ;
            int    c2BaseIndex, c2Index;
            int    cell2X, cell2Y, cell2Z;
            int    cell2MinXus, cell2MinYus, cell2MinZus;
            int    cell2Xplus, cell2Yplus, cell2Zplus;
            real8  x, y, z;
            Node_t *node;

            home->cell2QentArray[i].node = (Node_t *)NULL;
            home->cell2QentArray[i].next = -1;

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            x = node->x;
            y = node->y;
            z = node->z;

/*
 *          Get the image of the point closest to the center of the
 *          current domain.
 */
            PBCPOSITION(param, centerX, centerY, centerZ, &x, &y, &z);

            cell2X = (int) floor((x - cell2MinX) / cell2SizeX);
            cell2Y = (int) floor((y - cell2MinY) / cell2SizeY);
            cell2Z = (int) floor((z - cell2MinZ) / cell2SizeZ);

/*
 *          If a native node falls outside native cell2 area 
 *          force it into the nearest native cell2.
 */
            if (cell2X < minNativeCell2X) {
                cell2X = minNativeCell2X;
            }
            if (cell2X > maxNativeCell2X) {
                cell2X = maxNativeCell2X;
            }

            if (cell2Y < minNativeCell2Y) {
                cell2Y = minNativeCell2Y;
            }
            if (cell2Y > maxNativeCell2Y) {
                cell2Y = maxNativeCell2Y;
            }

            if (cell2Z < minNativeCell2Z) {
                cell2Z = minNativeCell2Z;
            }
            if (cell2Z > maxNativeCell2Z) {
                cell2Z = maxNativeCell2Z;
            }

/*
 *          Find minimum cell2 containing node in the domain range.
 *          This minimum cell will be adjusted to the domain range, rather
 *          than problem range, of cell2 indices.
 *
 *          X direction
 */
            cell2MinXus = cell2X - numBodyCell2X;
            cell2Xplus  = cell2X + numBodyCell2X;

            if (cell2MinXus >= minDomCell2X)
                baseX = cell2MinXus - minDomCell2X;
            else if (cell2X >= minDomCell2X && cell2X <= maxDomCell2X)
                baseX = cell2X - minDomCell2X;
            else if (cell2Xplus <= maxDomCell2X)
                baseX = cell2Xplus - minDomCell2X;
            else 
                continue;  /* node moved outside of ghost cells; ignore */

/*
 *          Y direction
 */
            cell2MinYus = cell2Y - numBodyCell2Y;
            cell2Yplus  = cell2Y + numBodyCell2Y;

            if (cell2MinYus >= minDomCell2Y)
                baseY = cell2MinYus - minDomCell2Y;
            else if (cell2Y >= minDomCell2Y && cell2Y <= maxDomCell2Y)
                baseY = cell2Y - minDomCell2Y;
            else if (cell2Yplus <= maxDomCell2Y)
                baseY = cell2Yplus - minDomCell2Y;
            else
                continue; /* node moved outside of ghost cells; ignore */

/*
 *          Z direction
 */
            cell2MinZus = cell2Z - numBodyCell2Z;
            cell2Zplus  = cell2Z + numBodyCell2Z;

            if (cell2MinZus >= minDomCell2Z)
                baseZ = cell2MinZus - minDomCell2Z;
            else if (cell2Z >= minDomCell2Z && cell2Z <= maxDomCell2Z)
                baseZ = cell2Z - minDomCell2Z;
            else if (cell2Zplus <= maxDomCell2Z)
                baseZ = cell2Zplus - minDomCell2Z;
            else
                continue; /* node moved outside of ghost cells; ignore */

/*
 *          Make cell2 index relative to the domain
 */
            cell2X -= minDomCell2X;
            cell2Y -= minDomCell2Y;
            cell2Z -= minDomCell2Z;

/*
 *          Queue the node on the minimum cell2 that it falls into, but save the
 *          ghost cell2 index in the node, since this is where it starts its
 *          neighbor search from.
 *
 *          Use a 'critical' section to prevent multiple threads from updating
 *          the cell2 stuff simultaneously
 */
            c2BaseIndex = EncodeCell2Idx(home, baseX, baseY, baseZ);
            c2Index = EncodeCell2Idx(home, cell2X, cell2Y, cell2Z);

            node->cell2Idx = c2Index;

#ifdef _OPENMP
#pragma omp critical (CRIT_SORT_FOR_COLLISION)
#endif
            {
                node->cell2QentIdx = nextQent;
                home->cell2QentArray[nextQent].node = node;
                home->cell2QentArray[nextQent].next = home->cell2[c2BaseIndex];
                home->cell2[c2BaseIndex] = nextQent;

                nextQent++;

            }  /* end "omp critical" section */
        }

/*
 *      Loop through the cell2s in the domain. If a cell2 is a higher image
 *      of a base cell2 also contained in the domain space, set the higher
 *      image's queue to point to the base image cell2.
 */
        for (i = 0; i < numCell2X; i++) {
            int xIndex, yIndex, zIndex, cell2Index, baseIndex;
        
            if (i >= numBodyCell2X) {
                xIndex = i - numBodyCell2X;
            } else {
                xIndex = i;
            }
        
            for (j = 0; j < numCell2Y; j++) {
        
                if (j >= numBodyCell2Y) {
                    yIndex = j - numBodyCell2Y;
                } else {
                    yIndex = j;
                }
        
                for (k = 0; k < numCell2Z; k++) {
        
                    if (k >= numBodyCell2Z) {
                        zIndex = k - numBodyCell2Z;
                    } else {
                        zIndex = k;
                    }
        
                    baseIndex = EncodeCell2Idx(home, xIndex, yIndex, zIndex);
                    cell2Index = EncodeCell2Idx(home, i, j, k);
                    home->cell2[cell2Index] = home->cell2[baseIndex];
                }
            }
        }
        
        return;
}
