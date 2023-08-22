/****************************************************************************
 *
 *      Module:      RSDecomp.c
 *      Description: Contains functions needed for both generating and
 *                   accessing a domain decomposition based on the original
 *                   recursive sectioning algorithm used from the early
 *                   days of ParaDis (aka DD3d).
 *
 *                   The basic idea of this algorithm is to perform a
 *                   domain decomposition over a 3-cycle period.  During
 *                   the first cycle, the entire problem space will be
 *                   sectioned along the X axis into <param->nXdoms> slabs
 *                   such that the computational cost of each slab will be
 *                   roughly equivalent.  The next cycle each slab will be
 *                   sectioned along the Y axis into <param->nYdoms> 
 *                   columns such that the computational cost of each
 *                   column in a particular slab is roughly equivalent.
 *                   On the third cycle, every column will be sectioned
 *                   along the Z axis into <param->nZdoms> chunks such that
 *                   the computational cost of every chunck within the
 *                   same column is roughly equivalent.
 *
 *
 *      Includes public functions:
 *          DLBalanceX()
 *          DLBalanceY()
 *          DLBalanceZ()
 *          FindRSDecompCoordDomain()
 *          FreeRSDecomp()
 *          GetAllRSDecompBounds()
 *          GetRSDecompLocalDomainBounds()
 *          GetRSDecompCellDomainList()
 *          ReadBinRSDecompBounds()
 *          ReadRSDecompBounds()
 *          UniformRSDecomp()
 *          WriteRSDecompBound()
 *          XPlotRSDecomp()
 *
 *      Includes private functions:
 *          DLBnewBounds()
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Decomp.h"
#include "RSDecomp.h"
#include "DisplayC.h"
#include "Restart.h"

#ifdef USE_HDF
#include <hdf5.h>
#endif

/*
 *      Define the minimum imbalance threshold determining if
 *      load-balance is required.  If imbalance is below this
 *      value, no rebalancing is done.
 */
#define DLBTHRESHOLD  0.0

/*
 *      Define the maximum fraction of the old partition size a
 *      boundary is permitted to shift in a single rebalance.
 */
#define MAX_SHIFT 0.05

/*
 *      Define an increment size by which the list of domains
 *      intersecting a cell will be increased when determining
 *      cell/domain intersection lists.
 */
#define DOMAIN_CNT_INCREMENT 20


/*---------------------------------------------------------------------------
 *
 *      Function:    GetRSDecompLocalDomainBounds
 *      Description: Extract the local domain boundaries from the
 *                   current domain decomposition and store them in
 *                   the <home> structure.  Makes life easier later
 *                   on in the code if we don't have to keep searching
 *                   the decomposition for the local bounds
 *
 *      Arguments:
 *          decomp  Pointer to current domain decomposition
 *
 *--------------------------------------------------------------------------*/
void GetRSDecompLocalDomainBounds(Home_t *home, RSDecomp_t *decomp)
{
        int       xIndex, yIndex, zIndex;

        DecodeDomainIdx(home, home->myDomain, &xIndex, &yIndex, &zIndex);

        home->domXmin = decomp->domBoundX[xIndex  ];
        home->domXmax = decomp->domBoundX[xIndex+1];

        home->domYmin = decomp->domBoundY[xIndex][yIndex  ];
        home->domYmax = decomp->domBoundY[xIndex][yIndex+1];

        home->domZmin = decomp->domBoundZ[xIndex][yIndex][zIndex  ];
        home->domZmax = decomp->domBoundZ[xIndex][yIndex][zIndex+1];

        return;
}


/*---------------------------------------------------------------------------
 *
 *    Function:     GetRSDecompCellDomainList
 *    Description:  Find out what domains intersect the specified cell.
 *                  Return a pointer to the list of domains to the caller
 *                  along with a count of the number of domains on the
 *                  returned list.
 *      Arguments:
 *          cellID   ID of the cell in question
 *          domCount Location at which to return to the caller, the
 *                   number of domain IDs placed in the <domList> array.
 *          domList  Location at which to return to the caller a pointer
 *                   to an array of the IDs of domains intersecting this
 *                   cell.  Note: the caller is responsible for freeing
 *                   this list when it is no longer necessary.
 *     
 *--------------------------------------------------------------------------*/
void GetRSDecompCellDomainList(Home_t *home, int cellID, int *domCount,
                              int **domList)
{
        int       xDom, yDom, zDom, domID;
        int       *domains, numDomains, maxDomains;
        int       xCell, yCell, zCell;
        real8     xCellSize, yCellSize, zCellSize;
        real8     xCellMin, yCellMin, zCellMin;
        real8     xCellMax, yCellMax, zCellMax;
        real8     *domBoundX;
        real8     **domBoundY;
        real8     ***domBoundZ;
        Param_t   *param;
        RSDecomp_t *decomp;
        Cell_t    *cell;


        param = home->param;
        decomp = (RSDecomp_t *)home->decomp;

        domBoundX = decomp->domBoundX;
        domBoundY = decomp->domBoundY;
        domBoundZ = decomp->domBoundZ;

/*
 *      Get the indices of the cell from the cell struct if we have
 *      it, otherwise calculate it.
 */
        cell = LookupCell(home, cellID);

        if (cell != (Cell_t *)NULL) {
            xCell = cell->xIndex;
            yCell = cell->yIndex;
            zCell = cell->zIndex;
        } else {
            DecodeCellIdx(home, cellID, &xCell, &yCell, &zCell);
        }

/*
 *      Get the min and max coordinates for this cell.
 *
 *      WARNING: The cell min and max coordinate values must be
 *      computed in the same manner or neighboring domains can end
 *      up with slightly different locations for a shared cell
 *      boundary (due to numeric round-off.)
 */
        xCellSize = (param->maxSideX - param->minSideX) / param->nXcells;
        yCellSize = (param->maxSideY - param->minSideY) / param->nYcells;
        zCellSize = (param->maxSideZ - param->minSideZ) / param->nZcells;

        xCellMin = rint((param->minSideX) + ((xCell-1) * xCellSize));
        yCellMin = rint((param->minSideY) + ((yCell-1) * yCellSize));
        zCellMin = rint((param->minSideZ) + ((zCell-1) * zCellSize));

        xCellMax = rint((param->minSideX) + (xCell * xCellSize));
        yCellMax = rint((param->minSideY) + (yCell * yCellSize));
        zCellMax = rint((param->minSideZ) + (zCell * zCellSize));

/*
 *      Allocate a list to hold the intersecting domains
 */
        numDomains = 0;
        maxDomains = DOMAIN_CNT_INCREMENT;
        domains = (int *)malloc(maxDomains * sizeof(int));

/*
 *      Loop through the x partitions of the problem. If an x
 *      partition doesn't overlap the x bounds of the cell, skip
 *      to next x partition. If there is overlap, loop through y
 *      partitions for that x partition looking for overlap with
 *      cell. If not, skip to next y partition. If overlap, loop
 *      through all z partitions for the current x and y partition,
 *      looking for overlap with cell If overlap in all dimensions,
 *      then that domain overlaps the cell, so save its encoded
 *      index in the cell's domains list
 */
        for (xDom = 0; xDom < param->nXdoms; xDom++) {

            if (domBoundX[xDom  ] >= xCellMax ||
                domBoundX[xDom+1] <= xCellMin) {
                continue;
            }

            for (yDom = 0; yDom < param->nYdoms; yDom++) {

                if (domBoundY[xDom][yDom  ] >= yCellMax ||
                    domBoundY[xDom][yDom+1] <= yCellMin) {
                    continue;
                }

                for (zDom = 0; zDom<param->nZdoms; zDom++) {

                    if (domBoundZ[xDom][yDom][zDom  ] >= zCellMax ||
                        domBoundZ[xDom][yDom][zDom+1] <= zCellMin) {
                        continue;
                    }

/*
 *                  If a domain survived all the above tests, then
 *                  it intersects the cell. Save its domain ID in
 *                  the domain list.
 */
                    domID = EncodeDomainIdx(home, xDom, yDom, zDom);

                    domains[numDomains++] = domID;

                    if (numDomains >= maxDomains) {
                        maxDomains += DOMAIN_CNT_INCREMENT;
                        domains = (int *)realloc(domains, maxDomains *
                                                 sizeof(int));
                    }
                }
            }
        }

/*
 *      Pass the domain list (and count) back to the caller.
 */
        *domCount = numDomains;
        *domList  = domains;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReadRSDecompBounds
 *      Description: Read the domain geometry and decomposition (if any)
 *                   from the nodal data file and return it to the
 *                   caller (if requested).
 *
 *      Arguments:
 *          fp          File pointer to the opened nodal data file;
 *                      should be positioned in the file such that
 *                      the next item in the file is the decomposition.
 *          numXDoms    Number of domains in X dimension of decomposition
 *                      contained in the file
 *          numYDoms    Number of domains in Y dimension of decomposition
 *                      contained in the file
 *          numZDoms    Number of domains in X dimension of decomposition
 *                      contained in the file
 *          saveDecomp  Flag indicating if decomposition is to be saved
 *                      and returned to the caller.  0 == don't save,
 *                      1 == save.
 *          oldDecomp   Location in which to return to the caller a
 *                      pointer to the old domain decomposition read
 *                      from the file (if necessary).
 *
 *------------------------------------------------------------------------*/
void ReadRSDecompBounds(void **filePtr, int numXDoms, int numYDoms,
                        int numZDoms, int saveDecomp, RSDecomp_t **oldDecomp)
{
        int       i, j, k;
        real8     *xBnd=0, **yBnd=0, ***zBnd=0;
        char      inLine[512];
        FILE      *fp = (FILE *)*filePtr;
        RSDecomp_t *decomp=0;

/*
 *      If we'll need to return the decompostion data to the caller
 *      allocate the needed arrays here.
 */
        if (saveDecomp) {

            *oldDecomp = (RSDecomp_t *)NULL;

            decomp = (RSDecomp_t *)malloc(sizeof(RSDecomp_t));

            xBnd = (real8 *) malloc((numXDoms+1)*sizeof(real8));
            yBnd = (real8 **) malloc(numXDoms*sizeof(real8 *));
            zBnd = (real8 ***) malloc(numXDoms*sizeof(real8 **));

            for (i = 0; i < numXDoms; i++) {

                yBnd[i] = (real8 *)malloc((numYDoms+1) * sizeof(real8));
                zBnd[i] = (real8 **)malloc(numYDoms * sizeof(real8));

                for (j = 0; j < numYDoms; j++) {

                    zBnd[i][j] =
                            (real8 *)malloc((numZDoms+1) * sizeof(real8));
                }
            }
        }
/*
 *      Read the old domain decomposition in and, if necessary, store
 *      store the boundaries for the caller.
 */
        for (i = 0; i < numXDoms; i++) {
            Getline(inLine, sizeof(inLine), fp);
            if (saveDecomp) sscanf(inLine, "%lf", &xBnd[i]);
            for (j = 0; j < numYDoms; j++) {
                Getline(inLine, sizeof(inLine), fp);
                if (saveDecomp) sscanf(inLine, "%lf", &yBnd[i][j]);
                for (k = 0; k < numZDoms; k++) {
                    Getline(inLine, sizeof(inLine), fp);
                    if (saveDecomp) sscanf(inLine, "%lf", &zBnd[i][j][k]);
                }
                Getline(inLine, sizeof(inLine), fp);
                if (saveDecomp) sscanf(inLine, "%lf", &zBnd[i][j][k]);
            }
            Getline(inLine, sizeof(inLine), fp);
            if (saveDecomp) sscanf(inLine, "%lf", &yBnd[i][j]);
        }
        Getline(inLine, sizeof(inLine), fp);
        if (saveDecomp) sscanf(inLine, "%lf", &xBnd[i]);

/*
 *      If the caller needs a copy of this old domain decomposition
 *      return him a pointer to it.
 */
        if (saveDecomp) {

            decomp->domBoundX = xBnd;
            decomp->domBoundY = yBnd;
            decomp->domBoundZ = zBnd;

            *oldDecomp = decomp;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReadBinRSDecompBounds
 *      Description: Read the domain decomposition from the HDF5 data
 *                   file and return it to the caller.
 *
 *      Arguments:
 *          fp          File pointer to the opened HDF5 data file;
 *          numXDoms    Number of domains in X dimension of decomposition
 *                      contained in the file
 *          numYDoms    Number of domains in Y dimension of decomposition
 *                      contained in the file
 *          numZDoms    Number of domains in X dimension of decomposition
 *                      contained in the file
 *          oldDecomp   Location in which to return to the caller a
 *                      pointer to the old domain decomposition read
 *                      from the file.
 *
 *------------------------------------------------------------------------*/
void ReadBinRSDecompBounds(void *filePtr, int numXDoms, int numYDoms,
                           int numZDoms, RSDecomp_t **oldDecomp)
{
#ifdef USE_HDF
        int        i, j, k, numValues, offset;
        int        status;
        real8      *xBnd, **yBnd, ***zBnd, *boundsBuf;
        hid_t      *fileID = (hid_t *)filePtr;
        RSDecomp_t *decomp;


/*
 *      Allocate arrays for the domain decomposition
 */
        numValues = numXDoms * (numYDoms * (numZDoms+1) + (numYDoms+1)) +
                    (numXDoms+1);
        boundsBuf = (real8 *)malloc(numValues * sizeof(real8));

        *oldDecomp = (RSDecomp_t *)NULL;

        decomp = (RSDecomp_t *)malloc(sizeof(RSDecomp_t));

        xBnd = (real8 *) malloc((numXDoms+1)*sizeof(real8));
        yBnd = (real8 **) malloc(numXDoms*sizeof(real8 *));
        zBnd = (real8 ***) malloc(numXDoms*sizeof(real8 **));

        for (i = 0; i < numXDoms; i++) {
            yBnd[i] = (real8 *)malloc((numYDoms+1) * sizeof(real8));
            zBnd[i] = (real8 **)malloc(numYDoms * sizeof(real8));
            for (j = 0; j < numYDoms; j++) {
                zBnd[i][j] = (real8 *)malloc((numZDoms+1) * sizeof(real8));
            }
        }

/*
 *      Read the old domain decomposition in and store
 *      the boundaries for the caller.
 */
        status = ReadHDFDataset(*fileID, "/decomposition", H5T_NATIVE_DOUBLE,
                                numValues, boundsBuf);
        if (status != 0) {
            Fatal("ReadBinRSDecompBounds: Error reading decomposition");
        }
        
        offset = 0;

        for (i = 0; i < numXDoms; i++) {
            xBnd[i] = boundsBuf[offset++];
            for (j = 0; j < numYDoms; j++) {
                yBnd[i][j] = boundsBuf[offset++];
                for (k = 0; k < numZDoms; k++) {
                    zBnd[i][j][k] = boundsBuf[offset++];
                }
                zBnd[i][j][k] = boundsBuf[offset++];
            }
            yBnd[i][j] = boundsBuf[offset++];
        }
        xBnd[i] = boundsBuf[offset++];

        decomp->domBoundX = xBnd;
        decomp->domBoundY = yBnd;
        decomp->domBoundZ = zBnd;

        *oldDecomp = decomp;

        free(boundsBuf);
#endif

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    WriteRSDecompBounds
 *      Description: Writes the domain geometry and boundaries for
 *                   the specified decomposition into the restart file.
 *
 *      Arguments:
 *          fp          open file descriptor for the restart file being
 *                      written
 *          decomp      pointer to current domain decomposition.
 *
 *------------------------------------------------------------------------*/
void WriteRSDecompBounds(Home_t *home, FILE *fp, RSDecomp_t *decomp)
{
        int     i, j, k, numXdoms, numYdoms, numZdoms;
        Param_t *param;

        param = home->param;

        numXdoms = param->nXdoms;
        numYdoms = param->nYdoms;
        numZdoms = param->nZdoms;

        fprintf(fp, " %11.4f\n", decomp->domBoundX[0]);

        for (i = 1; i <= numXdoms; i++) {
            fprintf(fp, "     %11.4f\n", decomp->domBoundY[i-1][0]);
            for (j = 1; j <= numYdoms; j++) {
                for (k = 0; k <= numZdoms; k++) {
                    fprintf(fp, "         %11.4f\n",
                            decomp->domBoundZ[i-1][j-1][k]);
                }
                fprintf(fp,"     %11.4f\n", decomp->domBoundY[i-1][j]);
            }
            fprintf(fp," %11.4f\n", decomp->domBoundX[i]);
        }
        fprintf(fp, "\n");

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetAllRSDecompBounds
 *      Description: Create a one-dimensional array containing all
 *                   boundaries defining the RS decomposition.
 *
 *      Arguments:
 *          decomp      pointer to current domain decomposition.
 *          bounds      Caller provided array in which to return
 *                      the domain boundaries.  NOTE: <bounds> array
 *                      must be large enough to hold all domain
 *                      boundaries or the results are undefined.
 *
 *------------------------------------------------------------------------*/
void GetAllRSDecompBounds(Home_t *home, RSDecomp_t *decomp, real8 *bounds)
{
        int i, j, k;
        int numXdoms, numYdoms, numZdoms;
        int offset = 0;

        numXdoms = home->param->nXdoms;
        numYdoms = home->param->nYdoms;
        numZdoms = home->param->nZdoms;

        bounds[offset++] = decomp->domBoundX[0];
        for (i = 1; i <= numXdoms; i++) {
            bounds[offset++] = decomp->domBoundY[i-1][0];
            for (j = 1; j <= numYdoms; j++) {
                for (k = 0; k <= numZdoms; k++) {
                    bounds[offset++] = decomp->domBoundZ[i-1][j-1][k];
                }
                bounds[offset++] = decomp->domBoundY[i-1][j];
            }
            bounds[offset++] = decomp->domBoundX[i];
        }

        return;
}


#ifndef NO_XWINDOW
/*---------------------------------------------------------------------------
 *
 *      Function:     XPlotRSDecomp
 *      Description:  Plots the domain boundaries for a recursive sectioning
 *                    decomposition in the active X-window display.
 *
 *      Arguments:
 *          xMin      Minimum permitted coordinate in the X dimension
 *          yMin      Minimum permitted coordinate in the Y dimension
 *          zMin      Minimum permitted coordinate in the Z dimension
 *          lMax      Length of the problem space in the largest dimension
 *          color     color to use for the lines defining the boundaries
 *          lineWidth self-explanatory
 *
 *-------------------------------------------------------------------------*/
void XPlotRSDecomp(Home_t *home, real8 xMin, real8 yMin, real8 zMin,
                  real8 lMax, int color, real8 lineWidth)
{
        int       i, j, k;
        int       numXDoms, numYDoms, numZDoms;
        real8     x, y, z;
        real8     x2, y2, z2;
        RSDecomp_t *decomp;
        Param_t   *param;

        param = home->param;
        decomp = (RSDecomp_t *)home->decomp;

        numXDoms = param->nXdoms;
        numYDoms = param->nYdoms;
        numZDoms = param->nZdoms;

/*
 *      Loop through boundaries and plot boundary lines.
 *
 *      Note: For plotting purposes, positions for the endpoints
 *      of the lines are scaled from absolute positions to
 *      values in the range -1 to +1.
 */
        for (i = 1; i <= numXDoms; i++) {

            x = decomp->domBoundX[i];
            x = (x-xMin) / lMax * 2 - 1;

            if (i < numXDoms) {
                y = -1;  z = -1;  y2 = -1;  z2 = 1;
                WinDrawLine(x,y,z,x,y2,z2,color,lineWidth/2,0);

                y = -1;  z = 1;  y2 = 1;  z2 = 1;
                WinDrawLine(x,y,z,x,y2,z2,color,lineWidth/2,0);

                y = 1;  z = 1;  y2 = 1;  z2 = -1;
                WinDrawLine(x,y,z,x,y2,z2,color,lineWidth/2,0);

                y = 1;  z = -1;  y2 = -1;  z2 = -1;
                WinDrawLine(x,y,z,x,y2,z2,color,lineWidth/2,0);
            }
            x2 = decomp->domBoundX[i-1];
            x2 = (x2-xMin) / lMax * 2 - 1;

            for (j = 1; j <= numYDoms; j++) {

                y = decomp->domBoundY[i-1][j];
                y = (y-yMin) / lMax * 2 - 1;

                if (j < numYDoms) {
                    z = -1;  z2 = 1;
                    WinDrawLine(x2,y,z,x,y,z,color, lineWidth/2,0);
                    WinDrawLine(x,y,z,x,y,z2,color, lineWidth/2,0);
                    WinDrawLine(x,y,z2,x2,y,z2,color, lineWidth/2,0);
                    WinDrawLine(x2,y,z2,x2,y,z,color, lineWidth/2,0);
                }

                y2 = decomp->domBoundY[i-1][j-1];
                y2 = (y2-yMin) / lMax * 2 - 1;
/*
                y2 = (y2-yMin) / Ly * 2 - 1;
*/

                for (k = 1; k < numZDoms; k++) {

                    z = decomp->domBoundZ[i-1][j-1][k];
                    z = (z-zMin) / lMax * 2 - 1;

                    WinDrawLine(x2,y,z,x,y,z,color, lineWidth/2,0);
                    WinDrawLine(x,y,z,x,y2,z,color, lineWidth/2,0);
                    WinDrawLine(x,y2,z,x2,y2,z,color, lineWidth/2,0);
                    WinDrawLine(x2,y2,z,x2,y,z,color, lineWidth/2,0);
                }
            }
        }

        return;
}
#endif


/*-------------------------------------------------------------------------
 *
 *      Function:    FreeRSDecomp
 *      Description: Free all arrays associated with the recursive
 *                   sectioning domain decomposition and zero the
 *                   array pointers.
 *
 *      Arguments:
 *          decomp   pointer to the domain decomposition to be freed
 *
 *------------------------------------------------------------------------*/
void FreeRSDecomp(Home_t *home, RSDecomp_t *decomp)
{
        int xDom, yDom;
        int numXDoms, numYDoms;

        numXDoms = home->param->nXdoms;
        numYDoms = home->param->nYdoms;

        if (decomp->domBoundX != (real8 *)NULL) {
            for (xDom = 0; xDom < numXDoms; xDom++) {
                for (yDom = 0; yDom < numYDoms; yDom++) {
                    free(decomp->domBoundZ[xDom][yDom]);
                    decomp->domBoundZ[xDom][yDom] = (real8 *)NULL;
                }
                free(decomp->domBoundZ[xDom]);
                decomp->domBoundZ[xDom] = (real8 **)NULL;
                free(decomp->domBoundY[xDom]);
                decomp->domBoundY[xDom] = (real8 *)NULL;
            }
            free(decomp->domBoundZ);
            decomp->domBoundZ = (real8 ***)NULL;
            free(decomp->domBoundY);
            decomp->domBoundY = (real8 **)NULL;
            free(decomp->domBoundX);
            decomp->domBoundX = (real8 *)NULL;
        }

        return;
}

/*---------------------------------------------------------------------------
 *
 *      Function:     FindRSDecompCoordDomain
 *      Description:  Determine the ID of the domain which contains the
 *                    specified coordinate.
 *
 *                    Note: coordinates outside the primary image are
 *                    treated as if they were within the closest domain
 *                    boundary.
 *
 *      Arguments:
 *          decomp    Pointer to decomposition data.
 *          x, y, z   Coordinates 
 *
 *-------------------------------------------------------------------------*/
int FindRSDecompCoordDomain(Home_t *home, RSDecomp_t *decomp, real8 x,
                           real8 y, real8 z)
{
        int     i, xDom, yDom, zDom, domID;
        Param_t *param;

        param = home->param;

/*
 *      Find the index in X dimension of the owning domain.
 */
        xDom = param->nXdoms - 1;

        for (i = 1; i < param->nXdoms; i++) {
            if (x < decomp->domBoundX[i]) {
                xDom = i - 1;
                break;
            }
        }

/*
 *      Find the index in Y dimension of the owning domain.
 */
        yDom = param->nYdoms - 1;

        for (i = 1; i < param->nYdoms; i++) {
            if (y < decomp->domBoundY[xDom][i]) {
                yDom = i - 1;
                break;
            }
        }

/*
 *      Find the index in Z dimension of the owning domain.
 */
        zDom = param->nZdoms - 1;

        for (i = 1; i < param->nZdoms; i++) {
            if (z < decomp->domBoundZ[xDom][yDom][i]) {
                zDom = i - 1;
                break;
            }
        }

        domID = EncodeDomainIdx(home, xDom, yDom, zDom);

        return(domID);
}


/*************************************************************************
 *
 *      Function:     DLBnewBounds
 *      Description:  Given the current computational weight (times[i])
 *                    for each region bounded by Bold[i] and Bold[i+1],
 *                    construct a weight density function wDensity, and
 *                    find a new partition Bnew such that the integrated
 *                    weight in each region Bnew[i] to Bnew[i+1] is as
 *                    close to possible as the target weight wTarget
 *                    for the region.
 *
 *                    Constrain Bnew[i] to lie inside the region bounded by:
 *
 *                        Bold[i] - MAX_SHIFT * (Bold[i] - Bold[i-1]) and
 *                        Bold[i] + MAX_SHIFT * (Bold[i+1] - Bold[i])
 *
 *      Arguments:
 *          Bold         Array of <numParts>+1 elements defining
 *                       the previous partition boundaries.
 *          Bnew         Array of <numParts>+1 elements in which
 *                       to return to the caller the new partition
 *                       boundaries.
 *          times        Array of <numParts> elements containing
 *                       computational weighting value associated with
 *                       each partition.
 *          numParts     Number of partitions 
 *          shiftLimit   Hard limit on the distance an individual
 *                       boundary may shift in a single rebalance.
 *
 *************************************************************************/
static void DLBnewBounds(real8 *Bold, real8 *Bnew, real8 *times, int numParts,
                         real8 shiftLimit)
{
        int   i, j;
        real8 minBound, maxBound;
        real8 wTotal, wRemaining, wTarget, *wDensity;
        real8 wDensitySum = 0.0;
/*
 *      Boundaries at lower and upper edges are static and don't
 *      move, so just copy them from the old to the new
 */
        Bnew[0]        = Bold[0] ;
        Bnew[numParts] = Bold[numParts] ;

/*
 *      Calculate the weight density for each of the old partitions
 */
        wDensity = (real8 *)malloc(numParts * sizeof(real8));

        for (i = 0; i < numParts; i++) {
            wDensity[i] = times[i] / (Bold[i+1]-Bold[i]);
            wDensitySum += wDensity[i];
        }

/*
 *      To avoid divide by zero.
 */
        if (fabs(wDensitySum) < 1.0e-8) {

            for (i = 1; i < numParts; i++) {
                Bnew[i] = Bold[i];
            }

            free(wDensity);

            return;
        }

/*
 *      For each partition, the lower boundary is already set, we
 *      just need to find the upper boundary.
 *
 *      The new partition <i> may overlap up to 3 of the old
 *      partitions; i-1, i, and i+1;
 */
        for (i = 0; i < (numParts-1); i++) {

/*
 *          Set outer limits on how far up or down the boundary
 *          may move base on the MAX_SHIFT factor and the hard
 *          limit <shiftLimit> indicated by the caller.
 */
            minBound = Bold[i+1] -
                       MIN(shiftLimit, (MAX_SHIFT * (Bold[i+1] - Bold[i])));

            maxBound = Bold[i+1] +
                       MIN(shiftLimit, (MAX_SHIFT * (Bold[i+2] - Bold[i+1])));

/*
 *          Calculate the target weight for the current partition.
 */
            wTotal = 0.0;

            for (j = 0; j < numParts+1; j++) {
                if (Bnew[i] > Bold[j]) {
                    continue;
                } else if (Bnew[i] < Bold[j]) {
                    wTotal += times[j-1] * ((Bold[j] - Bnew[i])/
                                             (Bold[j] - Bold[j-1]));
                }
                break;
            }

            for ( ; j < numParts; j++) {
                wTotal += times[j];
            }

            wTarget = wTotal / (numParts - i);

/*
 *          If the new partition overlaps the previous i-1 partition...
 */
            if (Bnew[i] < Bold[i]) {
                wTarget -= (Bold[i] - Bnew[i]) * wDensity[i-1];
                if (wTarget < 0.0) {
                    Bnew[i+1] = minBound;
                }

/*
 *              If it does not overlap the i+1 partition...
 */
                if (wTarget < times[i]) { 
                    Bnew[i+1] = Bold[i] + (wTarget / wDensity[i]);
                    if (Bnew[i+1] < minBound) {
                        Bnew[i+1] = minBound;
                    }
                }
/*
 *              Does overlap the i+1 partition...
 */
                else {
                    wTarget -= times[i];
                    Bnew[i+1] = Bold[i+1] + (wTarget / wDensity[i+1]);
                    if (Bnew[i+1] > maxBound) {
                        Bnew[i+1] = maxBound;
                    }
                }
            }
/*
 *          New partition does not overlap previous i-1 partition
 */
            else {
                wRemaining = wDensity[i] * (Bold[i+1] - Bnew[i]);
/*
 *              Does not overlap previous i+1 partition
 */
                if (wRemaining >= wTarget) {
                    Bnew[i+1] = Bnew[i] + (wTarget / wDensity[i]);
                    if (Bnew[i+1] < minBound) {
                        Bnew[i+1] = minBound;
                    }
                }
/*
 *              Overlaps previous i+1 partition
 */
                else {
                    wTarget -= wRemaining;
/*
 *                  To avoid divide by zero...
 */
                    if (wDensity[i+1] > 0.0) {
                        Bnew[i+1] = Bold[i+1] + (wTarget / wDensity[i+1]);
                    } else {
                        Bnew[i+1] = Bold[i+1];
                    }

                    if (Bnew[i+1] > maxBound) {
                        Bnew[i+1] = maxBound;
                    }
                }
            }
        }

        free(wDensity) ;

        return;
}


/*-----------------------------------------------------------------------
 *
 *      Function:     DLBalanceX
 *      Description:  Using the current load data, determine if load
 *                    imbalance warrants rebalancing (if so) and
 *                    calculate new boundaries in the X dimension.
 *
 *      Arguments:
 *          loadData  Array (1 element per domain) containing the
 *                    per task load data to be used for rebalancing.
 *                            
 *----------------------------------------------------------------------*/
void DLBalanceX(Home_t *home, real8 *loadData)
{
        int nXdoms, nYdoms, nZdoms, i, j, k, iDom;
        real8 maxVal, sumVal, avgVal, imbalance;
        real8 *valX, *Bnew;
        real8 cellSize, shiftLimit;
        RSDecomp_t *decomp;
        
        decomp = (RSDecomp_t *)home->decomp;

        nXdoms = home->param->nXdoms;
        nYdoms = home->param->nYdoms;
        nZdoms = home->param->nZdoms;
        
/*
 *      We never want a boundary to shift more than 1/2 the size
 *      of a cell in the dimension being balanced, so set an
 *      upper limit on the possible distance the bounds may shift.
 */
        cellSize = home->param->Lx / home->param->nXcells;
        shiftLimit = 0.49 * cellSize;

/*
 *      Sum the total values in each x-partition
 */
        valX = (real8 *) malloc (nXdoms * sizeof(real8));

        for (i = 0; i < nXdoms; i++) {
            valX[i] = 0.0;
            for (j = 0; j < nYdoms; j++) {
                for (k = 0; k < nZdoms; k++) {
                    iDom = EncodeDomainIdx (home, i, j, k);
                    valX[i] += loadData[iDom];
                }
            }
        }
        
/*
 *      See if load imbalance warrants a new domain decomposition
 *      in this dimension.
 */
        
        maxVal = 0.0;
        sumVal = 0.0;
        
        for (i = 0; i < nXdoms; i++) {
            sumVal += valX[i];
            if (valX[i] > maxVal) maxVal = valX[i];
        }

        avgVal = sumVal / nXdoms;
        
        if(avgVal==0) {
            imbalance=0;
            fprintf(stderr,"DLBalance: avgVal = 0 !\n");
        } else {
            imbalance = (maxVal - avgVal) / avgVal;
        }
           
/*
 *      If the imbalance is below the specified threshhold,
 *      we'll consider the problem already balanced and bail out
 */
        if (imbalance < DLBTHRESHOLD) {
            free(valX);
            return;
        }
        
/*
 *      Calculate new boundaries in the X dimension
 */
        Bnew = (real8 *)malloc((nXdoms+1) * sizeof (real8));

        DLBnewBounds(decomp->domBoundX, Bnew, valX, nXdoms, shiftLimit);
        
/*
 *      Replace old X domain boundaries with new ones
 */
        free(decomp->domBoundX);
        decomp->domBoundX = Bnew;
        
        free(valX);
        
        return;
}
        
        
/*-----------------------------------------------------------------------
 *
 *      Function:     DLBalanceY
 *      Description:  Using the current load data, determine if load
 *                    imbalance warrants rebalancing (if so) and
 *                    calculate new boundaries in the Y dimension.
 *
 *      Arguments:
 *          loadData  Array (1 element per domain) containing the
 *                    per task load data to be used for rebalancing.
 *
 *----------------------------------------------------------------------*/
void DLBalanceY(Home_t *home, real8 *loadData)
{
        int   nXdoms, nYdoms, nZdoms, i, j, k, iDom, need;
        real8 maxVal, sumVal, avgVal, imbalance;
        real8 *Bnew;
        real8 **valY;
        real8 cellSize, shiftLimit;
        RSDecomp_t *decomp;
        
        decomp = (RSDecomp_t *)home->decomp;

        nXdoms = home->param->nXdoms;
        nYdoms = home->param->nYdoms;
        nZdoms = home->param->nZdoms;
        
/*
 *      We never want a boundary to shift more than 1/2 the size
 *      of a cell in the dimension being balanced, so set an
 *      upper limit on the possible distance the bounds may shift.
 */
        cellSize = home->param->Ly / home->param->nYcells;
        shiftLimit = 0.49 * cellSize;

/*
 *      Sum the total values in each y-partition
 */
        valY = (real8 **) malloc (nXdoms * sizeof(real8 *));

        for (i = 0; i < nXdoms; i++) {
            valY[i] = (real8 *) malloc (nYdoms * sizeof(real8));
            for (j = 0; j < nYdoms; j++) {
                valY[i][j] = 0.0;
                for (k = 0; k < nZdoms; k++) {
                    iDom = EncodeDomainIdx (home, i, j, k);
                    valY[i][j] += loadData[iDom];
                }
            }
        }
        
/*
 *      See if load imbalance warrants a new domain decomposition
 *      in this dimension.  If any Y partition in any of the X
 *      partitions needs it, they all get it.
 */
        
        need = 0;

        for (i = 0; i < nXdoms; i++) {
            maxVal = 0.0;
            sumVal = 0.0;
            for (j = 0; j < nYdoms; j++) {
        
                sumVal += valY[i][j];
                if (valY[i][j] > maxVal) maxVal = valY[i][j];
        
            }
        
            avgVal = sumVal / nYdoms;
        
            if (avgVal == 0) 
                imbalance = 0;
            else
                imbalance = (maxVal - avgVal) / avgVal;
              
            if (imbalance >= DLBTHRESHOLD) need++;
        }

/*
 *      If the imbalance was below the specified threshhold,
 *      we consider the problem already balanced and bail out
 */
        if (!need) {
            for (i = 0; i < nXdoms; i++) {
                free (valY[i]);
            }
            free(valY);
            return;
        }
        
/*
 *      Calculate new Y boundaries for each X partition
 */
        Bnew = (real8 *) malloc ((nYdoms+1) * sizeof(real8));
        
        for (i = 0; i < nXdoms; i++) {
            DLBnewBounds(decomp->domBoundY[i], Bnew, valY[i],
                         nYdoms, shiftLimit);
        
/*
 *          Replace old Y domain boundaries with new ones
 */
            for (j = 0; j < nYdoms+1; j++) {
                decomp->domBoundY[i][j] = Bnew[j];
            }
        }
        

/*
 *      Free up the temporary arrays.
 */
        for (i = 0; i < nXdoms; i++) {
            free(valY[i]);
        }

        free(valY);
        free(Bnew);
        
        return;
}
        
        
/*-----------------------------------------------------------------------
 *
 *      Function:     DLBalanceZ
 *      Description:  Using the current load data, determine if load
 *                    imbalance warrants rebalancing (if so) and
 *                    calculate new boundaries in the Z dimension.
 *
 *      Arguments:
 *          loadData  Array (1 element per domain) containing the
 *                    per task load data to be used for rebalancing.
 *
 *----------------------------------------------------------------------*/
void DLBalanceZ(Home_t *home, real8 *loadData)
{
        int   nXdoms, nYdoms, nZdoms, i, j, k, need, startIdx;
        real8 maxVal, sumVal, avgVal, imbalance;
        real8 *Bnew;
        real8 ***valZ;
        real8 cellSize, shiftLimit;
        RSDecomp_t *decomp;
        
        decomp = (RSDecomp_t *)home->decomp;

        nXdoms = home->param->nXdoms;
        nYdoms = home->param->nYdoms;
        nZdoms = home->param->nZdoms;
        
/*
 *      We never want a boundary to shift more than 1/2 the size
 *      of a cell in the dimension being balanced, so set an
 *      upper limit on the possible distance the bounds may shift.
 */
        cellSize = home->param->Lz / home->param->nZcells;
        shiftLimit = 0.49 * cellSize;

/*
 *      Reorganize loadData into a 3D array, valZ, for convenience.
 *      IMPORTANT: Assumes domain # assignment in dimension order
 *      Z, Y, X
 */
        valZ = (real8 ***) malloc (nXdoms * sizeof (real8 **));

        for (i = 0; i < nXdoms; i++) {
            valZ[i] = (real8 **) malloc (nYdoms * sizeof (real8 *));
            for (j = 0; j < nYdoms; j++) {
                startIdx = EncodeDomainIdx (home, i, j, 0);
                valZ[i][j] = &loadData[startIdx];
            }
        }
        
/*
 *      See if load imbalance warrants a new domain decomposition
 *      in this dimension.  If any of the Z partition of any Y
 *      partition in any of the X partitions needs it, they all
 *      get it.
 */
        need = 0;

        for (i = 0; i < nXdoms; i++) {
            for (j = 0; j < nYdoms; j++) {
                maxVal = 0.0;
                sumVal = 0.0;

                for (k = 0; k < nZdoms; k++) {
                    sumVal += valZ[i][j][k];
                    if (valZ[i][j][k] > maxVal) maxVal = valZ[i][j][k];
                }

                avgVal = sumVal / nZdoms;
        
                if (avgVal == 0.0)
                    imbalance = 0.0;
                else
                    imbalance = (maxVal - avgVal) / avgVal;
        
                if (imbalance >= DLBTHRESHOLD) need++;
            }
        }
        
/*
 *      If the imbalance was below the specified threshhold,
 *      we consider the problem already balanced and bail out
 */
        if (!need) {
            for (i = 0; i < nXdoms; i++)  {
                free(valZ[i]);
            }
            free(valZ);
            return;
        }
        

/*
 *      Calculate new Z boundaries for each XY block
 */
        Bnew = (real8 *) malloc ((nZdoms+1) * sizeof (real8));
        
        for (i = 0; i < nXdoms; i++) {
            for (j = 0; j < nYdoms; j++) {
                DLBnewBounds(decomp->domBoundZ[i][j], Bnew, valZ[i][j],
                             nZdoms, shiftLimit);

/*
 *              Replace old Z boundaries with new ones
 */
                for (k = 0; k < nZdoms+1; k++) {
                    decomp->domBoundZ[i][j][k] = Bnew[k];
                }
            }
        }
        
/*
 *      Free up the temporary arrays.
 */
        for (i = 0; i < nXdoms; i++) {
            free(valZ[i]);
        }

        free(valZ);
        free(Bnew);
        
        return;
}


/*---------------------------------------------------------------------------
 *
 *	Function:	UniformRSDecomp
 *	Description:	Generate a new domain decomposition based solely
 *			on the domain geometry.  In essence, this will
 *			divide the problems space into a set of domains
 *			of identical size and shape.
 *
 *-------------------------------------------------------------------------*/
void UniformRSDecomp(Param_t *param, RSDecomp_t **uniDecomp)
{
        int      i, xDom, yDom, zDom;
        int      nXdoms, nYdoms, nZdoms;
        real8    Lx, Ly, Lz;
        real8    xDomLen, yDomLen, zDomLen;
        real8    *xCut, *yCut, *zCut;
        real8    *newXBounds, **newYBounds, ***newZBounds;
        RSDecomp_t *decomp;

/*
 *      Set up the boundaries for the uniform decomp
 */
        nXdoms = param->nXdoms;
        nYdoms = param->nYdoms;
        nZdoms = param->nZdoms;

        Lx = param->Lx;
        Ly = param->Ly;
        Lz = param->Lz;

        xDomLen = Lx / param->nXdoms;
        yDomLen = Ly / param->nYdoms;
        zDomLen = Lz / param->nZdoms;

        xCut = (real8 *) malloc((nXdoms+1) * sizeof(real8));
        yCut = (real8 *) malloc((nYdoms+1) * sizeof(real8));
        zCut = (real8 *) malloc((nZdoms+1) * sizeof(real8));

        for (i = 0; i < param->nXdoms; i++) {
            xCut[i] = rint(param->minSideX + i*xDomLen);
        }
        xCut[param->nXdoms] = param->maxSideX;

        for (i = 0; i < param->nYdoms; i++) {
            yCut[i] = rint(param->minSideY + i*yDomLen);
        }
        yCut[param->nYdoms] = param->maxSideY;

        for (i = 0; i < param->nZdoms; i++) {
            zCut[i] = rint(param->minSideZ + i*zDomLen);
        }
        zCut[param->nZdoms] = param->maxSideZ;

/*
 *      Create a new domain decomposition.
 */
        newXBounds = (real8 *)malloc((nXdoms + 1) * sizeof(real8));
        newYBounds = (real8 **)malloc(nXdoms * sizeof(real8 *));
        newZBounds = (real8 ***)malloc(nXdoms * sizeof(real8 **));

        for (xDom = 0; xDom < nXdoms; xDom++) {

            newXBounds[xDom] = xCut[xDom];
            newYBounds[xDom] = (real8 *)malloc((nYdoms+1)*sizeof(real8));
            newZBounds[xDom] = (real8 **)malloc(nYdoms*sizeof(real8 *));

            for (yDom = 0; yDom < nYdoms; yDom++) {

                newYBounds[xDom][yDom] = yCut[yDom];
                newZBounds[xDom][yDom] =
                        (real8 *)malloc((nZdoms+1)* sizeof(real8));

                for (zDom = 0; zDom < nZdoms; zDom++) {

                    newZBounds[xDom][yDom][zDom] = zCut[zDom];

                }
                newZBounds[xDom][yDom][zDom] = zCut[zDom];
            }
            newYBounds[xDom][yDom] = yCut[yDom];
        }
        newXBounds[xDom] = xCut[xDom];

/*
 *      Free up the temporary arrays and return the new decomposition
 *      to the caller.
 */
        free(xCut);
        free(yCut);
        free(zCut);

        decomp = (RSDecomp_t *)malloc(sizeof(RSDecomp_t));

        decomp->domBoundX = newXBounds;
        decomp->domBoundY = newYBounds;
        decomp->domBoundZ = newZBounds;

        *uniDecomp = decomp;

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:       BroadcastRSDecomp
 *      Description:    Have task zero broadcast a domain decomposition
 *                      to all other tasks.  The domain decomposition
 *                      in home->decomp will be created from this
 *                      broadcast data before control is returned
 *                      to the caller.
 *      Arguments:
 *          decomp   Pointer to decomposition to be broadcast from domain
 *                   zero to all other domains. Pointer is NULL or 
 *                   uninitialized on all tasks but zero.
 *
 *-----------------------------------------------------------------------*/
void BroadcastRSDecomp(Home_t *home, RSDecomp_t *decomp)
{
        int       numVals, bufIndex;
        int       xDom, yDom, zDom;
        int       numXDoms, numYDoms, numZDoms;
        real8     *bcastBuf;
        Param_t   *param;
        RSDecomp_t *inDecomp;

        param = home->param;

        numXDoms = param->nXdoms;
        numYDoms = param->nYdoms;
        numZDoms = param->nZdoms;

        bufIndex = 0;
        numVals = numXDoms * (1 + numYDoms * (1 + numZDoms));
        bcastBuf = (real8 *)calloc(1, numVals * sizeof(real8));

        if (home->myDomain == 0) {

            for (xDom = 0; xDom < numXDoms; xDom++) {

                bcastBuf[bufIndex++] = decomp->domBoundX[xDom];

                for (yDom = 0; yDom < numYDoms; yDom++) {

                    bcastBuf[bufIndex++] = decomp->domBoundY[xDom][yDom];

                    for (zDom = 0; zDom < numZDoms; zDom++) {

                        bcastBuf[bufIndex++] =
                                decomp->domBoundZ[xDom][yDom][zDom];
                    }
                }
            }
        }

#ifdef PARALLEL
        MPI_Bcast((real8 *)bcastBuf, numVals, MPI_DOUBLE, 0, MPI_COMM_WORLD);
#endif

/*
 *      Store the domain decomposition from task zero in the <home>
 *      structure on all domains.  Note: domain zero does this as
 *      well since it currently has the decomposition in temporary arrays.
 */
        inDecomp = (RSDecomp_t *)malloc(sizeof(RSDecomp_t));
        home->decomp = (void *)inDecomp;

        inDecomp->domBoundX =
                  (real8 *)malloc((param->nXdoms+1)*sizeof(real8));
        inDecomp->domBoundY =
                  (real8 **)malloc(param->nXdoms*sizeof(real8 *));
        inDecomp->domBoundZ =
                  (real8 ***)malloc(param->nXdoms*sizeof(real8 **));

        bufIndex = 0;

        for (xDom = 0; xDom < numXDoms; xDom++) {

            inDecomp->domBoundX[xDom] = bcastBuf[bufIndex++];

            inDecomp->domBoundY[xDom] =
                      (real8 *) malloc((numYDoms+1)*sizeof(real8));
            inDecomp->domBoundZ[xDom] =
                      (real8 **) malloc(numYDoms*sizeof(real8 *));

            for (yDom = 0; yDom < numYDoms; yDom++) {

                inDecomp->domBoundY[xDom][yDom] = bcastBuf[bufIndex++];

                inDecomp->domBoundZ[xDom][yDom] =
                          (real8 *)malloc((numZDoms+1)* sizeof(real8));

                for (zDom = 0; zDom < numZDoms; zDom++) {
                    inDecomp->domBoundZ[xDom][yDom][zDom] =
                            bcastBuf[bufIndex++];
                }

                inDecomp->domBoundZ[xDom][yDom][numZDoms] = param->maxSideZ;
            }

            inDecomp->domBoundY[xDom][numYDoms] = param->maxSideY;
        }

        inDecomp->domBoundX[numXDoms] = param->maxSideX;

        free(bcastBuf);

        return;
}
