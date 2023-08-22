/****************************************************************************
 *
 *      Module:      RBDecomp.c
 *      Description: Contains functions needed for both generating and
 *                   accessing a domain decomposition based on a recursive
 *                   bisection algorithm.
 *
 *                   The general idea with this algorithm is to
 *                   start with the entire problem space and subdivide
 *                   along the X, Y and/or Z dimensions (depending on the
 *                   number of domains per dimension) into octants,
 *                   quarters or halves such that the per-domain computational
 *                   cost of each subsection will be roughly the same.
 *                   The decomposition is then recursively applied to
 *                   each subsection until no further decomposition
 *                   is required.
 *
 *                   The decomposition lends itself nicely to a
 *                   hierarchical oct-tree, where each terminating "leaf"
 *                   in the tree (a subsection that is not decomposed
 *                   any further) represents a single domain.
 *
 *      Includes public functions:
 *          AllocRBDecomp()
 *          BroadcastRBDecomp()
 *          DomID2DecompID()
 *          ExchangeRBDecomp()
 *          FindRBDecompCoordDomain()
 *          FreeRBDecomp()
 *          GetAllRBDecompBounds()
 *          GetRBDecompLocalDomainBounds()
 *          GetRBDecompCellDomainList()
 *          PovPlotRBDecomp()
 *          RBCheckLoadBalance()
 *          RBDecomp()
 *          ReadRBDecompBounds()
 *          UniformRBDecomp()
 *          WriteRBDecompBounds()
 *          XPlotRBDecomp()
 *
 *      Includes private functions:
 *          DecompID2DomID()
 *          FindVolumeDomainList()
 *          GetBisection()
 *          GetDecompCnt()
 *          InitRBDecomp()
 *          SetAllRBDecompBounds()
 *
 ***************************************************************************/
#include <stdio.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Decomp.h"
#include "RBDecomp.h"
#include "Restart.h"
#include "DisplayC.h"

#ifdef USE_HDF
#include <hdf5.h>
#endif

#define MAX_BOUND_SHIFT_FACT 0.30

/*
 *      Each byte of a subdecomp ID identifies (at the corresponding
 *      level of the decomposition hierarchy) which octant the sub-
 *      partition is of the parent partition.  If we view the octants
 *      as a 3-dimensional structure we can identify any octant by
 *      a 3-bit value where each bit represents one of the X, Y or Z
 *      dimensions and is either a 1 or a zero.  Below we define
 *      the bit number corresponding to each dimension and the bit
 *      masks that can be used to isolate an individual one of the bits.
 */
#define X_BIT 0
#define Y_BIT 1
#define Z_BIT 2

#define X_BITMASK 0x01
#define Y_BITMASK 0x02
#define Z_BITMASK 0x04

/*
 *      Define an increment size by which the list of domains
 *      intersecting a cell will be increased when determining
 *      cell/domain intersection lists.
 */
#define DOMAIN_CNT_INCREMENT 20

static void GetBisection(real8 newMinCoord, real8 newMaxCoord, int numElements,
        real8 *load, real8 *bounds, real8 *coord, real8 minSideRatio);

static void SetAllRBDecompBounds(Param_t *param, RBDecomp_t *decomp,
        int level, real8 *globalBounds);


/*-------------------------------------------------------------------------
 *
 *      Function:    RBCheckLoadBalance
 *      Description: Use an imbalance threshold to choose the proper
 *                   level at which to rebalance the domain boundaries.
 *
 *      Arguments:
 *          decomp             Current domain decomposition
 *          loadData           Array containing the per-process load
 *                             data based on the decomposition in <decomp>
 *          imbalanceThreshold Indicates the maximum 'imbalance' permitted
 *                             at any level before rebalancing is required
 *                             for that level.  Value must be in the range
 *                             0.0 <= imbalanceThreshold <= 1.0
 *          currLevel          Current level in the decomposition hierarchy.
 *          rebalanceLevel     Location in which to return to the caller
 *                             the level at which rebalance is needed.
 *                             NOTE: the contents of this should be set
 *                             to indicate the finest decomposition level
 *                             by the caller prior to the calling this
 *                             function at the coarsest decomposition
 *                             level!
 *
 *      Returns:  The total load associated with the portion of the
 *                domain decomposition indicated by <decomp>
 *
 *------------------------------------------------------------------------*/
real8 RBCheckLoadBalance(Home_t *home, RBDecomp_t *decomp, real8 *loadData,
                         real8 imbalanceThreshold, int currLevel,
                         int *rebalanceLevel)
{
        int   octant, activeOctants;
        real8 minLoad, maxLoad, avgLoad, subpartLoad, imbalance;


        minLoad   = 1.0e+20;
        maxLoad   = 0.0;
        avgLoad   = 0.0;
        imbalance = 0.0;
        activeOctants = 0;
        decomp->totLoad = 0.0;

/*
 *      Loop through all active octants summing the total load
 *      from all octants, as well as identifying the minimum and
 *      maximum loads for any of the octants.
 */
        for (octant = 0; octant < 8; octant++) {

            if (decomp->subDecomp[octant] != (RBDecomp_t *)NULL) {

/*
 *              If the subpartition represents a single domain, just
 *              grab the load for that domain, otherwise recursively
 *              calculate the load for the subpartition.
 */
                if (decomp->subDecomp[octant]->domID >= 0) {
                    subpartLoad = loadData[decomp->subDecomp[octant]->domID];
                    decomp->subDecomp[octant]->totLoad = subpartLoad;
                } else {
                    subpartLoad = RBCheckLoadBalance(home,
                                                 decomp->subDecomp[octant],
                                                 loadData, imbalanceThreshold,
                                                 currLevel+1, rebalanceLevel);
                }

                if (subpartLoad > maxLoad) maxLoad = subpartLoad;
                if (subpartLoad < minLoad) minLoad = subpartLoad;

                decomp->totLoad += subpartLoad;
                activeOctants++;
            }
        }

/*
 *      Calculate the average load for each active octant and
 *      the load imbalance among the octants.  If the imbalance
 *      exceeds the provided threshold, select this level for
 *      rebalancing *unless* a coarser level (i.e. lower number
 *      level) has already been selected for rebalancing.
 */
        if ((avgLoad = decomp->totLoad / activeOctants) > 0.0) {
            imbalance = (maxLoad - avgLoad) / avgLoad;
        }

        if (imbalance > imbalanceThreshold) {
            *rebalanceLevel = MIN(*rebalanceLevel, currLevel);
        }

/*
 *      Return to the caller, the total load for all octants this
 *      decomposition substructure.
 */
        return(decomp->totLoad);
}


/*---------------------------------------------------------------------------
 *
 *      Function:    GetRBDecompLocalDomainBounds
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
void GetRBDecompLocalDomainBounds(Home_t *home, RBDecomp_t *decomp)
{
        int         octant, level = 0;
        RBDecomp_t  *decompPtr;
        static int  firstTime = 1;
        static char domDecompID[MAX_DECOMP_LVLS];

        if (firstTime) {
            firstTime = 0;
            DomID2DecompID(home, home->myDomain, domDecompID);
        }

/*
 *      The decompID identifies the unique branch of the decomposition
 *      tree associated with this domain, so we can just follow the
 *      branch all the way down.
 */
        decompPtr = decomp;

        while (decompPtr->domID < 0) {
            octant = domDecompID[++level];
            decompPtr = decompPtr->subDecomp[octant];
        }

        home->domXmin = decompPtr->cMin[X];
        home->domXmax = decompPtr->cMax[X];

        home->domYmin = decompPtr->cMin[Y];
        home->domYmax = decompPtr->cMax[Y];

        home->domZmin = decompPtr->cMin[Z];
        home->domZmax = decompPtr->cMax[Z];

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:    FindVolumeDomainList
 *      Description: Given the min/max coordinates bounding a volume
 *                   of the problem space, build a list of all domains
 *                   that intersect that volume.
 *                   
 *      Arguments:
 *          vMin           Minimum boundaries of the volume in question
 *          vMax           Maximum boundaries of the volume in question
 *          decomp         Pointer to current domain decomposition
 *          domainList     location in which to return the list of
 *                         intersecting domains to the caller.
 *          doaminCnt      location in which to return the number of
 *                         domains on <domainList> to the caller.
 *          domainListEnts pointer to integer containing the current
 *                         number of entries allocated in the domainList
 *                         array
 *
 *--------------------------------------------------------------------------*/
static void FindVolumeDomainList(real8 vMin[3], real8 vMax[3],
                                 RBDecomp_t *decomp, int **domainList,
                                 int *domainCnt, int *domainListEnts)
{
        int  octant;

        if (decomp == (RBDecomp_t *)NULL) {
            return;
        }

/*
 *      If this branch of the domain decomposition does not intersect
 *      the specified volume, no need to continue down this branch of
 *      the decomposition tree.
 */
        if ((decomp->cMin[X] >= vMax[X]) ||
            (decomp->cMin[Y] >= vMax[Y]) ||
            (decomp->cMin[Z] >= vMax[Z]) ||
            (decomp->cMax[X] <= vMin[X]) ||
            (decomp->cMax[Y] <= vMin[Y]) ||
            (decomp->cMax[Z] <= vMin[Z])) {
            return;
        }

/*
 *      If this is a leaf node of the decomposition tree (and hence
 *      associated with a single domain) add the domain ID to the list
 *      of domains that intersect the volume in question.  Otherwise
 *      recursively check each of the octants of this portion of the tree.
 */
        if (decomp->domID >= 0) {
            (*domainList)[*domainCnt] = decomp->domID;
            *domainCnt += 1;
            if (*domainCnt >= *domainListEnts) {
                *domainListEnts += DOMAIN_CNT_INCREMENT;
                *domainList = (int *)realloc(*domainList,
                                             *domainListEnts * sizeof(int));
            }
        } else {
            for (octant = 0; octant < 8; octant++) {
                FindVolumeDomainList(vMin, vMax, decomp->subDecomp[octant],
                                     domainList, domainCnt, domainListEnts);
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     GetRBDecompCellDomainList
 *      Description:  Find out what domains intersect the specified cell
 *                    return a pointer to the list of domains to the caller
 *                    along with a count of the number of domains on the
 *                    returned list.
 *
 *      Arguments:
 *          cellID     ID of the cell as returned by EncodeCellIdx()
 *          domainCnt  location in which to return the count of domains
 *                     intersecting the specified celll
 *          domainList location in which to return an array containing
 *                     the IDs of all domains intersecting the specified
 *                     cell.
 *
 *--------------------------------------------------------------------------*/
void GetRBDecompCellDomainList(Home_t *home, int cellID, int *domainCnt,
                               int **domainList)
{
        int        domainListEnts = 0;
        int        xCell, yCell, zCell;
        real8      xCellSize, yCellSize, zCellSize;
        real8      vMin[3], vMax[3];
        RBDecomp_t *decomp;
        Param_t    *param;
        Cell_t     *cell;

        param  = home->param;
        decomp = (RBDecomp_t *)home->decomp;

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
 *      Get the min and max coordinates for this cell...
 *
 *      WARNING: The cell min and max coordinate values must be
 *      computed in the same manner or neighboring domains can end
 *      up with slightly different locations for a shared cell
 *      boundary (due to numeric round-off.)
 */
        xCellSize = (param->maxSideX - param->minSideX) / param->nXcells;
        yCellSize = (param->maxSideY - param->minSideY) / param->nYcells;
        zCellSize = (param->maxSideZ - param->minSideZ) / param->nZcells;

        vMin[X] = rint((param->minSideX) + ((xCell-1) * xCellSize));
        vMin[Y] = rint((param->minSideY) + ((yCell-1) * yCellSize));
        vMin[Z] = rint((param->minSideZ) + ((zCell-1) * zCellSize));

        vMax[X] = rint((param->minSideX) + (xCell * xCellSize));
        vMax[Y] = rint((param->minSideY) + (yCell * yCellSize));
        vMax[Z] = rint((param->minSideZ) + (zCell * zCellSize));

/*
 *      Allocate a small initial array to hold the list of
 *      intersecting domains.  The size of this array will be
 *      increased as necessary if the domain count exceeds
 *      the array size during the search below.
 */
        *domainCnt = 0;
        domainListEnts = DOMAIN_CNT_INCREMENT;
        *domainList = (int *)malloc(domainListEnts * sizeof(int));

/*
 *      Recursively search the decomposition tree to find all
 *      domains intesecting the volume encompassed by the cell.
 */
        FindVolumeDomainList(vMin, vMax, decomp, domainList, domainCnt,
                             &domainListEnts);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       GetDecompCnt
 *      Description:    Used for calculating the total number of RBDecomp_t
 *                      structures that would be required to contain
 *                      the complete domain decomposition (all levels)
 *                      for the current domain count and geometry.
 *
 *      Arguments:
 *          xNumDoms,
 *          yNumDoms,
 *          zNumDoms    Number of domains (in each dimension) that are
 *                      assigned to this protion of the decomposition.
 *          totDecomps  Pointer to count of total decomposition structures
 *                      needed.  The value at this location will be
 *                      incremented with the number of structures needed
 *                      for this portion of the decompsition before
 *                      return to the caller.
 *
 *-------------------------------------------------------------------------*/
static void GetDecompCnt
(
   const    int  xNumDoms  ,   ///< number of domains (x)
   const    int  yNumDoms  ,   ///< number of domains (y)
   const    int  zNumDoms  ,   ///< number of domains (z)
   volatile int *totDecomps    ///< resulting domain count
)
{
        int  i;
        int  xDoms[2]={0,0}, yDoms[2]={0,0}, zDoms[2]={0,0};
        int  cut[3], subDecompUsed[8];

        *totDecomps += 1;

/*
 *      Figure out the dimensions (if any) in which this portion of the
 *      deocmposition must be further decomposed.
 */
        cut[X] = (xNumDoms > 1);
        cut[Y] = (yNumDoms > 1);
        cut[Z] = (zNumDoms > 1);

        if (cut[X] + cut[Y] + cut[Z] == 0) {
            return;
        }

        for (i = 0; i < 8; i++) {
            subDecompUsed[i] = 1;
        }

/*
 *      If we're not cutting in the X dimension, none of the
 *      right side octants will be used.  Also set the number of
 *      domains in the X dimension that will be encompassed by
 *      octants on each side of the X bisection (left and right
 *      side octants).
 *
 *      Then the same sort of thing for Y and Z dimensions.
 */
        if (!cut[X]) {
            subDecompUsed[LRF] = 0;           
            subDecompUsed[URF] = 0;           
            subDecompUsed[LRB] = 0;           
            subDecompUsed[URB] = 0;           
            xDoms[0] = xNumDoms;
        } else {
            xDoms[0] = xNumDoms >> 1;
            xDoms[1] = xNumDoms - xDoms[0];
        }

        if (!cut[Y]) {
            subDecompUsed[ULF] = 0;           
            subDecompUsed[URF] = 0;           
            subDecompUsed[ULB] = 0;           
            subDecompUsed[URB] = 0;           
            yDoms[0] = yNumDoms;
        } else {
            yDoms[0] = yNumDoms >> 1;
            yDoms[1] = yNumDoms - yDoms[0];
        }

        if (!cut[Z]) {
            subDecompUsed[LLB] = 0;           
            subDecompUsed[LRB] = 0;           
            subDecompUsed[ULB] = 0;           
            subDecompUsed[URB] = 0;           
            zDoms[0] = zNumDoms;
        } else {
            zDoms[0] = zNumDoms >> 1;
            zDoms[1] = zNumDoms - zDoms[0];
        }

/*
 *      Loop through all the octants that will be populated, incrementing
 *      the total decomp count with the number of portions into which the
 *      octant is sub-partitioned.
 */
        for (i = 0; i < 8; i++) {

            if (!subDecompUsed[i]) {
                continue;
            }

            switch (i) {
                case LLF: { GetDecompCnt(xDoms[0], yDoms[0], zDoms[0], totDecomps); break; }
                case LRF: { GetDecompCnt(xDoms[1], yDoms[0], zDoms[0], totDecomps); break; }
                case ULF: { GetDecompCnt(xDoms[0], yDoms[1], zDoms[0], totDecomps); break; }
                case URF: { GetDecompCnt(xDoms[1], yDoms[1], zDoms[0], totDecomps); break; }
                case LLB: { GetDecompCnt(xDoms[0], yDoms[0], zDoms[1], totDecomps); break; }
                case LRB: { GetDecompCnt(xDoms[1], yDoms[0], zDoms[1], totDecomps); break; }
                case ULB: { GetDecompCnt(xDoms[0], yDoms[1], zDoms[1], totDecomps); break; }
                case URB: { GetDecompCnt(xDoms[1], yDoms[1], zDoms[1], totDecomps); break; }
            }
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     GetAllRBDecompBounds
 *      Description:  Build a simple array containing the boundaries of all
 *                    domains.  Array will contain 6 values (min boundaries
 *                    and max boundaries) for every domain.
 *
 *      Arguments:
 *          decomp     pointer to the current domain decomposition.
 *          boundsBuf  Array into which to store the boundary data.
 *                     This array must be large enough to hold 6 real8
 *                     values per domain.
 *
 *--------------------------------------------------------------------------*/
void GetAllRBDecompBounds(Home_t *home, RBDecomp_t *decomp, real8 *boundsBuf)
{
        int        index, offset;
        static int totDecompCnt = -1;

        if (totDecompCnt < 0) {
            totDecompCnt = 0;
            GetDecompCnt(home->param->nXdoms, home->param->nYdoms,
                         home->param->nZdoms, &totDecompCnt);
        }

        for (index = 0; index < totDecompCnt; index++) {
            if (decomp[index].domID >= 0) {
                offset = decomp[index].domID * 6;
                boundsBuf[offset    ] = decomp[index].cMin[X];
                boundsBuf[offset + 1] = decomp[index].cMin[Y];
                boundsBuf[offset + 2] = decomp[index].cMin[Z];
                boundsBuf[offset + 3] = decomp[index].cMax[X];
                boundsBuf[offset + 4] = decomp[index].cMax[Y];
                boundsBuf[offset + 5] = decomp[index].cMax[Z];
            }
        }

        return;
}


/*------------------------------------------------------------------------
 *
 *      Function:       BroadcastRBDecomp
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
void BroadcastRBDecomp(Home_t *home, RBDecomp_t *decomp)
{
        int        numValues, level;
        real8      *boundsBuf;
        RBDecomp_t *inDecomp;

        level = 0;
        numValues = 6 * home->numDomains;
        boundsBuf = (real8 *)malloc(numValues * sizeof(real8));

/*
 *      Have domain 0 pack all the boundaries into a buffer
 *      for transmission
 */
        if (home->myDomain == 0) {
            GetAllRBDecompBounds(home, decomp, boundsBuf);
        }

#ifdef PARALLEL
        MPI_Bcast((real8 *)boundsBuf, numValues, MPI_DOUBLE,
                  0, MPI_COMM_WORLD);
#endif

/*
 *      All domains create a new domain decomposition tree and 
 *      initialize it using the transmitted domain boundaries.
 */
        AllocRBDecomp(home, (RBDecomp_t *)NULL, &inDecomp, ALLOC_NEW_DECOMP);
        SetAllRBDecompBounds(home->param, inDecomp, level, boundsBuf);
        free(boundsBuf);

        home->decomp = (void *)inDecomp;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReadRBDecompBounds
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
void ReadRBDecompBounds(Home_t *home, void **filePtr, int numXDoms,
                        int numYDoms, int numZDoms, int saveDecomp,
                        RBDecomp_t **oldDecomp)
{
        int   i, domID, offset, numDomains;
        int   level = 0;
        real8 xMin, yMin, zMin;
        real8 xMax, yMax, zMax;
        real8 *boundsBuf;
        char  inLine[256];
        FILE  *fp = (FILE *)*filePtr;

/*
 *      Number of domains in the decomposition the restart file
 *      may not be the same as the current domain count, so
 *      we need to be sure to read the proper number of entries
 *      from the restart file.
 */
        numDomains = numXDoms * numYDoms * numZDoms;

/*
 *      If we need to return the old domain decomposition from the
 *      restart file to the caller, pre-allocate the decomposition
 *      tree and an array to temporarily store the domain boundaries
 */
        if (saveDecomp) {
            AllocRBDecomp(home, (RBDecomp_t *)NULL, oldDecomp,
                          ALLOC_NEW_DECOMP);
            boundsBuf = (real8 *)malloc(6 * numDomains * sizeof(real8));
        }

/*
 *      Read all domain boundaries and store them in the temporary buffer
 */
        for (i = 0; i < numDomains; i++) {
            Getline(inLine, sizeof(inLine), fp);
            if (saveDecomp) {
                sscanf(inLine, "%d %lf %lf %lf %lf %lf %lf",
                       &domID, &xMin, &yMin, &zMin, &xMax, &yMax, &zMax);
                offset = domID * 6;
                boundsBuf[offset    ] = xMin;
                boundsBuf[offset + 1] = yMin;
                boundsBuf[offset + 2] = zMin;
                boundsBuf[offset + 3] = xMax;
                boundsBuf[offset + 4] = yMax;
                boundsBuf[offset + 5] = zMax;
            }
        }

/*
 *      If necessary, use the array of domain boundaries to
 *      initialize all the boundaries in the decomposition tree.
 */
        if (saveDecomp) {
            SetAllRBDecompBounds(home->param, *oldDecomp, level, boundsBuf);
            free(boundsBuf);
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    ReadBinRBDecompBounds
 *      Description: Read the domain decomposition (if any) from the HDF5
 *                   data file and return it to the caller.
 *
 *      Arguments:
 *          filePtr     File pointer to the opened HDF5 data file;
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
void ReadBinRBDecompBounds(Home_t *home, void *filePtr, int numXDoms,
                           int numYDoms, int numZDoms, RBDecomp_t **oldDecomp)
{
#ifdef USE_HDF
        int    status, numDomains;
        int    level = 0;
        real8  *boundsBuf;
        hid_t  *fileID = (hid_t *)filePtr;

        numDomains = numXDoms * numYDoms * numZDoms;

        AllocRBDecomp(home, (RBDecomp_t *)NULL, oldDecomp, ALLOC_NEW_DECOMP);
        boundsBuf = (real8 *)malloc(6 * numDomains * sizeof(real8));

        status = ReadHDFDataset(*fileID, (char *) "/decomposition", H5T_NATIVE_DOUBLE, numDomains * 6, boundsBuf);

        if (status != 0) {
            Fatal("ReadBinRBDecompBounds: Error reading decomposition");
        }

        SetAllRBDecompBounds(home->param, *oldDecomp, level, boundsBuf);
        free(boundsBuf);
#endif
        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    WriteRBDecompBounds
 *      Description: Writes the domain geometry and boundaries for
 *                   the specified decomposition into the restart file.
 *
 *      Arguments:
 *          fp          open file descriptor for the restart file being
 *                      written
 *          decomp      pointer to current domain decomposition.
 *          level       current level of domain decomposition being written.
 *
 *------------------------------------------------------------------------*/
void WriteRBDecompBounds(Home_t *home, FILE *fp, RBDecomp_t *decomp,
                         int level)
{
        int        index;
        static int totDecompCnt = -1;

/*
 *      Only when we're at the top level... write the current domain
 *      geometry and decomposition type to the restart file
 */
        if (level == 0) {
            fprintf(fp, "# Dom_ID  Minimum XYZ bounds   Maximum XYZ bounds\n");
            if (home->numDomains == 0) {
                return;
            }
        }

        if (decomp == (RBDecomp_t *)NULL) {
            return;
        }

        if (totDecompCnt < 0) {
            totDecompCnt = 0;
            GetDecompCnt(home->param->nXdoms, home->param->nYdoms,
                         home->param->nZdoms, &totDecompCnt);
        }

/*
 *      Loop through all the structures of the decomposition tree.  If
 *      the structure corresponds to a leaf node of the tree, print
 *      the associated domain boundaries.
 */
        for (index = 0; index < totDecompCnt; index++) {
            if (decomp[index].domID >= 0) {
                fprintf(fp,"  %d  %11.4f %11.4f %11.4f  %11.4f %11.4f %11.4f\n",
                        decomp[index].domID, decomp[index].cMin[X],
                        decomp[index].cMin[Y], decomp[index].cMin[Z],
                        decomp[index].cMax[X], decomp[index].cMax[Y],
                        decomp[index].cMax[Z]);
            }
        }

        return;
}


#ifndef NO_XWINDOW
/*---------------------------------------------------------------------------
 *
 *      Function:     XPlotRSDecomp
 *      Description:  Plots the domain boundaries for a recursive bisection
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
void XPlotRBDecomp(Home_t *home, RBDecomp_t *decomp, real8 xMin,
                   real8 yMin, real8 zMin, real8 lMax, int color,
                   real8 lineWidth)
{
        int        index;
        real8      x1, x2, y1, y2, z1, z2;
        static int totDecompCnt = -1;


        if (decomp == (RBDecomp_t *)NULL) {
            return;
        }

        if (totDecompCnt < 0) {
            totDecompCnt = 0;
            GetDecompCnt(home->param->nXdoms, home->param->nYdoms,
                         home->param->nZdoms, &totDecompCnt);
        }

        for (index = 0; index < totDecompCnt; index++) {
/*
 *          If this is a base decomposition for a single domain, plot
 *          the boundaries.
 *
 *          Note: For plotting purposes, positions for the endpoints
 *          of the lines are scaled from absolute positions to
 *          values in the range -1 to +1.
 */
            if (decomp[index].domID >= 0) {

                x2 = (decomp[index].cMax[X] - xMin) / lMax * 2 - 1;
                x1 = (decomp[index].cMin[X] - xMin) / lMax * 2 - 1;

                y2 = (decomp[index].cMax[Y] - yMin) / lMax * 2 - 1;
                y1 = (decomp[index].cMin[Y] - yMin) / lMax * 2 - 1;

                z2 = (decomp[index].cMax[Z] - zMin) / lMax * 2 - 1;
                z1 = (decomp[index].cMin[Z] - zMin) / lMax * 2 - 1;

                WinDrawLine(x1, y1, z1, x2, y1, z1, color, lineWidth/2, 0);
                WinDrawLine(x2, y1, z1, x2, y2, z1, color, lineWidth/2, 0);
                WinDrawLine(x2, y2, z1, x1, y2, z1, color, lineWidth/2, 0);
                WinDrawLine(x1, y2, z1, x1, y1, z1, color, lineWidth/2, 0);

                WinDrawLine(x1, y1, z1, x1, y1, z2, color, lineWidth/2, 0);
                WinDrawLine(x2, y1, z1, x2, y1, z2, color, lineWidth/2, 0);
                WinDrawLine(x1, y2, z1, x1, y2, z2, color, lineWidth/2, 0);
                WinDrawLine(x2, y2, z1, x2, y2, z2, color, lineWidth/2, 0);

                WinDrawLine(x1, y1, z2, x2, y1, z2, color, lineWidth/2, 0);
                WinDrawLine(x2, y1, z2, x2, y2, z2, color, lineWidth/2, 0);
                WinDrawLine(x2, y2, z2, x1, y2, z2, color, lineWidth/2, 0);
                WinDrawLine(x1, y2, z2, x1, y1, z2, color, lineWidth/2, 0);

            }
        }

        return;
}
#endif


/*-------------------------------------------------------------------------
 *
 *      Function:    FreeRBDecomp
 *      Description: Free all arrays associated with the recursive
 *                   bisection domain decomposition 
 *
 *                   Note: At this time, all the RBDecomp_t structures
 *                   required for a single decomposition tree are allocated
 *                   as a single block and the highest level of the
 *                   decomposition pointed to the beginning of the block. 
 *                   So, we can do a single free() to release the whole
 *                   structure.
 *
 *      Arguments:
 *          decomp   pointer to the domain decomposition to be freed
 *
 *------------------------------------------------------------------------*/
void FreeRBDecomp(RBDecomp_t *decomp)
{
        free(decomp);
}


/*---------------------------------------------------------------------------
 *
 *      Function:       FindRBDecompCoordDomain
 *      Description:    Determines the ID of the domain which contains
 *                      the specified coordinate.
 *
 *                      A coordinate directly on a domain boundary
 *                      will be associated with the first domain
 *                      found which shares that boundary.
 *
 *      Returns:   -1 if the coordinate is outside the problem space;
 *                 otherwise returns the ID of the domain owning the
 *                 specified coordinate.
 *
 *-------------------------------------------------------------------------*/
int FindRBDecompCoordDomain(RBDecomp_t *decomp, real8 x, real8 y, real8 z)
{
        int i, domID;

        if (decomp == (RBDecomp_t *)NULL) {
            return(-1);
        }

/*
 *      If the specified coordinate is not contained within
 *      this portion of the decomposition, just return without
 *      any additional checks.
 */
        if ((x <  decomp->cMin[X]) || (x > decomp->cMax[X]) ||
            (y <  decomp->cMin[Y]) || (y > decomp->cMax[Y]) ||
            (z <  decomp->cMin[Z]) || (z > decomp->cMax[Z])) {
            return(-1);
        }

/*
 *      Coordinates are contained within this portion of the
 *      the decomposition.  If this portiton is associated with
 *      a single domain, just return the domain ID, otherwise
 *      recursively search the octants of this partition to find
 *      the proper domain.
 */
        if (decomp->domID >= 0) {
            return(decomp->domID);
        }

        for (i = 0; i < 8; i++) {
            domID = FindRBDecompCoordDomain(decomp->subDecomp[i], x, y, z);
            if (domID >= 0) break;
        }

        return(domID);
}


/*-------------------------------------------------------------------------
 *
 *      Function:       DecompID2DomID
 *      Description:    Calculate the domain ID associated with 
 *                      the specified decomp ID.
 *
 *      Arguments:
 *          decompID    Character array containing the decomp ID
 *                      to be converted.  Each byte of the decompID
 *                      identifies the octant of the domposition
 *                      encompassing the domain at a particular level
 * 
 *      Returns:    The domain ID associated with the specified
 *                  decomp ID.
 *
 *-----------------------------------------------------------------------*/
static int DecompID2DomID(Home_t *home, char *decompID)
{
        int x, y, z, domID, octant, level;
        int xMaxLevel, yMaxLevel, zMaxLevel;
        int tmpDomCnt, domsOnLeft, domsOnRight;

        x = 0;
        y = 0;
        z = 0;

        xMaxLevel = home->xMaxLevel;
        yMaxLevel = home->yMaxLevel;
        zMaxLevel = home->zMaxLevel;

/*
 *      Find the domain's index in the X dimension
 */
        tmpDomCnt = home->param->nXdoms;

        for (level = 1; level <= xMaxLevel; level++) {
            domsOnLeft  = tmpDomCnt > 1 ? tmpDomCnt / 2 : 1;
            domsOnRight = tmpDomCnt - domsOnLeft;
            octant = decompID[level];
            if ((octant & X_BITMASK) >> X_BIT) {
                x += domsOnLeft;
                tmpDomCnt = domsOnRight;
            } else {
                tmpDomCnt = domsOnLeft;
            }
        }

/*
 *      Find the domain's index in the Y dimension
 */
        tmpDomCnt = home->param->nYdoms;

        for (level = 1; level <= yMaxLevel; level++) {
            domsOnLeft  = tmpDomCnt > 1 ? tmpDomCnt / 2 : 1;
            domsOnRight = tmpDomCnt - domsOnLeft;
            octant = decompID[level];
            if ((octant & Y_BITMASK) >> Y_BIT) {
                y += domsOnLeft;
                tmpDomCnt = domsOnRight;
            } else {
                tmpDomCnt = domsOnLeft;
            }
        }

/*
 *      Find the domain's index in the Z dimension
 */
        tmpDomCnt = home->param->nZdoms;

        for (level = 1; level <= zMaxLevel; level++) {
            domsOnLeft  = tmpDomCnt > 1 ? tmpDomCnt / 2 : 1;
            domsOnRight = tmpDomCnt - domsOnLeft;
            octant = decompID[level];
            if ((octant & Z_BITMASK) >> Z_BIT) {
                z += domsOnLeft;
                tmpDomCnt = domsOnRight;
            } else {
                tmpDomCnt = domsOnLeft;
            }
        }

        domID = EncodeDomainIdx(home, x, y, z);

        return(domID);
}


/*-------------------------------------------------------------------------
 *
 *      Function:       DomID2DecompID
 *      Description:    Calculate the decomp ID associated with 
 *                      the specified domain.
 *
 *      Arguments:
 *          domID       Domain ID to be converted.
 *          decompID    Character array in which the decomp ID
 *                      associated with the specified domain ID will
 *                      be returned to the caller.  Each byte of the
 *                      decompID identifies the octant of the domposition
 *                      encompassing the domain at a particular level
 * 
 *-----------------------------------------------------------------------*/
void DomID2DecompID(Home_t *home, int domID, char decompID[MAX_DECOMP_LVLS])
{
        int x, y, z, level;
        int xMaxLevel, yMaxLevel, zMaxLevel;
        int minDom, midDom, maxDom;
        int rightDomCount, leftDomCount, portionDomCount;

        xMaxLevel = home->xMaxLevel;
        yMaxLevel = home->yMaxLevel;
        zMaxLevel = home->zMaxLevel;

        DecodeDomainIdx(home, domID, &x, &y, &z);

        memset(decompID, 0, MAX_DECOMP_LVLS);

        minDom = 0;
        maxDom = home->param->nXdoms - 1;

        for (level = 1; level <= xMaxLevel; level++) {

            portionDomCount = maxDom - minDom + 1;
            rightDomCount = portionDomCount - (portionDomCount >> 1);
            leftDomCount  = portionDomCount - rightDomCount;
            midDom = minDom + (leftDomCount - 1);
            midDom = MIN(MAX(minDom, midDom), maxDom);
            
            if (x <= midDom) {
                maxDom = midDom;
            } else {
                decompID[level] |= 1 << X_BIT;
                minDom = midDom + 1;
            }
        }

        minDom = 0;
        maxDom = home->param->nYdoms - 1;

        for (level = 1; level <= yMaxLevel; level++) {

            portionDomCount = maxDom - minDom + 1;
            rightDomCount = portionDomCount - (portionDomCount >> 1);
            leftDomCount  = portionDomCount - rightDomCount;
            midDom = minDom + (leftDomCount - 1);
            midDom = MIN(MAX(minDom, midDom), maxDom);
            
            if (y <= midDom) {
                maxDom = midDom;
            } else {
                decompID[level] |= 1 << Y_BIT;
                minDom = midDom + 1;
            }
        }

        minDom = 0;
        maxDom = home->param->nZdoms - 1;

        for (level = 1; level <= zMaxLevel; level++) {

            portionDomCount = maxDom - minDom + 1;
            rightDomCount = portionDomCount - (portionDomCount >> 1);
            leftDomCount  = portionDomCount - rightDomCount;
            midDom = minDom + (leftDomCount - 1);
            midDom = MIN(MAX(minDom, midDom), maxDom);
            
            if (z <= midDom) {
                maxDom = midDom;
            } else {
                decompID[level] |= 1 << Z_BIT;
                minDom = midDom + 1;
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       GetBisection
 *      Description:    Calculate the coordinate of the plane that
 *                      bisects a bounded volume such that the ratio
 *                      of the load on the minimum side of the bisection
 *                      plane to the total load is <minSideratio>.
 *
 *      Arguments:
 *          newMinCoord
 *          newMaxCoord
 *          numElements   Number of values in the <load> array.
 *          load          Array of <numElements> values containing load data
 *                        associated with volumes defined by <bounds>
 *          bounds        Array of <numElements>+1 values defining volumes
 *                        for which there is load data.
 *          coord         Location in which to return to the caller the
 *                        coordinate of the bisecting plane.
 *          minSideRatio  portion of load to be accumulated on the minimum
 *                        side of the bisection coordinate.  Must be 
 *                        between 0.0 and 1.0.
 * 
 *-----------------------------------------------------------------------*/
void GetBisection(real8 newMinCoord, real8 newMaxCoord, int numElements,
                  real8 *load, real8 *bounds, real8 *coord, real8 minSideRatio)
{
        int   i;
        real8 neededLoad, tmpCoord, ratio;

        neededLoad = 0.0;

        for (i = 0; i < numElements; i++) {
            neededLoad += load[i];
        }

        neededLoad = neededLoad * minSideRatio;

/*
 *      Just a quick sanity check... 
 */
        if (neededLoad == 0.0) {
            *coord = bounds[0] + 0.5 * (bounds[numElements] - bounds[0]);
            return;
        }

        for (i = 0; i < numElements; i++) {
            if (neededLoad >= load[i]) {
                neededLoad -= load[i];
                tmpCoord = bounds[i+1];
            } else {
                tmpCoord = bounds[i] + (bounds[i+1] - bounds[i]) *
                                       (neededLoad / load[i]);
                break;
            }
        }

        ratio = (tmpCoord - bounds[0]) / (bounds[numElements] - bounds[0]);
        *coord = newMinCoord + (ratio * (newMaxCoord - newMinCoord));

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:       SetAllRBDecompBounds
 *      Description:    Set the boundaries for all remote domains and
 *                      all intermediate levels of the decomposition
 *                      tree based on the boundaries communicated
 *                      from the remote domains.
 *
 *      Arguments:
 *          decomp        Pointer to the highest level decomposition in the
 *                        decomposition tree.
 *          level         Current level in the decomposition tree.
 *          globalBounds  Array of domain boundaries (6 values per domain).
 *
 *-----------------------------------------------------------------------*/
static void SetAllRBDecompBounds(Param_t *param, RBDecomp_t *decomp,
                                 int level, real8 *globalBounds)
{
        int         octant, offset;
        RBDecomp_t *subpart;

        if (decomp == (RBDecomp_t *)NULL) {
            return;
        }

/*
 *      If this portion of the decomposition is associated with a single
 *      domain, set the boundaries.  Otherwise, recursively set the
 *      boundaries for each octant of this portion of the decomposition
 *      then set the boundaries for this intermediate level based on
 *      the boundaries of the sub-portions.
 */
        if (decomp->domID >= 0) {
            offset = decomp->domID * 6;
            decomp->cMin[X] = globalBounds[offset    ];
            decomp->cMin[Y] = globalBounds[offset + 1];
            decomp->cMin[Z] = globalBounds[offset + 2];
            decomp->cMax[X] = globalBounds[offset + 3];
            decomp->cMax[Y] = globalBounds[offset + 4];
            decomp->cMax[Z] = globalBounds[offset + 5];
        } else {

/*
 *          Yes, we do want to set the min to the max and the
 *          max to the min here!
 */
            decomp->cMin[X] = param->maxSideX;
            decomp->cMin[Y] = param->maxSideY;
            decomp->cMin[Z] = param->maxSideZ;

            decomp->cMax[X] = param->minSideX;
            decomp->cMax[Y] = param->minSideY;
            decomp->cMax[Z] = param->minSideZ;

            for (octant = 0; octant < 8; octant++) {

                subpart = decomp->subDecomp[octant];

                if (subpart == (RBDecomp_t *)NULL) {
                    continue;
                }

                SetAllRBDecompBounds(param, subpart, level+1, globalBounds);

                decomp->cMin[X] = MIN(decomp->cMin[X], subpart->cMin[X]);
                decomp->cMin[Y] = MIN(decomp->cMin[Y], subpart->cMin[Y]);
                decomp->cMin[Z] = MIN(decomp->cMin[Z], subpart->cMin[Z]);
                decomp->cMax[X] = MAX(decomp->cMax[X], subpart->cMax[X]);
                decomp->cMax[Y] = MAX(decomp->cMax[Y], subpart->cMax[Y]);
                decomp->cMax[Z] = MAX(decomp->cMax[Z], subpart->cMax[Z]);
            }
        }
       
        return;
}


#ifdef PARALLEL
/*---------------------------------------------------------------------------
 *
 *      Function:       ExchangeRBDecomp
 *      Description:    Distributes the locally calculated domain
 *                      boundaries among all the processors.
 *
 *      Arguments:
 *          decomp  Pointer to highest level of the decomposition tree.
 *
 *-------------------------------------------------------------------------*/
void ExchangeRBDecomp(Home_t *home, RBDecomp_t *decomp)
{
        int         numLocalValues, numGlobalValues, level, allocSize;
        real8       *localBounds, *globalBounds;
        Param_t     *param;

        param = home->param;

        numLocalValues = 6;
        numGlobalValues = home->numDomains * numLocalValues;
        allocSize = numGlobalValues * sizeof(real8);

        localBounds  = (real8 *)calloc(1, numLocalValues * sizeof(real8));
        globalBounds = (real8 *)calloc(1, allocSize);

/*
 *      Each domain sets its own domain boundaries in the array,
 *      with all other boundaries zero'ed.  When we sum all the
 *      boundary data with an MPI allreduce, each domain will
 *      have an array of all domain bounds.
 */
        localBounds[0] = home->domXmin;
        localBounds[1] = home->domYmin;
        localBounds[2] = home->domZmin;
        localBounds[3] = home->domXmax;
        localBounds[4] = home->domYmax;
        localBounds[5] = home->domZmax;

        MPI_Allgather(localBounds, numLocalValues, MPI_DOUBLE,
                      globalBounds, numLocalValues, MPI_DOUBLE,
                      MPI_COMM_WORLD);

/*
 *      Use the array of domain boundaries to set the boundaries in
 *      the decomposition tree for all the remote domains (and all
 *      intermediate levels of the decomposition tree).
 */
        level = 0;

        SetAllRBDecompBounds(param, decomp, level, globalBounds);

        free(localBounds);
        free(globalBounds);

        return;
}
#endif


/*---------------------------------------------------------------------------
 *
 *      Function:       InitRBDecomp
 *      Description:    Creates the base decomposition tree appropriate to
 *                      the current domain count and geometry.  On return
 *                      to the caller, the tree will have been built and
 *                      all items except the boundaries will have been
 *                      initialized.  Boundary initialization must be
 *                      done by a subsequent call to either RBDecomp() or
 *                      RBUniDecomp().
 *
 *      Arguments:
 *          decomp        Pointer to highest level of the decomposition tree.
 *          decompArray   Array of RBDecomp_t structures available for
 *                        use.
 *          decompsUsed   Pointer to number of RBDecomp_t structures
 *                        in <decompArray> that have already been assigned.
 *          level         Indicates the current level in the decomposition
 *                        tree.
 *
 *-------------------------------------------------------------------------*/
static void InitRBDecomp(Home_t *home, RBDecomp_t *decomp,
                         RBDecomp_t *decompArray, int *decompsUsed,
                         int level)
{
        int        i, cut[3], subDecompUsed[8];
        int        xDoms[2]={0,0}, yDoms[2]={0,0}, zDoms[2]={0,0};
        RBDecomp_t *subDecomp;
        Param_t    *param;

        param = home->param;

/*
 *      At top of the decomposition tree, explicitly set the
 *      boundaries to the edges of the problem space and
 *      indicate all domains are included in the decomposition.
 */
        if (level == 0) {
            decomp->cMin[X] = param->minSideX;
            decomp->cMin[Y] = param->minSideY;
            decomp->cMin[Z] = param->minSideZ;
            decomp->cMax[X] = param->maxSideX;
            decomp->cMax[Y] = param->maxSideY;
            decomp->cMax[Z] = param->maxSideZ;

            decomp->numDoms[X] = param->nXdoms;
            decomp->numDoms[Y] = param->nYdoms;
            decomp->numDoms[Z] = param->nZdoms;
        }

/*
 *      If this portion of the deomposition encompasses more than 1 domain
 *      in a dimension, it must be further decomposed in that dimension.
 *      Sub-decomposition will be into octants, quadrants halves or none
 *      at all.
 */
        cut[X] = (decomp->numDoms[X] > 1);
        cut[Y] = (decomp->numDoms[Y] > 1);
        cut[Z] = (decomp->numDoms[Z] > 1);

/*
 *      If there is no need to further decompose this portion, we're
 *      down to a leaf node of the tree associated with an individual
 *      domain, so set the associated domain ID and return.
 */
        if ((cut[X] + cut[Y] + cut[Z]) == 0) {
            decomp->domID = DecompID2DomID(home, decomp->decompID);
            return;
        }

        decomp->domID = -1;

/*
 *      To start off, assume the current branch of the tree will be
 *      further decomposed.  If that is not the case, these values will
 *      be reset as necessary.
 */
        for (i = 0; i < 8; i++) {
            subDecompUsed[i] = 1;
        }

/*
 *      If we're not cutting in the X dimension, none of the
 *      right side octants will be used.  Also set the number of
 *      domains in the X dimension that will be encompassed by
 *      octants on each side of the X bisection (left and right
 *      side octants).
 */
        if (!cut[X]) {
            subDecompUsed[LRF] = 0;           
            subDecompUsed[URF] = 0;           
            subDecompUsed[LRB] = 0;           
            subDecompUsed[URB] = 0;           
            xDoms[0] = decomp->numDoms[X];
        } else {
            xDoms[0] = decomp->numDoms[X] >> 1;
            xDoms[1] = decomp->numDoms[X] - xDoms[0];
        }

/*
 *      If we're not cutting in the Y dimension, none of the
 *      upper octants will be used.  Also set the number of
 *      domains in the Y dimension that will be encompassed by
 *      octants on each side of the Y bisection (lower and upper
 *      octants).
 */
        if (!cut[Y]) {
            subDecompUsed[ULF] = 0;           
            subDecompUsed[URF] = 0;           
            subDecompUsed[ULB] = 0;           
            subDecompUsed[URB] = 0;           
            yDoms[0] = decomp->numDoms[Y];
        } else {
            yDoms[0] = decomp->numDoms[Y] >> 1;
            yDoms[1] = decomp->numDoms[Y] - yDoms[0];
        }

/*
 *      If we're not cutting in the Z dimension, none of the
 *      back octants will be used.  Also set the number of
 *      domains in the Z dimension that will be encompassed by
 *      octants on each side of the Z bisection (front and back
 *      octants).
 */
        if (!cut[Z]) {
            subDecompUsed[LLB] = 0;           
            subDecompUsed[LRB] = 0;           
            subDecompUsed[ULB] = 0;           
            subDecompUsed[URB] = 0;           
            zDoms[0] = decomp->numDoms[Z];
        } else {
            zDoms[0] = decomp->numDoms[Z] >> 1;
            zDoms[1] = decomp->numDoms[Z] - zDoms[0];
        }

/*
 *      Loop through all the possible octants of this partition.
 *      For any octant that will be used, assign an RBDecomp_t
 *      structure to the octant, set the octant's subDecomp ID,
 *      and set the number of domains encompassed by each octant
 *      of the decomp.  all other initialization will be handled elsewhere.
 */
        for (i = 0; i < 8; i++) {

            if (!subDecompUsed[i]) {
                continue;
            }

            decomp->subDecomp[i] = &decompArray[*decompsUsed];
            *decompsUsed += 1;

            subDecomp = decomp->subDecomp[i];
/*
 *          For each octant used, we need to set the number of
 *          domains encompassed by the octant. 
 */
            switch (i) {
                case LLF:
                    subDecomp->numDoms[X] = xDoms[0];
                    subDecomp->numDoms[Y] = yDoms[0];
                    subDecomp->numDoms[Z] = zDoms[0];
                    break;
                case LRF:
                    subDecomp->numDoms[X] = xDoms[1];
                    subDecomp->numDoms[Y] = yDoms[0];
                    subDecomp->numDoms[Z] = zDoms[0];
                    break;
                case ULF:
                    subDecomp->numDoms[X] = xDoms[0];
                    subDecomp->numDoms[Y] = yDoms[1];
                    subDecomp->numDoms[Z] = zDoms[0];
                    break;
                case URF:
                    subDecomp->numDoms[X] = xDoms[1];
                    subDecomp->numDoms[Y] = yDoms[1];
                    subDecomp->numDoms[Z] = zDoms[0];
                    break;
                case LLB:
                    subDecomp->numDoms[X] = xDoms[0];
                    subDecomp->numDoms[Y] = yDoms[0];
                    subDecomp->numDoms[Z] = zDoms[1];
                    break;
                case LRB:
                    subDecomp->numDoms[X] = xDoms[1];
                    subDecomp->numDoms[Y] = yDoms[0];
                    subDecomp->numDoms[Z] = zDoms[1];
                    break;
                case ULB:
                    subDecomp->numDoms[X] = xDoms[0];
                    subDecomp->numDoms[Y] = yDoms[1];
                    subDecomp->numDoms[Z] = zDoms[1];
                    break;
                case URB:
                    subDecomp->numDoms[X] = xDoms[1];
                    subDecomp->numDoms[Y] = yDoms[1];
                    subDecomp->numDoms[Z] = zDoms[1];
                    break;
            }

/*
 *          Subpartition inherits a portion of its decomp ID from
 *          the parent partition; the last applicable byte identifies
 *          which octant the subdecomp is.
 */
            memcpy(decomp->subDecomp[i]->decompID,
                   decomp->decompID, level+1);
            decomp->subDecomp[i]->decompID[level+1] = i;

            InitRBDecomp(home, decomp->subDecomp[i], decompArray,
                          decompsUsed, level+1);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       AllocRBDecomp
 *      Description:    Allocates the tree needed for a full domain
 *                      decomposition, and either duplicates the structures
 *                      of a previously existing decomposition, or 
 *                      initializes the new decomposition from scratch.
 *
 *      Arguments:
 *          oldDecomp  Pointer to an existing decoposition.  Only
 *                     required if <allocType> is ALLOC_DUPLICATE_DECOMP.
 *          newDecomp  Location in which to return to the caller
 *                     a pointer to the highest level structure in the
 *                     decomposition tree.  The caller is responsible
 *                     for freeing this memory.
 *          allocType  Defines the desired behaviour.  Valid values are
 *                     ALLOC_NEW_DECOMP or ALLOC_DUPLICATE_DECOMP
 *
 *-------------------------------------------------------------------------*/
void AllocRBDecomp(Home_t *home, RBDecomp_t *oldDecomp, RBDecomp_t **newDecomp,
                   int allocType)
{
        int        i, j, level, decompsUsed, allocSize, index;
        RBDecomp_t *decompArray;
        static int totDecompCnt = 0;

/*
 *      Given the domain geometry, we know in advance how many
 *      RBDecomp_t structures will be required to hold the entire
 *      domain decomposition tree.  So preallocate the structures
 *      as a single block.
 */
        if (totDecompCnt == 0) {
            GetDecompCnt(home->param->nXdoms, home->param->nYdoms,
                         home->param->nZdoms, &totDecompCnt);
        }

        allocSize = totDecompCnt * sizeof(RBDecomp_t);
        decompArray = (RBDecomp_t *)calloc(1, allocSize);

/*
 *      Leave this check in.  At least 1 version of the gnu compiler
 *      has a bug that results in an allocSize of zero without something
 *      like this in place.
 */
        if (decompArray == (RBDecomp_t *)NULL) {
            printf("allocSize = %d\n", (int)allocSize);
            Fatal("AllocRBDecomp(): allocate error decompArray, "
                  "totDecompCnt %d, size %d", totDecompCnt, allocSize);
        }

/*
 *      The highest level of the decomposition will use the first
 *      element of the array of RBDecomp_t structures... which means
 *      we can free the entire decomposition tree by using the
 *      main decomposition pointer in a free() call.
 */
        *newDecomp = decompArray;
        decompsUsed = 1;
        level = 0;

/*
 *      If we just want to duplicate the structure of a previously
 *      initialized decomposition, do a wholesale copy of the 
 *      old to the new.
 */
        if (allocType == ALLOC_DUPLICATE_DECOMP) {

            if (oldDecomp == (RBDecomp_t *)NULL) {
                Fatal("AllocRBDecomp: Unable to copy NULL decomposition");
            }

            memcpy(decompArray, oldDecomp, allocSize);
/*
 *          Copying the subDecomp array within each structure copied
 *          pointers into the old decomposition.  Loop through the
 *          structures, do some pointer arithmetic to figure out
 *          indices of the structures in the old decomposition that
 *          are being referenced, then update the pointers to reference
 *          the corresponding structures in the duplicate copy.
 */
            for (i = 0; i < totDecompCnt; i++) {
                for (j = 0; j < 8; j++) {
                    if (decompArray[i].subDecomp[j] != (RBDecomp_t *)NULL) {
                        index = decompArray[i].subDecomp[j] - oldDecomp;
                        decompArray[i].subDecomp[j] = &decompArray[index];
                    }
                }
            }
        } else {
/*
 *          Recursively initialize the new decomposition tree, setting
 *          domain IDs, decompIDs, and pointers to sub-decompositions as
 *          necessary.
 */
            InitRBDecomp(home, *newDecomp, decompArray, &decompsUsed, level);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       RBDecomp
 *      Description:    Calculate the boundaries for the branch of
 *                      the decomposition tree terminating at the
 *                      local domain's decomposition based on the
 *                      current load (computation cost) distribution among
 *                      the processors.
 *
 *                      Once each domain has calculated its own boundaries
 *                      A call to ExchangeRBDecomp() can be used to
 *                      share those boundaries with all other processors
 *                      so all processors have the full decomposition.
 *
 *      Arguments:
 *          decomp       pointer to new decomposition tree which
 *                       has been pre-allocated and partially initialized
 *                       via AllocRBDecomp() and InitRBDecomp().
 *          oldDecomp    pointer to previous decomposition tree.
 *          domDecompID  ID of portion of decomposition associated with
 *                       the current domain
 *          currLevel    current level in the domain decomposition.
 *          startLevel   Level at which to start the decomposition.  The
 *                       caller may wish to do only a partial decomposition
 *                       in which case some number of levels of the
 *                       decompositiontree will simply have been copied
 *                       from a previous decomposition.  All levels less
 *                       than <startLevel> will simply be left as they are
 *
 *-------------------------------------------------------------------------*/
void RBDecomp(Home_t *home, RBDecomp_t *decomp, RBDecomp_t *oldDecomp,
              char *domDecompID, int currLevel, int startLevel)
{
        int          i, octant, thisSubpart=0;
        int          skipLevel = 0, restrictMovement = 0, doAllBoundaries = 0;
        int          numSlices;
        int          cut[3];
        real8        coord, oldBisection;
        real8        octantLoad[8], sliceLoad[2], oldBounds[3], minRatio[3];
        RBDecomp_t   *subDecomp, **subDecompList;


        if (decomp == (RBDecomp_t *)NULL) {
            return;
        }

        if (domDecompID == (char *)NULL) {
            doAllBoundaries = 1;
        } else {
            thisSubpart = domDecompID[currLevel+1];
        }

/*
 *      If we have the boundaries from the previous decomposition,
 *      We *may* want to restrict the motion of domain boundaries at
 *      the coarsest level at which we're doing a decomposition this
 *      cycle.
 *
 *      On the other hand, if we haven't yet reached the level at which
 *      to start the decomposition, set a flag so we skip this level --
 *      Domain boundaries at skipped levels have been copied from the
 *      the previous decomposition.
 */
        if (currLevel == startLevel) {
            if (oldDecomp != (RBDecomp_t *)NULL) restrictMovement = 1;
        } else if (currLevel < startLevel) {
            skipLevel = 1;
        }

/*
 *      Figure out the dimensions in which this partition is to be
 *      sub-partitioned.  Sub-partitioning will be into octants, quadrants
 *      halves or no subpartitions at all.
 */
        cut[X] = decomp->numDoms[X] > 1;
        cut[Y] = decomp->numDoms[Y] > 1;
        cut[Z] = decomp->numDoms[Z] > 1;

/*
 *      If there is no need to do further decomposition, we're
 *      down to an individual domain, so we can return to the caller.
 */
        if ((cut[X] + cut[Y] + cut[Z]) == 0) {
            return;
        }

/*
 *      When we cut the partion in a dimension, the two 'halves' may
 *      not be asigned the same number of procesors, so we need to 
 *      calculate the portion of the load that will be associated with
 *      the lower side of the bisection plane in each dimension.
 */
        minRatio[X] = (real8)(decomp->numDoms[X] / 2) /
                      (real8)decomp->numDoms[X];
        minRatio[Y] = (real8)(decomp->numDoms[Y] / 2) /
                      (real8)decomp->numDoms[Y];
        minRatio[Z] = (real8)(decomp->numDoms[Z] / 2) /
                      (real8)decomp->numDoms[Z];

        
/*
 *      The caller may have requested a partial decomposition rather
 *      than a full decomposition.  If that's the case, we may be able
 *      to simply skip this level (leaving the boundaries as they were
 *      in the previous decomposition) and move to the next level down.
 */
        if (!skipLevel) {
/*
 *          Set up the initial boundaries on all used subpartitions so
 *          they are identical to the current partition boundaries.
 *          As we determine the proper planes along which to subpartition
 *          the current partition, the subpartition boundaries will be 
 *          updated appropriately.
 */
            for (i = 0; i < 8; i++) {

                subDecomp = decomp->subDecomp[i];

                if (subDecomp != (RBDecomp_t *)NULL) {
    
                    subDecomp->cMin[X] = decomp->cMin[X];
                    subDecomp->cMin[Y] = decomp->cMin[Y];
                    subDecomp->cMin[Z] = decomp->cMin[Z];

                    subDecomp->cMax[X] = decomp->cMax[X];
                    subDecomp->cMax[Y] = decomp->cMax[Y];
                    subDecomp->cMax[Z] = decomp->cMax[Z];

                    octantLoad[i] = oldDecomp->subDecomp[i]->totLoad;
                } else {
                    octantLoad[i] = 0.0;
                }
            }

/*
 *          If the partition needs to be cut in the X dimension
 *          Find the coordinate of the plane along which the 
 *          partition needs to be cut in the X dimension and
 *          update the min/max coordinates for each subpartition
 *          that could be affected by slicing the partition in
 *          this dimension.
 */
            numSlices = 2;
            subDecompList = decomp->subDecomp;

            if (cut[X]) {
/*
 *              Given the load on each side of the current X bisection
 *              determine the coordinate for the new bisection.  At the
 *              coarsest level to be rebalanced this cycle, the distance
 *              a boundary may shift is restricted.  This helps prevent
 *              oscillation of boundaries.
 */
                sliceLoad[0] = octantLoad[LLF] + octantLoad[LLB] +
                               octantLoad[ULF] + octantLoad[ULB];
                sliceLoad[1] = octantLoad[LRF] + octantLoad[LRB] +
                               octantLoad[URF] + octantLoad[URB];
                oldBounds[0] = oldDecomp->cMin[X];
                oldBounds[1] = oldDecomp->subDecomp[LLF]->cMax[X];
                oldBounds[2] = oldDecomp->cMax[X];

                GetBisection(decomp->cMin[X], decomp->cMax[X], numSlices,
                             sliceLoad, oldBounds, &coord, minRatio[X]);

/*
 *              If we're restricting the domain boundary motion at the
 *              starting level, find where corresponding boundary was
 *              in the previous decomposition, and only move the new
 *              boundary a portion of the estimated distance it needs
 *              to move... helps limit boundary oscillation where
 *              density is low.
 */
                if (restrictMovement) {
                    oldBisection = oldBounds[1];
                    coord = oldBisection +
                            ((coord - oldBisection) * MAX_BOUND_SHIFT_FACT);
                }
/*
 *              If the current partition is not being subpartitioned
 *              in all 3 dimensions, not all of the octants will be
 *              have been allocated... We know certain octants must
 *              exist, but we have to verify the existence of the others
 *              before we try accessing them.  The same situation also
 *              applies further on when we cut the Y and Z dimension.
 */
                subDecompList[LLF]->cMax[X] = coord;
                subDecompList[LRF]->cMin[X] = coord;

                if (subDecompList[ULF]) subDecompList[ULF]->cMax[X] = coord;
                if (subDecompList[URF]) subDecompList[URF]->cMin[X] = coord;
                if (subDecompList[LLB]) subDecompList[LLB]->cMax[X] = coord;
                if (subDecompList[LRB]) subDecompList[LRB]->cMin[X] = coord;
                if (subDecompList[ULB]) subDecompList[ULB]->cMax[X] = coord;
                if (subDecompList[URB]) subDecompList[URB]->cMin[X] = coord;
            }

/*
 *          If the partition needs to be cut in the Y dimension
 *          determine the Y-coordinates of the planes along which
 *          to slice the (one or more) subpartitions formed after
 *          the partition has been sliced (possibly) in the X dimension,
 *          then update the min/max coordinates for each subpartition
 *          that could be affected by slicing the partition in this
 *          dimension.
 */
            if (cut[Y]) {

/*
 *              Unless we're calculating *all* the domain boundaries 
 *              we only need to get boundaries for the octant containing
 *              this domain... so only need to get the Y bisection of
 *              either the LLF or LRF octant but not both.
 */
                if (doAllBoundaries || ((thisSubpart & X_BITMASK) == 0)) {

                    sliceLoad[0] = octantLoad[LLF] + octantLoad[LLB];
                    sliceLoad[1] = octantLoad[ULF] + octantLoad[ULB];
                    oldBounds[0] = oldDecomp->cMin[Y];
                    oldBounds[1] = oldDecomp->subDecomp[LLF]->cMax[Y];
                    oldBounds[2] = oldDecomp->cMax[Y];

                    GetBisection(decomp->cMin[Y], decomp->cMax[Y], numSlices,
                                 sliceLoad, oldBounds, &coord, minRatio[Y]);

                    if (restrictMovement) {
                        oldBisection = oldBounds[1];
                        coord = oldBisection +
                                ((coord - oldBisection) * MAX_BOUND_SHIFT_FACT);
                    }

                    subDecompList[LLF]->cMax[Y] = coord;
                    subDecompList[ULF]->cMin[Y] = coord;

                    if (subDecompList[LLB]) subDecompList[LLB]->cMax[Y] = coord;
                    if (subDecompList[ULB]) subDecompList[ULB]->cMin[Y] = coord;
                }

/*
 *              Cut right side in Y
 */
                if ((subDecompList[LRF]) && (doAllBoundaries ||
                     ((thisSubpart & X_BITMASK) == X_BITMASK))) {

                    sliceLoad[0] = octantLoad[LRF] + octantLoad[LRB];
                    sliceLoad[1] = octantLoad[URF] + octantLoad[URB];
                    oldBounds[0] = oldDecomp->cMin[Y];
                    oldBounds[1] = oldDecomp->subDecomp[LRF]->cMax[Y];
                    oldBounds[2] = oldDecomp->cMax[Y];

                    GetBisection(decomp->cMin[Y], decomp->cMax[Y], numSlices,
                                 sliceLoad, oldBounds, &coord, minRatio[Y]);

                    if (restrictMovement) {
                        oldBisection = oldBounds[1];
                        coord = oldBisection +
                                ((coord - oldBisection) * MAX_BOUND_SHIFT_FACT);
                    }

                    subDecompList[LRF]->cMax[Y] = coord;
                    subDecompList[URF]->cMin[Y] = coord;

                    if (subDecompList[LRB]) subDecompList[LRB]->cMax[Y] = coord;
                    if (subDecompList[URB]) subDecompList[URB]->cMin[Y] = coord;
                }
            }  /* if (cut[Y]) */

/*
 *          If the partition needs to be cut in the Z dimension
 *          determine the z-coordinates of the planes along which
 *          to slice the (one or more) subpartitions formed after
 *          the partition has been sliced (possibly) in the X and
 *          Y dimensions, then update the min/max coordinates for
 *          each subpartition that could be affected by slicing
 *          the partition in this dimension.
 */
            if (cut[Z]) {

/*
 *              If we don't need *all* boundaries for the repartition utility
 *              we only need to get boundaries for the octant containing this
 *              domain... so only need to get the Z bisection of only
 *              one of the LLF, LRF, ULF, or URF octants
 */
                if (doAllBoundaries ||
                    (((thisSubpart & X_BITMASK) == 0) &&
                     ((thisSubpart & Y_BITMASK) == 0))) {

                    sliceLoad[0] = octantLoad[LLF];
                    sliceLoad[1] = octantLoad[LLB];
                    oldBounds[0] = oldDecomp->cMin[Z];
                    oldBounds[1] = oldDecomp->subDecomp[LLF]->cMax[Z];
                    oldBounds[2] = oldDecomp->cMax[Z];

                    GetBisection(decomp->cMin[Z], decomp->cMax[Z], numSlices,
                                 sliceLoad, oldBounds, &coord, minRatio[Z]);

                    if (restrictMovement) {
                        oldBisection = oldBounds[1];
                        coord = oldBisection +
                                ((coord - oldBisection) * MAX_BOUND_SHIFT_FACT);
                    }

                    subDecompList[LLF]->cMax[Z] = coord;
                    subDecompList[LLB]->cMin[Z] = coord;
                }

/*
 *              Cut lower right quadrant in Z
 */
                if ( subDecompList[LRF] && (doAllBoundaries || (((thisSubpart & X_BITMASK)==X_BITMASK) && ((thisSubpart & Y_BITMASK)==0))) ) 
                {
                    sliceLoad[0] = octantLoad[LRF];
                    sliceLoad[1] = octantLoad[LRB];
                    oldBounds[0] = oldDecomp->cMin[Z];
                    oldBounds[1] = oldDecomp->subDecomp[LRF]->cMax[Z];
                    oldBounds[2] = oldDecomp->cMax[Z];

                    GetBisection(decomp->cMin[Z], decomp->cMax[Z], numSlices,
                                 sliceLoad, oldBounds, &coord, minRatio[Z]);

                    if (restrictMovement) {
                        oldBisection = oldBounds[1];
                        coord = oldBisection +
                                ((coord - oldBisection) * MAX_BOUND_SHIFT_FACT);
                    }

                    subDecompList[LRF]->cMax[Z] = coord;
                    subDecompList[LRB]->cMin[Z] = coord;
                }

/*
 *              Cut upper left quadrant in Z
 */
                if ( subDecompList[ULF] && (doAllBoundaries || (((thisSubpart & X_BITMASK)==0) && ((thisSubpart & Y_BITMASK)==Y_BITMASK))) ) 
                {
                    sliceLoad[0] = octantLoad[ULF];
                    sliceLoad[1] = octantLoad[ULB];
                    oldBounds[0] = oldDecomp->cMin[Z];
                    oldBounds[1] = oldDecomp->subDecomp[ULF]->cMax[Z];
                    oldBounds[2] = oldDecomp->cMax[Z];

                    GetBisection(decomp->cMin[Z], decomp->cMax[Z], numSlices,
                                 sliceLoad, oldBounds, &coord, minRatio[Z]);

                    if (restrictMovement) {
                        oldBisection = oldBounds[1];
                        coord = oldBisection +
                                ((coord - oldBisection) * MAX_BOUND_SHIFT_FACT);
                    }

                    subDecompList[ULF]->cMax[Z] = coord;
                    subDecompList[ULB]->cMin[Z] = coord;
                }

/*
 *              Cut upper right quadrant in Z
 */
                if ( subDecompList[URF] && (doAllBoundaries || (((thisSubpart & X_BITMASK)==X_BITMASK) && ((thisSubpart & Y_BITMASK)==Y_BITMASK))) ) 
                {
                    sliceLoad[0] = octantLoad[URF];
                    sliceLoad[1] = octantLoad[URB];
                    oldBounds[0] = oldDecomp->cMin[Z];
                    oldBounds[1] = oldDecomp->subDecomp[URF]->cMax[Z];
                    oldBounds[2] = oldDecomp->cMax[Z];

                    GetBisection(decomp->cMin[Z], decomp->cMax[Z], numSlices,
                                 sliceLoad, oldBounds, &coord, minRatio[Z]);

                    if (restrictMovement) {
                        oldBisection = oldBounds[1];
                        coord = oldBisection +
                                ((coord - oldBisection) * MAX_BOUND_SHIFT_FACT);
                    }

                    subDecompList[URF]->cMax[Z] = coord;
                    subDecompList[URB]->cMin[Z] = coord;
                }
    
            }  /* if (cut[Z]) */

        }  /* if (!skipLevel) */

/*
 *      Boundaries are now properly set for all necessary octants of
 *      the current portion of the decomposition.  Now either continue
 *      with a full decomposition if we're generating a full decomposition
 *      for the repartitioning utility, or decompose just the octant of
 *      the tree with which current domain is associated.
 */
        if (doAllBoundaries) {
            for (octant = 0; octant < 8; octant++) {
                RBDecomp(home, decomp->subDecomp[octant],
                         oldDecomp->subDecomp[octant],
                         domDecompID, currLevel+1, startLevel);
            }
        } else {
            RBDecomp(home, decomp->subDecomp[thisSubpart],
                     oldDecomp->subDecomp[thisSubpart],
                     domDecompID, currLevel+1, startLevel);
        }


        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       UniformRBDecomp
 *      Description:    Generates a new uniform domain decomposition,
 *                      initializing domain boundaries at all levels
 *                      of a previously allocated decomposition tree.
 *                      This is a recursive process which continues
 *                      down to the single-domain partitions of the
 *                      decomposition tree.
 *      Arguments:
 *          decomp   pointer to new decomposition tree which
 *                   has been pre-allocated and partially initialized
 *                   via AllocRBDecomp() and InitRBDecomp().
 *          level    current level in the domain decomposition.
 *
 *-------------------------------------------------------------------------*/
void UniformRBDecomp(Param_t *param, RBDecomp_t *decomp, int level)
{
        int        i, octant;
        int        cut[3], subDecompUsed[8];
        real8      ratio;
        real8      Ctr[3];
        RBDecomp_t *subP;

/*
 *      Figure out the dimensions in which this subpartition must
 *      be further partitioned.  If no further partitioning is
 *      necessary, we can return to the caller.
 */
        cut[X] = (decomp->numDoms[X] > 1);
        cut[Y] = (decomp->numDoms[Y] > 1);
        cut[Z] = (decomp->numDoms[Z] > 1);

        if ((cut[X] + cut[Y] + cut[Z]) == 0) {
            return;
        }

/*
 *      Calculate the position in the volume composing this portion
 *      of the domain decomposition such that each octant contains
 *      roughly equivalent volume per domain.  The coordinates of
 *      the center will be used to define the planes for bisection
 *      in each required dimension.
 */
        for (i = 0; i < 3; i++) {
            if (cut[i]) {
                ratio = (real8)(decomp->numDoms[i]/2) /
                        (real8)decomp->numDoms[i];
            } else {
                ratio = 0.5;
            }
            Ctr[i] = decomp->cMin[i] + ratio *
                     (decomp->cMax[i] - decomp->cMin[i]);
        }

/*
 *      To start off assume all subpartitions will be
 *      further partitioned.  If that is not the case, 
 *      these values will be adjusted appropriately later.
 */
        for (octant = 0; octant < 8; octant++) {
            subDecompUsed[octant] = 1;
        }

/*
 *      If we're not cutting in the X dimension, none of the
 *      right side octants will be used
 */
        if (!cut[X]) {
            subDecompUsed[LRF] = 0;
            subDecompUsed[URF] = 0;
            subDecompUsed[LRB] = 0;
            subDecompUsed[URB] = 0;
        }

/*
 *      If we're not cutting in the Y dimension, none of the
 *      upper octants will be used
 */
        if (!cut[Y]) {
            subDecompUsed[ULF] = 0;
            subDecompUsed[URF] = 0;
            subDecompUsed[ULB] = 0;
            subDecompUsed[URB] = 0;
        }

/*
 *      If we're not cutting in the Z dimension, none of the
 *      back octants will be used
 */
        if (!cut[Z]) {
            subDecompUsed[LLB] = 0;
            subDecompUsed[LRB] = 0;
            subDecompUsed[ULB] = 0;
            subDecompUsed[URB] = 0;
        }


        for (octant = 0; octant < 8; octant++) {

            if (subDecompUsed[octant] == 0) {
                continue;
            }

/*
 *          Start each subpartition off with the same boundaries
 *          as the current partition then reset any boundaries
 *          that are affected by bisection of the current
 *          partition in one of more of the dimensions.
 */
            subP = decomp->subDecomp[octant];

            subP->cMin[X] = decomp->cMin[X];
            subP->cMax[X] = decomp->cMax[X];
            subP->cMin[Y] = decomp->cMin[Y];
            subP->cMax[Y] = decomp->cMax[Y];
            subP->cMin[Z] = decomp->cMin[Z];
            subP->cMax[Z] = decomp->cMax[Z];

            switch(octant) {
                case LLF:
                    if (cut[X]) subP->cMax[X] = Ctr[X];
                    if (cut[Y]) subP->cMax[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMax[Z] = Ctr[Z];
                    break;
                case LRF:
                    if (cut[X]) subP->cMin[X] = Ctr[X];
                    if (cut[Y]) subP->cMax[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMax[Z] = Ctr[Z];
                    break;
                case ULF:
                    if (cut[X]) subP->cMax[X] = Ctr[X];
                    if (cut[Y]) subP->cMin[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMax[Z] = Ctr[Z];
                    break;
                case URF:
                    if (cut[X]) subP->cMin[X] = Ctr[X];
                    if (cut[Y]) subP->cMin[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMax[Z] = Ctr[Z];
                    break;
                case LLB:
                    if (cut[X]) subP->cMax[X] = Ctr[X];
                    if (cut[Y]) subP->cMax[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMin[Z] = Ctr[Z];
                    break;
                case LRB:
                    if (cut[X]) subP->cMin[X] = Ctr[X];
                    if (cut[Y]) subP->cMax[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMin[Z] = Ctr[Z];
                    break;
                case ULB:
                    if (cut[X]) subP->cMax[X] = Ctr[X];
                    if (cut[Y]) subP->cMin[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMin[Z] = Ctr[Z];
                    break;
                case URB:
                    if (cut[X]) subP->cMin[X] = Ctr[X];
                    if (cut[Y]) subP->cMin[Y] = Ctr[Y];
                    if (cut[Z]) subP->cMin[Z] = Ctr[Z];
                    break;

            }

            UniformRBDecomp(param, subP, level+1);
        }

        return;
}
