/****************************************************************************
 *
 *      Module:      Decomp.c
 *      Description: Contains a set of generic functions used to create,
 *                   access, and manipulate the domain decomposition
 *                   regardless of the type of decomposition in use.
 *                   These generic functions will do any necessary
 *                   setup and invoke the appropriate functions to do
 *                   the real work based on the type of domain decomposition
 *                   currently in use. 
 *
 *      Includes public functions
 *          BroadcastDecomp()
 *          FindCoordDomain()
 *          FreeDecomp()
 *          GetAllDecompBounds()
 *          GetCellDomainList()
 *          GetLocalDomainBounds()
 *          ReadDecompBounds()
 *          Rebalance()
 *          UniformDecomp()
 *          WriteDecompBounds()
 *          XPlotDecomp()
 *
 *      Includes private functions
 *          DLBStats()
 *          GetLoadData()
 *
 ***************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Restart.h"
#include "Decomp.h"
#include "RSDecomp.h"
#include "RBDecomp.h"

/*
 *      Define some flags that toggle compile-time generation of
 *      some code for debugging the load-balancing algorithms.
 *      For each definition, 0 == OFF, 1 == ON.
 *
 *      DLB_DEBUG           Highest level flag enabling general debug code.
 *      DLB_PRINT_STATS     Enables message printed to stdout each DLB cycle
 *                          with load-balance stats such as min and max
 *                          per-process load, average load, and so on.  This
 *                          has no effect when DLB_DEBUG == 0.
 *      DLB_DUMP_DATAFILES  Enables creation of a file
 *                          (DLB_stats.<sequnce_num>) every DLB cycle with
 *                          the per-process load data.  This has no effect
 *                          when DLB_DEBUG == 0.
 */
#define DLB_DEBUG          0
#define DLB_PRINT_STATS    1
#define DLB_DUMP_DATAFILES 0


#if DLB_DEBUG
/*-----------------------------------------------------------------------
 *
 *      Function:     DLBStats
 *      Description:  This function is only provided for debugging
 *                    purposes.  It handles writing (to stdout and/or
 *                    disk) some data related to the current load
 *                    on all processors.  Only domain zero will
 *                    do any output here...
 *
 *      Arguments:
 *          loadData  Array containing current load data (wallclock
 *                    time or force calc counts) for each domain.
 *                    Array is assumed to have same number of elements
 *                    as the job's domain count.
 *
 *----------------------------------------------------------------------*/
static void DLBStats(Home_t *home, real8 *loadData)
{
        int        domain;
        real8      minVal, maxVal, totVal, avgVal, imbalance;
        static int seqNum = 1;

#if DLB_DUMP_DATAFILES
        char       statFileName[64];
        FILE       *statFP;

        if (home->myDomain != 0) {
            return;
        }

/*
 *      If enabled, dump the current load data for all processors
 *      to a disk file.
 */
        sprintf(statFileName, "DLB_stats.%d", seqNum);
        statFP = fopen(statFileName, "w");
        if (statFP == (FILE *)NULL) {
            printf("  *** Error %d opening file %s\n", errno, statFileName);
            return;
        }

        for (domain = 0; domain < home->numDomains; domain++) {
            fprintf(statFP, "%d  %lf\n", domain, loadData[domain]);
        }

        fclose(statFP);
#endif  // DLB_DUMP_DATAFILES

#if DLB_PRINT_STATS
        if (home->myDomain != 0) {
            return;
        }

/*
 *      If enabled, print to stdout a single line with some stats
 *      about the current load balance.
 */
        minVal = loadData[0];
        maxVal = 0.0;
        totVal = 0.0;
        avgVal = 0.0;

        for (domain = 0; domain < home->numDomains; domain++) {
            if (loadData[domain] < minVal) minVal = loadData[domain];
            if (loadData[domain] > maxVal) maxVal = loadData[domain];
            totVal += loadData[domain];
        }

        avgVal = totVal / home->numDomains;

        if (avgVal == 0.0) {
            imbalance = 0.0;
        } else {
            imbalance = (maxVal - avgVal) / avgVal;
        }

        printf(" +++ DLB %d: min=%lf max=%lf avg=%lf imbalance=%lf\n",
               seqNum, minVal, maxVal, avgVal, imbalance);
#endif  // DLB_PRINT_STATS

        seqNum++;

        return;
}
#endif  // DLB_DEBUG

/*------------------------------------------------------------------------
 *
 *      Function:       BroadcastDecomp
 *      Description:    Generic function to control broadcast of the
 *                      domain decomposition from task zero to all other
 *                      tasks.  This function will determine the type of
 *                      domain decomposition is use and invoke the
 *                      broadcast function appropriate to the decomposition.
 *
 *      Arguments:
 *          decomp   Pointer to decomposition to be broadcast from domain
 *                   zero to all other domains. Pointer is NULL or
 *                   uninitialized on all tasks but zero.
 *
 *-----------------------------------------------------------------------*/
void BroadcastDecomp(Home_t *home, void *decomp)
{
        switch (home->param->decompType) {
        case 1:
            BroadcastRSDecomp(home, (RSDecomp_t *)decomp);
            break;
        case 2:
            BroadcastRBDecomp(home, (RBDecomp_t *)decomp);
            break;
        }

/*
 *      Have each process copy its own local boundaries into the
 *      <home> structure... saves multiple lookups later on.
 */
        GetLocalDomainBounds(home, home->decomp);

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       UniformRBDecomp
 *      Description:    Generic function to generate a uniform domain
 *                      decomposition in which all domains will encompass
 *                      volumes of equivalent size and shape.  This function
 *                      will determine the type of domain decomposition in
 *                      use and invoke the appropriate decomposition function.
 *
 *      Arguments:
 *          decomp  Location in which to return to the caller a pointer
 *                  to a new uniform domain decomposition
 *
 *-------------------------------------------------------------------------*/
void UniformDecomp(Home_t *home, void **decomp)
{
        int        level = 0;
        Param_t    *param;

        param = home->param;

        switch (home->param->decompType) {
        case 1:
            UniformRSDecomp(param, (RSDecomp_t **)decomp);
            break;
        case 2:
/*
 *          Allocate and initialize the uniform decomposition
 *          but DO NOT do an exchange of boundaries here.  Task
 *          0 (which is doing the uniform decomposition) will
 *          broadcast the decomposition later.
 */
            AllocRBDecomp(home, (RBDecomp_t *)NULL, (RBDecomp_t **)decomp,
                          ALLOC_NEW_DECOMP);
            UniformRBDecomp(param, (RBDecomp_t *)*decomp, level);
            break;
        }
 
        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:       ReadDecompBounds
 *      Description:    Generic function to read the domain decomposition
 *                      from a previous restart file.  This function
 *                      will determine the type of domain decomposition
 *                      the file contains and invoke a the read function
 *                      appropriate to the decomposition type.
 *
 *                      NOTE: If the restart file does not contain a 
 *                      domain decomposition or the decomposition type
 *                      in the file does not match the decomposition
 *                      type for this run of the application, no
 *                      decomposition will be returned to the caller.
 *                      A uniform decomposition will have to be generated
 *                      elsewhere.
 *
 *      Arguments:
 *          filePtr        Handle to open file from which the decomposition
 *                         will be read.
 *          dobinRead      flag indicating that the data is being read
 *                         from an HDF5 restart file.
 *          newDecompType  Indicates the type of decomposition which will
 *                         be used for this run of the application.
 *          oldDecomp      location in which to return to the caller the
 *                         decomposition (if any) read from the file.
 *
 *-------------------------------------------------------------------------*/
void ReadDecompBounds(Home_t *home, void **filePtr, int doBinRead,
                      int newDecompType, void **oldDecomp)
{
        int     saveDecomp, oldDecompType;
        int     decompDomX, decompDomY, decompDomZ;
        Param_t *param;

        param = home->param;

        saveDecomp = 1;

        oldDecompType = param->dataDecompType;

        decompDomX = param->dataDecompGeometry[X];
        decompDomY = param->dataDecompGeometry[Y];
        decompDomZ = param->dataDecompGeometry[Z];

/*
 *      If the domain geometry was not provided, no domain decomposition
 *      will have been provided either, so just return to the caller.
 */
        if ((decompDomX * decompDomY * decompDomZ) == 0) {
            return;
        }

/*
 *      If the domain geometry in the restart file does not match
 *      the current domain geometry, or the type of domain
 *      decomposition contained in the restart file does not
 *      match the type of decomposition to be used for this run,
 *      there's no point in storing the decomposition for the
 *      caller -- we'll still need to read through the decomposition
 *      in order to get to the nodal data, though.
 */
        if ((decompDomX != home->param->nXdoms) ||
            (decompDomY != home->param->nYdoms) ||
            (decompDomZ != home->param->nZdoms) ||
            (oldDecompType != newDecompType)) {
            saveDecomp = 0;
        }

        if (oldDecomp == (void **)NULL) {
            saveDecomp = 0;
        }

/*
 *      If we're reading from a binary restart file, and we
 *      don't need to save the decomposition, we don't even need
 *      read the decomposition.
 */
        if ((doBinRead && saveDecomp == 0)) {
            return;
        }

/*
 *      Read the appropriate type of domain decomposition from 
 *      the restart file.
 */
        switch (oldDecompType) {
        case 1:
            if (doBinRead) {
                ReadBinRSDecompBounds((void *)filePtr, decompDomX, decompDomY,
                                      decompDomZ, (RSDecomp_t **)oldDecomp);
            } else {
                ReadRSDecompBounds(filePtr, decompDomX, decompDomY, decompDomZ,
                                   saveDecomp, (RSDecomp_t **)oldDecomp);
            }
            break;
        case 2:
            if (doBinRead) {
                ReadBinRBDecompBounds(home, (void *)filePtr, decompDomX, decompDomY,
                                      decompDomZ, (RBDecomp_t **)oldDecomp);
            } else {
                ReadRBDecompBounds(home, filePtr, decompDomX, decompDomY,
                               decompDomZ, saveDecomp,
                               (RBDecomp_t **)oldDecomp);
            }
            break;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    WriteDecompBounds
 *      Description: Generic function to write the domain decomposition
 *                   into a restart file.  This function will determine
 *                   the type of domain decomposition in use and invoke
 *                   the write function appropriate to the decomposition
 *                   type.
 *
 *      Arguments:
 *          fp          open file descriptor for the restart file being
 *                      written
 *          versionNum  version number of the nodal data file being written.
 *
 *------------------------------------------------------------------------*/
void WriteDecompBounds(Home_t *home, FILE *fp)
{
        int level = 0;

        switch(home->param->decompType) {
        case 1:
            WriteRSDecompBounds(home, fp, (RSDecomp_t *)home->decomp);
            break;
        case 2:
            WriteRBDecompBounds(home, fp, (RBDecomp_t *)home->decomp,
                                level);
            break;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetAllDecompBounds
 *      Description: Generic function to create a single dimension array
 *                   containing all the domain boundaries.  This function
 *                   will determine the type of domain decomposition in
 *                   use, calculate the necessary array size and invoke
 *                   the 'get' function appropriate to the decomposition
 *                   type.
 *
 *      Arguments:
 *          bounds      location in whcih to return to the caller the
 *                      pointer to the array of domain boundaries.
 *          numBounds   location in which to return to the caller the
 *                      number of elements in the <*bounds> array.
 *
 *------------------------------------------------------------------------*/
void GetAllDecompBounds(Home_t *home, real8 **bounds, int *numBounds)
{
        int xDoms, yDoms, zDoms;

        switch(home->param->decompType) {
        case 1:
            xDoms = home->param->nXdoms;
            yDoms = home->param->nYdoms;
            zDoms = home->param->nZdoms;
            *numBounds = xDoms * (yDoms * (zDoms+1) + (yDoms+1)) + (xDoms+1);
            *bounds = (real8 *)malloc(*numBounds * sizeof(real8));
            GetAllRSDecompBounds(home, (RSDecomp_t *)home->decomp, *bounds);
            break;
        case 2:
            *numBounds = 6 * home->numDomains;
            *bounds = (real8 *)malloc(*numBounds * sizeof(real8));
            GetAllRBDecompBounds(home, (RBDecomp_t *)home->decomp, *bounds);
            break;
        }

        return;
}


#ifndef NO_XWINDOW
/*-------------------------------------------------------------------------
 *
 *      Function:    XPlotDecomp
 *      Description: Generic function to plot the domain decomposition
 *                   in the X-Window display.  This function will 
 *                   invoke the plotting function appropriate to the
 *                   decomposition type.
 *
 *      Arguments:
 *          xMin    Minimum boundary of problem space in X dimension
 *          yMin    Minimum boundary of problem space in Y dimension
 *          zMin    Minimum boundary of problem space in Z dimension
 *          lMax    Length of problem space in largest dimension
 *          color   index of color to use for plotting domain boundaries
 *          lineWidth width of lines used to plot domain boundaries
 *
 *------------------------------------------------------------------------*/
void XPlotDecomp(Home_t *home, real8 xMin, real8 yMin, real8 zMin, real8 lMax,
                 int color, real8 lineWidth)
{

        switch(home->param->decompType) {
        case 1:
            XPlotRSDecomp(home, xMin, yMin, zMin, lMax, color, lineWidth);
            break;
        case 2:
            XPlotRBDecomp(home, (RBDecomp_t *)home->decomp, xMin, yMin, zMin,
                        lMax, color, lineWidth);
            break;
        }

        return;
}
#endif



/*-------------------------------------------------------------------------
 *
 *      Function:    FreeDecomp
 *      Description: Generic function to release memory associated with
 *                   a domain decomposition.  This function will 
 *                   invoke the memory freeing function appropriate
 *                   to the decomposition type.
 *
 *      Arguments:
 *          decomp   pointer to the domain decomposition to be freed
 *
 *------------------------------------------------------------------------*/
void FreeDecomp(Home_t *home, void *decomp)
{
        if (decomp != (void *)NULL) {

            switch(home->param->decompType) {
            case 1:
                FreeRSDecomp(home, (RSDecomp_t *)decomp);
                free(decomp);
                break;
            case 2:
                FreeRBDecomp((RBDecomp_t *)decomp);
                break;
            }
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetLocalDomainBounds
 *      Description: Generic function to search the domain decomposition
 *                   for the current domain's boundaries and store them
 *                   in the <home> structure (to save multiple lookups
 *                   later).
 *
 *      Arguments:
 *          decomp   pointer to the current domain decomposition
 *
 *------------------------------------------------------------------------*/
void GetLocalDomainBounds(Home_t *home, void *decomp)
{
        switch (home->param->decompType) {
        case 1:
            GetRSDecompLocalDomainBounds(home, (RSDecomp_t *)decomp);
            break;
        case 2:
            GetRBDecompLocalDomainBounds(home, (RBDecomp_t *)decomp);
            break;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetCellDomainList
 *      Description: Generic function to search the domain decomposition
 *                   for the list of domains intersecting the specified
 *                   cell.  This function will invoke the search function
 *                   appropriate to the type of domain decomposition
 *                   being used.
 *
 *      Arguments:
 *          cellID    ID of cell as returned by EncodeCellIdx().
 *          domCount  Location in which to return to caller the number of
 *                    domains intersecting the specified cell.
 *          domList   Location in which to return to caller the array
 *                    containing the IDs of all domains intersecting the
 *                    specified cell.
 *
 *------------------------------------------------------------------------*/
void GetCellDomainList(Home_t *home, int cellID, int *domCount, int **domList)
{

        switch (home->param->decompType) {
        case 1:
            GetRSDecompCellDomainList(home, cellID, domCount, domList);
            break;
        case 2:
            GetRBDecompCellDomainList(home, cellID, domCount, domList);
            break;
        }

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    GetLoadData
 *      Description: Obtains the per-process load information required
 *                   to do dynamic load-balancing across the processors.
 *
 *      Arguments:
 *          criteria  Identifies the type of load data to obtain.
 *                    DLB_USE_FORCECALC_COUNT = number of force calculations
 *                                              done this cycle
 *                    DLB_USE_WALLCLK_TIME    = wallclock time spent doing
 *                                              force calculations this
 *                                              cycle.
 *          loadData  Location in which to return to the caller a pointer
 *                    to the load data array.  The array consists of 1
 *                    value per process.
 *
 *------------------------------------------------------------------------*/
static void GetLoadData(Home_t *home, int criteria, real8 **loadData)
{
        real8  localLoad;
        real8  *globalLoadData;


        switch (criteria) {
            case DLB_USE_FORCECALC_COUNT:
                localLoad = (real8)home->cycleForceCalcCount;
                break;
            case DLB_USE_WALLCLK_TIME:
            default:
                localLoad = home->timers[CALC_FORCE].save;
                break;
        }

/*
 *      Caller is responsible for freeing this buffer!
 */
        globalLoadData = (real8 *)malloc(home->numDomains * sizeof(real8));

#ifdef PARALLEL
        MPI_Allgather(&localLoad, 1, MPI_DOUBLE, globalLoadData, 1,
                      MPI_DOUBLE, MPI_COMM_WORLD);
#else
        *globalLoadData = localLoad;
#endif

        *loadData = globalLoadData;

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    Rebalance
 *      Description: Generic function to dynamically rebalance the
 *                   work load by adjusting the domain decomposition.
 *                   This function will obtain updated per-process load
 *                   data which it will then pass on to the decomposition
 *                   (and related rebalancing) functions appropriate to
 *                   the type of domain decomposition being used.
 *
 *      Arguments:
 *          criteria  Identifies the type of per process load data to obtain.
 *                    DLB_USE_FORCECALC_COUNT = number of force calculations
 *                                              done this cycle
 *                    DLB_USE_WALLCLK_TIME    = wallclock time spent doing
 *                                              force calculations this
 *                                              cycle.
 *------------------------------------------------------------------------*/
void Rebalance(Home_t *home, int criteria)
{
        int         didDLB = 0, currLevel = 0;
        real8       imbalanceThreshold;
        real8       *loadData;
        Param_t     *param;
        RBDecomp_t  *newDecomp;
        static int  firstTime = 1, maxLevel = 0;
        static int  copyLevel = 0;
        static char domDecompID[MAX_DECOMP_LVLS];


        param = home->param;

/*
 *      Save the current time spent on force calculations during this
 *      cycle as a basis for determining load balance then start the
 *      timer for the load balancing.
 */
        TimerSave(home, CALC_FORCE);
        TimerSave(home, CELL_CHARGE);
        TimerStart(home, LOAD_BALANCE);

/*
 *      There's a number of conditions that determine whether to do
 *      load balancing this cycle... check 'em all.
 */
        if ((param->numDLBCycles > 0) || ((param->DLBfreq > 0) &&
            (((param->decompType == 1)&&(home->cycle%param->DLBfreq <= 3)) ||
             ((param->decompType == 2)&&(home->cycle%param->DLBfreq == 0))))) {

/*
 *          Decomposition requires the raw per-process load
 *          data, so get the data (and dump statistics about
 *          the current load imbalance if necessary).
 */
            GetLoadData(home, criteria, &loadData);
#if DLB_DEBUG
            DLBStats(home, loadData);
#endif

            switch (param->decompType) {
            case 1:
/*
 *              For recursive sectioning domain decompositions, load
 *              balancing occurs over a 3-cycle period; balance along
 *              X dimension the first cycle, along the Y dimension the
 *              second cycle and along Z the third cycle.
 */
                switch (home->cycle % param->DLBfreq) {
                case 0:
                    DLBalanceX(home, loadData);
                    break;
                case 1:
                    DLBalanceY(home, loadData);
                    break;
                case 2:
                    DLBalanceZ(home, loadData);
                    break;
                }

                GetLocalDomainBounds(home, home->decomp);
                didDLB = 1;
                break;

            case 2:
/*
 *              First time into this section we need to calculate
 *              the current domain's decomp ID for use elsewhere...
 */
                if (firstTime) {
                    firstTime = 0;
                    DomID2DecompID(home, home->myDomain, domDecompID);
                    maxLevel = MAX(home->xMaxLevel, home->yMaxLevel);
                    maxLevel = MAX(home->zMaxLevel, maxLevel);
                }

/*
 *              Use an imbalance threshold to choose the level to
 *              rebalance.  Also set the load for each branch of
 *              the decomposition tree based on the raw load data.
 */
                copyLevel = maxLevel;
                imbalanceThreshold = 0.05;

                (void) RBCheckLoadBalance(home, (RBDecomp_t *)home->decomp,
                                          loadData, imbalanceThreshold,
                                          currLevel, &copyLevel);

/*
 *              If copyLevel is >= the maximum level, there were no
 *              levels of the decomposition tree that had an imbalanced
 *              greater than the threshold and no rebalancing is needed
 */
                if (copyLevel < maxLevel) {
/*
 *                  Allocate a new decomposition tree that initially
 *                  duplicates the previous one.  We'll adjust boundaries
 *                  for whichever levels of the tree we need to later on.
 */
                    AllocRBDecomp(home, (RBDecomp_t *)home->decomp, &newDecomp,
                                  ALLOC_DUPLICATE_DECOMP);

/*                  Have each task calculate the boundaries for only the
 *                  branch of the tree terminating at the local domain's leaf.
 *                  Once that is done, an exchange is done so all domains
 *                  learn all the domain boundaries and initialize the
 *                  rest of the decomposition tree.
 */
                    RBDecomp(home, newDecomp, (RBDecomp_t *)home->decomp,
                             domDecompID, currLevel, copyLevel);

                    GetLocalDomainBounds(home, newDecomp);

#ifdef PARALLEL
                    ExchangeRBDecomp(home, newDecomp);
#endif
                    FreeRBDecomp((RBDecomp_t *)home->decomp);
                    home->decomp = (void *)newDecomp;
                    didDLB = 1;
                }

                break;
            }

/*
 *          If a rebalance was actually performed, then free
 *          the old cell and neighbor domain structures (since
 *          cell membership, neighbors, and neighboring domains are
 *          dependent on the domain boundaries) and reinitialize
 *          the structures.
 */
            free(loadData);

            if (didDLB) {
                DLBfreeOld(home);

                InitCellNatives(home);
                InitCellNeighbors(home);
                InitCellDomains(home);
                InitRemoteDomains(home);
            }
        }


        TimerStop(home, LOAD_BALANCE);

        return;
}


/*-------------------------------------------------------------------------
 *
 *      Function:    FindCoordDomain
 *      Description: Generic function to search the domain decomposition
 *                   to locate the domain 'owning' the specified
 *                   coordinates.  This function will shift the specified
 *                   coordinates into the primary image space (if
 *                   necessary) then invoke the search function appropriate
 *                   to the type of domain decomposition being used.
 *
 *      Arguments:
 *          updateCoords Flag indicating if coordinates outside the primary
 *                       image should be converted to the corresponding
 *                       coordinates.  Any non-zero value permits the
 *                       coordinate update.  Note: regardless of the
 *                       value of this flag, no coordinate update will
 *                       be done if periodic boundaries are not enabled.
 *          x, y, z      Pointers to coordinates.  Contents may be updated
 *                       depending o9n the value of <updateCoords>.
 *
 *------------------------------------------------------------------------*/
int FindCoordDomain(Home_t *home, int updateCoords,
                    real8 *x, real8 *y, real8 *z)
{
        int     domID;
        real8   newX, newY, newZ;
        real8   xMin, yMin, zMin;
        real8   xMax, yMax, zMax;
        real8   xAdjustment, yAdjustment, zAdjustment;
        Param_t *param;

        param = home->param;
        domID = home->myDomain;

        newX = *x;
        newY = *y;
        newZ = *z;

/*
 *      Find the min/max coordinates of the current domain.  If the
 *      specified position is contained within the current domain
 *      we don't need to go any further.
 */
        xMin = home->domXmin;
        yMin = home->domYmin;
        zMin = home->domZmin;

        xMax = home->domXmax;
        yMax = home->domYmax;
        zMax = home->domZmax;

        if (((newX >= xMin) && (newX < xMax)) &&
            ((newY >= yMin) && (newY < yMax)) &&
            ((newZ >= zMin) && (newZ < zMax))) {
            return(domID);
        }

/*
 *      We may need to shift coordinates outside the primary image
 *      back into the primary.  If periodic boundaries are enabled
 *      we define a non-zero adjustment value for the coordinates;
 *      for free surfaces, we leave the coordinates as they are.
 */
        xAdjustment = param->maxSideX - param->minSideX;
        yAdjustment = param->maxSideY - param->minSideY;
        zAdjustment = param->maxSideZ - param->minSideZ;

/*
 *      Adjustments become zero if periodic boundaries are not enabled.
 */
        xAdjustment *= (real8)(param->xBoundType == Periodic);
        yAdjustment *= (real8)(param->yBoundType == Periodic);
        zAdjustment *= (real8)(param->zBoundType == Periodic);


        if (newX < param->minSideX) {
            newX += xAdjustment;
        } else if (newX >= param->maxSideX) {
            newX -= xAdjustment;
        }

        if (newY < param->minSideY) {
            newY += yAdjustment;
        } else if (newY >= param->maxSideY) {
            newY -= yAdjustment;
        }

        if (newZ < param->minSideZ) {
            newZ += zAdjustment;
        } else if (newZ >= param->maxSideZ) {
            newZ -= zAdjustment;
        }

/*
 *      If necessary, return new coordinates to the caller.
 */
        if (updateCoords) {
            *x = newX;
            *y = newY;
            *z = newZ;
        }

/*
 *      Call the function appropriate to the current domain decomposition
 *      to identify the domain encompassing the specified coordinate.
 */
        switch (param->decompType) {
        case 1:
            domID = FindRSDecompCoordDomain(home, (RSDecomp_t *)home->decomp,
                                           newX, newY, newZ);
            break;
        case 2:
/*
 *          If periodic boundaries are not enabled, it *may* be possible
 *          that the coordinates are still outside the boundaries of the
 *          problem space.  For this type of decomposition, determining
 *          ownership of such coordinates would be problematic, so we'll
 *          shift the coordinates to the closest boundary (temporarily)
 *          in order to determine the owning domain.
 */
            if      (newX < param->minSideX)  newX = param->minSideX;
            else if (newX > param->maxSideX)  newX = param->maxSideX;

            if      (newY < param->minSideY)  newY = param->minSideY;
            else if (newY > param->maxSideY)  newY = param->maxSideY;

            if      (newZ < param->minSideZ)  newZ = param->minSideZ;
            else if (newZ > param->maxSideZ)  newZ = param->maxSideZ;

            domID = FindRBDecompCoordDomain((RBDecomp_t *)home->decomp,
                                            newX, newY, newZ);

            break;
        }

        return(domID);
}
