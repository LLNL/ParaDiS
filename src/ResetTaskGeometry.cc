/****************************************************************************
 * 
 *      Module:       ResetTaskGeometry.c
 *
 *      Description:  Contains a generic dispatch function and the
 *                    platform-specific functions to query the system
 *                    for the shape/size of the hardware partition and
 *                    reset the logical task geometry (if appropriate)
 *                    to match the hardware partition.
 *
 *      Includes public functions:
 *          ResetTaskGeometry()
 *
 *      Includes private functions:
 *          GetOptimal3DTaskGeometry()
 *          ResetTaskGeometry_BGP()
 *          ResetTaskGeometry_BGQ()
 *
 ****************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

#ifdef PARALLEL

#ifdef _BGP
#include <spi/kernel_interface.h>
#include <common/bgp_personality.h>
#include <common/bgp_personality_inlines.h>
#define TASK_MAP_ENV_VAR         "BG_MAPPING"
#define TASK_DISTRIBUTION_DEFLT  "XYZT"
#endif

#ifdef _BGQ
#include <limits.h>
#include <spi/include/kernel/location.h>
#endif

#define TASK_MAP_MISMATCH_USE_PHYSICAL 0
#define TASK_MAP_MISMATCH_USE_LOGICAL  1
#define TASK_MAP_MISMATCH_ABORT        2

#ifdef _BGQ
/*
 *      Define the number of coordinates needed to identify the
 *      location of a single task.  For BG/Q systems this corresponds
 *      to the 5 dimensions of the hardware partition and we add an
 *      extra coordinate to represent the task index within each node.
 *
 */
#define BGQ_NUM_COORD 6


/*---------------------------------------------------------------------------
 *
 *      Function:     GetOptimal3DTaskGeometry
 *
 *      Description:  For BG/Q systems, we may not know the physical geometry
 *                    of the 5D partition until run-time, which means the
 *                    logical domain decomposition specified by the user
 *                    may not match the hardware partition.  This function
 *                    will query the system to determine the shape/size
 *                    of the hardware partition and attempt to calculate
 *                    an optimal task geometry based on the total number
 *                    of tasks in the simulation, the shape of the physical
 *                    hardware partition and the order in which the tasks
 *                    are distributed over the dimensions of the partition.
 *
 *                    NOTES:
 *
 *                    1) This process must be called by ALL tasks in
 *                       the simulation since they all participate in sending
 *                       task 0 information it needs to analyze.  Only task
 *                       0 will calculate the optimal geometry.
 *
 *                    2) The optimal geometry will be returned to the caller
 *                       in the optimalTaskGeometry array of the Param_t
 *                       structure.  However, if the optimal geometry can
 *                       not be determined, the values in this array will each
 *                       be set to -1.
 *
 *                    3) This function MUST be called before the call
 *                       to ResetTaskGeometry().
 *
 *-------------------------------------------------------------------------*/
void GetOptimal3DTaskGeometry(Home_t *home, int optimalTaskGeometry[3])
{
        int           numTasks;
        int           localCoordinates[6];
        int           *coordinateList = (int *)NULL;
        Personality_t personality;


/*
 *      Task 0 needs to allocate a buffer to hold the aggregated list
 *      of task coordinates.
 */
        if (home->myDomain == 0) {
            int bufSize;

            numTasks = home->numDomains;
            bufSize = BGQ_NUM_COORD * numTasks * sizeof(int); 

            coordinateList = (int *)malloc(bufSize);
        }

/*
 *      Have each task query the OS to determine it's coordinates in
 *      the hardware partition allocated to the job.
 *
 *          coordinate[0] : A dimension
 *          coordinate[1] : B dimension
 *          coordinate[2] : C dimension
 *          coordinate[3] : D dimension
 *          coordinate[4] : E dimension
 *          coordinate[5] : Task index within a node
 */
        Kernel_GetPersonality(&personality, sizeof(Personality_t));

        localCoordinates[0] = personality.Network_Config.Acoord;
        localCoordinates[1] = personality.Network_Config.Bcoord;
        localCoordinates[2] = personality.Network_Config.Ccoord;
        localCoordinates[3] = personality.Network_Config.Dcoord;
        localCoordinates[4] = personality.Network_Config.Ecoord;

        localCoordinates[5] = Kernel_ProcessorID(); /* 0 - 63 for BGQ */

/*
 *      Ship all the task coordinates to task 0 which will analyze the
 *      data.
 */
        MPI_Gather(localCoordinates, BGQ_NUM_COORD, MPI_INT, 
                   coordinateList, BGQ_NUM_COORD, MPI_INT, 0, MPI_COMM_WORLD);


/*
 *      Task 0 now has to analyze all the tasks' coordinates to determine
 *        1) The task distribution among dimensions (i.e. which dimensions
 *           vary fastest, slowest)
 *        2) The minimum/maximum per-dimension coordinates in order to 
 *           define the physical shape of the partition
 */
        if (home->myDomain == 0) {
            int i;
            int cIndex, taskID;
            int lenSumMin;
            int computeNodeCount, tasksPerNode;
            int allocOrderIndex;
            int allocOrder[BGQ_NUM_COORD];
            int allocOrderFound[BGQ_NUM_COORD];
            int dimSize[BGQ_NUM_COORD];
            int minCoord[BGQ_NUM_COORD], maxCoord[BGQ_NUM_COORD];
            int optimal3DGeometry[3];

            allocOrderIndex = BGQ_NUM_COORD - 1;

/*
 *          Initially set the per-dimenson min/max coordinates to those
 *          of task 0, and other initialization.
 */
            for (cIndex = 0; cIndex < BGQ_NUM_COORD; cIndex++) {

                minCoord[cIndex] = coordinateList[cIndex];
                maxCoord[cIndex] = coordinateList[cIndex];

                allocOrderIndex = BGQ_NUM_COORD - 1;
                allocOrderFound[cIndex] = 0;
                allocOrder[cIndex] = -1;
            }

/*
 *          Now loop through coordinates of all tasks other than task 0
 */
            for (taskID = 1; taskID < numTasks; taskID++) {

                for (cIndex = 0; cIndex < BGQ_NUM_COORD; cIndex++) {
                    int offset = BGQ_NUM_COORD * taskID;
                    int offsetPrev = offset - BGQ_NUM_COORD;

/*
 *                  If necessary, adjust the min/max coordinates of the
 *                  partition based on the task coordinates
 */
                    if (coordinateList[offset+cIndex] < minCoord[cIndex]) {
                        minCoord[cIndex] = coordinateList[offset+cIndex];
                    }

                    if (coordinateList[offset+cIndex] > maxCoord[cIndex]) {
                        maxCoord[cIndex] = coordinateList[offset+cIndex];
                    }

/*
 *                  If the task coordinate in this dimension increased
 *                  and this is the first time we encounted an increase
 *                  in this dimension, store that info away for later.
 *
 *                  NOTE: The allocOrder[*] array stores the dimension
 *                        indices in order of the slowest to fastest
 *                        changing dimensions
 */
                    if (coordinateList[cIndex+offset] >
                        coordinateList[cIndex+offsetPrev]) {

                        if (allocOrderFound[cIndex]) {
                            continue;
                        }

                        allocOrder[allocOrderIndex--] = cIndex;
                        allocOrderFound[cIndex] = 1;
                    }

                }

            }  /* end for (taskID = 0; ... ) */

/*
 *          It is possible that the length in some dimensions is only
 *          1 in which case the code above did not fully fill in the
 *          allocOrder array.  For the dimensions of length 1, the
 *          order doesn't matter, but we need to finish initializing
 *          the array with valid numbers or we'll run into problems
 *          later on.
 */
            while (allocOrderIndex >= 0) {
                for (cIndex = 0; cIndex < BGQ_NUM_COORD; cIndex++) {
                    if (!allocOrderFound[cIndex]) {
                        allocOrder[allocOrderIndex--] = cIndex;
                    }
                }
            }

/*
 *          Calculate the size of the hardware partition from the
 *          min/max task coordinates.  The min/max coordinates
 *          cannot be used to calculate the number of tasks per
 *          node (dimension 6) since they are not necessarily sequential
 *          numbers so calculate that based on the node and task counts.
 */
            for (cIndex = 0; cIndex < (BGQ_NUM_COORD-1); cIndex++) {
                dimSize[cIndex] = maxCoord[cIndex] - minCoord[cIndex] + 1;
            }

            computeNodeCount = dimSize[0] * dimSize[1] * dimSize[2] *
                               dimSize[3] * dimSize[4];

            tasksPerNode = numTasks / computeNodeCount;

            dimSize[5] = tasksPerNode;

#if 1
            printf("Node count:         %d\n", computeNodeCount);
            printf("Task count:         %d\n", numTasks);
            printf("Tasks/node:         %d\n", tasksPerNode);
            printf("Task distribution:  ");

            for (i = 0; i < BGQ_NUM_COORD; i++) {
                printf("%c", (allocOrder[i] < 5 ? ('A'+allocOrder[i]) : 'T'));
            }
            printf("\n");

            printf("Physical geometry:  %d X %d X %d X %d X %d\n",
                   dimSize[0], dimSize[1], dimSize[2],
                   dimSize[3], dimSize[4]);
#endif
    
/*
 *          Since every partition on BG/Q is a 5D 'rectangle', the OS *may*
 *          allocate more nodes than requested in order to guarantee
 *          rectangularness.  This function does not currently have
 *          enough intelligence to figure out a reasonable task geometry
 *          in such a situation so set the task geometry to reflect
 *          this and return.
 */
            if ((tasksPerNode * computeNodeCount) != numTasks) {

                optimalTaskGeometry[X] = -1;
                optimalTaskGeometry[Y] = -1;
                optimalTaskGeometry[Z] = -1;

                return;
            }

/*
 *          Now we have to figure out how to reasonably map a 3D domain
 *          geometry onto the 5D hardware partition (or 6D when you consider
 *          the processors per node as an extra dimension).  We'll do
 *          this by grouping the dimensions into 3 groups, each of
 *          which maps to 1 of the {XYZ} dimensions of the ParaDiS
 *          task geometry.
 *
 *          The goal is to create a 3D geometry as close to cubic
 *          as possible.
 */

            lenSumMin = INT_MAX;

            for (i = 0; i < (BGQ_NUM_COORD - 2); i++) {
                int j, k;
                int lenSum;
                int tasksInGrp[3];
                int firstInGrp[3];
                int dimInGrp[3];

                firstInGrp[0] = 0;
                firstInGrp[1] = i + 1;

                tasksInGrp[0] = dimSize[allocOrder[0]];

                for (k = 1; k <= i; k++) {
                    tasksInGrp[0] *= dimSize[allocOrder[k]];
                }

                for (j = i + 1; j < (BGQ_NUM_COORD - 1); j++) {

                    firstInGrp[2] = j + 1;

                    tasksInGrp[1] = dimSize[allocOrder[firstInGrp[1]]];

                    for (k = firstInGrp[1]+1; k < firstInGrp[2]; k++) {
                        tasksInGrp[1] *= dimSize[allocOrder[k]];
                    }

                    tasksInGrp[2] = dimSize[allocOrder[firstInGrp[2]]];

                    for (k = firstInGrp[2]+1; k <= 5; k++) {
                        tasksInGrp[2] *= dimSize[allocOrder[k]];
                    }

                    dimInGrp[0] = i + 1;
                    dimInGrp[1] = j - dimInGrp[0] + 1;
                    dimInGrp[2] = 6 - (dimInGrp[0] + dimInGrp[1]);

/*
 *                  If this configuration is closer to cubic than anything
 *                  yet seen, it becomes the optimal configuration
 */
                    lenSum = tasksInGrp[0] + tasksInGrp[1] + tasksInGrp[2];

                    if (lenSum < lenSumMin) {
                        int optimalDimIn3DGrp[3];
/*
 *             FIX ME!  It is likely that some of the dimensions will
 *                      not be tori.  We should add some sort of penalty
 *                      score for each grouping of dimensions to take
 *                      into consideration the increased number of hops
 *                      communications will require for dimensions that
 *                      are not tori.  Then if two configurations are
 *                      found to be equally optimal so far as the
 *                      shape of the ParaDiS 3D task geometry, the
 *                      penalty score of each can ba assessed to choose
 *                      between them.
 */
                        lenSumMin = lenSum;

                        optimalDimIn3DGrp[0] = dimInGrp[0];
                        optimalDimIn3DGrp[1] = dimInGrp[1];
                        optimalDimIn3DGrp[2] = dimInGrp[2];

                        optimal3DGeometry[0] = tasksInGrp[0];
                        optimal3DGeometry[1] = tasksInGrp[1];
                        optimal3DGeometry[2] = tasksInGrp[2];
                    }
                }
            }

/*
 *          Save the optimal task geometry.  NOTE: Only domain 0 knows
 *          the optimal geometry for now.  Domain 0 will later determine
 *          if the task geometry needs to change, and will pass the
 *          selected task geometry to remote domains when the param
 *          sturct gets broadcast out.
 */
            VECTOR_COPY(optimalTaskGeometry, optimal3DGeometry);

/*
 *          Don't forget to free the temporary buffer before we're done.
 */
            free(coordinateList);

        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     ResetTaskGeometry_BGQ
 *      Description:  For BG/Q systems, we may not know the physical geometry
 *                    of the 5D partition until run-time, which means the
 *                    logical domain decomposition specified by the user
 *                    may not match the hardware partition.  Before this
 *                    function is called a call must be made to the
 *                    GetOptimal3DTaskGeometry() function.  Once that is
 *                    done, this function will verify whether an optimal
 *                    task geometry was obtained and based on the value
 *                    of the <taskMappingMode> control parameter determine
 *                    whether to reset the task geometry to match the
 *                    physical partition, do nothing, or abort.
 *
 *-------------------------------------------------------------------------*/
static void ResetTaskGeometry_BGQ(Home_t *home)
{
        int     haveOptimalGeometry, geometryIsMatch;
        Param_t *param;

        param = home->param;

/*
 *      If the user wants to use the logical task geometry specified
 *      in the control parameter file regardless of the shape of the
 *      physical hardwar partition, there's no need to do anything
 *      more.
 */
        if (param->taskMappingMode == TASK_MAP_MISMATCH_USE_LOGICAL) {
            return;
        }

/*
 *      Check if we were able to obtain an "optimal" task geometry
 *      based on the shape of the physical hardware partition.
 */
        haveOptimalGeometry = (param->optimalTaskGeometry[X] > 0) &&
                              (param->optimalTaskGeometry[Y] > 0) &&
                              (param->optimalTaskGeometry[Z] > 0);

        geometryIsMatch = (param->optimalTaskGeometry[X] == param->nXdoms) &&
                          (param->optimalTaskGeometry[Y] == param->nYdoms) &&
                          (param->optimalTaskGeometry[Z] == param->nZdoms);

/*
 *      If we were unable to determine an "optimal" task geometry earlier
 *      on, and the caller wants to abort if the task geometry he specified
 *      does not match the hardware partition, abort.
 */
        if (!haveOptimalGeometry) {
            switch(param->taskMappingMode) {
                case TASK_MAP_MISMATCH_ABORT:
                case TASK_MAP_MISMATCH_USE_PHYSICAL:

                    Fatal("Unable to determine optimal task geometry for "
                          "this partition.\n  Try setting the "
                          "<taskMappingMode> control parameter to 1.");
                    break;
            }
        } else if (geometryIsMatch) {
/*
 *          The user-specified task geometry already matches the optimal
 *          geometry, so nothing more to do.
 */
            return;
        } else {
/*
 *          The user-supplied geometry does not match the optimal geometry
 *          for this hardware partition, so either abort or reset the
 *          geometry.
 */
            switch(param->taskMappingMode) {
                case TASK_MAP_MISMATCH_ABORT:

                    Fatal("logical task geometry %dX%dX%d mismatch with the "
                          "physical task\n    geometry %dX%dX%d: Aborting!",
                           param->nXdoms, param->nYdoms, param->nZdoms,
                           param->optimalTaskGeometry[X],
                           param->optimalTaskGeometry[Y],
                           param->optimalTaskGeometry[Z]);

                    break;

                case TASK_MAP_MISMATCH_USE_PHYSICAL:

                    printf("Warning: logical domain geometry %d X %d X %d is "
                           "inconsistent with\n    physical geometry of "
                           "%d X %d X %d.\n    Resetting logical geometry to "
                           "be consistent with physical!\n",

                    param->nXdoms, param->nYdoms, param->nZdoms,
                    param->optimalTaskGeometry[X],
                    param->optimalTaskGeometry[Y],
                    param->optimalTaskGeometry[Z]);

                    param->nXdoms = param->optimalTaskGeometry[X];
                    param->nYdoms = param->optimalTaskGeometry[Y];
                    param->nZdoms = param->optimalTaskGeometry[Z];

                    break;
            }
        }

        return;
}
#endif  // _BGQ 

#ifdef _BGP
/*---------------------------------------------------------------------------
 *
 *      Function:     ResetTaskGeometry_BGP
 *      Description:  For BG/P systems, we may not know the physical geometry
 *                    of the 3D partition until run-time, which means the
 *                    logical domain decomposition specified by the user
 *                    may not match the hardware partition.  This function
 *                    will determine the geometry of the hardware partition,
 *                    and compare it to the user-specified domain
 *                    decomposition.
 *
 *                    If the two are consistent, all is okay.  If not, the
 *                    code will either ignore the difference, abort, or try to
 *                    reset the logical domain decomposition to be consistent
 *                    with the underlying hardware partition.  Which action is
 *                    taken depends on the value of the <taskMappingMode>
 *                    control file parameter.
 *
 *-------------------------------------------------------------------------*/
static void ResetTaskGeometry_BGP(Home_t *home)
{
        int  i, length, offset;
        int  taskCount, nodeCount, tasksPerNode, sequentialPerNode;
        int  allocationOrder[3];
        int  physicalGeom[3], physicalTasks[3], logicalGeom[3];
        char taskDistributionSpec[5], *strPtr;
        Param_t *param;
        _BGP_Personality_t personality;

        param = home->param;

/*
 *      If the user requested to use the logical task mapping regardless of
 *      the geometry of the physical partition on which we're running, then
 *      there's no reason to check the physical partition geometry.
 */
        if (param->taskMappingMode == TASK_MAP_MISMATCH_USE_LOGICAL) {
            return;
        }

        logicalGeom[X] = param->nXdoms;
        logicalGeom[Y] = param->nYdoms;
        logicalGeom[Z] = param->nZdoms;

/*
 *      Get the physical domain geoemtry
 */
        Kernel_GetPersonality(&personality, sizeof(personality));

        physicalGeom[X] = personality.Network_Config.Xnodes;
        physicalGeom[Y] = personality.Network_Config.Ynodes;
        physicalGeom[Z] = personality.Network_Config.Znodes;

/*
 *      Search for an environment variable specifying the method of MPI
 *      task mapping.  If not found, use the default of "XYZT".
 *
 *      Don't know if it's needed, but convert the string to upper case
 *      to simplify handling it.
 */
        strPtr = getenv(TASK_MAP_ENV_VAR);
        if (strPtr == (char *)NULL) {
            strcpy(taskDistributionSpec, TASK_DISTRIBUTION_DEFLT);
        } else {
            taskDistributionSpec[4] = 0;
            strncpy(taskDistributionSpec, strPtr, 4);
        }

        length = strlen(taskDistributionSpec);

        for (i = 0; i < length; i++) {
            taskDistributionSpec[i] = toupper((int)taskDistributionSpec[i]);
        }

/*
 *      Figure out if multiple tasks on a node are assigned sequentially
 *      and determine how tasks are allocate in each dimension (i.e. which
 *      dimension varies fastest, which slowest)
 */
        if (taskDistributionSpec[0] == 'T') {
            sequentialPerNode = 1;
            offset = 1;
        } else {
            sequentialPerNode = 0;
            offset = 0;
        }

        for (i = 0; i < 3; i++) {
            allocationOrder[i] = taskDistributionSpec[i+offset] - 'X';
        }

/*
 *      If the task count is greater than the node count, it means we're
 *      going to assign multiple tasks per node... so figure out how many
 *      tasks per node and set up the 'task geometry' array.  The physical
 *      task array is set to the number of tasks per dimension and adjusts
 *      the value for the fastest changing dimension to account for multiple
 *      tasks per node.
 */
        taskCount = home->numDomains;
        nodeCount = physicalGeom[X] * physicalGeom[Y] * physicalGeom[Z];

        if ((taskCount != nodeCount) && (sequentialPerNode == 0)) {
            Fatal("Can't figure out how to map logical geometry %dx%dx%d\n"
                  "to physical geometry %dx%dx%d using %s\n",
                  logicalGeom[X], logicalGeom[Y], logicalGeom[Z],
                  physicalGeom[X], physicalGeom[Y], physicalGeom[Z],
                  taskDistributionSpec);
        }

        physicalTasks[0] = physicalGeom[0];
        physicalTasks[1] = physicalGeom[1];
        physicalTasks[2] = physicalGeom[2];

        tasksPerNode = taskCount / nodeCount;
        physicalTasks[allocationOrder[0]] *= tasksPerNode;

/*
 *      If the physical geometry is a match for the logical, there's no
 *      problem and we're done.
 *
 *      Note: ParaDiS assigns domains in Z, Y, X order, so we have to
 *      compare the logical mappings for Z with physical mapping of
 *      the fastest changing dimension, and so on.
 */
        if ((logicalGeom[Z] == physicalTasks[allocationOrder[0]]) &&
            (logicalGeom[Y] == physicalTasks[allocationOrder[1]]) &&
            (logicalGeom[X] == physicalTasks[allocationOrder[2]])) {
            return;
        }

#if 1
        printf("Node count:         %d\n", nodeCount);
        printf("Task count:         %d\n", taskCount);
        printf("Task distribution:  %s\n", taskDistributionSpec);
        printf("Physical geometry:  %d X %d X %d\n",
               physicalGeom[X], physicalGeom[Y], physicalGeom[Z]);

        printf("Physical tasks:     %d X %d X %d\n",
               physicalTasks[X], physicalTasks[Y], physicalTasks[Z]);

        printf("Logical geometry:   %d X %d X %d\n",
               logicalGeom[X], logicalGeom[Y], logicalGeom[Z]);
#endif

/*
 *      Logical geometry does not match physical geometry.  If the user
 *      specified the logical mapping must be used, we left this function
 *      way before this point, which leaves us with two options:
 *          a) modify the logical mapping to match the physical mapping
 *          b) abort if the logical mapping does not match the physical
 */
        switch (param->taskMappingMode) {
            case TASK_MAP_MISMATCH_USE_PHYSICAL:
                printf("Warning: logical domain geometry %d X %d X %d is "
                       "inconsistent with\n    physical task geometry of"
                       "%d X %d X %d.\n    Resetting logical geometry to be "
                       "consistent with physical!\n",
                       logicalGeom[X], logicalGeom[Y], logicalGeom[Z],
                       physicalTasks[X], physicalTasks[Y], physicalTasks[Z]);
                logicalGeom[Z] = physicalTasks[allocationOrder[0]];
                logicalGeom[Y] = physicalTasks[allocationOrder[1]];
                logicalGeom[X] = physicalTasks[allocationOrder[2]];
                break;
            case TASK_MAP_MISMATCH_ABORT:
                Fatal("logical domain geometry %dx%dx%d mismatch with physical"
                      "task \ngeometry %dx%dx%d:  Aborting!\n",
                      logicalGeom[X], logicalGeom[Y], logicalGeom[Z],
                      physicalTasks[X], physicalTasks[Y], physicalTasks[Z]);
                break;
            default:
                Fatal("Unsupported taskMappingMode = %d!\n",
                      param->taskMappingMode);
                break;
        }

/*
 *      Be sure to reset the domain geometry before returning to the caller.
 *      This new domain geometry will then get distributed to remote domains
 *      when task 0 broadcasts the final param structure out.
 */
        param->nXdoms = logicalGeom[X];
        param->nYdoms = logicalGeom[Y];
        param->nZdoms = logicalGeom[Z];

        return;
}
#endif  // _BGP 

/*
 *      Simple dispatch function to invoke the ResetTaskGeometry_*()
 *      function appropriate to the current hardware.
 */
void ResetTaskGeometry(Home_t *home)
{

#ifdef _BGP
        ResetTaskGeometry_BGP(home);
#endif

#ifdef _BGQ
        ResetTaskGeometry_BGQ(home);
#endif

        return;
}

#endif  //  PARALLEL
