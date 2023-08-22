/*-------------------------------------------------------------------------
 *
 *      Module:       NodeTaskGroup.cc
 *      Description:  Contains functions used in identifying the group
 *                    of MPI tasks co-located on the same computer-node.
 *
 *                    NOTE: This is only needed when support for both
 *                    parallel execution and shared memory buffers has been
 *                    compiled into the executable.
 *
 *  Includes public functions:
 *      FreeNodeTaskGroups()
 *      GetNodeTaskGroups()
 *
 *-----------------------------------------------------------------------*/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

#if defined PARALLEL && defined SHARED_MEM
#include <sys/utsname.h>
#endif

#define MAX_NODE_NAME_LEN 64  // set maximum length of compute-node name


/*-------------------------------------------------------------------------
 *
 *  Function:    FreeNodeTaskGroups
 *  Description: Destroys the MPI groups and communicators used for
 *               collective communications within the group of tasks
 *               co-locate on the compute node, and the group of
 *               tasks that are identified as being the 'first' tasks
 *               on each node.
 *
 *               Since this communicator is only used when compiled
 *               for parallel execution with shared memory enabled, if
 *               either is not in use, this function does nothing.
 *
 *------------------------------------------------------------------------*/
void FreeNodeTaskGroups(Home_t *home)
{
#if defined PARALLEL && defined SHARED_MEM
    MPI_Comm_free(&home->commNode);
    MPI_Group_free(&home->groupOnNode);

/*
 *  The first task on each node is also part of an additional group/commuicator
 *  that needs to be be freed
 */
    if (home->isFirstTaskOnNode) {
        MPI_Comm_free(&home->commFirstOnNode);
        MPI_Group_free(&home->groupFirstOnNode);
    }
#endif

    return;
}


/*-------------------------------------------------------------------------
 *
 *  Function:    GetNodeTaskGroups
 *
 *  Description: This function will determine the list of MPI tasks
 *               that are co-located on each node, determine which task
 *               is the first (lowest task ID), set a flag in the
 *               home struct for the first task on the node, and define
 *               an MPI communicator to enable collective communications
 *               within the group of tasks on the node.
 *
 *               If the code has been compiled in serial mode, no
 *               communicator is defined, and the task will group
 *               contains only the single running task.
 *
 *               Currently, this information is only needed when support
 *               for shared memory buffers has been enabled.  Therefore
 *               when shared memory is not enabled, this function takes
 *               no action.
 *
 *------------------------------------------------------------------------*/
void GetNodeTaskGroups(Home_t *home)
{
#if defined PARALLEL && defined SHARED_MEM
    int              i, thisTask, numTasks, numComputeNodes;
    int              firstTaskOnNode;
    int              isFirstTaskOnNode;
    int              nodeNameLen;
    int              nodeTaskListSize, numTasksOnNode;
    int             *nodeTaskList;
    char             thisTaskNodeName[MAX_NODE_NAME_LEN];
    char           (*nodeNameList)[MAX_NODE_NAME_LEN];
    struct utsname   utsInfo;
    MPI_Group        groupWorld;

//
//  Do some initialization and set up some temporary arrays
//
    numTasks = home->numDomains;
    thisTask = home->myDomain;

    thisTaskNodeName[MAX_NODE_NAME_LEN-1] = 0;

    numTasksOnNode  = 0;
    firstTaskOnNode = 0;
    nodeTaskListSize = 64;

    nodeTaskList = (int *)malloc(nodeTaskListSize * sizeof(int));

    nodeNameList = (char (*)[MAX_NODE_NAME_LEN])malloc(numTasks *
                                                       MAX_NODE_NAME_LEN);

//
//  Get the name of this task's compute-node.  If the node name doesn't
//  fit in the buffer, any comparison to other node names would use
//  truncated strings which could result in incorrect matches.  To avoid
//  that, set the node name to a short unique value to force the task
//  to assume it is not on the same node as it's neighbors.
//
    (void)uname(&utsInfo);
    nodeNameLen = strlen(utsInfo.nodename);

    if ((nodeNameLen > MAX_NODE_NAME_LEN-1) || (nodeNameLen < 1)) {
        sprintf(thisTaskNodeName, "%d:node_name_error", thisTask);
    } else {
        strncpy(thisTaskNodeName, utsInfo.nodename, MAX_NODE_NAME_LEN-1);
    }

//
//  Gather the compute-node name for every task into a single array
//
    MPI_Allgather(thisTaskNodeName, MAX_NODE_NAME_LEN, MPI_CHAR,
                  nodeNameList, MAX_NODE_NAME_LEN, MPI_CHAR, MPI_COMM_WORLD);

//
//  Loop through the list of compute-node names for every task.  Save the
//  ID of each task running on the local compute-node, get a total count
//  of the tasks running on the local compute-node, and set a flag if the
//  current task is the first (lowest task ID) task on the compute-node.
//
    for (i = 0; i < numTasks; i++) {

        if (strcmp(thisTaskNodeName, nodeNameList[i]) == 0) {

            if (numTasksOnNode > nodeTaskListSize) {
                nodeTaskListSize += 64;
                nodeTaskList = (int *)realloc(nodeTaskList,
                                              nodeTaskListSize * sizeof(int));
            }

            nodeTaskList[numTasksOnNode++] = i;
        }
    }

    home->isFirstTaskOnNode = (thisTask == nodeTaskList[0]);

//
//  Build an MPI communicator specific to only the tasks on the
//  current compute-node.
//
    MPI_Comm_group(MPI_COMM_WORLD, &groupWorld);

    MPI_Group_incl(groupWorld, numTasksOnNode, nodeTaskList,
                   &home->groupOnNode);

    MPI_Comm_create(MPI_COMM_WORLD, home->groupOnNode,
                    &home->commNode);

//
//  Now we need to create a communicator that encompasses the
//  first MPI task on all compute nodes.
//
//  First all the tasks need the identify all "first-on-node" tasks.
//
    home->localReduceBuf[home->myDomain] = home->isFirstTaskOnNode;

    MPI_Allgather(&home->isFirstTaskOnNode, 1, MPI_INT, home->globalReduceBuf,
                  1, MPI_INT, MPI_COMM_WORLD);

//
//  Each task has an array of integer flags (1 per task) indicating if the
//  associated task is the first on a compute node.  Overwrite the array
//  with the indices of the "first-on-node" tasks.
//
    for (numComputeNodes = 0, i = 0; i < numTasks; i++) {
        if (home->globalReduceBuf[i] == 1) {
            home->globalReduceBuf[numComputeNodes++] = i;
        }
    }

//
//  Now we create a communicator that encompasses all the 
//  "first-on-node" tasks.  Re-use <groupWorld> which was
//  created above
//
    MPI_Group_incl(groupWorld, numComputeNodes, home->globalReduceBuf,
                   &home->groupFirstOnNode);

    MPI_Comm_create(MPI_COMM_WORLD, home->groupFirstOnNode,
                    &home->commFirstOnNode);

//
//  Free up the temporary arrays.  The newly created MPI groups and
//  communicators will be freed at program termination.
//
    free(nodeNameList);
    nodeNameList = (char (*)[MAX_NODE_NAME_LEN])NULL;

    free(nodeTaskList);
    nodeTaskList = (int *)NULL;

#else  // Either PARALLEL or SHARED_MEM is not defined

    home->isFirstTaskOnNode = 1;

#endif  // end if defined PARALLEL && defined SHARED_MEM

    return;
}
