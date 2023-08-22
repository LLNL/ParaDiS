/***************************************************************************
 *
 *      Function:    GetNewNativeNode
 *      Description: Pop a free Node_t from the free node Queue, and assign
 *                   it a free tag. NOTE: This routine does not allocate
 *                   any arrays specific to the number of neighbors. The
 *                   caller must use AllocNodeArms() to do this.
 *
 **************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"
#include "QueueOps.h"

Node_t *GetNewNativeNode(Home_t *home)
{
        int     newIdx;
        Node_t *newNode;

        newNode = PopFreeNodeQ(home);
        newIdx  = GetFreeNodeTag(home);

        home->nodeKeys[newIdx] = newNode;

        newNode->myTag.domainID = home->myDomain;
        newNode->myTag.index    = newIdx;
        newNode->cellIdx        = -1;
        newNode->cell2Idx       = -1;
        newNode->cell2QentIdx   = -1;
/*
 *      Explicitly zero out velocity so we don't end up
 *      with garbage when we are calculating the velocity
 *      delta between timesteps.
 */
        newNode->vX = 0.0;
        newNode->vY = 0.0;
        newNode->vZ = 0.0;

        newNode->oldvX = 0.0;
        newNode->oldvY = 0.0;
        newNode->oldvZ = 0.0;

        newNode->flags = 0;

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
        newNode->multiNodeLife = 0;
#endif

#ifdef _ARLFEM
        newNode->surfaceNorm[0] = 0.0;
        newNode->surfaceNorm[1] = 0.0;
        newNode->surfaceNorm[2] = 0.0;
#endif

        return(newNode);
}
