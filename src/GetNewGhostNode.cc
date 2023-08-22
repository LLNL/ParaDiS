/***************************************************************************
 *
 *  Function    : GetNewGhostNode
 *  Description : Get a free node, and assign it to the specified domain 
 *                and index. Extend the RemoteDomain's nodeKeys table, if
 *                necessary. Also, queue the node onto the ghost node queue
 *
 ***************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "QueueOps.h"
#include "Util.h"

Node_t *GetNewGhostNode(Home_t *home, int domain, int index)
{
        int            i;
        Node_t         *newNode;
        RemoteDomain_t *remDom;

        newNode = PopFreeNodeQ (home);
        remDom  = home->remoteDomainKeys[domain];

        if (!remDom) {
            Fatal("GetNewGhostNode: domain %d is not a neighbor of domain %d",
                  domain, home->myDomain);
        }

/*
 *      If index is beyond current limit on remDom->nodeKeys, expand it big
 *      enough to include index, and initialize new elements to zero
 */
        if (index >= remDom->maxTagIndex) {
            remDom->nodeKeys = (Node_t **)realloc(remDom->nodeKeys, 
                                                  (index+1)*sizeof(Node_t *));

        for (i = remDom->maxTagIndex; i < index; i++)
            remDom->nodeKeys[i] = 0;
            remDom->maxTagIndex = index + 1;
        }

        remDom->nodeKeys[index] = newNode;
        newNode->myTag.domainID = domain;
        newNode->myTag.index    = index;
        newNode->cellIdx        = -1;
        newNode->cell2Idx       = -1;
        newNode->cell2QentIdx   = -1;

        newNode->vX = 0.0;
        newNode->vY = 0.0;
        newNode->vZ = 0.0;

        newNode->oldvX = 0.0;
        newNode->oldvY = 0.0;
        newNode->oldvZ = 0.0;

        newNode->flags = 0;

        PushGhostNodeQ(home, newNode);

#ifdef _ARLFEM
        newNode->surfaceNorm[0] = 0.0;
        newNode->surfaceNorm[1] = 0.0;
        newNode->surfaceNorm[2] = 0.0;
#endif

        return(newNode);
}
