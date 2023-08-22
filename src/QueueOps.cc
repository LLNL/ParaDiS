/*****************************************************************************
 *
 *  Module      : QueueOps
 *  Description : Contains the various queue operations on the native, ghost,
 *                and free queues.
 *
 ****************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "ParadisThread.h"
#include "Node.h"
#include "QueueOps.h"

/*---------------------------------------------------------------------------
 *
 * Function    : PopFreeNodeQ
 * Description : dequeue the node at the top of the freeNodeQ and return it
 *               to the caller. If there are no free nodes, allocate a new
 *               block of nodes, queue them all on the freeNodeQ, and then
 *               pop the first one
 *
 *---------------------------------------------------------------------------*/
Node_t *PopFreeNodeQ(Home_t *home)
{
	int		i;
	NodeBlock_t	*nodeBlock;
	Node_t		*currNode, *freeNode;

	if (home->freeNodeQ == 0) {

		nodeBlock = (NodeBlock_t *) malloc(sizeof(NodeBlock_t));
		nodeBlock->next = home->nodeBlockQ;
		home->nodeBlockQ = nodeBlock;

		nodeBlock->nodes = (Node_t *)calloc(1, NODE_BLOCK_COUNT *
						    sizeof(Node_t));

		currNode = nodeBlock->nodes;
		home->freeNodeQ = currNode;

		for (i = 1; i < NODE_BLOCK_COUNT; i++) {
			INIT_LOCK(&currNode->nodeLock);
			currNode->next = currNode + 1;
			currNode++;
		}

		currNode->next = 0;    /* last free node */
		INIT_LOCK(&currNode->nodeLock);
		home->lastFreeNode = currNode;
	}

/*
 *	Dequeue the first free node and return it to the caller
 */
	freeNode = home->freeNodeQ;
	home->freeNodeQ = freeNode->next;

	return(freeNode);
}


/*--------------------------------------------------------------------------
 *
 *	Function:	PushFreeNodeQ
 *	Description:	Push a node onto the top of the free queue
 *
 *-------------------------------------------------------------------------*/
void PushFreeNodeQ(Home_t *home, Node_t *node)
{
	node->next = home->freeNodeQ;
	home->freeNodeQ = node;

	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	PushNativeNodeQ
 *	Description:	Push a node onto the top of the native queue
 *
 *-------------------------------------------------------------------------*/
void PushNativeNodeQ(Home_t *home, Node_t *node)
{
	node->next = home->nativeNodeQ;
	home->nativeNodeQ = node;

	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	PushGhostNodeQ
 *      Description:    Push a node onto the top of the ghost queue, then
 *                      place it onto the ghost node array as well.  We
 *                      maintain the information in two ways because the
 *                      original ghost node queue was good for insertions
 *                      and deletions, but was not good for loops through
 *                      all the ghost nodes (couldn't be threaded).
 *
 *                      NOTE: This function accesses global ghost node
 *                      data structures without any locking, so the
 *                      caller must insure the following:
 *
 *                        1) Only 1 thread may call this function at a time
 *                        2) None of the threads will be executing a portion
 *                           of the code that loops through the ghost nodes
 *                           while another thread is updating the ghost node
 *                           data in this function!
 *
 *-------------------------------------------------------------------------*/
void PushGhostNodeQ(Home_t *home, Node_t *node)
{

/*
 *      First update the queue
 */
        if (home->ghostNodeQ == 0) home->lastGhostNode = node;
        node->next = home->ghostNodeQ;
        home->ghostNodeQ = node;

/*
 *      Now update the array and count of ghost node pointers.  Increase
 *      the array size if needed.
 *
 *      Note: We never reduce the size of this array, so it will be rare
 *            that expensive calls to reallocate the array will have
 *            to be done.
 */
        if (home->ghostNodeCount >= home->ghostNodeListSize) {
            int size;
            home->ghostNodeListSize += 1000;
            size = home->ghostNodeListSize * sizeof(Node_t *);
            home->ghostNodeList = (Node_t **)realloc(home->ghostNodeList, size);
        }
        home->ghostNodeList[home->ghostNodeCount] = node;
        home->ghostNodeCount += 1;

	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	RecycleGhostNodes
 *	Description:	Put all the nodes in the ghostNodeQ back onto
 *			the freeNodeQ and reset the ghostNodeQ to empty
 *
 *-------------------------------------------------------------------------*/
void RecycleGhostNodes(Home_t *home)
{
	if (home->ghostNodeQ == NULL) return;  /* nothing to do */

	if (home->freeNodeQ == NULL) {
		home->freeNodeQ = home->ghostNodeQ;
	} else {
		home->lastFreeNode->next = home->ghostNodeQ;
	}

	home->lastFreeNode = home->lastGhostNode;
	home->lastGhostNode = NULL;
	home->ghostNodeQ = NULL;

        home->ghostNodeCount = 0;
	return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	RemoveNodeFromCellQ
 *	Description:	Remove the specified node from the queue for
 *                      the cell containing the node.  Typically only
 *                      needed when deleting a node during topological
 *                      changes.
 *
 *-------------------------------------------------------------------------*/
void RemoveNodeFromCellQ(Home_t *home, Node_t *node)
{
        Cell_t *cell;
        Node_t *prevNode = (Node_t *)NULL;
        Node_t *nextNode;

        if (node == (Node_t *)NULL) return;
        if (node->cellIdx < 0) return;

        cell = LookupCell(home, node->cellIdx);

        if (node == cell->nodeQ) {
            cell->nodeQ = node->nextInCell;
            cell->nodeCount--;
        } else {
            nextNode = cell->nodeQ;
            while (nextNode != (Node_t *)NULL) {
                if (nextNode == node) {
                    prevNode->nextInCell = node->nextInCell;
                    cell->nodeCount--;
                    return;
                }
                prevNode = nextNode;
                nextNode = nextNode->nextInCell;
            }
        }

        return;
}


/*--------------------------------------------------------------------------
 *
 *	Function:	RemoveNodeFromCell2Q
 *	Description:	Remove the specified node from the queue for
 *                      the cell2 containing the node.  Typically only
 *                      needed when deleting a node during topological
 *                      changes.
 *
 *-------------------------------------------------------------------------*/
void RemoveNodeFromCell2Q(Home_t *home, Node_t *node)
{

	if (home->cell2QentArray == (C2Qent_t *)NULL) {
		return;
	}

        if (node == (Node_t *)NULL) {
		return;
	}

        if ((node->cell2Idx < 0) || (node->cell2QentIdx < 0)) {
		return;
	}

        home->cell2QentArray[node->cell2QentIdx].node = (Node_t *)NULL;
        home->cell2QentArray[node->cell2QentIdx].next = -1;

        return;
}
