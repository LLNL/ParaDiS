/**************************************************************************
 *
 *  Function    : RemoveNode
 *  Description : Unlink a local node from its two neighbors, 
 *                return it to the free Queue, and recycle the node index
 *                If node is not local, just unlink it and zero its 
 *                nodeKeys entry.
 *
 **************************************************************************/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"

/*
 *  Modified by Wei Cai 9/7/2001, so that the neighbor nodes
 *  do not necessarily occupy the nbrTag[] slots sequentially
 */
   
void RemoveNode(Home_t *home, Node_t *node, int Log)
{
    int             domain, index;
    Node_t         *nbr1, *nbr2;
    RemoteDomain_t *remDom;


    domain = node->myTag.domainID;
    index = node->myTag.index;

/*
 *   It can happen that a node is pinned and has two neighbors.
 *   Prevent removal of such a node.
 */
    if (HAS_ANY_OF_CONSTRAINTS(node->constraint, PINNED_NODE)) {
       return;
    }  


    if (Log) {
        AddOpRemoveNode(home, &node->myTag);
    }

/*
 * For local nodes, all the nodes arms should have been removed
 * prior to the invocation of this function, so just free the
 * node structure and return.
 */
    if ((domain == home->myDomain) && (node->numNbrs == 0)) {
        FreeNode(home, index);
        return;
    }

/*
 *  There are situations where this function is invoked (via FixRemesh())
 *  to remove a ghost node which whose arms have not been removed.  In
 *  this case, just change the connectivity around it and zero
 *  out its remDom->nodeKeys entry. The node struct itself will be
 *  freed up with the rest of the ghosts at the end of the cycle
 *  (another reason why the ghost queue is not a reliable list of active
 *  ghost nodes)
 */
    if (domain != home->myDomain) {

        if (node->numNbrs == 0) {
            return;
        }
           
        if (node->numNbrs == 2) {

            nbr1 = GetNeighborNode(home, node, 0);
            nbr2 = GetNeighborNode(home, node, 1);

/*
 *          Since MeshCoarsen() does not notify neighboring domains of
 *          operations involving only nodes which have neither off-domain
 *          connections nor connections to other nodes with off-domain
 *          connections, it is possible to have a situation in which a
 *          request is received to remove a ghost node which the current
 *          domain still believes to have connected segments.  In this
 *          case, it is safe to simply cut whatever connections it can
 *          and go on.
 */
#if 0
            if ((nbr1 == (Node_t *)NULL) ||
                (nbr2 == (Node_t *)NULL)) {
                printf("Task %d: RemoveNode() - node (%d,%d), "
                       "nbr1 (%d,%d) ptr %s, nbr2 (%d,%d) ptr %s\n",
                       home->myDomain,
                       node->myTag.domainID, node->myTag.index,
                       node->nbrTag[0].domainID,
                       node->nbrTag[0].index,
                       (nbr1 == (Node_t *)NULL ? "NULL" : "OKAY"),
                       node->nbrTag[1].domainID,
                       node->nbrTag[1].index,
                       (nbr2 == (Node_t *)NULL ? "NULL" : "OKAY"));

                if ((nbr1 != (Node_t *)NULL) ||
                    (nbr2 != (Node_t *)NULL)) {
                    Fatal("Task %d: Neighbor not found at %s "
                          "line %d\n", home->myDomain,
                          __FILE__, __LINE__);
                }
            } 
#endif

            if (nbr1 != (Node_t *)NULL) {
                ChangeConnection(home, nbr1, &node->myTag,
                                 &node->nbrTag[1], 0);
            }

            if (nbr2 != (Node_t *)NULL) {
                ChangeConnection(home, nbr2, &node->myTag,
                                 &node->nbrTag[0], 0);
            }

            remDom = home->remoteDomainKeys[domain];
            remDom->nodeKeys[index] = 0;
            return;
        }

        if ((node->numNbrs != 0) && (node->numNbrs != 2)) {
/*
 *          Wei Cai, 11/20/2003 We should be careful here in principle
 *          this should not be allowed but in EventHandle.c we call
 *          RemoveNode to remove a redundant node4, all of whose arms
 *          have been directed to node2 already in this case we really
 *          should not rediect the arms of the neighbors of node4, instead
 *          we should simply do FreeNode(node4).  Here we are lucky that
 *          things work out because when calling ChangeConnection, nbr1
 *          and nbr2 already do not have node4 as their neighbors
 */
            /*Fatal("RemoveNode not with 2 arms!\n");*/
            return;
        }
    }

    Fatal("RemoveNode: Removing local node with non-zero arm count!\n");

    return;
}
