/****************************************************************************
 *
 *  Function:    FixRemesh
 *  Description: This function is called by each domain after
 *               a set of topology changes has been done (i.e.
 *               remesh, collision handling or SplitMultiNodes()).
 *               The function uses the list of operations provided
 *               by each domain to update its local topology with
 *               the necessary changes.
 *
 *               The opList received from each neighbor has been
 *               stored in remDom->inBuf, with length in
 *               remDom->inBufLen, by CommSendRemsh.
 *
 *****************************************************************************/
#include <stdio.h>
#include <string.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"


void FixRemesh(Home_t *home)
{
    int      idom;
    Param_t *param;

/*
 *  Loop through the neighbor domains
 */
    param = home->param;

    for (idom = 0; idom < home->remoteDomainCount; idom++) {
        int             offSet, domIdx, recvOpListLen;
        char           *recvOpListBuf;
        Node_t         *node1, *node2;
        RemoteDomain_t *remDom;

        domIdx = home->remoteDomains[idom];
        remDom = home->remoteDomainKeys[domIdx];

        recvOpListBuf = (char *)remDom->inBuf;
        recvOpListLen = remDom->inBufLen;

/*
 *      Loop through the operations from this neighbor
 */
        offSet = 0;

        while (offSet < recvOpListLen) {
            int opType;

            memcpy(&opType, &recvOpListBuf[offSet], sizeof(int));

            if (opType == REM_OP_CHANGE_CONN) {
                int               opLen = sizeof(RemOpChangeConn_t);
                RemOpChangeConn_t opBuf;

/*
 *              Change the node1->node2 link so it becomes
 *              a node1->node3 link instead.
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != NULL) {
                    ChangeConnection (home, node1, &opBuf.tag2, &opBuf.tag3, 0);
                }

            } else if (opType == REM_OP_INSERT_ARM) {
                int              opLen = sizeof(RemOpInsertArm_t);
                RemOpInsertArm_t opBuf;

/*
 *              Add an arm from node1 to node2
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != NULL) {
                    InsertArm(home, node1, &opBuf.tag2,
                              opBuf.burg[0], opBuf.burg[1], opBuf.burg[2],
                              opBuf.plane[0],opBuf.plane[1],opBuf.plane[2], 0);
                }

            } else if (opType == REM_OP_SPLIT_NODE) {
                int              opLen = sizeof(RemOpSplitNode_t);
                real8            cellSize[3], cellIndex[3];
                RemOpSplitNode_t opBuf;

/*
 *              Create a ghost node with the necessary domain
 *              and index, and set the coordinates and velocity
 *              appropriately.  Note: the node initially has
 *              no arms, they will be added later.  
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node2 = GetNewGhostNode(home,
                                        opBuf.tag2.domainID,
                                        opBuf.tag2.index);
                FreeNodeArms(node2);

                SET_CONSTRAINTS(node2->constraint, UNCONSTRAINED);

                node2->x = opBuf.newCoord[0];
                node2->y = opBuf.newCoord[1];
                node2->z = opBuf.newCoord[2];

                node2->oldx = opBuf.newCoord[0];
                node2->oldy = opBuf.newCoord[1];
                node2->oldz = opBuf.newCoord[2];

                node2->vX = opBuf.newVelocity[0];
                node2->vY = opBuf.newVelocity[1];
                node2->vZ = opBuf.newVelocity[2];

/*
 *              Insure the node is assigned to the proper cell
 */
                cellSize[X] = param->Lx / param->nXcells;
                cellSize[Y] = param->Ly / param->nYcells;
                cellSize[Z] = param->Lz / param->nZcells;

                cellIndex[X] = (int)(node2->x-param->minSideX) /
                               (int)cellSize[X];

                cellIndex[Y] = (int)(node2->y-param->minSideY) /
                               (int)cellSize[Y];

                cellIndex[Z] = (int)(node2->z-param->minSideZ) /
                               (int)cellSize[Z];

                cellIndex[X] = MAX(0, cellIndex[X]);
                cellIndex[Y] = MAX(0, cellIndex[Y]);
                cellIndex[Z] = MAX(0, cellIndex[Z]);

                cellIndex[X] = MIN(param->nXcells,cellIndex[X]);
                cellIndex[Y] = MIN(param->nYcells,cellIndex[Y]);
                cellIndex[Z] = MIN(param->nZcells,cellIndex[Z]);

                cellIndex[X]++;
                cellIndex[Y]++;
                cellIndex[Z]++;

                node2->cellIdx = EncodeCellIdx(home,
                                               cellIndex[X],
                                               cellIndex[Y],
                                               cellIndex[Z]);

            } else if (opType == REM_OP_MARK_FORCE_OBSOLETE) {
                int  opLen = sizeof(RemOpMarkForceObsolete_t);
                RemOpMarkForceObsolete_t opBuf;

/*
 *              Mark a node so that new force and velocity
 *              values will be computed for the node.  Probably
 *              because topology changes in another domain have
 *              affected this node.
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != (Node_t *)NULL) {
                    node1->flags |= NODE_RESET_FORCES;
                }

            } else if (opType == REM_OP_RESET_SEG_FORCE1) {
                int                   opLen = sizeof(RemOpResetSegForce1_t);
                RemOpResetSegForce1_t opBuf;

/*
 *              Reset the forces at node1 for the specified segment
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != (Node_t *)NULL) {

                    ResetSegForces(home, node1, &opBuf.tag2,
                                   opBuf.f1[0], opBuf.f1[1], opBuf.f1[2], 0);
                }

            } else if (opType == REM_OP_RESET_SEG_FORCE2) {
                int                   opLen = sizeof(RemOpResetSegForce2_t);
                RemOpResetSegForce2_t opBuf;

/*
 *              Reset the forces at node1 and node2 for the specified segment
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != (Node_t *)NULL) {

                    ResetSegForces2(home, node1, &opBuf.tag2,
                                    opBuf.f1[0], opBuf.f1[1], opBuf.f1[2],
                                    opBuf.f2[0], opBuf.f2[1], opBuf.f2[2], 0);
                }

            } else if (opType == REM_OP_RESET_COORD) {
                int               opLen = sizeof(RemOpResetCoord_t);
                RemOpResetCoord_t opBuf;

/*
 *              Reset the coordinates of the specified node.
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != (Node_t *)NULL) {

                    node1->x = opBuf.newCoord[0];
                    node1->y = opBuf.newCoord[1];
                    node1->z = opBuf.newCoord[2];
                }

            } else if (opType == REM_OP_REMOVE_NODE) {
                int               opLen = sizeof(RemOpRemoveNode_t);
                RemOpRemoveNode_t opBuf;

/*
 *              Remove the specified node
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != (Node_t *)NULL) {
                    RemoveNode(home, node1, 0);
                }

            } else if (opType == REM_OP_CHANGE_BURG) {
                int               opLen = sizeof(RemOpChangeBurg_t);
                RemOpChangeBurg_t opBuf;

/*
 *              Change the burgers vector of a  node's arm
 *              to the given burgers vector.  If the new
 *              burgers vector is zero, the arm will be
 *              removed.
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                node1 = GetNodeFromTag(home, opBuf.tag1);

                if (node1 != NULL) {

                    ChangeArmBurg(home, node1, &opBuf.tag2, opBuf.newBurg[0],
                                  opBuf.newBurg[1], opBuf.newBurg[2],
                                  opBuf.newPlane[0], opBuf.newPlane[1],
                                  opBuf.newPlane[2], 0, DEL_SEG_NONE);
                }

            } else if (opType == REM_OP_RESET_PLANE) {
                int               opLen = sizeof(RemOpResetPlane_t);
                RemOpResetPlane_t opBuf;

/*
 *              Reset the glide plane for a specific segment.
 */
                memcpy(&opBuf, &recvOpListBuf[offSet], opLen);
                offSet += opLen;

                ResetGlidePlane(home, opBuf.newPlane, &opBuf.tag1,
                                &opBuf.tag2, 0);

            } else {
                    Fatal("FixRemesh : invalid Remesh operation");
            }

        }  // end while (offSet < recvOpListLen)

/*
 *      We're done with this domain's opList; release it
 */
        free(remDom->inBuf);
    }

    return;
}
