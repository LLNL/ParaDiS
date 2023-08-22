/*****************************************************************************
 *
 *      Module:       SMNEvalSplitForces
 *
 *      Description:  This function examines all nodes with at least
 *                    four arms and decides if the node should be split
 *                    and some of the node's arms moved to a new node.
 *                    If it is determined this is necessary, the function
 *                    will invoke the lower level SplitNode() function
 *                    providing the list of arms to be moved to the new
 *                    node.
 *
 *                    NOTE: All topology changes due to this function
 *                    must be communicate to remote domains.  Currently
 *                    these topology changes are combined with operations
 *                    from the collision handling and distributed to
 *                    remote daomins via CommSendRemesh()/FixRemesh()
 *                    calls in ParadisStep().
 *
 *****************************************************************************/

#include "mpi_portability.h"

#include "Home.h"
#include "SMN.h"


static int dbgDom;

/*
 *      Use a predfined of list of the unique ways it is possible to
 *      split 2 or more arms from nodes of up to 15 arms.  We predefine
 *      these values since they are constant and there's no point in
 *      recomputing the values every time we need them.
 */
extern int POSSIBLE_SPLITS[16];


/*---------------------------------------------------------------------------
 *
 *      Function:     SMNEvalSplitForces
 *
 *      Parameters:
 *          nodeCount  IN:  Number of local multi-nodes in <nodeInfo> array
 *          nodeInfo   IN:  Array of <nodeCount> structs containing
 *                          information about each local multi-node to
 *                          be evaluated when doing parallel multi-node
 *                          splitting.
 *
 *-------------------------------------------------------------------------*/
void SMNEvalSplitForces(Home_t *home, int nodeCount, SMN_Info_t *nodeInfo)
{
        int     i;
        real8   splitDist;
        real8   eps;
        Param_t *param;


#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom   = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom   = -1;
#endif
        param     = home->param;

        splitDist = param->rann * 2.0;
        eps = 1.0e-12;

/*
 *      Loop through all native multinodes
 */
        for (i = 0; i < nodeCount; i++) {
            int    k;
            int    isGhostMultiNode;
            int    segID, setID, numNbrs;
            int    totalSets;
            int    bkupNodeCount;
            int    *armList;
            int    **armSets;
            real8  origPos[3], origVel[3];
            Node_t *bkupNodeList;
            Node_t *node;

            node = nodeInfo[i].multiNode;

            isGhostMultiNode = (node->myTag.domainID != home->myDomain);
            
/*
 *          Make a backup copy of the multinode and each of its neighbors
 */
            numNbrs = node->numNbrs;

            bkupNodeCount = numNbrs + 1;
            bkupNodeList = (Node_t *)malloc(bkupNodeCount * sizeof(Node_t));

            BackupNode(home, node, &bkupNodeList[0]);

            for (segID = 0; segID < node->numNbrs; segID++) {
                Node_t *nbrNode;

                nbrNode = GetNodeFromTag(home, node->nbrTag[segID]);
                BackupNode(home, nbrNode, &bkupNodeList[segID+1]);
            }

#ifdef ESHELBY
            SegPartListClear(home);
#endif

/*
 *          Get the total number of ways in which this node can
 *          be split and build arrays of flags we'll use later
 *          for selecting which segments should be split off the
 *          node wich each corresponding split.
 */
            SMNGetArmSets(numNbrs, &totalSets, &armSets);

/*
 *          Allocate the storage for the partial forces for
 *          all possible ways to split this multinode
 */
            nodeInfo[i].numSplits = totalSets;
            nodeInfo[i].status = 0;
            nodeInfo[i].splitInfo =
                    (SMN_Split_t *)calloc(1, totalSets * sizeof(SMN_Split_t));

/*
 *          Make a copy of the original node position and velocity for later
 */
            origPos[X] = node->x;
            origPos[Y] = node->y;
            origPos[Z] = node->z;

            origVel[X] = node->vX;
            origVel[Y] = node->vY;
            origVel[Z] = node->vZ;

            armList = (int *)calloc(1, (numNbrs - 2) * sizeof(int));

/*
 *          For each possible way to split the node, attempt the split,
 *          re-evaluate the velocity of the two nodes, calculate the
 *          portion of the forces for which this domain is responsible,
 *          store the partial forces, and merge the nodes back into one.
 */
            for (setID = 0; setID < totalSets; setID++) {
                int    armIndex;
                int    splitStatus, mobStatus, mergeStatus;
                int    tmpRepositionNode1, tmpRepositionNode2;
                int    repositionBothNodes;
                real8  vd1, vd2;
                Node_t *splitNode1, *splitNode2;
                Node_t *origNode, *newNode, *mergedNode;
                MobArgs_t mobArgs1, mobArgs2;

/*
 *              Attempt to split the node.  If the split fails,
 *              just move on and check the next possible way of
 *              splitting the node, but if the split succeeded
 *              we need to evaluate the results then restore
 *              the nodes to the original unsplit state.
 */
                for (armIndex = 0, k = numNbrs - 1; k >= 0; k--) {
                    if (armSets[setID][k] == 1) {
                        armList[armIndex++] = k;
                    }
                }

/*
 *              IMPORTANT! Splitting a ghost node rather than a local node
 *                         requires special handling, so if were dealing
 *                         with a ghost node, use a different 'split'
 *                         function.
 *
 *              Note: for this test, both the node being split
 *              and the resulting new node will intially retain the
 *              location and velocity of the node being split.
 *              This may lead to a zero length segment between the
 *              two nodes, but the mobility function should be okay
 *              with that.  Also, this operation MUST NOT be
 *              communicated to remote domains from processing!
 *
 *              Note: on success, the arms selected to be split
 *              will be connected to splitNode2 and all unselected
 *              arms will be attached to splitNode1.
 */
                if (isGhostMultiNode) {
                    splitStatus = SplitGhostNode(home, OPCLASS_SEPARATION,
                                                 node, origPos, origPos,
                                                 origVel, origVel,
                                                 armSets[setID][numNbrs],
                                                 armList, &splitNode1,
                                                 &splitNode2, 0);
                } else {
                    splitStatus = SplitNode(home, OPCLASS_SEPARATION,
                                            node, origPos, origPos,
                                            origVel, origVel,
                                            armSets[setID][numNbrs],
                                            armList, 0, &splitNode1,
                                            &splitNode2, 0);
                }

                if (splitStatus != SPLIT_SUCCESS) {
                    setID = totalSets;
                    nodeInfo[i].status = 1;
                    continue;
                }

/*
 *              For the initial split, the two nodes involved
 *              both remain at the same location as the orginal
 *              node.  This means that the forces on all segments
 *              attached to the nodes should be unchanged from
 *              the corresponding segments attached to the original
 *              node; this means we can avoid the expense of
 *              calculating forces here, and just copy the
 *              corresponding forces from the backup copies of
 *              the original nodes...
 */
                GetForcesFromBkup(home, splitNode1, bkupNodeList);
                GetForcesFromBkup(home, splitNode2, bkupNodeList);

#ifdef ESHELBY
/*
 *              If we're using special mobility for segments
 *              intersecting eshelby inclusions, we need to
 *              generate the segment/particle intersection list
 *              (and clear it after calculating velocities)
 */
                FindNodePartIntersects(home, splitNode1);
                FindNodePartIntersects(home, splitNode2);

                SegPartListSort(home);

#endif
                mobStatus  = EvaluateMobility(home, splitNode1, &mobArgs1);

                mobStatus  = EvaluateMobility(home, splitNode2, &mobArgs2);

#ifdef ESHELBY
                SegPartListClear(home);
#endif

                if (mobStatus != 0) {
                    nodeInfo[i].status = 2;
                    setID = totalSets;
                    continue;
                }


/*
 *              When the split creates a new junction segment, the
 *              segment is initially zero length but this segment
 *              will affect the force, velocity and direction of
 *              the two attached nodes.  So, we need to try to
 *              determine the optimal junction direction and adjust
 *              the force/velocity of the two nodes appropriately.
 */
                AdjustJuncNodeForceAndVel(home, splitNode1, &mobArgs1);

                vd1 = (splitNode1->vX * splitNode1->vX) +
                      (splitNode1->vY * splitNode1->vY) +
                      (splitNode1->vZ * splitNode1->vZ);

                AdjustJuncNodeForceAndVel(home, splitNode2, &mobArgs2);

                vd2 = (splitNode2->vX * splitNode2->vX) +
                      (splitNode2->vY * splitNode2->vY) +
                      (splitNode2->vZ * splitNode2->vZ);

/*
 *              This configuration should only be considered as
 *              a viable option if the velocity of the nodes
 *              after the split are not near zero.
 */
                if ((fabs(vd1) > eps) || (fabs(vd2) > eps)) {
                    int   numSegs;
                    real8 dirx, diry, dirz;
                    real8 invvNorm;
                    real8 distFactor;
                    real8 minSplitDist;
                    real8 vel1Dotvel2, vd1Sqrt, vd2Sqrt;
                    SegNodes_t *segList;
                    SMN_Force_t *node1Forces, *node2Forces;

                    vel1Dotvel2 = (splitNode1->vX * splitNode2->vX) +
                                  (splitNode1->vY * splitNode2->vY) +
                                  (splitNode1->vZ * splitNode2->vZ);

                    vd1Sqrt = sqrt(vd1);
                    vd2Sqrt = sqrt(vd2);

                    if (fabs(vel1Dotvel2/vd1Sqrt/vd2Sqrt+1.0) < 0.01) {
/*
 *                      In this case, the velocity of the two nodes
 *                      is nearly equal and opposite, so we'll treat
 *                      them as *exactly* equal and opposite and
 *                      later reposition both nodes.
 */
                        tmpRepositionNode1 = 1;
                        tmpRepositionNode2 = 1;

                        if (vd1 > vd2) {

                            splitNode2->vX = -splitNode1->vX;
                            splitNode2->vY = -splitNode1->vY;
                            splitNode2->vZ = -splitNode1->vZ;

                            vd2 = (splitNode2->vX * splitNode2->vX) +
                                  (splitNode2->vY * splitNode2->vY) +
                                  (splitNode2->vZ * splitNode2->vZ);

                        } else {

                            splitNode1->vX = -splitNode2->vX;
                            splitNode1->vY = -splitNode2->vY;
                            splitNode1->vZ = -splitNode2->vZ;

                            vd1 = (splitNode1->vX * splitNode1->vX) +
                                  (splitNode1->vY * splitNode1->vY) +
                                  (splitNode1->vZ * splitNode1->vZ);

                        }

                    }  /* end if (fabs(vel1Dotvel2 / ... ) < 0.01) */

                    if (vd1 > vd2) {

                        tmpRepositionNode1 = 1;
                        invvNorm = 1.0 / sqrt(vd1);

                        dirx = -splitNode1->vX * invvNorm;
                        diry = -splitNode1->vY * invvNorm;
                        dirz = -splitNode1->vZ * invvNorm;

                    } else {

                        tmpRepositionNode2 = 1;
                        invvNorm = 1.0 / sqrt(vd2);

                        dirx = splitNode2->vX * invvNorm;
                        diry = splitNode2->vY * invvNorm;
                        dirz = splitNode2->vZ * invvNorm;
                    }

/*
 *                  We initially want to separate the two nodes
 *                  by a distance of <splitDist>.  If we're
 *                  moving both nodes, we'll move them each half
 *                  that distance, if we're only moving one node
 *                  we move it the full distance.
 */
                    repositionBothNodes = tmpRepositionNode1 *
                                          tmpRepositionNode2;

                    if (repositionBothNodes) {
                        minSplitDist = splitDist * 0.5;
                    } else {
                        minSplitDist = splitDist;
                    }

                    distFactor = minSplitDist * tmpRepositionNode1;

                    splitNode1->x -= dirx * distFactor;
                    splitNode1->y -= diry * distFactor;
                    splitNode1->z -= dirz * distFactor;

                    FoldBox(param, &splitNode1->x, &splitNode1->y,
                            &splitNode1->z);

                    distFactor = minSplitDist * tmpRepositionNode2;

                    splitNode2->x += dirx * distFactor;
                    splitNode2->y += diry * distFactor;
                    splitNode2->z += dirz * distFactor;

                    FoldBox(param, &splitNode2->x, &splitNode2->y,
                            &splitNode2->z);

/*
 *                  Evaluate this domain's contribution to the forces on
 *                  the segments of nodes involved in the split and store
 *                  the forces for later use.
 */
                    node1Forces = &nodeInfo[i].splitInfo[setID].node1Force;
                    node2Forces = &nodeInfo[i].splitInfo[setID].node2Force;

/*
 *                  Build the list of locally owned segments and use this
 *                  list to calculate this domain's contribution to the
 *                  forces on the segments of nodes involved in the split,
 *                  and store the forces for later use.
 */
                    BuildSMNLocalSegList(home, &segList, &numSegs);

                    SMNPartialForce(home, splitNode1, numSegs, segList,
                                    node1Forces);
                    SMNPartialForce(home, splitNode2, numSegs, segList,
                                    node2Forces);

                    if (segList != (SegNodes_t *)NULL) {
                        free(segList);
                        segList = (SegNodes_t *)NULL;
                    }

                    nodeInfo[i].splitInfo[setID].node1Segs=splitNode1->numNbrs;
                    nodeInfo[i].splitInfo[setID].node2Segs=splitNode2->numNbrs;
                }

/*
 *              Merge the new node back into the original node then restore
 *              the original node and it's neighbors to the pristine unsplit
 *              state.  (Need to do a restore to insure the ordering of arms
 *              since the arrays of splitting possibilites use the
 *              original node's arm indices, plus we get back all the original
 *              force/velocity/etc values.)
 *
 *              Assumes that MergeNode() will merge <newNode> into <origNode>!
 */
                if (node == splitNode1) {
                    origNode = splitNode1;
                    newNode  = splitNode2;
                } else {
                    newNode  = splitNode1;
                    origNode = splitNode2;
                }

/*
 *              This is a kludge to handle pinned multinodes.  If the first
 *              node is pinned, move the second node to the original position,
 *              pin it with the same constraints as were on the original 
 *              node and unpin the first node so the merge will take place
 *              as expected.
 *
 *              NOTE: For now, a node that is pinned in ANY dimension is
 *                    treated here as if it is pinned in ALL dimensions.
 */
                if (HAS_ANY_OF_CONSTRAINTS(newNode->constraint, PINNED_NODE)) {
                    origNode->x = origPos[X];
                    origNode->y = origPos[Y];
                    origNode->z = origPos[Z];
                    ADD_CONSTRAINTS(origNode->constraint,
                                    newNode->constraint & PINNED_NODE);
                    REMOVE_CONSTRAINTS(newNode->constraint, PINNED_NODE);
                }


                if (isGhostMultiNode) {
                    MergeGhostNode(home, OPCLASS_SEPARATION, newNode, origNode,
                                   origPos, &mergedNode, &mergeStatus);
                } else {
                    MergeNode(home, OPCLASS_SEPARATION, newNode, origNode,
                              origPos, &mergedNode, &mergeStatus, 0);
                }

/*
 *              If the merge failed, something's screwy!  Abort.
 */
                if ((mergeStatus & MERGE_SUCCESS) == 0) {
                    Fatal("Unable to merge (%d,%d) into (%d,%d)"
                          "during separation test",
                          newNode->myTag.domainID, newNode->myTag.index,
                          origNode->myTag.domainID, origNode->myTag.index);
                }

/*
 *              Restore nodes to the original state
 */
                for (k = 0; k < bkupNodeCount; k++) {
                    Node_t *tmpNode;

                    tmpNode = GetNodeFromTag(home, bkupNodeList[k].myTag);
                    RestoreNode(home, tmpNode, &bkupNodeList[k]);
                }

            }  /* for (setID = 0; setID < totalSets; ...) */

            for (k = 0; k < totalSets; k++) {
                free(armSets[k]);
            }

            free(armSets);
            free(armList);

            armList = (int *)NULL;

/*
 *          Need to free the arrays and node structures used
 *          above for backup/restore or we end up with a memory leak.
 */
            for (k = 0; k < bkupNodeCount; k++) {
                FreeNodeArrays(&bkupNodeList[k]);
            }

            free(bkupNodeList);
            bkupNodeList = (Node_t *)NULL;

        }  /* loop over multinodes in list*/

        return;
}
