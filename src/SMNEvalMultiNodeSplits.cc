/*****************************************************************************
 *
 *      Module:       SMNEvalMultiNodeSplits.c
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
 *                    remote domains via CommSendRemesh()/FixRemesh()
 *                    calls in ParadisStep().
 *
 *      Included functions:
 *
 *****************************************************************************/
#include <stdio.h>

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
 *      Function:     SMNEvalMultiNodeSplits
 *
 *      Parameters:
 *          nodeCount  IN:  Number of local multi-nodes in <nodeInfo> array
 *          nodeInfo   IN:  Array of <nodeCount> structs containing
 *                          information about each local multi-node to
 *                          be evaluated when doing parallel multi-node
 *                          splitting.
 *
 *-------------------------------------------------------------------------*/
void SMNEvalMultiNodeSplits(Home_t *home, int nodeCount, SMN_Info_t *nodeInfo)
{
        int     listIndex;
        int     numNbrs;
        int     globalOp;
        int     bkupNodeCount;
        int     *armList, *armList2;
        real8   vNoise;
        real8   splitDist;
        real8   eps;
        Node_t  *bkupNodeList;
        Param_t *param;
        SegForce_t *SegForces;


#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom   = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom   = -1;
#endif
        param     = home->param;

        vNoise    = param->splitMultiNodeAlpha * param->rTol / param->deltaTT;
        splitDist = param->rann * 2.0;

        armList  = (int *)NULL;
        armList2 = (int *)NULL;
        SegForces  = (SegForce_t *)NULL;

        globalOp = 1;
        eps      = 1.0e-12;

/*
 *      Allocate a temporary array in which to store segment forces while
 *      evaluating possible node splits.  We'll allocate the array large
 *      enough to handle the split of a node with the maximum number of
 *      segments.
 */
        SegForces = (SegForce_t *)malloc((MAX_NBRS + 2) * sizeof(SegForce_t));

/*
 *      Loop through all native multinodes
 */
        for (listIndex = 0; listIndex < nodeCount; listIndex++) {
            int    k;
            int    origNumNbrs;
            int    segID, setID, setIndex;
            int    armIndex, armIndex2;
            int    segIndex;
            int    colliding;
            int    splitStatus, mergeStatus, mobStatus;
            int    totalSets;
            int    repositionNode1, repositionNode2;
            int    repositionBothNodes;
            int    **armSets;
            real8  powerMax, powerTest;
            real8  vd1, vd2;
            real8  distFactor, minSplitDist;
            real8  dirx, diry, dirz;
            real8  origPos[3], origVel[3];
            real8  vel1[3], vel2[3];
            Node_t *node, *nbrNode, *origNode, *newNode;
            Node_t *splitNode1, *splitNode2, *mergedNode;

            real8  splitPos1[3] = { 0.0, 0.0, 0.0 };
            real8  splitPos2[3] = { 0.0, 0.0, 0.0 };

            node = nodeInfo[listIndex].multiNode;

/*
 *          If any of the remote domains encountered a problem doing the
 *          multinpde split, we won't have all the info needed to determine
 *          whether to split the node or not, so just skip this multinode
 */
            if (nodeInfo[listIndex].status != 0) {
                continue;
            }

/*
 *          If we haven't already done so, allocate a temporary array
 *          in which to store segment forces while evaluating possible
 *          node splits.  We'll allocate the array large enough to
 *          handle the split of a 15-node.  If we need an array larger
 *          than that, we've probably got other problems.
 */
            if (SegForces == (SegForce_t *)NULL) {
                SegForces = (SegForce_t *)malloc(17 * sizeof(SegForce_t));
            }

/*
 *          Make a backup copy of the multinode and each of its neighbors
 */
            numNbrs = node->numNbrs;
            origNumNbrs = numNbrs;

            bkupNodeCount = node->numNbrs + 1;
            bkupNodeList = (Node_t *)malloc(bkupNodeCount * sizeof(Node_t));

            BackupNode(home, node, &bkupNodeList[0]);

            for (segID = 0; segID < node->numNbrs; segID++) {
                nbrNode = GetNodeFromTag(home, node->nbrTag[segID]);
                BackupNode(home, nbrNode, &bkupNodeList[segID+1]);
            }

#ifdef ESHELBY
            SegPartListClear(home);
#endif

/*
 *          The determination of which (if any) arm-splitting
 *          possibilities to use when breaking apart a node is
 *          based on the maximum energy release rate.  Use the
 *          original unsplit node as the baseline condition
 *          for the energy comparison.
 */
            powerMax = ((node->fX * node->vX) +
                        (node->fY * node->vY) +
                        (node->fZ * node->vZ));

/*
 *          Get the total number of ways in which this node can
 *          be split and build arrays of flags we'll use later
 *          for selecting which segments should be split off the
 *          node wich each corresponding split.
 */
            SMNGetArmSets(numNbrs, &totalSets, &armSets);

/*
 *          Make a copy of the original node position and velocity for later
 */
            colliding = 0;
            setIndex = -1;

            origPos[X] = node->x;
            origPos[Y] = node->y;
            origPos[Z] = node->z;

            origVel[X] = node->vX;
            origVel[Y] = node->vY;
            origVel[Z] = node->vZ;

            armList = (int *)calloc(1, (numNbrs - 2) * sizeof(int));
            armList2 = (int *)calloc(1, (numNbrs - 2) * sizeof(int));

            for (setID = 0; setID < totalSets; setID++) {
                int   tmpSegCount1, tmpSegCount2;
                int   tmpRepositionNode1, tmpRepositionNode2;
                real8 vDiffx, vDiffy, vDiffz;
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
                splitStatus = SplitNode(home, OPCLASS_SEPARATION,
                                        node, origPos, origPos,
                                        origVel, origVel,
                                        armSets[setID][numNbrs],
                                        armList, 0, &splitNode1,
                                        &splitNode2, 0);

                if (splitStatus != SPLIT_SUCCESS) {
                    setID = totalSets;
                    nodeInfo[listIndex].status = 1;
                    continue;
                }

/*
 *              Make sure the node segment counts match...
 *
 *              Just a temporary sanity check we can remove once the code
 *              is debugged and working correctly
 */
                tmpSegCount1 = nodeInfo[listIndex].splitInfo[setID].node1Segs;
                tmpSegCount2 = nodeInfo[listIndex].splitInfo[setID].node2Segs;

                if ((tmpSegCount1 != 0) &&
                    (tmpSegCount1 != splitNode1->numNbrs)) {
                    Fatal("Mismatch on node1 segment count");
                }

                if ((tmpSegCount2 != 0) &&
                    (tmpSegCount2 != splitNode2->numNbrs)) {
                    Fatal("Mismatch on node2 segment count");
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

                mobStatus |= EvaluateMobility(home, splitNode2, &mobArgs2);

#ifdef ESHELBY
                SegPartListClear(home);
#endif

/*
 *              After the split is done, reposition the nodes and then
 *              we can update the forces with the forces generated
 *              earlier.
 *
 *              Note: If the direction of velocity changes,
 *              don't do the split
 */
                if (mobStatus == 0) {
/*
 *                  When the split creates a new junction segment, the
 *                  segment is initially zero length but this segment
 *                  will affect the force, velocity and direction of
 *                  the two attached nodes.  So, we need to try to
 *                  determine the optimal junction direction and adjust
 *                  the force/velocity of the two nodes appropriately.
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
 *                  This configuration should only be considered as
 *                  a viable option if the velocity of the nodes
 *                  after the split are not near zero.
 */
                    if ((fabs(vd1) > eps) || (fabs(vd2) > eps)) {
                        real8 vel1Dotvel2;
                        real8 vd1Sqrt, vd2Sqrt;
                        real8 invvNorm;
                        SMN_Force_t *node1Forces, *node2Forces;

                        vel1Dotvel2 = (splitNode1->vX * splitNode2->vX) +
                                      (splitNode1->vY * splitNode2->vY) +
                                      (splitNode1->vZ * splitNode2->vZ);

                        vd1Sqrt = sqrt(vd1);
                        vd2Sqrt = sqrt(vd2);

                        if (fabs(vel1Dotvel2/vd1Sqrt/vd2Sqrt+1.0) < 0.01) {
/*
 *                          In this case, the velocity of the two nodes
 *                          is nearly equal and opposite, so we'll treat
 *                          them as *exactly* equal and opposite and
 *                          later reposition both nodes.
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
 *                      We initially want to separate the two nodes
 *                      by a distance of <splitDist>.  If we're
 *                      moving both nodes, we'll move them each half
 *                      that distance, if we're only moving one node
 *                      we move it the full distance.
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
 *                      The force components from all domains were previously
 *                      summed, so copy them into the nodes' force arrays
 */
                        node1Forces = &nodeInfo[listIndex].splitInfo[setID].node1Force;
                        node2Forces = &nodeInfo[listIndex].splitInfo[setID].node2Force;

                        CopyForcesToNode(home, splitNode1, node1Forces);
                        CopyForcesToNode(home, splitNode2, node2Forces);

/*
 *                      When the original node was split above, both new
 *                      nodes were left at the original position.  If a
 *                      new segment was created connecting the nodes, it
 *                      intially zero length and therefore had no glide-
 *                      plane, but the mobility fuynction can handle zero-
 *                      length segs.  After calculating velocities, though
 *                      we've shifted the nodes making the connecting
 *                      segment (if any) non-zero length, so now we have
 *                      to explicitly set the glide plane of that seg
 *                      because some of the mobility functions depend on
 *                      the glide plane being set for non-zero length segs
 *                          
 *                      Added by Wei Cai, 2009/09/07
 */
                        armIndex = GetArmID(splitNode1, splitNode2);

                        if (armIndex >= 0 ) {
                            real8 burg[3], dirvec[3], nplane[3];

/*
 *                          Calculate glide plane normal as
 *                          n = cross ( burg(x,y,z), dir(x,y,z))
 */
                            burg[X] = splitNode1->burgX[armIndex];
                            burg[Y] = splitNode1->burgY[armIndex];
                            burg[Z] = splitNode1->burgZ[armIndex];

                            dirvec[X] = dirx;
                            dirvec[Y] = diry;
                            dirvec[Z] = dirz;

                            FindPreciseGlidePlane(home, burg, dirvec,
                                                  nplane, home->param->allowFuzzyGlidePlanes);

/*
 *                          If the new segment is screw the glide
 *                          plane is undefined.  If that's the case, we
 *                          need to randomly select a glide plane
 *                          appropriate to the burgers vector.
 */
                            if (DotProduct(nplane, nplane) < 1.0e-03) {
                                PickScrewGlidePlane(home, burg, nplane);
                            }

                            Normalize(&nplane[X], &nplane[Y], &nplane[Z]);

                            splitNode1->nx[armIndex] = nplane[X];
                            splitNode1->ny[armIndex] = nplane[Y];
                            splitNode1->nz[armIndex] = nplane[Z];

                            armIndex = GetArmID(splitNode2, splitNode1);

                            splitNode2->nx[armIndex] = nplane[X];
                            splitNode2->ny[armIndex] = nplane[Y];
                            splitNode2->nz[armIndex] = nplane[Z];
                        }

#ifdef ESHELBY
/*
 *                      If we're using special mobility for segments
 *                      intersecting eshelby inclusions, we need to
 *                      generate the segment/particle intersection list
 *                      (and clear it after calculating velocities)
 */
                        FindNodePartIntersects(home, splitNode1);
                        FindNodePartIntersects(home, splitNode2);

                        SegPartListSort(home);

#endif

                        mobStatus  = EvaluateMobility(home, splitNode1,
                                                      &mobArgs1);

                        mobStatus |= EvaluateMobility(home, splitNode2,
                                                      &mobArgs2);

#ifdef ESHELBY
                        SegPartListClear(home);
#endif

                        vDiffx = splitNode2->vX - splitNode1->vX;
                        vDiffy = splitNode2->vY - splitNode1->vY;
                        vDiffz = splitNode2->vZ - splitNode1->vZ;

                    } else {
/*
 *                      Both nodes have near-zero velocity after the split, 
 *                      so don't choose this split
 */
                        vDiffx = -1.0;
                        vDiffy = -1.0;
                        vDiffz = -1.0;

                        dirx = 1.0;
                        diry = 1.0;
                        dirz = 1.0;
                    }

                }  /* if (mobStatus == 0) */

                if (mobStatus != 0) {
/*
 *                  If the mobility function was unable to converge on
 *                  velocities for either node in their new positions,
 *                  we can't really evaluate the node split properly.
 *                  so skip evaluation of any other possible splits for this
 *                  node this cycle.
 */
                    setIndex = -1;
                    setID = totalSets;

                } else if (((vDiffx*dirx)+(vDiffy*diry)+(vDiffz*dirz)) > 0) {

                    powerTest = ((splitNode1->fX * splitNode1->vX) +
                                 (splitNode1->fY * splitNode1->vY) +
                                 (splitNode1->fZ * splitNode1->vZ)) +
                                ((splitNode2->fX * splitNode2->vX) +
                                 (splitNode2->fY * splitNode2->vY) +
                                 (splitNode2->fZ * splitNode2->vZ)) -
                                vNoise * (
                                    sqrt(splitNode1->fX * splitNode1->fX +
                                         splitNode1->fY * splitNode1->fY +
                                         splitNode1->fZ * splitNode1->fZ)+
                                    sqrt(splitNode2->fX * splitNode2->fX +
                                         splitNode2->fY * splitNode2->fY +
                                         splitNode2->fZ * splitNode2->fZ));

/*
 *                  If this potential node split would result in the
 *                  highest energey release, save enough info to 
 *                  perform this split later.
 */
                    if ((powerTest - powerMax) > eps) {

                        vel1[X] = splitNode1->vX;
                        vel1[Y] = splitNode1->vY;
                        vel1[Z] = splitNode1->vZ;

                        vel2[X] = splitNode2->vX;
                        vel2[Y] = splitNode2->vY;
                        vel2[Z] = splitNode2->vZ;

                        splitPos1[X] = splitNode1->x;
                        splitPos1[Y] = splitNode1->y;
                        splitPos1[Z] = splitNode1->z;

                        splitPos2[X] = splitNode2->x;
                        splitPos2[Y] = splitNode2->y;
                        splitPos2[Z] = splitNode2->z;

                        repositionNode1 = tmpRepositionNode1;
                        repositionNode2 = tmpRepositionNode2;

                        powerMax = powerTest;
                        setIndex = setID;

/*
 *                      Save the forces on all segments of the two nodes
 *                      resulting from this split.  If this configuration is
 *                      the one chosen for the final split, we can then use
 *                      these preserved force values rather than go to the
 *                      expense of recomputing full forces on the nodes yet
 *                      once again.
 */
                        segIndex = 0;
                        SaveSegForces(home, splitNode1, SegForces, &segIndex);
                        SaveSegForces(home, splitNode2, SegForces, &segIndex);

/*
 *                      For 4-nodes only check if the split would likely
 *                      result in an immediate collision reversing the effects
 *                      of the split.  If this is the case, we'll skip the
 *                      split to avoid flicker that might affect the timestep.
 *                      We don't want to form many-armed nodes, though, so if
 *                      we've got more than 4 arms, don't use this criteria to
 *                      restrict a node split.
 */       
                        if (origNumNbrs == 4) {
                            colliding = CheckCollisionConditions(home,
                                                                 splitNode2,
                                                                 splitNode1);
                        } else {
                            colliding = 0;
                        }
                    }
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
 *              pin it in place and unpin the first node so the merge will
 *              take place as expected.
 *
 *              NOTE: For now, a node that is pinned in ANY dimension is
 *                    treated here as if it is pinned in ALL dimensions
 */
                if (newNode->constraint == PINNED_NODE) {
                    origNode->x = origPos[X];
                    origNode->y = origPos[Y];
                    origNode->z = origPos[Z];
                    REMOVE_CONSTRAINTS(newNode->constraint, PINNED_NODE);
                    ADD_CONSTRAINTS(origNode->constraint, PINNED_NODE);
                }

                MergeNode(home, OPCLASS_SEPARATION, newNode, origNode,
                          origPos, &mergedNode, &mergeStatus, 0);

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

/*
 *          If it is necessary (and possible) to split the node,
 *          armSet[setIndex] points to the array indicating which
 *          arms are to be split off and which are to remain with
 *          the original node.
 *
 *          However, if configuration resulting after the split
 *          would result in an imminent collision which might
 *          possibly restore the original configuration, we won't
 *          do the split.
 */
            if (colliding) {
                setIndex = -1;
            }
                
            if (setIndex >= 0) {
                int   q;
                int   tmpArm;
                int   armSets2;
                real8 normsplit;
                real8 invnormsplit, dsPx, dsPy, dsPz;
                real8 junctionLen1, junctionLen2;
                real8 pos1[3], pos2[3];

                junctionLen1 = 0.0;
                junctionLen2 = 0.0;

/*
 *              Set up the list of arms to be passed to SplitNode.
 *              Note: the arm list MUST be in reverse order (i.e.
 *              highest arm index to lowest) These operations are
 *              always considered *global* so other domains will be
 *              notified of these ops.
 *
 *              In some cases we also need the list (armList2) of
 *              arms NOT being split off, so create that here as well.
 */
                armIndex = 0;
                armIndex2 = 0;
                armSets2 = numNbrs - armSets[setIndex][numNbrs];

                for (k = numNbrs - 1; k >= 0; k--) {
                    if (armSets[setIndex][k] == 1) {
                        armList[armIndex++] = k;
                    } else {
                        armList2[armIndex2++] = k;
                    }
                }

/*
 *              Get direction from saved data
 *              By definition, dir is (x2-x1)/|x2-x1|
 */
                dsPx = splitPos2[X] - splitPos1[X];
                dsPy = splitPos2[Y] - splitPos1[Y];
                dsPz = splitPos2[Z] - splitPos1[Z];

                ZImage(param, &dsPx, &dsPy, &dsPz);

                normsplit = dsPx*dsPx + dsPy*dsPy + dsPz*dsPz;
                invnormsplit = 1.0 / sqrt(normsplit);

                dirx = dsPx * invnormsplit;
                diry = dsPy * invnormsplit;
                dirz = dsPz * invnormsplit;

/*
 *              If we are restricted to creating junctions of a
 *              uniform length, move 1 node the uniform distance or
 *              both nodes half that distance.  Otherwise, use a line
 *              tension approximation to get the optimal distance to
 *              move the nodes along the junction direction.
 *
 *              By convention, splitNode2 has the armList, here pos2 has
 *              the arms list. pos1 has the other arms.
 */
                repositionBothNodes = repositionNode1 * repositionNode2;

                if (repositionBothNodes) {
                    minSplitDist = splitDist * 0.5;
                } else {
                    minSplitDist = splitDist;
                }

                distFactor = 0;

                if (repositionNode1) {

                    if (param->useUniformJunctionLen) {

                        distFactor = minSplitDist;

                    } else {

                        junctionLen1 =
                                CalcJunctionLen(home, node, armSets2, armList2,
                                                -dirx, -diry, -dirz,
                                                minSplitDist,
                                                repositionBothNodes);

                        distFactor = MAX(minSplitDist, junctionLen1);
                    }
                }

                pos1[X] = node->x - dirx * distFactor;
                pos1[Y] = node->y - diry * distFactor;
                pos1[Z] = node->z - dirz * distFactor;

                FoldBox(param, &pos1[X], &pos1[Y], &pos1[Z]);

                distFactor = 0;

                if (repositionNode2) {

                    if (param->useUniformJunctionLen) {

                        distFactor = minSplitDist;

                    } else {
                        junctionLen2 =
                                CalcJunctionLen(home, node,
                                                armSets[setIndex][numNbrs],
                                                armList, dirx, diry, dirz,
                                                minSplitDist,
                                                repositionBothNodes);

                        distFactor = MAX(minSplitDist, junctionLen2);
                    }
                }

                pos2[X] = node->x + dirx * distFactor;
                pos2[Y] = node->y + diry * distFactor;
                pos2[Z] = node->z + dirz * distFactor;
                

                FoldBox(param, &pos2[X], &pos2[Y], &pos2[Z]);

/*
 *             Note: after the split, the arms selected to be split
 *             will be connected to splitNode2 which will be
 *             positioned at pos2 with velocity vel2.  Node
 *             splitNode1 will be attached to the arms not selected
 *             to be moved, and will have the position and velocity
 *             indicated by pos1 and vel1.
 */
               splitStatus = SplitNode(home, OPCLASS_SEPARATION,
                                       node, pos1, pos2, vel1, vel2,
                                       armSets[setIndex][numNbrs],
                                       armList, globalOp,
                                       &splitNode1, &splitNode2, 0);
/*
 *             If the split failed, clean up then go back and
 *             continue looking for more multinodes to split.
 */
               if (splitStatus != SPLIT_SUCCESS) {
                   free(armList);
                   armList = (int *)NULL;
                   for (k = 0; k < totalSets; k++) free(armSets[k]);
                   free(armSets);
                   free(armList2);
                   continue;
               }

#ifdef DEBUG_LOG_MULTI_NODE_SPLITS
               splitNode1->multiNodeLife = 0;
               splitNode2->multiNodeLife = 0;
#endif

/*
 *             The forces for all segments attached to the two nodes were
 *             preserved during the evaluation step above, use those
 *             preserved values to reset the forces on these nodes rather
 *             than explicitly recompute the forces.
 *
 *             WARNING! Make sure GetSavedSegForces() is called
 *             for the two nodes in the same order SaveSegForces()
 *             was called for the two nodes.
 */
               segIndex = 0;
               GetSavedSegForces(home, splitNode1, SegForces, &segIndex);
               GetSavedSegForces(home, splitNode2, SegForces, &segIndex);

/*
 *             Nodal velocities passed into SplitNode() should be okay,
 *             but we do need to call ResetSegForces2() to distribute
 *             the new forces out to the remote domains.
 */
               if (node == splitNode1) {
                   origNode = splitNode1;
                   newNode  = splitNode2;
               } else {
                   newNode  = splitNode1;
                   origNode = splitNode2;
               }

               for (q = 0; q < origNode->numNbrs; q++) {
                   nbrNode = GetNodeFromTag(home, origNode->nbrTag[q]);
                   tmpArm = GetArmID(nbrNode, origNode);
                   ResetSegForces2(home, origNode, &nbrNode->myTag,
                                   origNode->armfx[q],
                                   origNode->armfy[q],
                                   origNode->armfz[q],
                                   nbrNode->armfx[tmpArm],
                                   nbrNode->armfy[tmpArm],
                                   nbrNode->armfz[tmpArm], globalOp);
               }

               for (q = 0; q < newNode->numNbrs; q++) {
                   nbrNode = GetNodeFromTag(home, newNode->nbrTag[q]);
                   tmpArm = GetArmID(nbrNode, newNode);
                   ResetSegForces2(home, newNode, &nbrNode->myTag,
                                   newNode->armfx[q],
                                   newNode->armfy[q],
                                   newNode->armfz[q],
                                   nbrNode->armfx[tmpArm],
                                   nbrNode->armfy[tmpArm],
                                   nbrNode->armfz[tmpArm], globalOp);
               }

/*
 *             Mark both nodes involved in the split as 'exempt'
 *             from subsequent collisions this time step. 
 */
               origNode->flags |= NO_COLLISIONS;
               newNode->flags  |= NO_COLLISIONS;

/*
 *             When debugging, dump some info on topological
 *             changes taking place and the nodes involved
 */
#ifdef DEBUG_TOPOLOGY_CHANGES
                if ((dbgDom < 0)||(dbgDom == home->myDomain)) {
                    printf("  Separate: (%d,%d) from (%d,%d)\n",
                           newNode->myTag.domainID, newNode->myTag.index,
                           origNode->myTag.domainID, origNode->myTag.index);
                    PrintNode(origNode);
                    PrintNode(newNode);
                }
#endif
            }  /* if (setIndex >= 0) */

            for (k = 0; k < totalSets; k++) {
                free(armSets[k]);
            }

            free(armSets);
            free(armList);
            free(armList2);

            armList = (int *)NULL;
            armList2 = (int *)NULL;

/*
 *          Need to free the arrays and node structures used
 *          above for backup/restore or we end up with a memory leak.
 */
            for (k = 0; k < bkupNodeCount; k++) {
                FreeNodeArrays(&bkupNodeList[k]);
            }

            free(bkupNodeList);
            bkupNodeList = (Node_t *)NULL;

        }  /* loop over local multinodes */

        if (SegForces != (SegForce_t *)NULL) {
            free(SegForces);
        }

        return;
}
