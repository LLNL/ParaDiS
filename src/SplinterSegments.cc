/*-------------------------------------------------------------------------
 *
 *      Module:       SplinterSegments.c
 *
 *      Description:  This module contains functions specific to the
 *                    process of 'splintering' a segment into two
 *                    segments of a different burgers vector.
 *
 *                    The splintering process actually begins with a
 *                    node with exactly two segments, each with the
 *                    same glide plane and burgers vector where the
 *                    burgers vector is one that may be decomposed into
 *                    a pair of identical burgers vectors (that sum
 *                    to the original burgers vector) and are mutually
 *                    repulsive, resulting in a lower energy configuration
 *                    than the original configuration.
 *
 *                    The burgers vectors of the original node's segments
 *                    are converted to the new burgers vector, a new
 *                    node is cerated and given duplicates of the original
 *                    node's segments, and both nodes are repositioned
 *                    slightly.
 *                    
 *      Includes public functions
 *          SplinterSegments()
 *
 *      Includes private functions
 *          IsSplinterableBurg()
 *
 *-----------------------------------------------------------------------*/
#include <stdio.h>

#include "mpi_portability.h"

#include "Home.h"

static int dbgDom;


/*-------------------------------------------------------------------------
 *
 *      Function:    IsSplinterableBurg
 *
 *      Description: Determine if segments of the specified burgers vector 
 *                   may be decomposed (or 'splintered') into two
 *                   mutually replusive segments of a different burgers
 *                   vectors, and if so, return to the caller the burgers
 *                   vector into which the segments it can be 'splintered'
 *
 *      Parameters:
 *          IN:  burgOrig  Burgers vector to look up in the list of
 *                         splinterable burgers vectors.  This burgers
 *                         vector is specified in the laboratory frame
 *                         (if different from the crystalographic frame).
 *
 *          OUT: burgNew   Burgers vector into which <burgOrig> may be
 *                         splintered.  This burgers vector will be
 *                         specified in the laboratory frame (if different
 *                         from the crystalographic frame).  If <origBurg>
 *                         is not splinterable, this vector is zeroed.
 *
 *      Returns:  0 if segments with this burgers vector may not splinter
 *                1 if they are permitted to do so.
 *
 *------------------------------------------------------------------------*/
static int IsSplinterableBurg(Home_t *home, real8 burgOrig[3], real8 burgNew[3])
{
        int      index;
        int      useLabFrame;
        int      doSplinter;
        Param_t *param = home->param;
        real8    bOrigCrystal[3], bNewCrystal[3];


        VECTOR_ZERO(burgNew);

/*
 *      For now, this stuff only applies to HCP materials
 */
        switch(home->param->materialType) {
            case MAT_TYPE_HCP:
                break;

            default:
                return(0);
        }

/*
 *      The burgers vector's on the list of burgers vectors that can
 *      splinter are in the crystalographic frame, so if necessary
 *      rotate the incoming burgers vector from the lab to crystal frame.
 */
        useLabFrame = param->useLabFrame;

        if (useLabFrame) {
            Matrix33Vector3Multiply(param->rotMatrixInverse,
                                    burgOrig, bOrigCrystal);
        } else {
            VECTOR_COPY(bOrigCrystal, burgOrig);
        }

/*
 *      Check if the provided burgers vector is one of the splinterable ones.
 */
        doSplinter = 0;

        for (index = 0; index < home->burgData.numSplinterableBurgs; index++) {
            real8 burgSign, refBurgDotOrig;
            real8 *newBurg;
            real8 *refBurg = home->burgData.splinterableBurgList[index].refBurg;

            refBurgDotOrig = DotProduct(refBurg, bOrigCrystal);

            burgSign = ((refBurgDotOrig < 0.0) ? -1.0 : 1.0);

/*
 *          If the burgers vector may be splintered, save the burgers vector
 *          into which it may be spintered for the caller.  Note: the provided
 *          burgers vector may match the reference burgers vector but be
 *          of the opposite sign.  If so, we have to reverse the sign on
 *          the returned burgers vector.
 */
            if ((fabs(Normal(refBurg) - Normal(bOrigCrystal)) < 1e-4) &&
                ((fabs(refBurg[X] - burgSign * bOrigCrystal[X]) < 1.e-4) &&
                 (fabs(refBurg[Y] - burgSign * bOrigCrystal[Y]) < 1.e-4) &&
                 (fabs(refBurg[Z] - burgSign * bOrigCrystal[Z]) < 1.e-4))) {

                newBurg = home->burgData.splinterableBurgList[index].newBurg;

                bNewCrystal[X] = newBurg[X] * burgSign;
                bNewCrystal[Y] = newBurg[Y] * burgSign;
                bNewCrystal[Z] = newBurg[Z] * burgSign;

                doSplinter = 1;
                break;
            }
        }

/*
 *      If needed, rotate the new burgers vector back to the laboratory
 *      frame for the caller.
 */
        if (doSplinter) {
            if (useLabFrame) {
                Matrix33Vector3Multiply(param->rotMatrix, bNewCrystal, burgNew);
            } else {
                VECTOR_COPY(burgNew, bNewCrystal);
            }
        }

        return(doSplinter);
}


/*-------------------------------------------------------------------------
 *
 *      Function:    SplinterSegments
 *
 *      Description: Main funtion that loops through all the nodes
 *                   looking for nodes with segments that should
 *                   be splintered.
 *
 *------------------------------------------------------------------------*/
void SplinterSegments(Home_t *home)
{
        int     globalOp;
        int     nodeIndex;
        real8   splitDist, eps;
        Param_t *param;


#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom   = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom   = -1;
#endif

        param = home->param;

/*
 *      If burgers vector splitting has not been enabled via the
 *      control file parameter, or it is not supposed to be done during
 *      the current cycle, no need to do anything.
 */
        if (param->splinterSegFrequency == 0) {
            return;
        }

        if ((home->cycle % param->splinterSegFrequency) != 0) {
            return;
        }

/*
 *      If the current material type does not produce segments that
 *      must be splintered into 2 segments with a different burgers
 *      vector, just return
 */
        switch(param->materialType) {
            case MAT_TYPE_HCP:
                break;

            default:
                return;
        }

        globalOp = 1;
        eps = 1.0e-6;
        splitDist = param->rann * 2.0 + eps;

/*
 *      Loop through all the local nodes looking for any discretization
 *      nodes whose segments' burgers vectors match a burgers vector
 *      that must be splintered.
 */
        for (nodeIndex = 0; nodeIndex < home->newNodeKeyPtr; nodeIndex++) {
            real8  burgSeg1[3], burgSeg2[3], burgCross[3], burgNew[3];
            real8  planeSeg1[3], planeSeg2[3], planeCross[3];
            real8  nodePos[3], nbrPos[3];
            real8  lineDir[3], splitDir[3];
            real8  newPos1[3], newPos2[3], newVel[3];
            Node_t *node, *neighborNode1, *neighborNode2, *newNode;

            node = home->nodeKeys[nodeIndex];

            if (node == (Node_t *)NULL) {
                continue;
            }

/*
 *          For now, only unconstrained 2-nodes are even considered
 */
            if ((node->numNbrs != 2) ||
                HAS_ANY_OF_CONSTRAINTS(node->constraint, CONSTRAINT_ANY))
            {
                continue;
            }

            neighborNode1 = GetNeighborNode(home, node, 0);
            neighborNode2 = GetNeighborNode(home, node, 1);

/*
 *          If the node's segments have different burgers vectors
 *          or glide planes, the node is not a pure discretization
 *          node, so skip it.
 */
            burgSeg1[X] = node->burgX[0];
            burgSeg1[Y] = node->burgY[0];
            burgSeg1[Z] = node->burgZ[0];

            burgSeg2[X] = node->burgX[1];
            burgSeg2[Y] = node->burgY[1];
            burgSeg2[Z] = node->burgZ[1];

            cross(burgSeg1, burgSeg2, burgCross);

            if (fabs(DotProduct(burgCross, burgCross)) > 1.0e-6) {
                continue;
            }

            planeSeg1[X] = node->nx[0];
            planeSeg1[Y] = node->ny[0];
            planeSeg1[Z] = node->nz[0];

            planeSeg2[X] = node->nx[1];
            planeSeg2[Y] = node->ny[1];
            planeSeg2[Z] = node->nz[1];

            cross(planeSeg1, planeSeg2, planeCross);

            if (fabs(DotProduct(planeCross, planeCross)) > 1.0e-6) {
                continue;
            }

/*
 *          Based on the burgers vector, check if the segment should
 *          be split into two, and if so, what the new burgers vector
 *          should be.
 */
            if (!IsSplinterableBurg(home, burgSeg1, burgNew)) {
                continue;
            }

/*
 *          In order to perform splinter operations in parallel
 *          in different domains, we have to abide by node/segment
 *          ownership rules...
 *
 *          If either of the node's neighbor's is in a remote domain,
 *          the splinter operation can only be done if both of the
 *          segments are owned by the local node.  This restriction
 *          means that the neighboring nodes have segments they
 *          do not own, which guarantees the neighboring nodes
 *          can not be deleted by the remote domain.  This guarantee
 *          is needed because this process will update the neighbor
 *          nodes, but since the only changes are to segments owned
 *          by this node or to add new segments (no deletions of any
 *          sort) we should be able to get away with this.
 *
 *          If all three nodes are owned by the same domain, we
 *          can go ahead and do the splinter operation regardless of
 *          which nodes own the segments.
 *
 *          WARNING:  Any attempt to thread this function MUST prevent
 *                    multiple threads within the same domain from
 *                    acting on overlapping groups of nodes!
 *                    threads 
 */
            if ((neighborNode1->myTag.domainID != node->myTag.domainID) ||
                (neighborNode2->myTag.domainID != node->myTag.domainID)) {

                if ((NodeOwnsSeg(home, node, neighborNode1) == 0) ||
                    (NodeOwnsSeg(home, node, neighborNode2) == 0)) {
                    continue;
                }
            }

/*
 *          First we need to change the burgers vector of the original
 *          two segments.  Note: the ChangeArmBurg() call only alters
 *          the segment burgers vector in the first specified node.
 *          It must be called twice for each segment to make the
 *          necessary changes in all affected nodes.
 */
            ChangeArmBurg(home, node, &neighborNode1->myTag,
                          burgNew[X], burgNew[Y], burgNew[Z],
                          planeSeg1[X], planeSeg1[Y], planeSeg1[Z],
                          globalOp, DEL_SEG_NONE);

            ChangeArmBurg(home, neighborNode1, &node->myTag,
                          -burgNew[X], -burgNew[Y], -burgNew[Z],
                          planeSeg1[X], planeSeg1[Y], planeSeg1[Z],
                          globalOp, DEL_SEG_NONE);

            ChangeArmBurg(home, node, &neighborNode2->myTag,
                          -burgNew[X], -burgNew[Y], -burgNew[Z],
                          planeSeg2[X], planeSeg2[Y], planeSeg2[Z],
                          globalOp, DEL_SEG_NONE);

            ChangeArmBurg(home, neighborNode2, &node->myTag,
                          burgNew[X], burgNew[Y], burgNew[Z],
                          planeSeg2[X], planeSeg2[Y], planeSeg2[Z],
                          globalOp, DEL_SEG_NONE);

/*
 *          Now calculate the locations at which the two resulting
 *          nodes will be positioned.
 */
            nodePos[X] = node->x;
            nodePos[Y] = node->y;
            nodePos[Z] = node->z;

            nbrPos[X] = neighborNode1->x;
            nbrPos[Y] = neighborNode1->y;
            nbrPos[Z] = neighborNode1->z;

            PBCPOSITION(param, nodePos[X], nodePos[Y], nodePos[Z],
                        &nbrPos[X], &nbrPos[Y], &nbrPos[Z]);

            lineDir[X] = nbrPos[X] - nodePos[X];
            lineDir[Y] = nbrPos[Y] - nodePos[Y];
            lineDir[Z] = nbrPos[Z] - nodePos[Z];

            cross(planeSeg1, lineDir, splitDir);
            NormalizeVec(splitDir);

/*
 *          New position for original node
 */
            newPos1[X] = nodePos[X] + (splitDir[X] * splitDist);
            newPos1[Y] = nodePos[Y] + (splitDir[Y] * splitDist);
            newPos1[Z] = nodePos[Z] + (splitDir[Z] * splitDist);

            FoldBox(param, &newPos1[X], &newPos1[Y], &newPos1[Z]);

            RepositionNode(home, newPos1, &node->myTag, globalOp);

/*
 *          newPos2 is where the newly created node will be placed
 */
            newPos2[X] = nodePos[X] - (splitDir[X] * splitDist);
            newPos2[Y] = nodePos[Y] - (splitDir[Y] * splitDist);
            newPos2[Z] = nodePos[Z] - (splitDir[Z] * splitDist);

            FoldBox(param, &newPos2[X], &newPos2[Y], &newPos2[Z]);

/*
 *          Allocate a new node, place it at the desired location
 */
            newNode = GetNewNativeNode(home);
            FreeNodeArms(newNode);

            newNode->native = 1;
            SET_CONSTRAINTS(newNode->constraint, UNCONSTRAINED);

            newNode->x = newPos2[X];
            newNode->y = newPos2[Y];
            newNode->z = newPos2[Z];

            newNode->oldvX = 0.0;
            newNode->oldvY = 0.0;
            newNode->oldvZ = 0.0;

            newNode->vX = 0.0;
            newNode->vY = 0.0;
            newNode->vZ = 0.0;

            VECTOR_ZERO(newVel);

            AssignNodeToCell(home, newNode);

/*
 *          If this is a 'global' topological change that remote domains
 *          need to know about, add the operations to the list that will
 *          be sent to the neighboring domains.
 */
            if (globalOp) {

                AddOpSplitNode(home, &node->myTag, &newNode->myTag,
                               newPos2, newVel);

            }  /* end if (globalOp) */

/*
 *          Insert segments between the new node and the two
 *          neighbors of the original node.  As with the ChangeArmBurg()
 *          calls above, this only adds the segment in the first specified
 *          node.  It must be called twice for each segment to make the
 *          necessary changes in all affected nodes.
 */
            InsertArm(home, newNode, &neighborNode1->myTag,
                      burgNew[X], burgNew[Y], burgNew[Z],
                      planeSeg1[X], planeSeg1[Y], planeSeg1[Z], globalOp);
                      
            InsertArm(home, neighborNode1, &newNode->myTag, 
                      -burgNew[X], -burgNew[Y], -burgNew[Z],
                      planeSeg1[X], planeSeg1[Y], planeSeg1[Z], globalOp);
                      
            InsertArm(home, newNode, &neighborNode2->myTag,
                      -burgNew[X], -burgNew[Y], -burgNew[Z],
                      planeSeg2[X], planeSeg2[Y], planeSeg2[Z], globalOp);
                      
            InsertArm(home, neighborNode2, &newNode->myTag, 
                      burgNew[X], burgNew[Y], burgNew[Z],
                      planeSeg2[X], planeSeg2[Y], planeSeg2[Z], globalOp);
                      
/*
 *          The topology changes will significantly affect the force
 *          and velocity of the four nodes, so flag them to indicate
 *          they need new force/vel calcs at the beginning of the next
 *          timestep.  NOTE: it would be unwise to do any further
 *          topological changes to these nodes this cycle if the
 *          changes require valid nodal force/vel values.
 */
            MarkNodeForceObsolete(home, node);
            MarkNodeForceObsolete(home, newNode);
            MarkNodeForceObsolete(home, neighborNode1);
            MarkNodeForceObsolete(home, neighborNode2);

#ifdef DEBUG_TOPOLOGY_CHANGES
            if ((dbgDom < 0) || (dbgDom == home->myDomain)) {
                printf("Splinter: (%d,%d) into (%d,%d), old burg %e %e %e "
                       "new burg %e %e %e\n",
                       node->myTag.domainID, node->myTag.index,
                       newNode->myTag.domainID, newNode->myTag.index,
                       burgSeg1[X], burgSeg1[Y], burgSeg1[Z],
                       burgNew[X], burgNew[Y], burgNew[Z]);
                PrintNode(node);
                PrintNode(newNode);
                PrintNode(neighborNode1);
                PrintNode(neighborNode2);
            }
#endif

/*
 *          The following checks have been added for debugging only
 *          and should be removed once the splintering code has been verified
 *          to perform properly without resulting in inconsistencies when
 *          used in conjunction with the other topology changing functions.
 */
            if ((CheckNodeBurgConservation(home, node)    == 0)       ||
                (CheckNodeBurgConservation(home, newNode) == 0)       ||
                (CheckNodeBurgConservation(home, neighborNode1) == 0) ||
                (CheckNodeBurgConservation(home, neighborNode2) == 0)) {
                printf("Splinter: (%d,%d) into (%d,%d), old burg %e %e %e "
                       "new burg %e %e %e\n",
                       node->myTag.domainID, node->myTag.index,
                       newNode->myTag.domainID, newNode->myTag.index,
                       burgSeg1[X], burgSeg1[Y], burgSeg1[Z],
                       burgNew[X], burgNew[Y], burgNew[Z]);
                PrintNode(node);
                PrintNode(newNode);
                PrintNode(neighborNode1);
                PrintNode(neighborNode2);
                Fatal("Burgers vector not conserved after splinter operation");
            }

        }  /* end for (nodeIndex = 0; ...) */

        return;
}
