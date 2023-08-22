/*****************************************************************************
 *
 *	Module:		SplitSurfaceNodes.c
 *	Description:	This module contains functions needed to split
 *                      (if necessary) a surface node into multiple nodes.
 *                      All segment and node ownership rules are the
 *                      same as for other operations (see Topology.c).
 *
 *	Included functions:
 *          GetSplitSets()
 *          SplitSurfaceNodes()
 *          TestSurfaceSplit()
 *
 *****************************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <sys/param.h>

#include "mpi_portability.h"

#include "Home.h"
#include "Util.h"


#ifdef _ARLFEM

#define MAX_ARMS 9

static int dbgDom;

/*
 *      Predefine the number of ways it is possible to split
 *      apart a surface node of up to 9 segments.
 *
 *      Some ways of splitting a node may be logically equivalent
 *      to other ways (i.e. for a 5-node, splitting off arms
 *      1 & 2 is logically equivalent to spliting off arms
 *      3, 4 & 5).  At times we want to treat these logically
 *      eqivalent cases as if they are exactly the same, but in
 *      other cases we need to treat them as if they are in fact
 *      unique.  Hence, we have two sets of counts of possible
 *      splits.  The first <maxSplitsEquiv> assumes equivalent
 *      cases are combined as a single possibility, while the
 *      other <maxSplitsNoEquiv> treats every case as unique.
 */
int maxSplitsEquiv[MAX_ARMS+1]   = { 0, 0, 2, 4, 8,  16, 32, 64,  128, 256 };
int maxSplitsNoEquiv[MAX_ARMS+1] = { 0, 0, 3, 7, 15, 31, 63, 127, 255, 511 };


/**************************************************************************
 *
 *      Function:       GetSplitSet
 *      Description:    Build a set of possible ways to split a surface
 *                      node with the specified number of segments. 
 *
 *      Arguments:
 *          nodeSegCnt   Number of segments attached to the node we
 *                       want to split.
 *          includeDups  If set to a non-zero value, the returned list
 *                       will contain *all* split possibilities, including
 *                       splits that may be logically equivalent to others
 *                       in the list (See comments by <maxSplitsEquiv>
 *                       and <maxSplitsNoEquiv> above.
 *          splitSets    Multidimensional array in which to store the
 *                       set of ways to split the node. 
 *
 * FIX ME! Describe contents of <splitSets>
 *
 *************************************************************************/
static void GetSplitSets(int nodeSegCnt, int includeDups,
                         int splitSets[][MAX_ARMS+1])
{
        int  i, index, num2Split, setCnt, done;
        int  minIndex, minSplitNum;
        int  segID[MAX_ARMS];


        minSplitNum = includeDups ? 1 : (nodeSegCnt+1) / 2;

/*
 *      Loop though the possible number of segments to split
 *      (i.e 2 <= num2Split <= nodeSegCnt)
 */
        setCnt = 0;

        for (num2Split = nodeSegCnt; num2Split >= minSplitNum; num2Split--) {

/*
 *          We'll use segID[0] thru segID[num2Split-1] to hold the
 *          list of segments to try to split from the original node.
 */
            for (index = 0; index < num2Split; index++) {
                segID[index] = index;
            }

            done     = 0;
            index    = num2Split - 1;

            if (includeDups) {
                minIndex = 0;
            } else {
                minIndex = (num2Split << 1 == nodeSegCnt) ? 1 : 0;
                minIndex = MIN(index, minIndex);
            }

/*
 *          The initialization of <segID[]> above sets the array to
 *          the set of IDs for the first valid group of segments that
 *          can be split off.  We won't exit the loop until we've 
 *          cycled through all possible valid sets.
 *
 *          For example, if we're dealing with a 5-node and we're on
 *          the loop where num2Split is 3, this loop would in this order 
 *          find the following sets of segment IDs:
 *
 *          segID[0]      segID[1]      segID[2]    segID[3] thru segID[N]
 *              0            1             2          Not Applicable
 *              0            1             3          Not Applicable
 *              0            1             4          Not Applicable
 *              0            2             3          Not Applicable
 *              0            2             4          Not Applicable
 *              0            3             4          Not Applicable
 *              1            2             3          Not Applicable
 *              1            2             4          Not Applicable
 *              1            3             4          Not Applicable
 *              2            3             4          Not Applicable
 */
            while (!done) {
/*
 *              If segID[index] has been incremented past the final segment
 *              ID that is valid (for this position of the array), we need
 *              to move back one element in the segID array, increment the
 *              segment ID at that position, and continue the loop from
 *              there.
 */
                if (((includeDups == 0) && (index == 0) &&
                     (num2Split << 1 == nodeSegCnt) && (segID[index] > 0)) ||
                    (segID[index] > nodeSegCnt - (num2Split - index))) {
/*
 *                  If we've exhausted all possible values for the first seg
 *                  ID, we've been through all possible ways to split this
 *                  number of segments from the node, so break out of the
 *                  while loop. 
 */
                    if (index == minIndex) {
                        done = 1;
                    } else {
                        segID[--index] += 1;
                    }

                } else if (index < (num2Split-1)) {
/*
 *                  We have a valid segment ID in the current array position
 *                  but we're not at the end of the array so move to the next
 *                  position and set the segment ID for that one to 
 *                  appropriate value and continue
 */
                    segID[index+1] = segID[index]+1;
                    index++;

                } else {
/*
 *                  segID[0] thru segID[num2Split-1] now contains a list of
 *                  segment IDs that could be split from the original node
 */

                    for (i = 0; i < num2Split; i++) {
                        splitSets[setCnt][i] = segID[i];
                    }
                    splitSets[setCnt][MAX_ARMS] = num2Split;

                    setCnt++;
/*
 *                  Bump up the segment ID in the current position for the
 *                  next loop iteration.
 */
                    segID[index] += 1;
                }

            }  /* while (1) */
        }

        return;
}


/**************************************************************************
 *
 *      Function:       TestSurfaceSplit
 *      Description:    Attempt to split the given segments of a surface
 *                      node to a new node, evaluate the energy release
 *                      from the operation, then merge the split nodes
 *                      and restore the incoming node to its original state.
 *
 *      Arguments:
 *          node         Pointer to node to be split
 *          segCnt       Number of segments to be split from the node.
 *          segList      Array containing indices of segments to be split from
 *                       the node.  Array is at least <segCnt> elements long.
 *          splitFlags   Additional flags to pass to SplitNode().
 *          bkupNodeCnt  Number of nodes in the <bkupNodeList> array.
 *          bkupNodeList Array in which the original node (and neighbor
 *                       node) data has been preserved.
 *          v1           Array in which to preserve velocity of first
 *                       node if this split is the optimal split.
 *          v2           Array in which to preserve velocity of second
 *                       node if this split is the optimal split.
 *          powerMax     Pointer to Baseline energy release used to determine
 *                       if the node-split is optimal.  If this split results
 *                       in a higher energy release, <powerMax> will be 
 *                       updated for the caller.
 *          errStatus    Location in which to return an error status to the
 *                       caller.  0 == no error, non-zero indicates the
 *                       split could not be done.
 *          isSelecetd   Location of flag returned to caller.  Set to 1 if
 *                       this split is the optimal split, 0 otherwise.
 *
 *************************************************************************/
static void TestSurfaceSplit(Home_t *home, Node_t *node, int segCnt,
                             int *segList, int splitFlags, int bkupNodeCount,
                             Node_t *bkupNodeList, real8 v1[3], real8 v2[3],
                             real8 *powerMax, int *errStatus, int *isSelected)
{
        int     i, armIndex;
        int     splitStatus, mergeStatus, mobError;
        real8   eps, splitDist, minDist, powerTest, vNoise;
        real8   vd1, vd2;
        real8   dirx, diry, dirz, invvNorm;
        real8   vDiffx, vDiffy, vDiffz;
        real8   origPos[3], origVel[3];
        Node_t  *splitNode1, *splitNode2;
        Node_t  *origNode, *newNode, *mergedNode, *tmpNode;
        Param_t *param;
        MobArgs_t mobArgs;


        param = home->param;

        vNoise   = param->splitMultiNodeAlpha * param->rTol / param->deltaTT;
        minDist  = param->rann * 2.0;

        splitDist = minDist;
        eps = 1.0e-12;

        *errStatus = 0;
        *isSelected = 0;

        origPos[X] = node->x;
        origPos[Y] = node->y;
        origPos[Z] = node->z;

        origVel[X] = node->vX;
        origVel[Y] = node->vY;
        origVel[Z] = node->vZ;

/*
 *      Note: for this test, both the node being split and the resulting
 *      new node will intially retain the location and velocity of the
 *      node being split.  This may lead to a zero length segment between
 *      the two nodes, but the mobility function should be okay with that.
 *      Also, this operation MUST NOT be communicated to remote domains
 *      for processing!
 *
 *      Note: on success, the arms selected to be split
 *      will be connected to splitNode2 and all unselected
 *      arms will be attached to splitNode1.
 */
        splitStatus = SplitNode(home, OPCLASS_SEPARATION, node, origPos,
                                origPos, origVel, origVel, segCnt,
                                segList, 0, &splitNode1, &splitNode2,
                                splitFlags);

/*
 *      If the split could not be done, return now.
 */
        if (splitStatus != SPLIT_SUCCESS) {
            *errStatus = 1;
            return;
        }

#if 0
/*
 *      Calculate forces and velocities for the two nodes after the split
 */
        SetOneNodeForce(home, splitNode1);
        SetOneNodeForce(home, splitNode2);
#endif
        GetForcesFromBkup(home, splitNode1, bkupNodeList);
        GetForcesFromBkup(home, splitNode2, bkupNodeList);

/*
 *      If we're using special mobility for segments intersecting eshelby
 *      inclusions, we need to generate the segment/particle intersection list
 *      (and clear it after calculating velocities)
 */
        FindNodePartIntersects(home, splitNode1);
        FindNodePartIntersects(home, splitNode2);
        SegPartListSort(home);

        mobError  = EvaluateMobility(home, splitNode1, &mobArgs);

        mobError |= EvaluateMobility(home, splitNode2, &mobArgs);

        SegPartListClear(home);

/*
 *      After the split is done, reposition the nodes and again calculate
 *      forces/velocities.  If the direction of velocity changes, don't
 *      do the split, but only if we didn't have any errors returned from
 *      the mobility function.
 */
        if (mobError == 0) {
            vd1 = splitNode1->vX*splitNode1->vX +
                  splitNode1->vY*splitNode1->vY +
                  splitNode1->vZ*splitNode1->vZ;

            vd2 = splitNode2->vX*splitNode2->vX +
                  splitNode2->vY*splitNode2->vY +
                  splitNode2->vZ*splitNode2->vZ;

/*
 *          If, after the split, the nodes have zero velocities, this
 *          configuration should not be considered as an option.
 */
            if ((fabs(vd1) > eps) || (fabs(vd2) > eps)) {

                if (vd1 > vd2) {

                    invvNorm = 1.0 / sqrt(vd1);

                    dirx = -splitNode1->vX * invvNorm;
                    diry = -splitNode1->vY * invvNorm;
                    dirz = -splitNode1->vZ * invvNorm;

                    splitNode1->x -= (splitDist*(1.0+eps)*dirx);
                    splitNode1->y -= (splitDist*(1.0+eps)*diry);
                    splitNode1->z -= (splitDist*(1.0+eps)*dirz);

                    FoldBox(param, &splitNode1->x, &splitNode1->y,
                            &splitNode1->z);

                } else {

                    invvNorm = 1.0 / sqrt(vd2);

                    dirx = splitNode2->vX * invvNorm;
                    diry = splitNode2->vY * invvNorm;
                    dirz = splitNode2->vZ * invvNorm;

                    splitNode2->x += (splitDist*(1.0+eps)*dirx);
                    splitNode2->y += (splitDist*(1.0+eps)*diry);
                    splitNode2->z += (splitDist*(1.0+eps)*dirz);

                    FoldBox(param, &splitNode2->x, &splitNode2->y,
                            &splitNode2->z);
                }

                SetOneNodeForce(home, splitNode1);
                SetOneNodeForce(home, splitNode2);

/*
 *              When the original node was split above, both new nodes were
 *              left at the original position.  If a new segment was
 *              created connecting the nodes, it is intially zero length and
 *              therefore had no glide-plane, but the mobility function
 *              can handle zero-length segs.  After calculating velocities,
 *              though we've shifted the nodes making the connecting segment
 *              (if any) non-zero length, so now we have to explicitly set
 *              the glide plane of that seg because some of the mobility
 *              functions depend on the glide plane being set for non-zero
 *              length segs
 *
 *              From Wei Cai, 2009/09/07
 */
                armIndex = GetArmID(splitNode1, splitNode2);

                if (armIndex >= 0 ) {
                    real8 burg[3], dirvec[3], nplane[3];

/*
 *                  Calculate glide plane normal as
 *                  n = cross ( burg(x,y,z), dir(x,y,z))
 */
                    burg[X] = node->burgX[armIndex];
                    burg[Y] = node->burgY[armIndex];
                    burg[Z] = node->burgZ[armIndex];

                    dirvec[X] = dirx;
                    dirvec[Y] = diry;
                    dirvec[Z] = dirz;

                    FindPreciseGlidePlane(home, burg, dirvec, nplane,
                                          home->param->allowFuzzyGlidePlanes);

/*
 *                  If the new segment is screw the glide plane is
 *                  undefined.  If that's the case, we need to randomly
 *                  select a glide plane appropriate to the burgers vector.
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

/*
 *              In case eshelby inclusions are being simulated, the
 *              segment/particle intersection list may have been generated
 *              during force calcs. We need to sort the list now, and clear
 *              it when we're done with it.
 */
                SegPartListSort(home);

                mobError  = EvaluateMobility(home, splitNode1, &mobArgs);

                mobError |= EvaluateMobility(home, splitNode2, &mobArgs);

                SegPartListClear(home);

                vDiffx = splitNode2->vX - splitNode1->vX;
                vDiffy = splitNode2->vY - splitNode1->vY;
                vDiffz = splitNode2->vZ - splitNode1->vZ;

            } else {

                vDiffx = -1.0;
                vDiffy = -1.0;
                vDiffz = -1.0;
                dirx = 1.0;
                diry = 1.0;
                dirz = 1.0;
            }
        } /* if (mobError == 0) */

        if (mobError != 0) {
/*
 *          If the mobility function was unable to converge on
 *          velocities for either node in their new positions,
 *          we can't really evaluate the node split properly.
 *          Given that, regardless of what other split
 *          possibilities were evaluated for this multi-node,
 *          set flags to avoid splitting this node and skip
 *          evaluation of any other possible splits for this
 *          node this cycle.
 */
            *errStatus = 1;

        } else if (((vDiffx*dirx) + (vDiffy*diry) + (vDiffz*dirz)) > 0) {

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
 *          If this potential node split would result in the highest
 *          energey release so far, save a little info in case we want
 *          to perform this split later.
 */
            if ((powerTest - (*powerMax)) > eps) {

                v1[X] = splitNode1->vX;
                v1[Y] = splitNode1->vY;
                v1[Z] = splitNode1->vZ;

                v2[X] = splitNode2->vX;
                v2[Y] = splitNode2->vY;
                v2[Z] = splitNode2->vZ;

                *isSelected = 1;
                *powerMax = powerTest;

            }
        }

/*
 *      Merge the new node back into the original node
 *      then restore the original node and it's neighbors
 *      to the pristine unsplit state.  (Need to do a 
 *      restore to insure the ordering of arms since
 *      the arrays of splitting possibilites use the
 *      original node's arm indices, plus we get back
 *      all the original force/velocity/etc values.)
 *
 *      Assumes that MergeNode() will merge <newNode>
 *      into <origNode>!
 */
        if (node == splitNode1) {
            origNode = splitNode1;
            newNode  = splitNode2;
        } else {
            newNode  = splitNode1;
            origNode = splitNode2;
        }

/*
 *      Since SplitNode() can play tricks with which node it
 *      attaches the split segments to and switch the  surface properties
 *      from the original node to the new node, we need to add a little
 *      kludge here to be certain we can merge the new node back into
 *      the original and still preserve the proper surface constraint.
 */
        if (HAS_NONE_OF_CONSTRAINTS(origNode->constraint, SURFACE_NODE)) {

            origNode->x = origPos[X];
            origNode->y = origPos[Y];
            origNode->z = origPos[Z];

            ADD_CONSTRAINTS(origNode->constraint, SURFACE_NODE);

            origNode->surfaceNorm[0] = newNode->surfaceNorm[0];
            origNode->surfaceNorm[1] = newNode->surfaceNorm[1];
            origNode->surfaceNorm[2] = newNode->surfaceNorm[2];
        }

        MergeNode(home, OPCLASS_SEPARATION, newNode, origNode,
                  origPos, &mergedNode, &mergeStatus, 0);

/*
 *      If the merge failed, something's screwy!  Abort.
 */
        if ((mergeStatus & MERGE_SUCCESS) == 0) {
            Fatal("TestSurfaceSplit: Unable to merge (%d,%d) into (%d,%d)"
                  "during separation test", newNode->myTag.domainID,
                  newNode->myTag.index, origNode->myTag.domainID,
                  origNode->myTag.index);
        }

/*
 *      Restore nodes to the original state
 */
        for (i = 0; i < bkupNodeCount; i++) {
            tmpNode = GetNodeFromTag(home, bkupNodeList[i].myTag);
            RestoreNode(home, tmpNode, &bkupNodeList[i]);
        }

        return;
}


/*---------------------------------------------------------------------------
 *
 *      Function:     SplitSurfaceNodes
 *
 *      Description:  Identify all surface nodes that are candidates to
 *                    be split into multiple nodes, evaluate the possible
 *                    ways to split such nodes, and split the nodes if
 *                    necessary.  Similar in function to the
 *                    SplitMultiNodes() routine in Topology.c
 *
 *      Note: Since this function performs operations that must be
 *            distributed to remote domains, it should be preceeded by
 *            a call to ClearOpList() and followed by calls to
 *            CommSendRemesh() and FixRemesh().
 *
 *-------------------------------------------------------------------------*/
void SplitSurfaceNodes(Home_t *home)
{
        int     i, j, k;
        int     doSplit, skipSplit, splitStatus, selectedFlags;
        int     isSelected, errStatus, splitFlags;
        int     nbrs, bkupNodeCount, globalOp, tmpArm;
        int     numSets, setListSize;
        int     selectedSet[MAX_ARMS+1];
        int     (*armSet)[MAX_ARMS+1];
        real8   eps;
        real8   invvNorm, splitDist, minDist;
        real8   dirx, diry, dirz;
        real8   vd1, vd2, powerMax;
        real8   v1[3], v2[3];
        real8   pos1[3], pos2[3];
        Node_t  *node, *nbrNode, *origNode, *newNode;
        Node_t  *splitNode1, *splitNode2;
        Node_t  *bkupNodeList;
        Param_t *param;


/*
 *      Only do surface node splitting on cycles that are multiples
 *      of the <splitMultiNodeFreq> control file parameter.
 */
        if (home->cycle % home->param->splitMultiNodeFreq) {
            return;
        }

#ifdef DEBUG_TOPOLOGY_DOMAIN
        dbgDom   = DEBUG_TOPOLOGY_DOMAIN;
#else
        dbgDom   = -1;
#endif
        param    = home->param;
        minDist  = param->rann * 2.0;

        eps      = 1.0e-12;
        globalOp = 1;
        splitDist = minDist;

/*
 *      In case we're simulating eshelby inclusions and using a different
 *      mobility for segments intersecting them, we need to re-initialize
 *      the segment/particle intersection list here.  Any other portion of
 *      code using the list *should* have cleared it when the list was no
 *      longer needed, but just for safety's sake we do it here.
 */
        SegPartListClear(home);

/*
 *      Loop through all native nodes; any surface node with at least 2 arms
 *      of differing burgers vectors is a candidate for being split.
 */
        for (i = 0; i < home->newNodeKeyPtr; i++) {

            if ((node = home->nodeKeys[i]) == (Node_t *)NULL) {
                continue;
            }

            if ((node->numNbrs < 2) ||
                HAS_NONE_OF_CONSTRAINTS(node->constraint, SURFACE_NODE))
            {
                continue;
            }

            if ((node->numNbrs == 2) &&
                ((fabs(node->burgX[0] + node->burgX[1]) < 1.0e-6) &&
                 (fabs(node->burgY[0] + node->burgY[1]) < 1.0e-6) &&
                 (fabs(node->burgZ[0] + node->burgZ[1]) < 1.0e-6))) {
                continue;
            }

/*
 *          We've got a hard-coded limit on the maximum number of arms
 *          we can deal with... if that's exceeded, abort with an error
 */
            nbrs = node->numNbrs;

            if (nbrs > MAX_ARMS) {
                Fatal("SplitSurfaceNodes: Node (%d,%d) has %d arms and "
                      "cannot be split!", node->myTag.domainID,
                       node->myTag.index, nbrs);
            }

/*
 *          Preserve the original state of the node being evaluated
 *          for separation as well as the state of all its neighbor
 *          nodes.
 */
            bkupNodeCount = node->numNbrs + 1;
            bkupNodeList = (Node_t *)malloc(bkupNodeCount * sizeof(Node_t));

            BackupNode(home, node, &bkupNodeList[0]);

            for (j = 0; j < node->numNbrs; j++) {
                nbrNode = GetNodeFromTag(home, node->nbrTag[j]);
                BackupNode(home, nbrNode, &bkupNodeList[j+1]);
            }

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
 *          For each of the node-splitting possibilities, attempt
 *          the node split, evaluate the forces and velocities of
 *          the two nodes after the split, and determine which 
 *          splitting possibility leads to the greatest energy
 *          release; that's the split that should ultimately be
 *          done.
 *
 *          First test the unique split possibilities constraining the
 *          new node to the surface
 */
            splitFlags = SPLIT_DUP_SURFACE_PROP;

            numSets  = maxSplitsEquiv[nbrs];
            setListSize = (MAX_ARMS+1) * numSets * sizeof(int);
            armSet   = (int (*)[MAX_ARMS+1])calloc(1, setListSize);

            GetSplitSets(nbrs, 0, armSet);

            doSplit = 0;
            skipSplit = 0;
            selectedFlags = 0;

            for (j = 0; j < numSets; j++) {

                TestSurfaceSplit(home, node, armSet[j][MAX_ARMS],
                                 armSet[j], splitFlags, bkupNodeCount,
                                 bkupNodeList, v1, v2, &powerMax, &errStatus,
                                 &isSelected);

                if (errStatus) {
                    skipSplit = 1;
                    doSplit = 0;
                    break;
                }

                if (isSelected) {
                    doSplit = 1;
                    selectedFlags = splitFlags;
                    for (k = 0; k <= MAX_ARMS; k++) {
                        selectedSet[k] = armSet[j][k];
                    }
                }
            }

            free(armSet);

/*
 *          Next test all possible splits without constraining the new node
 *          to the surface
 */
            splitFlags = 0;

            numSets  = maxSplitsNoEquiv[nbrs];
            setListSize = (MAX_ARMS+1) * numSets * sizeof(int);
            armSet  = (int (*)[MAX_ARMS+1])calloc(1, setListSize);

            GetSplitSets(nbrs, 1, armSet);

            for (j = 0; (j < numSets) && (skipSplit == 0); j++) {

                TestSurfaceSplit(home, node, armSet[j][MAX_ARMS], 
                                 armSet[j], splitFlags, bkupNodeCount,
                                 bkupNodeList, v1, v2, &powerMax, &errStatus,
                                 &isSelected);

                if (errStatus) {
                    skipSplit = 1;
                    doSplit = 0;
                    break;
                }

                if (isSelected) {
                    doSplit = 1;
                    selectedFlags = splitFlags;
                    for (k = 0; k <= MAX_ARMS; k++) {
                        selectedSet[k] = armSet[j][k];
                    }
                }
            }

            free(armSet);
            armSet = NULL;

/*
 *          If we've determined the node should be split, perform
 *          the split in the optimal manner.
 */
            if (doSplit) {
/*
 *              The node with the lower velocity will remain
 *              at the original location and node with the
 *              higher velocity will be moved.
 */
                vd1 = v1[X] * v1[X] +
                      v1[Y] * v1[Y] +
                      v1[Z] * v1[Z];

                vd2 = v2[X] * v2[X] +
                      v2[Y] * v2[Y] +
                      v2[Z] * v2[Z];

                if (vd1 > vd2) {
                    invvNorm = 1.0 / sqrt(vd1);
                    dirx = v1[X] * invvNorm;
                    diry = v1[Y] * invvNorm;
                    dirz = v1[Z] * invvNorm;
                    pos2[X] = node->x;
                    pos2[Y] = node->y;
                    pos2[Z] = node->z;
                    pos1[X] = node->x+(splitDist*(1.0+eps)*dirx);
                    pos1[Y] = node->y+(splitDist*(1.0+eps)*diry);
                    pos1[Z] = node->z+(splitDist*(1.0+eps)*dirz);
                    FoldBox(param, &pos1[X], &pos1[Y], &pos1[Z]);
                } else {
                    invvNorm = 1.0 / sqrt(vd2);
                    dirx = v2[X] * invvNorm;
                    diry = v2[Y] * invvNorm;
                    dirz = v2[Z] * invvNorm;
                    pos1[X] = node->x;
                    pos1[Y] = node->y;
                    pos1[Z] = node->z;
                    pos2[X] = node->x+(splitDist*(1.0+eps)*dirx);
                    pos2[Y] = node->y+(splitDist*(1.0+eps)*diry);
                    pos2[Z] = node->z+(splitDist*(1.0+eps)*dirz);
                    FoldBox(param, &pos2[X], &pos2[Y], &pos2[Z]);
                }

/*
 *              Split the node. This operation is always considered
 *              *global* so other domains will be notified of this.
 *
 *              Note: after the split, the arms selected to be split
 *              will be connected to splitNode2 which will be
 *              positioned at pos2 with velocity vel2.  Node
 *              splitNode1 will be attached to the arms not selected
 *              to be moved, and will have the position and velocity
 *              indicated by pos1 and vel1.
 */
                splitStatus = SplitNode(home, OPCLASS_SEPARATION,
                                        node, pos1, pos2, v1, v2,
                                        selectedSet[MAX_ARMS],
                                        selectedSet, globalOp,
                                        &splitNode1, &splitNode2,
                                        selectedFlags);

/*
 *              If the split failed, clean up then go back and
 *              continue looking for more multinodes to split.
 */
                if (splitStatus != SPLIT_SUCCESS) {
                    for (k = 0; k < bkupNodeCount; k++) {
                        FreeNodeArrays(&bkupNodeList[k]);
                    }

                    free(bkupNodeList);
                    bkupNodeList = (Node_t *)NULL;

                    continue;
                }

/*
 *              Velocities passed to split should be okay, but we need to
 *              set the force of the nodes involved in the split, then call
 *              ResetSegForces2() to distribute the new forces out to the
 *              remotedomains.
 */
/*
 * FIX ME!  We should change this code to save forces during the split testing
 *          and use them here rather than recalculate them here!
 */
                SetOneNodeForce(home, splitNode1);
                SetOneNodeForce(home, splitNode2);

                SegPartListClear(home);

                if (node == splitNode1) {
                    origNode = splitNode1;
                    newNode  = splitNode2;
                } else {
                    newNode  = splitNode1;
                    origNode = splitNode2;
                }

                for (j = 0; j < origNode->numNbrs; j++) {
                    nbrNode = GetNodeFromTag(home, origNode->nbrTag[j]);
                    tmpArm = GetArmID(nbrNode, origNode);
                    ResetSegForces2(home, origNode, &nbrNode->myTag,
                                    origNode->armfx[j],
                                    origNode->armfy[j],
                                    origNode->armfz[j],
                                    nbrNode->armfx[tmpArm],
                                    nbrNode->armfy[tmpArm],
                                    nbrNode->armfz[tmpArm], globalOp);
                }

                for (j = 0; j < newNode->numNbrs; j++) {
                    nbrNode = GetNodeFromTag(home, newNode->nbrTag[j]);
                    tmpArm = GetArmID(nbrNode, newNode);
                    ResetSegForces2(home, newNode, &nbrNode->myTag,
                                    newNode->armfx[j],
                                    newNode->armfy[j],
                                    newNode->armfz[j],
                                    nbrNode->armfx[tmpArm],
                                    nbrNode->armfy[tmpArm],
                                    nbrNode->armfz[tmpArm], globalOp);
                }

/*
 *              Mark both nodes involved in the split as 'exempt'
 *              from subsequent collisions this time step. 
 */
                origNode->flags |= NO_COLLISIONS;
                newNode->flags  |= NO_COLLISIONS;

#ifdef DEBUG_TOPOLOGY_CHANGES
/*
 *              When debugging, dump some info on topological
 *              changes taking place and the nodes involved
 */
                if ((dbgDom < 0)||(dbgDom == home->myDomain)) {
                    printf("  SurfaceSplit: (%d,%d) from (%d,%d)\n",
                           newNode->myTag.domainID, newNode->myTag.index,
                           origNode->myTag.domainID, origNode->myTag.index);
                    PrintNode(origNode);
                    PrintNode(newNode);
                }
#endif
            }  /* if (doSplit) */

/*
 *          Need to free the arrays and node structures used
 *          above for backup/restore or we end up with a memory
 *          leak.
 */
            for (k = 0; k < bkupNodeCount; k++) {
                FreeNodeArrays(&bkupNodeList[k]);
            }

            free(bkupNodeList);
            bkupNodeList = (Node_t *)NULL;

        }  /* loop over nodes */

        return;
}
#endif  /* ifdef _ARLFEM */
